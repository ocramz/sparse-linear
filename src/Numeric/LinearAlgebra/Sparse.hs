{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse where

import Control.Monad (when)
import Control.Monad.ST (ST, runST)
import Data.Function (fix)
import Data.Ord (comparing)
import qualified Data.Vector.Algorithms.Intro as Intro
import Data.Vector.Unboxed (Unbox, Vector)
import qualified Data.Vector.Unboxed as V
import Data.Vector.Unboxed.Mutable (STVector)
import qualified Data.Vector.Unboxed.Mutable as MV
import GHC.Stack (errorWithStackTrace)

data Orient = Row | Col

data Matrix (or :: Orient) a = Matrix
    { nmajor :: {-# UNPACK #-} !Int -- positive
    , nminor :: {-# UNPACK #-} !Int -- positive
    , majors :: !(Vector Int) -- length == nmajor + 1; head == 0
    , minors :: !(Vector Int) -- length == nz; all (< nminor)
    , entries :: !(Vector a) -- length == nz
    }

nnz :: Unbox a => Matrix or a -> Int
nnz = V.length . entries
{-# INLINE nnz #-}

class OrientOf or where
    orientOf :: Matrix or a -> Orient

instance OrientOf Row where
    orientOf = \_ -> Row
    {-# INLINE orientOf #-}

instance OrientOf Col where
    orientOf = \_ -> Col
    {-# INLINE orientOf #-}

dedup :: (Num a, Unbox a) => STVector s (Int, Int, a) -> ST s ()
dedup v = do
    let len = MV.length v
        loop prevR prevC prevA prevN nextN
          | nextN < len = do
              (r, c, a) <- MV.read v nextN
              if r == prevR && c == prevC
                then loop r c (a + prevA) prevN (nextN + 1)
                else do
                    MV.write v prevN (prevR, prevC, prevA)
                    loop r c a (prevN + 1) (nextN + 1)
          | otherwise = return ()
    if len > 0
      then do
          (r0, c0, a0) <- MV.read v 0
          loop r0 c0 a0 0 1
      else return ()
{-# INLINE dedup #-}

compress
  :: (OrientOf or, Num a, Unbox a)
  => Int
  -> Int
  -> Vector (Int, Int, a)
  -> Matrix or a
compress = \nr nc coords -> fix $ \compressed ->
    let orient = orientOf compressed
        nmajor = case orient of { Row -> nr; Col -> nc }
        nminor = case orient of { Row -> nc; Col -> nr }

        comparison = case orient of
          Row -> comparing $ \(r, c, _) -> (r, c)
          Col -> comparing $ \(r, c, _) -> (c, r)
        sorted = V.modify (\v -> do { Intro.sortBy comparison v; dedup v }) coords
        (rows, cols, entries) = V.unzip3 sorted

        minors = case orient of { Row -> cols; Col -> rows }
        count ixs = V.create $ do
            v <- MV.replicate 0 (nmajor + 1)
            V.forM_ ixs $ \ix ->
                if ix >= nmajor
                  then errorWithStackTrace "major index out of bounds"
                  else MV.write v (ix + 1) . (+ 1) =<< MV.read v (ix + 1)
            return v
        majors = case orient of
          Row -> count rows
          Col -> count cols
    in Matrix{..}
{-# INLINE compress #-}

slice :: Unbox a => Int -> Matrix or a -> Vector (Int, a)
slice n Matrix{..}
  | n >= nmajor = errorWithStackTrace "major index out of bounds"
  | otherwise =
      let start = majors V.! n
          end = majors V.! (n + 1)
          len = end - start
          minors' = V.slice start len minors
          entries' = V.slice start len entries
      in V.zip minors' entries'
{-# INLINE slice #-}

dim :: OrientOf or => Matrix or a -> Int
dim mat@Matrix{..} =
    case orientOf mat of
      Row -> nmajor
      Col -> nminor
{-# INLINE dim #-}

codim :: OrientOf or => Matrix or a -> Int
codim mat@Matrix{..} =
    case orientOf mat of
      Col -> nmajor
      Row -> nminor
{-# INLINE codim #-}

forI :: Monad m => Int -> (Int -> m ()) -> m ()
forI end act = forI_go 0 where
  forI_go n
    | n < end = act n >> forI_go (n + 1)
    | otherwise = return ()
{-# INLINE forI #-}

zipWithI
  :: (Unbox a, Unbox b, Unbox c)
  => (a -> b -> c) -> Vector (Int, a) -> Vector b -> Vector c
zipWithI f ias bs = V.map (\(i, a) -> f a (bs V.! i)) ias
{-# INLINE zipWithI #-}

mulV
  :: (Num a, Unbox a)
  => Matrix Row a
  -> Vector a
  -> Vector a
mulV mat@Matrix{..} xs
  | V.length xs /= dim mat =
      errorWithStackTrace "matrix dimension does not match vector"
  | otherwise = V.create $ do
      ys <- MV.replicate (V.length xs) 0
      forI nmajor $ \r -> do
          let ax = V.sum $ zipWithI (*) (slice r mat) xs
          MV.write ys r . (+ ax) =<< MV.read ys r
      return ys
{-# INLINE mulV #-}

fromDiag
  :: (Unbox a)
  => Vector a
  -> Matrix or a
fromDiag = \entries ->
    let nmajor = V.length entries
        nminor = nmajor
        majors = V.enumFromN 0 (nmajor + 1)
        minors = V.enumFromN 0 nminor
    in Matrix{..}
{-# INLINE fromDiag #-}

scale :: (Num a, Unbox a) => a -> Matrix or a -> Matrix or a
scale = \sc mat -> mat { entries = V.map (* sc) $ entries mat }
{-# INLINE scale #-}

add :: (Num a, Unbox a) => Matrix or a -> Matrix or a -> Matrix or a
add matA matB
  | nmajor matA /= nmajor matB =
      errorWithStackTrace "matrix major dimensions do not agree"
  | nminor matA /= nminor matB =
      errorWithStackTrace "matrix minor dimensions do not agree"
  | nmajor matA <= 0 =
      errorWithStackTrace "matrix major dimension must be positive"
  | nminor matA <= 0 =
      errorWithStackTrace "matrix minor dimension must be positive"
  | otherwise = runST $ do
      let ixls = minors matA
          xls = entries matA
          ixrs = minors matB
          xrs = entries matB
          next (j, _, ia, iaEnd, ib, ibEnd, _) =
              let ixl = ixls V.! ia
                  ixr = ixrs V.! ib
                  xl = xls V.! ia
                  xr = xrs V.! ib
                  (ia', ib', entry) = case compare ixl ixr of
                    EQ -> (ia + 1, ib + 1, (ixl, xl + xr))
                    LT -> (ia + 1, ib, (ixl, xl))
                    GT -> (ia, ib + 1, (ixr, xr))
              in case (ia < iaEnd, ib < ibEnd) of
                (True, True) ->
                    Just (j, False, ia', iaEnd, ib', ibEnd, entry)
                (False, True) ->
                    Just (j, False, ia, iaEnd, ib + 1, ibEnd, (ixr, xr))
                (True, False) ->
                    Just (j, False, ia + 1, iaEnd, ib, ibEnd, (ixl, xl))
                (False, False) ->
                    let j' = j + 1
                        iaEnd' = (majors matA) V.! (j' + 1)
                        ibEnd' = (majors matB) V.! (j' + 1)
                    in if j' < nmajor matA
                      then Just (j', True, ia', iaEnd', ib', ibEnd', entry)
                      else Nothing
          nnzUB = nnz matA + nnz matB
          nmajor_ = nmajor matA
          nminor_ = nminor matA
      majors__ <- MV.new (nmajor_ + 1)
      minors__ <- MV.new nnzUB
      entries__ <- MV.new nnzUB
      let loop ix = \case
            Nothing -> return ix
            Just it@(j, wj, _, _, _, _, (mnr, x)) ->  do
                when wj $ MV.write majors__ j ix
                MV.write minors__ ix mnr
                MV.write entries__ ix x
                loop (ix + 1) (next it)
      nnz_ <- loop 0 $ next (-1, True, -1, -1, -1, -1, (0, 0))
      MV.write majors__ nmajor_ nnz_
      majors_ <- V.unsafeFreeze majors__
      minors_ <- V.unsafeFreeze $ MV.slice 0 nnz_ minors__
      entries_ <- V.unsafeFreeze $ MV.slice 0 nnz_ entries__
      return $ Matrix
          { nmajor = nmajor_
          , nminor = nminor_
          , majors = majors_
          , minors = minors_
          , entries = entries_
          }
