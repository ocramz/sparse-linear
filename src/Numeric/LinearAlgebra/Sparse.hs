{-# LANGUAGE DataKinds #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse where

import Control.Monad.ST (ST)
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
