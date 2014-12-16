{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse where

import Data.Foldable (forM_)
import Data.Ord (comparing)
import qualified Data.Vector.Algorithms.Intro as Intro
import Data.Vector.Unboxed (Unbox, Vector)
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as MV
import GHC.Stack (errorWithStackTrace)

data Orient = Row | Col

data Matrix (or :: Orient) a
    = Matrix {-# UNPACK #-} !Int !(Vector Int) !(Vector (Int, a))

class OrientOf or where
    orientOf :: Matrix or a -> Orient

instance OrientOf Row where
    orientOf = \_ -> Row
    {-# INLINE orientOf #-}

instance OrientOf Col where
    orientOf = \_ -> Col
    {-# INLINE orientOf #-}

sortCoords :: Unbox a => Orient -> Vector (Int, Int, a) -> Vector (Int, Int, a)
sortCoords = \case
    Row -> V.modify $ Intro.sortBy $ comparing $ \(r, c, _) -> (r, c)
    Col -> V.modify $ Intro.sortBy $ comparing $ \(r, c, _) -> (c, r)
{-# INLINE sortCoords #-}

countMajors :: Int -> Vector Int -> Vector Int
countMajors dim_ ixs = V.create $ do
    v <- MV.replicate 0 (dim_ + 1)
    V.forM_ ixs $ \ix ->
        if ix >= dim_
          then errorWithStackTrace "index out of bounds"
          else MV.write v (ix + 1) . (+ 1) =<< MV.read v (ix + 1)
    return v
{-# INLINE countMajors #-}

compress
  :: (OrientOf or, Unbox a)
  => Int
  -> Int
  -> Vector (Int, Int, a)
  -> Matrix or a
compress = \r c coords ->
    let orient = orientOf compressed
        (rows, cols, entries) = V.unzip3 $ sortCoords orient coords
        minors = case orient of { Row -> cols; Col -> rows }
        majors = case orient of
          Row -> countMajors r rows
          Col -> countMajors c cols
        minorDim = case orient of { Row -> c; Col -> r }
        compressed = Matrix minorDim majors $ V.zip minors entries
    in compressed
{-# INLINE compress #-}

slices :: Unbox a => Matrix or a -> [(Int, Vector (Int, a))]
slices = \(Matrix _ rows entries) -> do
    i <- enumFromTo 0 $ V.length rows - 1
    let start = rows V.! i
        end = rows V.! (i + 1)
    return $ (i, V.slice start (end - start) entries)
{-# INLINE slices #-}

dim :: OrientOf or => Matrix or a -> Int
dim = \mat@(Matrix minorDim majors _) ->
    case orientOf mat of
      Row -> V.length majors - 1
      Col -> minorDim
{-# INLINE dim #-}

codim :: OrientOf or => Matrix or a -> Int
codim = \mat@(Matrix minorDim majors _) ->
    case orientOf mat of
      Col -> V.length majors - 1
      Row -> minorDim
{-# INLINE codim #-}

mulV
  :: (Num a, Unbox a)
  => Matrix Row a
  -> Vector a
  -> Vector a
mulV = \mat xs ->
    if V.length xs /= dim mat
      then errorWithStackTrace "matrix dimension does not match vector"
      else V.create $ do
          ys <- MV.replicate (V.length xs) 0
          forM_ (slices mat) $ \(r, as) -> do
              y <- MV.read ys r
              let ax = V.sum $ flip V.map as $ \(c, a) -> a * (xs V.! c)
              MV.write ys r $ ax + y
          return ys
{-# INLINE mulV #-}
