{-# LANGUAGE DataKinds #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE TypeFamilies #-}

module Numeric.LinearAlgebra.Sparse where

import Data.Complex
import Data.Ord (comparing)
import qualified Data.Vector.Algorithms.Intro as Intro
import Data.Vector.Unboxed (Unbox, Vector)
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as MV
import GHC.Stack (errorWithStackTrace)

data Orient = Row | Col

class Linear a where
    data family Matrix (or :: Orient) a
    compress :: OrientOf or => Int -> Int -> Vector (Int, Int, a) -> Matrix or a

sortCoords :: Unbox a => Orient -> Vector (Int, Int, a) -> Vector (Int, Int, a)
sortCoords = \case
    Row -> V.modify $ Intro.sortBy $ comparing $ \(r, c, _) -> (r, c)
    Col -> V.modify $ Intro.sortBy $ comparing $ \(r, c, _) -> (c, r)
{-# INLINE sortCoords #-}

countMajors :: Int -> Vector Int -> Vector Int
countMajors dim ixs = V.create $ do
    v <- MV.replicate 0 (dim + 1)
    V.forM_ ixs $ \ix ->
        if ix >= dim
          then errorWithStackTrace "index out of bounds"
          else MV.write v (ix + 1) . (+ 1) =<< MV.read v (ix + 1)
    return v
{-# INLINE countMajors #-}

instance Linear Double where
    data Matrix or Double
        = MatrixD {-# UNPACK #-} !Int !(Vector Int) !(Vector Int) !(Vector Double)

    compress = \r c coords ->
        let orient = orientOf compressed
            (rows, cols, entries) = V.unzip3 $ sortCoords orient coords
            minors = case orient of { Row -> cols; Col -> rows }
            majors = case orient of
              Row -> countMajors r rows
              Col -> countMajors c cols
            minorDim = case orient of { Row -> c; Col -> r }
            compressed = MatrixD minorDim majors minors entries
        in compressed
    {-# INLINE compress #-}

instance (RealFloat a, Unbox a) => Linear (Complex a) where
    data Matrix or (Complex a)
        = MatrixZ {-# UNPACK #-} !Int !(Vector Int) !(Vector Int) !(Vector a) !(Vector a)

    compress = \r c coords ->
        let orient = orientOf compressed
            (rows, cols, entries) = V.unzip3 $ sortCoords orient coords
            minors = case orient of { Row -> cols; Col -> rows }
            majors = case orient of
              Row -> countMajors r rows
              Col -> countMajors c cols
            minorDim = case orient of { Row -> c; Col -> r }
            reals = V.map realPart entries
            imags = V.map imagPart entries
            compressed = MatrixZ minorDim majors minors reals imags
        in compressed
    {-# INLINE compress #-}

class OrientOf or where
    orientOf :: Matrix or a -> Orient

instance OrientOf Row where
    orientOf = \_ -> Row
    {-# INLINE orientOf #-}

instance OrientOf Col where
    orientOf = \_ -> Col
    {-# INLINE orientOf #-}
