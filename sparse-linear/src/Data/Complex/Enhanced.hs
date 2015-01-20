{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeFamilies #-}

module Data.Complex.Enhanced
       ( RealOf, ComplexOf, IsReal(..), IsImag(..)
       , Complex(..)
       ) where

import Data.Complex
import Foreign.Storable.Complex ()

type family RealOf a where
  RealOf (Complex a) = a
  RealOf a = a

type family ComplexOf a where
  ComplexOf (Complex a) = Complex a
  ComplexOf a = a

class IsReal a where
  real :: Num (RealOf a) => RealOf a -> a
  conj :: a -> a
  mag :: a -> RealOf a

class IsImag a where
  imag :: Num (RealOf a) => RealOf a -> a

instance IsReal Double where
  {-# INLINE real #-}
  {-# INLINE conj #-}
  {-# INLINE mag #-}
  real = id
  conj = id
  mag = abs

instance IsReal (Complex Double) where
  {-# INLINE real #-}
  {-# INLINE conj #-}
  {-# INLINE mag #-}
  real = (:+ 0)
  conj = conjugate
  mag = magnitude

instance IsImag (Complex Double) where
  {-# INLINE imag #-}
  imag = (0 :+)