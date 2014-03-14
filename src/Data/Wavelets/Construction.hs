{-# LANGUAGE OverloadedStrings #-}

module Data.Wavelets.Construction where 

import Data.Wavelets 
import Numeric.Statistics
import Data.Vector hiding (map)
-- import Linear


-- | Vector wrapper for dwt so that when I am ready to turn the whole thing vectorized I can without too much hassle

vdwt :: (Num t) =>  Int -> WaveletFilter t -> WaveletPacker t c -> Vector t -> c
vdwt i wft wptc vt = dwt i wft wptc (toList vt)



