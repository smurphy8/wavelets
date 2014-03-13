{-# LANGUAGE NoImplicitPrelude,OverloadedStrings,BangPatterns #-}
module Main where
import Data.Wavelets
import Prelude

import System.IO



main = do
  ctrl <- return $ runWavelet testWaveletData
  seq (length.head $ ctrl) (print "done")





{-| The test waveletData below was transformed into several of the result dataFiles |-}

-- 3 sinusoids added together to make an interesting data set that is easy to understand 
testWaveletData :: [Double]
testWaveletData = [ (sin (pi*t*2))+ (sin (pi * t * 15) + (sin (pi * t * 0.002)))| t <- [1 .. 100000] ]


runWavelet !d = dwt 1 haar wp_separate d



