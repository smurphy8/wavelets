module Data.Wavelets.ReconstructionSpec (main, spec) where

import Test.Hspec

import Data.Wavelets.Reconstruction
import Data.Wavelets
import System.IO
import qualified Data.Vector.Storable as V

{-| The test waveletData below was transformed into several of the result dataFiles |-}

-- 3 sinusoids added together to make an interesting data set that is easy to understand 
testWaveletData :: [Double]
testWaveletData = [ ((sin (pi*t*2))+ (sin (pi * t * 2.1) + (sin (pi * t * 2.002))))* 12 | t <- [1 .. 1000] ]

impulse = replicate 499 1
          ++ replicate 1 100 
          ++ replicate 500 1

impulseAvg = (1.3956047904191924)

waveletHaar_packer_separate_testStub  :: IO [[Double]]
waveletHaar_packer_separate_testStub = do
  (read `fmap` readFile "./test/Data/haar_separate.tst" )


testWaveletHaar_PackerSeparate = dwt 12 haar wp_separate impulse

compareWaveletHaarResults = do
  let rslt = testWaveletHaar_PackerSeparate
  ctrl <- waveletHaar_packer_separate_testStub
  return $ (length rslt ) == (length ctrl)


testReconstructTimeSeries  = reconstructTimeSeries (12-n) haar wp_separate $ drop n testWaveletHaar_PackerSeparate
  where n = 3

main :: IO ()
main = do
  haar_separate_test_data <- waveletHaar_packer_separate_testStub 
  hspec $ spec

-- | Have to bring in test data from a file to test this  
spec :: Spec
spec  = do
  describe "reconstructTimeSeries" $ do
    it "shouldReturn a scaled wavelet" $ do 
      let tstData = testWaveletData
      tst <- compareWaveletHaarResults
      tst `shouldBe` True


