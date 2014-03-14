module Data.WaveletsSpec (main, spec) where

import Test.Hspec

import Data.Wavelets
import System.IO



{-| The test waveletData below was transformed into several of the result dataFiles |-}

-- 3 sinusoids added together to make an interesting data set that is easy to understand 
testWaveletData :: [Double]
testWaveletData = [ ((sin (pi*t*2))+ (sin (pi * t * 2.1) + (sin (pi * t * 2.002))))* 12 | t <- [1 .. 10000] ]


waveletHaar_packer_separate_testStub  :: IO [[Double]]
waveletHaar_packer_separate_testStub = do
  (read `fmap` readFile "./test/Data/haar_separate.tst" )


testWaveletHaar_PackerSeparate = dwt 10 haar wp_separate testWaveletData

compareWaveletHaarResults = do
  let rslt = testWaveletHaar_PackerSeparate
  ctrl <- waveletHaar_packer_separate_testStub
  return $ (length rslt ) == (length ctrl)



testIdwtSyth = idwtsynth 3 haar wp_separate $  drop 4 testWaveletHaar_PackerSeparate


main :: IO ()
main = do
  haar_separate_test_data <- waveletHaar_packer_separate_testStub 
  hspec $ spec

-- | Have to bring in test data from a file to test this  
spec :: Spec
spec  = do
  describe "dwt 1 haar wp_separate testWaveletData" $ do
    it "should return a dwt of the same length as the test file" $ do
      let tstData = testWaveletData
      tst <- compareWaveletHaarResults
      tst `shouldBe` True

