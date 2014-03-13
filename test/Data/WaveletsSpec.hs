module Data.WaveletsSpec (main, spec) where


import Data.Wavelets
import System.IO
import Test.Hspec



{-| The test waveletData below was transformed into several of the result dataFiles |-}

-- 3 sinusoids added together to make an interesting data set that is easy to understand 
testWaveletData :: [Double]
testWaveletData = [ (sin (pi*t*2))+ (sin (pi * t * 15) + (sin (pi * t * 0.002)))| t <- [1 .. 1000] ]


waveletHaar_packer_separate_testStub  :: IO [[Double]]
waveletHaar_packer_separate_testStub = do
  (read `fmap` readFile "./test/Data/haar_separate.tst" )


testWaveletHaar_PackerSeparate = dwt 1 haar wp_separate testWaveletData

compareWaveletHaarResults = do
  let rslt = testWaveletHaar_PackerSeparate
  ctrl <- waveletHaar_packer_separate_testStub
  return $ (length rslt ) == (length ctrl)


main :: IO ()
main = do
  haar_separate_test_data <- waveletHaar_packer_separate_testStub 
  hspec $ spec




-- | Have to bring in test data from a file to test this  
spec :: Spec
spec  = do
  describe "someFunction" $ do
    it "should work fine" $ do
      let tstData = testWaveletData
      tst <- compareWaveletHaarResults
      tst `shouldBe` True

