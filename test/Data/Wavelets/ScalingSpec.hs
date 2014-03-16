module Data.Wavelets.ScalingSpec (main, spec) where

import Test.Hspec
import Data.Wavelets.Reconstruction
import Data.Wavelets.Scaling
import Data.Wavelets
import System.IO
import Numeric.Statistics
import qualified Data.Vector.Storable as V

{-| The test waveletData below was transformed into several of the result dataFiles |-}

-- 3 sinusoids added together to make an interesting data set that is easy to understand 
testWaveletData :: [Double]
testWaveletData = [ ((sin (pi*t*2))+ (sin (pi * t * 2.1) + (sin (pi * t * 2.002))))* 12 | t <- [1 .. 1000] ]

impulse = replicate 499 1
          ++ replicate 1 200
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


main :: IO ()
main = do
  haar_separate_test_data <- waveletHaar_packer_separate_testStub 
  hspec $ spec

-- | Have to bring in test data from a file to test this  
spec :: Spec
spec  = do
  describe "absAvgItrFcn" $ do
    it "should return the average when folded over a series" $ do 
      let tstData = replicate 10 1.3
          (avg,tot) = foldl (\(a,t) b -> absAvgItrFcn a b (t+1) ) (0,0) tstData
      avg `shouldBe` 1.3
  describe "absMaxItrFcn" $ do
    it "should return the maximum when folded overa  series" $ do 
      let tstData = (replicate 10 1.3) ++ [100] :: [Double] 
          mx = foldl (\a b -> absMaxItrFcn a b )  0 tstData
      mx `shouldBe` 100
  describe "minMaxAndAverage" $ do
    it "should return the average, max , length of a given vector of doubles" $ do 
      let tstData = testVectorData
          (mn,mx,av,tot) = minMaxAndAverage tstData
      mx `shouldBe` 100



testVectorData = V.fromList impulse 



testReconstructTimeSeries  = reconstructTimeSeries (12-n) haar wp_separate $ drop n testWaveletHaar_PackerSeparate
  where n = 3



--OSF (SeriesFactors {seriesMin = 1.0, seriesMax = 100.0, seriesMean = 1.099, seriesCount = 1000})
testComputeOldSeriesmatrix = OSF . computeSeriesFactors $ testVectorData


-- NSF (SeriesFactors {seriesMin = 2.8284271247461934, seriesMax = 37.830212793480314, seriesMean = 3.108441410096067, seriesCount = 125})
testComputeNewSeriesmatrix = NSF . computeSeriesFactors $ testReconstructTimeSeries


-- SM (2><1)
--  [ -7.374250048383981
--  , 2.8382461515606137 ]
testComputeScalingMatrix = computeScalingMatrix testComputeNewSeriesmatrix testComputeOldSeriesmatrix


testApplyScalingMatrix = applyScalingMatrix testComputeScalingMatrix testReconstructTimeSeries

testGetNewScalingFactors = computeSeriesFactors testApplyScalingMatrix
