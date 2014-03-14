{-# LANGUAGE OverloadedStrings #-}

module Data.Wavelets.Reconstruction where 
import Prelude hiding (map,maximum,minimum)
import Data.Wavelets 
import Numeric.Statistics
import Statistics.Sample
import Data.Vector.Storable -- hiding (map)
-- import Linear

-- |explicitly name a Fractional and Ord instance on the idwt function
fidwt ::  Int -> WaveletFilter Double -> WaveletPacker Double c -> c -> [Double]
fidwt = idwt


-- | adding a max value so that you can scale versus the other max
reconstructTimeSeries ::  Double -> Int -> WaveletFilter Double -> WaveletPacker Double c -> c -> Vector Double
reconstructTimeSeries tavg i wft wptc c = scaledReconstruction 
    where rslt = fidwt i wft wptc c
          vRslt = fromList rslt 
--          (avg,max) = map abs vRslt
--          avg     = mean absRslt          
          scaledReconstruction = vRslt



-- | compute the average of the absolute value of a foldable type by setting it up like 
-- foldl (\(a,t) b -> absAvgItrFcn a b (t+1) ) (0,0) foldableThing
absAvgItrFcn ::(Num a, Fractional a)=> a -> a -> a -> (a,a)
absAvgItrFcn a b t = avgtuple
   where 
     diff = (a - abs(b))
     div  = abs(diff/t)
     tot  = a + div 
     avgtuple = (tot,t)


-- | compute the maximum absolute value of a foldable type by setting it up like 
-- foldl (\mx b -> absMaxItrFcn mx b ) 0 foldableThing
absMaxItrFcn a b = max a (abs b)


-- | compute the max and average of a Vector in one pass
maxAndAverage :: Vector Double -> (Double, Double , Double )
maxAndAverage vl = foldl' maxAndAverageItrFcn seed vl 
    where
      seed = (0,0,0)
      maxAndAverageItrFcn (avg, mx, tot) i = maxAndAverageItrFcn' avg mx tot i
      maxAndAverageItrFcn' avg mx tot i = let
          (avg',tot') = absAvgItrFcn avg i (tot + 1)
          mx' = absMaxItrFcn mx i
          in (avg', mx', tot')
