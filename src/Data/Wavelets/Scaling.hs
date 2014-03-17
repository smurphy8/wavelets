{-# LANGUAGE OverloadedStrings, DeriveGeneric #-}

module Data.Wavelets.Scaling where 
import Prelude hiding (length,replicate,map,maximum,minimum)
-- import Data.Wavelets 
import Numeric.Statistics
import Statistics.Function
import Statistics.Sample
import Data.Vector.Storable -- hiding (map)
import Data.Packed.Matrix
-- import Linear

{-| Scale factors are designed to make it possible to use a scaled wavelet to reconstruct an approximation 
    of the original data series. ordinary least squares is used so that translation is available for approximation,
    in addition to scaling. 


|-}


data SeriesFactors = SeriesFactors { 
      seriesMin :: Double ,
      seriesMax :: Double ,
      seriesMean :: Double ,
      seriesCount  :: Int } deriving (Show,Generic)

computeSeriesFactors :: Vector Double -> SeriesFactors
computeSeriesFactors v = SeriesFactors mn mx avg tot
    where (mn, mx,avg, tot) = minMaxAndAverage v 
 

newtype OldSeriesFactors = OSF SeriesFactors  deriving (Show)    
newtype NewSeriesFactors = NSF SeriesFactors  deriving (Show)  
newtype NewSeriesMatrix = NSM (Matrix Double) deriving (Show)   
newtype OldSeriesMatrix = OSM (Matrix Double) deriving (Show)   

-- | x scaling matrix is B in Y = XB + e matrix formulation of ordinary least squares

newtype ScalingMatrix = SM (Matrix Double) deriving (Show)

type SeriesMatrix = Matrix Double 


-- ScaleMatrix is the scaling matrix for the least squares problem
computeSeriesMatrix :: SeriesFactors -> SeriesMatrix
computeSeriesMatrix (SeriesFactors mn mx amean _ ) = (fromColumns vlst)
    where onesList = replicate 3 1
          nsf' = fromList [mn,mx,amean]
          vlst = [onesList,nsf']

computeNewSeriesMatrix :: NewSeriesFactors -> NewSeriesMatrix
computeNewSeriesMatrix (NSF sf) = NSM $ computeSeriesMatrix sf


computeSeriesResult :: SeriesFactors -> SeriesMatrix
computeSeriesResult (SeriesFactors mn mx avg _ ) = (fromColumns vlst)
    where nsf' = fromList [mn,mx,avg]
          vlst = [nsf']

computeOldSeriesMatrix :: OldSeriesFactors -> OldSeriesMatrix
computeOldSeriesMatrix (OSF sf) = OSM $ computeSeriesResult sf


-- | ordinary least squares estimation for the multivariate model
--   Y = X B + e        rows are observations, columns are elements
--   mean e = 0, cov e = kronecker s I

computeScalingMatrix :: NewSeriesFactors -> OldSeriesFactors -> ScalingMatrix
computeScalingMatrix nsf osf = let (OSM _Y) = computeOldSeriesMatrix osf
                                   (NSM _X) = computeNewSeriesMatrix nsf 
                                   (_B,_ ,_ ) = ols _X _Y 
                               in SM _B
    


-- | compute the average of the absolute value of a foldable type by setting it up like 
-- foldl (\(a,t) b -> absAvgItrFcn a b (t+1) ) (0,0) foldableThing
absAvgItrFcn ::(Num a, Fractional a)=> a -> a -> a -> (a,a)
absAvgItrFcn a b t = avgtuple
   where 
     diff = (a - abs(b))
     dv  = abs(diff/t)
     tot  = a + dv 
     avgtuple = (tot,t)


-- | compute the maximum absolute value of a foldable type by setting it up like 
-- foldl (\mx b -> absMaxItrFcn mx b ) 0 foldableThing
absMaxItrFcn :: (Num a, Ord a) => a -> a -> a
absMaxItrFcn a b = max a (abs b)


-- | compute the max min and average of a Vector in one pass 
-- (avg, max, min, tot)
minMaxAndAverage :: Vector Double -> (Double, Double ,Double , Int )
minMaxAndAverage v = (mn,mx,av,l)
    where
      av = mean v
      (mn, mx) = minMax v
      l = length v



applyScalingMatrix :: ScalingMatrix -> Vector Double -> Vector Double 
applyScalingMatrix (SM sm) v = map scaleFcn v
    where tCoef = sm @@> (0,0)
          sCoef = sm @@> (1,0)
          scaleFcn x = sCoef * x + tCoef 


    
