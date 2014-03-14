module Data.Wavelets (
    -- * Data types
    WaveletFilter(MakeWaveletFilter), WaveletPacker(MakeWaveletPacker),
    -- * Wavelet types
    haar, daubechies, daubechies_la, daubechies_bl, coiflet,
    -- * Creating custom wavelets
    wavelet, wavelet_s, verifywavelet,
    -- * Wavelet coefficient packers
    wp_interleave, wp_separate,
    -- * Transform functions
    dwt, idwt, idwtsynth
    )


where


-- | Wavelet filter type containing the wavelet linear filter, its
-- corresponding smoothing filter and their common length.
data WaveletFilter t = MakeWaveletFilter [t] [t] Int


instance (Show t) => Show (WaveletFilter t)
  where show (MakeWaveletFilter w s l) = show (w, s, l)


-- | Create wavelet object from wavelet filter coefficients.
wavelet :: Num t => [t] -> WaveletFilter t
wavelet l
  | odd len = error "wavelet: Length of wavelet must be even"
  | otherwise = MakeWaveletFilter l (scalingfromwavelet l) len
  where len = length l


-- | Create wavelet object from smoothing filter coefficients.
wavelet_s :: Num t => [t] -> WaveletFilter t
wavelet_s l
  | odd len = error "wavelet_s: Length of smoothing filter must be even"
  | otherwise = MakeWaveletFilter (waveletfromscaling l) l len
  where len = length l


-- | Compute coefficients of scaling (smoothing) filter from wavelet (detail)
-- filter.
scalingfromwavelet :: Num t => [t] -> [t]
scalingfromwavelet l = w2s l [] 1
  where
    w2s [] r sig = r
    w2s (h:t) r sig = w2s t (sig*h : r) (-sig)


-- | Compute coefficients of wavelet (detail) filter from scaling (smoothing)
-- filter.
waveletfromscaling :: Num t => [t] -> [t]
waveletfromscaling l = s2w l [] (-1)
  where
    s2w [] r sig = r
    s2w (h:t) r sig = s2w t (sig*h : r) (-sig)


-- | Verify the orthonormality relation of a wavelet filter.  Actually the
-- judgement has still to be made by the calling context, as the relation will
-- not be fulfilled exactly for reasons of numerical errors.  This function
-- merely computes the RMS deviation of the normality condition (first return
-- value) and orthogonality condition (second).

verifywavelet :: (Floating t, Eq t) => WaveletFilter t -> (t, t)
verifywavelet wavelet
  | len /= wlen || len /= slen = error "verifywavelet: Explicit length inconsistent with length of filters"
  | len == 0 = error "verifywavelet: Wavelet length must be non-zero"
  | odd len = error "verifywavelet: Length of wavelet must be even"
  | otherwise = (diagdev, offdev)
  where
    MakeWaveletFilter wfilt sfilt len = wavelet
    wlen = length wfilt
    slen = length sfilt
    allshifts filt = take (len `div` 2) $ iterate rot2left filt
    matrix = allshifts wfilt ++ allshifts sfilt
    diagdev = rms $ map (1-) $ zipWith scalarprod matrix matrix
    offdev = rms $ concat $ map (maphead scalarprod) $ take (len - 2) $ iterate tail matrix


-- | Map an operation the first argument of which is the head of the list onto
-- the tail of the list.
maphead :: (t -> t -> r) -> [t] -> [r]
maphead _ [] = []
maphead f (h:t) = map (f h) t


-- | Rotate a list two places to the left by removing and appending the first
-- two elements.
rot2left :: [t] -> [t]
rot2left (h:t) = t1 ++ [h, h1]
  where (h1:t1) = t


-- | Root of mean square of a list of numbers.
rms :: (Floating t, Eq t) => [t] -> t
rms l
  | s == 0 = 0
  | otherwise = sqrt $ s / fromIntegral len
  where
    len = length l
    s = sum $ map (\x -> x * x) l


-- | Haar wavelet, the simplest wavelet, which has length 2.
haar :: Floating t => WaveletFilter t
haar = wavelet [sqrt_2, -sqrt_2]


-- | Square root of 2, used in Haar wavelet.
sqrt_2 :: Floating t => t
sqrt_2 = sqrt 0.5


-- | Daubechies wavelets in even sizes from 4 to 20.
daubechies :: Floating t => Int -> WaveletFilter t
daubechies n
  | n < 4 || n > 20 = error "Daubechies wavelets are only available for lengths 4 to 20"
  | odd n = error "Wavelet lengths must be even (in daubechies function)"
  | otherwise = wavelet_s (daubechies_scaling !! (n `div` 2 - 2))

-- | Table of smoothing filters of Daubechies wavelets.
daubechies_scaling :: Floating t => [[t]]
daubechies_scaling = [
  [ 0.4829629131445341, 0.8365163037378079, 0.2241438680420134,
    -0.1294095225512604 ],
  [ 0.3326705529500826, 0.8068915093110926, 0.4598775021184916,
    -0.1350110200102546, -0.0854412738820267, 0.0352262918857095 ],
  [ 0.2303778133074431, 0.7148465705484058, 0.6308807679358788,
    -0.0279837694166834, -0.1870348117179132, 0.0308413818353661,
    0.0328830116666778, -0.0105974017850021 ],
  [ 0.1601023979741930, 0.6038292697971898, 0.7243085284377729,
    0.1384281459013204, -0.2422948870663824, -0.0322448695846381,
    0.0775714938400459, -0.0062414902127983, -0.0125807519990820,
    0.0033357252854738 ],
  [ 0.1115407433501094, 0.4946238903984530, 0.7511339080210954,
    0.3152503517091980, -0.2262646939654399, -0.1297668675672624,
    0.0975016055873224, 0.0275228655303053, -0.0315820393174862,
    0.0005538422011614, 0.0047772575109455, -0.0010773010853085 ],
  [ 0.0778520540850081, 0.3965393194819136, 0.7291320908462368,
    0.4697822874052154, -0.1439060039285293, -0.2240361849938538,
    0.0713092192668312, 0.0806126091510820, -0.0380299369350125,
    -0.0165745416306664, 0.0125509985560993, 0.0004295779729214,
    -0.0018016407040474, 0.0003537137999745 ],
  [ 0.0544158422431049, 0.3128715909143031, 0.6756307362972904,
    0.5853546836541907, -0.0158291052563816, -0.2840155429615702,
    0.0004724845739124, 0.1287474266204837, -0.0173693010018083,
    -0.0440882539307952, 0.0139810279173995, 0.0087460940474061,
    -0.0048703529934518, -0.0003917403733770, 0.0006754494064506,
    -0.0001174767841248 ],
  [ 0.0380779473638791, 0.2438346746125939, 0.6048231236901156,
    0.6572880780512955, 0.1331973858249927, -0.2932737832791761,
    -0.0968407832229524, 0.1485407493381306, 0.0307256814793395,
    -0.0676328290613302, 0.0002509471148340, 0.0223616621236805,
    -0.0047232047577520, -0.0042815036824636, 0.0018476468830564,
    0.0002303857635232, -0.0002519631889427, 0.0000393473203163 ],
  [ 0.0266700579005546, 0.1881768000776863, 0.5272011889317202,
    0.6884590394536250, 0.2811723436606485, -0.2498464243272283,
    -0.1959462743773399, 0.1273693403357890, 0.0930573646035802,
    -0.0713941471663697, -0.0294575368218480, 0.0332126740593703,
    0.0036065535669880, -0.0107331754833036, 0.0013953517470692,
    0.0019924052951930, -0.0006858566949566, -0.0001164668551285,
    0.0000935886703202, -0.0000132642028945  ]
  ]


-- | Least-asymmetric Daubechies wavelets in even sizes from 4 to 20.  These
-- transforms are a good compromise between good phase behaviour, separation of
-- frequency bands and lack of artifacts.  The order of the coefficients has
-- been reversed as recommended by Percival and Walden.
daubechies_la :: Floating t => Int -> WaveletFilter t
daubechies_la n
  | n < 4 || n > 20 = error "L. a. Daubechies wavelets are only available for lengths 4 to 20"
  | odd n = error "Wavelet lengths must be even (in daubechies_la function)"
  | otherwise = wavelet_s (daubechies_la_scaling !! (n `div` 2 - 2))

-- | Table of smoothing filters of least-asymmetric Daubechies wavelets.
daubechies_la_scaling :: Floating t => [[t]]
daubechies_la_scaling = [
  [ 0.4829629131445341, 0.8365163037378079, 0.2241438680420134,
    -0.1294095225512604 ],
  [ 0.3326705529500826, 0.8068915093110926, 0.4598775021184916,
    -0.1350110200102546, -0.0854412738820267, 0.0352262918857095 ],
  [ -0.0757657147893407, -0.0296355276459541, 0.4976186676324578,
    0.8037387518052163, 0.2978577956055422, -0.0992195435769354,
    -0.0126039672622612, 0.0322231006040713 ],
  [ 0.0195388827353869, -0.0211018340249298, -0.1753280899081075,
    0.0166021057644243, 0.6339789634569490, 0.7234076904038076,
    0.1993975339769955, -0.0391342493025834, 0.0295194909260734,
    0.0273330683451645 ],
  [ 0.0154041093273377, 0.0034907120843304, -0.1179901111484105,
    -0.0483117425859981, 0.4910559419276396, 0.7876411410287941,
    0.3379294217282401, -0.0726375227866000, -0.0210602925126954,
    0.0447249017707482, 0.0017677118643983, -0.0078007083247650 ],
  [ 0.0102681767084968, 0.0040102448717033, -0.1078082377036168,
    -0.1400472404427030, 0.2886296317509833, 0.7677643170045710,
    0.5361019170907720, 0.0174412550871099, -0.0495528349370410,
    0.0678926935015971, 0.0305155131659062, -0.0126363034031526,
    -0.0010473848889657, 0.0026818145681164 ],
  [ -0.0033824159513594, -0.0005421323316355, 0.0316950878103452,
    0.0076074873252848, -0.1432942383510542, -0.0612733590679088,
    0.4813596512592012, 0.7771857516997478, 0.3644418948359564,
    -0.0519458381078751, -0.0272190299168137, 0.0491371796734768,
    0.0038087520140601, -0.0149522583367926, -0.0003029205145516,
    0.0018899503329007 ],
  [ 0.0010694900326538, -0.0004731544985879, -0.0102640640276849,
    0.0088592674935117, 0.0620777893027638, -0.0182337707798257,
    -0.1915508312964873, 0.0352724880359345, 0.6173384491413523,
    0.7178970827642257, 0.2387609146074182, -0.0545689584305765,
    0.0005834627463312, 0.0302248788579895, -0.0115282102079848,
    -0.0132719677815332, 0.0006197808890549, 0.0014009155255716 ],
  [ 0.0007701598091030, 0.0000956326707837, -0.0086412992759401,
    -0.0014653825833465, 0.0459272392237649, 0.0116098939129724,
    -0.1594942788575307, -0.0708805358108615, 0.4716906668426588,
    0.7695100370143388, 0.3838267612253823, -0.0355367403054689,
    -0.0319900568281631, 0.0499949720791560, 0.0057649120455518,
    -0.0203549398039460, -0.0008043589345370, 0.0045931735836703,
    0.0000570360843390, -0.0004593294205481 ]
  ]


-- | Best localised Daubechies wavelets in even sizes from 4 to 20.  These
-- transforms are optimised to be nearly zero-phase.  The order of the
-- coefficients has been reversed as recommended by Percival and Walden.
daubechies_bl :: Floating t => Int -> WaveletFilter t
daubechies_bl n
  | n < 4 || n > 20 = error "B. l. Daubechies wavelets are only available for lengths 4 to 20"
  | odd n = error "Wavelet lengths must be even (in daubechies_bl function)"
  | n == 14 = wavelet_s daubechies_bl_scaling_14
  | n == 18 = wavelet_s daubechies_bl_scaling_18
  | n == 20 = wavelet_s daubechies_bl_scaling_20
  | otherwise = wavelet_s (daubechies_la_scaling !! (n `div` 2 - 2))

-- | Table of smoothing filters of best localised Daubechies wavelets.
daubechies_bl_scaling_14 :: Floating t => [t]
daubechies_bl_scaling_14 = [
    0.0120154192834842, 0.0172133762994439, -0.0649080035533744,
    -0.0641312898189170, 0.3602184608985549, 0.7819215932965554,
    0.4836109156937821, -0.0568044768822707, -0.1010109208664125,
    0.0447423494687405, 0.0204642075778225, -0.0181266051311065,
    -0.0032832978473081, 0.0022918339541009 ]
daubechies_bl_scaling_18 :: Floating t => [t]
daubechies_bl_scaling_18 = [
    0.0002594576266544, -0.0006273974067728, -0.0019161070047557,
    0.0059845525181721, 0.0040676562965785, -0.0295361433733604,
    -0.0002189514157348, 0.0856124017265279, -0.0211480310688774,
    -0.1432929759396520, 0.2337782900224977, 0.7374707619933686,
    0.5926551374433956, 0.0805670008868546, -0.1143343069619310,
    -0.0348460237698368, 0.0139636362487191, 0.0057746045512475 ]
daubechies_bl_scaling_20 :: Floating t => [t]
daubechies_bl_scaling_20 = [
    0.0008625782242896, 0.0007154205305517, -0.0070567640909701,
    0.0005956827305406, 0.0496861265075979, 0.0262403647054251,
    -0.1215521061578162, -0.0150192395413644, 0.5137098728334054,
    0.7669548365010849, 0.3402160135110789, -0.0878787107378667,
    -0.0670899071680668, 0.0338423550064691, -0.0008687519578684,
    -0.0230054612862905, -0.0011404297773324, 0.0050716491945793,
    0.0003401492622332, -0.0004101159165852 ]


-- | Coifman wavelets in mutiple-of-6 sizes from 6 to 30.  The order of the
-- coefficients has been reversed as recommended by Percival and Walden.
coiflet :: Floating t => Int -> WaveletFilter t
coiflet n
  | n < 6 || n > 30 || (n `div` 6) * 6 /= n = error "Coifman wavelets are only available for lengths that are multiples of 6 from 6 to 30"
  | otherwise = wavelet_s (coiflet_scaling !! (n `div` 6 - 1))

-- | Table of smoothing filters of Coifman wavelets.
coiflet_scaling :: Floating t => [[t]]
coiflet_scaling = [
  [ -0.0156557285289848, -0.0727326213410511, 0.3848648565381134,
    0.8525720416423900, 0.3378976709511590, -0.0727322757411889 ],
  [ -0.0007205494453679, -0.0018232088707116, 0.0056114348194211,
    0.0236801719464464, -0.0594344186467388, -0.0764885990786692,
    0.4170051844236707, 0.8127236354493977, 0.3861100668229939,
    -0.0673725547222826, -0.0414649367819558, 0.0163873364635998 ],
  [ -0.0000345997728362, -0.0000709833031381, 0.0004662169601129,
    0.0011175187708906, -0.0025745176887502, -0.0090079761366615,
    0.0158805448636158, 0.0345550275730615, -0.0823019271068856,
    -0.0717998216193117, 0.4284834763776168, 0.7937772226256169,
    0.4051769024096150, -0.0611233900026726, -0.0657719112818552,
    0.0234526961418362, 0.0077825964273254, -0.0037935128644910 ],
  [ -0.0000017849850031, -0.0000032596802369, 0.0000312298758654,
    0.0000623390344610, -0.0002599745524878, -0.0005890207562444,
    0.0012665619292991, 0.0037514361572790, -0.0056582866866115,
    -0.0152117315279485, 0.0250822618448678, 0.0393344271233433,
    -0.0962204420340021, -0.0666274742634348, 0.4343860564915321,
    0.7822389309206135, 0.4153084070304910, -0.0560773133167630,
    -0.0812666996808907, 0.0266823001560570, 0.0160689439647787,
    -0.0073461663276432, -0.0016294920126020, 0.0008923136685824 ],
  [ -0.0000000951765727, -0.0000001674428858, 0.0000020637618516,
    0.0000037346551755, -0.0000213150268122, -0.0000413404322768,
    0.0001405411497166, 0.0003022595818445, -0.0006381313431115,
    -0.0016628637021860, 0.0024333732129107, 0.0067641854487565,
    -0.0091642311634348, -0.0197617789446276, 0.0326835742705106,
    0.0412892087544753, -0.1055742087143175, -0.0620359639693546,
    0.4379916262173834, 0.7742896037334738, 0.4215662067346898,
    -0.0520431631816557, -0.0919200105692549, 0.0281680289738655,
    0.0234081567882734, -0.0101311175209033, -0.0041593587818186,
    0.0021782363583355, 0.0003585896879330, -0.0002120808398259 ]
  ]



-- | Delay a series of numbers by prepending zeros as needed.  Negative delays
-- are realised by dropping leading elements.
delay :: (Num t) => Int -> [t] -> [t]
delay d
  | d <= 0 = drop (-d)
  | otherwise = (++) (replicate d 0)


-- | A data type for storing wavelet coefficients after performing the
-- transform.  The first two parameters are a start value and a function for
-- adding the list of wavelet coefficients from a stage of the transform to the
-- stored set.  This function will be applied to the smoothing coefficient list
-- and then to the detail coefficients with rising frequency.  The other two
-- type parameters are functions for unpacking the stored coefficients and
-- recreating the lists of coefficients for each stage.  The third should
-- return the detail coefficients from highest to lowest frequencies with each
-- successive call.  The fourth should return the smoothing coefficients.
data WaveletPacker t c = MakeWaveletPacker c ([t] -> c -> c) (c -> ([t], c)) (c -> [t])


-- | This packer interleaves wavelet coefficients in one list so that the n-th
-- level coefficient occurs every 2^n entries starting with the 2^(n-1)th
-- (1-based) (n in [1 .. order of the transform] starting with the highest
-- frequency; the smoothing coefficient takes the remaining place at the end of
-- each block of 2^order coefficients).
wp_interleave :: Num t => WaveletPacker t [t]
wp_interleave = MakeWaveletPacker [] interleave_pack interleave_unpack id
interleave_pack :: [t] -> [t] -> [t]
interleave_pack l1 [] = l1
interleave_pack [] l2 = l2
interleave_pack (h1:t1) (h2:t2) = h1 : h2 : interleave_pack t1 t2
interleave_unpack [] = ([], [])
interleave_unpack (h:t) = (h : fst rest, h2 : snd rest)
  where
    (h2:t2) = t
    rest = interleave_unpack t2


-- | This packer stores the wavelet coefficients corresponding to different
-- frequencies / levels of detail in separate lists.  The highest-frequency
-- coefficients are in the first list.  Note that the number of coefficients
-- describing a given stretch of the original time series differs by a factor
-- of two between neighbouring coefficient lists.
wp_separate :: Num t => WaveletPacker t [[t]]
wp_separate = MakeWaveletPacker [] (:) separate_unpack head
separate_unpack :: [[t]] -> ([t], [[t]])
separate_unpack (h:t) = (h, t)


-- | Scalar product of two vectors; or sum of products of two lists of numbers.
scalarprod :: Num t => [t] -> [t] -> t
scalarprod l1 l2 = sum $ zipWith (*) l1 l2


-- | The number of zeros to prepend to a time series so that the filter just
-- overlaps with the data and its middle is at an odd position relative to the
-- start of the data.
databias :: Num t => WaveletFilter t -> Int
databias wavelet = len `div` 2 - 1 + 2 * (len `div` 4)
   where MakeWaveletFilter _ _ len = wavelet


-- | One stage of the discrete wavelet transform.  The first argument is the
-- type of wavelet.  The second argument is the original data series at the
-- first stage, or successively higher-scale smoothing coefficients.  The
-- returned pair contains the resulting detail and smoothing coefficients.
dwtstage :: Num t => WaveletFilter t -> [t] -> ([t], [t])
dwtstage wavelet lowscalcoeffs = (waveletcoeffs, scalingcoeffs)
  where
    MakeWaveletFilter wfilt sfilt len = wavelet
    timeseries = takeWhile (not . null) $ map (flip delay lowscalcoeffs) [0,-2..]
    waveletcoeffs = map (scalarprod wfilt) timeseries
    scalingcoeffs = map (scalarprod sfilt) timeseries


-- | Discrete wavelet transform.  The first argument is the number of stages of
-- the transform.  The result will be (stages + 1) sets of coefficients,
-- including one set of smoothed values.  The second argument is the type of
-- wavelet, the third the way to wrap the resulting coefficients.  The last
-- argument is the series of data.
--
-- The data are interpreted as part of an infinitely long series, not (as is
-- sometimes done) a set of blocks with implied cyclical boundary conditions.
-- To allow reversing the transform, enough zeros are prepended to the data so
-- that the first coefficients from the last stage corresponds to minimal
-- overlap between the filter and its input data (i.e. one or two values).
dwt :: Num t => Int -> WaveletFilter t -> WaveletPacker t c -> [t] -> c
dwt order wavelet wrapper indata
  | order <= 0 = error "dwt: order has to be at least 1"
  | otherwise = foldr wrap_step wrap_init $ getcoeffs bands
  where
    stage = dwtstage wavelet
    MakeWaveletPacker wrap_init wrap_step _ _ = wrapper
    init = ([], delay totalbias indata)
    totalbias = (2 ^ order - 1) * (databias wavelet)
    bands = take order $ tail $ iterate (dwtstage wavelet . snd) init
    getcoeffs :: [(a, a)] -> [a]
    getcoeffs l
      | length l == 1 = [fst h, snd h]
      | otherwise = fst h : getcoeffs t
      where h:t = l


-- | Scalar multiplication: multiply each element of a list with a constant.
scalarmult :: Num t => t -> [t] -> [t]
scalarmult s = map (s *)

 
-- | Similar to zipWith (+), but the result is always as long as the second
-- argument.  "Missing" elements of the first list are taken as zero.
asymzipadd :: Num t => [t] -> [t] -> [t]
asymzipadd _ [] = []
asymzipadd [] l = l
asymzipadd (h1:t1) (h2:t2) = (h1 + h2) : asymzipadd t1 t2


-- | One stage of the inverse discrete wavelet transform.  The first argument
-- is the type of wavelet, the other two arguments are the detail and smoothing
-- coefficients of the current stage.  The return value is the smoothing
-- coefficients one level more detailed, or (for the last stage) the original
-- data.
idwtstage :: Num t => WaveletFilter t -> [t] -> [t] -> [t]
idwtstage wavelet waveletcoeffs scalingcoeffs = zipWith (+) wavepart scalepart
  where
    MakeWaveletFilter wfilt sfilt _ = wavelet
    wavepart = superimpose $ map (flip scalarmult wfilt) waveletcoeffs
    scalepart = superimpose $ map (flip scalarmult sfilt) scalingcoeffs
    superimpose :: Num t1 => [[t1]] -> [t1]
    superimpose [] = []
    superimpose (h:t) = h1 : h2 : asymzipadd t2 (superimpose t)
      where h1:t1 = h
            h2:t2 = t1


-- | Inverse discrete wavelet transform suitable for synthesis.  The first
-- three arguments are the same as for dwt: The number of stages of the
-- transform, the wavelet type and the coefficient packer.  The last argument
-- is the result of the transform with the same packer.  The result contains
-- the zero values that dwt has prepended to create a reversible result (which
-- may in fact  differ slightly from zero for numerical reasons).  This allows
-- to use this function for synthesis (where the first values will usually be
-- non-zero).  If you do not want them, use idwt instead.
idwtsynth :: Num t => Int -> WaveletFilter t -> WaveletPacker t c -> c -> [t]
idwtsynth order wavelet wrapper coeffs
  | order <= 0 = error "idwtsynth: order has to be at least 1"
  | otherwise = foldr (idwtstage wavelet) lowestscale unpacked
  where
    MakeWaveletPacker _ _ unwrap_step unwrap_last = wrapper
    unpacked_pre = take order $ tail $ iterate (unwrap_step . snd) ([0], coeffs)
    -- unpacked list is in order from high to low frequencies
    lowestscale = unwrap_last $ snd $ last unpacked_pre
    unpacked = map fst unpacked_pre


-- | Inverse discrete wavelet transform reversing the transform of the dwt
-- function.  Arguments and operation are the same as for idwtsynth, only the
-- zero values prepended to the data by dwt are removed.  One approximate zero
-- value will remain after the original data for wavelets of lengths 4*n+2,
-- because the result of the inverse transform (before removing the prefix) is
-- always even in length.
idwt :: Num t => Int -> WaveletFilter t -> WaveletPacker t c -> c -> [t]
idwt order wavelet wrapper coeffs = delay (-totalbias) synthresult
  where
    totalbias = (2 ^ order - 1) * (databias wavelet)
    synthresult = idwtsynth order wavelet wrapper coeffs
