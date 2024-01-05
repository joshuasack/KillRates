--- ***********************
--- **** InVivo module ****
--- ***********************
--- ---------------------------------------------------------
--- **** Brief description of the main procedure inVivo: ****
--- ---------------------------------------------------------
---
--- > run a simulation with cell number given by constant parameters for timeStepsInVivo time steps
--- 
--- effects (in vivo)
--- > write position data in the dataPosition\ folder for every stage: beginning of run (time step 0), and 
---                                                                    end of each time step from 1 on)
--- > write distance data in the dataPosition\ folder for every stage
--- > write cell counts in the dataCounts\ folder for each run: 
---     each row is a time step (starting with 1)
---     each column corresponds to a type of cell same ordering as in the composition data type:
---       TumLoMHC TumHiMHC FCa FCu FCr1 FCr2 PCa PCu PCr1 PCr2 FNa FNu FNr1 FNr2 PNa PNu PNr1 PNr2
--- > write average (among all runs) counts in dataCounts\ 
---     each row is a time step (starting with 1)
---     with 6 columns corresponding to:
---     TumLoMHC TumHiMHC CTL_Pre_Recognize CTL_Post_Recognize NK_Cell_Pre_Recognize KN-Cell_Post_Recognize  
--- > Write various statistics for each run in the dataStats\ folder 
---     with explanation printed to make it easier for us to read (but harder for MatLab to read)
--- > Write average statistics for each run in the dataStats\ folder
---     with explanation printed too 
--- > print various run-time progress
--- > print the total number of tumor cells lysed
--- > print the maximum number of effector cells that have been in the same grid element as a tumor cell after a movement phase.
---
--- ------------------------------
--- **** Position Data Format ****
--- ------------------------------
--- 
--- Each position data file is for a particular time step in a particular run
--- It consists of an array with:
--- 18 columns (corresponding to the components of the composition data type)
---     TumLoMHC TumHiMHC FCa FCu FCr1 FCr2 PCa PCu PCr1 PCr2 FNa FNu FNr1 FNr2 PNa PNu PNr1 PNr2
--- gridSize rows (corresponding to the positions: the first gridDim correspond to the first row)
---
--- The value of each entry is the number of cells (whose type is determined by the column) in the appropriate position (determined by the row)

module InVivo where

import System.Random
import Data.List
import Data.Sequence
import System.Directory
import System.IO
import Biofunctions

--- ---------------------
--- Constant parameters
--- ---------------------

--- Used throughout the program
 
iTumNumInVivo  = 10000 --- :: Num a => a - initial number of tumor cells (in vitro: to match Seki, this should be about 10,000)
chanceLoMHCvivo = 0.2 --- :: Num a => a --- the chance a tumor cell has low MHC expression

--- Used inVivoSimulation
iFCNumInVivo = 50000 --- :: Integral a => a - initial number of fasL CTLs in in vivo experiment
iPCNumInVivo = 150000 --- :: Integral a => a - initial number of perforin CTLs in in vivo experiment
iFNNumInVivo = 0 --- :: Integral a => a - initial number of fasL NKLs in in vivo experiment
iPNNumInVivo = 0 --- :: Integral a => a - initial number of perforin NKLs in in vivo experiment
repsInVivo = 5 --- :: Integral a => a - number of times the experiment is repeated with the same settings


--- Used in inVivoTimeStep
timeStepsInVivo = 200 --- :: Num a => a - the number of time steps for an in vivo simulation
speedBoundInVivo = 4 --- :: Num a => a - the upper limit of the number of movement iterations in a time step
probTumRep = 0.04 --- :: Fractional a => a - probability of a tumor cell completing the process of replication in a single timestep
avFNadd = 1 --- :: Integral a => a - average number of fasL NKLs to be added at the end of each (in vivo) time step
avPNadd = 10 --- :: Integral a => a - average number of perforin NKLs to be be added at the end of each (in vivo) time step
addRange = 5 --- :: Integral a => a - half the range of numbers of each type of cell to add
addRangeFN = minimum [addRange, avFNadd]
addRangePN = minimum [addRange, avPNadd]
addMultFC = 2 --- :: Integral a => a - the number of tumor cells lysed within a certain period of time multiplied by *this number* gives the average number of FasCTLs to add during a timestep
addMultPC = 6 --- :: Integral a => a - the number of tumor cells lysed within a certain period of time multiplied by *this number* gives the average number of Perforin CTLs to add during a timestep
ctlProdTime = 4 --- :: Integral a => a - the number of timesteps it takes to produce and recruit a new ctl
memory = 8 --- :: Integral a => a - the number of tumor cells lysed as many as "memory" time steps ago can influence the number of new CTL's produced and recruited


--- in positionLst and availablePos
maxTumNumVivo = 1 :: Num a => a --  the maximum number of tumor cells that can be in a grid element


--- Used in createDistSeq
distDepth = 3 -- :: Integral a => a -- (NOT MUCH FLEXIBILITY TO CHANG THIS) maximum number of grid elements a CTL can detect the proximity of a tumor cell 

--- Concerning distribution of effector cells
distributeTypeOrig = 1 -- :: Int --- 1 for avoid tumor; 2 for edge of tumor; 3 for edge of grid, 4 for uniformly throughout grid
distributeTypeAdd = 1 -- :: Int --- 1 for avoid tumor; 2 for edge of tumor; 3 for edge of grid, 4 for uniformly throughout grid


--- -----------------------
--- main procedure
--- -----------------------

inVivo :: IO ()
--- Effects
---   Write position data (a file for each exponent and each point in time), average percent lysis data, and variance data
---   Print progress
---   Print he maximum number of leucocytes in the same grid element as a tumor cells after a movement phase 
inVivo = do 
--  gen <- getStdGen  --- If we want "real" randomization
  let gen = (mkStdGen 100) --- If we want to fix a "seed"
  inVivoSimulation gen 
  
  

inVivoSimulation :: StdGen -> IO ()
--- Inputs
---   gen :: StdGen - a random generator
--- Global constants used:
---   iFNNumInVivo :: Integral a => a - initial number of fasL NKLs in in vivo experiment
---   iPNNumInVivo :: Integral a => a - initial number of perforin NKLs in in vivo experiment
---   iFCNumInVivo :: Integral a => a - initial number of fasL CTLs in in vivo experiment
---   iPCNumInVivo :: Integral a => a - initial number of perforin CTLs in in vivo experiment
--- Effects:
---   Write position data (a file for each run and point in time)
---   Write distance data (a file for each run and point in time)
---   Write detailed stats data (a file for each run)
---   Write counts data (to more easily be read by Matlab)
---   Write average statistics
---   Write average counts (to more easily be read by Matlab)
---   Print progress
inVivoSimulation gen = do
  createDirectoryIfMissing True "dataPosition"  -- True could be changed to False.  It is only relevant if we specify a parent directory
  createDirectoryIfMissing True "dataDistance"  -- True could be changed to False.  It is only relevant if we specify a parent directory
  createDirectoryIfMissing True "dataStats"  -- True could be changed to False.  It is only relevant if we specify a parent directory
  createDirectoryIfMissing True "dataCounts"
  let inVivoSimRec :: StdGen -> Int -> [[Int]] -> [[Int]] -> [[Composition]] -> IO ([[Int]], [[Int]], [[Composition]])
      --- genrec - random generator 
      --- repsrec - the number of runs left
      --- lyslstlst - a list of lists: one list for each run of the cumulative lysis by the end of each time step
      --- tumlstlst - a list of lists: one list for each run of number of the cumulative tumor cells added by the beginning of each time step
      --- countslstlst - one list for each run of the grid composition of cell counts
      inVivoSimRec genrec repsrec lyslstlst tumlstlst countslstlst = do
       let repcur = repsInVivo-repsrec+1   --- current repetition number (this the "repcur-th" run of the experiment)
       let emptyGrid = Data.Sequence.replicate gridSize (Comp 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)   --- create a data sequence of eighteen-tuples
       let distSeq = Data.Sequence.replicate gridSize (distDepth +1) --- create data sequence of 4's for tumor distance grid
       putStrLn $ "Starting to fill tumor positions."
       let fstPos = (div gridSize 2) + (div gridDim 2)
           (randn, newgen) = randomR (0::Double, 1::Double) genrec  --- determine whether first tumor cell has high MHC
           !oneTumGrid = if randn > chanceLoMHCvivo  --- place the first tumor cell near the middle of the grid
                         then update fstPos (adjustTumHiMHC (\x->x+1) (index emptyGrid fstPos)) emptyGrid
                         else update fstPos (adjustTumLoMHC (\x->x+1) (index emptyGrid fstPos)) emptyGrid
           !oneDistSeq = addTumDistSeq distSeq fstPos
           frstAvailLst = if maxTumNumVivo == 1       --- for where the first round of tumor cells can go after the initial cell
                          then nbrPos fstPos
                          else fstPos:(nbrPos fstPos)
           !(withTumGrid, nDistSeq, genFC) = addTumsContiguous (iTumNumInVivo-1) oneTumGrid oneDistSeq frstAvailLst (Prelude.length frstAvailLst) newgen 
       putStrLn $ "Finished filling in tumor cells."
       let !(withFCgrid, genPC) 
             | distributeTypeOrig <= 1 = distributeEffAvoidTum withTumGrid nDistSeq iFCNumInVivo (adjustFCa (\x -> x+1)) genFC
             | distributeTypeOrig <= 2 = distributeEffTumEdge withTumGrid nDistSeq iFCNumInVivo (adjustFCa (\x -> x+1)) genFC
             | distributeTypeOrig <= 3 = distributeEffEdge withTumGrid iFCNumInVivo (adjustFCa (\x -> x+1)) genFC
             | otherwise = distributeEff withTumGrid iFCNumInVivo (adjustFCa (\x -> x+1)) genFC
       putStrLn $ "Finished filling in Fas CTLs."
       let !(withPCgrid, genFN) 
             | distributeTypeOrig <= 1 = distributeEffAvoidTum withFCgrid nDistSeq iPCNumInVivo (adjustPCa (\x -> x+1)) genPC
             | distributeTypeOrig <= 2 = distributeEffTumEdge withFCgrid nDistSeq iPCNumInVivo (adjustPCa (\x -> x+1)) genPC
             | distributeTypeOrig <= 3 = distributeEffEdge withFCgrid iPCNumInVivo (adjustPCa (\x -> x+1)) genPC
             | otherwise = distributeEff withFCgrid iPCNumInVivo (adjustPCa (\x -> x+1)) genPC
       putStrLn $ "Finished filling in perforin CTLs."
       let !(withFNgrid, genPN) 
             | distributeTypeOrig <= 1 = distributeEffAvoidTum withPCgrid nDistSeq iFNNumInVivo (adjustFNa (\x -> x+1)) genFN
             | distributeTypeOrig <= 2 = distributeEffTumEdge withPCgrid nDistSeq iFNNumInVivo (adjustFNa (\x -> x+1)) genFN
             | distributeTypeOrig <= 3 = distributeEffEdge withPCgrid iFNNumInVivo (adjustFNa (\x -> x+1)) genFN
             | otherwise = distributeEff withPCgrid iFNNumInVivo (adjustFNa (\x -> x+1)) genFN
       putStrLn $ "Finished filling in Fas NKLs."
       let !(withAllgrid, genTS) 
             | distributeTypeOrig <= 1 = distributeEffAvoidTum withFNgrid nDistSeq iPNNumInVivo (adjustPNa (\x -> x+1)) genPN
             | distributeTypeOrig <= 2 = distributeEffTumEdge withFNgrid nDistSeq iPNNumInVivo (adjustPNa (\x -> x+1)) genPN
             | distributeTypeOrig <= 3 = distributeEffEdge withFNgrid iPNNumInVivo (adjustPNa (\x -> x+1)) genPN
             | otherwise = distributeEff withFNgrid iPNNumInVivo (adjustPNa (\x -> x+1)) genPN
       putStrLn $ "Finished filling in perforin NKLs."
       --- Write position and distance data
       if repcur == 1
       then do
         handlePos <- openFile ("dataPosition/inVivoRun" ++ (show repcur) ++ "Time" ++ (show 0) ++ "Pos.dat") WriteMode
         writePosData withAllgrid 0 handlePos
         putStrLn $ "Finished writing position data for run " ++ (show repcur) ++ " time " ++ (show 0) ++ "."
         hClose handlePos
         handleDist <- openFile ("dataDistance/inVivoRun" ++ (show repcur) ++ "Time" ++ (show 0) ++ "Dist.dat") WriteMode
         writeDistData nDistSeq 0 handleDist
         hClose handleDist
       else return ()
       --- Write stats data
       let iEffNum = iFNNumInVivo + iPNNumInVivo + iFCNumInVivo + iPCNumInVivo
       handleStats <- openFile ("dataStats/inVivoRun" ++ (show repcur) ++ "Stats.dat") WriteMode
       hPutStrLn handleStats $ "Time-step: 0"
       hPutStrLn handleStats $ "  The original number of tumor cells in run " ++ (show repcur) ++ " is: " ++ (show iTumNumInVivo)
       hPutStrLn handleStats $ "  The original number of effector cells in run " ++ (show repcur) ++ " is: " ++ (show iEffNum) 
       hPutStrLn handleStats $ "  The original effector cell composition in run " ++ (show repcur) ++ " is: " ++ (show iFCNumInVivo) ++ " FC's, " ++ (show iPCNumInVivo) ++ " PC's, " ++ (show iFNNumInVivo) ++ " FN's, " ++ (show iPNNumInVivo) ++ " PN's"
       (lysisCountList, lysisCumList, tumAddList, tumCumList, countsList, tumNum, effNum, wCellCount, nHandleStats, genNxt) <- inVivoTimeStep withAllgrid nDistSeq repcur 0 [] [] [] [] [] iTumNumInVivo iEffNum 0 handleStats genTS 
       hPutStrLn nHandleStats $ "The total number of tumor cells lysed in run " ++ (show repcur) ++ " is " ++ (show $ sum lysisCountList) ++ "."
       hPutStrLn nHandleStats $ "The total number of tumor cells accumulated in run " ++ (show repcur) ++ " is " ++ (show $ iTumNumInVivo + (sum tumAddList)) ++ "."
       hPutStrLn nHandleStats $ "The total number of tumor cells remaining in run " ++ (show repcur) ++ " is " ++ (show tumNum) ++ "."
       hPutStrLn nHandleStats $ "The total number of effector cells remaining in run " ++ (show repcur) ++ " is " ++ (show effNum) ++ "." 
       hPutStrLn nHandleStats $ "The maximum number of effector cells in the same grid element as a tumor cells after a movement phase in run " ++ (show repcur) ++ " is " ++ (show wCellCount) ++ "."
       hClose handleStats
       --- write cell counts to file
       handleCountsData <- openFile ("dataCounts/inVivoRun" ++ (show repcur) ++ "Counts.dat") WriteMode
       writeCounts handleCountsData countsList
       hClose handleCountsData
       --- recursion 
       if repsrec > 1
       then inVivoSimRec genNxt (repsrec-1) (lysisCumList:lyslstlst) (tumCumList:tumlstlst) (countsList:countslstlst)
       else return ((lysisCumList:lyslstlst), (tumCumList:tumlstlst), (countsList:countslstlst))
  (lyslstlsttot, tumlstlsttot, countslstlsttot) <- inVivoSimRec gen repsInVivo [] [] []
  handleAveStats <- openFile ("dataStats/aveStats.dat") WriteMode
  handleAveCounts <- openFile ("dataCounts/aveCounts.dat") WriteMode
  writeAveStats handleAveStats handleAveCounts lyslstlsttot tumlstlsttot countslstlsttot
  hClose handleAveStats      
  hClose handleAveCounts

writeCounts :: Handle -> [Composition] -> IO ()
--- Inputs:
---   handle :: Handle - a handle for InVivoRun" ++ (show run) ++ "Counts.dat
---   compLst :: [Composition] - a list of grid compositions of cell counts, one for each time step (last time step first)
--- Effect:
---   write grid compositions for each time step in handle (any of files inVivoRun1Stats.dat, inVivoRun2Stats.dat, etc)
writeCounts handle compLst = 
    if (not (Data.List.null compLst))
    then do
      writeCounts handle (tail compLst)
      let thiscomp = (head compLst)
      let (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = thiscomp
      hPutStrLn handle $ (show tumLoMHC) ++ " " ++ (show tumHiMHC) ++ " " ++ (show fca) ++ " " ++ (show fcu) ++ " " ++ (show fc1) ++ " " ++ (show fc2) ++ " " ++ (show pca) ++ " " ++ (show pcu) ++ " " ++ (show pc1) ++ " " ++ (show pc2) ++ " " ++ (show fna) ++ " " ++ (show fnu) ++ " " ++ (show fn1) ++ " " ++ (show fn2) ++ " " ++ (show pna) ++ " " ++ (show pnu) ++ " " ++ (show pn1) ++ " " ++ (show pn2)
    else return ()

writeAveStats :: Handle -> Handle -> [[Int]] -> [[Int]] -> [[Composition]] -> IO ()  --- Write average lysis data
--- Inputs:
---   handleStats :: Handle - a handle for aveStats.dat
---   handleCounts :: Handle - a handle for counts
---   lysLstLst :: [[Int]] - a list for each run (last run appearing first) of a list of cumulative lysis at each time step (with last time step appearing first)
---   tumLstLst :: [[Int]] - a list for each run (last run appearing first) of a list of cumulative tumors cells added to the grid (with last time step appearing first)
---   countsLstLst :: [[Composition]] - a list for each run (last run appearing first) of a list grid compositions of cell counts
--- Effects:
---   Write average statistics, variance, and standard deviations in handleStats (file \dataStats\aveStats.dat) in a way that is easy for people to read
---   Write average statistics and standard deviations in file for handleCounts (file \dataConuts\aveCounts.dat) in such a way that may be easy for matLab to read
writeAveStats handleStats handleCounts lysLstLst tumLstLst countsLstLst = 
    if (not (Data.List.null $ head lysLstLst))
    then do
      writeAveStats handleStats handleCounts (map tail lysLstLst) (map tail tumLstLst) (map tail countsLstLst)
      let lysList = (map head lysLstLst)
          tumList = (map head tumLstLst)
          countsList = (map head countsLstLst)
          aveLys = (fromIntegral $ sum lysList)/(fromIntegral $ Prelude.length lysList)
          aveTum = (fromIntegral $ sum tumList)/(fromIntegral $ Prelude.length tumList)
          aplo = (fromIntegral 100)*aveLys/(fromIntegral iTumNumInVivo)  --- average lysis divided by original tumor count
          aplc = (fromIntegral 100)*aveLys/aveTum  --- average lysis divided by average cumulative tumor count (tumors accumulated by beginning of time step)
          varo = let lnth = (fromIntegral $ Prelude.length lysList)  --- sample variance of percent lysis original
                  in if lnth == 1
                     then 0
                     else (sum $ map (\x -> (100*(fromIntegral x)/(fromIntegral iTumNumInVivo) - aplo)^2) lysList)/(lnth-1)
          varc = let lnth = (fromIntegral $ Prelude.length lysList)  --- sample variance of percent lysis original
                  in if lnth == 1
                     then 0
                     else (sum $ Data.List.zipWith (\x y -> (100*(fromIntegral x)/(fromIntegral y) - aplc)^2) lysList tumList)/(lnth - 1)
          aveTumLoMHC = (fromIntegral $ sum (map (\x -> (countTumLoMHC x)) countsList))/(fromIntegral $ Prelude.length countsList)
          varTumLoMHC = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
                         in if lnth == 1
                            then 0
                            else (sum $ map (\x -> ((fromIntegral $ countTumLoMHC x) - aveTumLoMHC)^2) countsList)/(lnth-1)
          aveTumHiMHC = (fromIntegral $ sum (map (\x -> (countTumHiMHC x)) countsList))/(fromIntegral $ Prelude.length countsList)
          varTumHiMHC = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
                         in if lnth == 1
                            then 0
                            else (sum $ map (\x -> ((fromIntegral $ countTumHiMHC x) - aveTumHiMHC)^2) countsList)/(lnth-1)  
          aveCtlPreRec = (fromIntegral $ sum (map (\x -> (countFCa x) + (countFCu x) + (countPCa x) + (countPCu x)) countsList))/(fromIntegral $ Prelude.length countsList) 
          aveCtlPostRec = (fromIntegral $ sum (map (\x -> (countFC1 x) + (countFC2 x) + (countPC1 x) + (countPC2 x)) countsList))/(fromIntegral $ Prelude.length countsList) 
          aveNklPreRec = (fromIntegral $ sum (map (\x -> (countFNa x) + (countFNu x) + (countPNa x) + (countPNu x)) countsList))/(fromIntegral $ Prelude.length countsList) 
          aveNklPostRec = (fromIntegral $ sum (map (\x -> (countFN1 x) + (countFN2 x) + (countPN1 x) + (countPN2 x)) countsList))/(fromIntegral $ Prelude.length countsList) 
          varCtlPreRec = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
              in if lnth == 1
                 then 0
                 else (sum $ map (\x -> ((fromIntegral  ((countFCa x) + (countFCu x) + (countPCa x) + (countPCu x))) - aveCtlPreRec)^2) countsList)/(lnth-1) 
          varCtlPostRec = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
              in if lnth == 1
                 then 0
                 else (sum $ map (\x -> ((fromIntegral  ((countFC1 x) + (countFC2 x) + (countPC1 x) + (countPC2 x))) - aveCtlPostRec)^2) countsList)/(lnth-1)
          varNklPreRec = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
              in if lnth == 1
                 then 0
                 else (sum $ map (\x -> ((fromIntegral  ((countFNa x) + (countFNu x) + (countPNa x) + (countPNu x))) - aveNklPreRec)^2) countsList)/(lnth-1)                 
          varNklPostRec = let lnth = (fromIntegral $ Prelude.length countsList)  --- sample variance
              in if lnth == 1
                 then 0
                 else (sum $ map (\x -> ((fromIntegral  ((countFN1 x) + (countFN2 x) + (countPN1 x) + (countPN2 x))) - aveNklPostRec)^2) countsList)/(lnth-1)                
      hPutStrLn handleStats $ " Time-step: " ++ (show (Prelude.length (head lysLstLst)))
      hPutStrLn handleStats $ "   Average percent lysis wrt original tumor count: " ++ (show aplo)
      hPutStrLn handleStats $ "   Average percent lysis wrt number of tumors accumulated by beginning of time step: " ++ (show aplc)
      hPutStrLn handleStats $ "   Variance of percent lysis wrt original tumor count: " ++ (show varo) 
      hPutStrLn handleStats $ "   Variance of percent lysis wrt tumors accumulated by beginning of time step: " ++ (show varc)
      hPutStrLn handleStats $ "   Average: loMHC hiMHC ctlprerec ctlpostrec nklprerec nkl postrec: " ++ (show aveTumLoMHC) ++ ", " ++ (show aveTumHiMHC) ++ ", " ++ (show aveCtlPreRec) ++ ", " ++ (show aveCtlPostRec) ++ ", " ++ (show aveNklPreRec) ++ ", " ++ (show aveNklPostRec)
      hPutStrLn handleStats $ "   Variance: loMHC hiMHC ctlprerec ctlpostrec nklprerec nkl postrec: " ++ (show varTumLoMHC) ++ ", " ++ (show varTumHiMHC) ++ ", " ++ (show varCtlPreRec) ++ ", " ++ (show varCtlPostRec) ++ ", " ++ (show varNklPreRec) ++ ", " ++ (show varNklPostRec)
      hPutStrLn handleStats$ "   Standard deviation: loMHC hiMHC ctlprerec ctlpostrec nklprerec nklpostrec: " ++ (show $ sqrt varTumLoMHC) ++ ", " ++ (show $ sqrt varTumHiMHC) ++ ", " ++ (show $ sqrt varCtlPreRec) ++ ", " ++ (show $ sqrt varCtlPostRec) ++ ", " ++ (show $ sqrt varNklPreRec) ++ ", " ++ (show $ sqrt varNklPostRec)
      hPutStrLn handleCounts $ (show aveTumLoMHC) ++ " " ++ (show aveTumHiMHC) ++ " " ++ (show aveCtlPreRec) ++ " " ++ (show aveCtlPostRec) ++ " " ++ (show aveNklPreRec) ++ " " ++ (show aveNklPostRec)
      hPutStrLn handleCounts  $ (show $ sqrt varTumLoMHC) ++ " " ++ (show $ sqrt varTumHiMHC) ++ " " ++ (show $ sqrt varCtlPreRec) ++ " " ++ (show $ sqrt varCtlPostRec) ++ " " ++ (show $ sqrt varNklPreRec) ++ " " ++ (show $ sqrt varNklPostRec)
    else return ()
         
 
inVivoTimeStep :: Seq Composition -> Seq Int -> Int -> Int -> [Int] -> [Int] -> [Int] -> [Int] -> [Composition] -> Int -> Int -> Int -> Handle -> StdGen -> IO ([Int], [Int], [Int], [Int], [Composition], Int, Int, Int, Handle, StdGen)          
--- Inputs:
---   grid :: Seq Composition - the input grid
---   distSeq :: Seq Int - distance sequence
---   run :: Int - the run number (1 if it is the 1st, 2 if it is the 2nd)
---   time :: Int - how many time steps have already ellapsed
---   lysCountLst :: [Int] - a list of the number of tumor cells lysed in each time step
---   lysCumLst :: [Int] - a list of the total number of tumor cells lysed by the end of the time step
---   tumAddLst :: [Int] - a list of the number of tumor cells added at each time step
---   tumCumLst :: [Int] - a list of the total number of tumor cells added by the beginning of the time step (including the initial tumor count)
---   countslst :: [Composition] - a list of grid (rather than grid element) compositions of cell counts, one for each time step
---   tumNum :: Int - the number of tumor cells
---   wCellNum :: Int - the number of effector cells
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   statsHandle :: Handle - a handle to put lysis counts and effector cell counts at each stage
---   gen :: StdGen - a random generator
--- Global constants:
---   timeStepsInVivo :: Num a => a - the number of time steps for an in vivo simulation
---   probTumRep :: Fractional a => a - probability of a tumor cell completing the process of replication in a single timestep
---   avFNadd :: Integral a => a - average number of fasL NKLs to be added at the end of each (in vivo) time step
---   avPNadd :: Integral a => a - average number of perforin NKLs to be be added at the end of each (in vivo) time step
---   addRange :: Integral a => a - half the range of numbers of each type of cell to add
---   addRangeFN = minimum [addRange, avFNadd]
---   addRangePN = minimum [addRange, avPNadd]
---   addMultFC :: Integral a => a - the number of tumor cells lysed within a certain period of time multiplied by *this number* gives the average number of FasCTLs to add during a timestep
---   addMultPC :: Integral a => a - the number of tumor cells lysed within a certain period of time multiplied by *this number* gives the average number of Perforin CTLs to add during a timestep
---   ctlProdTime :: Integral a => a - the number of timesteps it takes to produce and recruit a new ctl
---   memory :: Integral a => a - the number of tumor cells lysed as many as "memory" time steps ago can influence the number of new CTL's produced and recruited  
--- Outputs:
---   :: [Int] - an updated list of the number of tumor cells lysed in each time step
---   :: [Int] - an updated list of cumulative lysis at end of each time step
---   :: [Int] - an updated list of number of tumor cells added in each time step
---   :: [Int] - an updated list of the cumulative number of tumor cells added by the beginning of each-time step (including time 0 in sum)
---   :: [Composition] - an update list of grid (not grid element) compositions of cell counts
---   :: Int - the number of tumor cells at the end of the time step
---   :: Int - the number of effector cells at the end of the time step
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   :: StdGen - an updated randam generator
--- Effects:
---   Write position data corresponding to the grid at the end of this time step
---   Write distance data corresponding to the grid at the end of this time step
---   Print program runtime progress
inVivoTimeStep grid distSeq run time lysCountLst lysCumLst tumAddLst tumCumLst countslst tumNum wCellNum wCellCount statsHandle gen = 
    if time == timeStepsInVivo  --- time starts at 0 and then increments with each call
    then return (lysCountLst, lysCumLst, tumAddLst, tumCumLst, countslst, tumNum, wCellNum, wCellCount, statsHandle, gen)
    else do
      -- movement
      let (times,timesgen) = randomR (1 :: Int, speedBoundInVivo :: Int) gen
      let !(movedGrid, genLys) = moveCellsInVivo grid distSeq times timesgen 
      putStrLn $ "Finished a movement phase.  " -- ++ (show gen)
      -- lysis
      let (nGrid, distSeqAftLys, nLysCount, nwCellRemCount, nwCellCount, nGen) = inVivoLysing movedGrid distSeq 0 0 0 wCellCount genLys
      putStrLn "Finished a lysis phase.  "
      -- new tumor cells
      let availablePosAftLys = newAvailableLst nGrid distSeqAftLys
      let (numNewTums,genRep) = randomBinom (tumNum-nLysCount) (probTumRep * ((fromIntegral 1) - (fromIntegral $ tumNum-nLysCount)/(fromIntegral gridSize))) nGen
      putStrLn $ "The number of new tumor cells is " ++ (show numNewTums)
      let !(withRepTumsGrid,nDistSeq,genFCrcr) = addTumsContiguous numNewTums nGrid distSeqAftLys availablePosAftLys (Prelude.length availablePosAftLys) genRep
      putStrLn "Finished tumor replication."
      -- new effector cells
      let avFCadd = addMultFC * (sum $ Data.List.take (memory - ctlProdTime) $ Data.List.drop (ctlProdTime) lysCountLst)
          avPCadd = addMultPC * (sum $ Data.List.take (memory - ctlProdTime) $ Data.List.drop (ctlProdTime) lysCountLst)
          addRangeFC = minimum [avFCadd, addRange]
          addRangePC = minimum [avPCadd, addRange] 
          (newFCnum, genPCrcr) = randomR (avFCadd - addRangeFC::Int, avFCadd + addRangeFC::Int) genFCrcr
          (newPCnum, genFNrcr) = randomR (avPCadd - addRangePC::Int, avPCadd + addRangePC::Int) genPCrcr
          (newFNnum, genPNrcr) = randomR (avFNadd - addRangeFN ::Int, avFNadd + addRangeFN) genFNrcr
          (newPNnum, genFCadd) = randomR (avPNadd - addRangePN ::Int, avPNadd + addRangePN) genPNrcr
      putStrLn $ "The number of new FC, new PC, new FN, new PN are " ++ (show newFCnum) ++ ", " ++ (show newPCnum) ++ ", " ++ (show newFNnum) ++ ", " ++ (show newPNnum)
      let !(withFCgrid, genPCadd) 
               | distributeTypeAdd <= 1 = distributeEffAvoidTum withRepTumsGrid nDistSeq newFCnum (adjustFCa (\x -> x+1)) genFCadd
               | distributeTypeAdd <= 2 = distributeEffTumEdge withRepTumsGrid nDistSeq newFCnum (adjustFCa (\x -> x+1)) genFCadd
               | distributeTypeAdd <= 3 = distributeEffEdge withRepTumsGrid newFCnum (adjustFCa (\x -> x+1)) genFCadd
               | otherwise = distributeEff withRepTumsGrid newFCnum (adjustFCa (\x -> x+1)) genFCadd
          !(withPCgrid, genFNadd)
               | distributeTypeAdd <= 1 = distributeEffAvoidTum withFCgrid nDistSeq newPCnum (adjustPCa (\x -> x+1)) genPCadd
               | distributeTypeAdd <= 2 = distributeEffTumEdge withFCgrid nDistSeq newPCnum (adjustPCa (\x -> x+1)) genPCadd
               | distributeTypeAdd <= 3 = distributeEffEdge withFCgrid newPCnum (adjustPCa (\x -> x+1)) genPCadd
               | otherwise = distributeEff withFCgrid newPCnum (adjustPCa (\x -> x+1)) genPCadd
          !(withFNgrid, genPNadd) 
               | distributeTypeAdd <= 1 = distributeEffAvoidTum withPCgrid nDistSeq newFNnum (adjustFNa (\x -> x+1)) genFNadd
               | distributeTypeAdd <= 2 = distributeEffTumEdge withPCgrid nDistSeq newFNnum (adjustFNa (\x -> x+1)) genFNadd
               | distributeTypeAdd <= 3 = distributeEffEdge withPCgrid newFNnum (adjustFNa (\x -> x+1)) genFNadd
               | otherwise = distributeEff withPCgrid newFNnum (adjustFNa (\x -> x+1)) genFNadd
          !(withAllgrid, genRet) 
               | distributeTypeAdd <= 1 = distributeEffAvoidTum withFNgrid nDistSeq newPNnum (adjustPNa (\x -> x+1)) genPNadd 
               | distributeTypeAdd <= 2 = distributeEffTumEdge withFNgrid nDistSeq newPNnum (adjustPNa (\x -> x+1)) genPNadd
               | distributeTypeAdd <= 3 = distributeEffEdge withFNgrid newPNnum (adjustPNa (\x -> x+1)) genPNadd
               | otherwise = distributeEff withFNgrid newPNnum (adjustPNa (\x -> x+1)) genPNadd
      putStrLn "Finished distributing effector cells."
      -- write to files
      let nTumNum = tumNum - nLysCount + numNewTums
      let nEffNum = wCellNum - nwCellRemCount + (newFCnum + newPCnum + newFNnum + newPNnum)
      --- write position and distance data
      if run == 1
      then do 
        handlePos <- openFile ("dataPosition/inVivoRun" ++ (show run) ++ "Time" ++ (show $ time+1) ++ "Pos.dat") WriteMode
        writePosData withAllgrid 0 handlePos
--      putStrLn $ "Finished writing position data for time " ++ (show $ time+1) ++ ".\n"
        hClose handlePos
        handleDist <- openFile ("dataDistance/inVivoRun" ++ (show run) ++ "Time" ++ (show $ time+1) ++ "Dist.dat") WriteMode
        writeDistData nDistSeq 0 handleDist
        hClose handleDist
      else return ()
      -- write to stats file
      hPutStrLn statsHandle $ ("Run " ++ (show run) ++ " time-step: " ++ (show $ time+1))
      hPutStrLn statsHandle $ ("   New tumor count: " ++ (show nTumNum))
      hPutStrLn statsHandle $ ("   Tumor cells cumulated by beginning of time step " ++ (show $ iTumNumInVivo + (sum tumAddLst)))
      hPutStrLn statsHandle $ ("  Lysis: " ++ (show $ nLysCount))
      hPutStrLn statsHandle $ ("  Accumulated lysis: " ++ (show $ sum (nLysCount:lysCountLst)))
      hPutStrLn statsHandle $ ("  Number of tumor cells added: " ++ (show numNewTums))
      hPutStrLn statsHandle $ ("   New effector cell count: " ++ (show nEffNum))
      hPutStrLn statsHandle $ ("  Number of effector cells becoming ineffective: " ++ (show nwCellRemCount))
      hPutStrLn statsHandle $ ("  New effector cells (total): " ++ (show (newFCnum + newPCnum + newFNnum + newPNnum)))
      hPutStrLn statsHandle $ ("  New effector cells by type: " ++ (show newFCnum) ++ " FC's, " ++ (show newPCnum) ++ " PC's, " ++ (show newFNnum) ++ " FN's, " ++ (show newPNnum) ++ " PN's")--      hPutStrLn statsHandle $ (" Max effector cells in single gridpoint: " ++ (show nwCellCount))
      putStrLn $ "Finished Run " ++ (show run) ++ " time " ++ (show $ time+1) ++ ".\n"
      -- recursion
      let newcounts = seqToCounts withAllgrid  -- updated grid composition of cell counts
      inVivoTimeStep withAllgrid nDistSeq run (time+1) (nLysCount:lysCountLst) ((sum (nLysCount:lysCountLst)):lysCumLst) (numNewTums:tumAddLst) ((iTumNumInVivo + sum tumAddLst):tumCumLst) (newcounts:countslst) nTumNum nEffNum nwCellCount statsHandle genRet
                                      
seqToCounts :: Seq Composition -> Composition
--- Inputs:
---   grid
--- Output:
---   countsComp
seqToCounts grid = seqToCountsRec 0 (Comp 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) where
  seqToCountsRec :: Int -> Composition -> Composition
  --- Inputs
  ---   pos :: Int - position in grid to add to counts
  ---   countsrec :: Composition - composition of the current tallies
  --- Output:
  ---   :: Composition - composition of cell counts
  seqToCountsRec pos countsrec =
    if pos < gridSize
    then
      let compHere = index grid pos 
          (Comp tumLoMHCorig tumHiMHCorig fcAorig fcUorig fc1orig fc2orig pcAorig pcUorig pc1orig pc2orig fnAorig fnUorig fn1orig fn2orig pnAorig pnUorig pn1orig pn2orig) = compHere
          (Comp tumLoMHCTal tumHiMHCTal fcATal fcUTal fc1Tal fc2Tal pcATal pcUTal pc1Tal pc2Tal fnATal fnUTal fn1Tal fn2Tal pnATal pnUTal pn1Tal pn2Tal) = countsrec
          countsnew = (Comp (tumLoMHCorig+tumLoMHCTal) (tumHiMHCorig+tumHiMHCTal) (fcAorig+fcATal) (fcUorig+fcUTal) (fc1orig+fc1Tal) (fc2orig+fc2Tal) (pcAorig+pcATal) (pcUorig+pcUTal) (pc1orig+pc1Tal) (pc2orig+pc2Tal) (fnAorig+fnATal) (fnUorig+fnUTal) (fn1orig+fn1Tal) (fn2orig+fn2Tal) (pnAorig+pnATal) (pnUorig+pnUTal) (pn1orig+pn1Tal) (pn2orig+pn2Tal))
        in seqToCountsRec (pos+1) countsnew
    else countsrec
          



                                                                        
moveCellsInVivo :: Seq Composition -> Seq Int -> Int -> StdGen -> (Seq Composition, StdGen)
--- Inputs:
---   gridOrig :: Seq Composition - an input grid of cells to move
---   distSeq :: Seq Int - distance sequence (indicating for each position how far to nearest tumor cell)
---   times :: Int - the number of movement iterations left
---   gen :: StdGen - a random generator
--- Global constants:
---   gridSize 
--- Outputs:
---   :: Seq Composition - an updated grid with CTL's having moved
---   :: StdGen - an updated random generator
moveCellsInVivo gridOrig distSeq times gen = if times == 0 then (gridOrig, gen) else moveCellsRec gridOrig 0 gen where 
    moveCellsRec :: Seq Composition -> Int -> StdGen -> (Seq Composition, StdGen)
    moveCellsRec gridRec pos genRec = 
        if pos == gridSize  --- pos starts at 0 and increments with each call
        then moveCellsInVivo gridRec distSeq (times-1) genRec
        else let compOrig = index gridOrig pos 
                 (Comp tumLoMHCorig tumHiMHCorig fcAorig fcUorig fc1orig fc2orig pcAorig pcUorig pc1orig pc2orig fnAorig fnUorig fn1orig fn2orig pnAorig pnUorig pn1orig pn2orig) = compOrig 
             in if ((countEff compOrig) == 0)
                then moveCellsRec gridRec (pos+1) genRec
                else if (countTum compOrig >0)
                     then let (nkNotMove, genRec1) = randomR (0::Int,1::Int) genRec  --- Coin toss
                            in if (nkNotMove >0) || (times > 1) --- NK cells can only move during last movement phase of time step
                               then moveCellsRec gridRec (pos+1) genRec1
                               else let compCur = index gridRec pos 
                                        (Comp tumLoMHCcur tumHiMHCcur fcAcur fcUcur fc1cur fc2cur pcAcur pcUcur pc1cur pc2cur fnAcur fnUcur fn1cur fn2cur pnAcur pnUcur pn1cur pn2cur) = compCur
                                        compNew = (Comp tumLoMHCcur tumHiMHCcur fcAcur fcUcur fc1cur fc2cur pcAcur pcUcur pc1cur pc2cur fnAcur (fnUcur - fnUorig) fn1cur (fn2cur - fn2orig) pnAcur (pnUcur - pnUorig) pn1cur (pn2cur - pn2orig))
                                        !gridNew = update pos compNew gridRec -- all original NK cells move out
                                        fcLst = []    --- no CTLs move in presence of tumor cell
                                        pcLst = []    --- no CTLs move in presence of tumor cell 
                                        nbrs = nbrPos pos --- finding positions for nkl that move out
                                        nbrNum = Prelude.length nbrs 
                                        (fnLst, genPN) = randomRsGen (0::Int, nbrNum-1::Int) ((countFNu compOrig)+(countFN2 compOrig)) [] genRec1
                                        (pnLst, genRet) = randomRsGen (0::Int, nbrNum-1::Int) ((countPNu compOrig)+(countPN2 compOrig)) [] genPN
                                      in moveCellsSubRec (nbrNum-1) nbrs pos (fnLst, pnLst, fcLst, pcLst) gridNew genRet
                     else let compCur = index gridRec pos   --- in absence of tumor cell
                              (Comp tumLoMHCcur tumHiMHCcur fcAcur fcUcur fc1cur fc2cur pcAcur pcUcur pc1cur pc2cur fnAcur fnUcur fn1cur fn2cur pnAcur pnUcur pn1cur pn2cur) = compCur
                              compNew = (Comp tumLoMHCcur tumHiMHCcur (fcAcur - fcAorig) (fcUcur - fcUorig) (fc1cur - fc1orig) (fc2cur - fc2orig) (pcAcur - pcAorig) (pcUcur - pcUorig) (pc1cur - pc1orig) (pc2cur - pc2orig) (fnAcur - fnAorig) (fnUcur - fnUorig) (fn1cur - fn1orig) (fn2cur - fn2orig) (pnAcur - pnAorig) (pnUcur - pnUorig) (pn1cur - pn1orig) (pn2cur - pn2orig))
                              !gridNew = update pos compNew gridRec -- all original leukocytes move out
                              (fcLst, pcLst, genFN) = if (countFC compOrig)==0 && (countPC compOrig) == 0  --- create list of positions for CTLs to move to
                                                      then ([],[], genRec)    --- In case with no ctls, no need to look for ctl attracting neighbors
                                                      else let nbrDist = map (\x-> index distSeq x) (nbrPos pos) --- create list of distances for each neighbor
                                                               ctlNbrs =  (minimum nbrDist) `elemIndices` nbrDist --- list indices for the minimum distance value
                                                               ctlnbrsNum = Data.List.length ctlNbrs --- number of neighbors whose positions have minimum distance
                                                               (fcIndices, genPC) = randomRsGen (0::Int, ctlnbrsNum-1::Int) (countFC compOrig) [] genRec
                                                               (pcIndices, genrt) = randomRsGen (0::Int, ctlnbrsNum-1::Int) (countPC compOrig) [] genPC 
                                                               fcLst = map (\x-> ctlNbrs!!x) fcIndices --- list of indices of positions to move to
                                                               pcLst = map (\x-> ctlNbrs!!x) pcIndices
                                                     in (fcLst, pcLst, genrt)
                              nbrs = nbrPos pos --- finding positions for nkl that move out
                              nbrNum = Prelude.length nbrs 
                              (fnLst, genPN) = randomRsGen (0::Int, nbrNum-1::Int) (countFN compOrig) [] genFN
                              (pnLst, genRet) = randomRsGen (0::Int, nbrNum-1::Int) (countPN compOrig) [] genPN
                            in moveCellsSubRec (nbrNum-1) nbrs pos (fnLst, pnLst, fcLst, pcLst) gridNew genRet
    moveCellsSubRec :: Int -> [Int] -> Int -> ([Int], [Int], [Int], [Int]) -> Seq Composition -> StdGen -> (Seq Composition, StdGen)
    moveCellsSubRec indx posLst posHere (fnSubRec, pnSubRec, fcSubRec, pcSubRec) gridSubRec genSubRec = 
        if indx >= 0  --- indx starts with nbrNum-1 and is decremented each call 
        then let !(fnPart, fnRest) = Data.List.partition (==indx) fnSubRec
                 !(pnPart, pnRest) = Data.List.partition (==indx) pnSubRec
                 !(fcPart, fcRest) = Data.List.partition (==indx) fcSubRec
                 !(pcPart, pcRest) = Data.List.partition (==indx) pcSubRec
                 !posThere = posLst!!indx
                 (Comp tumLoMHCthere tumHiMHCthere fcAthere fcUthere fc1there fc2there pcAthere pcUthere pc1there pc2there fnAthere fnUthere fn1there fn2there pnAthere pnUthere pn1there pn2there) = index gridSubRec posThere
                 !fnAnew = fnAthere + (Prelude.length fnPart)
                 !pnAnew = pnAthere + (Prelude.length pnPart)
                 !fcAnew = fcAthere + (Prelude.length fcPart)
                 !pcAnew = pcAthere + (Prelude.length pcPart)
                 !gridNew = update posThere (Comp tumLoMHCthere tumHiMHCthere fcAnew fcUthere fc1there fc2there pcAnew pcUthere pc1there pc2there fnAnew fnUthere fn1there fn2there pnAnew pnUthere pn1there pn2there) gridSubRec
             in moveCellsSubRec (indx-1) posLst posHere (fnRest, pnRest, fcRest, pcRest) gridNew genSubRec
        else moveCellsRec gridSubRec (posHere+1) genSubRec



inVivoLysing :: Seq Composition -> Seq Int -> Int -> Int -> Int -> Int -> StdGen -> (Seq Composition, Seq Int, Int, Int, Int, StdGen)
--- Inputs:
---   grid :: Seq Composition - a position sequence of cell compositions
---   distSeq :: Seq Int - a distance sequence
---   pos :: Int - the position of the grid we are focusing on
---   lysCount :: Int - a count of the number of tumor cells that have been lysed so far
---   wCellRemCount :: Int - the number of effector cells removed this time step so far
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   gen :: StdGen - a random generator
--- Global constants
---   gridSize :: Num a => a - the number of grid elements
---   chanceCTLrg  :: Fractional a => a - the chance a CTL recognizes a tumor cell
---   chanceNKLrg  :: Fractional a => a - the chance a NKL recognizes a tumor cell 
--- Outputs
---   :: Seq Composition - updated position sequence of cell compositions
---   :: Seq Int - updated distance sequence
---   :: Int - tumor lysis count
---   :: Int - the number of effector cells removed this time step so far
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   :: StdGen - updated random generator
inVivoLysing grid distSeq pos lysCount wCellRemCount wCellCount gen = 
    if pos == gridSize  --- pos starts at 0 and is incremented each call
    then (grid, distSeq, lysCount, wCellRemCount, wCellCount, gen)
    else
        let origComp = index grid pos
            (Comp tumLoMHC tumHiMHC fcA fcU fc1 fc2 pcA pcU pc1 pc2 fnA fnU fn1 fn2 pnA pnU pn1 pn2) = origComp
            !wCellCountRet = max wCellCount (countEff origComp) --- for statistics to print at end of program
            --- to divide effector cells up according to tumor cells
            tumNum = tumLoMHC + tumHiMHC
            effNum = fcA + fcU + fc1 + fc2 + pcA + pcU + pc1 + pc2 + fnA + fnU + fn1 + fn2 + pnA + pnU + pn1 + pn2
        in if (tumNum == 0) || (effNum == 0)
           then inVivoLysing grid distSeq (pos+1) lysCount wCellRemCount wCellCountRet gen
           else let (fcAlst, genfcUlst) = randomRsGen (0::Int, tumNum-1::Int) fcA [] gen
                    (fcUlst, genfc1lst) = randomRsGen (0::Int, tumNum-1::Int) fcU [] genfcUlst
                    (fc1lst, genfc2lst) = randomRsGen (0::Int, tumNum-1::Int) fc1 [] genfc1lst
                    (fc2lst, genpcAlst) = randomRsGen (0::Int, tumNum-1::Int) fc2 [] genfc2lst
                    (pcAlst, genpcUlst) = randomRsGen (0::Int, tumNum-1::Int) pcA [] genpcAlst
                    (pcUlst, genpc1lst) = randomRsGen (0::Int, tumNum-1::Int) pcU [] genpcUlst
                    (pc1lst, genpc2lst) = randomRsGen (0::Int, tumNum-1::Int) pc1 [] genpc1lst
                    (pc2lst, genfnAlst) = randomRsGen (0::Int, tumNum-1::Int) pc2 [] genpc2lst
                    (fnAlst, genfnUlst) = randomRsGen (0::Int, tumNum-1::Int) fnA [] genfnAlst
                    (fnUlst, genfn1lst) = randomRsGen (0::Int, tumNum-1::Int) fnU [] genfnUlst
                    (fn1lst, genfn2lst) = randomRsGen (0::Int, tumNum-1::Int) fn1 [] genfn1lst
                    (fn2lst, genpnAlst) = randomRsGen (0::Int, tumNum-1::Int) fn2 [] genfn2lst
                    (pnAlst, genpnUlst) = randomRsGen (0::Int, tumNum-1::Int) pnA [] genpnAlst
                    (pnUlst, genpn1lst) = randomRsGen (0::Int, tumNum-1::Int) pnU [] genpnUlst
                    (pn1lst, genpn2lst) = randomRsGen (0::Int, tumNum-1::Int) pn1 [] genpn1lst
                    (pn2lst, gen2rec) = randomRsGen (0::Int, tumNum-1::Int) pn2 [] genpn2lst
                    interactRec :: Composition -> Int -> Int -> StdGen -> (Composition, Int, StdGen)
                    interactRec compAcc effRemAcc num genInt = 
                        if num == tumNum   --- num starts at 0 and is incremented each call
                        then (compAcc, effRemAcc, genInt)
                        else let tumLoMHCRec = if num<tumLoMHC
                                               then 1
                                               else 0
                                 tumHiMHCRec = if num>=tumLoMHC
                                               then 1
                                               else 0
                                 fcArec = Prelude.length (Prelude.filter (== num) fcAlst)
                                 fcUrec = Prelude.length (Prelude.filter (== num) fcUlst)
                                 fc1rec = Prelude.length (Prelude.filter (== num) fc1lst)
                                 fc2rec = Prelude.length (Prelude.filter (== num) fc2lst)
                                 pcArec = Prelude.length (Prelude.filter (== num) pcAlst)
                                 pcUrec = Prelude.length (Prelude.filter (== num) pcUlst)
                                 pc1rec = Prelude.length (Prelude.filter (== num) pc1lst)
                                 pc2rec = Prelude.length (Prelude.filter (== num) pc2lst)
                                 fnArec = Prelude.length (Prelude.filter (== num) fnAlst)
                                 fnUrec = Prelude.length (Prelude.filter (== num) fnUlst)
                                 fn1rec = Prelude.length (Prelude.filter (== num) fn1lst)
                                 fn2rec = Prelude.length (Prelude.filter (== num) fn2lst)
                                 pnArec = Prelude.length (Prelude.filter (== num) pnAlst)
                                 pnUrec = Prelude.length (Prelude.filter (== num) pnUlst)
                                 pn1rec = Prelude.length (Prelude.filter (== num) pn1lst)
                                 pn2rec = Prelude.length (Prelude.filter (== num) pn2lst)
                                 compRec = (Comp tumLoMHCRec tumHiMHCRec fcArec fcUrec fc1rec fc2rec pcArec pcUrec pc1rec pc2rec fnArec fnUrec fn1rec fn2rec pnArec pnUrec pn1rec pn2rec)
                                 (compRecNew, effRemNew, genIntRet) = eachInteract compRec genInt
                                 compAccRet = sumComp compAcc compRecNew
                                 effRemRet = effRemAcc + effRemNew
                              in interactRec compAccRet effRemRet (num+1) genIntRet
                    compA = (Comp 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0) 
                    (newComp, effRem, genRet) = interactRec compA 0 0 gen2rec      
                    !gridRet = update pos newComp grid
                    tumLys = tumNum - (countTum newComp)
                    lysCountRet = lysCount + tumLys
                    wCellRemCountRet = wCellRemCount + effRem
                    !distSeqRet = if ((countTum newComp) == 0) && ((countTum origComp) /= 0)
                                  then delTumDistSeq distSeq pos
                                  else distSeq
                in inVivoLysing gridRet distSeqRet (pos+1) lysCountRet wCellRemCountRet wCellCountRet genRet


newAvailableLst :: Seq Composition -> Seq Int -> [Int]
--- Inputs:
---   grid :: Seq Composition - the grid before lysis
---   distSeq :: Seq Int - distance sequence
--- Outputs:
---   :: [Int] - a new list of available tumor cells 
newAvailableLst grid distSeq = nLstRec 0 [] where
    nLstRec pos lst = 
        if pos == gridSize --- pos starts at 0 and is incremented each call
        then lst
        else let dist = index distSeq pos
              in if dist == 1 || (dist == 0 && (countTum (index grid pos)) < maxTumNumVivo) --- either one away from tumor cell or in grid point that has a tumor cell but that is not at maximum capacity
                 then nLstRec (pos+1) (pos:lst)
                 else nLstRec (pos+1) lst 


--------------------------------------------
--- procedures for adding cells to the grid
--------------------------------------------

--- This procedure distributes effector cells uniformly over grid points not occupied by a tumor cell
distributeEffAvoidTum :: Seq Composition -> Seq Int -> Int -> (Composition -> Composition) -> StdGen -> (Seq Composition, StdGen)
distributeEffAvoidTum grid distSeq num fn gen = distributeRec grid num gen where
    (availablePos, availableNum) = available 0 0 []
    available :: Int -> Int -> [Int] -> ([Int],Int)
    available pos num lst = 
        if pos == gridSize
        then (lst, num)
        else if index distSeq pos > 0
             then available (pos+1) (num+1) (pos:lst)
             else available (pos+1) num lst    
    distributeRec gridRec numRec genRec = 
        if numRec == 0
        then (gridRec, genRec)
        else let (newPos, newGen) = if availableNum == 0
                                    then randomR (0 :: Int, gridSize-1::Int) genRec
                                    else let (indx, nGen) = randomR (0 :: Int, availableNum-1::Int) genRec
                                             nPos = availablePos!!indx
                                          in (nPos, nGen)
                 cellLst = index gridRec newPos  --- the composition at position newPos
                 !newGrid = update newPos (fn cellLst) gridRec --- Increment the designated component of the composition at location newPos
              in distributeRec newGrid (numRec-1) newGen


--- Input
---   grid :: Seq Composition - a position sequence of cell compositions
---   distSeq :: Seq Int - a list of distances to tumor cells
---   num :: Int - the number of cells to add (must be non-negative)
---   fn :: (Composition -> Composition) - a function that increments an effector cell component of the composition
---   gen :: StdGen - a random generator
distributeEffTumEdge :: Seq Composition -> Seq Int -> Int -> (Composition -> Composition) -> StdGen -> (Seq Composition, StdGen)
distributeEffTumEdge grid distSeq num fn gen = distributeRec grid num gen where
    (availablePos, availableNum) = available 0 0 []
    available :: Int -> Int -> [Int] -> ([Int],Int)
    available pos num lst = 
        if pos == gridSize
        then (lst, num)
        else if index distSeq pos == 1
             then available (pos+1) (num+1) (pos:lst)
             else available (pos+1) num lst    
    distributeRec gridRec numRec genRec = 
        if numRec == 0
        then (gridRec, genRec)
        else let (newPos, newGen) = if availableNum == 0
                                    then randomR (0 :: Int, gridSize-1::Int) gen
                                    else let (indx, nGen) = randomR (0 :: Int, availableNum-1::Int) gen
                                             nPos = availablePos!!indx
                                          in (nPos, nGen)
                 cellLst = index gridRec newPos  --- the composition at position newPos
                 !newGrid = update newPos (fn cellLst) gridRec --- Increment the designated component of the composition at location newPos
              in distributeRec newGrid (numRec-1) newGen

distributeEffEdge :: Seq Composition -> Int -> (Composition -> Composition) -> StdGen -> (Seq Composition, StdGen)
--- Input
---   grid :: Seq Composition - a position sequence of cell compositions
---   num :: Int - the number of cells to add (must be non-negative)
---   fn :: (Composition -> Composition) - a function that increments an effector cell component of the composition
---   gen :: StdGen - a random generator
--- Global Constants
---   gridDim :: Num a => a - square root of the total number of grid elements
--- Output
---   :: Seq Composition - an updated sequence of compositions with the numbers of a given effector cell added to the grid 
---   :: StdGen - an updated random generator
distributeEffEdge grid num fn gen = 
    if num == 0
    then (grid, gen)
    else let edgeSize = 4*gridDim -4
             (indx, newGen) = randomR (0 :: Int, edgeSize-1 :: Int) gen
             newPos = if indx <= gridDim
                      then indx  --- First gridDim+1 indices coincide with the position
                      else if indx <= 3*gridDim-4
                           then let row = 1 + div (indx-gridDim) 2 --- for example, if the index is gridDim, then the row is 1 (below the top row) and column is 0
                                    column = if rem (indx-gridDim) 2 == 0
                                             then 0
                                             else gridDim-1
                                in gridDim * row + column
                           else gridSize - (edgeSize-indx)
             cellLst = index grid newPos  --- the composition at position newPos
             !newGrid = update newPos (fn cellLst) grid --- Increment the designated component of the composition at location newPos
   in distributeEffEdge newGrid (num - 1) fn newGen
   
addTumsContiguous :: Int -> Seq Composition -> Seq Int -> [Int] -> Int -> StdGen -> (Seq Composition, Seq Int, StdGen)
--- Inputs
---   num :: Int - the number of (additional) tumor cells to add to the grid
---   grid :: Seq Composition - the grid
---   distSeq :: Seq Int - a list of distances to tumor cells
---   availablePos :: [Int] - a list of available positions for tumor cells
---   availableNum :: Int - this should be the length of availablePos
---   gen :: StdGen - a random generator
--- Outputs
---   :: Seq Composition - the updated grid
---   :: Seq Int - updated distance sequence
---   :: StdGen - an updated random generator
addTumsContiguous num grid distSeq availablePos availableNum gen = 
    if num == 0 || (availableNum == 0)
    then (grid, distSeq, gen) --- output after last call of procedure
    else let (idx, newGen) = randomR (0 :: Int, availableNum-1::Int) gen
             newPos = availablePos!!idx
             (randn, newerGen) = randomR (0::Double, 1::Double) newGen  --- determine whether first tumor cell has high MHC
             !newGrid = if randn > chanceLoMHCvivo
                        then update newPos (adjustTumHiMHC (\x->x+1) (index grid newPos)) grid
                        else update newPos (adjustTumLoMHC (\x->x+1) (index grid newPos)) grid
             !nDistSeq = if 0 == countTum (index grid newPos) --- if the position originally didn't have any tumor cells
                         then addTumDistSeq distSeq newPos
                         else distSeq
             (adjAvailablePos,adjAvailableNum) = 
                             if (countTum $ index newGrid newPos) == maxTumNumVivo
                             then ([x | x <- availablePos, x /= newPos],availableNum-1)
                             else (availablePos,availableNum)
             nbrs = nbrPos newPos
             toAdd = [x | x<- nbrs \\ adjAvailablePos, (countTum $ index newGrid x) < maxTumNumVivo] 
             !nAvailableNum = adjAvailableNum + (Prelude.length toAdd)
             !nAvailablePos = toAdd++adjAvailablePos
         in addTumsContiguous (num-1) newGrid nDistSeq nAvailablePos nAvailableNum newerGen

-------------------------------------------
--- Procedures for updating position sequence
-------------------------------------------

addTumDistSeq :: Seq Int -> Int -> Seq Int
--- Input
---   seq : Seq Int - an input position sequence that has not been updated since the last complete update
---   pos :: Int - position where tumor is added  
--- Output
---   :: Seq Int - an updated distance sequence of distances from tumor cells  
addTumDistSeq seq pos = addTumDistSeqRec nSeq 1 (nbrPos pos) where
    !nSeq = update pos 0 seq
    addTumDistSeqRec :: Seq Int -> Int -> [Int] -> Seq Int
    addTumDistSeqRec seqRec dist list 
        | (Prelude.null list) && (dist == 3) = seqRec --- output of whole procedure
        | (Prelude.null list) && (dist == 2) = addTumDistSeqRec seqRec 3 (nbrPos3away pos)
        | (Prelude.null list) && (dist == 1) = addTumDistSeqRec seqRec 2 (nbrPos2away pos)
        | otherwise = let here = head list
                          val = index seqRec here
                          newVal | val == (distDepth+1) = dist
                                 | otherwise = (min val dist)
                          !nSeqRec = update here newVal seqRec
                       in addTumDistSeqRec nSeqRec dist (tail list)

delTumDistSeq :: Seq Int -> Int -> Seq Int             
--- Input
---   distSeq :: Seq Int - an input position sequence that has not been updated since the last complete update    
---   pos :: Int - position where tumor is deleted
--- Global constants
---   distDepth :: Integral a => a - (fixed to 3) the maximum number of grid elements a CTL can detect the proximity of a tumor cell
--- Output
---   :: Seq Int - an updated distance sequence of distances from tumor cells                  
delTumDistSeq distSeq pos = delTumDistSeqClear nSeq allPos where
    !nSeq = update pos (distDepth+1) distSeq --- (distDepth + 1) is equal to 4
    allPos = pos:(nbrPos pos)++(nbrPos2away pos)++(nbrPos3away pos) --- positions to be updated
    delTumDistSeqClear seqRec lst
      | (Prelude.null lst) = delTumDistSeqUpdate seqRec allPos distDepth --- distDepth is equal to 3
      | otherwise = 
          let here = (head lst)
              val = index seqRec here
              newVal | val==0 = 0
                     | otherwise = (distDepth +1)
              !nSeqRec = update here newVal seqRec
           in delTumDistSeqClear nSeqRec (tail lst)


delTumDistSeqUpdate :: Seq Int -> [Int] -> Int -> Seq Int                          
--- Input
---   distSeq :: Seq Int - input partially updated distance sequence
---   posLst :: [Int] - list of positions that have yet to be updated
---   num :: Int - starts off at 3 (distDepth) and iterates down each time
--- Output
---   :: Seq Int - an updated distance sequence of distances from tumor cell
delTumDistSeqUpdate distSeq posLst num --- num starts off at 3 (distDepth) and iterates down each time
    | (num==0) = distSeq  --- distance sequence output by procedure after all calls
    | otherwise = delTumDistSeqUpdate nDistSeq posLst (num-1) where
        !nDistSeq = updatePositions distSeq posLst
        updatePositions :: Seq Int -> [Int] -> Seq Int
        updatePositions distSeqRec posLstRec   --- update all positions in distSeqRec with respect to distSeq
            | (Prelude.null posLstRec) = distSeqRec --- output of updatePositions
            | otherwise = let here = head posLstRec
                              newVal = createNewVal here
                              !nDistSeq = update here newVal distSeqRec --- the only place where the data sequence changes
                              in updatePositions nDistSeq (tail posLstRec)
        createNewVal :: Int -> Int
        createNewVal here =             --- create new value for "here" based on "distSeq"
           let nearVals = Data.List.filter (<=distDepth) (map (\x -> index distSeq x) (nbrPos here))
            in if Prelude.null nearVals --- if neighbors all have value 4
               then index distSeq here --- then keep value of here the same
               else minimum [index distSeq here, (minimum nearVals) + 1] 
               
               
writeDistData :: Seq Int  -> Int -> Handle -> IO()
--- Inputs:
---   posSeq :: Seq Composition - a sequence of distances to nearest tumor cell
---   pos :: Int - position to be writtin at this stage
---   handle :: Handle - a handle for the file to write data in.
--- Effects:
---   Write a line of distance data
writeDistData seq pos handle = do
  if pos < gridSize
  then do
    let val = index seq pos 
    hPutStrLn handle $ (show val) 
    writeDistData seq (pos+1) handle
  else return ()
