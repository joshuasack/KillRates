--- ************************
--- **** InVitro module ****
--- ************************
--- -----------------------------------------------------------
--- **** Brief description of the main procedure inVitro: ****
--- -----------------------------------------------------------
---
---  > A number of types of experiments are run (how many is determined by parameter dataSize), each with a different tumor to effector cell ratio (by fixing the number of tumor cells and adjusting the number of effector cells)
---  > for each experiment type, there are reps many assays 
---  > the number of time steps of each assay is set to the parameter timeStepsInVitro
---   
--- effects
---  > write position data for every stage (beginning and end of each time step) of the the last run of each experiment type
---  > write, for each experiment type, the percent lysis of each run, the average percent lysis, and the variance of the percent lysis
---  > print progress throughout the program
---  > print the maximum number of effector cells that have been in the same grid element as a tumor cell after a movement phase.
---
--- ------------------------------
--- **** Position Data Format ****
--- ------------------------------
--- 
--- Each position data file given the position data for the last run of an experiment type, and consists of an array with:
--- 18 columns, corresponding to the components of the composition data type, i.e. from left to right:
---    > low MHC expressing tumor, high MHC tumor, 
---    > fasCTLa (just arrived in the grid element or is not colocated with a tumor cell), 
---    > FasCTLu (previously colocated with a tumor cell, but still not recognizing it),
---    > fasCTL1 (recognized tumor cell 1 time step ago), 
---    > fasCTL2 (recognized tumor cell 2 time steps ago), 
---    > perCTLa, perCTLu, perCTL1, perCTL2, 
---    > fasNKLa, fasNKLu, fasNKL1, fasNKL2, 
---    > perNKLa, perNKLu, perNKL1, and perNKL2
--- gridSize many rows (corresponding to the positions: the first gridDim correspond to the first row)
---
--- The value of each entry is the number of cells (whose type is determined by the column) in the appropriate position (determined by the row)
---
--- --------------------------------
--- **** statistics file format ****
--- --------------------------------
---
--- For each experiment type, there are the following lines:
---    > a line explaining the experiment type (the exponent of a given base)
---    > a line giving the average percent lysis
---    > a line giving the variance of percent lysis
---    > a line, for each run, giving the percent lysis of that run
--- dataSize many such lists, each distinguished by their first line (giving exponent)
---
--- ***********************************

module InVitro where

import System.Random
import Data.List
import Data.Sequence
import System.IO
import System.Directory
import Biofunctions



--- ---------------------
--- Constant parameters
--- ---------------------

dataSize = 1  --- :: Num a => a - number of effector/tumor ratios to plot 
iExponent = 1 --- :: Num a => a - floor(iTumNumInVitro*base^^iExponent) will be the initial number of effector cells.

--- Used throughout the program
iTumNumInVitro = 10000 --- :: Num a => a --- initial number of tumor cells (in vitro: to match Seki, this should be about 10,000 divided by 96 wells which is approximately 100)
chanceLoMHCvitro = 0.2 --- :: Num a => a --- the chance a tumor cell has low MHC

--- Used inVitroExpTypes
fractionFC = 0.25 --- :: Fractional a => a - The approximate fraction of effector cells that are fasL expressing CTLs
fractionPC = 0.75  --- :: Fractional a => a - The approximate fraction of effector cells that are perforin expressing CTLs
fractionFN = 0 --- :: Fractional a => a - The approximate fraction of effector cells that are fasL expressing NKLs
fractionPN = 0 --- :: Fractional a => a - The approximate fraction of effector cells that are perforin expressing NKLs
base = 20 --- :: Num a => a - The factor by which we increase each number of CTLs or NKLs for each data point.
reps = 5 --- :: Num a => a - The number of assays run with the same setting

--- Used in inVitroTimeSteps
timeStepsInVitro = 8 --- :: Num a => a - the number of timesteps for each assay (3 for "8 hour assay" and 6 for "18 hour assay")
speedBoundInVitro = 4 --- :: Num a => a - the upper limit of the number of movement iterations in a time step

--- in positionLst and availablePos
maxTumNumVitro = 1  --- :: Num a => a -  the maximum number of tumor cells that can be in a grid element

--- in distributeTums
aveClumpSize = 75 --- :: Num a => a - the expected number of tumor cells in each cluster

--- ------------------
--- main procedure
--- ------------------

inVitro :: IO () 
--- Global contstants
---   dataSize :: Num a => a - number of effector/tumor ratios to plot 
---   iExponent :: Num a => a - floor(base^iExponent*iTumNumInVitro) will be the initial number of effector cells.
--- Effects
---   Write position data (a file for each exponent and each point in time) and lysis statistics
---   Print progress
---   Print he maximum number of leucocytes in the same grid element as a tumor cells after a movement phase 
inVitro = do 
--  gen <- getStdGen  --- If we want "real" randomization
    let gen = (mkStdGen 100) --- If we want to fix a "seed"
    createDirectoryIfMissing True "data"  -- True could be changed to False.  It is only relevant if we specify a parent directory
    statsDataHandle <- openFile "data/inVitroStats.dat" WriteMode
    tableDataHandle <- openFile "data/inVitroTableData.dat" WriteMode
    wCellCount <- inVitroExpTypes iExponent dataSize 0 statsDataHandle tableDataHandle gen
    putStrLn $ "The maximum number of leucocytes in the same grid element as a tumor cells after a movement phase is " ++ (show wCellCount) ++ "."
    hClose statsDataHandle
    hClose tableDataHandle

inVitroExpTypes :: Int -> Int -> Int -> Handle -> Handle -> StdGen -> IO Int
--- Input
---   exponent :: Int - iTumNumInVitro * base^^exponent is the number of effector cells (modulo some rounding down)
---   left :: Int - the number of ratios left to plot
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   statsDataHandle :: Handle - a handle for the file to write percent lysis data in
---   tableDataHandle :: Handle - a handle for the file to write table of output for plotting purposes
---   gen :: StdGen - a random regerator
--- Global constants used
---   fractionFC :: Fractional a => a - The approximate fraction of effector cells that are fasL expressing CTLs
---   fractionPC :: Fractional a => a - The approximate fraction of effector cells that are perforin expressing CTLs
---   fractionFN :: Fractional a => a - The approximate fraction of effector cells that are fasL expressing NKLs
---   fractionPN :: Fractional a => a - The approximate fraction of effector cells that are perforin expressing NKLs
---   base :: Num a => a - The factor by which we increase each number of CTLs or NKLs for each data point.
---   reps :: Num a => a -  number of assays run with the same setting
--- Output:
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
--- Effects:
---   Write position data (a file for each exponent and point of time) and lysis statistics
---   Print progress
inVitroExpTypes exponent left wCellCount statsDataHandle tableDataHandle gen = 
    if left == 0
    then return wCellCount
    else do
      let iFCNum = floor $ fractionFC*((fromIntegral base)^^exponent * (fromIntegral iTumNumInVitro)) 
          iPCNum = floor $ fractionPC*((fromIntegral base)^^exponent * (fromIntegral iTumNumInVitro)) 
          iFNNum = floor $ fractionFN*((fromIntegral base)^^exponent * (fromIntegral iTumNumInVitro))
          iPNNum = floor $ fractionPN*((fromIntegral base)^^exponent * (fromIntegral iTumNumInVitro))
      (lysListList, remListList, speedList, nwCellCount, newGen) <- assayReps iFCNum iPCNum iFNNum iPNNum exponent reps [] [] [] wCellCount gen --- run assays
      hPutStrLn statsDataHandle $ "Exponent of base " ++ (show base) ++ " for E:T ratio: " ++ (show exponent)
      writeData statsDataHandle tableDataHandle lysListList remListList speedList
      inVitroExpTypes (exponent+1) (left-1) nwCellCount statsDataHandle tableDataHandle newGen

writeData :: Handle -> Handle -> [[Int]] -> [[Int]] -> [Int] -> IO ()
--- Input
---   handle :: Handle - a handle for the lysis statistics file
---   tablehandle :: Handle - a handle for the table of lysis file
---   lysLstLst :: [[Int]] - a list of lists of lysis counts (there are rep-many lists; each list in the list gives the lysis counts of a given assay with the last time step at the head of the list and the first time step at the tail of the list)
---   remLstLst :: [[Int]] - a list of lists of effector cell removal counts
---   speedLst :: [Int] - a list of speeds, with the speed of the latest run at the head.
writeData handle tablehandle lysLstLst remLstLst speedLst = 
    if (not (Data.List.null $ head lysLstLst))
    then do
        writeData handle tablehandle (map tail lysLstLst) (map tail remLstLst) speedLst
        let lysList = (map head lysLstLst)
            remList = (map head remLstLst)
            aveLys = (fromIntegral $ sum lysList)/(fromIntegral $ Prelude.length lysList)
            aveRem = (fromIntegral $ sum remList)/(fromIntegral $ Prelude.length remList)
            apl = (fromIntegral 100)*aveLys/(fromIntegral iTumNumInVitro)
            var = let lnth = (fromIntegral $ Prelude.length lysList)  --- sample variance of percent lysis
                   in if lnth == 1
                      then 0
                      else (sum $ map (\x -> (100*(fromIntegral x)/(fromIntegral iTumNumInVitro) - apl)^2) lysList)/(lnth-1)
        hPutStrLn handle $ " Time-step: " ++ (show (Prelude.length (head lysLstLst)))
        hPutStrLn handle $ "   Average percent lysis: " ++ (show apl)
        hPutStrLn handle $ "   Variance of percent lysis: " ++ (show var) 
        hPutStrLn handle $ "   Average number of effector cells becoming ineffective since beginning of experiment: " ++ (show aveRem)
        writePL handle lysList speedLst
        hPutStrLn tablehandle $ (show (Prelude.length (head lysLstLst))) ++  "," ++ (show apl) ++ ","
    else return ()

writePL :: Handle -> [Int] -> [Int] -> IO () --- Write "percent lysis"
--- Input
---   handle :: Handle - a handle for the lysis statistics file 
---   lyslist :: [Int] - a list of actual tumor lysis, one entry for each run
---   speedlist :: [Int] - a list of speeds
--- Global constants
---   reps :: Num a => a - The number of assays run with the same setting
--- Effects
---   Write in the stats file reps-many lines, each giving the percent lysis far the corresponding run (or assay)
writePL handle lyslist speedlist = 
    if (Prelude.length lyslist) > 0
    then do hPutStrLn handle $ "     Run " ++ (show (reps - (Prelude.length lyslist)+1)) ++ " with speed " ++ (show (head speedlist)) ++ " has percent lysis: " ++ (show $ 100*(fromIntegral $ head lyslist)/(Prelude.fromIntegral iTumNumInVitro))
            writePL handle (tail lyslist) (tail speedlist)
    else return ()
                            

assayReps ::  Int -> Int -> Int -> Int -> Int -> Int -> [[Int]] -> [[Int]] -> [Int] -> Int -> StdGen -> IO ([[Int]], [[Int]], [Int], Int, StdGen)
--- Inputs
---   iFCNum :: Int - initial number of Fas CTL's 
---   iPCNum :: Int - initial number of perforin CTL's 
---   iFNNum :: Int - initial number of Fas NKL's 
---   iPNNum :: Int - initial number of perforin NKL's
---   exponent :: Int - there will be base^exponent times as many effector cells as tumor cells
---   repsleft :: Int - number of assays of the same type to perform
---   lysLstLst :: [[Int]] - a list of lists of cumulative lysis counts, with the list at the head being the most recent run
---   remLstLst :: [[Int]] - a list of lists of cumulative effector cell removal counts, with the list at the head being the most recent run
---   speedLst :: [Int] - a list of speeds, with the speed at the head being the most recent run
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   gen :: StdGen - a random generator
--- Global constants used
---   gridSize :: Num a => a - the number of grid elements (used for making the position arrays)
--- Outputs
---   :: [[Int]] - a list of lists of lysis counts (each list in the list gives the lysis counts of a given assay with the last time step at the head of the list and the first time step at the tail of the list)
---   :: [[Int]] - a list of lists of effector cell removal counts
---   :: [Int] - a list of speeds
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   :: StdGen - an updated random generator
--- Effects:
---   Write position data (a file for each exponent and point of time)
---   Print progress
assayReps iFCNum iPCNum iFNNum iPNNum exponent repsleft lysLstLst remLstLst speedLst wCellCount gen = 
    if repsleft == 0
    then do return (lysLstLst,remLstLst,speedLst,wCellCount,gen)
    else do
      let emptySeq = Data.Sequence.replicate gridSize (Comp 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0)   --- create a data sequence of eighteen-tuples
          !(withTumSeq, genFC) = distributeTums iTumNumInVitro emptySeq (maxTumNumVitro * gridSize) gen
          posLst = [x | x <- [0..gridSize-1], (countTum $ index withTumSeq x) > 0]
      putStrLn $ "Finished filling in tumor cells."
      let !(withFCseq, genPC) = distributeEff withTumSeq iFCNum (adjustFCa (\x -> x+1)) genFC
      putStrLn $ "Finished filling in Fas CTLs."
      let !(withPCseq, genFN) = distributeEff withFCseq iPCNum (adjustPCa (\x -> x+1)) genPC
      putStrLn $ "Finished filling in perforin CTLs."
      let !(withFNseq, genPN) = distributeEff withPCseq iFNNum (adjustFNa (\x -> x+1)) genFN
      putStrLn $ "Finished filling in Fas NKLs."
      let !(withAllseq, genSpeed) = distributeEff withFNseq iPNNum (adjustPNa (\x -> x+1)) genPN
      putStrLn $ "Finished filling in perforin NKLs."
      let (speed, genTS) = randomR (1 :: Int, speedBoundInVitro :: Int) genSpeed
      if repsleft == 1
      then do
        handle <- openFile ("data/inVitroPosExp" ++ (show exponent) ++ "Time" ++ (show 0) ++ ".dat") WriteMode
        writePosData withAllseq 0 handle
        putStrLn $ "Finished writing position data for assay with exponent " ++ (show exponent) ++ " at time " ++ (show 0) ++ "."
        hClose handle
        (lysisLst, wCellRemLst, wCellCountNew, genNxt) <- inVitroTimeStep True withAllseq exponent (nub posLst) 0 [] [] wCellCount speed genTS 
        putStrLn $ "With speed " ++ (show speed) ++ ", the final lysis is " ++ (show (head lysisLst)) ++ "."
        (lysListRet,remListRet,speedListRet,wCellCountRet, genRet) <- assayReps iFCNum iPCNum iFNNum iPNNum exponent (repsleft - 1) (lysisLst:lysLstLst) (wCellRemLst:remLstLst) (speed:speedLst) wCellCountNew genNxt
        return (lysListRet,remListRet, speedListRet, wCellCountRet, genRet)
      else do 
        (lysisLst, wCellRemLst, wCellCountNew, genNxt) <- inVitroTimeStep False withAllseq exponent (nub posLst) 0 [] [] wCellCount speed genTS
        putStrLn $ "With speed " ++ (show speed) ++ ", the final lysis is " ++ (show (head lysisLst)) ++ "."
        (lysListRet,remListRet, speedListRet, wCellCountRet, genRet) <- assayReps iFCNum iPCNum iFNNum iPNNum exponent (repsleft - 1) (lysisLst:lysLstLst) (wCellRemLst:remLstLst) (speed:speedLst) wCellCountNew genNxt
        return (lysListRet, remListRet, speedListRet, wCellCountRet, genRet)


inVitroTimeStep :: Bool ->  Seq Composition -> Int -> [Int] -> Int -> [Int] -> [Int] -> Int -> Int -> StdGen -> IO ([Int], [Int], Int, StdGen)
--- Input 
---   write :: Bool - if true, then return a real position list (this is the assay of a given type), otherwise do not modify position list (leave it empty)
---   seq :: Seq Composition - a position sequence of cell compositions
---   exponent :: Int - the exponent for the assay used (for creating a file name and printing status)
---   tumPosLst :: [Int] - a list of positions where there is at least one tumor cell (so position is repeated)
---   time :: Int - the number of timesteps already completed (this value is 0 if this is the first time step)
---   lysCountLst :: [Int] - list of tallies up how many tumor cells have been lysed
---   wCellRemLst :: [Int] - list of numbers of effecter cells that become ineffective
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   speed :: Int - the number of movement phases in each time step
---   gen :: StdGen - random generator
--- Global constants used:
---   timeStepsInVitro :: Num a => a - the number of timesteps for each assay (3 for "8 hour assay" and 6 for "18 hour assay")
--- Output
---   :: [Int] - list of tumor cell lysis counts
---   :: [Int] - list of numbers of effecter cells that become ineffective
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   :: StdGen - an updated random generator
--- Effects (only if "write" is true):
---   Write position data (for the given exponent and time step)
---   Print progress on writing position data 
inVitroTimeStep write seq exponent tumPosLst time lysCountLst wCellRemLst wCellCount speed gen = 
    if time == timeStepsInVitro
    then return (lysCountLst, wCellRemLst, wCellCount, gen)
    else do
      let !(movedSeq, genLys) = moveCellsInVitro seq speed gen
      putStr $ "Finished the " ++ (show speed) ++ " movement phases for this time step.  "
      let !(nSeq, ntumPosLst, nLysCount, nwCellRemCount, nwCellCount, nGen) = inVitroLysing movedSeq tumPosLst (if (Data.List.null lysCountLst) then 0 else (head lysCountLst)) 0 wCellCount genLys
      putStrLn "Finished a lysis phase.  "
      if write
      then do
        handle <- openFile ("data/inVitroPosExp" ++ (show exponent) ++ "Time" ++ (show $ time+1) ++ ".dat") WriteMode
        writePosData nSeq 0 handle
        putStrLn $ "Finished writing position data for assay with exponent " ++ (show exponent) ++ " at time " ++ (show $ time+1) ++ "."
        hClose handle
        (lysCountRet, wCellRemCountRet, wCellCountRet, genRet) <- inVitroTimeStep write nSeq exponent ntumPosLst (time+1) (nLysCount:lysCountLst) (nwCellRemCount:wCellRemLst) nwCellCount speed nGen
        return (lysCountRet, wCellRemCountRet, wCellCountRet, genRet)
      else do
        (lysCountRet, wCellRemCountRet, wCellCountRet, genRet) <- inVitroTimeStep write nSeq exponent ntumPosLst (time+1) (nLysCount:lysCountLst) (nwCellRemCount:wCellRemLst) nwCellCount speed nGen
        return (lysCountRet, wCellRemCountRet, wCellCountRet, genRet)
          
moveCellsInVitro :: Seq Composition -> Int -> StdGen -> (Seq Composition, StdGen)
--- Inputs:
---   seqOrig :: Seq Composition - the position sequence of cell compositions from the beginning of the lysis stage of a time step
---   times :: Int - the number of movement iterations left
---   gen :: StdGen - a random generator
--- Global constants:
---   gridSize :: Num a => a - The number of positions on the grid
--- Outputs:  
---   :: Seq Composition - an updated position sequence of cell compositions
---   :: StdGen - an updated random generator
moveCellsInVitro seqOrig times gen = if times == 0 then (seqOrig, gen) else moveCellsRec seqOrig 0 gen
    where moveCellsRec :: Seq Composition -> Int -> StdGen -> (Seq Composition, StdGen)
          moveCellsRec seqRec pos genRec =
              if pos == gridSize
              then moveCellsInVitro seqRec (times-1) genRec 
              else let compOrig = index seqOrig pos 
                       (Comp tumLoMHCorig tumHiMHCorig fcAorig fcUorig fc1orig fc2orig pcAorig pcUorig pc1orig pc2orig fnAorig fnUorig fn1orig fn2orig pnAorig pnUorig pn1orig pn2orig) = compOrig --- defined here to be used in prep for subrecursive routines
                    in if ((countTum compOrig) > 0) || ((countEff compOrig) == 0)
                       then moveCellsRec seqRec (pos+1) genRec 
                           else let compCur = index seqRec pos 
                                    (Comp tumLoMHCcur tumHiMHCcur fcAcur fcUcur fc1cur fc2cur pcAcur pcUcur pc1cur pc2cur fnAcur fnUcur fn1cur fn2cur pnAcur pnUcur pn1cur pn2cur) = compCur
                                    compNew = (Comp tumLoMHCcur tumHiMHCcur (fcAcur - fcAorig) (fcUcur - fcUorig) (fc1cur - fc1orig) (fc2cur - fc2orig) (pcAcur - pcAorig) (pcUcur - pcUorig) (pc1cur - pc1orig) (pc2cur - pc2orig) (fnAcur - fnAorig) (fnUcur - fnUorig) (fn1cur - fn1orig) (fn2cur - fn2orig) (pnAcur - pnAorig) (pnUcur - pnUorig) (pn1cur - pn1orig) (pn2cur - pn2orig))
                                    !seqNew = update pos compNew seqRec -- all original white cells move out
                                    nbrs = nbrPos pos
                                    nbrNum = Prelude.length nbrs
                                    (fcLst,genPC) = randomRsGen (0 :: Int, nbrNum -1 :: Int) (countFC compOrig) [] genRec  -- determine where cells go
                                    (pcLst,genFN) = randomRsGen (0 :: Int, nbrNum -1 :: Int) (countPC compOrig) [] genPC
                                    (fnLst,genPN) = randomRsGen (0 :: Int, nbrNum -1 :: Int) (countFN compOrig) [] genFN
                                    (pnLst,genRet) = randomRsGen (0 :: Int, nbrNum -1 :: Int) (countPN compOrig) [] genPN
                                 in moveCellsSubRec (nbrNum-1) nbrs pos (fcLst, pcLst, fnLst, pnLst) seqNew genRet
          moveCellsSubRec :: Int -> [Int] -> Int -> ([Int], [Int], [Int], [Int]) -> Seq Composition -> StdGen -> (Seq Composition, StdGen)
          --- Input:
          ---   indx :: Int - the particular position to add cells to
          ---   nbrs :: [Int] - a list of positions to add cells to
          ---   posHere :: Int - the position whence cells are moving
          ---   fcSubRec :: [Int] - a list of positions to move fas CTL's to
          ---   seqSubRec :: Seq Composition - the grid with which to add cells to
          ---   genSubRec :: StdGen - a random generator
          moveCellsSubRec indx posLst posHere (fcSubRec, pcSubRec, fnSubRec, pnSubRec) seqSubRec genSubRec = 
              if indx >= 0
              then let (fcPart, fcRest) = Data.List.partition (==indx) fcSubRec
                       (pcPart, pcRest) = Data.List.partition (==indx) pcSubRec
                       (fnPart, fnRest) = Data.List.partition (==indx) fnSubRec
                       (pnPart, pnRest) = Data.List.partition (==indx) pnSubRec
                       posThere = posLst!!indx
                       (Comp tumLoMHCthere tumHiMHCthere fcAthere fcUthere fc1there fc2there pcAthere pcUthere pc1there pc2there fnAthere fnUthere fn1there fn2there pnAthere pnUthere pn1there pn2there) = index seqSubRec posThere
                       fcAnew = fcAthere + (Prelude.length fcPart)
                       pcAnew = pcAthere + (Prelude.length pcPart)
                       fnAnew = fnAthere + (Prelude.length fnPart)
                       pnAnew = pnAthere + (Prelude.length pnPart)
                       !seqNew = update posThere (Comp tumLoMHCthere tumHiMHCthere fcAnew fcUthere fc1there fc2there pcAnew pcUthere pc1there pc2there fnAnew fnUthere fn1there fn2there pnAnew pnUthere pn1there pn2there) seqSubRec
                   in moveCellsSubRec (indx-1) posLst posHere (fcRest, pcRest, fnRest, pnRest) seqNew genSubRec
              else moveCellsRec seqSubRec (posHere+1) genSubRec

inVitroLysing :: Seq Composition -> [Int] -> Int -> Int -> Int -> StdGen -> (Seq Composition, [Int], Int, Int, Int, StdGen)
--- Inputs:
---   grid :: Seq Composition - a position sequence of cell compositions
---   tumPosLst :: [Int] - a list of remaining positions to consider for cell interaction
---   lysCount :: Int - a count of the number of tumor cells that have been lysed so far
---   wCellRemCount :: Int - the number of effector cells removed this time step so far
---   wCellCount :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   gen :: StdGen - a random generator
--- Global constants
---   gridSize :: Num a => a - the number of grid elements
---   maxEffNum :: Num a => a - the maximum numbers of cells that can recognize a tumor cell
---   chanceNKLrgLoMHC = 0.7 :: Fractional a => a - the chance an NKL recognizes a tumor cell with low MHC
---   chanceNKLrgHiMHC = 0.3 :: Fractional a => a - the chance an NKL recognizes a tumor cell with high MHC
---   chanceCTLrgLoMHC = 0.1 :: Fractional a => a - the chance a CTL recognizes a tumor cell with low MHC
---   chanceCTLrgHiMHC = 0.8 :: Fractional a => a - the chance a CTL recognizes a tumor cell with high MHC
---   chanceFasNotUsedup :: Fractional a => a - the chance that the cell's supply of FasL is not exhausted
---   chancePerNotUsedup :: Fractional a => a - the chance that the cell's supply of perforin is not exhausted
--- Outputs
---   :: Seq Composition - updated position sequence of cell compositions
---   :: [Int] - updated part of tumor position list
---   :: Int - tumor lysis count
---   :: Int - the number of effector cells removed this time step so far
---   :: Int - the maximum number of leucocytes in the same grid element as a tumor cells after a movement phase so far
---   :: StdGen - updated random generator
inVitroLysing grid tumPosLst lysCount wCellRemCount wCellCount gen = 
    if Data.List.null tumPosLst
    then (grid, tumPosLst, lysCount, wCellRemCount, wCellCount, gen)
    else
        let pos = head tumPosLst
            origComp = index grid pos
            (Comp tumLoMHC tumHiMHC fcA fcU fc1 fc2 pcA pcU pc1 pc2 fnA fnU fn1 fn2 pnA pnU pn1 pn2) = origComp
            !wCellCountRet = max wCellCount (countEff origComp) --- for statistics to print at end of program
            --- to divide effector cells up according to tumor cells 
            tumNum = tumLoMHC + tumHiMHC
            effNum = fcA + fcU + fc1 + fc2 + pcA + pcU + pc1 + pc2 + fnA + fnU + fn1 + fn2 + pnA + pnU + pn1 + pn2           
        in if effNum == 0
           then let (gridFinal, tLstFinal, lysCountFinal, wCellRemCountFinal, wCellCountFinal, genFinal) = inVitroLysing grid (tail tumPosLst) lysCount wCellRemCount wCellCountRet gen
                in (gridFinal, (head tumPosLst):tLstFinal, lysCountFinal, wCellRemCountFinal, wCellCountFinal, genFinal)
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
                    (newComp,effRem,genRet) = interactRec compA 0 0 gen2rec      
                    !gridRet = update pos newComp grid
                    tumLys = tumNum - (countTum newComp)
                    lysCountRet = lysCount + tumLys
                    wCellRemCountRet = wCellRemCount + effRem
                in if (countTum newComp) == 0    --- whether to adjust the tumor position list
                   then inVitroLysing gridRet (tail tumPosLst) lysCountRet wCellRemCountRet wCellCountRet genRet
                   else let (gridFinal, tLstFinal, lysCountFinal, wCellRemCountFinal, wCellCountFinal, genFinal) = inVitroLysing gridRet (tail tumPosLst) lysCountRet wCellRemCountRet wCellCountRet genRet
                        in (gridFinal, (head tumPosLst):tLstFinal, lysCountFinal, wCellRemCountFinal, wCellCountFinal, genFinal)  


------------------------------------------------
--- Procedures for adding tumor cells in clumps
------------------------------------------------

distributeTums :: Int -> Seq Composition -> Int -> StdGen -> (Seq Composition, StdGen)
--- Inputs
---   num :: Int - the number of tumor cells to place
---   grid :: Seq Composition - the grid
---   unfilled :: Int - number of slots (out of maxTumNumVitro*gridSize) that are still unfilled
---   gen :: StdGen - random number generator
--- Global constants
---   maxTumNumVitro :: Num a => a - the maximum number of tumor cells in a grid element
---   gridSize :: Num a => a - the size of the grid
--- Outputs
---   :: Seq Composition - updated grid, with tumors placed in clumps
---   :: StdGen - updated random generator
distributeTums num grid unfilled gen
    | num == 0  = (grid, gen)
    | num == 1  = let (nthPos, newgen) = randomR (0 :: Int, unfilled-1 :: Int) gen
                      newPoint = findInx grid nthPos
                      (randn, newergen) = randomR (0::Double, 1::Double) newgen
                      !newGrid = if randn > chanceLoMHCvitro
                                 then update newPoint (adjustTumHiMHC (\x->x+1) (index grid newPoint)) grid
                                 else update newPoint (adjustTumLoMHC (\x->x+1) (index grid newPoint)) grid
                   in (newGrid, newergen)
    | otherwise = let (almostThisClumpSize, interGen) = randomPoisson aveClumpSize (num-1) gen
                      thisClumpSize = almostThisClumpSize + 1
                      (nthPos, interGen2) = randomR (0 :: Int, unfilled-1 :: Int) interGen
                      thisCenterPoint = findInx grid nthPos
                      (randn, interGen3) = randomR (0::Double, 1::Double) interGen2
                      !interGrid = if randn > chanceLoMHCvitro
                                   then update thisCenterPoint (adjustTumHiMHC (\x->x+1) (index grid thisCenterPoint)) grid
                                   else update thisCenterPoint (adjustTumLoMHC (\x->x+1) (index grid thisCenterPoint)) grid
                      newunfilled = unfilled-1
                      nbrs = nbrPos thisCenterPoint
                      availablePos = [x | x<- thisCenterPoint:nbrs, (countTum $ index interGrid x) < maxTumNumVitro]
                      availableNum = Prelude.length(availablePos)
                      (newgrid, newGen) = addTumsClump (thisClumpSize-1) interGrid availablePos availableNum newunfilled interGen3
                   in distributeTums (num - thisClumpSize) newgrid (unfilled - thisClumpSize) newGen

findInx :: Seq Composition -> Int -> Int
--- Inputs
---   grid :: Seq Composition - the grid
---   nthPos :: Int - the position ranging from 0 to (maxTumNumVitro * gridSize - currentTumNum), where currentTumNum is the current number of tumor cells.  This is not a grid position, but a position of available slots for tumor cells
--- Global constants
---   maxTumNumVitro :: Num a => a - the maximum number of tumor cells in a grid element
--- Output
---   :: Int - the position on the grid that a new tumor cell can be placed
findInx grid nthPos = findInxrec 0 nthPos
    where findInxrec thisPos target =
           let currentTums = (countTum $ index grid thisPos)
               available = maxTumNumVitro - currentTums
            in if target + 1 <= available
               then thisPos
               else findInxrec (thisPos + 1) (target - available)


addTumsClump :: Int -> Seq Composition -> [Int] -> Int -> Int -> StdGen -> (Seq Composition, StdGen)
--- Inputs
---   num :: Int - the number of (additional) tumor cells to add to the grid
---   grid :: Seq Composition - the grid
---   availablePos :: [Int] - a list of available positions for tumor cells
---   availableNum :: Int - this should be the length of availablePos
---   unfilled :: Int - number of slots (out of maxTumNumVitro*gridSize) that are still unfilled
---   gen :: StdGen - a random generator
--- Outputs
---   :: Seq Composition - the updated grid
---   :: Int - the number of tumor cells not successfully entered onto grid (due to lack of space locally)
---   :: [Int] - an updated list of positions of tumor cell already distributed 
---   :: StdGen - an updated random generator
addTumsClump num grid availablePos availableNum unfilled gen = 
    if (num == 0) || (availableNum == 0)
    then distributeTums num grid unfilled gen --- (grid, gen) -- (grid, num, posLst, gen)
    else let (idx, newGen) = randomR (0 :: Int, availableNum-1::Int) gen
             newPos = availablePos!!idx
             (randn, newGen2) = randomR (0::Double, 1::Double) newGen
             !newGrid = if randn > chanceLoMHCvitro
                        then update newPos (adjustTumHiMHC (\x->x+1) (index grid newPos)) grid
                        else update newPos (adjustTumLoMHC (\x->x+1) (index grid newPos)) grid
             newunfilled = unfilled - 1
             (adjAvailablePos,adjAvailableNum) = if (countTum $ index newGrid newPos) == maxTumNumVitro
                             then ([x | x <- availablePos, x /= newPos],availableNum-1)
                             else (availablePos,availableNum)
             nbrs = nbrPos newPos
             toAdd = [x | x<- nbrs \\ adjAvailablePos, (countTum $ index newGrid x) < maxTumNumVitro]
             !nAvailableNum = adjAvailableNum + (Prelude.length toAdd)
             !nAvailablePos = toAdd++adjAvailablePos
          in addTumsClump (num-1) newGrid nAvailablePos nAvailableNum newunfilled newGen2


