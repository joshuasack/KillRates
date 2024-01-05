--- *****************************
--- **** Biofunctions module ****
--- *****************************
--- ------------------------------------------
--- **** Brief description of the module: ****
--- ------------------------------------------
---
--- > This module contains types and procedures used by both the in vivo and in vitro programs:
--- > Defines data type Composition (representing the comoposition of a grid element)
--- > Provides basic procedures randomness (randomBinom, randomPoisson, randomRsGen)
--- > Provides procedure for cell distribution (distributeEff)
--- > Provides supporting procedures for cell interaction, recognition, and lysis (eachInteract, effCellRec, probFasKill, probPerKill)
--- > Provides spatial procedures for movement and adding to grid (nbrPos, nbrPos2away, nbrPos3away)
--- > Introduces basic functions on Compositions (counting certain cells and adjusting composition)
--- > Provides supporting function writePosData
--- > Further explanation and relevant documentation is provided with the code of the procedures

module Biofunctions where

import System.Random
import Data.List
import Data.Sequence
import System.IO

--- -----------------------
--- definition of data type
--- ----------------------- 

type TumLoMHC = Int -- tumor cell with low MHC expression
type TumHiMHC = Int -- tumor cell with high MHC expression
type FCa = Int -- CTL expressing FasL that was not previously co-located with a tumor cell
type FCu = Int -- CTL expressing FasL that has not yet recognized a tumor cell given at least 1 timestep of co-location 
type FCr1 = Int -- CTL expressing Fasl that recognized tumor cell last timestep
type FCr2 = Int -- CTL expressing FasL that recognized tumor cell two timesteps ago
type PCa = Int -- CTL expressing perforin that was not previously co-located with a tumor cell
type PCu = Int -- CTL expressing perforin that has not yet recognized a tumor cell given at least 1 timestep of co-location 
type PCr1 = Int -- CTL expressing perforin that recognized tumor cell last timestep
type PCr2 = Int -- CTL expressing perforin that recognized tumor cell two timesteps ago
type FNa = Int -- NKL expressing FasL that was not previously co-located with a tumor cell
type FNu = Int -- NKL expressing FasL that has not yet recognized a tumor cell given at least 1 timestep of co-location 
type FNr1 = Int -- NKL expressing FasL that recognized tumor cell last timestep
type FNr2 = Int -- NKL expressing FasL that recognized tumor cell two timesteps ago
type PNa = Int -- NKL expressing perforin that was not previously co-located with a tumor cell
type PNu = Int -- KNKL expressing perforin that has not yet recognized a tumor cell given at least 1 timestep of co-location 
type PNr1 = Int -- NKL expressing perforin that recognized tumor cell last timestep
type PNr2 = Int -- NKL expressing perforin that recognized tumor cell two timesteps ago
data Composition = Comp TumLoMHC TumHiMHC FCa FCu FCr1 FCr2 PCa PCu PCr1 PCr2 FNa FNu FNr1 FNr2 PNa PNu PNr1 PNr2 

--- Used throughout the program
gridDim = 300 --- :: Num a => a - cell interaction takes place in a gridDim by gridDim matrix (to match Seki, use 300)
gridSize = gridDim*gridDim 

--- in probFasKill
scale = 5 --- :: Num a => a - the smallest positive number of FasL expressing cells that make the kill rate 0.

--- in lysing
-------- chance of recognition
chanceNKLrgLoMHCa = 0.02 --- :: Fractional a => a - the chance an NKL recognizes a tumor cell with low MHC upon arrival in tumor grid point
chanceNKLrgHiMHCa = 0.01 --- :: Fractional a => a - the chance an NKL recognizes a tumor cell with high MHC upon arrival in tumor grid point
chanceCTLrgLoMHCa = 0.01 --- :: Fractional a => a - the chance a CTL recognizes a tumor cell with low MHC upon arrival in tumor grid point
chanceCTLrgHiMHCa = 0.03 --- :: Fractional a => a - the chance a CTL recognizes a tumor cell with high MHC upon arrival in tumor grid point
chanceNKLrgLoMHCu = 0.5 --- :: Fractional a => a - the chance an NKL recognizes a tumor cell with low MHC having been co-located with tumor cell for at least one turn
chanceNKLrgHiMHCu = 0.1 --- :: Fractional a => a - the chance an NKL recognizes a tumor cell with high MHC having been co-located with tumor cell for at least one turn
chanceCTLrgLoMHCu = 0.1 --- :: Fractional a => a - the chance a CTL recognizes a tumor cell with low MHC having been co-located with tumor cell for at least one turn
chanceCTLrgHiMHCu = 0.5 --- :: Fractional a => a - the chance a CTL recognizes a tumor cell with high MHC having been co-located with tumor cell for at least one turn
-------- chance of FasL involvement
chanceFasInv0 = 0.05 --- FasL cell's chance of involvement in lysis immediately after recognition (perhaps 5%)
chanceFasInv1 = 0.1 --- FasL cell's chance of involvement in lysis after 1 round of recognition (perhaps 30%)
chanceFasInv2 = 0.5 --- FasL cell's chance of involvement in lysis after 2 rounds of recognition (perhaps 100%)
--------  chance of lysis
chancePerLoMHC0 = 0.0 --- Perforin chance of lysis immediately after recognition low MHC 
chancePerHiMHC0 = 0.0 --- Perforin chance of lysis immediately after recognition high MHC 
chancePerLoMHC1 = 0.04 --- Perforin chance of lysis after 1 round low MHC 
chancePerHiMHC1 = 0.06 --- Perforin chance of lysis after 1 round high MHC
chancePerLoMHC2 = 0.15 --- Perforin chance of lysis after 2 rounds low MHC 
chancePerHiMHC2 = 0.2 --- Perforin chance of lysis after 2 rounds high MHC

--- Used in eachInteract
maxEffNum = 1000 --- :: Num a => a - the maximum numbers of cells that can recognize a tumor cell
chanceFasUsedup = 0.0 --- :: Fractional a => a - the chance that the cell's supply of FasL is exhausted
chancePerUsedup0 = 0.0 --- :: Fractional a => a - the chance that a cell's supply of perforin is exhausted right after recognizing a new tumor cell
chancePerUsedup1 = 0.0 --- :: Fractional a = a - the chance that a cell's supply of perforin is exhausted 1 time step after recognizing a tumor cell 
chancePerUsedup2 = 0.3 --- :: Fractional a => a - the chance that the cell's supply of perforin is exhausted 2 or more time steps after recognizing a tumor cell


-------------------------------------------
--- Some procedures regarding randomization
-------------------------------------------

        
randomBinom :: Int -> Double -> StdGen -> (Int, StdGen)
--- Inputs
---   n :: (Integral a) => a - the total number of trials
---   p :: Double - the probability of a success for each trial
--- Outputs
---   :: (Integral a) => a - a randomly chosen number between 0 and n distributed according to the binomial distribution with n trials and probability p of success 
randomBinom n p gen = 
    let randomBinomRec :: Int -> Int -> StdGen -> (Int, StdGen)
        randomBinomRec nRec valRec genRec = 
            if nRec==0
            then (valRec, genRec)
            else let (randn, newGen) = randomR (0::Double, 1::Double) genRec
                 in if randn == 0 || randn == 1
                    then randomBinomRec nRec valRec newGen
                    else if p < randn
                         then randomBinomRec (nRec-1) valRec newGen
                         else randomBinomRec (nRec-1) (valRec+1) newGen
    in randomBinomRec n 0 gen

-- randomPoissonIntegral :: Int -> Int -> StdGen -> (Int, StdGen)
--- Inputs
---   e :: (Integral a) => a - the expected value of the distribution
---   m :: (Integral b) => b - the maximum value that can be generated
--- Outputs
---   :: (Integral a) => a - a randomly chosen number between 0 and m according to the Poisson distrubution with mean e
randomPoisson mu mx gen = (min result mx, ngen)
    where (result, ngen) = psn mu 0 1 gen
          psn :: Int -> Int -> Double -> StdGen -> (Int, StdGen)
          psn muleft k p g =
              let (rval, newgen) = randomR (0::Double, 1::Double) g
                  interp = rval * p
                in if (interp < exp(1)) && (muleft >0)
                   then if muleft > 500
                        then let newp = interp * exp(500)
                               in if newp>1 then psn (muleft-500) (k+1) newp newgen else (k, newgen)
                        else let newp = (interp::Double) * exp(fromIntegral(muleft))
                               in if newp>1 then psn (-1) (k+1) newp newgen else (k,newgen)
                   else if interp > 1 then psn muleft (k+1) interp newgen else (k, newgen)

randomRsGen :: (Int, Int) -> Int -> [Int] -> StdGen -> ([Int],StdGen)
--- Inputs:
---   lo :: Int - lower bound
---   hi :: Int - upper bound
---   num :: Int - number of randomly generated values to take
---   lst :: [Int] - list af values so far
---   gen :: StdGen - random generator
--- Outputs:
---   :: [Int] - a list of randomly generated integers (each ranging from 'lo' to 'hi')
---   :: StdGen - an update random generator
randomRsGen (lo, hi) num lst gen = 
    if num > 0
    then let (rnum,newGen) = randomR (lo, hi) gen
         in randomRsGen (lo, hi) (num-1) (rnum:lst) newGen           
    else (lst,gen)

         
--------------------------------
--- Regarding distributing cells
---------------------------------

distributeEff :: Seq Composition -> Int -> (Composition -> Composition) -> StdGen -> (Seq Composition, StdGen)
--- Input
---   seq :: Seq Composition - a senquence of cell data indexed by their positions
---   num :: Int - the number of cells to add (must be non-negative)
---   fn :: (Composition -> Composition) - a function that increments an effector cell component of the composition
---   gen :: StdGen - a random generator
--- Global Constants
---   gridSize :: Num a => a - total number of grid elements
--- Output
---   :: Seq Composition - an updated sequence of compositions with the numbers of a given effector cell added to the grid 
---   :: StdGen - an updated random generator
distributeEff seq num fn gen = 
    if num == 0
    then (seq, gen)
    else let (newPos, newGen) = randomR (0 :: Int, gridSize-1 :: Int) gen
             cellLst = index seq newPos  --- the composition at position newPos
             !newSeq = update newPos (fn cellLst) seq --- Increment the fas CTL component of the composition at location newPos
         in distributeEff newSeq (num - 1) fn newGen


---------------------
--- Cell interaction
---------------------


eachInteract :: Composition -> StdGen -> (Composition, Int, StdGen)
--- Inputs:
---   origComp :: Composition - a composition of tumor and immune cells in a portion of a grid point interaction takes place
---   gen :: StdGen - a random generator
--- Outputs:
---   :: Composition - an updated compositions resulting from the lysis that takes place in a portion of a grid point
---   :: StdGen - an updated random generator
eachInteract origComp gen = 
    let (Comp tumLoMHC tumHiMHC fcA fcU fc1 fc2 pcA pcU pc1 pc2 fnA fnU fn1 fn2 pnA pnU pn1 pn2) = origComp
        ---
        --- first determine how many effector cells interact with the tumor cell
        ---
        (fcRptl,pcRptl,fnRptl,pnRptl,genPtl) = effCellRec (fcA+fcU) (pcA+pcU) (fnA+fnU) (pnA+pnU) (max 0 (maxEffNum - (fc1+fc2+pc1+pc2+fn1+fn2+pn1+pn2))) gen --- the largest number of new cells of each type that can interact (potentially recognize a tumor cell (Ptl for "potential"))
        (fcUptl, pcUptl, fnUptl, pnUptl) = ((min fcRptl fcU), (min pcRptl pcU), (min fnRptl fnU), (min pnRptl pnU)) --- cells that have been around longer have priority
        (fcAptl, pcAptl, fnAptl, pnAptl) = ((fcRptl-fcUptl), (pcRptl - pcUptl), (fnRptl - fnUptl), (pnRptl - pnUptl)) --- remaining cells will be those that just arrived
        hiMHC = (tumHiMHC>0) --- a boolean indicating whether or not the tumor has high MHC
        ---
        --- next determine how many of each type of effector cell actually recognize tumor cell
        ---
        (fcArec, genNewFC) = randomBinom fcAptl (if hiMHC then chanceCTLrgHiMHCa else chanceCTLrgLoMHCa) genPtl
        (pcArec, genNewPC) = randomBinom pcAptl (if hiMHC then chanceCTLrgHiMHCa else chanceCTLrgLoMHCa) genNewFC
        (fnArec, genNewFN) = randomBinom fnAptl (if hiMHC then chanceNKLrgHiMHCa else chanceNKLrgLoMHCa) genNewPC
        (pnArec, genFCptlU) = randomBinom pnAptl (if hiMHC then chanceNKLrgHiMHCa else chanceNKLrgLoMHCa) genNewPC
        (fcUrec, genPCptlU) = randomBinom fcUptl (if hiMHC then chanceCTLrgHiMHCu else chanceCTLrgLoMHCu) genFCptlU
        (pcUrec, genFNptlU) = randomBinom pcUptl (if hiMHC then chanceCTLrgHiMHCu else chanceCTLrgLoMHCu) genPCptlU
        (fnUrec, genPNptlU) = randomBinom fnUptl (if hiMHC then chanceNKLrgHiMHCu else chanceNKLrgLoMHCu) genFNptlU
        (pnUrec, genFC0inv) = randomBinom pnUptl (if hiMHC then chanceNKLrgHiMHCu else chanceNKLrgLoMHCu) genPNptlU
        ---
        --- next determine how many of each type of effector cell expresses fasL this time step
        --- 
        (fc0inv,genFC1inv) = randomBinom (fcArec+fcUrec) chanceFasInv0 genFC0inv
        (fc1inv,genFC2inv) = randomBinom fc1 chanceFasInv1 genFC1inv
        (fc2inv,genFN0inv) = randomBinom fc2 chanceFasInv2 genFC2inv
        (fn0inv,genFN1inv) = randomBinom (fnArec+fnUrec) chanceFasInv0 genFN0inv
        (fn1inv,genFN2inv) = randomBinom fn1 chanceFasInv1 genFN1inv
        (fn2inv,genFas2) = randomBinom fn2 chanceFasInv2 genFN2inv
        fasLinv = fc0inv + fc1inv + fc2inv + fn0inv + fn1inv + fn2inv
        ---
        --- determine the number of tumor cells lysed
        --- 
        tnum = if hiMHC then tumHiMHC else tumLoMHC   --- if there are hiMHC tumor cells here, then just count those
        (tumByFas, genPer0) = randomBinom tnum (probFasKill fasLinv) genFas2
        (tumByPer0, genPer1) = randomBinom (tnum - tumByFas) (probPerKill (pcArec + pcUrec + pnArec + pnUrec) 0 hiMHC) genPer0  
        (tumByPer1, genPer2) = randomBinom (tnum - tumByFas-tumByPer0) (probPerKill (pc1 + pn1) 1 hiMHC) genPer1
        (tumByPer2, genFC0sv) = randomBinom (tnum - tumByFas-tumByPer0 -tumByPer1) (probPerKill (pc2+pn2) 2 hiMHC) genPer2
        --- 
        --- effector cells becoming ineffective
        --- 
        (fc0Rm, genFC1sv) = randomBinom fc0inv chanceFasUsedup genFC0sv
        (fc1Rm, genFC2sv) = randomBinom fc1inv chanceFasUsedup genFC1sv
        (fc2Rm, genPC0sv) = randomBinom fc2inv chanceFasUsedup genFC2sv
        (pc0Rm, genPC1sv) = randomBinom (pcArec+pcUrec) chancePerUsedup0 genPC0sv
        (pc1Rm, genPC2sv) = randomBinom pc1 chancePerUsedup1 genPC1sv
        (pc2Rm, genFN0sv) = randomBinom pc2 chancePerUsedup2 genPC2sv            
        (fn0Rm, genFN1sv) = randomBinom fn0inv chanceFasUsedup genFN0sv
        (fn1Rm, genFN2sv) = randomBinom fn1inv chanceFasUsedup genFN1sv
        (fn2Rm, genPN0sv) = randomBinom fn2inv chanceFasUsedup genFN2sv
        (pn0Rm, genPN1sv) = randomBinom (pnArec+pnUrec) chancePerUsedup0 genPN0sv
        (pn1Rm, genPN2sv) = randomBinom pn1 chancePerUsedup1 genPN1sv
        (pn2Rm, genRet) = randomBinom pn2 chancePerUsedup2 genPN2sv
        effRem = fc0Rm + fc1Rm + fc2Rm + pc0Rm + pc1Rm + pc2Rm + fn0Rm + fn1Rm + fn2Rm + pn0Rm + pn1Rm + pn2Rm -- for stats
        ---
        --- record remaining values to return
        --- 
        tumLys = tumByFas + tumByPer0 + tumByPer1 + tumByPer2
        nTumLoMHC = if hiMHC then tumLoMHC else tumLoMHC - tumLys
        nTumHiMHC = if hiMHC then tumHiMHC - tumLys else tumHiMHC
     in if tumLys == 0   --- if the tumor was not lysed, then improve the status of the effector cells not removed, otherwise consider the cells no longer colocated with a tumor cell
        then let nfcA = 0 
                 npcA = 0
                 nfnA = 0
                 npnA = 0
                 nfcU = fcU+fcA-fcUrec-fcArec
                 npcU = pcU+pcA-pcUrec-pcArec
                 nfnU = fnU+fnA-fnUrec-fnArec
                 npnU = pnU+pnA-pnUrec-pnArec
                 nfc1 = fcArec+fcUrec-fc0Rm
                 npc1 = pcArec+pcUrec-pc0Rm
                 nfn1 = fnArec+fnUrec-fn0Rm
                 npn1 = pnArec+pnUrec-pn0Rm
                 nfc2 = fc2+fc1-fc2Rm-fc1Rm
                 npc2 = pc2+pc1-pc2Rm-pc1Rm
                 nfn2 = fn2+fn1-fn2Rm-fn1Rm
                 npn2 = pn2+pn1-pn2Rm-pn1Rm
                 newComp = Comp nTumLoMHC nTumHiMHC nfcA nfcU nfc1 nfc2 npcA npcU npc1 npc2 nfnA nfnU nfn1 nfn2 npnA npnU npn1 npn2
              in (newComp, effRem, genRet)
        else let nfcA = fcA+fcU+fc1+fc2 - fc0Rm - fc1Rm - fc2Rm
                 npcA = pcA+pcU+pc1+pc2 - pc0Rm - pc1Rm - pc2Rm
                 nfnA = fnA+fnU+fn1+fn2 - fn0Rm - fn1Rm - fn2Rm
                 npnA = pnA+pnU+pn1+pn2 - pn0Rm - pn1Rm - pn2Rm
                 nfcU = 0
                 npcU = 0
                 nfnU = 0
                 npnU = 0
                 nfc1 = 0
                 npc1 = 0
                 nfn1 = 0
                 npn1 = 0
                 nfc2 = 0
                 npc2 = 0
                 nfn2 = 0
                 npn2 = 0
                 newComp = Comp nTumLoMHC nTumHiMHC nfcA nfcU nfc1 nfc2 npcA npcU npc1 npc2 nfnA nfnU nfn1 nfn2 npnA npnU npn1 npn2
              in (newComp, effRem, genRet)



effCellRec :: Int -> Int -> Int -> Int -> Int -> StdGen -> (Int, Int, Int, Int, StdGen)
--- Inputs:
---   fcNmr :: Int - the number of Fas CTLs that could (without cap) recognize a tumor cell
---   pcNmr :: Int - the number of perforin CTLs that could (without cap) recognize a tumor cell
---   fnNmr :: Int - the number of Fas NKLs that could (without cap) recognize a tumor cell
---   pnNmr :: Int - the number of perforin NKLs that could (without cap) recognize a tumor cell
---   max :: Int - maximum number of remaining cells that can recognize a tumor cell
---   eGen :: StdGen - a random generator
--- Output:
---   :: Int - the actual number of Fas CTLs that could recognize a tumor cell
---   :: Int - the actual number of perforin CTLs that could recognize a tumor cell
---   :: Int - the actual number of Fas NKLs that could recognize a tumor cell
---   :: Int - the actual number of perforin NKLs that could recognize a tumor cell
---   :: StdGen - an updated random generator
effCellRec fcNmr pcNmr fnNmr pnNmr max eGen = 
    let effSum = fcNmr + pcNmr + fnNmr + pnNmr
    in if effSum>max
       then let (randnmr, nGen) = randomR (0::Int, effSum-1 :: Int) eGen
            in if randnmr < fcNmr
               then effCellRec (fcNmr-1) pcNmr fnNmr pnNmr max nGen
               else if randnmr < fcNmr+pcNmr
                    then effCellRec fcNmr (pcNmr-1) fnNmr pnNmr max nGen
                    else if randnmr < fcNmr+pcNmr+fnNmr
                         then effCellRec fcNmr pcNmr (fnNmr-1) pnNmr max nGen
                         else effCellRec fcNmr pcNmr fnNmr (pnNmr-1) max nGen
       else (fcNmr, pcNmr, fnNmr, pnNmr, eGen)  

probFasKill :: (Fractional a) => Int -> a
--- Input
---   num :: Int - The number of FasL effector cells relevant to a given tumor cell
--- Global constants used
---   scale :: Num a => a - the smallest positive number of FasL expressing cells that make the kill rate 0.
--- Output
---   :: (Fractional a) => a - The probability that a given cell will be lyced 
probFasKill num = if num>scale
                then 0
                else let x = fromIntegral num
                         s = fromIntegral scale
                in -1*x*(x-s)*(x-3*s)^2/(6*(s^4))

probPerKill :: (Integral a) => a -> a -> Bool -> Double
--- Input 
---   num :: (Integral a) => a - The number of perforin effector cells that have recognized tumor cell "age" times steps ago
---   age :: (Integral a) => a - The number of number of time steps ago the "num" effector cells have recognized a tumor cell (ranges over 0, 1, 2, with 0 representing cells that just recognized tumor cell, and 2 representing cells that recognized tumor cell at least 2 time steps age)
---   hiMHC :: Bool - attains value 1 if tumor cell expresses high level of MHC, and 0 otherwise
--- Global constant used
---   chancePerLoMHC0 :: (Fractional a) => a - Perforin chance of lysis immediately after recognition low MHC 
---   chancePerHiMHC0 :: (Fractional a) => a - Perforin chance of lysis immediately after recognition high MHC 
---   chancePerLoMHC1 :: (Fractional a) => a - Perforin chance of lysis after 1 round low MHC 
---   chancePerHiMHC1 :: (Fractional a) => a - Perforin chance of lysis after 1 round high MHC
---   chancePerLoMHC2 :: (Fractional a) => a - Perforin chance of lysis after 2 rounds low MHC 
---   chancePerHiMHC2 :: (Fractional a) => a - Perforin chance of lysis after 2 rounds high MHC
--- Output
---   :: Double - the probability that a given cell is lysed via perforin
probPerKill num age hiMHC
    | ((age == 0) && (not hiMHC)) = 1-(1-chancePerLoMHC0)^num
    | ((age == 0) && hiMHC)       = 1-(1-chancePerHiMHC0)^num
    | ((age == 1) && (not hiMHC)) = 1-(1-chancePerLoMHC1)^num
    | ((age == 1) && hiMHC)       = 1-(1-chancePerHiMHC1)^num
    | ((age == 2) && (not hiMHC)) = 1-(1-chancePerLoMHC2)^num
    | otherwise                   = 1-(1-chancePerHiMHC2)^num


-----------------------
--- spatial procedures
-----------------------

nbrPos :: Int -> [Int]
--- input
---   i :: Int input integer representing a position in the array of size arraySize
--- global constant used:
---   gridDim :: Num a => a - (cell interaction takes place in a gridDim by gridDim matrix)
--- output
---   [Int] - an array of positions of all the immediate neighbors of the input position
nbrPos i = 
    let 
        x = (rem i gridDim)  --- remainder
        y = (div i gridDim)  --- quotient
    in
      if (even y)  --- (for example, the first line is y==0)
      then 
          if y>0  --- (where y is not the first line)
          then 
              if y<gridDim-1 
              then 
                  if x>0 
                  then 
                      if x<gridDim-1
                      then [i-gridDim-1,i-gridDim,i-1,i+1,i+gridDim-1,i+gridDim]
                      else [i-gridDim-1,i-gridDim,i-1,i+gridDim-1,i+gridDim] --- (case where x==gridDim-1) 
                  else [i-gridDim,i+1,i+gridDim] --- (case where x==0)
              else --- (case where y==gridDim-1)
                  if x>0
                  then 
                      if x<gridDim-1
                      then [i-gridDim-1,i-gridDim,i-1,i+1]
                      else [i-gridDim-1,i-gridDim,i-1] --- (case where x==gridDim-1)
                  else [i-gridDim,i+1] --- (case where x==0)
          else --- (case where y==0)
              if x>0
              then if x<gridDim-1
                   then [i-1,i+1,i+gridDim-1,i+gridDim]
                   else [i-1,i+gridDim-1,i+gridDim] --- (case where x==gridDim-1)
              else [i+1,i+gridDim] --- (case where x==0)
      else --- (case (odd x), note that y cannot be the first line)
          if y<gridDim-1
          then 
              if x>0
              then 
                  if x<gridDim-1
                  then [i-gridDim+1,i-gridDim,i-1,i+1,i+gridDim,i+gridDim+1]
                  else [i-gridDim,i-1,i+gridDim] --- (case where x== gridDim-1)
              else [i-gridDim+1,i-gridDim,i+1,i+gridDim,i+gridDim+1] --- (case where x==0)
          else --- (case where y==gridDim-1)
              if x>0
              then 
                  if x<gridDim-1
                  then [i-gridDim+1,i-gridDim,i-1,i+1]
                  else [i-gridDim,i-1] --- (case where x==gridDim-1)
              else [i-gridDim+1,i-gridDim,i+1]

nbrPos2away :: Int -> [Int]
--- input
---   i :: Int input integer representing a position in the array of size arraySize
--- global constant used:
---   gridDim :: Num a => a - (cell interaction takes place in a gridDim by gridDim matrix)
--- output
---   [Int] - an array of positions of all the immediate neighbors of the input position
nbrPos2away i = 
    let 
      x = (rem i gridDim)  --- remainder
      y = (div i gridDim)  --- quotient
      iPosLst = if (even y) --- such as y==0, the even rows are shifted left and odd rows are shifted right
                then [(x-1,y-2),(x,y-2),(x+1,y-2),(x-2,y-1),(x+1,y-1),(x-2,y),(x+2,y),(x-2,y+1),(x+1,y+1),(x-1,y+2),(x,y+2),(x+1,y+2)]
                else --- (odd y)
                 [(x-1,y-2),(x,y-2),(x+1,y-2),(x-1,y-1),(x+2,y-1),(x-2,y),(x+2,y),(x-1,y+1),(x+2,y+1),(x-1,y+2),(x,y+2),(x+1,y+2)]
      f (a,b) = ((a>=0) && (a<gridDim) && (b>=0) && (b<gridDim))  --- to filter out points that don't fit on grid
      g (a,b) = b*gridDim + a     --- to convert back to single coordinate form
      fPosLst = Data.List.filter f iPosLst
     in (map g fPosLst)



nbrPos3away :: Int -> [Int]
--- input
---   i :: Int input integer representing a position in the array of size arraySize
--- global constant used:
---   gridDim :: Num a => a - (cell interaction takes place in a gridDim by gridDim matrix)
--- output
---   [Int] - an array of positions of all the immediate neighbors of the input position
nbrPos3away i = 
    let 
      x = (rem i gridDim)  --- remainder
      y = (div i gridDim)  --- quotient
      iPosLst = if (even y) --- such as y==0, the even rows are shifted left and odd rows are shifted right
                then [(x-2,y-3),(x-1,y-3),(x,y-3),(x+1,y-3),(x-2,y-2),(x+2,y-2),(x-3,y-1),(x+2,y-1),(x-3,y),(x+3,y),(x-3,y+1),(x+2,y+1),(x-2,y+2),(x+2,y+2),(x-2,y+3),(x-1,y+3),(x,y+3),(x+1,y+3)]
                else --- (odd y)
                 [(x-1,y-3),(x,y-3),(x+1,y-3),(x+2,y-3),(x-2,y-2),(x+2,y-2),(x-2,y-1),(x+3,y-1),(x-3,y),(x+3,y),(x-2,y+1),(x+3,y+1),(x-2,y+2),(x+2,y+2),(x-1,y+3),(x,y+3),(x+1,y+3),(x+2,y+3)]
      f (a,b) = ((a>=0) && (a<gridDim) && (b>=0) && (b<gridDim))  --- to filter out points that don't fit on grid
      g (a,b) = b*gridDim + a     --- to convert back to single coordinate form
      fPosLst = Data.List.filter f iPosLst
     in (map g fPosLst)



------------------------------
--- Functions on Compositions 
------------------------------

countTumLoMHC :: Composition -> Int
--- Input:
---   (Comp TumLoMHC TumHiMHC FCa FCu FCr1 FCr2 PCa PCu PCr1 PCr2 FNa FNu FNr1 FNr2 PNa PNu PNr1 PNr2) :: Composition - the cell composition at a position of the grid
--- Output:
---   tum :: Int - the number of tumor cells at that position of the grid
countTumLoMHC (Comp tumLoMHC _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _) = tumLoMHC

countTumHiMHC :: Composition -> Int
countTumHiMHC (Comp _ tumHiMHC _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _) = tumHiMHC

countFCa :: Composition -> Int
countFCa (Comp _ _ fca _ _ _ _ _ _ _ _ _ _ _ _ _ _ _) = fca

countFCu :: Composition -> Int
countFCu (Comp _ _ _ fcu _ _ _ _ _ _ _ _ _ _ _ _ _ _) = fcu

countFC1 :: Composition -> Int
countFC1 (Comp _ _ _ _ fc1 _ _ _ _ _ _ _ _ _ _ _ _ _) = fc1

countFC2 :: Composition -> Int
countFC2 (Comp _ _ _ _ _ fc2 _ _ _ _ _ _ _ _ _ _ _ _) = fc2

countPCa :: Composition -> Int
countPCa (Comp _ _ _ _ _ _ pca _ _ _ _ _ _ _ _ _ _ _) = pca

countPCu :: Composition -> Int
countPCu (Comp _ _ _ _ _ _ _ pcu _ _ _ _ _ _ _ _ _ _) = pcu

countPC1 :: Composition -> Int
countPC1 (Comp _ _ _ _ _ _ _ _ pc1 _ _ _ _ _ _ _ _ _) = pc1

countPC2 :: Composition -> Int
countPC2 (Comp _ _ _ _ _ _ _ _ _ pc2 _ _ _ _ _ _ _ _) = pc2

countFNa :: Composition -> Int
countFNa (Comp _ _ _ _ _ _ _ _ _ _ fna _ _ _ _ _ _ _) = fna

countFNu :: Composition -> Int
countFNu (Comp _ _ _ _ _ _ _ _ _ _ _ fnu _ _ _ _ _ _) = fnu

countFN1 :: Composition -> Int
countFN1 (Comp _ _ _ _ _ _ _ _ _ _ _ _ fn1 _ _ _ _ _) = fn1

countFN2 :: Composition -> Int
countFN2 (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ fn2 _ _ _ _) = fn2

countPNa :: Composition -> Int
countPNa (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ _ pna _ _ _) = pna

countPNu :: Composition -> Int
countPNu (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ pnu _ _) = pnu

countPN1 :: Composition -> Int
countPN1 (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ pn1 _) = pn1

countPN2 :: Composition -> Int
countPN2 (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ pn2) = pn2

countTum :: Composition -> Int
countTum comp = 
    (countTumLoMHC comp) + (countTumHiMHC comp)
    
countEff :: Composition -> Int
countEff (Comp _ _ fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
    fca + fcu + fc1 + fc2 + pca + pcu + pc1 + pc2 + fna + fnu + fn1 + fn2 + pna + pnu + pn1 + pn2

countFC :: Composition -> Int
countFC (Comp _ _ fca fcu fc1 fc2 _ _ _ _ _ _ _ _ _ _ _ _) =
    fca + fcu + fc1 + fc2

countPC :: Composition -> Int
countPC (Comp _ _ _ _ _ _ pca pcu pc1 pc2 _ _ _ _ _ _ _ _) = 
    pca + pcu + pc1 + pc2

countFN :: Composition -> Int
countFN (Comp _ _ _ _ _ _ _ _ _ _ fna fnu fn1 fn2 _ _ _ _) = 
    fna + fnu + fn1 + fn2

countPN :: Composition -> Int
countPN (Comp _ _ _ _ _ _ _ _ _ _ _ _ _ _ pna pnu pn1 pn2) = 
    pna + pnu + pn1 + pn2

adjustTumLoMHC :: (Int -> Int) -> Composition -> Composition
--- Inputs:
---   f :: (Int -> Int) - a function used to adjust the number of low MHC tumor cells
---   Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) :: Composition - the cell composition at a position of the grid
--- Output: 
---  :: Composition - a new composition that is the same as the old, except the number of tumor cells has been changed according the function f
adjustTumLoMHC f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp (f tumLoMHC) tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustTumHiMHC :: (Int -> Int) -> Composition -> Composition
adjustTumHiMHC f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC (f tumHiMHC) fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustFCa :: (Int -> Int) -> Composition -> Composition
adjustFCa f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
  (Comp tumLoMHC tumHiMHC (f fca) fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) 

adjustFCu :: (Int -> Int) -> Composition -> Composition
adjustFCu f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
  (Comp tumLoMHC tumHiMHC fna (f fcu) fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) 

adjustFC1 :: (Int -> Int) -> Composition -> Composition
adjustFC1 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) =  
   (Comp tumLoMHC tumHiMHC fca fcu (f fc1) fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustFC2 :: (Int -> Int) -> Composition -> Composition
adjustFC2 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 (f fc2) pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustPCa :: (Int -> Int) -> Composition -> Composition
adjustPCa f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 (f pca) pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)
   
adjustPCu :: (Int -> Int) -> Composition -> Composition
adjustPCu f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca (f pcu) pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustPC1 :: (Int -> Int) -> Composition -> Composition
adjustPC1 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu (f pc1) pc2 fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustPC2 :: (Int -> Int) -> Composition -> Composition
adjustPC2 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 (f pc2) fna fnu fn1 fn2 pna pnu pn1 pn2)

adjustFNa :: (Int -> Int) -> Composition -> Composition
adjustFNa f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 (f fna) fnu fn1 fn2 pna pnu pn1 pn2)
   
adjustFNu :: (Int -> Int) -> Composition -> Composition
adjustFNu f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna(f fnu) fn1 fn2 pna pnu pn1 pn2)

adjustFN1 :: (Int -> Int) -> Composition -> Composition
adjustFN1 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca pcu fc1 fc2 pca pcu pc1 pc2 fna pnu (f fn1) fn2 pna pnu pn1 pn2)

adjustFN2 :: (Int -> Int) -> Composition -> Composition
adjustFN2 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 (f fn2) pna pnu pn1 pn2)

adjustPNa :: (Int -> Int) -> Composition -> Composition
adjustPNa f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 (f pna) pnu pn1 pn2)

adjustPNu :: (Int -> Int) -> Composition -> Composition
adjustPNu f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna (f pnu) pn1 pn2)
   
adjustPN1 :: (Int -> Int) -> Composition -> Composition
adjustPN1 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) =
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu (f pn1) pn2)

adjustPN2 :: (Int -> Int) -> Composition -> Composition
adjustPN2 f (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = 
   (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 (f pn2))

sumComp :: Composition -> Composition -> Composition
sumComp (Comp tumLoMHC_1 tumHiMHC_1 fca_1 fcu_1 fc1_1 fc2_1 pca_1 pcu_1 pc1_1 pc2_1 fna_1 fnu_1 fn1_1 fn2_1 pna_1 pnu_1 pn1_1 pn2_1) (Comp tumLoMHC_2 tumHiMHC_2 fca_2 fcu_2 fc1_2 fc2_2 pca_2 pcu_2 pc1_2 pc2_2 fna_2 fnu_2 fn1_2 fn2_2 pna_2 pnu_2 pn1_2 pn2_2) =
 (Comp (tumLoMHC_1 + tumLoMHC_2) (tumHiMHC_1+tumHiMHC_2) (fca_1+fca_2) (fcu_1 + fcu_2) (fc1_1 + fc1_2) (fc2_1 + fc2_2) (pca_1+pca_2) (pcu_1 + pcu_2) (pc1_1 + pc1_2) (pc2_1 + pc2_2) (fna_1+fna_2) (fnu_1 + fnu_2) (fn1_1 + fn1_2) (fn2_1 + fn2_2) (pna_1+pna_2) (pnu_1 + pnu_2) (pn1_1 + pn1_2) (pn2_1 + pn2_2))

--------------------------
--- Writing Position Data 
--------------------------


writePosData :: Seq Composition  -> Int -> Handle -> IO()
--- Inputs:
---   posSeq :: Seq Composition - a sequence of position compositions
---   pos :: Int - position to be writtin at this point
---   handle :: Handle - a handle for the file to write data in.
--- Effects:
---   Write a line of position data (corresponding to the cells at a particular position)
writePosData seq pos handle = do
  if pos < gridSize
  then do
    let compHere = index seq pos 
    let (Comp tumLoMHC tumHiMHC fca fcu fc1 fc2 pca pcu pc1 pc2 fna fnu fn1 fn2 pna pnu pn1 pn2) = compHere
    hPutStrLn handle $ (show tumLoMHC) ++ " " ++ (show tumHiMHC) ++ " " ++ (show fca) ++ " " ++ (show fcu) ++ " " ++ (show fc1) ++ " " ++ (show fc2) ++ " " ++ (show pca) ++ " " ++ (show pcu) ++ " " ++ (show pc1) ++ " " ++ (show pc2) ++ " " ++ (show fna) ++ " " ++ (show fnu) ++ " " ++ (show fn1) ++ " " ++ (show fn2) ++ " " ++ (show pna) ++ " " ++ (show pnu) ++ " " ++ (show pn1) ++ " " ++ (show pn2)
    writePosData seq (pos+1) handle
  else return ()