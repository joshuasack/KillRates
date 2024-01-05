--- ************************
--- **** main module ****
--- ************************
---
--- To switch between between in vivo and in vitro program, all you need to do is add and remove comments.
---
--- If you want to run the in vivo program, leave the program as is:
--- make sure that the in vivo lines "import InVivo" and "main = do inVivo" are uncommented
--- and the in vitro lines "import InVitro" and "main = do inVitro" are commented out using double dashes
---
--- If you want to run the in vitro program, 
--- comment out the in vivo lines and uncomment the in vitro lines:
--- make sure that the in vitro lines "import InVitro" and "main = do inVitro" are uncommented
--- and the in vivo lines "import InVivo" and "main = do inVivo" are commented out using double dashes

module Main where

import InVivo
--import InVitro

main :: IO ()
main = do inVivo
--main = do inVitro