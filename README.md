analysis
========

These files must be renamed and added into the `PWGLF/NUCLEX/Exotica/Sexaquark` directory of [**AliPhysics**](https://github.com/alisw/AliPhysics):

* `sexaquark/AddSexaquark.C`
* `sexaquark/AliAnalysisTaskSexaquark.C`
* `sexaquark/AliAnalysisTaskSexaquark.C`
* `sexaquark/runAnalysis.C`
* (and respective `macros/`)

These files must be added into the `PWGLF/RESONANCES/???` directory of [**AliPhysics**](https://github.com/alisw/AliPhysics):

* `lambda-1520/AddLambda1520Lpipi.C`
* `lambda-1520/AliAnalysisTaskLambda1520Lpipi.C`
* `lambda-1520/AliAnalysisTaskLambda1520Lpipi.C`
* `lambda-1520/runAnalysis.C`
* (and respective `macros/`)

## Testing

* **Local test**

  1. Enter the environment:
     ```
     alienv enter AliRoot/latest
     ```
  3. Download OCDB file via:
     ```
     alien.py cp alien:///alice/sim/2020/LHC20e3/OCDB/297595/OCDBsim.root file://OCDBsim.root
     ```
  4. Execute:
     ```
     $ALIDPG_ROOT/bin/aliroot_dpgsim.sh --run 297595 --mode sim --uid 1 --nevents 1 --generator PWGLF:Hijing_Sexaquark:A1.8 --simulation SimulationDefaultIonTail --system Pb-Pb --detector NoAD
     ```
