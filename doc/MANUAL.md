MANUAL FOR USER
===============

(date: 12/December/2022)

The output files of this analysis are organized as follows:

```
<analysis-dir>/output/<sim-set>/<rn>/SexaquarkResults_<v0-option>_<index>.root
<analysis-dir>/output/<sim-set>/<rn>/AnalysisResults_<v0-option>_<index>.root
```

where:

* `<analysis-dir>` : main analysis directory, `/home/ceres/borquez/work/analysis`
* `<sim-set>` : could be `data`, `only_bkg`, `only_V0s`, `signal+bkg`
* `<rn>` : run number
* `<v0-option>` : `TrueV0s`, `CustomV0s`, `OfficialV0s`.
* `<index>` : 3-digit number

`SexaquarkResults.root` contain a `TTree` called `"Sexaquarks"`, which contains the information of the formed anti-sexaquark candidates and both of the V0s that form it, assuming the **anti-sexaquark + neutron -> neutral kaon short + anti-lambda** reaction channel.

`AnalysisResults.root` contain a `TList` called `"Trees"`, which contains a `TTree` called `"Events"`, which contains the information of every generated event: the MC particles, the reconstructed tracks and the V0s formed with the chosen `<v0-option>`.

There are certain branches in the `"Sexaquarks"` tree (within a `SexaquarkResults_<v0-option>_<index>.root` file) that index to a position in the vectors of the `"Events"` tree (within the analogous `AnalysisResults_<v0-option>_<index>.root`).

## How to read a `SexaquarkResults` file
---

The structure of the `Sexaquarks` tree is a **flat format**, i.e., each entry (or row) corresponds to a different anti-sexaquark candidate.

The branches of this tree, that contain the information of the **reconstructed candidate**, are:

* `event`    : event number (entry number in the analogous `AnalysisResults_<v0-option>_<index>.root` file)
* `Idx_V0A`  : index of first daughter (in the V0s vector at `AnalysisResults_<v0-option>_<index>.root`)
* `Idx_V0B`  : index of second daughter (in the V0s vector at `AnalysisResults_<v0-option>_<index>.root`)
* `E`        : energy
* `Px`       : x-component of the momentum
* `Py`       : y-component of the momentum
* `Pz`       : z-component of the momentum
* `M`        : inv. mass, **reconstructed by assuming nucleon at rest**. When kinematics don't allow the subtraction of the neutron invariant mass in the conservation of energy (resulting in a negative value inside the square root), the mass is arbitrarily set to -1.
* `X`        : x-coordinate of secondary vertex
* `Y`        : y-coordinate of secondary vertex
* `Z`        : z-coordinate of secondary vertex
* `DCA`      : distance of closest approach between the V0s that form the anti-sexaquark, after being propagated to a common vertex
* `isSignal` : `kTRUE` if signal, `kFALSE` if background

The remaining branches of this tree, that contain the information of the **V0s** (`V0A` and `V0B`) that form the current candidate on that row, are:

* `V0<A,B>_Idx_Pos`       : index of positive daughter (in the tracks vector at `AnalysisResults_<v0-option>_<index>.root`)
* `V0<A,B>_Idx_Neg`       : index of negative daughter (in the tracks vector at `AnalysisResults_<v0-option>_<index>.root`)
* `V0<A,B>_Px`            : x-component of the momentum
* `V0<A,B>_Py`            : y-component of the momentum
* `V0<A,B>_Pz`            : z-component of the momentum
* `V0<A,B>_X`             : x-coordinate of V0
* `V0<A,B>_Y`             : y-coordinate of V0
* `V0<A,B>_Z`             : z-coordinate of V0
* `V0<A,B>_Pos_Px`        : x-component of positive track momentum at V0 position
* `V0<A,B>_Pos_Py`        : y-component of positive track momentum at V0 position
* `V0<A,B>_Pos_Pz`        : z-component of positive track momentum at V0 position
* `V0<A,B>_Neg_Px`        : x-component of negative track momentum at V0 position
* `V0<A,B>_Neg_Py`        : y-component of negative track momentum at V0 position
* `V0<A,B>_Neg_Pz`        : z-component of negative track momentum at V0 position
* `V0<A,B>_isSignal`      : `kTRUE` if signal, `kFALSE` if background
* `V0<A,B>_E_asK0`        : energy of V0, assuming it's a K0
* `V0<A,B>_E_asAL`        : energy of V0, assuming it's an anti-lambda
* `V0<A,B>_couldBeK0`     : `kTRUE` if K0 candidate, `kFALSE` if not
* `V0<A,B>_couldBeAL`     : `kTRUE` if anti-lambda candidate, `kFALSE` if not
* `V0<A,B>_onFlyStatus`   : `kTRUE` if comes from the on-the-fly finder, `kFALSE` if not
* `V0<A,B>_Chi2`          : chi2 value from the Kalman Filter (?)
* `V0<A,B>_DCA_Daughters` : distance of closest approach between daughters
* `V0<A,B>_IP_wrtPV`      : impact parameter w.r.t. Primary Vertex
* `V0<A,B>_CPA_wrtPV`     : cosine of pointing angle w.r.t. Primary Vertex
* `V0<A,B>_ArmAlpha`      : Armenteros-Podolanski variable alpha
* `V0<A,B>_ArmPt`         : Armenteros-Podolanski variable Pt
* `V0<A,B>_DecayLength`   : distance between PV and V0

You can also access information of the positive and negative daughters of each V0, with the following branches:

* `V0<A,B>_<Pos,Neg>_isDuplicate` : `kTRUE` if the particle is a duplicate particle, `kFALSE` if it was the best match with its respective MC
* `V0<A,B>_<Pos,Neg>_Rec_Px`      : x-component of the reconstructed momentum of the track
* `V0<A,B>_<Pos,Neg>_Rec_Py`      : y-component of the reconstructed momentum of the track
* `V0<A,B>_<Pos,Neg>_Rec_Pz`      : z-component of the reconstructed momentum of the track

## How to read an `AnalysisResults` file
---

The structure of the `Events` tree is in a **vector format**, which you can load from the `TList` called `Trees`.

Variables and vectors with MC particle information:

 * `N_MCGen`     : (int) number of MC particles, which corresponds to the size of each vector in this paragraph
 * `MC_Px`       : (vector of float) x-component of true momentum
 * `MC_Py`       : (vector of float) y-component of true momentum
 * `MC_Pz`       : (vector of float) z-component of true momentum
 * `MC_X`        : (vector of float) x-coordinate of generation vertex
 * `MC_Y`        : (vector of float) y-coordinate of generation vertex
 * `MC_Z`        : (vector of float) z-coordinate of generation vertex
 * `MC_PID`      : (vector of int) PDG code
 * `MC_Mother`   : (vector of int) index of mother, 0 if anti-sexaquark, -1 if particle is a primary, -2 if mother not relevant
 * `MC_FirstDau` : (vector of int) index of first daughter, -1 if no daughters, -2 if daughter is not relevant
 * `MC_LastDau`  : (vector of int) index of last daughter, -1 if no daughters, -2 if daughter is not relevant
 * `MC_Gen`      : (vector of int) generation (to diff. between daughters and grandaughters) (only valid for signal particles)
 * `MC_Status`   : (vector of int) MC status code
 * `MC_isSignal` : (vector of bool) kTRUE if it belongs to anti-sexaquark signal, kFALSE if background

Variables and vectors of MC Rec. Tracks:

* `N_MCRec`          : (int) number of MC reconstructed tracks, which corresponds to the size of each vector in this paragraph
* `Idx_True`         : (vector of int) index of true MC particle
* `Rec_Px`           : (vector of float) x-component of reconstructed momentum
* `Rec_Py`           : (vector of float) y-component of reconstructed momentum
* `Rec_Pz`           : (vector of float) z-component of reconstructed momentum
* `Rec_Charge`       : (vector of short) measured charge
* `Rec_NSigmaPion`   : (vector of float) absolute value of likeness to be a charged pion (closer to 0, the most likely)
* `Rec_NSigmaKaon`   : (vector of float) absolute value of likeness to be a charged kaon
* `Rec_NSigmaProton` : (vector of float) absolute value of likeness to be a proton or anti-proton
* `Rec_isDuplicate`  : (vector of bool) kTRUE if track is a duplicate, kFALSE if not
* `Rec_isSignal`     : (vector of bool) kTRUE if it belongs to anti-sexaquark signal, kFALSE if background

Variables and vectors of V0s:

* `N_V0s`            : (int) number of formed V0s, which corresponds to the size of each vector in this paragraph
* `Idx_Pos`          : (vector of int) index of positive daughter
* `Idx_Neg`          : (vector of int) index of negative daughter
* `V0_Px`            : (vector of float) x-component of V0 momentum
* `V0_Py`            : (vector of float) y-component of V0 momentum
* `V0_Pz`            : (vector of float) z-component of V0 momentum
* `V0_X`             : (vector of float) x-coordinate of V0
* `V0_Y`             : (vector of float) y-coordinate of V0
* `V0_Z`             : (vector of float) z-coordinate of V0
* `Pos_Px`           : (vector of float) x-component of positive track momentum at V0 position
* `Pos_Py`           : (vector of float) y-component of positive track momentum at V0 position
* `Pos_Pz`           : (vector of float) z-component of positive track momentum at V0 position
* `Neg_Px`           : (vector of float) x-component of negative track momentum at V0 position
* `Neg_Py`           : (vector of float) y-component of negative track momentum at V0 position
* `Neg_Pz`           : (vector of float) z-component of negative track momentum at V0 position
* `V0_isSignal`      : (vector of bool) kTRUE if signal, kFALSE if background
* `V0_E_asK0`        : (vector of float) energy of V0, assuming it's a K0
* `V0_E_asAL`        : (vector of float) energy of V0, assuming it's an anti-lambda
* `V0_couldBeK0`     : (vector of bool) kTRUE if K0 candidate, kFALSE if not
* `V0_couldBeAL`     : (vector of bool) kTRUE if anti-lambda candidate, kFALSE if not
* `V0_onFlyStatus`   : (vector of bool) kTRUE if comes from the on-the-fly finder, kFALSE if not
* `V0_Chi2`          : (vector of float) chi2 value from the Kalman Filter (?)
* `V0_DCA_Daughters` : (vector of float) distance of closest approach between daughters
* `V0_IP_wrtPV`      : (vector of float) impact parameter w.r.t. Primary Vertex
* `V0_CPA_wrtPV`     : (vector of float) cosine of pointing angle w.r.t. Primary Vertex
* `V0_ArmAlpha`      : (vector of float) Armenteros-Podolanski variable alpha
* `V0_ArmPt`         : (vector of float) Armenteros-Podolanski variable Pt
* `V0_DecayLength`   : (vector of float) distance between PV and V0

Using C++ and ROOT, you can follow the example of `/home/ceres/borquez/work/analysis/SexaquarkFinder.cxx` and `/home/ceres/borquez/work/analysis/macros/include/TreeFunctions.hxx` on how to read the `Events` vectors and variables. It would be also possible to read it using Python and Uproot, but that's work in progress...
