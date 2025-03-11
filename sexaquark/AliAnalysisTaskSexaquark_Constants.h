#ifndef TASKSEXAQUARK_CONSTANTS_H
#define TASKSEXAQUARK_CONSTANTS_H

#include <cstddef>

#include "RtypesCore.h"

/** Constants **/

namespace PdgCode {
const Short_t AntiLambda = -3122;
const Short_t Lambda = 3122;
const Short_t KaonZeroShort = 310;
const Short_t AntiNeutron = -2112;
const Short_t Neutron = 2112;
const Short_t AntiProton = -2212;
const Short_t Proton = 2212;
const Short_t NegKaon = -321;
const Short_t PosKaon = 321;
const Short_t PiMinus = -211;
const Short_t PiPlus = 211;
}  // namespace PdgCode

namespace SexaConst {
const size_t PV_CovMatrix_Size = 6;
}  // namespace SexaConst

/** Cuts **/

namespace SexaCuts {

namespace Event {
const Double_t AbsMax_PV_Zv = 12.;
}  // namespace Event

namespace Track {
const Double_t Min_Pt = 1E-2;  // 10 MeV
const Double_t Max_Pt = 1E2;   // 100 GeV
const Double_t AbsMax_PID_NSigma = 3.;
const Double_t AbsMax_Eta = 1.;
const UShort_t Min_NTPCClusters = 50;
const Double_t Max_Chi2PerNTPCClusters = 2.;
const Bool_t TurnedOn_StatusCuts = kTRUE;
const Bool_t TurnedOn_RejectKinks = kFALSE;
const Double_t AbsMin_DCAxy_wrtPV = 2.;
}  // namespace Track

namespace Lambda {
/* kinematics */
const Double_t Min_Pt = 2.;
const Double_t Min_Mass = 1.1;
const Double_t Max_Mass = 1.13;
const Double_t AbsMax_Eta = 0.9;
const Double_t Min_CPAwrtPV = 0.35;
const Double_t Max_CPAwrtPV = 0.9;
const Double_t Min_DCAwrtPV = 40.;
const Double_t AbsMax_ArmQtOverAlpha = 0.2;
/* geometric */
const Double_t AbsMax_Zv = 50.;
const Double_t Min_Radius = 45.;
const Double_t Max_Radius = 140.;
const Double_t Max_DCAnegV0 = 0.2;
const Double_t Max_DCAposV0 = 0.15;
const Double_t Max_DCAbtwDau = 0.3;
}  // namespace Lambda

namespace KaonZeroShort {
/* kinematics */
const Double_t Min_Pt = 1.5;
const Double_t Min_Mass = 0.475;
const Double_t Max_Mass = 0.525;
const Double_t AbsMax_Eta = 0.8;
const Double_t Min_CPAwrtPV = 0.35;
const Double_t Max_CPAwrtPV = 0.9;
const Double_t Min_DCAwrtPV = 20.;
/* geometric */
const Double_t AbsMax_Zv = 50.;
const Double_t Min_Radius = 20.;
const Double_t Max_Radius = 180.;
const Double_t Max_DCAnegV0 = 0.2;
const Double_t Max_DCAposV0 = 0.2;
const Double_t Max_DCAbtwDau = 0.2;
}  // namespace KaonZeroShort

namespace ChannelA {
/* kinematics-dependent */
const Double_t AbsMax_Rapidity = 0.7;
const Double_t Min_CPAwrtPV = 0.9;
const Double_t Max_CPAwrtPV = 1.;
const Double_t Min_MassAsDecay = 0.;
const Double_t Max_MassAsDecay = 6.;
/* geometry-exclusive */
const Double_t Min_Radius = 20.;
const Double_t Max_Radius = 150.;
const Double_t Max_DCALaSV = 10.;
const Double_t Max_DCALaNegSV = 10.;
const Double_t Max_DCALaPosSV = 10.;
const Double_t Max_DCAK0SV = 10.;
const Double_t Max_DCAK0NegSV = 10.;
const Double_t Max_DCAK0PosSV = 10.;
const Double_t Max_DCAbtwV0s = 10.;
const Double_t Max_DecayLengthLa = 100.;
const Double_t Max_DecayLengthK0 = 100.;
}  // namespace ChannelA

namespace ChannelD {
/* kinematics-dependent */
const Double_t AbsMax_Rapidity = 0.6;
const Double_t Min_CPAwrtPV = 0.99;
const Double_t Max_CPAwrtPV = 1.;
/* geometry-exclusive */
const Double_t Min_Radius = 60.;
const Double_t Max_Radius = 170.;
const Double_t Max_DCALaSV = 2.;
const Double_t Max_DCALaNegSV = 2.;
const Double_t Max_DCALaPosSV = 4.;
const Double_t Max_DCAKaSV = 1.5;
const Double_t Max_DCAKaLa = 1.;
const Double_t Max_DCALaNegKa = 0.5;
const Double_t Max_DCALaPosKa = 2.;
}  // namespace ChannelD

}  // namespace SexaCuts

#endif  // TASKSEXAQUARK_CONSTANTS_H
