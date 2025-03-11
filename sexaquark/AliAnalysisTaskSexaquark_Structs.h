#ifndef TASKSEXAQUARK_STRUCTS_H
#define TASKSEXAQUARK_STRUCTS_H

#include "RtypesCore.h"

#include "Math/Point3D.h"
#include "Math/Vector4D.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFParticle.h"

using PxPyPzMVector = ROOT::Math::PxPyPzMVector;
using XYZPoint = ROOT::Math::XYZPoint;

/** Structs **/

struct MC_V0 {
    Bool_t is_signal{};
    UInt_t reaction_id{};
    Bool_t is_hybrid{};
};

struct KF_V0 {
    KF_V0() : kf(), kf_neg(), kf_pos(), lv(), lv_neg(), lv_pos(), v3() {}

    KFParticle kf, kf_neg, kf_pos;
    PxPyPzMVector lv, lv_neg, lv_pos;
    XYZPoint v3;

    Int_t idx{}, idx_neg{}, idx_pos{};

    Double_t impar_neg[2]{};
    Double_t param_d_neg_v0{};
    Float_t dca_neg_v0{};
    Float_t dcaxy_neg_v0{};

    Double_t impar_pos[2]{};
    Double_t param_d_pos_v0{};
    Float_t dca_pos_v0{};
    Float_t dcaxy_pos_v0{};

    Double_t cpa_wrt_pv{};
    Double_t dca_wrt_pv{};
    Double_t arm_qt{};
    Double_t arm_alpha{};
    Float_t dca_btw_dau{};
};

#endif  // TASKSEXAQUARK_STRUCTS_H
