#ifndef ALIANALYSISTASKSEXAQUARK_H
#define ALIANALYSISTASKSEXAQUARK_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "TArray.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TGrid.h"
#include "TH1.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "AliAnalysisManager.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDv0.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFParticle.h"
#include "KFVertex.h"

#include "AliAnalysisTaskSexaquark_Constants.h"

using namespace ROOT::Math::VectorUtil;
using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;
using PxPyPzMVector = ROOT::Math::PxPyPzMVector;
using PxPyPzEVector = ROOT::Math::PxPyPzEVector;

class AliPIDResponse;
class KFParticle;
class KFVertex;

class AliAnalysisTaskSexaquark : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskSexaquark();
    AliAnalysisTaskSexaquark(const char* name);
    virtual ~AliAnalysisTaskSexaquark();

    /* Settings ~ stored in Analysis Manager */
    void Initialize(Bool_t is_mc, Bool_t is_signal_mc);
    void DefineCuts_Tracks();
    void DefineCuts_V0s();
    void DefineCuts_TypeA();
    void DefineCuts_TypeD();

    /* Main ~ executed at runtime */
    void UserCreateOutputObjects();
    Bool_t UserNotify();
    void UserExec(Option_t* option);
    void Terminate(Option_t* option) {}

    /* Tree */
    void AssociateBranches_Events();
    void AssociateBranches_Injected();
    void AssociateBranches_V0s();
    void AssociateBranches_TypeA();
    void AssociateBranches_TypeD();

    /* Events */
    Bool_t ProcessEvent();
    Bool_t PassesEventSelection();

    /* MC Particles */
    void ProcessMCParticles();

    /* Injected */
    void ProcessInjected();
    void BringSignalLogs();
    Bool_t LoadSignalLogs();
    void ClearSignalLogs();

    /* Tracks */
    void ProcessTracks();
    Bool_t PassesTrackSelection(const AliESDtrack* track);

    /* V0s */
    void KF_FindV0s(Short_t pdg_code_v0, Short_t pdg_code_neg, Short_t pdg_code_pos);
    Bool_t PassesV0CutsAs(Short_t pdg_code_v0, const KFParticle& kf_v0, const KFParticle& kf_neg, const KFParticle& kf_pos,
                          const PxPyPzMVector& lv_v0, const PxPyPzMVector& lv_neg, const PxPyPzMVector& lv_pos);
    void StoreV0As(Short_t pdg_code_v0, Int_t esd_idx_neg, Int_t esd_idx_pos, const KFParticle& kf_v0, const KFParticle& kf_neg,
                   const KFParticle& kf_pos, const PxPyPzMVector& lv_v0, const PxPyPzMVector& lv_neg, const PxPyPzMVector& lv_pos);

    /* Sexaquarks */
    /* -- Channel A : AntiSexaquark Neutron -> AntiLambda KaonZeroShort */
    void KF_FindSexaquarks_TypeA(Short_t pdg_struck_nucleon, const std::vector<Short_t>& pdg_reaction_products);
    Bool_t PassesSexaquarkCuts_TypeA(const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa,
                                     const PxPyPzEVector& lv_sexa_asdecay,                                                  //
                                     const KFParticle& kf_v0a, const KFParticle& kf_v0a_neg, const KFParticle& kf_v0a_pos,  //
                                     const KFParticle& kf_v0b, const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos);
    void StoreSexaquark_TypeA(Int_t idx_v0a, Int_t idx_v0b, const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa, const KFParticle& kf_v0a,
                              const PxPyPzEVector& lv_v0a, const KFParticle& kf_v0a_neg, const KFParticle& kf_v0a_pos, const KFParticle& kf_v0b,
                              const PxPyPzEVector& lv_v0b, const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos);
    /* -- Channel D : AntiSexaquark Proton -> AntiLambda K+ */
    void KF_FindSexaquarks_TypeD(Short_t pdg_struck_nucleon, const std::vector<Short_t>& pdg_reaction_products);
    Bool_t PassesSexaquarkCuts_TypeD(const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa, const KFParticle& kf_v0, const KFParticle& kf_v0_neg,
                                     const KFParticle& kf_v0_pos, const KFParticle& kf_ka);
    void StoreSexaquark_TypeD(Int_t idx_v0, Int_t esd_idx_ka, const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa, const KFParticle& kf_v0,
                              const PxPyPzEVector& lv_v0, const KFParticle& kf_v0_neg, const KFParticle& kf_v0_pos, const KFParticle& kf_ka,
                              const PxPyPzEVector& lv_ka);

    /* Utilities */
    /* -- Math */
    inline Double_t CosinePointingAngle(const XYZVector& mom, const XYZPoint& pos, const XYZPoint& ref) { return TMath::Cos(Angle(mom, pos - ref)); }
    inline Double_t ArmenterosAlpha(const XYZVector& mom_v0, const XYZVector& mom_neg, const XYZVector& mom_pos) {
        Double_t lQlNeg = mom_neg.Dot(mom_v0) / mom_v0.R();
        Double_t lQlPos = mom_pos.Dot(mom_v0) / mom_v0.R();
        if (TMath::Abs(lQlPos + lQlNeg) < 1E-6) return 2.;  // protection
        return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
    }
    inline Double_t ArmenterosQt(const XYZVector& mom_v0, const XYZVector& mom_neg) { return Perp(mom_v0, mom_neg); }
    inline Double_t LinePointDCA(const XYZVector& mom, const XYZPoint& pos, const XYZPoint& ref) { return (ref - pos).Cross(mom).R() / mom.R(); }
    /* -- Kalman Filter */
    KFParticle CreateKFParticle(const AliExternalTrackParam* track_param, Double_t mass, Int_t charge);
    KFVertex CreateKFVertex(const AliVVertex* vertex);
    /* -- Containers */
    void ClearBranches_Injected();
    void ClearBranches_V0s();
    void ClearBranches_TypeA();
    void ClearBranches_TypeD();
    void ClearContainers();

   private:
    /* Settings ~ stored in Analysis Manager ~ all persistent */
    Bool_t fIsMC;        // kTRUE if MC simulation, kFALSE if data
    Bool_t fIsSignalMC;  // kTRUE to read and load signal logs

    /* AliRoot Objects */
    AliPIDResponse* fPIDResponse;  //! pid response object
    /* -- MC */
    AliMCEvent* fMC;                      //! MC event
    const AliVVertex* fMC_PrimaryVertex;  //! MC gen. (or true) primary vertex
    /* -- Event */
    AliESDEvent* fESD;                   //! reconstructed event
    const AliESDVertex* fPrimaryVertex;  //! primary vertex
    XYZPoint v3_pv;                      //! primary vertex (optimization: load once as XYZPoint)
    AliEventCuts fEventCuts;             //! event cuts

    /* Signal Logs */
    TString fAliEnPath;              //! loaded in `UserNotify()`
    TString fSignalLog_NewBasename;  //!

    /* Utilities */
    std::unordered_map<UInt_t, std::vector<UInt_t>> fReactionID_;     //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fSexaquark_Px_;  //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fSexaquark_Py_;  //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fSexaquark_Pz_;  //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fNucleon_Px_;    //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fNucleon_Py_;    //! key: event_n
    std::unordered_map<UInt_t, std::vector<Float_t>> fNucleon_Pz_;    //! key: event_n

    /* Tree */
    TTree* fOutputTree;  //!
    /* -- Event properties */
    UInt_t fRunNumber;                                                      //! run number
    UInt_t fDirNumber;                                                      //! directory number
    UInt_t fDirNumberB;                                                     //! directory number
    UInt_t fEventNumber;                                                    //! event number
    Float_t fCentrality;                                                    //! centrality percentile
    Float_t fMagneticField;                                                 //! magnetic field
    Float_t tEvents_MC_PV_Xv;                                               //!
    Float_t tEvents_MC_PV_Yv;                                               //!
    Float_t tEvents_MC_PV_Zv;                                               //!
    Float_t tEvents_PV_Xv;                                                  //!
    Float_t tEvents_PV_Yv;                                                  //!
    Float_t tEvents_PV_Zv;                                                  //!
    std::array<Float_t, SexaConst::PV_CovMatrix_Size> tEvent_PV_CovMatrix;  //!
    Int_t tEvents_NTracks;                                                  //!
    Int_t tEvents_NTPCClusters;                                             //!
    /* -- KF */
    KFVertex kf_pv;  //!
    /* -- Signal Reaction properties */
    UShort_t tInjected_Nucleon_PdgCode;          //!
    Double_t tInjected_Mass;                     //!
    std::vector<UShort_t> tInjected_ReactionID;  //!
    std::vector<Float_t> tInjected_Px;           //!
    std::vector<Float_t> tInjected_Py;           //!
    std::vector<Float_t> tInjected_Pz;           //!
    std::vector<Float_t> tInjected_Nucleon_Px;   //!
    std::vector<Float_t> tInjected_Nucleon_Py;   //!
    std::vector<Float_t> tInjected_Nucleon_Pz;   //!
    std::vector<Float_t> tInjected_Xv;           //!
    std::vector<Float_t> tInjected_Yv;           //!
    std::vector<Float_t> tInjected_Zv;           //!
    std::vector<Float_t> tInjected_Post_Px;      //!
    std::vector<Float_t> tInjected_Post_Py;      //!
    std::vector<Float_t> tInjected_Post_Pz;      //!
    /* -- V0s properties */
    std::map<Short_t, std::vector<Int_t>> tV0_Idx;          //!
    std::map<Short_t, std::vector<Float_t>> tV0_Px;         //!
    std::map<Short_t, std::vector<Float_t>> tV0_Py;         //!
    std::map<Short_t, std::vector<Float_t>> tV0_Pz;         //!
    std::map<Short_t, std::vector<Float_t>> tV0_E;          //!
    std::map<Short_t, std::vector<Float_t>> tV0_Xv;         //!
    std::map<Short_t, std::vector<Float_t>> tV0_Yv;         //!
    std::map<Short_t, std::vector<Float_t>> tV0_Zv;         //!
    std::map<Short_t, std::vector<Int_t>> tV0_Neg_EsdIdx;   //!
    std::map<Short_t, std::vector<Float_t>> tV0_Neg_Px;     //!
    std::map<Short_t, std::vector<Float_t>> tV0_Neg_Py;     //!
    std::map<Short_t, std::vector<Float_t>> tV0_Neg_Pz;     //!
    std::map<Short_t, std::vector<Int_t>> tV0_Pos_EsdIdx;   //!
    std::map<Short_t, std::vector<Float_t>> tV0_Pos_Px;     //!
    std::map<Short_t, std::vector<Float_t>> tV0_Pos_Py;     //!
    std::map<Short_t, std::vector<Float_t>> tV0_Pos_Pz;     //!
    std::map<Short_t, std::vector<Float_t>> tV0_DCAnegV0;   //!
    std::map<Short_t, std::vector<Float_t>> tV0_DCAposV0;   //!
    std::map<Short_t, std::vector<Float_t>> tV0_DCAbtwDau;  //!
    /* -- V0s true information */
    std::map<Short_t, std::vector<Int_t>> tV0_McIdx;        //!
    std::map<Short_t, std::vector<Int_t>> tV0_PdgCode;      //!
    std::map<Short_t, std::vector<Bool_t>> tV0_IsSignal;    //!
    std::map<Short_t, std::vector<UInt_t>> tV0_ReactionID;  //!
    std::map<Short_t, std::vector<Bool_t>> tV0_IsHybrid;    //!
    /* -- Sexaquarks properties (Channel A) */
    std::vector<Float_t> tTypeA_Px;               //!
    std::vector<Float_t> tTypeA_Py;               //!
    std::vector<Float_t> tTypeA_Pz;               //!
    std::vector<Float_t> tTypeA_E;                //!
    std::vector<Float_t> tTypeA_E_asDecay;        //!
    std::vector<Float_t> tTypeA_Xv;               //!
    std::vector<Float_t> tTypeA_Yv;               //!
    std::vector<Float_t> tTypeA_Zv;               //!
    std::vector<Int_t> tTypeA_V0a_Idx;            //!
    std::vector<Float_t> tTypeA_V0a_Px;           //!
    std::vector<Float_t> tTypeA_V0a_Py;           //!
    std::vector<Float_t> tTypeA_V0a_Pz;           //!
    std::vector<Float_t> tTypeA_V0a_E;            //!
    std::vector<Float_t> tTypeA_V0a_DecayLength;  //!
    std::vector<Float_t> tTypeA_DCAV0aSV;         //!
    std::vector<Float_t> tTypeA_DCAV0aNegSV;      //!
    std::vector<Float_t> tTypeA_DCAV0aPosSV;      //!
    std::vector<Int_t> tTypeA_V0b_Idx;            //!
    std::vector<Float_t> tTypeA_V0b_Px;           //!
    std::vector<Float_t> tTypeA_V0b_Py;           //!
    std::vector<Float_t> tTypeA_V0b_Pz;           //!
    std::vector<Float_t> tTypeA_V0b_E;            //!
    std::vector<Float_t> tTypeA_V0b_DecayLength;  //!
    std::vector<Float_t> tTypeA_DCAV0bSV;         //!
    std::vector<Float_t> tTypeA_DCAV0bNegSV;      //!
    std::vector<Float_t> tTypeA_DCAV0bPosSV;      //!
    std::vector<Float_t> tTypeA_DCAbtwV0s;        //!
    /* -- Sexaquarks true info (Channel A) */
    std::vector<Bool_t> tTypeA_IsSignal;    //!
    std::vector<UInt_t> tTypeA_ReactionID;  //!
    std::vector<Bool_t> tTypeA_IsHybrid;    //!
    /* -- Sexaquarks properties (Channel D) */
    std::vector<Float_t> tTypeD_Px;              //!
    std::vector<Float_t> tTypeD_Py;              //!
    std::vector<Float_t> tTypeD_Pz;              //!
    std::vector<Float_t> tTypeD_E;               //!
    std::vector<Float_t> tTypeD_E_asDecay;       //!
    std::vector<Float_t> tTypeD_Xv;              //! secondary x-vertex
    std::vector<Float_t> tTypeD_Yv;              //! secondary y-vertex
    std::vector<Float_t> tTypeD_Zv;              //! secondary z-vertex
    std::vector<Int_t> tTypeD_V0_Idx;            //!
    std::vector<Float_t> tTypeD_V0_Px;           //!
    std::vector<Float_t> tTypeD_V0_Py;           //!
    std::vector<Float_t> tTypeD_V0_Pz;           //!
    std::vector<Float_t> tTypeD_V0_E;            //!
    std::vector<Float_t> tTypeD_V0_DecayLength;  //!
    std::vector<Float_t> tTypeD_DCAV0SV;         //!
    std::vector<Float_t> tTypeD_DCAV0NegSV;      //!
    std::vector<Float_t> tTypeD_DCAV0PosSV;      //!
    std::vector<Float_t> tTypeD_DCAV0NegKa;      //!
    std::vector<Float_t> tTypeD_DCAV0PosKa;      //!
    std::vector<Int_t> tTypeD_Ka_EsdIdx;         //!
    std::vector<Float_t> tTypeD_Ka_Px;           //!
    std::vector<Float_t> tTypeD_Ka_Py;           //!
    std::vector<Float_t> tTypeD_Ka_Pz;           //!
    std::vector<Float_t> tTypeD_DCAKaSV;         //!
    std::vector<Float_t> tTypeD_DCAKaV0;         //!
    /* -- Sexaquarks true info (Channel D) */
    std::vector<Bool_t> tTypeD_IsSignal;    //!
    std::vector<UInt_t> tTypeD_ReactionID;  //!
    std::vector<Bool_t> tTypeD_IsHybrid;    //!

    /* Containers */
    /* -- key: `mc_idx` */
    std::unordered_map<Int_t, Int_t> fMC_PdgCode_;       //!
    std::unordered_map<Int_t, Int_t> fMC_Mother_McIdx_;  //!
    std::unordered_map<Int_t, Bool_t> fMC_IsSignal_;     //!
    std::unordered_map<Int_t, UInt_t> fMC_ReactionID_;   //!
    /* -- key: `reaction_id` */
    std::unordered_map<UInt_t, std::vector<Int_t>> fReactionProducts_McIdx_;  //!
    /* -- key: `esdIdx` */
    std::unordered_map<Int_t, Int_t> fLinked_McIdx_;  //!
    /* -- looped over in `KalmanV0Finder()` and Sexaquark Finders */
    std::vector<Int_t> fAntiProton_EsdIdx;  //!
    std::vector<Int_t> fProton_EsdIdx;      //!
    std::vector<Int_t> fNegKaon_EsdIdx;     //!
    std::vector<Int_t> fPosKaon_EsdIdx;     //!
    std::vector<Int_t> fPiMinus_EsdIdx;     //!
    std::vector<Int_t> fPiPlus_EsdIdx;      //!
    /*  */
    std::vector<KFParticle> kfAntiLambdas;     //!
    std::vector<KFParticle> kfKaonsZeroShort;  //!

    AliAnalysisTaskSexaquark(const AliAnalysisTaskSexaquark&);             // not implemented
    AliAnalysisTaskSexaquark& operator=(const AliAnalysisTaskSexaquark&);  // not implemented

    /// \cond CLASSDEF
    ClassDef(AliAnalysisTaskSexaquark, 78);  // = number of persistent members
    /// \endcond
};

#endif
