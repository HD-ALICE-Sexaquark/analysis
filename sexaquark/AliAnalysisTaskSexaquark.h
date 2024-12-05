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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Math/Point3D.h"   // because "TVector3.h" is deprecated
#include "Math/Vector3D.h"  // because "TVector3.h" is deprecated
#include "Math/Vector4D.h"  // because "TLorentzVector.h" is deprecated
#include "Math/VectorUtil.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

using namespace ROOT;

class AliPIDResponse;
class KFParticle;
class KFVertex;

/*
 * Auxiliary class to make use of protected function `KFParticleBase::GetMeasurement()`
 * (Copied from `/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.h`)
 */
class KFParticleMother : public KFParticle {
   public:
    Bool_t CheckDaughter(KFParticle daughter) {
        Float_t m[8], mV[36], D[3][3];
        return KFParticleBase::GetMeasurement(daughter, m, mV, D);
    }
};

class AliAnalysisTaskSexaquark : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskSexaquark();
    AliAnalysisTaskSexaquark(const char* name);
    virtual ~AliAnalysisTaskSexaquark();
    virtual void Terminate(Option_t* option) { return; }

    /* Settings ~ stored in Analysis Manager */

    void IsMC(Bool_t IsMC, Bool_t IsSignalMC = kFALSE) {
        fIsMC = IsMC;
        fIsSignalMC = IsSignalMC;
    };
    void SetSourceOfV0s(TString SourceOfV0s) { fSourceOfV0s = SourceOfV0s; };
    void DoQA(Bool_t DoQA) { fDoQA = DoQA; };
    void Initialize();

    /* Main ~ executed at runtime */

    virtual void UserCreateOutputObjects();
    virtual Bool_t UserNotify();
    virtual void UserExec(Option_t* option);
    void EndOfEvent();

    /* Trees */

    void PrepareEventsBranches();
    void PrepareInjectedBranches();
    void PrepareMCParticlesBranches();
    void PrepareTracksBranches();
    void PrepareV0sBranches();
    void PrepareCommonSexaquarkBranches(TTree* Tree_Sexaquarks);
    void PrepareSexaquarksBranches();

    void ClearEventsBranches();
    void ClearInjectedBranches();
    void ClearMCParticlesBranches();
    void ClearTracksBranches();
    void ClearV0sBranches();
    void ClearSexaquarksBranches();

    /* Histograms */

    void PrepareQAHistograms();
    void PrepareTracksHistograms();
    void PrepareV0sHistograms();
    void PrepareSexaquarksHistograms();

    /* Cuts */

    void DefineTracksCuts(TString cuts_option);
    void DefineV0Cuts(TString cuts_option);
    void DefineSexaquarkCuts(TString cuts_option);

    /* MC Generated */

    void ProcessMCGen();
    void ProcessSignalReactions();

    /* Tracks */

    Bool_t PassesEventSelection();
    void ProcessTracks();
    Bool_t PassesTrackSelection(AliESDtrack* track, Bool_t isSecondary, Bool_t isSignal);
    Bool_t PassesSpeciesSelection(AliESDtrack* track, Int_t esdPdgCode, Bool_t isSecondary, Bool_t isSignal);

    /* V0s, K+K+ pairs, Anti-Sexaquarks -- That Require True Info */

    void ProcessFindableV0s();
    void ProcessFindableSexaquarks();

    /* V0s -- ALICE V0 Finders */

    void GetV0sFromESD(Bool_t onTheFly);
    void CustomV0Finder(Int_t pdgV0);
    Bool_t PassesV0CutsAs(Int_t pdgV0, AliESDv0* esdV0, AliESDtrack* esdNegDaughter, AliESDtrack* esdPosDaughter, Bool_t isTrue, Bool_t isSecondary,
                          Bool_t isSignal);

    /* Sexaquark -- Geometrical Finder */

    /** `AntiSexaquark Neutron -> AntiLambda KaonZeroShort` **/
    /** `Sexaquark AntiNeutron -> Lambda KaonZeroShort` **/

    void GeoSexaquarkFinder_TypeA(TString ReactionChannel);
    Bool_t PassesSexaquarkCuts_TypeA(TString ReactionChannel, Math::XYZPoint SecondaryVertex, Math::PxPyPzEVector lvSexaquark, AliESDv0 esdLambda,
                                     AliESDv0 esdKaonZeroShort, AliESDtrack* esdLambdaNegDau, AliESDtrack* esdLambdaPosDau,
                                     AliESDtrack* esdKaonZeroShortNegDau, AliESDtrack* esdKaonZeroShortPosDau, Bool_t isSignal);

    /* V0s -- ALICE V0 Finders */

    void KalmanV0Finder(Int_t pdgNegDaughter, Int_t pdgPosDaughter, Int_t pdgV0 = -1);
    Bool_t PassesV0CutsAs(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfNegDaughter, KFParticle kfPosDaughter, Math::PxPyPzEVector lvV0,
                          Math::PxPyPzEVector lvNegDaughter, Math::PxPyPzEVector lvPosDaughter, Bool_t isTrue, Bool_t isSecondary, Bool_t isSignal);
    void StoreV0As(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfNegDaughter, KFParticle kfPosDaughter, Math::PxPyPzEVector lvV0,
                   Math::PxPyPzEVector lvNegDaughter, Math::PxPyPzEVector lvPosDaughter, Int_t esdIdxNegDau, Int_t esdIdxPosDau, Int_t mcIdxV0,
                   Int_t mcPdgCode, Bool_t isTrue, Bool_t isSecondary, Bool_t isSignal, Int_t ReactionID, Bool_t isHybrid);

    /* Sexaquark -- Kalman Filter */

    /** `AntiSexaquark Neutron -> AntiLambda KaonZeroShort` **/
    /** `Sexaquark AntiNeutron -> Lambda KaonZeroShort` **/

    void KalmanSexaquarkFinder_TypeA(TString ReactionChannel);
    Bool_t PassesSexaquarkCuts_TypeA(TString ReactionChannel, KFParticleMother kfSexaquark, Math::PxPyPzEVector lvSexaquark, KFParticle kfLambda,
                                     KFParticle kfKaonZeroShort, KFParticle kfLambdaNegDau, KFParticle kfLambdaPosDau,
                                     KFParticle kfKaonZeroShortNegDau, KFParticle kfKaonZeroShortPosDau, Bool_t isSignal);
    void StoreSexaquark_TypeA(TString ReactionChannel, KFParticleMother kfSexaquark, KFParticleMother kfKaonZeroShort, KFParticleMother kfLambda,
                              KFParticle kfKaonZeroShortNegDau, KFParticle kfKaonZeroShortPosDau, KFParticle kfLambdaNegDau,
                              KFParticle kfLambdaPosDau, Math::PxPyPzEVector lvSexaquark, Math::PxPyPzEVector lvKaonZeroShort,
                              Math::PxPyPzEVector lvLambda, Int_t idxKaonZeroShort, Int_t esdIdxNeg_KaonZeroShort, Int_t esdIdxPos_KaonZeroShort,
                              Int_t idxLambda, Int_t esdIdxNeg_Lambda, Int_t esdIdxPos_Lambda, Bool_t isSignal, Int_t ReactionID, Bool_t isHybrid);

    /** `AntiSexaquark Proton -> AntiLambda K+` **/
    /** `AntiSexaquark Proton -> AntiLambda K+ pi- pi+` **/
    /** `Sexaquark AntiProton -> Lambda K-` **/
    /** `Sexaquark AntiProton -> Lambda K- pi- pi+` **/

    void KalmanSexaquarkFinder_TypeDE(TString ReactionChannel);
    Bool_t PassesSexaquarkCuts_TypeD(TString ReactionChannel, KFParticleMother kfSexaquark, Math::PxPyPzEVector lvSexaquark, KFParticle kfLambda,
                                     KFParticle kfKaon, KFParticle kfLambdaNegDau, KFParticle kfLambdaPosDau, Bool_t isSignal);
    Bool_t PassesSexaquarkCuts_TypeE(TString ReactionChannel, KFParticleMother kfSexaquark, Math::PxPyPzEVector lvSexaquark, KFParticle kfLambda,
                                     KFParticle kfLambdaNegDau, KFParticle kfLambdaPosDau, KFParticle kfKaon, KFParticle kfPiMinus,
                                     KFParticle kfPiPlus, Bool_t isSignal);
    void StoreSexaquark_TypeD(TString ReactionChannel, KFParticleMother kfSexaquark, KFParticleMother kfLambda, KFParticle kfLambdaNegDau,
                              KFParticle kfLambdaPosDau, KFParticle kfKaon, Math::PxPyPzEVector lvSexaquark, Math::PxPyPzEVector lvLambda,
                              Math::PxPyPzEVector lvKaon, Int_t idxLambda, Int_t esdIdxNeg_Lambda, Int_t esdIdxPos_Lambda, Int_t esdIdxKaon,
                              Bool_t isSignal, Int_t ReactionID, Bool_t isHybrid);
    void StoreSexaquark_TypeE(TString ReactionChannel, KFParticleMother kfSexaquark, KFParticleMother kfLambda, KFParticle kfLambdaNegDau,
                              KFParticle kfLambdaPosDau, KFParticle kfKaon, KFParticle kfPiMinus, KFParticle kfPiPlus,
                              Math::PxPyPzEVector lvSexaquark, Math::PxPyPzEVector lvLambda, Math::PxPyPzEVector lvKaon,
                              Math::PxPyPzEVector lvPiMinus, Math::PxPyPzEVector lvPiPlus, Int_t idxLambda, Int_t esdIdxNeg_Lambda,
                              Int_t esdIdxPos_Lambda, Int_t esdIdxKaon, Int_t idxPionPair, Int_t esdIdxPiMinus, Int_t esdIdxPiPlus, Bool_t isSignal,
                              Int_t ReactionID, Bool_t isHybrid);

    /** `AntiSexaquark Proton -> K+ K+ X` **/
    /** `Sexaquark AntiProton -> K- K- X` **/

    void KalmanSexaquarkFinder_TypeH(TString ReactionChannel);
    Bool_t PassesSexaquarkCuts_TypeH(TString ReactionChannel, KFParticleMother kfKaonPair, KFParticle kfKaonA, KFParticle kfKaonB,
                                     Math::PxPyPzEVector lvKaonPair, Math::PxPyPzEVector lvKaonA, Math::PxPyPzEVector lvKaonB, Bool_t isSignal);
    void StoreSexaquark_TypeH(TString ReactionChannel, KFParticleMother kfKaonPair, KFParticle kfKaonA, KFParticle kfKaonB,
                              Math::PxPyPzEVector lvKaonPair, Math::PxPyPzEVector lvKaonA, Math::PxPyPzEVector lvKaonB, Int_t esdIdxKaonA,
                              Int_t esdIdxKaonB, Bool_t isSignal, Int_t ReactionID, Bool_t isHybrid);

    /* Utilities */
    void ClearContainers();
    Bool_t Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn, Double_t* lPreprocessxp,
                       const Double_t b);
    void GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b);
    void Evaluate(const Double_t* h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]);
    Float_t MCRec_GetImpactParameter(AliESDtrack* track);
    Double_t CosinePointingAngle(Math::PxPyPzEVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);
    Double_t CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);
    Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                             Double_t Pos_Py, Double_t Pos_Pz);
    Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz);
    Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z, Double_t refPointX,
                          Double_t refPointY, Double_t refPointZ);
    Double_t SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]);
    Double_t SquaredDistanceBetweenLines(const Double_t* t, Double_t pos0[], Double_t dir0[], Double_t pos1[], Double_t dir1[]);
    Double_t Calculate_TwoLinesDCA_v1(Math::XYZPoint v3_pos0, Math::XYZVector v3_dir0, Math::XYZPoint v3_pos1, Math::XYZVector v3_dir1,
                                      Math::XYZPoint& PCA0, Math::XYZPoint& PCA1);
    Double_t Calculate_TwoLinesDCA_v2(Math::XYZPoint v3_pos0, Math::XYZVector v3_dir0, Math::XYZPoint v3_pos1, Math::XYZVector v3_dir1,
                                      Math::XYZPoint& PCA0, Math::XYZPoint& PCA1);

    /* Kalman Filter Utilities */
    KFParticle CreateKFParticle(AliESDtrack& track, Double_t mass, Int_t charge, Bool_t inner = kTRUE);
    KFVertex CreateKFVertex(const AliVVertex& vertex);
    KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Double_t massThis, Int_t chargeThis);
    void GetDStoParticleBz(float Bz, const KFParticleBase& p, const KFParticleBase& q, float dS[2], float dsdr[4][6], const float* param1 = 0,
                           const float* param2 = 0) const;

    /* External Files */
    Bool_t ReadSignalLogs();

   private:
    /* Settings ~ stored in Analysis Manager ~ all persistent */
    Bool_t fIsMC;          // kTRUE if MC simulation, kFALSE if data
    Bool_t fIsSignalMC;    // kTRUE to read and load signal logs
    TString fSourceOfV0s;  // choose V0 finder: "kalman", "custom", "on-the-fly", "offline"
    Bool_t fDoQA;          // kTRUE to do QA

    /* AliRoot Objects */
    AliMCEvent* fMC;                    //! MC event
    AliVVertex* fMC_PrimaryVertex;      //! MC gen. (or true) primary vertex
    Math::XYZPoint v3MC_PrimaryVertex;  //!
    AliESDEvent* fESD;                  //! reconstructed event
    AliESDVertex* fPrimaryVertex;       //! primary vertex
    KFVertex kfPrimaryVertex;           //!
    Math::XYZPoint v3PrimaryVertex;     //!
    AliPIDResponse* fPIDResponse;       //! pid response object
    AliEventCuts fEventCuts;            //! event cuts
    Double_t fMagneticField;            //! magnetic field
    Int_t fRunNumber;                   //! run number
    Float_t fDirNumber;                 //! directory number
    Int_t fEventNumber;                 //! event number
    Float_t fCentrality;                //! centrality percentile

    /* ROOT Objects */
    TDatabasePDG fPDG;                    //!
    TList* fList_QA_Hists;                //!
    TList* fList_AntiProton_Hists;        //!
    TList* fList_Proton_Hists;            //!
    TList* fList_NegKaon_Hists;           //!
    TList* fList_PosKaon_Hists;           //!
    TList* fList_PiMinus_Hists;           //!
    TList* fList_PiPlus_Hists;            //!
    TList* fList_AntiLambda_Hists;        //!
    TList* fList_Lambda_Hists;            //!
    TList* fList_KaonZeroShort_Hists;     //!
    TList* fList_PionPair_Hists;          //!
    TList* fList_Sexaquark_ALK0_Hists;    //!
    TList* fList_Sexaquark_ALPK_Hists;    //!
    TList* fList_Sexaquark_ALPKPP_Hists;  //!
    TList* fList_Sexaquark_PKPKX_Hists;   //!
    TList* fList_Sexaquark_LK0_Hists;     //!
    TList* fList_Sexaquark_LNK_Hists;     //!
    TList* fList_Sexaquark_LNKPP_Hists;   //!
    TList* fList_Sexaquark_NKNKX_Hists;   //!

    /* External Files */
    TString fAliEnPath;        //! loaded in `UserNotify()`
    Bool_t fIsFirstEvent;      //!
    TString fReactionChannel;  //! derived from `fAliEnPath` in `UserNotify()`

    /* Filled at `Initialize()` ~ persistent */
    Double_t kMass_Neutron;                                        //
    Double_t kMass_Proton;                                         //
    Double_t kMass_Kaon;                                           //
    Double_t kMass_Pion;                                           //
    std::unordered_map<Int_t, Int_t> getNegPdgCode_fromV0PdgCode;  // key: `V0 pdg code`, value: `negative daughter pdg code`
    std::unordered_map<Int_t, Int_t> getPosPdgCode_fromV0PdgCode;  // key: `V0 pdg code`, value: `positive daughter pdg code`

    /*** Trees ***/

    TTree* fTree_Events;             //!
    TTree* fTree_Injected;           //!
    TTree* fTree_MC;                 //!
    TTree* fTree_Tracks;             //!
    TTree* fTree_V0s;                //!
    TTree* fTree_Sexaquarks_ALK0;    //!
    TTree* fTree_Sexaquarks_ALPK;    //!
    TTree* fTree_Sexaquarks_ALPKPP;  //!
    TTree* fTree_Sexaquarks_PKPKX;   //!
    TTree* fTree_Sexaquarks_LK0;     //!
    TTree* fTree_Sexaquarks_LNK;     //!
    TTree* fTree_Sexaquarks_LNKPP;   //!
    TTree* fTree_Sexaquarks_NKNKX;   //!

    /*** Branches ***/

    /** Events **/
    Float_t tEvents_PV_TrueXv;            //!
    Float_t tEvents_PV_TrueYv;            //!
    Float_t tEvents_PV_TrueZv;            //!
    Bool_t tEvents_IsGenPileup;           //!
    Bool_t tEvents_IsSBCPileup;           //!
    Float_t tEvents_PV_RecXv;             //!
    Float_t tEvents_PV_RecYv;             //!
    Float_t tEvents_PV_RecZv;             //!
    Int_t tEvents_PV_NContributors;       //!
    Float_t tEvents_PV_ZvErr_FromSPD;     //!
    Float_t tEvents_PV_ZvErr_FromTracks;  //!
    Float_t tEvents_PV_Zv_FromSPD;        //!
    Float_t tEvents_PV_Zv_FromTracks;     //!
    Float_t tEvents_PV_Dispersion;        //!
    Int_t tEvents_NTracks;                //!
    Int_t tEvents_NTPCClusters;           //!
    Bool_t tEvents_IsMB;                  //!
    Bool_t tEvents_IsHighMultV0;          //!
    Bool_t tEvents_IsHighMultSPD;         //!
    Bool_t tEvents_IsCentral;             //!
    Bool_t tEvents_IsSemiCentral;         //!

    /** Injected Sexaquarks **/
    Int_t tInjected_RunNumber;         //!
    Int_t tInjected_DirNumber;         //!
    Int_t tInjected_EventNumber;       //!
    Int_t tInjected_ReactionID;        //!
    Float_t tInjected_Px;              //!
    Float_t tInjected_Py;              //!
    Float_t tInjected_Pz;              //!
    Float_t tInjected_Mass;            //!
    Int_t tInjected_Nucleon_PdgCode;   //!
    Float_t tInjected_Nucleon_Px;      //!
    Float_t tInjected_Nucleon_Py;      //!
    Float_t tInjected_Nucleon_Pz;      //!
    UInt_t tInjected_ReactionChannel;  //! char -> ASCII -> uint

    /** MC Particles **/
    Int_t tMC_Idx;           //!
    Int_t tMC_PdgCode;       //!
    Int_t tMC_Idx_Mother;    //!
    Int_t tMC_NDaughters;    //!
    Int_t tMC_Idx_FirstDau;  //!
    Int_t tMC_Idx_LastDau;   //!
    Float_t tMC_Px;          //!
    Float_t tMC_Py;          //!
    Float_t tMC_Pz;          //!
    Float_t tMC_Xv;          //! origin x-vertex
    Float_t tMC_Yv;          //! origin y-vertex
    Float_t tMC_Zv;          //! origin z-vertex
    UInt_t tMC_Status;       //!
    Bool_t tMC_IsOOBPileup;  //!
    Bool_t tMC_IsSecondary;  //!
    Bool_t tMC_IsSignal;     //!
    Int_t tMC_ReactionID;    //!

    /** Tracks **/
    Int_t tTrack_Idx;              //!
    Float_t tTrack_Px;             //! inner parametrization
    Float_t tTrack_Py;             //! inner parametrization
    Float_t tTrack_Pz;             //! inner parametrization
    Short_t tTrack_Charge;         //!
    Float_t tTrack_NSigmaPion;     //!
    Float_t tTrack_NSigmaKaon;     //!
    Float_t tTrack_NSigmaProton;   //!
    TBits tTrack_TPCFitMap;        //!
    TBits tTrack_TPCClusterMap;    //!
    TBits tTrack_TPCSharedMap;     //!
    Bool_t tTrack_IsKinkDaughter;  //!
    Bool_t tTrack_ITSin;           //!
    Bool_t tTrack_ITSout;          //!
    Bool_t tTrack_ITSrefit;        //!
    Bool_t tTrack_ITSpid;          //!
    Bool_t tTrack_TPCin;           //!
    Bool_t tTrack_TPCout;          //!
    Bool_t tTrack_TPCrefit;        //!
    Bool_t tTrack_TPCpid;          //!
    Int_t tTrack_Idx_True;         //!
    Int_t tTrack_True_PdgCode;     //!
    Bool_t tTrack_IsSecondary;     //!
    Bool_t tTrack_IsSignal;        //!
    Int_t tTrack_ReactionID;       //!

    /** V0s **/
    Int_t tV0_Idx;           //!
    Int_t tV0_Idx_Neg;       //!
    Int_t tV0_Idx_Pos;       //!
    Int_t tV0_PID;           //!
    Float_t tV0_Px;          //!
    Float_t tV0_Py;          //!
    Float_t tV0_Pz;          //!
    Float_t tV0_E;           //!
    Float_t tV0_Xv;          //! V0 x-vertex
    Float_t tV0_Yv;          //! V0 y-vertex
    Float_t tV0_Zv;          //! V0 z-vertex
    Float_t tV0_Neg_Px;      //!
    Float_t tV0_Neg_Py;      //!
    Float_t tV0_Neg_Pz;      //!
    Float_t tV0_Pos_Px;      //!
    Float_t tV0_Pos_Py;      //!
    Float_t tV0_Pos_Pz;      //!
    Int_t tV0_Idx_True;      //!
    Int_t tV0_True_PdgCode;  //!
    Bool_t tV0_IsSecondary;  //!
    Bool_t tV0_IsSignal;     //!
    Int_t tV0_ReactionID;    //!
    Bool_t tV0_IsHybrid;     //!

    /** (Anti-)Sexaquark Candidates **/
    Float_t tSexaquark_Px;          //!
    Float_t tSexaquark_Py;          //!
    Float_t tSexaquark_Pz;          //!
    Float_t tSexaquark_E;           //!
    Float_t tSexaquark_E_asDecay;   //!
    Float_t tSexaquark_Xv;          //! secondary x-vertex
    Float_t tSexaquark_Yv;          //! secondary y-vertex
    Float_t tSexaquark_Zv;          //! secondary z-vertex
    Float_t tSexaquark_DistFromPV;  //!
    Float_t tSexaquark_CPAwrtPV;    //!
    Float_t tSexaquark_DCAwrtPV;    //!
    Float_t tSexaquark_Chi2ndf;     //!
    /* True Information */
    Bool_t tSexaquark_IsSignal;   //!
    Int_t tSexaquark_ReactionID;  //!
    Bool_t tSexaquark_IsHybrid;   //!
    /* Channels "A"+"D"+"E" */
    Int_t tSexaquark_Idx_Lambda;            //!
    Int_t tSexaquark_Idx_Lambda_Neg;        //!
    Int_t tSexaquark_Idx_Lambda_Pos;        //!
    Float_t tSexaquark_Lambda_DecayLength;  //!
    Float_t tSexaquark_DCALaSV;             //!
    Float_t tSexaquark_DCALaNegSV;          //!
    Float_t tSexaquark_DCALaPosSV;          //!
    /* Channels "A"+"D"+"H" */
    Float_t tSexaquark_OpeningAngle;  //!
    /* Channel "A" */
    Int_t tSexaquarkA_Idx_K0S;            //!
    Int_t tSexaquarkA_Idx_K0S_Neg;        //!
    Int_t tSexaquarkA_Idx_K0S_Pos;        //!
    Float_t tSexaquarkA_K0S_DecayLength;  //!
    Float_t tSexaquarkA_DCAK0SV;          //!
    Float_t tSexaquarkA_DCAK0NegSV;       //!
    Float_t tSexaquarkA_DCAK0PosSV;       //!
    Float_t tSexaquarkA_DCAbtwV0s;        //!
    /* Channel "D" */
    Int_t tSexaquarkD_Idx_Kaon;      //!
    Float_t tSexaquarkD_DCAKaSV;     //!
    Float_t tSexaquarkD_DCAKaLa;     //!
    Float_t tSexaquarkD_DCALaNegKa;  //!
    Float_t tSexaquarkD_DCALaPosKa;  //!
    /* Channel "E" */
    Int_t tSexaquarkE_Idx_Kaon;     //!
    Int_t tSexaquarkE_Idx_PP;       //!
    Int_t tSexaquarkE_Idx_PiMinus;  //!
    Int_t tSexaquarkE_Idx_PiPlus;   //!
    Float_t tSexaquarkE_DCAKaSV;    //!
    Float_t tSexaquarkE_DCAKaLa;    //!
    Float_t tSexaquarkE_DCApmSV;    //!
    Float_t tSexaquarkE_DCAppSV;    //!
    Float_t tSexaquarkE_DCApmLa;    //!
    Float_t tSexaquarkE_DCAppLa;    //!
    Float_t tSexaquarkE_DCApmKa;    //!
    Float_t tSexaquarkE_DCAppKa;    //!
    /* Channel "H" */
    Int_t tKaonPair_Idx_KaonA;   //!
    Int_t tKaonPair_Idx_KaonB;   //!
    Float_t tKaonPair_DCAbtwKK;  //!
    Float_t tKaonPair_DCAkaSV;   //!
    Float_t tKaonPair_DCAkbSV;   //!

    /* Histograms */

    /** QA Histograms **/

    /*
    fHist_QA["Tracks_Bookkeep"]:
     - `0`    : MC tracks
     - `10`   : SECONDARY MC tracks
     - `20`   : SIGNAL MC tracks
     - `30+`  : effect of cuts on found
     - `50`   : found tracks that passed all cuts
     - `90+`  : effect of cuts on SECONDARY tracks
     - `110`  : found SECONDARY tracks that passed all cuts
     - `120+` : effect of cuts on SIGNAL tracks
     - `140`  : found SIGNAL tracks that passed all cuts
    */
    std::map<TString, TH1F*> fHist_QA;    //! key: `property`
    std::map<TString, TH2F*> f2DHist_QA;  //! key: `hist name`

    /** Tracks Histograms **/

    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_Tracks;  //! key: `stage, set, pdg code, property`
    std::map<std::tuple<TString, Int_t>, TH2F*> f2DHist_Tracks_TPCsignal;        //! key: `set, pdg code`

    /** V0s Histograms **/

    /*
     - `0+`   : MC V0s           ( 1: relevant channels)
     - `10+`  : SECONDARY MC V0s (11: relevant channels)
     - `20+`  : SIGNAL MC V0s    (21: relevant channels)
     - `30`   : findable true V0s
     - `40`   : findable true SECONDARY V0s
     - `50`   : findable true SIGNAL V0s
     - `60+`  : effect of cuts on found
     - `80`   : found V0s that passed all cuts
     - `90+`  : effect of cuts on TRUE V0s
     - `110`  : found TRUE V0s that passed all cuts
     - `120+` : effect of cuts on SECONDARY V0s
     - `140`  : found SECONDARY V0s that passed all cuts
     - `150+` : effect of cuts on SIGNAL V0s
     - `170`  : found SIGNAL V0s that passed all cuts
    */
    std::unordered_map<Int_t, TH1F*> fHist_V0s_Bookkeep;                           //! key: `V0 pdg code`
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_V0s;       //! key: `stage, set, V0 pdg code, property`
    std::map<std::tuple<TString, Int_t>, TH2F*> f2DHist_V0s_ArmenterosPodolanski;  //! key: `set, V0 pdg code`

    /** (Anti-)Sexaquark Histograms **/

    /*
     - `0`   : MC Gen.
     - `10`  : findable
     - `20+` : effect of cuts on found candidates
     - `50`  : found candidates that passed all cuts
     - `60+` : effect of cuts on SIGNAL candidates
     - `90`  : found SIGNAL candidates that passed all cuts
    */
    std::map<TString, TH1D*> fHist_Sexaquarks_Bookkeep;                                //! key: `reaction channel`
    std::map<std::tuple<TString, TString, TString, TString>, TH1D*> fHist_Sexaquarks;  //! key: `reaction channel, stage, set, property`

    /*** Containers -- vectors and hash tables ***/
    /* Note: when adding a new one, don't forget to clear it on `ClearContainers()` */

    /* filled at `ProcessMCGen()` */
    std::unordered_map<Int_t, Int_t> getPdgCode_fromMcIdx;  //! key: `mcIdx`

    std::unordered_map<Int_t, Bool_t> isMcIdxSignal;     //! key: `mcIdx`
    std::unordered_map<Int_t, Bool_t> isMcIdxSecondary;  //! key: `mcIdx`

    std::unordered_map<Int_t, Int_t> getReactionID_fromMcIdx;                   //! key: `mcIdx`
    std::unordered_map<Int_t, std::vector<Int_t>> getMcIndices_fromReactionID;  //! key: `ReactionID`

    std::unordered_map<Int_t, Bool_t> doesMcIdxHaveMother;      //! key: `mcIdx`
    std::unordered_map<Int_t, Int_t> getMotherMcIdx_fromMcIdx;  //! key: `mcIdx`
    std::unordered_map<Int_t, Int_t> getNegDauMcIdx_fromMcIdx;  //! key: `mcIdx`
    std::unordered_map<Int_t, Int_t> getPosDauMcIdx_fromMcIdx;  //! key: `mcIdx`

    /* filled at `ProcessTracks()` */
    std::unordered_map<Int_t, Int_t> getMcIdx_fromEsdIdx;  //! key: `esdIdx`

    /* used in `GetV0sFromESD()` */
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsProton;      //! key: `esdIdx`
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsAntiProton;  //! key: `esdIdx`
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsPosKaon;     //! key: `esdIdx`
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsNegKaon;     //! key: `esdIdx`
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsPiPlus;      //! key: `esdIdx`
    std::unordered_map<Int_t, Bool_t> isEsdIdxSelectedAsPiMinus;     //! key: `esdIdx`

    /* used in `ProcessFindableSexaquarks()` */
    std::unordered_map<Int_t, std::vector<Int_t>> getEsdIndices_fromReactionID;  //! key: `ReactionID`

    /* used in `ProcessFindableV0s()` */
    std::vector<Int_t> mcIndicesOfTrueV0s;                       //!
    std::unordered_map<Int_t, Int_t> getNegDauEsdIdx_fromMcIdx;  //! key: `mcIdx`
    std::unordered_map<Int_t, Int_t> getPosDauEsdIdx_fromMcIdx;  //! key: `mcIdx`

    /* looped over in `KalmanV0Finder()` and Sexaquark Finders */
    std::vector<Int_t> esdIndicesOfAntiProtonTracks;  //!
    std::vector<Int_t> esdIndicesOfProtonTracks;      //!
    std::vector<Int_t> esdIndicesOfNegKaonTracks;     //!
    std::vector<Int_t> esdIndicesOfPosKaonTracks;     //!
    std::vector<Int_t> esdIndicesOfPiMinusTracks;     //!
    std::vector<Int_t> esdIndicesOfPiPlusTracks;      //!

    /* filled at `GetV0sFromESD()` and `[Custom|Kalman]V0Finder()` */
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromLambdaIdx;         //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromLambdaIdx;         //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromAntiLambdaIdx;     //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromAntiLambdaIdx;     //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromKaonZeroShortIdx;  //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromKaonZeroShortIdx;  //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromPionPairIdx;       //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromPionPairIdx;       //!

    std::vector<AliESDv0> esdLambdas;         //!
    std::vector<AliESDv0> esdAntiLambdas;     //!
    std::vector<AliESDv0> esdKaonsZeroShort;  //!

    std::vector<KFParticleMother> kfLambdas;         //!
    std::vector<KFParticleMother> kfAntiLambdas;     //!
    std::vector<KFParticleMother> kfKaonsZeroShort;  //!
    std::vector<KFParticleMother> kfPionPairs;       //!

    /*** Cuts ~ persistent, because they are set on `Initialize()` ***/

    /** Tracks **/
    Float_t kMax_NSigma_Pion;                          //
    Float_t kMax_NSigma_Kaon;                          //
    Float_t kMax_NSigma_Proton;                        //
    Float_t kMax_Track_Eta;                            //
    Float_t kMin_Track_NTPCClusters;                   //
    Float_t kMax_Track_Chi2PerNTPCClusters;            //
    Bool_t kTurnedOn_Track_StatusCuts;                 //
    Bool_t kTurnedOn_Track_RejectKinks;                //
    Float_t kMin_Track_DCA_wrtPV;                      //
    Float_t kMin_Track_DCAxy_wrtPV;                    //
    Float_t kMin_Track_DCAz_wrtPV;                     //
    std::unordered_map<Int_t, Float_t> kMin_Track_Pt;  // key: `pdg code`

    /** V0s **/
    /** Note: the key of these maps has to be a positive value **/
    std::unordered_map<UInt_t, Float_t> kMin_V0_Mass;            // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_Mass;            // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_Eta;             // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_ArmPtOverAlpha;  // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMin_V0_Pt;              // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMin_V0_Radius;          // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMin_V0_DistFromPV;      // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_DistFromPV;      // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMin_V0_CPAwrtPV;        // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_CPAwrtPV;        // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMin_V0_DCAwrtPV;        // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_DCAbtwDau;       // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_DCAnegV0;        // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_DCAposV0;        // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_ImprvDCAbtwDau;  // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_Chi2;            // key: `pdg code`
    std::unordered_map<UInt_t, Float_t> kMax_V0_Chi2ndf;         // key: `pdg code`

    /** (Anti-)Sexaquark Candidates **/
    Float_t kMin_Sexaquark_Pt;          //
    Float_t kMax_Sexaquark_Eta;         //
    Float_t kMax_Sexaquark_Rapidity;    //
    Float_t kMin_Sexaquark_DistFromPV;  //
    Float_t kMin_Sexaquark_Radius;      //
    Float_t kMin_Sexaquark_CPAwrtPV;    //
    Float_t kMax_Sexaquark_CPAwrtPV;    //
    Float_t kMin_Sexaquark_DCAwrtPV;    //
    Float_t kMax_Sexaquark_DCAwrtPV;    //
    Float_t kMax_Sexaquark_Chi2ndf;     //
    /* Channels "A"+"D"+"E" */
    Float_t kMax_Sexaquark_DCALaSV;     //
    Float_t kMax_Sexaquark_DCALaNegSV;  //
    Float_t kMax_Sexaquark_DCALaPosSV;  //
    /* Channel "A" */
    Float_t kMax_SexaquarkA_DCAK0SV;     //
    Float_t kMax_SexaquarkA_DCAK0NegSV;  //
    Float_t kMax_SexaquarkA_DCAK0PosSV;  //
    Float_t kMax_SexaquarkA_DCAbtwV0s;   //
    /* Channel "D" */
    Float_t kMax_SexaquarkD_DCAKaSV;     //
    Float_t kMax_SexaquarkD_DCAKaLa;     //
    Float_t kMax_SexaquarkD_DCALaNegKa;  //
    Float_t kMax_SexaquarkD_DCALaPosKa;  //
    /* Channel "E" */
    Float_t kMax_SexaquarkE_DCAKaSV;  //
    Float_t kMax_SexaquarkE_DCAKaLa;  //
    Float_t kMax_SexaquarkE_DCApmSV;  //
    Float_t kMax_SexaquarkE_DCAppSV;  //
    Float_t kMax_SexaquarkE_DCApmLa;  //
    Float_t kMax_SexaquarkE_DCAppLa;  //
    Float_t kMax_SexaquarkE_DCApmKa;  //
    Float_t kMax_SexaquarkE_DCAppKa;  //
    /* Channel "H" */
    Float_t kMax_KaonPair_DCAbtwKK;  //
    Float_t kMax_KaonPair_DCAkaSV;   //
    Float_t kMax_KaonPair_DCAkbSV;   //

    AliAnalysisTaskSexaquark(const AliAnalysisTaskSexaquark&);             // not implemented
    AliAnalysisTaskSexaquark& operator=(const AliAnalysisTaskSexaquark&);  // not implemented

    /// \cond CLASSDEF
    ClassDef(AliAnalysisTaskSexaquark, 78);  // = number of persistent members
    /// \endcond
};

#endif
