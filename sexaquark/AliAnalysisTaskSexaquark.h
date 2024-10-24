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
#include <vector>

#include "TArray.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TGrid.h"
#include "TH1.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

class AliPIDResponse;
class KFParticle;
class KFVertex;

/*
 Auxiliary class to make use of protected function KFParticleBase::GetMeasurement()
 (Copied from `/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.h`)
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
    void IsMC(Bool_t IsMC) { fIsMC = IsMC; };
    void SetSourceOfV0s(TString SourceOfV0s) { fSourceOfV0s = SourceOfV0s; };
    void SetReactionID(Char_t ReactionID) { fReactionID = ReactionID; };
    void DoQA(Bool_t DoQA) { fDoQA = DoQA; };
    void StopAfter(TString StopAfter) { fStopAfter = StopAfter; };
    void ReadSignalLogs(Bool_t ReadSignalLogs) { fReadSignalLogs = ReadSignalLogs; };
    void Initialize();

    /* Main ~ executed at runtime */
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    void EndOfEvent();
    virtual Bool_t UserNotify();

    /* Trees */
    void PrepareEventsBranches();
    void PrepareInjectedBranches();
    void PrepareMCParticlesBranches();
    void PrepareTracksBranches();
    void PrepareV0sBranches();
    void PrepareSexaquarksBranches();
    void PrepareKaonPairsBranches();

    void ClearEventsBranches();
    void ClearInjectedBranches();
    void ClearMCParticlesBranches();
    void ClearTracksBranches();
    void ClearV0sBranches();
    void ClearSexaquarksBranches();
    void ClearKaonPairsBranches();

    /* Histograms */
    void PrepareQAHistograms();
    void PrepareTracksHistograms();
    void PrepareV0Histograms();
    void PrepareSexaquarkHistograms();
    void PreparePosKaonPairHistograms();

    /* Cuts */
    void DefineTracksCuts(TString cuts_option);
    void DefineV0Cuts(TString cuts_option);
    void DefineSexaquarkCuts(TString cuts_option);
    void DefinePosKaonPairCuts(TString cuts_option);

    /* MC Generated */
    void ProcessMCGen();
    void ProcessSignalInteractions();

    /* Tracks */
    Bool_t PassesEventSelection();
    void ProcessTracks();
    Bool_t PassesTrackSelectionAs(Int_t HypothesisPID, Bool_t is_true, Bool_t is_secondary, Bool_t is_signal, AliESDtrack* track);
    void PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode);

    /* V0s, Anti-Sexaquarks, K+K+ pairs -- That Require True Info */
    void ProcessFindableV0s();
    void ProcessFindablePionPairs();
    void ProcessFindableSexaquarks();
    void ProcessFindablePosKaonPairs();

    /* V0s */
    void GetV0sFromESD(Bool_t onTheFly);
    void CustomV0Finder(Int_t pdgV0);
    void KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos);
    Bool_t PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, AliESDv0* v0, AliESDtrack* neg_track,
                          AliESDtrack* pos_track);
    Bool_t PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, KFParticleMother kfV0, KFParticle kfDaughterNeg,
                          KFParticle kfDaughterPos, TLorentzVector lvV0, TLorentzVector lvTrackNeg, TLorentzVector lvTrackPos);

    /* Sexaquark */
    void GeoSexaquarkFinder_ChannelA();
    void GeoSexaquarkFinder_ChannelD();
    void GeoSexaquarkFinder_ChannelE();
    void KalmanSexaquarkFinder_ChannelA();
    void KalmanSexaquarkFinder_ChannelD();
    void KalmanSexaquarkFinder_ChannelE();
    Bool_t PassesSexaquarkCuts_ChannelA(Bool_t isSignal, TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark, AliESDv0 esdAntiLambda,
                                        AliESDv0 esdKaonZeroShort, AliESDtrack* esdAntiLambdaNeg, AliESDtrack* esdAntiLambdaPos,
                                        AliESDtrack* esdKaonZeroShortNeg, AliESDtrack* esdKaonZeroShortPos);
    Bool_t PassesSexaquarkCuts_ChannelA(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark, KFParticle kfAntiLambda,
                                        KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos,
                                        KFParticle kfKaonZeroShortNeg, KFParticle kfKaonZeroShortPos);
    // Bool_t PassesSexaquarkCuts_ChannelD(TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark, AliESDv0 esdAntiLambda, AliESDtrack* esdPosKaon,
    // AliESDtrack* esdAntiLambdaNeg, AliESDtrack* esdAntiLambdaPos);
    Bool_t PassesSexaquarkCuts_ChannelD(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark, KFParticle kfAntiLambda,
                                        KFParticle kfPosKaon, KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos);
    // Bool_t PassesSexaquarkCuts_ChannelE(TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark, AliESDv0 esdAntiLambda, AliESDtrack* esdPosKaon,
    // AliESDtrack* esdAntiLambdaNeg, AliESDtrack* esdAntiLambdaPos);
    Bool_t PassesSexaquarkCuts_ChannelE(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark, KFParticle kfAntiLambda,
                                        KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos, KFParticle kfPosKaon, KFParticle kfPiMinus,
                                        KFParticle kfPiPlus);

    /* Channel H */
    void KalmanPosKaonPairFinder();
    Bool_t PassesPosKaonPairCuts(Bool_t isSignal, KFParticleMother kfPosKaonPair, KFParticle kfPosKaonA, KFParticle kfPosKaonB,
                                 TLorentzVector lvPosKaonPair, TLorentzVector lvPosKaonA, TLorentzVector lvPosKaonB);

    /* Utilities */
    void ClearContainers();
    Bool_t Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn, Double_t* lPreprocessxp,
                       const Double_t b);
    void GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b);
    void Evaluate(const Double_t* h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]);
    Float_t MCRec_GetImpactParameter(AliESDtrack* track);
    Double_t CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
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
    Double_t Calculate_TwoLinesDCA_v1(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0, TVector3& PCA1);
    Double_t Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0, TVector3& PCA1);

    /* Kalman Filter Utilities */
    KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge);
    KFVertex CreateKFVertex(const AliVVertex& vertex);
    KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis);
    void GetDStoParticleBz(float Bz, const KFParticleBase& p, const KFParticleBase& q, float dS[2], float dsdr[4][6], const float* param1 = 0,
                           const float* param2 = 0) const;

    /* External Files */
    Bool_t LoadLogsIntoTree();

   private:
    /* Settings ~ stored in Analysis Manager ~ all persistent */
    Bool_t fIsMC;            // kTRUE if MC simulation, kFALSE if data
    TString fSourceOfV0s;    // choose V0 finder: "kalman", "custom", "on-the-fly", "offline"
    Char_t fReactionID;      // reaction channel identifier, could be: 'A', 'D', 'E', 'H'
    Bool_t fDoQA;            // kTRUE to do QA
    TString fStopAfter;      // "MC", "Tracks", "Findables", "V0s" (default: "")
    Bool_t fReadSignalLogs;  // kTRUE to read and load signal logs

    /* AliRoot Objects */
    AliMCEvent* fMC;                //! MC event
    AliVVertex* fMC_PrimaryVertex;  //! MC gen. (or true) primary vertex
    AliESDEvent* fESD;              //! reconstructed event
    AliESDVertex* fPrimaryVertex;   //! primary vertex
    AliPIDResponse* fPIDResponse;   //! pid response object
    AliEventCuts fEventCuts;        //! event cuts
    Double_t fMagneticField;        //! magnetic field
    Int_t fRunNumber;               //! run number
    Int_t fDirNumber;               //! directory number
    Int_t fEventNumber;             //! event number
    Float_t fCentrality;            //! centrality percentile

    /* ROOT Objects */
    TDatabasePDG fPDG;                //!
    TList* fList_QA_Hists;            //!
    TList* fList_Tracks_Hists;        //!
    TList* fList_V0s_Hists;           //!
    TList* fList_Sexaquarks_Hists;    //!
    TList* fList_PosKaonPairs_Hists;  //!

    /* External Files */
    TString fAliEnPath;    //!
    Bool_t fIsFirstEvent;  //!

    /* Filled at Initialize() ~ persistent */
    std::vector<Int_t> fProductsPDG;                               // PDG codes if the immediate products of the reaction
    std::vector<Int_t> fFinalStateProductsPDG;                     // PDG codes of the final state products
    Int_t fStruckNucleonPDG;                                       // PDG code of the struck nucleon, could be: 2112 o 2212
    std::unordered_map<Int_t, Int_t> getNegPdgCode_fromV0PdgCode;  // PDG code of the negative daughter from the V0 PDG code
    std::unordered_map<Int_t, Int_t> getPosPdgCode_fromV0PdgCode;  // PDG code of the positive daughter from the V0 PDG code

    /*** Trees ***/
    TTree* fTree_Events;      //!
    TTree* fTree_Injected;    //!
    TTree* fTree_MC;          //!
    TTree* fTree_Tracks;      //!
    TTree* fTree_V0s;         //!
    TTree* fTree_Sexaquarks;  //!
    TTree* fTree_KaonPairs;   //!

    /** Branches **/

    Float_t tEvents_PV_TrueXv;     //!
    Float_t tEvents_PV_TrueYv;     //!
    Float_t tEvents_PV_TrueZv;     //!
    Float_t tEvents_PV_RecXv;      //!
    Float_t tEvents_PV_RecYv;      //!
    Float_t tEvents_PV_RecZv;      //!
    UInt_t tEvents_NTracks;        //!
    Bool_t tEvents_IsMB;           //!
    Bool_t tEvents_IsHighMultV0;   //!
    Bool_t tEvents_IsHighMultSPD;  //!
    Bool_t tEvents_IsCentral;      //!
    Bool_t tEvents_IsSemiCentral;  //!

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

    Int_t tMC_Idx;           //!
    Int_t tMC_PdgCode;       //!
    Int_t tMC_Idx_Mother;    //!
    Int_t tMC_NDaughters;    //!
    Int_t tMC_Idx_FirstDau;  //!
    Int_t tMC_Idx_LastDau;   //!
    Float_t tMC_Px;          //!
    Float_t tMC_Py;          //!
    Float_t tMC_Pz;          //!
    Float_t tMC_Xv_i;        //!
    Float_t tMC_Yv_i;        //!
    Float_t tMC_Zv_i;        //!
    UInt_t tMC_Status;       //!
    Bool_t tMC_IsSecondary;  //!
    Bool_t tMC_IsSignal;     //!
    Int_t tMC_ReactionID;    //!

    Int_t tTrack_Idx;             //!
    Int_t tTrack_Idx_True;        //!
    Float_t tTrack_Px;            //!
    Float_t tTrack_Py;            //!
    Float_t tTrack_Pz;            //!
    Short_t tTrack_Charge;        //!
    Float_t tTrack_NSigmaPion;    //!
    Float_t tTrack_NSigmaKaon;    //!
    Float_t tTrack_NSigmaProton;  //!
    Bool_t tTrack_IsTrue;         //!
    Bool_t tTrack_IsSecondary;    //!
    Bool_t tTrack_IsSignal;       //!
    Int_t tTrack_ReactionID;      //!

    Int_t tV0_Idx;           //!
    Int_t tV0_Idx_Pos;       //!
    Int_t tV0_Idx_Neg;       //!
    Int_t tV0_Idx_True;      //!
    Int_t tV0_PID;           //!
    Float_t tV0_Px;          //!
    Float_t tV0_Py;          //!
    Float_t tV0_Pz;          //!
    Float_t tV0_E;           //!
    Float_t tV0_Xv_f;        //!
    Float_t tV0_Yv_f;        //!
    Float_t tV0_Zv_f;        //!
    Float_t tV0_Neg_Px;      //!
    Float_t tV0_Neg_Py;      //!
    Float_t tV0_Neg_Pz;      //!
    Float_t tV0_Pos_Px;      //!
    Float_t tV0_Pos_Py;      //!
    Float_t tV0_Pos_Pz;      //!
    Bool_t tV0_IsTrue;       //!
    Bool_t tV0_IsSecondary;  //!
    Bool_t tV0_IsSignal;     //!
    Int_t tV0_ReactionID;    //!

    Int_t tSexaquark_Idx_AL;         //!
    Int_t tSexaquark_Idx_AL_Neg;     //!
    Int_t tSexaquark_Idx_AL_Pos;     //!
    Float_t tSexaquark_Px;           //!
    Float_t tSexaquark_Py;           //!
    Float_t tSexaquark_Pz;           //!
    Float_t tSexaquark_E;            //!
    Float_t tSexaquark_SV_Xv;        //!
    Float_t tSexaquark_SV_Yv;        //!
    Float_t tSexaquark_SV_Zv;        //!
    Float_t tSexaquark_DecayLength;  //!
    Float_t tSexaquark_CPAwrtPV;     //!
    Float_t tSexaquark_DCAwrtPV;     //!
    Float_t tSexaquark_Chi2ndf;      //!
    Bool_t tSexaquark_IsHybrid;      //!
    Bool_t tSexaquark_IsSignal;      //!
    Int_t tSexaquark_ReactionID;     //!

    Int_t tSexaquarkA_Idx_K0S;        //!
    Int_t tSexaquarkA_Idx_K0S_Neg;    //!
    Int_t tSexaquarkA_Idx_K0S_Pos;    //!
    Float_t tSexaquarkA_DCAbtwV0s;    //!
    Float_t tSexaquarkA_DCAV0aSV;     //!
    Float_t tSexaquarkA_DCAV0bSV;     //!
    Float_t tSexaquarkA_DCAV0anegSV;  //!
    Float_t tSexaquarkA_DCAV0aposSV;  //!
    Float_t tSexaquarkA_DCAV0bnegSV;  //!
    Float_t tSexaquarkA_DCAV0bposSV;  //!

    Int_t tSexaquarkD_Idx_PosKaon;   //!
    Float_t tSexaquarkD_DCAV0Ba;     //!
    Float_t tSexaquarkD_DCAV0negBa;  //!
    Float_t tSexaquarkD_DCAV0posBa;  //!
    Float_t tSexaquarkD_DCAV0SV;     //!
    Float_t tSexaquarkD_DCABaSV;     //!
    Float_t tSexaquarkD_DCAV0negSV;  //!
    Float_t tSexaquarkD_DCAV0posSV;  //!

    Int_t tSexaquarkE_Idx_PP;        //!
    Int_t tSexaquarkE_Idx_PP_Neg;    //!
    Int_t tSexaquarkE_Idx_PP_Pos;    //!
    Int_t tSexaquarkE_Idx_PosKaon;   //!
    Float_t tSexaquarkE_DCAV0SV;     //!
    Float_t tSexaquarkE_DCAV0negSV;  //!
    Float_t tSexaquarkE_DCAV0posSV;  //!
    Float_t tSexaquarkE_DCApkSV;     //!
    Float_t tSexaquarkE_DCApmSV;     //!
    Float_t tSexaquarkE_DCAppSV;     //!
    Float_t tSexaquarkE_DCApkV0;     //!
    Float_t tSexaquarkE_DCApmV0;     //!
    Float_t tSexaquarkE_DCAppV0;     //!
    Float_t tSexaquarkE_DCApmPK;     //!
    Float_t tSexaquarkE_DCAppPK;     //!

    Int_t tKaonPair_Idx;             //!
    Int_t tKaonPair_Idx_KaonA;       //!
    Int_t tKaonPair_Idx_KaonB;       //!
    Float_t tKaonPair_Px;            //!
    Float_t tKaonPair_Py;            //!
    Float_t tKaonPair_Pz;            //!
    Float_t tKaonPair_E;             //!
    Float_t tKaonPair_SV_Xv;         //!
    Float_t tKaonPair_SV_Yv;         //!
    Float_t tKaonPair_SV_Zv;         //!
    Float_t tKaonPair_DecayLength;   //!
    Float_t tKaonPair_DCAwrtPV;      //!
    Float_t tKaonPair_DCAbtwKaons;   //!
    Float_t tKaonPair_DCAkaSV;       //!
    Float_t tKaonPair_DCAkbSV;       //!
    Float_t tKaonPair_OpeningAngle;  //!
    Float_t tKaonPair_Chi2ndf;       //!
    Bool_t tKaonPair_IsHybrid;       //!
    Bool_t tKaonPair_IsSignal;       //!
    Int_t tKaonPair_ReactionID;      //!

    /*** Histograms ***/

    /** QA Histograms **/

    std::map<TString, TH1F*> fHist_QA;    //!
    std::map<TString, TH2F*> f2DHist_QA;  //!

    /** Tracks Histograms **/

    /*
     key: `pdg code`
     - `0`    : MC tracks
     - `10`   : SECONDARY MC tracks
     - `20`   : SIGNAL MC tracks
     - `30+`  : effect of cuts on found
     - `50`   : found tracks that passed all cuts
     - `60+`  : effect of cuts on TRUE tracks
     - `80`   : found TRUE tracks that passed all cuts
     - `90+`  : effect of cuts on SECONDARY tracks
     - `110`  : found SECONDARY tracks that passed all cuts
     - `120+` : effect of cuts on SIGNAL tracks
     - `140`  : found SIGNAL tracks that passed all cuts
    */
    std::unordered_map<Int_t, TH1F*> fHist_Tracks_Bookkeep;                      //!
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_Tracks;  //! key: `stage, set, pdg, property`
    std::map<TString, TH2F*> f2DHist_Tracks_TPCsignal;                           //! key: `set`, value: histogram

    /** V0s Histograms **/

    /*
     key: `V0 pdg code`
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
    std::unordered_map<Int_t, TH1F*> fHist_V0s_Bookkeep;                      //!
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_V0s;  //! key: `stage, set, pdg, property`
    std::map<TString, TH2F*> f2DHist_V0s_ArmenterosPodolanski;                //! key: `set`

    /** Anti-Sexaquark Histograms **/

    /*
     - `0`   : MC Gen.
     - `10`  : findable
     - `20+` : effect of cuts on found candidates
     - `50`  : found candidates that passed all cuts
     - `60+` : effect of cuts on SIGNAL candidates
     - `90`  : found SIGNAL candidates that passed all cuts
    */
    TH1D* fHist_AntiSexaquarks_Bookkeep;                                          //!
    std::map<std::tuple<TString, TString, TString>, TH1D*> fHist_AntiSexaquarks;  //! key: `stage, set, property`

    /** Channel H / Pos. Kaon Pairs Histograms **/

    /*
     - `0`   : MC Gen.
     - `10`  : findable
     - `20+` : effect of cuts on found candidates
     - `50`  : found candidates that passed all cuts
     - `60+` : effect of cuts on SIGNAL candidates
     - `90`  : found SIGNAL candidates that passed all cuts
    */
    TH1D* fHist_PosKaonPairs_Bookkeep;                                          //!
    std::map<std::tuple<TString, TString, TString>, TH1D*> fHist_PosKaonPairs;  //! key: `stage, set, property`

    /*** Containers -- vectors and hash tables ***/
    // NOTE: when adding a new one, don't forget to clear it on ClearContainers()!

    /* filled at `ProcessMCGen()` */
    std::unordered_map<Int_t, Int_t> getPdgCode_fromMcIdx;  //!

    std::unordered_map<Int_t, Bool_t> isMcIdxSignal;     //!
    std::unordered_map<Int_t, Bool_t> isMcIdxSecondary;  //!

    std::unordered_map<Int_t, UInt_t> getReactionIdx_fromMcIdx;               //!
    std::unordered_map<UInt_t, std::vector<Int_t>> getMcIdx_fromReactionIdx;  //!
    std::unordered_map<UInt_t, Float_t> getPt_fromReactionIdx;                //!
    std::unordered_map<UInt_t, Float_t> getRadius_fromReactionIdx;            //!

    std::unordered_map<Int_t, Bool_t> doesMcIdxHaveMother;      //!
    std::unordered_map<Int_t, Int_t> getMotherMcIdx_fromMcIdx;  //!
    std::unordered_map<Int_t, Int_t> getNegDauMcIdx_fromMcIdx;  //!
    std::unordered_map<Int_t, Int_t> getPosDauMcIdx_fromMcIdx;  //!

    std::vector<Int_t> mcIndicesOfTrueV0s;  //!

    /* filled at `ProcessTracks()` */
    std::unordered_map<Int_t, Int_t> getMcIdx_fromEsdIdx;                       //!
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsProton;      //! necessary for `GetV0sFromESD()`
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsAntiProton;  //!
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsPosKaon;     //!
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsNegKaon;     //!
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsPiPlus;      //!
    std::unordered_map<Int_t, Bool_t> esdIdxPassesTrackSelection_AsPiMinus;     //!
    std::unordered_map<Int_t, Int_t> getEsdIdx_fromMcIdx;                       //!

    std::unordered_map<Int_t, Bool_t> isEsdIdxSignal;                          //!
    std::unordered_map<Int_t, UInt_t> getReactionIdx_fromEsdIdx;               //!
    std::unordered_map<UInt_t, std::vector<Int_t>> getEsdIdx_fromReactionIdx;  //!

    std::unordered_map<Int_t, Bool_t> doesEsdIdxHaveMother;      //!
    std::unordered_map<Int_t, Int_t> getMotherMcIdx_fromEsdIdx;  //!
    std::unordered_map<Int_t, Int_t> getNegDauEsdIdx_fromMcIdx;  //!
    std::unordered_map<Int_t, Int_t> getPosDauEsdIdx_fromMcIdx;  //!

    std::vector<Int_t> esdIndicesOfAntiProtonTracks;  //!
    std::vector<Int_t> esdIndicesOfPosKaonTracks;     //!
    std::vector<Int_t> esdIndicesOfPiMinusTracks;     //!
    std::vector<Int_t> esdIndicesOfPiPlusTracks;      //!

    /* filled at `GetV0sFromESD() and [Custom|Kalman]V0Finder()` */
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromAntiLambdaIdx;     //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromAntiLambdaIdx;     //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromKaonZeroShortIdx;  //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromKaonZeroShortIdx;  //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromPionPairIdx;       //!
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromPionPairIdx;       //!

    std::vector<AliESDv0> esdAntiLambdas;     //!
    std::vector<AliESDv0> esdKaonsZeroShort;  //!

    std::vector<KFParticleMother> kfAntiLambdas;     //!
    std::vector<KFParticleMother> kfKaonsZeroShort;  //!
    std::vector<KFParticleMother> kfPionPairs;       //!

    /*** Cuts ~ persistent, because they are set on Initialize() ***/

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
    std::unordered_map<Int_t, Float_t> kMin_Track_Pt;  //

    /** V0s **/
    std::unordered_map<Int_t, Float_t> kMin_V0_Mass;            //
    std::unordered_map<Int_t, Float_t> kMax_V0_Mass;            //
    std::unordered_map<Int_t, Float_t> kMax_V0_Eta;             //
    std::unordered_map<Int_t, Float_t> kMax_V0_ArmPtOverAlpha;  //
    std::unordered_map<Int_t, Float_t> kMin_V0_Pt;              //
    std::unordered_map<Int_t, Float_t> kMin_V0_Radius;          //
    std::unordered_map<Int_t, Float_t> kMin_V0_DecayLength;     //
    std::unordered_map<Int_t, Float_t> kMax_V0_DecayLength;     //
    std::unordered_map<Int_t, Float_t> kMin_V0_CPAwrtPV;        //
    std::unordered_map<Int_t, Float_t> kMax_V0_CPAwrtPV;        //
    std::unordered_map<Int_t, Float_t> kMin_V0_DCAwrtPV;        //
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAbtwDau;       //
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAnegV0;        //
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAposV0;        //
    std::unordered_map<Int_t, Float_t> kMax_V0_ImprvDCAbtwDau;  //
    std::unordered_map<Int_t, Float_t> kMax_V0_Chi2;            //
    std::unordered_map<Int_t, Float_t> kMax_V0_Chi2ndf;         //

    /** Anti-Sexaquark Candidates **/
    Float_t kMin_Sexa_Mass;         //
    Float_t kMax_Sexa_Mass;         //
    Float_t kMin_Sexa_Pt;           //
    Float_t kMax_Sexa_Eta;          //
    Float_t kMax_Sexa_Rapidity;     //
    Float_t kMin_Sexa_DecayLength;  //
    Float_t kMin_Sexa_Radius;       //
    Float_t kMin_Sexa_CPAwrtPV;     //
    Float_t kMax_Sexa_CPAwrtPV;     //
    Float_t kMin_Sexa_DCAwrtPV;     //
    Float_t kMax_Sexa_DCAwrtPV;     //
    Float_t kMax_Sexa_Chi2ndf;      //
    /* Channel A */
    Float_t kMax_Sexa_DCAbtwV0s;    //
    Float_t kMax_Sexa_DCAv0aSV;     //
    Float_t kMax_Sexa_DCAv0bSV;     //
    Float_t kMax_Sexa_DCAv0anegSV;  //
    Float_t kMax_Sexa_DCAv0aposSV;  //
    Float_t kMax_Sexa_DCAv0bnegSV;  //
    Float_t kMax_Sexa_DCAv0bposSV;  //
    /* Channel D */
    Float_t kMax_Sexa_DCAbaSV;     //
    Float_t kMax_Sexa_DCAv0ba;     //
    Float_t kMax_Sexa_DCAv0negba;  //
    Float_t kMax_Sexa_DCAv0posba;  //
    /* Channel E */
    Float_t kMax_Sexa_DCAv0SV;     //
    Float_t kMax_Sexa_DCAv0negSV;  //
    Float_t kMax_Sexa_DCAv0posSV;  //
    Float_t kMax_Sexa_DCApkSV;     //
    Float_t kMax_Sexa_DCApmSV;     //
    Float_t kMax_Sexa_DCAppSV;     //
    Float_t kMax_Sexa_DCAv0pk;     //
    Float_t kMax_Sexa_DCAv0pm;     //
    Float_t kMax_Sexa_DCAv0pp;     //
    Float_t kMax_Sexa_DCApkpm;     //
    Float_t kMax_Sexa_DCApkpp;     //

    /** Secondary K+K+ Pairs **/
    Float_t kMin_PosKaonPair_Radius;       //
    Float_t kMin_PosKaonPair_Pt;           //
    Float_t kMin_PosKaonPair_DCAwrtPV;     //
    Float_t kMax_PosKaonPair_DCAbtwKaons;  //
    Float_t kMax_PosKaonPair_DCApkaSV;     //
    Float_t kMax_PosKaonPair_DCApkbSV;     //
    Float_t kMax_PosKaonPair_Chi2ndf;      //

    AliAnalysisTaskSexaquark(const AliAnalysisTaskSexaquark&);             // not implemented
    AliAnalysisTaskSexaquark& operator=(const AliAnalysisTaskSexaquark&);  // not implemented

    /// \cond CLASSDEF
    ClassDef(AliAnalysisTaskSexaquark, 81);
    /// \endcond
};

#endif
