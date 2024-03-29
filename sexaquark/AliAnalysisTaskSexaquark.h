#ifndef AliAnalysisTaskSexaquark_H
#define AliAnalysisTaskSexaquark_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "TArray.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TH1.h"
#include "TH1F.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliPIDResponse.h"

#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliVVertex.h"

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

struct V0_tt {
    //
    // Default common structure for the V0s
    // (since I want to study 3 alternatives, let's keep the same format for each)
    //
    Int_t Idx_Pos;          // index of positive track
    Int_t Idx_Neg;          // index of negative track
    Float_t Px;             // x-component of V0 momentum
    Float_t Py;             // y-component of V0 momentum
    Float_t Pz;             // z-component of V0 momentum
    Float_t X;              // x-coordinate of V0
    Float_t Y;              // y-coordinate of V0
    Float_t Z;              // z-coordinate of V0
    Float_t Pos_Px;         // x-component of positive track momentum at V0 position
    Float_t Pos_Py;         // y-component of positive track momentum at V0 position
    Float_t Pos_Pz;         // z-component of positive track momentum at V0 position
    Float_t Neg_Px;         // x-component of negative track momentum at V0 position
    Float_t Neg_Py;         // y-component of negative track momentum at V0 position
    Float_t Neg_Pz;         // z-component of negative track momentum at V0 position
    Bool_t isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    Float_t E_asK0;         // energy of V0, assuming it's a K0
    Float_t E_asAL;         // energy of V0, assuming it's an anti-lambda
    Bool_t couldBeK0;       // kTRUE if K0 candidate, kFALSE if not
    Bool_t couldBeAL;       // kTRUE if anti-lambda candidate, kFALSE if not
    Float_t Chi2;           // chi2 value from the Kalman Filter (?)
    Float_t DCA_Daughters;  // distance of closest approach between daughters
    Float_t IP_wrtPV;       // impact parameter w.r.t. Primary Vertex
    Float_t CPA_wrtPV;      // cosine of pointing angle w.r.t. Primary Vertex
    Float_t ArmAlpha;       // Armenteros-Podolanski variable alpha
    Float_t ArmPt;          // Armenteros-Podolanski variable Pt
    Float_t DecayLength;    // distance between PV and V0
    Float_t DCA_wrtPV;      // distance of closest approach to Primary Vertex
};

struct Event_tt {
    //
    // This is the structure of the output TTree
    //

    /* MC particles */

    Int_t N_MCGen;                          // number of MC particles
    std::vector<Float_t> MC_Px;             // x-component of true momentum
    std::vector<Float_t> MC_Py;             // y-component of true momentum
    std::vector<Float_t> MC_Pz;             // z-component of true momentum
    std::vector<Float_t> MC_X;              // initial x-coordinate of generation vertex
    std::vector<Float_t> MC_Y;              // initial y-coordinate of generation vertex
    std::vector<Float_t> MC_Z;              // initial z-coordinate of generation vertex
    std::vector<Float_t> MC_Xf;             // final x-coordinate of generation vertex
    std::vector<Float_t> MC_Yf;             // final y-coordinate of generation vertex
    std::vector<Float_t> MC_Zf;             // final z-coordinate of generation vertex
    std::vector<Int_t> MC_PID;              // PDG code
    std::vector<Int_t> MC_Mother;           // index of mother
    std::vector<Int_t> MC_PID_Mother;       // PID of mother
    std::vector<Int_t> MC_PID_GrandMother;  // PID of grandmother
    std::vector<Int_t> MC_NDaughters;       // number of daughters
    std::vector<Int_t> MC_FirstDau;         // index of first daughter
    std::vector<Int_t> MC_LastDau;          // index of last daughter
    std::vector<Int_t> MC_Status;           // MC status code
    std::vector<Bool_t> MC_isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background

    /* MC Rec. Tracks */

    Int_t N_MCRec;                          // number of MC reconstructed tracks
    std::vector<Int_t> Idx_True;            // index of true MC particle
    std::vector<Float_t> Rec_Px;            // x-component of reconstructed momentum
    std::vector<Float_t> Rec_Py;            // y-component of reconstructed momentum
    std::vector<Float_t> Rec_Pz;            // z-component of reconstructed momentum
    std::vector<Short_t> Rec_Charge;        // measured charge
    std::vector<Float_t> Rec_IP_wrtPV;      // transverse impact parameter w.r.t. Primary Vertex
    std::vector<Float_t> Rec_NSigmaPion;    // absolute value of likeness to be a charged pion (closer to 0, the most likely)
    std::vector<Float_t> Rec_NSigmaProton;  // absolute value of likeness to be a proton (closer to 0, the most likely)
    std::vector<Short_t> Rec_NClustersTPC;  // number of clusters assigned in the TPC
    std::vector<Bool_t> Rec_isDuplicate;    // kTRUE if track is a duplicate, kFALSE if not (v1: compared to true particle)
    std::vector<Bool_t> Rec_isSimilar;      // kTRUE if track is a duplicate, kFALSE if not (v2: compared to similar particles)
    std::vector<Bool_t> Rec_isSignal;       // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    std::vector<Float_t> Rec_HelixParam0;   // helix param. 0
    std::vector<Float_t> Rec_HelixParam1;   // helix param. 1
    std::vector<Float_t> Rec_HelixParam2;   // helix param. 2
    std::vector<Float_t> Rec_HelixParam3;   // helix param. 3
    std::vector<Float_t> Rec_HelixParam4;   // helix param. 4
    std::vector<Float_t> Rec_HelixParam5;   // helix param. 5

    /* V0s */

    Int_t N_V0s;                            // number of formed V0s
    std::vector<Int_t> Idx_Pos;             // index of positive daughter
    std::vector<Int_t> Idx_Neg;             // index of negative daughter
    std::vector<Float_t> V0_Px;             // x-component of V0 momentum
    std::vector<Float_t> V0_Py;             // y-component of V0 momentum
    std::vector<Float_t> V0_Pz;             // z-component of V0 momentum
    std::vector<Float_t> V0_X;              // x-coordinate of V0
    std::vector<Float_t> V0_Y;              // y-coordinate of V0
    std::vector<Float_t> V0_Z;              // z-coordinate of V0
    std::vector<Float_t> Pos_Px;            // x-component of positive track momentum at V0 position
    std::vector<Float_t> Pos_Py;            // y-component of positive track momentum at V0 position
    std::vector<Float_t> Pos_Pz;            // z-component of positive track momentum at V0 position
    std::vector<Float_t> Neg_Px;            // x-component of negative track momentum at V0 position
    std::vector<Float_t> Neg_Py;            // y-component of negative track momentum at V0 position
    std::vector<Float_t> Neg_Pz;            // z-component of negative track momentum at V0 position
    std::vector<Bool_t> V0_isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    std::vector<Float_t> V0_E_asK0;         // energy of V0, assuming it's a K0
    std::vector<Float_t> V0_E_asAL;         // energy of V0, assuming it's an anti-lambda
    std::vector<Bool_t> V0_couldBeK0;       // kTRUE if K0 candidate, kFALSE if not
    std::vector<Bool_t> V0_couldBeAL;       // kTRUE if anti-lambda candidate, kFALSE if not
    std::vector<Float_t> V0_Chi2;           // chi2 value from the Kalman Filter (?)
    std::vector<Float_t> V0_DCA_Daughters;  // distance of closest approach between daughters
    std::vector<Float_t> V0_IP_wrtPV;       // impact parameter w.r.t. Primary Vertex
    std::vector<Float_t> V0_CPA_wrtPV;      // cosine of pointing angle w.r.t. Primary Vertex
    std::vector<Float_t> V0_ArmAlpha;       // Armenteros-Podolanski variable alpha
    std::vector<Float_t> V0_ArmPt;          // Armenteros-Podolanski variable Pt
    std::vector<Float_t> V0_DecayLength;    // distance between PV and V0
    std::vector<Float_t> V0_DCA_wrtPV;      // distance of closest approach to Primary Vertex
};

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
    AliAnalysisTaskSexaquark(const char* name, Bool_t IsMC, TString SourceOfV0s, TString SimulationSet);
    virtual ~AliAnalysisTaskSexaquark();
    virtual void Terminate(Option_t* option) { return; }

    /* Initialization */
    virtual void UserCreateOutputObjects();
    void CheckForInputErrors();
    void Initialize();
    void PrepareTracksHistograms();
    void PrepareV0Histograms();
    void PrepareAntiSexaquarkHistograms();

    /* Main */
    virtual void UserExec(Option_t* option);

    /* Cuts */
    void DefineTracksCuts(TString cuts_option);
    void DefineV0Cuts(TString cuts_option);
    void DefineSexaquarkCuts(TString cuts_option);

    /* MC Generated */
    void ProcessMCGen();
    void GetDaughtersInfo(AliMCParticle* mcPart, Int_t& mcIdxNegDaughter, Int_t& mcIdxPosDaughter, TVector3& decay_vertex);
    void ProcessSignalInteractions();

    /* Tracks */
    void ProcessTracks();
    Bool_t PassesTrackSelection(AliESDtrack* track);
    void PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode);

    /* V0s and Anti-Sexaquarks -- That Require True Info */
    void ProcessFindableV0s();
    void ProcessFindableSexaquarks();

    /* V0s */
    void GetV0sFromESD(Bool_t onTheFly);
    void CustomV0Finder(Int_t pdgV0);
    void KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos);
    Bool_t PassesV0Cuts(Int_t pdgV0, AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track);
    Bool_t PassesV0Cuts(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfDaughterNeg, KFParticle kfDaughterPos, TLorentzVector lvV0,
                        TLorentzVector lvTrackNeg, TLorentzVector lvTrackPos);

    /* Sexaquark */
    void GeoSexaquarkFinder_ChannelA();
    // void GeoSexaquarkFinder_ChannelD();
    // void GeoSexaquarkFinder_ChannelE();
    void KalmanSexaquarkFinder_ChannelA();
    // void KalmanSexaquarkFinder_ChannelD();
    // void KalmanSexaquarkFinder_ChannelE();
    Bool_t PassesSexaquarkCuts_ChannelA(TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark, AliESDv0 esdAntiLambda, AliESDv0 esdKaonZeroShort,
                                        AliESDtrack* esdAntiLambdaNeg, AliESDtrack* esdAntiLambdaPos, AliESDtrack* esdKaonZeroShortNeg,
                                        AliESDtrack* esdKaonZeroShortPos);
    Bool_t PassesSexaquarkCuts_ChannelA(KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark, KFParticle kfAntiLambda,
                                        KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos,
                                        KFParticle kfKaonZeroShortNeg, KFParticle kfKaonZeroShortPos);
    // Bool_t PassesSexaquarkCuts_ChannelD();
    // Bool_t PassesSexaquarkCuts_ChannelD();
    // Bool_t PassesSexaquarkCuts_ChannelE();
    // Bool_t PassesSexaquarkCuts_ChannelE();

    /* Utilities */
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

    /* Tree Operations */
    Event_tt fEvent;
    void SetBranches();
    void MCGen_PushBack(Int_t evt_mc);
    void MCRec_PushBack(Int_t evt_track);
    void V0_PushBack(Int_t evt_v0);
    void FillTree();
    void ClearEvent();

    /* Kalman Filter Utilities */
    KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge);
    KFVertex CreateKFVertex(const AliVVertex& vertex);
    KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis);

   private:
    /* Input options */
    Bool_t fIsMC;                    //
    TString fSourceOfV0s;            // choose V0 finder: "offline", "on-the-fly", "custom", "kalman"
    TString fSimulationSet;          // <ReactionID><InjectedSexaquarkMass>
    Char_t fReactionID;              // could be: 'A', 'D', 'E', 'H'
    Float_t fInjectedSexaquarkMass;  // (in GeV/c^2), could be: 1.73, 1.8, 1.87, 1.94, 2.01
    // reaction channel, could be:
    // - 'A' = "AntiSexaquark,N->AntiLambda,K0S"
    // - 'D' = "AntiSexaquark,P->AntiLambda,K+"
    // - 'E' = "AntiSexaquark,P->AntiLambda,K+,pi-,pi+"
    // - 'H' = "AntiSexaquark,P->AntiProton,K+,K+,pi0"
    TString fReactionChannel;
    std::vector<Int_t> fProductsPDG;            //
    std::vector<Int_t> fFinalStateProductsPDG;  //
    Int_t fStruckNucleonPDG;                    // PDG Code of the struck nucleon, could ber: 2112 o 2212

    /* AliRoot objects */
    AliMCEvent* fMC;                // MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent* fESD;              // reconstructed event
    AliESDVertex* fPrimaryVertex;   // primary vertex
    AliPIDResponse* fPIDResponse;   // pid response object
    Double_t fMagneticField;        // magnetic field
    std::unordered_map<Int_t, Int_t> getNegPdgCode_fromV0PdgCode;
    std::unordered_map<Int_t, Int_t> getPosPdgCode_fromV0PdgCode;

    /* ROOT objects */
    TDatabasePDG fPDG;
    TList* fOutputListOfTrees;
    TTree* fTree;
    TList* fOutputListOfHists;

    /** Tracks Histograms **/

    /*
     key: `pdg code`
     - `0` : MC tracks
     - `1` : secondary MC tracks
     - `2` : signal MC tracks
     (PENDING)
     -- PENDING: in a sim.set of reactionID='A', signal K+ appear where they shouldn't exist, because of signal particles that were mis-id
     -- Note: should I only, besides if they're signal or not, use their true id? I think so, because it's a bookkeeping of TRUE particles
    */
    std::unordered_map<Int_t, TH1F*> fHist_Tracks_Bookkeep;
    // key: `stage, set, pdg, property`, value: histogram
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_Tracks;
    // key: `set`, value: histogram
    std::map<TString, TH2F*> f2DHist_Tracks_TPCsignal;

    /** V0s Histograms **/

    /*
     key: `V0 pdg code`
     - `0` : MC V0s
     - `1` : MC V0s w/ charged decay products
     - `2` : secondary MC V0s
     - `3` : secondary MC V0s w/ charged decay products
     - `4` : signal MC V0s
     - `5` : signal MC V0s w/ charged decay products
     - `6` : findable true V0s
     - `7` : findable true secondary V0s
     - `8` : findable true signal V0s
     - `30` : found V0s
     - `31` : found true V0s
     - `32` : found true secondary V0s
     - `33` : found true signal V0s
    */
    std::unordered_map<Int_t, TH1F*> fHist_V0s_Bookkeep;
    // key: `stage, set, pdg, property`, value: histogram
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_V0s;
    // key: `set`, value: histogram
    std::map<TString, TH2F*> f2DHist_V0s_ArmenterosPodolanski;

    /** Anti-Sexaquark Histograms **/
    /*
     - `0`    : MC Gen.
     - `1`    : Findable
     - `2-18` : Effect of cuts
     - `19`   : Found
     - `20`   : Found signal
    */
    TH1F* fHist_AntiSexaquarks_Bookkeep;
    // key: `stage, set, property`, value: histogram
    std::map<std::tuple<TString, TString, TString>, TH1F*> fHist_AntiSexaquarks;

    /** Containers -- vectors and hash tables **/

    /* filled at `ProcessMCGen()` */
    std::unordered_map<Int_t, Int_t> getPdgCode_fromMcIdx;  //

    std::unordered_map<Int_t, Bool_t> isMcIdxSignal;     //
    std::unordered_map<Int_t, Bool_t> isMcIdxSecondary;  //

    std::unordered_map<Int_t, UInt_t> getReactionIdx_fromMcIdx;               //
    std::unordered_map<UInt_t, std::vector<Int_t>> getMcIdx_fromReactionIdx;  //

    std::unordered_map<Int_t, Bool_t> doesMcIdxHaveMother;      //
    std::unordered_map<Int_t, Int_t> getMotherMcIdx_fromMcIdx;  //
    std::unordered_map<Int_t, Int_t> getNegDauMcIdx_fromMcIdx;  //
    std::unordered_map<Int_t, Int_t> getPosDauMcIdx_fromMcIdx;  //

    std::vector<Int_t> mcIndicesOfTrueV0s;  //

    /* filled at `ProcessTracks()` */
    std::unordered_map<Int_t, Int_t> getMcIdx_fromEsdIdx;            //
    std::unordered_map<Int_t, Bool_t> doesEsdIdxPassTrackSelection;  //
    std::unordered_map<Int_t, Int_t> getEsdIdx_fromMcIdx;            //

    std::unordered_map<Int_t, Bool_t> isEsdIdxSignal;                          //
    std::unordered_map<Int_t, UInt_t> getReactionIdx_fromEsdIdx;               //
    std::unordered_map<UInt_t, std::vector<Int_t>> getEsdIdx_fromReactionIdx;  //

    std::unordered_map<Int_t, Bool_t> doesEsdIdxHaveMother;      //
    std::unordered_map<Int_t, Int_t> getMotherMcIdx_fromEsdIdx;  //
    std::unordered_map<Int_t, Int_t> getNegDauEsdIdx_fromMcIdx;  //
    std::unordered_map<Int_t, Int_t> getPosDauEsdIdx_fromMcIdx;  //

    std::vector<Int_t> esdIndicesOfAntiProtonTracks;  //
    std::vector<Int_t> esdIndicesOfPosKaonTracks;     //
    std::vector<Int_t> esdIndicesOfPiMinusTracks;     //
    std::vector<Int_t> esdIndicesOfPiPlusTracks;      //

    /* filled at `GetV0sFromESD() and [Custom|Kalman]V0Finder()` */
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromAntiLambdaIdx;     //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromAntiLambdaIdx;     //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromKaonZeroShortIdx;  //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromKaonZeroShortIdx;  //

    std::vector<AliESDv0> esdAntiLambdas;     //
    std::vector<AliESDv0> esdKaonsZeroShort;  //

    std::vector<KFParticleMother> kfAntiLambdas;     //
    std::vector<KFParticleMother> kfKaonsZeroShort;  //

    void ClearContainers();

    /*** Cuts ***/

    /* Track selection */
    Float_t kMax_NSigma_Pion;
    Float_t kMax_NSigma_Kaon;
    Float_t kMax_NSigma_Proton;
    Float_t kMax_Track_Eta;
    Float_t kMin_Track_NTPCClusters;
    Float_t kMax_Track_Chi2PerNTPCClusters;
    Bool_t kTurnedOn_Track_StatusCuts;
    Bool_t kTurnedOn_Track_RejectKinks;
    Float_t kMin_Track_DCAwrtPV;
    std::unordered_map<Int_t, Float_t> kMin_Track_Pt;

    /* V0s */
    std::unordered_map<Int_t, Float_t> kMin_V0_Mass;
    std::unordered_map<Int_t, Float_t> kMax_V0_Mass;
    std::unordered_map<Int_t, Float_t> kMax_V0_Eta;
    std::unordered_map<Int_t, Float_t> kMax_V0_ArmPtOverAlpha;
    std::unordered_map<Int_t, Float_t> kMin_V0_Pt;
    std::unordered_map<Int_t, Float_t> kMin_V0_Radius;
    std::unordered_map<Int_t, Float_t> kMin_V0_DecayLength;
    std::unordered_map<Int_t, Float_t> kMax_V0_DecayLength;
    std::unordered_map<Int_t, Float_t> kMin_V0_CPAwrtPV;
    std::unordered_map<Int_t, Float_t> kMax_V0_CPAwrtPV;
    std::unordered_map<Int_t, Float_t> kMin_V0_DCAwrtPV;
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAbtwDau;
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAnegV0;
    std::unordered_map<Int_t, Float_t> kMax_V0_DCAposV0;
    std::unordered_map<Int_t, Float_t> kMax_V0_ImprvDCAbtwDau;
    std::unordered_map<Int_t, Float_t> kMax_V0_Chi2;
    std::unordered_map<Int_t, Float_t> kMax_V0_Chi2ndf;

    /* Anti-Sexaquark candidates */
    Float_t kMin_Sexa_Mass, kMax_Sexa_Mass;
    Float_t kMin_Sexa_Pt;
    Float_t kMax_Sexa_Eta;
    Float_t kMax_Sexa_Rapidity;
    Float_t kMin_Sexa_DecayLength;
    Float_t kMin_Sexa_Radius;
    Float_t kMin_Sexa_CPAwrtPV, kMax_Sexa_CPAwrtPV;
    Float_t kMin_Sexa_DCAwrtPV, kMax_Sexa_DCAwrtPV;
    Float_t kMax_Sexa_DCAbtwV0s;
    Float_t kMax_Sexa_DCAv0aSV;
    Float_t kMax_Sexa_DCAv0bSV;
    Float_t kMax_Sexa_DCAv0anegSV;
    Float_t kMax_Sexa_DCAv0aposSV;
    Float_t kMax_Sexa_DCAv0bnegSV;
    Float_t kMax_Sexa_DCAv0bposSV;
    Float_t kMax_Sexa_Chi2ndf;

    AliAnalysisTaskSexaquark(const AliAnalysisTaskSexaquark&);             // not implemented
    AliAnalysisTaskSexaquark& operator=(const AliAnalysisTaskSexaquark&);  // not implemented

    ClassDef(AliAnalysisTaskSexaquark, 1);
};

#endif
