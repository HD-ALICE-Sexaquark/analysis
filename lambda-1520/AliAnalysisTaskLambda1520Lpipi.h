#ifndef AliAnalysisTaskLambda1520Lpipi_h
#define AliAnalysisTaskLambda1520Lpipi_h

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "TArray.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TList.h"
#include "TLorentzVector.h"
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

#define HomogeneousField  // homogenous field in z direction, required by KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFVertex.h"

class AliPIDResponse;
class KFParticle;
class KFVertex;

/*
 Auxiliary class to make use of protected function KFParticleBase::GetMeasurement()
 [copied from /PWGLF/.../AliAnalysisTaskDoubleHypNucTree.h]
*/
class KFParticleMother : public KFParticle {
   public:
    Bool_t CheckDaughter(KFParticle daughter) {
        Float_t m[8], mV[36], D[3][3];
        if (KFParticleBase::GetMeasurement(daughter, m, mV, D)) return kTRUE;
        return kFALSE;
    }
};

class AliAnalysisTaskLambda1520Lpipi : public AliAnalysisTaskSE {

   public:
    AliAnalysisTaskLambda1520Lpipi();
    AliAnalysisTaskLambda1520Lpipi(const char* name, Bool_t IsMC);
    virtual ~AliAnalysisTaskLambda1520Lpipi();
    virtual void Terminate(Option_t* option) { return; }

    /* Initialization */
    virtual void UserCreateOutputObjects();
    void CheckForInputErrors();
    void Initialize();
    void PrepareTracksHistograms();
    void PrepareV0Histograms();
    void PrepareLambda1520Histograms();

    /* Main */
    virtual void UserExec(Option_t* option);

    /* Cuts */
    void DefineTracksCuts(TString cuts_option);
    void DefineV0Cuts(TString cuts_option);
    void DefineLambda1520Cuts(TString cuts_option);

    /* MC Generated */
    void ProcessMCGen();
    void GetDaughtersInfo(AliMCParticle* mcPart, Int_t& mcIdxNeutralDaughter, Int_t& mcIdxNegDaughter, Int_t& mcIdxPosDaughter,
                          TVector3& decayVertex);

    /* Tracks */
    void ProcessTracks();
    void PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode);
    Bool_t PassesTrackSelection(AliESDtrack* track);

    /* V0s and Lambdas1520 -- That Require True Info */
    void ProcessFindableV0s();
    void ProcessFindableLambdas1520();

    /* V0s */
    void KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos);
    Bool_t PassesV0Cuts(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfDaughterNeg, KFParticle kfDaughterPos, TLorentzVector lvV0,
                        TLorentzVector lvTrackNeg, TLorentzVector lvTrackPos);

    /* Lambda1520 */
    void KalmanLambda1520Finder();
    Bool_t PassesLambda1520Cuts(KFParticleMother kfLambda1520, TLorentzVector lvLambda1520, KFParticle kfAntiLambda, KFParticle kfKaonZeroShort,
                                KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg, KFParticle kfKaonZeroShortPos);

    /* Utilities */
    Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z, Double_t PV_X, Double_t PV_Y,
                          Double_t PV_Z);
    Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                             Double_t Pos_Py, Double_t Pos_Pz);
    Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Dau_Px, Double_t Dau_Py, Double_t Dau_Pz);
    Double_t CosinePointingAngle(TLorentzVector particles_mom, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);
    Double_t CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);

    /* Kalman Filter Utilities */
    KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge);
    KFVertex CreateKFVertex(const AliVVertex& vertex);
    KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis);

   private:
    /* Input options */
    Bool_t fIsMC;

    /* AliRoot objects */
    AliMCEvent* fMC;                // corresponding MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent* fESD;              // input event
    AliESDVertex* fPrimaryVertex;   // primary vertex
    AliPIDResponse* fPIDResponse;   // pid response object
    Double_t fMagneticField;        // magnetic field

    /* ROOT objects */
    TDatabasePDG fPDG;
    TList* fOutputListOfHists;

    /* Tracks Histograms */
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

    /*****************/
    /* V0 Histograms */
    /*
     key: `V0 pdg code`
     - `0` : MC V0s
     - `1` : MC V0s w/ charged decay products
     - `4` : signal MC V0s
     - `5` : signal MC V0s w/ charged decay products
     - `6` : findable true V0s
     - `7` : findable true secondary V0s
     - `8` : findable true signal V0s
     - `9` : found V0s
     - `10` : found true V0s
     - `11` : found true secondary V0s
     - `12` : found true signal V0s
    */
    std::unordered_map<Int_t, TH1F*> fHist_V0s_Bookkeep;
    // key: `stage, set, pdg, property`, value: histogram
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_V0s;
    // key: `set`, value: histogram
    std::map<TString, TH2F*> f2DHist_V0s_ArmenterosPodolanski;

    /********************************/
    /* (Anti)Lambdas1520 Histograms */
    /*
     key: `Lambda(1520) pdg code`
     - 0: n events
     - 1: n MC (anti)lambda(1520)
     - 2: n MC (anti)lambda(1520) -> lambda pi pi
     - 3: n MC (anti)lambda(1520) -> lambda pi+ pi-
    */
    std::unordered_map<Int_t, TH1F*> fHist_Lambdas1520_Bookkeep;
    // key: `stage, set, pdg code, property`, value: histogram
    std::map<std::tuple<TString, TString, Int_t, TString>, TH1F*> fHist_Lambdas1520;

    /**************/
    /* Containers */
    std::unordered_map<Int_t, Bool_t> isMcIdxSignal;  //

    std::unordered_map<Int_t, Int_t> getPdgCode_fromMcIdx;                 //
    std::unordered_map<Int_t, UInt_t> getLStarIdx_fromMcIdx;               //
    std::unordered_map<UInt_t, std::vector<Int_t>> getMcIdx_fromLStarIdx;  //

    std::unordered_map<Int_t, UInt_t> getLStarIdx_fromEsdIdx;               //
    std::unordered_map<UInt_t, std::vector<Int_t>> getEsdIdx_fromLStarIdx;  //

    std::vector<Int_t> esdIndicesOfAntiProtonTracks;  //
    std::vector<Int_t> esdIndicesOfPiMinusTracks;     //
    std::vector<Int_t> esdIndicesOfPiPlusTracks;      //

    std::unordered_map<Int_t, Int_t> getMcIdx_fromEsdIdx;  //

    std::unordered_map<Int_t, Bool_t> isEsdIdxSignal;  //

    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromLambdaIdx;      //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromLambdaIdx;      //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromAntiLambdaIdx;  //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromAntiLambdaIdx;  //

    std::vector<Bool_t> isLambdaIdxATrueV0;      //
    std::vector<Bool_t> isAntiLambdaIdxATrueV0;  //

    std::unordered_map<Int_t, Int_t> getMcIdx_fromLambdaIdx;      // must protect first with `isAntiLambdaIdxTrueV0`
    std::unordered_map<Int_t, Int_t> getMcIdx_fromAntiLambdaIdx;  // must protect first with `isAntiLambdaIdxTrueV0`

    std::vector<KFParticleMother> kfLambdas;      //
    std::vector<KFParticleMother> kfAntiLambdas;  //

    void ClearContainers();

    /** Cuts **/

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

    /* Lambda1520 candidates */
    Float_t kMin_LStar_Mass, kMax_LStar_Mass;
    Float_t kMin_LStar_Pt;
    Float_t kMax_LStar_Eta;
    Float_t kMax_LStar_Rapidity;
    Float_t kMin_LStar_DecayLength;
    Float_t kMin_LStar_Radius;
    Float_t kMin_LStar_CPAwrtPV, kMax_LStar_CPAwrtPV;
    Float_t kMin_LStar_DCAwrtPV, kMax_LStar_DCAwrtPV;
    Float_t kMax_LStar_DCAbtwV0s;
    Float_t kMax_LStar_DCAv0aSV;
    Float_t kMax_LStar_DCAv0bSV;
    Float_t kMax_LStar_DCAv0anegSV;
    Float_t kMax_LStar_DCAv0aposSV;
    Float_t kMax_LStar_DCAv0bnegSV;
    Float_t kMax_LStar_DCAv0bposSV;
    Float_t kMax_LStar_Chi2ndf;

    AliAnalysisTaskLambda1520Lpipi(const AliAnalysisTaskLambda1520Lpipi&);             // not implemented
    AliAnalysisTaskLambda1520Lpipi& operator=(const AliAnalysisTaskLambda1520Lpipi&);  // not implemented

    ClassDef(AliAnalysisTaskLambda1520Lpipi, 1);
};

#endif
