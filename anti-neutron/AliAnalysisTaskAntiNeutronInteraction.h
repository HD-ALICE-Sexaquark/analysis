#ifndef AliAnalysisTaskAntiNeutronInteraction_H
#define AliAnalysisTaskAntiNeutronInteraction_H

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

class AliAnalysisTaskAntiNeutronInteraction : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskAntiNeutronInteraction();
    AliAnalysisTaskAntiNeutronInteraction(const char* name, Bool_t IsMC);
    virtual ~AliAnalysisTaskAntiNeutronInteraction();

   public:
    virtual void UserCreateOutputObjects();
    void PrepareTracksHistograms();
    void PrepareAntiNeutronHistograms();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option) { return; }

   public:
    /* Initialization */
    void Initialize();

   public:
    /* MC Generated */
    void ProcessMCGen();

   public:
    /* Cuts */
    void DefineTracksCuts(TString cuts_option);
    void DefineAntiNeutronCuts(TString cuts_option);
    /*  */
    Bool_t PassesAntiNeutronCuts(KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark, KFParticle kfAntiLambda,
                                 KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg,
                                 KFParticle kfKaonZeroShortPos);

   public:
    /* Tracks */
    void ProcessTracks();
    Bool_t PassesTrackSelection(AliESDtrack* track);
    void PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode);

   public:
    /* Anti-Neutrons -- That Require True Info */
    void ProcessFindableAntiNeutrons();

   public:
    /* Anti-Neutrons -- Kalman Filter */
    void KalmanAntiNeutronFinder();

   public:
    /* Utilities */
    Double_t CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);
    Double_t CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX, Double_t refPointY,
                                 Double_t refPointZ);
    Double_t Calculate_ArmAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz, Double_t Pos_Px,
                                Double_t Pos_Py, Double_t Pos_Pz);
    Double_t Calculate_ArmPt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz);
    Double_t LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z, Double_t refPointX,
                          Double_t refPointY, Double_t refPointZ);
    Double_t SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]);
    Double_t SquaredDistanceBetweenLines(const Double_t* t, Double_t pos0[], Double_t dir0[], Double_t pos1[], Double_t dir1[]);
    Double_t Calculate_TwoLinesDCA_v1(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0, TVector3& PCA1);
    Double_t Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0, TVector3& PCA1);

   public:
    /* Kalman Filter Utilities */
    KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge);
    KFVertex CreateKFVertex(const AliVVertex& vertex);
    KFParticle TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis);

   private:
    /* Input options */
    Bool_t fIsMC;                               //
    std::vector<Int_t> fProductsPDG;            //
    std::vector<Int_t> fFinalStateProductsPDG;  //
    Int_t fStruckNucleonPDG;                    // PDG Code of the struck nucleon, could ber: 2112 o 2212

   private:
    /* AliRoot objects */
    AliMCEvent* fMC;                // MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent* fESD;              // reconstructed event
    AliPIDResponse* fPIDResponse;   // pid response object
    AliESDVertex* fPrimaryVertex;   // primary vertex
    Double_t fMagneticField;        // magnetic field

   private:
    /* ROOT objects */
    TDatabasePDG fPDG;
    TList* fOutputListOfHists;

   private:
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

   private:
    /* Anti-Neutron Histograms */
    // key: `stage, set, property`, value: histogram
    std::map<std::tuple<TString, TString, TString>, TH1F*> fHist_AntiNeutrons;

   private:
    /* Anti-Neutron Histograms */
    TH1F* fHist_All_AntiNeutron_Pt;           //!
    TH1F* fHist_Signal_AntiNeutron_Pt;        //!
    TH1F* fHist_Signal_AntiNeutron_Radius;    //!
    TH1F* fHist_True_AntiNeutron_NDaughters;  //!
    TH1F* fHist_True_AntiNeutron_Bookkeep;    //!

   private:
    /* Containers -- vectors and hash tables */

    std::unordered_map<Int_t, Bool_t> isMcIdxSignal;        //
    std::unordered_map<Int_t, Int_t> getPdgCode_fromMcIdx;  //

    std::unordered_map<Int_t, UInt_t> getReactionIdx_fromEsdIdx;               //
    std::unordered_map<UInt_t, std::vector<Int_t>> getEsdIdx_fromReactionIdx;  //

    std::unordered_map<Int_t, Bool_t> doesEsdIdxPassTrackSelection;  //

    std::vector<Int_t> esdIndicesOfAntiProtonTracks;  //
    std::vector<Int_t> esdIndicesOfPosKaonTracks;     //
    std::vector<Int_t> esdIndicesOfPiMinusTracks;     //
    std::vector<Int_t> esdIndicesOfPiPlusTracks;      //

    std::unordered_map<Int_t, Int_t> getMcIdx_fromEsdIdx;  //

    std::unordered_map<Int_t, Bool_t> isEsdIdxSignal;  //

    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromAntiLambdaIdx;     //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromAntiLambdaIdx;     //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfNegDau_fromKaonZeroShortIdx;  //
    std::unordered_map<Int_t, Int_t> getEsdIdxOfPosDau_fromKaonZeroShortIdx;  //

    void ClearContainers();

   private:
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

    /* Anti-Sexaquark candidates */
    Float_t kMin_AN_Mass, kMax_AN_Mass;
    Float_t kMin_AN_Pt;
    Float_t kMax_AN_Eta;
    Float_t kMax_AN_Rapidity;
    Float_t kMin_AN_DecayLength;
    Float_t kMin_AN_Radius;
    Float_t kMin_AN_CPAwrtPV, kMax_AN_CPAwrtPV;
    Float_t kMin_AN_DCAwrtPV, kMax_AN_DCAwrtPV;
    Float_t kMax_AN_DCAbtwV0s;
    Float_t kMax_AN_DCAv0aSV;
    Float_t kMax_AN_DCAv0bSV;
    Float_t kMax_AN_DCAv0anegSV;
    Float_t kMax_AN_DCAv0aposSV;
    Float_t kMax_AN_DCAv0bnegSV;
    Float_t kMax_AN_DCAv0bposSV;
    Float_t kMax_AN_Chi2ndf;

    AliAnalysisTaskAntiNeutronInteraction(const AliAnalysisTaskAntiNeutronInteraction&);             // not implemented
    AliAnalysisTaskAntiNeutronInteraction& operator=(const AliAnalysisTaskAntiNeutronInteraction&);  // not implemented

    ClassDef(AliAnalysisTaskAntiNeutronInteraction, 1);
};

#endif
