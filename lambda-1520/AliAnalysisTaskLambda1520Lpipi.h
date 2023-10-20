#ifndef AliAnalysisTaskLambda1520Lpipi_h
#define AliAnalysisTaskLambda1520Lpipi_h

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;

class AliAnalysisTaskLambda1520Lpipi : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskLambda1520Lpipi();
  AliAnalysisTaskLambda1520Lpipi(const char* name, Bool_t IsMC);
  virtual ~AliAnalysisTaskLambda1520Lpipi();

 public:
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option) { return; }

 public:
  virtual void CheckForInputErrors();

  void DefineCuts(TString cuts_option);

  virtual void ProcessMCGen(std::set<Int_t>&, std::set<Int_t>&);
  virtual void ProcessTracks(std::set<Int_t> Indices_MCGen_FS_Signal,                                        // (pending)
                             std::vector<Int_t>& idxPiPlusTracks, std::vector<Int_t>& idxPiMinusTracks,      //
                             std::vector<Int_t>& idxKaonPlusTracks, std::vector<Int_t>& idxKaonMinusTracks,  //
                             std::vector<Int_t>& idxProtonTracks, std::vector<Int_t>& idxAntiProtonTracks);
  virtual void ReconstructV0s_KF(std::vector<Int_t> idxNegativeTracks,  //
                                 std::vector<Int_t> idxPositiveTracks,  //
                                 Int_t pdgTrackNeg, Int_t pdgTrackPos,  //
                                 std::vector<KFParticle>& kfV0s,        //
                                 std::vector<std::vector<Int_t>>& idxDaughters);
  virtual void Lambda1520Finder(std::vector<KFParticle> kfFirstV0s,                    //
                                std::vector<std::vector<Int_t>> idxFirstV0Daughters,   //
                                std::vector<KFParticle> kfSecondV0s,                   //
                                std::vector<std::vector<Int_t>> idxSecondV0Daughters,  //
                                std::vector<Int_t> pdgDaughters,                       //
                                std::vector<KFParticle>& kfLambdas1520);
  void SexaquarkFinder_ChannelA(std::vector<KFParticle> kfFirstV0s,                    //
                                std::vector<std::vector<Int_t>> idxFirstV0Daughters,   //
                                std::vector<KFParticle> kfSecondV0s,                   //
                                std::vector<std::vector<Int_t>> idxSecondV0Daughters,  //
                                std::vector<Int_t> pdgDaughters,                       //
                                std::vector<KFParticle>& kfAntiSexaquarks);
  virtual void ProcessOfficialV0s(std::vector<Int_t>);

  virtual void ProcessTrueV0s();

  // Utilities
  virtual Float_t MCRec_GetImpactParameter(AliESDtrack* track);

  Double_t Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                  Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                  Double_t PV_X, Double_t PV_Y, Double_t PV_Z);
  Bool_t Preoptimize(const AliExternalTrackParam*, AliExternalTrackParam*, Double_t*, Double_t*, const Double_t);
  void GetHelixCenter(const AliExternalTrackParam*, Double_t[2], const Double_t);
  KFParticle CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge);
  KFVertex CreateKFVertex(const AliVVertex& vertex);
  Double_t ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
                           Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz,  //
                           Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz);
  Double_t ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                        Double_t Dau_Px, Double_t Dau_Py, Double_t Dau_Pz);
  Double_t CosinePointingAngle(TLorentzVector particles_mom, Double_t X, Double_t Y, Double_t Z,  //
                               Double_t refPointX, Double_t refPointY, Double_t refPointZ);

 private:
  /* Input options */
  Bool_t fIsMC;

 private:
  /* Cuts */
  Float_t kMaxNSigma_Pion;
  Float_t kMaxNSigma_Kaon;
  Float_t kMaxNSigma_Proton;

 private:
  /* AliRoot objects */
  AliMCEvent* fMC;               // corresponding MC event
  AliESDEvent* fESD;             // input event
  AliPIDResponse* fPIDResponse;  // pid response object
  AliESDVertex* fPrimaryVertex;  // primary vertex
  Double_t fMagneticField;       // magnetic field

  TList* fOutputListOfTrees;
  TList* fOutputListOfHists;

 private:
  /* ROOT objects */
  TDatabasePDG fPDG;
  TTree* fTree;
  // Histogram to keep count of
  // - 0: n events
  // - 1: n gen lambda(1520)
  // - 2: n gen lambda(1520) -> lambda pi pi
  // - 3: n gen lambda(1520) -> lambda pi+ pi-
  // - 4: n gen anti-lambda(1520)
  // - 5: n gen anti-lambda(1520) -> anti-lambda pi pi
  // - 6: n gen anti-lambda(1520) -> anti-lambda pi+ pi-
  TH1I* fHist_Bookkeeper;
  TH1F* fHist_TrueLambda_Pt;
  TH1F* fHist_TrueLambda_Pz;
  TH1F* fHist_TrueLambda_Radius;
  TH1F* fHist_TrueAntiLambda_Pt;
  TH1F* fHist_TrueAntiLambda_Pz;
  TH1F* fHist_TrueAntiLambda_Radius;

 private:
  AliAnalysisTaskLambda1520Lpipi(const AliAnalysisTaskLambda1520Lpipi&);             // not implemented
  AliAnalysisTaskLambda1520Lpipi& operator=(const AliAnalysisTaskLambda1520Lpipi&);  // not implemented

  ClassDef(AliAnalysisTaskLambda1520Lpipi, 1);
};

#endif
