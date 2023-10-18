#ifndef AliAnalysisTaskLambda1520Lpipi_h
#define AliAnalysisTaskLambda1520Lpipi_h

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;

class AliAnalysisTaskLambda1520Lpipi : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskLambda1520Lpipi()
      : AliAnalysisTaskSE(),
        fIsMC(0),
        fOutputListOfTrees(0),
        fOutputListOfHists(0),
        kMassNeutron(0),
        kMassProton(0),
        kMassPion(0),
        kMassKaon(0),
        fMC(0),
        fESD(0),
        fPIDResponse(0),
        fPrimaryVertex(0),
        fTree(0),
        fHist_Bookkeeper(0),
        fHist_TrueLambda_Pt(0),
        fHist_TrueLambda_Pz(0),
        fHist_TrueLambda_Radius(0),
        fHist_TrueAntiLambda_Pt(0),
        fHist_TrueAntiLambda_Pz(0),
        fHist_TrueAntiLambda_Radius(0) {}
  AliAnalysisTaskLambda1520Lpipi(const char* name, Bool_t IsMC)
      : AliAnalysisTaskSE(name),
        fIsMC(IsMC),
        fOutputListOfTrees(0),
        fOutputListOfHists(0),
        kMassNeutron(0),
        kMassProton(0),
        kMassPion(0),
        kMassKaon(0),
        fMC(0),
        fESD(0),
        fPIDResponse(0),
        fPrimaryVertex(0),
        fTree(0),
        fHist_Bookkeeper(0),
        fHist_TrueLambda_Pt(0),
        fHist_TrueLambda_Pz(0),
        fHist_TrueLambda_Radius(0),
        fHist_TrueAntiLambda_Pt(0),
        fHist_TrueAntiLambda_Pz(0),
        fHist_TrueAntiLambda_Radius(0) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    CheckForInputErrors();
  }
  ~AliAnalysisTaskLambda1520Lpipi() {
    if (fOutputListOfTrees) {
      delete fOutputListOfTrees;
    }
    if (fOutputListOfHists) {
      delete fOutputListOfHists;
    }
  }

 public:
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option) { return; }

 public:
  virtual void InitPDGMasses();
  virtual void CheckForInputErrors();

  virtual void ProcessMCGen(std::set<Int_t>&, std::set<Int_t>&);
  virtual void ProcessMCRec(std::set<Int_t>, std::vector<Int_t>&);
  virtual void ProcessTrueV0s_KF(std::vector<Int_t>, std::map<Int_t, Int_t>&, std::map<Int_t, Int_t>&, std::vector<AliKFParticle>&,
                                 std::map<Int_t, std::vector<Int_t>>&, std::vector<AliKFParticle>&);
  virtual void Lambda1520Finder(std::vector<AliKFParticle>, std::map<Int_t, std::vector<Int_t>>, std::vector<AliKFParticle>,
                                std::vector<AliKFParticle>&);
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

 private:
  // (input options)
  Bool_t fIsMC;

  Double_t kMassNeutron;
  Double_t kMassProton;
  Double_t kMassPion;
  Double_t kMassKaon;

 private:
  // AliRoot objects
  AliMCEvent* fMC;               // corresponding MC event
  AliESDEvent* fESD;             // input event
  AliPIDResponse* fPIDResponse;  // pid response object
  AliESDVertex* fPrimaryVertex;  // primary vertex

  TList* fOutputListOfTrees;
  TList* fOutputListOfHists;

 private:
  TTree* fTree;

 private:
  /** Histogram to keep count of
  - 0: n events
  - 1: n gen lambda(1520)
  - 2: n gen lambda(1520) -> lambda pi pi
  - 3: n gen lambda(1520) -> lambda pi+ pi-
  - 4: n gen anti-lambda(1520)
  - 5: n gen anti-lambda(1520) -> anti-lambda pi pi
  - 6: n gen anti-lambda(1520) -> anti-lambda pi+ pi-
  */
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
