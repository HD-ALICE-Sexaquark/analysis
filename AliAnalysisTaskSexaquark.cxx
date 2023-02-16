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

#include "AliAnalysisTaskSexaquark.h"

#define N_SIGMA 3.0

class AliAnalysisTaskSexaquark;

//_____________________________________________________________________________
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fHasSexaquark(0),
      fSourceOfV0s(),
      fSimulationSet(0),
      fOutputListOfHists(0),
      fOutputListOfTrees(0),
      kMassNeutron(0),
      kMassProton(0),
      kMassPion(0),
      kMassKaon(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0),
      fMap_DuplicatedTracks(),
      fIndex_NonRelevantMC(),
      fHist_PDG(0),
      // fHist_Kp_P(0),
      // fHist_Kp_Pt(0),
      // fHist_Pip_P(0),
      // fHist_Pip_Pt(0),
      // fHist_Pim_P(0),
      // fHist_Pim_Pt(0),
      // fHist_AntiP_P(0),
      // fHist_AntiP_Pt(0),
      // fHist_K0_P(0),
      // fHist_K0_Pt(0),
      // fHist_K0_M(0),
      // fHist_AntiL_P(0),
      // fHist_AntiL_Pt(0),
      // fHist_AntiL_M(0),
      // fHist_Sexaquark_P(0),
      // fHist_Sexaquark_Pt(0),
      // fHist_Sexaquark_M(0),
      fHistMC_PDG(0),
      // fHistMC_Kp_P(0),
      // fHistMC_Kp_Pt(0),
      // fHistMC_Pip_P(0),
      // fHistMC_Pip_Pt(0),
      // fHistMC_Pim_P(0),
      // fHistMC_Pim_Pt(0),
      // fHistMC_AntiP_P(0),
      // fHistMC_AntiP_Pt(0),
      // fHistMC_AntiN_P(0),
      // fHistMC_AntiN_Pt(0),
      // fHistMC_Gamma_P(0),
      // fHistMC_Gamma_Pt(0),
      // fHistMC_AntiL_P(0),
      // fHistMC_AntiL_Pt(0),
      // fHistMC_AntiL_M(0),
      // fHistMC_AntiL_CPA(0),
      // fHistMC_K0_P(0),
      // fHistMC_K0_Pt(0),
      // fHistMC_K0_M(0),
      // fHistMC_K0_CPA(0),
      // fHistMC_Sexaquark_P(0),
      // fHistMC_Sexaquark_Pt(0),
      // fHistMC_Sexaquark_M(0),
      fTree(0),
      fIndexer_fMC_to_fEvent(),
      fIndexer_fTrack_to_fEvent(),
      fSexaquarkDaughters(),
      fEvent() {
  //
  // default constructor
  // this is used by root for IO purposes, it needs to remain empty
  //
}

//_____________________________________________________________________________
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name, Bool_t IsMC, Bool_t HasSexaquark, TString SourceOfV0s,
                                                   char SimulationSet)
    : AliAnalysisTaskSE(name),
      fSimulationSet(SimulationSet),
      fIsMC(IsMC),
      fSourceOfV0s(SourceOfV0s),
      fHasSexaquark(HasSexaquark),
      fOutputListOfHists(0),
      fOutputListOfTrees(0),
      kMassNeutron(0),
      kMassProton(0),
      kMassPion(0),
      kMassKaon(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0),
      fMap_DuplicatedTracks(),
      fIndex_NonRelevantMC(),
      fHist_PDG(0),
      // fHist_Kp_P(0),
      // fHist_Kp_Pt(0),
      // fHist_Pip_P(0),
      // fHist_Pip_Pt(0),
      // fHist_Pim_P(0),
      // fHist_Pim_Pt(0),
      // fHist_AntiP_P(0),
      // fHist_AntiP_Pt(0),
      // fHist_K0_P(0),
      // fHist_K0_Pt(0),
      // fHist_K0_M(0),
      // fHist_AntiL_P(0),
      // fHist_AntiL_Pt(0),
      // fHist_AntiL_M(0),
      // fHist_Sexaquark_P(0),
      // fHist_Sexaquark_Pt(0),
      // fHist_Sexaquark_M(0),
      fHistMC_PDG(0),
      // fHistMC_Kp_P(0),
      // fHistMC_Kp_Pt(0),
      // fHistMC_Pip_P(0),
      // fHistMC_Pip_Pt(0),
      // fHistMC_Pim_P(0),
      // fHistMC_Pim_Pt(0),
      // fHistMC_AntiP_P(0),
      // fHistMC_AntiP_Pt(0),
      // fHistMC_AntiN_P(0),
      // fHistMC_AntiN_Pt(0),
      // fHistMC_Gamma_P(0),
      // fHistMC_Gamma_Pt(0),
      // fHistMC_AntiL_P(0),
      // fHistMC_AntiL_Pt(0),
      // fHistMC_AntiL_M(0),
      // fHistMC_AntiL_CPA(0),
      // fHistMC_K0_P(0),
      // fHistMC_K0_Pt(0),
      // fHistMC_K0_M(0),
      // fHistMC_K0_CPA(0),
      // fHistMC_Sexaquark_P(0),
      // fHistMC_Sexaquark_Pt(0),
      // fHistMC_Sexaquark_M(0),
      fTree(0),
      fIndexer_fMC_to_fEvent(),
      fIndexer_fTrack_to_fEvent(),
      fSexaquarkDaughters(),
      fEvent() {
  //
  // constructor
  //
  DefineInput(0, TChain::Class());  // define the input of the analysis: in this case we take a 'chain' of events
                                    // this chain is created by the analysis manager, so no need to worry about it,
                                    // it does its work automatically
  DefineOutput(1, TList::Class());  // define the output of the analysis: in this case it's a list of histograms
                                    // you can add more output objects by calling DefineOutput(2, classname::Class())
                                    // if you add more output objects, make sure to call PostData for all of them, and to
                                    // make changes to your AddTask macro!
  DefineOutput(2, TList::Class());

  CheckForInputErrors();
}

//_____________________________________________________________________________
AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() {
  //
  // destructor
  //
  // at the end of your task, it is deleted from memory by calling this function
  if (fOutputListOfHists) {
    delete fOutputListOfHists;
  }
  if (fOutputListOfTrees) {
    delete fOutputListOfTrees;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {
  //
  // create output objects
  // this function is called ONCE at the start of your analysis (RUNTIME)
  //

  // add the PID routine
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  if (!man) {
    AliFatal("ERROR: AliAnalysisManager couldn't be found.");
  }

  AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
  if (!inputHandler) {
    AliFatal("ERROR: AliESDInputHandler couldn't be found.");
  }

  fPIDResponse = inputHandler->GetPIDResponse();

  /*** Initializing ***/

  AliInfo("Initializing AliAnalysisTaskSexaquark...");
  AliInfo("INPUT OPTIONS:");
  AliInfoF(" >> IsMC           = %i", (Int_t)fIsMC);
  AliInfoF(" >> HasSexaquark   = %i", (Int_t)fHasSexaquark);
  AliInfoF(" >> Source of V0s  = %s", fSourceOfV0s.Data());
  AliInfoF(" >> Simulation Set = %c", fSimulationSet);
  InitPDGMasses();

  // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written
  // to the output file

  /*** List of Histograms ***/

  fOutputListOfHists = new TList();
  fOutputListOfHists->SetOwner(kTRUE);  // the list is owner of all objects it contains and will delete them if requested

  // MC Rec. Hists
  fHist_PDG = new TH1F("Hist_PDG", "Hist_PDG", 400, -3500, 500);
  fOutputListOfHists->Add(fHist_PDG);
  /*
  fHist_Kp_P = new TH1F("Hist_Kp_P", "Hist_Kp_P", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Kp_P);

  fHist_Kp_Pt = new TH1F("Hist_Kp_Pt", "Hist_Kp_Pt", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Kp_Pt);

  fHist_Pip_P = new TH1F("Hist_Pip_P", "Hist_Pip_P", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Pip_P);

  fHist_Pip_Pt = new TH1F("Hist_Pip_Pt", "Hist_Pip_Pt", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Pip_Pt);

  fHist_Pim_P = new TH1F("Hist_Pim_P", "Hist_Pim_P", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Pim_P);

  fHist_Pim_Pt = new TH1F("Hist_Pim_Pt", "Hist_Pim_Pt", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Pim_Pt);

  fHist_AntiP_P = new TH1F("Hist_AntiP_P", "Hist_AntiP_P", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_AntiP_P);

  fHist_AntiP_Pt = new TH1F("Hist_AntiP_Pt", "Hist_AntiP_Pt", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_AntiP_Pt);

  fHist_K0_P = new TH1F("Hist_K0_P", "Hist_K0_P", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_K0_P);

  fHist_K0_Pt = new TH1F("Hist_K0_Pt", "Hist_K0_Pt", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_K0_Pt);

  fHist_K0_M = new TH1F("Hist_K0_M", "Hist_K0_M", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_K0_M);

  fHist_AntiL_P = new TH1F("Hist_AntiL_P", "Hist_AntiL_P", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_AntiL_P);

  fHist_AntiL_Pt = new TH1F("Hist_AntiL_Pt", "Hist_AntiL_Pt", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_AntiL_Pt);

  fHist_AntiL_M = new TH1F("Hist_AntiL_M", "Hist_AntiL_M", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_AntiL_M);

  fHist_Sexaquark_P = new TH1F("Hist_Sexaquark_P", "Hist_Sexaquark_P", 8, 0., 8.);
  fOutputListOfHists->Add(fHist_Sexaquark_P);

  fHist_Sexaquark_Pt = new TH1F("Hist_Sexaquark_Pt", "Hist_Sexaquark_Pt", 5, 0., 5.);
  fOutputListOfHists->Add(fHist_Sexaquark_Pt);

  fHist_Sexaquark_M = new TH1F("Hist_Sexaquark_M", "Hist_Sexaquark_M", 100, 0., 5.);
  fOutputListOfHists->Add(fHist_Sexaquark_M);
  */

  // MC Gen.

  fHistMC_PDG = new TH1F("HistMC_PDG", "PDG Code of Generated Particles (after all decays)", 400, -3500, 500);
  fOutputListOfHists->Add(fHistMC_PDG);
  /*
  fHistMC_Kp_P = new TH1F("HistMC_Kp_P", "Momentum of Generated K^{+}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Kp_P);

  fHistMC_Kp_Pt = new TH1F("HistMC_Kp_Pt", "Pt of Generated K^{+}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Kp_Pt);

  fHistMC_Pip_P = new TH1F("HistMC_Pip_P", "Momentum of Generated #pi^{+}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Pip_P);

  fHistMC_Pip_Pt = new TH1F("HistMC_Pip_Pt", "Pt of Generated #pi^{+}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Pip_Pt);

  fHistMC_Pim_P = new TH1F("HistMC_Pim_P", "Momentum of Generated #pi^{-}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Pim_P);

  fHistMC_Pim_Pt = new TH1F("HistMC_Pim_Pt", "Pt of Generated #pi^{-}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Pim_Pt);

  fHistMC_AntiP_P = new TH1F("HistMC_AntiP_P", "Momentum of Generated anti-protons", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiP_P);

  fHistMC_AntiP_Pt = new TH1F("HistMC_AntiP_Pt", "Pt of Generated anti-protons", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiP_Pt);

  fHistMC_AntiN_P = new TH1F("HistMC_AntiN_P", "Momentum of Generated anti-neutrons", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiN_P);

  fHistMC_AntiN_Pt = new TH1F("HistMC_AntiN_Pt", "Pt of Generated anti-neutrons", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiN_Pt);

  fHistMC_Gamma_P = new TH1F("HistMC_Gamma_P", "Momentum of Generated #gamma", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Gamma_P);

  fHistMC_Gamma_Pt = new TH1F("HistMC_Gamma_Pt", "Pt of Generated #gamma", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Gamma_Pt);

  fHistMC_AntiL_P = new TH1F("HistMC_AntiL_P", "Momentum of Generated anti-#Lambda", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiL_P);

  fHistMC_AntiL_Pt = new TH1F("HistMC_AntiL_Pt", "Pt of Generated anti-#Lambda", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiL_Pt);

  fHistMC_AntiL_M = new TH1F("HistMC_AntiL_M", "Inv. Mass of Generated anti-#Lambda", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_AntiL_M);

  fHistMC_AntiL_CPA = new TH1F("HistMC_AntiL_CPA", "CPA Generated #bar{#Lambda}", 200, -1., 1.);
  fOutputListOfHists->Add(fHistMC_AntiL_CPA);

  fHistMC_K0_P = new TH1F("HistMC_K0_P", "Momentum of Generated K^{0}_{S}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_K0_P);

  fHistMC_K0_Pt = new TH1F("HistMC_K0_Pt", "Pt of Generated K^{0}_{S}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_K0_Pt);

  fHistMC_K0_M = new TH1F("HistMC_K0_M", "Inv. Mass of Generated K^{0}_{S}", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_K0_M);

  fHistMC_K0_CPA = new TH1F("HistMC_K0_CPA", "CPA Generated K^{0}_{S}", 200, -1., 1.);
  fOutputListOfHists->Add(fHistMC_K0_CPA);

  fHistMC_Sexaquark_P = new TH1F("HistMC_Sexaquark_P", "Momentum of Generated Sexaquark", 8, 0., 8.);
  fOutputListOfHists->Add(fHistMC_Sexaquark_P);

  fHistMC_Sexaquark_Pt = new TH1F("HistMC_Sexaquark_Pt", "Pt of Generated Sexaquark", 5, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Sexaquark_Pt);

  fHistMC_Sexaquark_M = new TH1F("HistMC_Sexaquark_M", "Inv. Mass of Generated Sexaquark", 100, 0., 5.);
  fOutputListOfHists->Add(fHistMC_Sexaquark_M);
  */

  PostData(1, fOutputListOfHists);  // postdata will notify the analysis manager of changes / updates to the
                                    // fOutputList object. the manager will in the end take care of writing your output to file
                                    // so it needs to know what's in the output

  /*** List of Trees ***/

  fOutputListOfTrees = new TList();
  fOutputListOfTrees->SetOwner(kTRUE);

  fTree = new TTree("Events", "Tree of Events");
  SetBranches();
  fOutputListOfTrees->Add(fTree);

  PostData(2, fOutputListOfTrees);
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::InitPDGMasses() {
  //
  // Set mass-at-rest of different particles
  //
  kMassNeutron = 0.939565;
  AliInfoF("kMassNeutron = %.6f", kMassNeutron);

  kMassProton = 0.938272;
  AliInfoF("kMassProton = %.6f", kMassProton);

  kMassPion = 0.139570;
  AliInfoF("kMassPion = %.6f", kMassPion);

  kMassKaon = 0.493677;
  AliInfoF("kMassKaon = %.6f", kMassKaon);
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {
  //
  // the body of the class, this function is called per each event
  // the manager will take care of reading the events from file, and with the static function InputEvent()
  // you have access to the current event
  // once you return from UserExec(), the manager will retrieve the next event from the chain
  //

  // get MC event
  fMC = MCEvent();

  // (1) process MC tree
  if (fMC) {
    MCGen_FirstProcess();
  } else {
    AliFatal("ERROR: AliMCEvent couldn't be found.");
  }

  // load reconstructed event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());

  // (2) process reconstructed tracks
  if (fESD) {
    MCRec_FirstProcess();
  } else {
    AliFatal("ERROR: AliESDEvent couldn't be found.");
  }

  // (3) add non-relevant MC particles, that were rec. as relevant
  if (fMC) {
    MCGen_SecondProcess();
  }

  // (4) label duplicated tracks
  if (fESD) {
    MCRec_SecondProcess();
  }

  // (5) assign daughters indices
  if (fMC) {
    MCGen_ThirdProcess();
  }

  // (6) V0 finder
  if (fSourceOfV0s == "official") {
    V0_ProcessESD();
  }
  if (fSourceOfV0s == "custom") {
    V0_ProcessCustom();
  }
  if (fSourceOfV0s == "true") {
    V0_ProcessTrue();
  }

  // (tree operations)
  // after all variables have been assigned, fill tree
  fTree->Fill();

  // stream the results the analysis of this event to the output manager
  // which will take care of writing it to a file
  PostData(1, fOutputListOfHists);
  PostData(2, fOutputListOfTrees);

  /*** End of Event: Clear STD Structures ***/

  // clear maps
  fMap_DuplicatedTracks.clear();

  // clear indexers
  fIndexer_fMC_to_fEvent.clear();
  fIndexer_fTrack_to_fEvent.clear();

  // (artificial sexaquark)
  fSexaquarkDaughters.clear();

  // (non-relevant, but rec. as relevant)
  fIndex_NonRelevantMC.clear();

  // (tree operations)
  // clear tree vars
  ClearEvent();
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::CheckForInputErrors() {
  //
  // Search for possible contradictions
  //
  if (fIsMC == kFALSE) {
    if (fHasSexaquark == kTRUE) {
      AliFatal("ERROR: data cannot have a MC anti-sexaquark inserted.");
    }
    if (fSourceOfV0s == "true") {
      AliFatal("ERROR: data cannot have access to true MC information.");
    }
  }
  if (fSourceOfV0s == "true" || fSourceOfV0s == "official" || fSourceOfV0s == "custom") {
    // nothing
  } else {
    AliFatal("ERROR: source of V0s must be true (available only for sim), official or custom.");
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_FirstProcess() {
  //
  // Search for the relevant MC particles in a single event
  //

  // define variables & objects
  AliMCParticle* mcPart;
  Bool_t isRelevant;
  Bool_t isSignal;
  Int_t generation;

  // (artificial anti-sexaquark)
  TVector3 antilambda;
  TVector3 neutral_kaon;

  // (debug)
  // AliInfo("MONTE CARLO TREE");
  // AliInfo("");
  // AliInfo(" Index    PID Mother Status  Label    Gen       Px       Py       Pz       Xv       Yv       Zv");

  // loop over MC gen. particles
  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

    mcPart = (AliMCParticle*)fMC->GetTrack(i);

    isRelevant = TMath::Abs(mcPart->PdgCode()) == 211;      // is neg. / pos. pion
    isRelevant = isRelevant || mcPart->PdgCode() == 321;    // is pos. kaon
    isRelevant = isRelevant || mcPart->PdgCode() == 310;    // is neutral kaon short
    isRelevant = isRelevant || mcPart->PdgCode() == -3122;  // is anti-lambda
    isRelevant = isRelevant || mcPart->PdgCode() == -2212;  // is anti-proton
    isRelevant = isRelevant || mcPart->PdgCode() == -3312;  // is Xi plus

    // IMPORTANT: mcPart->Label() != index... the real index is i!!
    if (isRelevant) {

      // verify source and generation of this particle
      isSignal = kFALSE;
      generation = 3;
      MCGen_VerifySourceAndGeneration(mcPart, generation, isSignal);

      // (artificial anti-sexaquark)
      if (isSignal && generation == 1) {
        if (mcPart->PdgCode() == -3122) {
          antilambda.SetXYZ(mcPart->Px(), mcPart->Py(), mcPart->Pz());
        }
        if (mcPart->PdgCode() == 310) {
          neutral_kaon.SetXYZ(mcPart->Px(), mcPart->Py(), mcPart->Pz());
        }
      }

      // (debug)
      // AliInfoF("%6i %6i %6i %6i %6i %6i %8.4f %8.4f %8.4f %8.3f %8.3f %8.3f", i, mcPart->PdgCode(), mcPart->GetMother(),
      //  mcPart->MCStatusCode(), isSignal, generation, mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->Xv(), mcPart->Yv(),
      //  mcPart->Zv());

      // (hist operations)
      // FillMCHistograms_Products(mcPart);

      // (tree operations)
      // before pushing back, it's important to correctly index this
      // the "+1" is because of the additional (artificial sexaquark)
      fIndexer_fMC_to_fEvent[i] = fEvent.N_MCGen + 1 * (Int_t)fHasSexaquark;
      MCGen_PushBack(mcPart, generation, isSignal);
    }
  }  // end of loop over MC particles

  // (artificial anti-sexaquark)
  if (fHasSexaquark) {
    MCGen_InsertAntiSexaquark(antilambda, neutral_kaon);
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_SecondProcess() {
  //
  // Add the non-relevant MC particles that were reconstructed as relevant
  //

  // define variables & objects
  AliMCParticle* mcPart;
  Bool_t isSignal;
  Int_t generation;

  // (debug)
  // AliInfo("MONTE CARLO TREE");
  // AliInfo("");
  // AliInfo(" Index    PID Mother Status  Label    Gen       Px       Py       Pz       Xv       Yv       Zv");

  // loop over indices of non-relevant MC particles
  for (Int_t index : fIndex_NonRelevantMC) {

    // load MC particle
    mcPart = (AliMCParticle*)fMC->GetTrack(index);

    // verify source and generation of this particle
    isSignal = kFALSE;
    generation = 3;
    MCGen_VerifySourceAndGeneration(mcPart, generation, isSignal);

    // (debug)
    // AliInfoF("%6i %6i %6i %6i %6i %6i %8.4f %8.4f %8.4f %8.3f %8.3f %8.3f", index, mcPart->PdgCode(), mcPart->GetMother(),
    //  mcPart->MCStatusCode(), isSignal, generation, mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->Xv(), mcPart->Yv(),
    //  mcPart->Zv());

    // (tree operations)
    // (artificial anti-sexaquark)
    // in this case, the +1 was already applied in "MC_InsertAntiSexaquark()"
    fIndexer_fMC_to_fEvent[index] = fEvent.N_MCGen;
    MCGen_PushBack(mcPart, generation, isSignal);
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_ThirdProcess() {
  //
  // Solve indexing to daughters
  // -2 : daughter is not relevant
  // -1 : particle has no daughter
  //  X : particle has a relevant daughter, located at index X
  //

  // define indices & objects
  Int_t fMC_index;
  Int_t fEvent_index;
  AliMCParticle* mcPart;

  Int_t idx_first_daughter;
  Int_t idx_last_daughter;

  // resize vectors
  fEvent.MC_FirstDau.resize(fEvent.MC_Mother.size());
  fEvent.MC_LastDau.resize(fEvent.MC_Mother.size());

  // (debug)
  // AliInfoF("first_dau.size = %i, last_dau.size = %i", (Int_t)fEvent.MC_FirstDau.size(), (Int_t)fEvent.MC_FirstDau.size());
  // AliInfo("   fMC_index fEvent_index    first_dau     last_dau");

  // loop over map
  for (auto it = fIndexer_fMC_to_fEvent.begin(); it != fIndexer_fMC_to_fEvent.end(); it++) {

    fMC_index = it->first;
    fEvent_index = it->second;
    mcPart = (AliMCParticle*)fMC->GetTrack(fMC_index);

    // default value: particle has no daughters
    idx_first_daughter = -1;
    idx_last_daughter = -1;

    // check if particle has daughters
    if (mcPart->GetNDaughters()) {
      // has daughters, but they weren't preselected in the fEvent struct
      idx_first_daughter = -2;
      idx_last_daughter = -2;

      // has daughters and they were preselected in the fEvent struct
      if (fIndexer_fMC_to_fEvent.count(mcPart->GetDaughterFirst())) {
        idx_first_daughter = fIndexer_fMC_to_fEvent[mcPart->GetDaughterFirst()];
      }
      if (fIndexer_fMC_to_fEvent.count(mcPart->GetDaughterLast())) {
        idx_last_daughter = fIndexer_fMC_to_fEvent[mcPart->GetDaughterLast()];
      }
    }

    // (tree operations)
    fEvent.MC_FirstDau[fEvent_index] = idx_first_daughter;
    fEvent.MC_LastDau[fEvent_index] = idx_last_daughter;

    // (debug)
    // AliInfoF("%12i %12i %12i %12i", fMC_index, fEvent_index, idx_first_daughter, idx_last_daughter);
  }  // end of loop
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_VerifySourceAndGeneration(AliMCParticle* mcPart, Int_t& generation, Bool_t& label) {
  //
  // Determine if a certain MC particle comes from a generated sexaquark-nucleon interaction
  // (make sure to run this function for relevant particles first!)
  // PENDING: what to do with the G reaction channel
  //

  // (1) define conditions
  Bool_t isSignal = mcPart->MCStatusCode() == 6;
  Bool_t isPrimary = mcPart->GetMother() == -1;

  // assign generation and leave
  if (isPrimary) {
    label = isSignal;
    generation = 1;
    return;
  }  // closure

  // (2) if it's not primary, it must have a mother
  AliMCParticle* mcMother = (AliMCParticle*)fMC->GetTrack(mcPart->GetMother());

  Bool_t isMotherSignal = mcMother->MCStatusCode() == 6;
  Bool_t isMotherPrimary = mcMother->GetMother() == -1;

  if (isMotherPrimary) {
    label = isMotherSignal;
    generation = 2;
    return;
  }

  // default values
  label = kFALSE;
  generation = 3;
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCRec_FirstProcess() {
  //
  // (First) Loop over the reconstructed tracks in a single event
  // Store tracks into the (duplicated tracks) map
  // IMPORTANT: don't trust track->GetID()
  //

  // define variables & objects
  AliESDtrack* track;
  Float_t track_PID[3];
  Bool_t isRecPIDRelevant;

  Int_t link_label;
  AliMCParticle* link;
  Bool_t isGenPIDRelevant;

  // (debug)
  // AliInfo("RECONSTRUCTED TRACKS");
  // AliInfo("");
  // AliInfo(" Index   Link Charge nSigma(    pion     kaon   proton)       Px       Py       Pz");
  // AliInfoF("Number of tracks in this event: %i", fESD->GetNumberOfTracks());

  // loop over MC rec. tracks
  for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {

    track = static_cast<AliESDtrack*>(fESD->GetTrack(i));

    // get reconstructed PID
    track_PID[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    track_PID[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
    track_PID[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    isRecPIDRelevant = track_PID[0] < N_SIGMA || track_PID[1] < N_SIGMA || track_PID[2] < N_SIGMA;

    // get linked MC particle
    link_label = TMath::Abs(track->GetLabel());
    link = (AliMCParticle*)fMC->GetTrack(link_label);

    // get true PID
    isGenPIDRelevant = TMath::Abs(link->PdgCode()) == 211;            // is neg. / pos. pion
    isGenPIDRelevant = isGenPIDRelevant || link->PdgCode() == 321;    // is pos. kaon
    isGenPIDRelevant = isGenPIDRelevant || link->PdgCode() == -2212;  // is anti-proton

    if (isRecPIDRelevant || isGenPIDRelevant) {

      // (debug)
      // AliInfoF("%6i %6i %6i        %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f", i, link_label, track->Charge(), track_PID[0], track_PID[1],
      //  track_PID[2], track->Px(), track->Py(), track->Pz());

      // (duplicated tracks)
      // add up the reconstructed indices that belong-to/come-from/are-matched-to a same MC index
      fMap_DuplicatedTracks[link_label].push_back(i);

      // (non-relevant, but rec. as relevant)
      // to be used by "MCGen_SecondProcess()"
      if (!isGenPIDRelevant && isRecPIDRelevant) {
        fIndex_NonRelevantMC.insert(link_label);
      }
    }
  }  // end of loop over tracks
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCRec_SecondProcess() {
  //
  // (Second) Loop over the (duplicated tracks) map
  // If a MC particle has been reconstructed more than once
  // Tag the best particle as "non-duplicate", and the rest as "duplicates"
  //

  // define objects & variables
  Int_t idx_MCGen;
  AliMCParticle* mcPart;
  TVector3 p_MCGen;

  std::vector<Int_t> indices_MCRec;
  AliESDtrack* track;
  TVector3 p_MCRec;

  Int_t idx_best;
  Double_t min_angle = TMath::Pi();
  Double_t delta_angle;

  // loop over (duplicated tracks) map
  for (auto it = fMap_DuplicatedTracks.begin(); it != fMap_DuplicatedTracks.end(); it++) {

    idx_MCGen = it->first;
    indices_MCRec = it->second;

    // has this MC index been called more than once?
    if ((Int_t)indices_MCRec.size() > 1) {

      // load MC gen. particle
      mcPart = (AliMCParticle*)fMC->GetTrack(idx_MCGen);
      p_MCGen.SetXYZ(mcPart->Px(), mcPart->Py(), mcPart->Pz());

      // reset minimum angle
      min_angle = TMath::Pi();

      // loop (1), look for best candidate
      for (Int_t idx_MCRec : indices_MCRec) {

        // load MC rec. track
        track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_MCRec));
        p_MCRec.SetXYZ(track->Px(), track->Py(), track->Pz());

        delta_angle = p_MCGen.Angle(p_MCRec);

        if (delta_angle < min_angle) {
          idx_best = idx_MCRec;
          min_angle = delta_angle;
        }
      }

      // loop (2), push all tracks accordingly
      // now that we know which one is the best
      for (Int_t idx_MCRec : indices_MCRec) {

        // load MC rec. track
        track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_MCRec));

        // (tree operations)
        // before pushing back, it's important to correctly index this
        fIndexer_fTrack_to_fEvent[idx_MCRec] = fEvent.N_MCRec;
        // "idx_MCRec != idx_best" == kTRUE, when duplicate
        // "idx_MCRec != idx_best" == kFALSE, when original/best
        MCRec_PushBack(track, idx_MCRec != idx_best);
      }
    } else {
      // if the MC particle has not been called more than once,
      // then, load the unique track and push it back accordingly
      Int_t idx_MCRec = indices_MCRec[0];
      track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_MCRec));

      // (tree operations)
      // before pushing back, it's important to correctly index this
      fIndexer_fTrack_to_fEvent[idx_MCRec] = fEvent.N_MCRec;
      // as kFALSE, because is an only track
      MCRec_PushBack(track, kFALSE);
    }
  }  // end of loop over (duplicated tracks) map
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::V0_ProcessESD() {
  //
  // Loop over all the found V0s in the event
  //

  // define variables & objects
  AliESDv0* V0;

  AliESDtrack* neg_track;
  AliESDtrack* pos_track;

  // (debug)
  AliMCParticle* mc_neg;
  AliMCParticle* mc_pos;

  V0_tt common_format;

  Double_t V0_X, V0_Y, V0_Z;
  Double_t N_Px, N_Py, N_Pz;
  Double_t P_Px, P_Py, P_Pz;

  Float_t nsigma_pos_pion;
  Float_t nsigma_neg_pion;
  Float_t nsigma_antiproton;

  Double_t pos_energy;
  Double_t neg_energy_asK0;
  Double_t neg_energy_asAL;

  Double_t decay_length;
  Double_t dca_to_pv;

  // get primary vertex of this event
  Double_t PV[3];
  const AliESDVertex* prim_vertex = fESD->GetPrimaryVertex();
  prim_vertex->GetXYZ(PV);

  // (debug)
  AliInfo("V0s FROM ESD");
  AliInfoF("-- original number of V0s in this event = %i --", fESD->GetNumberOfV0s());
  AliInfo(" Index  N_Idx  P_Idx  N_PID  P_PID   V0_X   V0_Y   V0_Z        CPA");

  // loop over V0s
  for (Int_t i = 0; i < fESD->GetNumberOfV0s(); i++) {

    V0 = fESD->GetV0(i);

    // (cut) if any daughter (neg. or pos.) is not in the indexer
    // => the particle was gen. non-relevant & rec. non-relevant
    if (fIndexer_fTrack_to_fEvent.count(V0->GetPindex()) == 0 || fIndexer_fTrack_to_fEvent.count(V0->GetNindex()) == 0) {
      // count in std::map, returns the number of elements matching specific key
      continue;
    }

    // get coordinates of V0
    V0->GetXYZ(V0_X, V0_Y, V0_Z);

    // load rec. tracks
    neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(V0->GetNindex()));
    pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(V0->GetPindex()));

    // load red. tracks param
    /*
    AliExternalTrackParam pos_track_param(*pos_track);
    AliExternalTrackParam neg_track_param(*neg_track);
    */

    // (debug) load MC particles
    mc_neg = (AliMCParticle*)fMC->GetTrack(TMath::Abs(neg_track->GetLabel()));
    mc_pos = (AliMCParticle*)fMC->GetTrack(TMath::Abs(pos_track->GetLabel()));

    // get NSigmas for PID
    nsigma_pos_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion));
    nsigma_neg_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kPion));
    nsigma_antiproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kProton));

    // get momentum of negative and positive components at V0 vertex
    V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
    V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

    // calculate energy
    // for positive daughter, assuming it's a positive pion
    pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + kMassPion * kMassPion);
    // for negative daughter, assuming it's K0 -> pi+ pi-
    neg_energy_asK0 = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassPion * kMassPion);
    // for negative daughter, assuming it's anti-lambda -> pi+ anti-proton
    neg_energy_asAL = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassProton * kMassProton);

    // (debug) remove on-the-fly V0s
    if (V0->GetOnFlyStatus()) {
      continue;
    }

    // (debug)
    TString level = "-> L4";
    if (V0->GetV0CosineOfPointingAngle() > 0.998) {
      level = "-> L1";
    } else if (V0->GetV0CosineOfPointingAngle() > 0.99) {
      level = "-> L2";
    } else if (V0->GetV0CosineOfPointingAngle() > 0.9) {
      level = "-> L3";
    }

    // (debug)
    AliInfoF("%6i %6i %6i %6i %6i %6.2f %6.2f %6.2f %10.6f %s",  //
             i,                                                  //
             V0->GetNindex(),                                    //
             V0->GetPindex(),                                    //
             mc_neg->PdgCode(),                                  //
             mc_pos->PdgCode(),                                  //
             V0_X, V0_Y, V0_Z,                                   //
             V0->GetV0CosineOfPointingAngle(),                   //
             level.Data());

    // geometry
    decay_length = TMath::Sqrt(TMath::Power(V0_X - PV[0], 2) + TMath::Power(V0_Y - PV[1], 2) + TMath::Power(V0_Z - PV[2], 2));
    dca_to_pv = Calculate_LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]);

    // (tree operations)
    // fill common format and push_back
    common_format.Idx_Pos = V0->GetPindex();
    common_format.Idx_Neg = V0->GetNindex();
    common_format.Px = V0->Px();
    common_format.Py = V0->Py();
    common_format.Pz = V0->Pz();
    common_format.X = V0_X;
    common_format.Y = V0_Y;
    common_format.Z = V0_Z;
    common_format.Pos_Px = P_Px;
    common_format.Pos_Py = P_Py;
    common_format.Pos_Pz = P_Pz;
    common_format.Neg_Px = N_Px;
    common_format.Neg_Py = N_Py;
    common_format.Neg_Pz = N_Pz;
    common_format.E_asK0 = pos_energy + neg_energy_asK0;
    common_format.E_asAL = pos_energy + neg_energy_asAL;
    common_format.couldBeK0 = nsigma_pos_pion < N_SIGMA && nsigma_neg_pion < N_SIGMA;
    common_format.couldBeAL = nsigma_pos_pion < N_SIGMA && nsigma_antiproton < N_SIGMA;
    common_format.onFlyStatus = V0->GetOnFlyStatus();
    common_format.Chi2 = V0->GetChi2V0();
    common_format.DCA_Daughters = V0->GetDcaV0Daughters();
    common_format.IP_wrtPV = V0->GetD(PV[0], PV[1], PV[2]);
    common_format.CPA_wrtPV = V0->GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
    common_format.ArmAlpha = V0->AlphaV0();
    common_format.ArmPt = V0->PtArmV0();
    common_format.DecayLength = decay_length;
    common_format.DCA_wrtPV = dca_to_pv;
    common_format.isPrimary = dca_to_pv < 2.;  // Fabio's condition

    V0_PushBack(common_format);
  }  // end of loop over V0s
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::V0_ProcessCustom() {
  //
  // hola
  // [copied and modified from AliRoot/STEER/ESD/AliV0vertexer.cxx]
  //

  // define variables & objects
  AliESDtrack* esdTrack;
  AliESDtrack* ntrk;
  AliESDtrack* ptrk;

  // (debug)
  AliMCParticle* mc_neg;
  AliMCParticle* mc_pos;

  V0_tt common_format;

  Int_t nentr = fESD->GetNumberOfTracks();
  TArrayI neg(nentr);
  TArrayI pos(nentr);

  Double_t V0_X, V0_Y, V0_Z;
  Double_t N_Px, N_Py, N_Pz;
  Double_t P_Px, P_Py, P_Pz;

  Float_t nsigma_pos_pion;
  Float_t nsigma_neg_pion;
  Float_t nsigma_antiproton;

  Double_t pos_energy;
  Double_t neg_energy_asK0;
  Double_t neg_energy_asAL;
  Double_t mass_asK0;
  Double_t mass_asAL;

  Double_t radius;

  Double_t decay_length;
  Double_t dca_to_pv;

  // define cuts
  Double_t COV0F_DNmin = 0.1;            // (d=0.1) min imp parameter for the negative daughter (depends on PV!)
  Double_t COV0F_DPmin = 0.1;            // (d=0.1) min imp parameter for the positive daughter (depends on PV!)
  Double_t COV0F_Chi2max = 33.;          // (d=33.) max chi2
  Double_t COV0F_DCAmax = 1.5;           // (d=1.) max DCA between the daughter tracks
  Double_t COV0F_CPAmax = 0.98;          // (d=0.998, and it's a minimum instead) max cosine of V0's pointing angle (depends on PV!)
  Double_t COV0F_Rmin = 50.;             // (d=0.9) min radius of the V0 radius
  Double_t COV0F_Rmax = 200.;            // (d=100) max radius of the V0 radius
  Double_t COV0F_Etamax_V0 = 1.0;        // (d=1.) max eta of V0
  Double_t COV0F_DCAmin_V0_wrtPV = 5.;   // min DCA between V0 and PV
  Double_t COV0F_Etamax_NegDau = 1.5;    // max eta of negative daughter
  Double_t COV0F_Etamax_PosDau = 1.5;    // max eta of positive daughter
  Int_t COV0F_NClustersTPCmin_Dau = 50;  // minimun number of clusters in the TPC
  Double_t COV0F_Ptmin = 1.0;            // min Pt of V0

  // get primary vertex of this event
  Double_t PV[3];
  const AliESDVertex* prim_vertex = fESD->GetPrimaryVertex();
  prim_vertex->GetXYZ(PV);

  // (debug)
  AliInfo("V0s FROM CUSTOM OFFLINE V0 FINDER");
  AliInfo("");

  Double_t b = fESD->GetMagneticField();

  if (nentr < 2) {
    return;
  }

  // (addition)
  Bool_t UseImprovedFinding = kTRUE;

  // (addition)
  // [based on AliV0HypSel::AccountBField(b) and AliV0HypSel::SetBFieldCoef(v)]
  // account effect of B-field on pT resolution, ignoring the fact that the V0 mass resolution
  // is only partially determined by the prongs pT resolution
  const Float_t kNomField = 5.00668e+00;
  Float_t babs = TMath::Abs(b);
  Float_t BFieldCoef;
  if (babs > 1e-3) {
    BFieldCoef = kNomField / babs > 0. ? kNomField / babs : 1.0;
  }

  // (addition)
  const Int_t nHypSel = 6;
  TString v0_hyp_name[nHypSel] = {"gamma", "K0", "Lambda", "antiLambda", "HyperTriton", "antiHyperTriton"};
  Float_t v0_hyp_m0[nHypSel] = {0.5486e-3, 139.570e-3, 938.272e-3, 139.570e-3, 2.8092, 139.570e-3};
  Float_t v0_hyp_m1[nHypSel] = {0.5486e-3, 139.570e-3, 139.570e-3, 938.272e-3, 139.570e-3, 2.8092};
  Float_t v0_hyp_mass[nHypSel] = {1.099e-3, 497.7e-3, 1115.683e-3, 1115.683e-3, 2.992, 2.992};
  Float_t v0_hyp_sigma[nHypSel] = {0.001, 0.003, 0.001, 0.001, 0.0025, 0.0025};
  Float_t v0_hyp_cf0[nHypSel] = {20, 20, 20, 20, 14, 14};
  Float_t v0_hyp_cf1[nHypSel] = {0.6, 0.07, 0.07, 0.07, 0.07, 0.07};
  Float_t v0_hyp_nsig[nHypSel] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  Float_t v0_hyp_margin[nHypSel] = {0.0, 0.5, 0.5, 0.5, 0.5, 0.5};

  // number of charged tracks found
  Int_t nneg = 0;
  Int_t npos = 0;

  // first loop: find neg and pos tracks
  for (Int_t i = 0; i < nentr; i++) {

    esdTrack = fESD->GetTrack(i);

    // (cut) if particle is not in the indexer
    // => the particle was gen. non-relevant & rec. non-relevant
    if (fIndexer_fTrack_to_fEvent.count(i) == 0) {
      // std::map.count() returns the number of elements matching a specific key
      continue;
    }

    // (cut) on status
    ULong64_t status = esdTrack->GetStatus();
    if ((status & AliESDtrack::kTPCrefit) == 0) {
      continue;
    }

    // (cut) min. number of TPC clusters
    if ((Int_t)esdTrack->GetTPCncls() < COV0F_NClustersTPCmin_Dau) {
      continue;
    }

    // cut on transverse impact parameter
    Double_t d = 0.0;
    if (esdTrack->GetSign() < 0.) {
      d = esdTrack->GetD(PV[0], PV[1], b);
      if (TMath::Abs(d) < COV0F_DNmin) {
        continue;
      }
      if (TMath::Abs(d) > COV0F_Rmax) {
        continue;
      }
      neg[nneg++] = i;
    } else {
      d = esdTrack->GetD(PV[0], PV[1], b);
      if (TMath::Abs(d) < COV0F_DPmin) {
        continue;
      }
      if (TMath::Abs(d) > COV0F_Rmax) {
        continue;
      }
      pos[npos++] = i;
    }
  }  // end of first loop

  // (debug)
  /*
  AliInfo(" Index  N_Idx  P_Idx  N_PID  P_PID   V0_X   V0_Y   V0_Z        CPA");
  */

  // second loop: test all neg-pos pairs
  for (Int_t i = 0; i < nneg; i++) {

    Int_t nidx = neg[i];
    ntrk = fESD->GetTrack(nidx);

    for (Int_t k = 0; k < npos; k++) {

      Int_t pidx = pos[k];
      ptrk = fESD->GetTrack(pidx);

      // (debug) load MC particles
      /*
      mc_neg = (AliMCParticle*)fMC->GetTrack(TMath::Abs(ntrk->GetLabel()));
      mc_pos = (AliMCParticle*)fMC->GetTrack(TMath::Abs(ptrk->GetLabel()));
      */

      // (cut) -1. < eta(pi+) < 1.
      if (TMath::Abs(ptrk->Eta()) > COV0F_Etamax_PosDau) {
        continue;
      }

      // (cut) -1.5 < eta(pi-,anti-p) < 1.5
      if (TMath::Abs(ntrk->Eta()) > COV0F_Etamax_NegDau) {
        continue;
      }

      AliExternalTrackParam nt(*ntrk), pt(*ptrk), *ntp = &nt, *ptp = &pt;

      Double_t xn, xp;
      if (UseImprovedFinding && Preoptimize(ntp, ptp, &xn, &xp, b)) {
        // Move tracks to a better position if that helps
        nt.PropagateTo(xn, b);
        pt.PropagateTo(xp, b);
      }

      // (cut) DCA between daughters
      Double_t dca = ntp->GetDCA(ptp, b, xn, xp);
      if (dca > COV0F_DCAmax) {
        continue;
      }

      // idk what's this...
      if ((xn + xp) > 2 * COV0F_Rmax) {
        continue;
      }
      if ((xn + xp) < 2 * COV0F_Rmin) {
        continue;
      }

      nt.PropagateTo(xn, b);
      pt.PropagateTo(xp, b);

      // important object!
      AliESDv0 vertex(nt, nidx, pt, pidx);

      if (UseImprovedFinding) {
      vertex.Refit();  // imp pos + cov mat
      }

      // (cut) chi2 ~ quality of fit
      if (vertex.GetChi2V0() > COV0F_Chi2max) {
        continue;
      }

      vertex.SetDcaV0Daughters(dca);

      // (cut) CPA
      Float_t cpa = vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
      vertex.SetV0CosineOfPointingAngle(cpa);
      if (cpa > COV0F_CPAmax) {
        continue;
      }

      // (cut) radius(V0)
      Double_t x = vertex.Xv();
      Double_t y = vertex.Yv();
      Double_t v0_radius = TMath::Sqrt(x * x + y * y);
      if (v0_radius < COV0F_Rmin) {
        continue;
      }
      if (v0_radius > COV0F_Rmax) {
        continue;
      }

      // load coordinates and momentum of each component of the V0
      vertex.GetXYZ(V0_X, V0_Y, V0_Z);

      // (cut) dca(V0,PV)
      dca_to_pv = Calculate_LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]);
      if (dca_to_pv < COV0F_DCAmin_V0_wrtPV) {
          continue;
        }

      // (addition)
      if (nHypSel) {  // do we select particular hypotheses?
        Bool_t reject = kTRUE;
        Float_t v0_pt = vertex.Pt();
        for (int hh = 0; hh < nHypSel; hh++) {
          Double_t v0_eff_mass = vertex.GetEffMassExplicit(v0_hyp_m0[hh], v0_hyp_m1[hh]);
          Float_t v0_mass_margin = v0_hyp_nsig[hh] * BFieldCoef * v0_hyp_sigma[hh] * (v0_hyp_cf0[hh] + v0_pt * v0_hyp_cf1[hh])  //
                                   + v0_hyp_margin[hh];
          if (TMath::Abs(v0_eff_mass - v0_hyp_mass[hh]) < v0_mass_margin) {
            reject = kFALSE;
            break;
          }
        }
        if (reject) {
          continue;
        }
      }

      // (addition) (cut) remove v0s "not contributing to physics"
      if (TMath::Abs(vertex.Eta()) > COV0F_Etamax) {
        continue;
      }

      // vertex.ChangeMassHypothesis(kK0Short);   // borquez: this was here...

      // store V0
      // (removed) we don't need to store in the ESD, since we store it in the tree
      // event->AddV0(&vertex);

      // get momentum of negative and positive components at V0 vertex
      vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
      vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

      // (debug)
      /*
      TString level = "-> L4";
      if (cpa > 0.998) {
        level = "-> L1";
      } else if (cpa > 0.99) {
        level = "-> L2";
      } else if (cpa > 0.9) {
        level = "-> L3";
      }
      */

      // (debug)
      /*
      AliInfoF("%6i %6i %6i %6i %6i %6.2f %6.2f %6.2f %10.6f %s",  //
               nvtx,                                               //
               vertex.GetNindex(),                                 //
               vertex.GetPindex(),                                 //
               mc_neg->PdgCode(),                                  //
               mc_pos->PdgCode(),                                  //
               V0_X, V0_Y, V0_Z,                                   //
               cpa,                                                //
               level.Data());
      */

      // (cut) on Pt
      if (TMath::Sqrt(vertex.Px() * vertex.Px() + vertex.Py() * vertex.Py()) < COV0F_Ptmin) {
        continue;
      }

      // get NSigmas for PID
      nsigma_pos_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ptrk, AliPID::kPion));
      nsigma_neg_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kPion));
      nsigma_antiproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(ntrk, AliPID::kProton));

      // calculate energy
      // for positive daughter, assume it's a positive pion
      pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + kMassPion * kMassPion);
      // for negative daughter, assuming it's K0 -> pi+ pi-
      neg_energy_asK0 = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassPion * kMassPion);
      // for negative daughter, assuming it's anti-lambda -> pi+ anti-proton
      neg_energy_asAL = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassProton * kMassProton);

      // calc. mass
      mass_asK0 = TMath::Sqrt((pos_energy + neg_energy_asK0) * (pos_energy + neg_energy_asK0) -  //
                              vertex.Px() * vertex.Px() - vertex.Py() * vertex.Py() - vertex.Pz() * vertex.Pz());
      mass_asAL = TMath::Sqrt((pos_energy + neg_energy_asAL) * (pos_energy + neg_energy_asAL) -  //
                              vertex.Px() * vertex.Px() - vertex.Py() * vertex.Py() - vertex.Pz() * vertex.Pz());

      // (cut) on mass
      if (nsigma_pos_pion < N_SIGMA && nsigma_neg_pion < N_SIGMA && (mass_asK0 < 0.485 || mass_asK0 > 0.515)) {
        continue;
      }
      if (nsigma_pos_pion < N_SIGMA && nsigma_antiproton < N_SIGMA && (mass_asAL < 1.1 || mass_asAL > 1.13)) {
        continue;
      }

      // geometry
      decay_length = TMath::Sqrt(TMath::Power(V0_X - PV[0], 2) + TMath::Power(V0_Y - PV[1], 2) + TMath::Power(V0_Z - PV[2], 2));

      // (tree operations)
      // fill common format and push_back
      common_format.Idx_Pos = vertex.GetPindex();
      common_format.Idx_Neg = vertex.GetNindex();
      common_format.Px = vertex.Px();
      common_format.Py = vertex.Py();
      common_format.Pz = vertex.Pz();
      common_format.X = V0_X;
      common_format.Y = V0_Y;
      common_format.Z = V0_Z;
      common_format.Pos_Px = P_Px;
      common_format.Pos_Py = P_Py;
      common_format.Pos_Pz = P_Pz;
      common_format.Neg_Px = N_Px;
      common_format.Neg_Py = N_Py;
      common_format.Neg_Pz = N_Pz;
      common_format.E_asK0 = pos_energy + neg_energy_asK0;
      common_format.E_asAL = pos_energy + neg_energy_asAL;
      common_format.couldBeK0 = nsigma_pos_pion < N_SIGMA && nsigma_neg_pion < N_SIGMA;
      common_format.couldBeAL = nsigma_pos_pion < N_SIGMA && nsigma_antiproton < N_SIGMA;
      common_format.onFlyStatus = kFALSE;  // custom V0 finder cannot be on-the-fly
      common_format.Chi2 = vertex.GetChi2V0();
      common_format.DCA_Daughters = vertex.GetDcaV0Daughters();
      common_format.IP_wrtPV = vertex.GetD(PV[0], PV[1], PV[2]);
      common_format.CPA_wrtPV = vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
      common_format.ArmAlpha = vertex.AlphaV0();
      common_format.ArmPt = vertex.PtArmV0();
      common_format.DecayLength = decay_length;
      common_format.DCA_wrtPV = dca_to_pv;
      common_format.isPrimary = dca_to_pv < 2.;  // Fabio's condition

      V0_PushBack(common_format);
    }  // end of loop over pos. tracks
  }    // end of loop over neg. tracks
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::V0_ProcessTrue() {
  //
  // Alternative to ProcessV0s()
  // we will store all lambdas and kaons with no ambiguity
  // - no dependence on PID, nor ambiguity on measurements
  // - only dependence: if they were reconstructed or not
  //

  // define variables & objects
  AliMCParticle* mc_mother;
  AliMCParticle* mc_firstdau;
  AliMCParticle* mc_lastdau;
  AliMCParticle* mc_pos;
  AliMCParticle* mc_neg;

  Int_t idx_pos_mc;
  Int_t idx_neg_mc;
  Int_t mother_idx;

  Bool_t antiLambdaNonPi0Channel;
  Bool_t neutralKaonNonPi0Channel;

  Double_t pos_energy;
  Double_t neg_energy_asK0;
  Double_t neg_energy_asAL;

  V0_tt common_format;

  // get primary vertex of this event
  Double_t PV[3];
  const AliVVertex* prim_vertex = fMC->GetPrimaryVertex();
  PV[0] = prim_vertex->GetX();
  PV[1] = prim_vertex->GetY();
  PV[2] = prim_vertex->GetZ();

  // (debug)
  AliInfo("V0s FROM TRUE");
  AliInfo("");
  AliInfo("   idx_pos    idx_neg  idx_pos_true  idx_neg_true    pid_pos    pid_neg idx_mc_mother pid_mother");

  // loop over MC gen. particles
  for (Int_t i = 0; i < fMC->GetNumberOfTracks(); i++) {

    mc_mother = (AliMCParticle*)fMC->GetTrack(i);

    // (cut) on PID, must be K0s or anti-lambda
    if (mc_mother->PdgCode() != 310 && mc_mother->PdgCode() != -3122) {
      continue;
    }

    // (cut) on number of daughters
    if (mc_mother->GetNDaughters() != 2) {
      continue;
    }

    // get the daughters
    mc_firstdau = (AliMCParticle*)fMC->GetTrack(mc_mother->GetDaughterFirst());
    mc_lastdau = (AliMCParticle*)fMC->GetTrack(mc_mother->GetDaughterLast());

    antiLambdaNonPi0Channel = (mc_firstdau->PdgCode() == -2212 && mc_lastdau->PdgCode() == 211) ||
                              (mc_firstdau->PdgCode() == 211 && mc_lastdau->PdgCode() == -2212);

    neutralKaonNonPi0Channel = (mc_firstdau->PdgCode() == -211 && mc_lastdau->PdgCode() == 211) ||
                               (mc_firstdau->PdgCode() == 211 && mc_lastdau->PdgCode() == -211);

    // check for non-pi0 decay channel
    if (!antiLambdaNonPi0Channel && !neutralKaonNonPi0Channel) {
      continue;
    }

    // check if both daughters were reconstructed
    if (fMap_DuplicatedTracks.count(mc_mother->GetDaughterFirst()) == 0 || fMap_DuplicatedTracks.count(mc_mother->GetDaughterLast()) == 0) {
      // reminder: count in std::map, returns the number of elements matching specific key
      continue;
    }

    // identify daughters as positive/negative
    if (mc_firstdau->PdgCode() > 0 && mc_lastdau->PdgCode() < 0) {
      idx_pos_mc = mc_mother->GetDaughterFirst();
      idx_neg_mc = mc_mother->GetDaughterLast();
    } else {
      idx_pos_mc = mc_mother->GetDaughterLast();
      idx_neg_mc = mc_mother->GetDaughterFirst();
    }

    // do combinations
    for (Int_t idx_pos_track : fMap_DuplicatedTracks[idx_pos_mc]) {
      for (Int_t idx_neg_track : fMap_DuplicatedTracks[idx_neg_mc]) {

        // load MC particle depending if positive/negative
        mc_pos = (AliMCParticle*)fMC->GetTrack(idx_pos_mc);
        mc_neg = (AliMCParticle*)fMC->GetTrack(idx_neg_mc);

        // (debug)
        AliInfoF("%10i %10i %13i %13i %10i %10i %13i %10i", idx_pos_track, idx_neg_track, idx_pos_mc, idx_neg_mc, mc_pos->PdgCode(),
                 mc_neg->PdgCode(), i, mc_mother->PdgCode());

        // calculate energy (assuming positive pion)
        pos_energy =
            TMath::Sqrt(mc_pos->Px() * mc_pos->Px() + mc_pos->Py() * mc_pos->Py() + mc_pos->Pz() * mc_pos->Pz() + kMassPion * kMassPion);
        // assuming K0 -> pi+ pi-
        neg_energy_asK0 =
            TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() + kMassPion * kMassPion);
        // assuming anti-lambda -> pi+ anti-proton
        neg_energy_asAL = TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() +
                                      kMassProton * kMassProton);

        // (tree operations)
        // fill common format and push_back
        common_format.Idx_Pos = idx_pos_track;
        common_format.Idx_Neg = idx_neg_track;
        common_format.Px = mc_mother->Px();
        common_format.Py = mc_mother->Py();
        common_format.Pz = mc_mother->Pz();
        common_format.X = mc_pos->Xv();  // use the origin of the daughters
        common_format.Y = mc_pos->Yv();
        common_format.Z = mc_pos->Zv();
        common_format.Pos_Px = mc_pos->Px();
        common_format.Pos_Py = mc_pos->Py();
        common_format.Pos_Pz = mc_pos->Pz();
        common_format.Neg_Px = mc_neg->Px();
        common_format.Neg_Py = mc_neg->Py();
        common_format.Neg_Pz = mc_neg->Pz();
        common_format.E_asK0 = pos_energy + neg_energy_asK0;
        common_format.E_asAL = pos_energy + neg_energy_asAL;
        common_format.couldBeK0 = mc_mother->PdgCode() == 310;
        common_format.couldBeAL = mc_mother->PdgCode() == -3122;
        common_format.onFlyStatus = kFALSE;  // true V0s don't come from any V0 finder
        common_format.Chi2 = 1.0;            // no KF was used for this
        common_format.DCA_Daughters = 0.;    // PENDING! too complicated to calculate using true information...
        common_format.IP_wrtPV = 0.;         // PENDING! still not sure about the definition of this one...
        common_format.CPA_wrtPV = Calculate_CPA(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),      //
                                                mc_pos->Xv(), mc_pos->Yv(), mc_pos->Zv(),               //
                                                PV[0], PV[1], PV[2]);                                   //
        common_format.ArmAlpha = Calculate_ArmAlpha(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),  //
                                                    mc_pos->Px(), mc_pos->Py(), mc_pos->Pz(),           //
                                                    mc_neg->Px(), mc_neg->Py(), mc_neg->Pz());          //
        common_format.ArmPt = Calculate_ArmPt(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),        //
                                              mc_neg->Px(), mc_neg->Py(), mc_neg->Pz());                //
        common_format.DecayLength = TMath::Sqrt(TMath::Power(mc_pos->Xv() - PV[0], 2) + TMath::Power(mc_pos->Yv() - PV[1], 2) +
                                                TMath::Power(mc_pos->Zv() - PV[2], 2));
        common_format.DCA_wrtPV = Calculate_LinePointDCA(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),  //
                                                         mc_pos->Xv(), mc_pos->Yv(), mc_pos->Zv(),           //
                                                         PV[0], PV[1], PV[2]);
        common_format.isPrimary = mc_mother->GetMother() == -1 && mc_mother->MCStatusCode() != 6;

        V0_PushBack(common_format);
      }  // end of loop over negative tracks
    }    // end of loop over positive tracks

  }  // end of loop over MC gen. particles
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::Terminate(Option_t*) {
  //
  // terminate
  // called at the END of the analysis (when all events are processed)
  //
  return;
}

/*** MATHEMATICAL FUNCTIONS ***/

//________________________________________________________________________
Double_t AliAnalysisTaskSexaquark::Calculate_CPA(Double_t Px, Double_t Py, Double_t Pz,  //
                                                 Double_t X, Double_t Y, Double_t Z,     //
                                                 Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
  //
  // calculates the cosine of the pointing angle
  // of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
  // [based on AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()]
  //

  // vector between the reference point and the particle's vertex
  Double_t deltaPos[3];
  deltaPos[0] = X - refPointX;
  deltaPos[1] = Y - refPointY;
  deltaPos[2] = Z - refPointZ;

  Double_t mom2 = Px * Px + Py * Py + Pz * Pz;
  Double_t deltaPos2 = deltaPos[0] * deltaPos[0] + deltaPos[1] * deltaPos[1] + deltaPos[2] * deltaPos[2];

  // (protection)
  if (deltaPos2 == 0. || mom2 == 0.) {
    return -2.;
  }

  return (deltaPos[0] * Px + deltaPos[1] * Py + deltaPos[2] * Pz) / TMath::Sqrt(mom2 * deltaPos2);
}

//________________________________________________________________________
Double_t AliAnalysisTaskSexaquark::Calculate_ArmAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
                                                      Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz,  //
                                                      Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
  //
  // This gives the Armenteros-Podolanski alpha
  // [based on AliRoot/STEER/ESD/AliESDv0::AlphaV0()]
  //
  TVector3 momTot(V0_Px, V0_Py, V0_Pz);
  TVector3 momPos(Pos_Px, Pos_Py, Pos_Pz);
  TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

  Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
  Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

  // (protection)
  if (lQlPos + lQlNeg == 0.) {
    return 2;
  }  // closure
  return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

//________________________________________________________________________
Double_t AliAnalysisTaskSexaquark::Calculate_ArmPt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                   Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
  //
  // This gives the Armenteros-Podolanski ptarm
  // [based on AliRoot/STEER/ESD/AliESDv0::PtArmV0()]
  //
  TVector3 momTot(V0_Px, V0_Py, V0_Pz);
  TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

  return momNeg.Perp(momTot);
}

//________________________________________________________________________
Double_t AliAnalysisTaskSexaquark::Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                Double_t PV_X, Double_t PV_Y, Double_t PV_Z) {
  //
  // Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
  // This function stores the point of closest approach, and returns its distance to the PV.
  //
  TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
  TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
  TVector3 PrimVertex(PV_X, PV_Y, PV_Z);

  TVector3 CrossProduct = (PrimVertex - V0Vertex).Cross(V0Momentum);

  return CrossProduct.Mag() / V0Momentum.Mag();
}

/*** CUSTOM OFFLINE V0 FINDER ***/

//________________________________________________
Bool_t AliAnalysisTaskSexaquark::Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn,
                                             Double_t* lPreprocessxp, const Double_t b) {
  //
  // This function pre-optimizes a two-track pair in the XY plane
  // and provides two X values for the tracks if successful
  // [copied from AliRoot/STEER/ESD/AliV0vertexer.cxx]
  //
  Double_t lMinimumX = -3.0;
  Double_t lMaximumX = 300.0;

  Double_t nhelix[6], phelix[6];
  nt->GetHelixParameters(nhelix, b);
  pt->GetHelixParameters(phelix, b);
  Double_t lNegCenterR[2], lPosCenterR[2];

  // Negative track parameters in XY
  GetHelixCenter(nt, lNegCenterR, b);
  Double_t xNegCenter = lNegCenterR[0];
  Double_t yNegCenter = lNegCenterR[1];
  Double_t NegRadius = TMath::Abs(1. / nhelix[4]);

  // Positive track parameters in XY
  GetHelixCenter(pt, lPosCenterR, b);
  Double_t xPosCenter = lPosCenterR[0];
  Double_t yPosCenter = lPosCenterR[1];
  Double_t PosRadius = TMath::Abs(1. / phelix[4]);

  // Define convenient coordinate system
  // Logical zero: position of negative center
  Double_t ux = xPosCenter - xNegCenter;
  Double_t uy = yPosCenter - yNegCenter;

  // Check center-to-center distance
  Double_t lDist = TMath::Sqrt(TMath::Power(xNegCenter - xPosCenter, 2) + TMath::Power(yNegCenter - yPosCenter, 2));
  // Normalize ux, uz to unit vector
  ux /= lDist;
  uy /= lDist;

  // Calculate perpendicular vector (normalized)
  Double_t vx = -uy;
  Double_t vy = +ux;

  Double_t lPreprocessDCAxy = 1e+3;  // define outside scope
  *lPreprocessxp = pt->GetX();       // start at current location
  *lPreprocessxn = nt->GetX();       // start at current location

  //============================================================
  // Pre-optimization in the XY plane: cases considered here
  //============================================================
  //
  //  Case 1: Circles do not touch, centers far away
  //          (D > R1 + R2)
  //
  //  Case 2: Circles touch, centers at reasonable distance wrt D
  //          (D < R1 + R2) && (D > |R1-R2|)
  //
  //  Case 3: Circles do not touch, one inside the other
  //          (D < |R1-R2|)
  //
  //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
  //  to be a problem with unlike-sign charged tracks): brute
  //  force minimization takes place in any case
  //
  //============================================================

  //______________________
  // CASE 1
  if ((lDist > NegRadius + PosRadius)) {
    //================================================================
    // Case 1: distance bigger than sum of radii ("gamma-like")
    //        re-position tracks along the center-to-center axis
    // Re-position negative track
    Double_t xNegOptPosition = xNegCenter + NegRadius * ux;
    Double_t yNegOptPosition = yNegCenter + NegRadius * uy;
    Double_t csNeg = TMath::Cos(nt->GetAlpha());
    Double_t snNeg = TMath::Sin(nt->GetAlpha());
    Double_t xThisNeg = xNegOptPosition * csNeg + yNegOptPosition * snNeg;

    // Re-position positive track
    Double_t xPosOptPosition = xPosCenter - PosRadius * ux;
    Double_t yPosOptPosition = yPosCenter - PosRadius * uy;
    Double_t csPos = TMath::Cos(pt->GetAlpha());
    Double_t snPos = TMath::Sin(pt->GetAlpha());
    Double_t xThisPos = xPosOptPosition * csPos + yPosOptPosition * snPos;

    if (xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX) {
      Double_t lCase1NegR[3] = {0.}, lCase1PosR[3] = {0.};
      if (nt->GetXYZAt(xThisNeg, b, lCase1NegR) && pt->GetXYZAt(xThisPos, b, lCase1PosR)) {
        lPreprocessDCAxy = TMath::Sqrt(TMath::Power(lCase1NegR[0] - lCase1PosR[0], 2) + TMath::Power(lCase1NegR[1] - lCase1PosR[1], 2) +
                                       TMath::Power(lCase1NegR[2] - lCase1PosR[2], 2));
        // Pass coordinates
        if (lPreprocessDCAxy < 999) {
          *lPreprocessxp = xThisPos;
          *lPreprocessxn = xThisNeg;
        }
      }
    }
    //================================================================
  }

  //______________________
  // CASE 2
  if ((lDist > TMath::Abs(NegRadius - PosRadius)) && (lDist < NegRadius + PosRadius)) {
    //================================================================
    // Case 2: distance smaller than sum of radii (cowboy/sailor configs)

    // Calculate coordinate for radical line
    Double_t lRadical = (lDist * lDist - PosRadius * PosRadius + NegRadius * NegRadius) / (2 * lDist);

    // Calculate absolute displacement from center-to-center axis
    Double_t lDisplace = (0.5 / lDist) * TMath::Sqrt((-lDist + PosRadius - NegRadius) * (-lDist - PosRadius + NegRadius) *
                                                     (-lDist + PosRadius + NegRadius) * (lDist + PosRadius + NegRadius));

    // 3D distances in the two cases studied (prefer smallest)
    Double_t lCase2aDCA = 1e+3;
    Double_t lCase2bDCA = 1e+3;

    // 2 cases: positive and negative displacement
    Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
    Double_t csNeg, snNeg, csPos, snPos;
    Double_t xThisNeg[2], xThisPos[2];

    csNeg = TMath::Cos(nt->GetAlpha());
    snNeg = TMath::Sin(nt->GetAlpha());
    csPos = TMath::Cos(pt->GetAlpha());
    snPos = TMath::Sin(pt->GetAlpha());

    // Case 2a: Positive displacement along v vector
    // Re-position negative track
    xNegOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
    yNegOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
    xThisNeg[0] = xNegOptPosition[0] * csNeg + yNegOptPosition[0] * snNeg;
    // Re-position positive track
    xPosOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
    yPosOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
    xThisPos[0] = xPosOptPosition[0] * csPos + yPosOptPosition[0] * snPos;

    // Case 2b: Negative displacement along v vector
    // Re-position negative track
    xNegOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
    yNegOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
    xThisNeg[1] = xNegOptPosition[1] * csNeg + yNegOptPosition[1] * snNeg;
    // Re-position positive track
    xPosOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
    yPosOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
    xThisPos[1] = xPosOptPosition[1] * csPos + yPosOptPosition[1] * snPos;

    // Test the two cases, please

    // Case 2a
    if (xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX) {
      Double_t lCase2aNegR[3] = {0.}, lCase2aPosR[3] = {0.};
      if (nt->GetXYZAt(xThisNeg[0], b, lCase2aNegR) && pt->GetXYZAt(xThisPos[0], b, lCase2aPosR)) {
        lCase2aDCA = TMath::Sqrt(TMath::Power(lCase2aNegR[0] - lCase2aPosR[0], 2) + TMath::Power(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                                 TMath::Power(lCase2aNegR[2] - lCase2aPosR[2], 2));
      }
    }

    // Case 2b
    if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX) {
      Double_t lCase2bNegR[3] = {0.}, lCase2bPosR[3] = {0.};
      if (nt->GetXYZAt(xThisNeg[1], b, lCase2bNegR) && pt->GetXYZAt(xThisPos[1], b, lCase2bPosR)) {
        lCase2bDCA = TMath::Sqrt(TMath::Power(lCase2bNegR[0] - lCase2bPosR[0], 2) + TMath::Power(lCase2bNegR[1] - lCase2bPosR[1], 2) +
                                 TMath::Power(lCase2bNegR[2] - lCase2bPosR[2], 2));
      }
    }
    // Minor detail: all things being equal, prefer closest X
    Double_t lCase2aSumX = xThisPos[0] + xThisNeg[0];
    Double_t lCase2bSumX = xThisPos[1] + xThisNeg[1];

    Double_t lDCAxySmallestR = lCase2aDCA;
    Double_t lxpSmallestR = xThisPos[0];
    Double_t lxnSmallestR = xThisNeg[0];

    Double_t lDCAxyLargestR = lCase2bDCA;
    Double_t lxpLargestR = xThisPos[1];
    Double_t lxnLargestR = xThisNeg[1];

    if (lCase2bSumX + 1e-6 < lCase2aSumX) {
      lDCAxySmallestR = lCase2bDCA;
      lxpSmallestR = xThisPos[1];
      lxnSmallestR = xThisNeg[1];
      lDCAxyLargestR = lCase2aDCA;
      lxpLargestR = xThisPos[0];
      lxnLargestR = xThisNeg[0];
    }

    // Pass conclusion to lPreprocess variables, please
    lPreprocessDCAxy = lDCAxySmallestR;
    *lPreprocessxp = lxpSmallestR;
    *lPreprocessxn = lxnSmallestR;
    if (lDCAxyLargestR + 1e-6 < lDCAxySmallestR) {  // beware epsilon: numerical calculations are unstable here
      lPreprocessDCAxy = lDCAxyLargestR;
      *lPreprocessxp = lxpLargestR;
      *lPreprocessxn = lxnLargestR;
    }
    // Protection against something too crazy, please
    if (lPreprocessDCAxy > 999) {
      *lPreprocessxp = pt->GetX();  // start at current location
      *lPreprocessxn = nt->GetX();  // start at current location
    }
  }
  // End of preprocessing stage!
  // at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
  Bool_t lWorked = kFALSE;
  if (lPreprocessDCAxy < 999) {
    lWorked = kTRUE;
  }
  return lWorked;
}

//________________________________________________
void AliAnalysisTaskSexaquark::GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b) {
  //
  // Copied from AliV0ReaderV1::GetHelixCenter
  // Get Center of the helix track parametrization
  // (based on AliRoot/STEER/ESD/AliV0vertexer.cxx)
  //

  Int_t charge = track->Charge();

  Double_t helix[6];
  track->GetHelixParameters(helix, b);

  Double_t xpos = helix[5];
  Double_t ypos = helix[0];
  Double_t radius = TMath::Abs(1. / helix[4]);
  Double_t phi = helix[2];
  if (phi < 0) {
    phi = phi + 2 * TMath::Pi();
  }

  phi -= TMath::Pi() / 2.;
  Double_t xpoint = radius * TMath::Cos(phi);
  Double_t ypoint = radius * TMath::Sin(phi);

  if (b < 0 && charge > 0) {
    xpoint = -xpoint;
    ypoint = -ypoint;
  }

  if (b > 0 && charge < 0) {
    xpoint = -xpoint;
    ypoint = -ypoint;
  }

  center[0] = xpos + xpoint;
  center[1] = ypos + ypoint;
  return;
}

/*** TREE OPERATIONS ***/

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::SetBranches() {
  //
  // Define fTree's branches, and load them into fEvent (which is an Event_tt object)
  //
  /* MC particles */
  fTree->Branch("N_MCGen", &fEvent.N_MCGen);
  fTree->Branch("MC_Px", &fEvent.MC_Px);
  fTree->Branch("MC_Py", &fEvent.MC_Py);
  fTree->Branch("MC_Pz", &fEvent.MC_Pz);
  fTree->Branch("MC_X", &fEvent.MC_X);
  fTree->Branch("MC_Y", &fEvent.MC_Y);
  fTree->Branch("MC_Z", &fEvent.MC_Z);
  fTree->Branch("MC_PID", &fEvent.MC_PID);
  fTree->Branch("MC_Mother", &fEvent.MC_Mother);
  fTree->Branch("MC_FirstDau", &fEvent.MC_FirstDau);
  fTree->Branch("MC_LastDau", &fEvent.MC_LastDau);
  fTree->Branch("MC_Gen", &fEvent.MC_Gen);
  fTree->Branch("MC_Status", &fEvent.MC_Status);
  fTree->Branch("MC_isSignal", &fEvent.MC_isSignal);
  /* MC Rec. Tracks */
  fTree->Branch("N_MCRec", &fEvent.N_MCRec);
  fTree->Branch("Idx_True", &fEvent.Idx_True);
  fTree->Branch("Rec_Px", &fEvent.Rec_Px);
  fTree->Branch("Rec_Py", &fEvent.Rec_Py);
  fTree->Branch("Rec_Pz", &fEvent.Rec_Pz);
  fTree->Branch("Rec_Charge", &fEvent.Rec_Charge);
  fTree->Branch("Rec_NSigmaPion", &fEvent.Rec_NSigmaPion);
  fTree->Branch("Rec_NSigmaProton", &fEvent.Rec_NSigmaProton);
  fTree->Branch("Rec_NClustersTPC", &fEvent.Rec_NClustersTPC);
  fTree->Branch("Rec_isDuplicate", &fEvent.Rec_isDuplicate);
  fTree->Branch("Rec_isSignal", &fEvent.Rec_isSignal);
  /* V0s */
  fTree->Branch("N_V0s", &fEvent.N_V0s);
  fTree->Branch("Idx_Pos", &fEvent.Idx_Pos);
  fTree->Branch("Idx_Neg", &fEvent.Idx_Neg);
  fTree->Branch("V0_Px", &fEvent.V0_Px);
  fTree->Branch("V0_Py", &fEvent.V0_Py);
  fTree->Branch("V0_Pz", &fEvent.V0_Pz);
  fTree->Branch("V0_X", &fEvent.V0_X);
  fTree->Branch("V0_Y", &fEvent.V0_Y);
  fTree->Branch("V0_Z", &fEvent.V0_Z);
  fTree->Branch("Pos_Px", &fEvent.Pos_Px);
  fTree->Branch("Pos_Py", &fEvent.Pos_Py);
  fTree->Branch("Pos_Pz", &fEvent.Pos_Pz);
  fTree->Branch("Neg_Px", &fEvent.Neg_Px);
  fTree->Branch("Neg_Py", &fEvent.Neg_Py);
  fTree->Branch("Neg_Pz", &fEvent.Neg_Pz);
  fTree->Branch("V0_isSignal", &fEvent.V0_isSignal);
  fTree->Branch("V0_E_asK0", &fEvent.V0_E_asK0);
  fTree->Branch("V0_E_asAL", &fEvent.V0_E_asAL);
  fTree->Branch("V0_couldBeK0", &fEvent.V0_couldBeK0);
  fTree->Branch("V0_couldBeAL", &fEvent.V0_couldBeAL);
  fTree->Branch("V0_onFlyStatus", &fEvent.V0_onFlyStatus);
  fTree->Branch("V0_Chi2", &fEvent.V0_Chi2);
  fTree->Branch("V0_DCA_Daughters", &fEvent.V0_DCA_Daughters);
  fTree->Branch("V0_IP_wrtPV", &fEvent.V0_IP_wrtPV);
  fTree->Branch("V0_CPA_wrtPV", &fEvent.V0_CPA_wrtPV);
  fTree->Branch("V0_ArmAlpha", &fEvent.V0_ArmAlpha);
  fTree->Branch("V0_ArmPt", &fEvent.V0_ArmPt);
  fTree->Branch("V0_DecayLength", &fEvent.V0_DecayLength);
  fTree->Branch("V0_DCA_wrtPV", &fEvent.V0_DCA_wrtPV);
  fTree->Branch("V0_isPrimary", &fEvent.V0_isPrimary);
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_PushBack(AliMCParticle* mcPart, Int_t generation, Bool_t isSignal) {
  //
  // Add one MC generated particle into the Event_tt object
  //
  // (1) prepare
  // (a) solve indexing to mother particle
  // -2 : mother is not relevant <-> mother outside the vector
  // -1 : particle has no mother <-> particle is primary
  //  0 : particle is primary && particle is signal, i.e., comes from a sexaquark
  //  X : particle has a relevant mother, located at index X
  Int_t idx_mother = -1;
  if (mcPart->GetMother() != -1) {
    // if mother is relevant <-> mother has been pushed before within the vector
    // ==> indexing attempt is correct!
    if (fIndexer_fMC_to_fEvent.count(mcPart->GetMother())) {
      // then, given that particle has a mother, load it
      idx_mother = fIndexer_fMC_to_fEvent[mcPart->GetMother()];
    } else {
      idx_mother = -2;
    }
  }  // closure: idx_mother stays at -1, because
  //   // particle has no mother <-> particle is primary
  // (artificial anti-sexaquark)
  // particle is primary && particle is signal
  // <->
  // comes from a sexaquark, that will be added on the first position
  if (fHasSexaquark && mcPart->MCStatusCode() == 6) {
    idx_mother = 0;
    // (artificial anti-sexaquark) add it to the daughters of the sexaquark
    // the "+1" is because of the inserted sexaquark
    fSexaquarkDaughters.push_back(fEvent.N_MCGen + 1);
  }
  // (2) assign
  fEvent.MC_Px.push_back(mcPart->Px());
  fEvent.MC_Py.push_back(mcPart->Py());
  fEvent.MC_Pz.push_back(mcPart->Pz());
  fEvent.MC_X.push_back(mcPart->Xv());
  fEvent.MC_Y.push_back(mcPart->Yv());
  fEvent.MC_Z.push_back(mcPart->Zv());
  fEvent.MC_PID.push_back(mcPart->PdgCode());
  fEvent.MC_Mother.push_back(idx_mother);
  fEvent.MC_Gen.push_back(generation);
  fEvent.MC_Status.push_back(mcPart->MCStatusCode());
  fEvent.MC_isSignal.push_back(isSignal);
  fEvent.N_MCGen++;
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCGen_InsertAntiSexaquark(TVector3 antilambda, TVector3 neutral_kaon) {
  //
  // Insert at first position an artificial MC generated anti-sexaquark into the Event_tt object
  //
  // (1) prepare
  // solve daughters indexing, add 1 because of the inserted sexaquark
  Int_t idx_first_daughter = fSexaquarkDaughters[0];
  Int_t idx_last_daughter = fSexaquarkDaughters[(Int_t)fSexaquarkDaughters.size() - 1];
  // (2) assign
  fEvent.MC_Px.insert(fEvent.MC_Px.begin(), antilambda.Px() + neutral_kaon.Px());
  fEvent.MC_Py.insert(fEvent.MC_Py.begin(), antilambda.Py() + neutral_kaon.Py());
  fEvent.MC_Pz.insert(fEvent.MC_Pz.begin(), antilambda.Pz() + neutral_kaon.Pz());
  fEvent.MC_X.insert(fEvent.MC_X.begin(), 0.);
  fEvent.MC_Y.insert(fEvent.MC_Y.begin(), 0.);
  fEvent.MC_Z.insert(fEvent.MC_Z.begin(), 0.);
  fEvent.MC_PID.insert(fEvent.MC_PID.begin(), -900000020);
  fEvent.MC_Mother.insert(fEvent.MC_Mother.begin(), -1);
  fEvent.MC_FirstDau.insert(fEvent.MC_FirstDau.begin(), idx_first_daughter);
  fEvent.MC_LastDau.insert(fEvent.MC_LastDau.begin(), idx_last_daughter);
  fEvent.MC_Gen.insert(fEvent.MC_Gen.begin(), 0);
  fEvent.MC_Status.insert(fEvent.MC_Status.begin(), 6);
  fEvent.MC_isSignal.insert(fEvent.MC_isSignal.begin(), kTRUE);
  fEvent.N_MCGen++;
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::MCRec_PushBack(AliESDtrack* track, Bool_t isDuplicate) {
  //
  // Add one more MC reconstructed track into the Event_tt object
  //
  // (1) prepare
  // check if rec. track is signal or not
  Float_t track_PID[2];
  track_PID[0] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  track_PID[1] = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  // get linked MC particle
  Int_t link_label = TMath::Abs(track->GetLabel());
  AliMCParticle* link = (AliMCParticle*)fMC->GetTrack(link_label);
  Int_t generation = 3;  // dummy int
  Bool_t isSignal = kFALSE;
  MCGen_VerifySourceAndGeneration(link, generation, isSignal);
  // (2) assign
  fEvent.Idx_True.push_back(fIndexer_fMC_to_fEvent[link_label]);
  fEvent.Rec_Px.push_back(track->Px());
  fEvent.Rec_Py.push_back(track->Py());
  fEvent.Rec_Pz.push_back(track->Pz());
  fEvent.Rec_Charge.push_back(track->Charge());
  fEvent.Rec_NSigmaPion.push_back(track_PID[0]);
  fEvent.Rec_NSigmaProton.push_back(track_PID[1]);
  fEvent.Rec_NClustersTPC.push_back(track->GetTPCncls());
  fEvent.Rec_isDuplicate.push_back(isDuplicate);
  fEvent.Rec_isSignal.push_back(isSignal);
  fEvent.N_MCRec++;
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::V0_PushBack(V0_tt common_format) {
  //
  // Add one more V0 into the Event_tt object
  //
  // (1) prepare
  // let's check if V0 is signal or not
  // define dummy int
  Int_t generation = 3;
  // positive component
  AliESDtrack* pos_track = fESD->GetTrack(common_format.Idx_Pos);
  Int_t pos_label = TMath::Abs(pos_track->GetLabel());
  AliMCParticle* pos = (AliMCParticle*)fMC->GetTrack(pos_label);
  Bool_t pos_isSignal = kFALSE;
  MCGen_VerifySourceAndGeneration(pos, generation, pos_isSignal);
  // negative component
  AliESDtrack* neg_track = fESD->GetTrack(common_format.Idx_Neg);
  Int_t neg_label = TMath::Abs(neg_track->GetLabel());
  AliMCParticle* neg = (AliMCParticle*)fMC->GetTrack(neg_label);
  Bool_t neg_isSignal = kFALSE;
  MCGen_VerifySourceAndGeneration(neg, generation, neg_isSignal);
  // get outcome
  Bool_t isSignal = pos_isSignal && neg_isSignal;
  // (2) assign
  fEvent.Idx_Pos.push_back(fIndexer_fTrack_to_fEvent[common_format.Idx_Pos]);
  fEvent.Idx_Neg.push_back(fIndexer_fTrack_to_fEvent[common_format.Idx_Neg]);
  fEvent.V0_Px.push_back(common_format.Px);
  fEvent.V0_Py.push_back(common_format.Py);
  fEvent.V0_Pz.push_back(common_format.Pz);
  fEvent.V0_X.push_back(common_format.X);
  fEvent.V0_Y.push_back(common_format.Y);
  fEvent.V0_Z.push_back(common_format.Z);
  fEvent.Pos_Px.push_back(common_format.Pos_Px);
  fEvent.Pos_Py.push_back(common_format.Pos_Py);
  fEvent.Pos_Pz.push_back(common_format.Pos_Pz);
  fEvent.Neg_Px.push_back(common_format.Neg_Px);
  fEvent.Neg_Py.push_back(common_format.Neg_Py);
  fEvent.Neg_Pz.push_back(common_format.Neg_Pz);
  fEvent.V0_isSignal.push_back(isSignal);
  fEvent.V0_E_asK0.push_back(common_format.E_asK0);
  fEvent.V0_E_asAL.push_back(common_format.E_asAL);
  fEvent.V0_couldBeK0.push_back(common_format.couldBeK0);
  fEvent.V0_couldBeAL.push_back(common_format.couldBeAL);
  fEvent.V0_onFlyStatus.push_back(common_format.onFlyStatus);
  fEvent.V0_Chi2.push_back(common_format.Chi2);
  fEvent.V0_DCA_Daughters.push_back(common_format.DCA_Daughters);
  fEvent.V0_IP_wrtPV.push_back(common_format.IP_wrtPV);
  fEvent.V0_CPA_wrtPV.push_back(common_format.CPA_wrtPV);
  fEvent.V0_ArmAlpha.push_back(common_format.ArmAlpha);
  fEvent.V0_ArmPt.push_back(common_format.ArmPt);
  fEvent.V0_DecayLength.push_back(common_format.DecayLength);
  fEvent.V0_DCA_wrtPV.push_back(common_format.DCA_wrtPV);
  fEvent.V0_isPrimary.push_back(common_format.isPrimary);
  fEvent.N_V0s++;
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::ClearEvent() {
  //
  // Clear all indices and vectors loaded in the Event_tt object
  //
  /* MC particles */
  fEvent.N_MCGen = 0;
  fEvent.MC_Px.clear();
  fEvent.MC_Py.clear();
  fEvent.MC_Pz.clear();
  fEvent.MC_X.clear();
  fEvent.MC_Y.clear();
  fEvent.MC_Z.clear();
  fEvent.MC_PID.clear();
  fEvent.MC_Mother.clear();
  fEvent.MC_FirstDau.clear();
  fEvent.MC_LastDau.clear();
  fEvent.MC_Gen.clear();
  fEvent.MC_Status.clear();
  fEvent.MC_isSignal.clear();
  /* MC Rec. Tracks */
  fEvent.N_MCRec = 0;
  fEvent.Idx_True.clear();
  fEvent.Rec_Px.clear();
  fEvent.Rec_Py.clear();
  fEvent.Rec_Pz.clear();
  fEvent.Rec_Charge.clear();
  fEvent.Rec_NSigmaPion.clear();
  fEvent.Rec_NSigmaProton.clear();
  fEvent.Rec_NClustersTPC.clear();
  fEvent.Rec_isDuplicate.clear();
  fEvent.Rec_isSignal.clear();
  /* V0s */
  fEvent.N_V0s = 0;
  fEvent.Idx_Pos.clear();
  fEvent.Idx_Neg.clear();
  fEvent.V0_Px.clear();
  fEvent.V0_Py.clear();
  fEvent.V0_Pz.clear();
  fEvent.V0_X.clear();
  fEvent.V0_Y.clear();
  fEvent.V0_Z.clear();
  fEvent.Pos_Px.clear();
  fEvent.Pos_Py.clear();
  fEvent.Pos_Pz.clear();
  fEvent.Neg_Px.clear();
  fEvent.Neg_Py.clear();
  fEvent.Neg_Pz.clear();
  fEvent.V0_isSignal.clear();
  fEvent.V0_E_asK0.clear();
  fEvent.V0_E_asAL.clear();
  fEvent.V0_couldBeK0.clear();
  fEvent.V0_couldBeAL.clear();
  fEvent.V0_onFlyStatus.clear();
  fEvent.V0_Chi2.clear();
  fEvent.V0_DCA_Daughters.clear();
  fEvent.V0_IP_wrtPV.clear();
  fEvent.V0_CPA_wrtPV.clear();
  fEvent.V0_ArmAlpha.clear();
  fEvent.V0_ArmPt.clear();
  fEvent.V0_DecayLength.clear();
  fEvent.V0_DCA_wrtPV.clear();
  fEvent.V0_isPrimary.clear();
}

/*** HIST OPERATIONS ***/

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::FillMCHistograms_Products(AliMCParticle* mcPart) {
  //
  // Fill the histograms that represent the generated MC final-state products
  //
  // fill PDG histogram
  fHistMC_PDG->Fill(mcPart->PdgCode());
  /*
  if (mcPart->PdgCode() == 321) {
    // Positive kaon
    fHistMC_Kp_P->Fill(mcPart->P());
    fHistMC_Kp_Pt->Fill(mcPart->Pt());
  } else if (mcPart->PdgCode() == 211) {
    // Positive pion
    fHistMC_Pip_P->Fill(mcPart->P());
    fHistMC_Pip_Pt->Fill(mcPart->Pt());
  } else if (mcPart->PdgCode() == -211) {
    // Negative pion
    fHistMC_Pim_P->Fill(mcPart->P());
    fHistMC_Pim_Pt->Fill(mcPart->Pt());
  } else if (mcPart->PdgCode() == -2212) {
    // Antiproton
    fHistMC_AntiP_P->Fill(mcPart->P());
    fHistMC_AntiP_Pt->Fill(mcPart->Pt());
  } else if (mcPart->PdgCode() == -2112) {
    // Antineutron
    fHistMC_AntiN_P->Fill(mcPart->P());
    fHistMC_AntiN_Pt->Fill(mcPart->Pt());
  } else if (mcPart->PdgCode() == 22) {
    // Gamma
    fHistMC_Gamma_P->Fill(mcPart->P());
    fHistMC_Gamma_Pt->Fill(mcPart->Pt());
  }
  */
}

//_____________________________________________________________________________
void AliAnalysisTaskSexaquark::FillMCHistograms_Sexaquark(std::vector<Int_t> fsProducts, Int_t fsCharge) {
  //
  // Fill the histograms that represent the generated sexaquarks
  // (at difference with FillHistograms_Sexaquark(), this function is general-purpose,
  // if MCComesFromSexaquark() works correctly, this one should work for every decay channel)
  //

  Double_t mcNucleon_M = kMassProton;
  if (fsCharge == 0) {
    mcNucleon_M = kMassNeutron;
  }

  // define nucleon at rest
  Double_t mcNucleon_Px = 0;
  Double_t mcNucleon_Py = 0;
  Double_t mcNucleon_Pz = 0;
  Double_t mcNucleon_E = mcNucleon_M;

  // define sexaquark kinematics
  Double_t mcSexaquark_Px = 0;
  Double_t mcSexaquark_Py = 0;
  Double_t mcSexaquark_Pz = 0;
  Double_t mcSexaquark_E = 0;

  // add four-momentum of each final-state particle
  AliMCParticle* fsPart;
  for (Int_t i = 0; i < (Int_t)fsProducts.size(); i++) {
    fsPart = (AliMCParticle*)fMC->GetTrack(fsProducts[i]);
    mcSexaquark_Px += fsPart->Px();
    mcSexaquark_Py += fsPart->Py();
    mcSexaquark_Pz += fsPart->Pz();
    mcSexaquark_E += fsPart->E();
  }

  // remove four-momentum of rest target
  mcSexaquark_Px -= 0.;
  mcSexaquark_Py -= 0.;
  mcSexaquark_Pz -= 0.;
  mcSexaquark_E -= mcNucleon_E;

  // reconstruct momentum of the sexaquark and fill histogram
  Double_t mcSexaquark_P = TMath::Sqrt(TMath::Power(mcSexaquark_Px, 2) + TMath::Power(mcSexaquark_Py, 2) + TMath::Power(mcSexaquark_Pz, 2));
  // fHistMC_Sexaquark_P->Fill(mcSexaquark_P);

  // reconstruct transverse-momentum of the sexaquark and fill histogram
  Double_t mcSexaquark_Pt = TMath::Sqrt(TMath::Power(mcSexaquark_Px, 2) + TMath::Power(mcSexaquark_Py, 2));
  // fHistMC_Sexaquark_Pt->Fill(mcSexaquark_Pt);

  // reconstruct mass of the sexaquark and fill histogram
  Double_t mcSexaquark_M;
  if (TMath::Power(mcSexaquark_E, 2) - TMath::Power(mcSexaquark_P, 2) >= 0) {
    mcSexaquark_M = TMath::Sqrt(TMath::Power(mcSexaquark_E, 2) - TMath::Power(mcSexaquark_P, 2));
    // fHistMC_Sexaquark_M->Fill(mcSexaquark_M);
    AliInfoF("M_S = %.6f", mcSexaquark_M);
  } else {
    AliInfo("M_S is OFF-SHELL");
  }
}
