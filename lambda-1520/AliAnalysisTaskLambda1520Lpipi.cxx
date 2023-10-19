#define HomogeneousField  // homogenous field in z direction, required by KFParticle

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

#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFVertex.h"

#include "AliAnalysisTaskLambda1520Lpipi.h"

#define DEBUG_MODE  // 1 : MC, 2 : Rec., 3 : V0s

class AliAnalysisTaskLambda1520Lpipi;

AliAnalysisTaskLambda1520Lpipi::AliAnalysisTaskLambda1520Lpipi()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      kMaxNSigma_Pion(0.),
      kMaxNSigma_Kaon(0.),
      kMaxNSigma_Proton(0.),
      fPDG(),
      fTree(0),
      fHist_Bookkeeper(0),
      fHist_TrueLambda_Pt(0),
      fHist_TrueLambda_Pz(0),
      fHist_TrueLambda_Radius(0),
      fHist_TrueAntiLambda_Pt(0),
      fHist_TrueAntiLambda_Pz(0),
      fHist_TrueAntiLambda_Radius(0) {}

AliAnalysisTaskLambda1520Lpipi::AliAnalysisTaskLambda1520Lpipi(const char* name, Bool_t IsMC)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      kMaxNSigma_Pion(0.),
      kMaxNSigma_Kaon(0.),
      kMaxNSigma_Proton(0.),
      fPDG(),
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

/*
 Destructor
 */
AliAnalysisTaskLambda1520Lpipi::~AliAnalysisTaskLambda1520Lpipi() {
  if (fOutputListOfTrees) {
    delete fOutputListOfTrees;
  }
  if (fOutputListOfHists) {
    delete fOutputListOfHists;
  }
}

/*
 Create output objects
 This function is called ONCE at the start of your analysis (RUNTIME)
 */
void AliAnalysisTaskLambda1520Lpipi::UserCreateOutputObjects() {

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

  AliInfo("Initializing AliAnalysisTaskLambda1520Lpipi...");
  AliInfo("INPUT OPTIONS:");
  AliInfoF(" >> IsMC           = %i", (Int_t)fIsMC);

  /*** Trees ***/

  fOutputListOfTrees = new TList();
  fOutputListOfTrees->SetOwner(kTRUE);

  fTree = new TTree("Events", "Tree of Events");
  // SetBranches();
  fOutputListOfTrees->Add(fTree);

  PostData(1, fOutputListOfTrees);

  /*** Histograms ***/

  fOutputListOfHists = new TList();
  fOutputListOfHists->SetOwner(kTRUE);

  fHist_Bookkeeper = new TH1I("Bookkeeper", "Bookkeeper", 7, 0., 7.);
  fOutputListOfHists->Add(fHist_Bookkeeper);

  fHist_TrueLambda_Pt = new TH1F("TrueLambda_Pt", "TrueLambda_Pt", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_TrueLambda_Pt);
  fHist_TrueLambda_Pz = new TH1F("TrueLambda_Pz", "TrueLambda_Pz", 100, 0., 20.);
  fOutputListOfHists->Add(fHist_TrueLambda_Pz);
  fHist_TrueLambda_Radius = new TH1F("TrueLambda_Radius", "TrueLambda_Radius", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_TrueLambda_Radius);
  fHist_TrueAntiLambda_Pt = new TH1F("TrueAntiLambda_Pt", "TrueAntiLambda_Pt", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_TrueAntiLambda_Pt);
  fHist_TrueAntiLambda_Pz = new TH1F("TrueAntiLambda_Pz", "TrueAntiLambda_Pz", 100, 0., 20.);
  fOutputListOfHists->Add(fHist_TrueAntiLambda_Pz);
  fHist_TrueAntiLambda_Radius = new TH1F("TrueAntiLambda_Radius", "TrueAntiLambda_Radius", 100, 0., 10.);
  fOutputListOfHists->Add(fHist_TrueAntiLambda_Radius);

  PostData(2, fOutputListOfHists);
}

/*
 Define cuts, depending on type of analysis.
 */
void AliAnalysisTaskLambda1520Lpipi::DefineCuts(TString cuts_option) {
  if (cuts_option == "Lambda1520->Lambda,pi,pi//Standard") {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    kMaxDCA_T1_T2 = 2.;  // cm
  }
  if (cuts_option == "AntiSexaquark,N->AntiLambda,K0S//Standard") {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    kMaxDCA_T1_T2 = 2.;  // cm
  }
  if (cuts_option == "AntiSexaquark,P->AntiLambda,K+//Standard") {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    kMaxDCA_T1_T2 = 2.;  // cm
  }
  if (cuts_option == "AntiSexaquark,P->AntiLambda,K+,pi-,pi+//Standard") {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    kMaxDCA_T1_T2 = 2.;  // cm
  }
  if (cuts_option == "AntiSexaquark,P->X,K+,K+//Standard") {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    kMaxDCA_T1_T2 = 2.;  // cm
  }
}

/*
 Main function, called per each event.
 */
void AliAnalysisTaskLambda1520Lpipi::UserExec(Option_t*) {

  // TString mode = "Lambda1520->Lambda,pi,pi//Standard";  // (pending) make this a global/input option...
  // (pending) don't forget that I'm currently asking only for signal particles...
  TString mode = "AntiSexaquark,N->AntiLambda,K0S//Standard";  // (pending) make this a global/input option...
  DefineCuts(mode);

  // containers
  std::set<Int_t> Indices_MCGen_FS;
  std::set<Int_t> Indices_MCGen_FS_Signal;  // FS: final-state

  std::vector<Int_t> idxPiPlusTracks, idxPiMinusTracks;
  std::vector<Int_t> idxKaonPlusTracks, idxKaonMinusTracks;
  std::vector<Int_t> idxProtonTracks, idxAntiProtonTracks;

  std::vector<KFParticle> kfLambdas;
  std::vector<std::vector<Int_t>> idxLambdaDaughters;
  std::vector<KFParticle> kfAntiLambdas;
  std::vector<std::vector<Int_t>> idxAntiLambdaDaughters;
  std::vector<KFParticle> kfPionPairs;
  std::vector<std::vector<Int_t>> idxPionPairDaughters;
  std::vector<KFParticle> kfNeutralShortKaons;
  std::vector<std::vector<Int_t>> idxNeutralShortKaonsDaughters;

  std::vector<Int_t> pdgTracks;

  std::vector<KFParticle> kfLambdas1520;
  std::vector<KFParticle> kfAntiLambdas1520;
  std::vector<KFParticle> kfAntiSexaquarks;

  // load MC generated event
  fMC = MCEvent();

  if (!fMC) {
    AliFatal("ERROR: AliMCEvent couldn't be found.");
  }

  // load reconstructed event
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!fESD) {
    AliFatal("ERROR: AliESDEvent couldn't be found.");
  }

  // load magnetic field
  fMagneticField = fESD->GetMagneticField();

  // initialize KFparticle
  KFParticle::SetField(fMagneticField);

  // load primary vertex
  // [note: const_cast removes the const qualifier from the pointer]
  fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());
#ifdef DEBUG_MODE
  AliInfoF("REC. PRIMARY VERTEX = (%8.3f, %8.3f, %8.3f)", fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
#endif

  ProcessMCGen(Indices_MCGen_FS, Indices_MCGen_FS_Signal);

  ProcessTracks(Indices_MCGen_FS_Signal,                //
                idxPiPlusTracks, idxPiMinusTracks,      //
                idxKaonPlusTracks, idxKaonMinusTracks,  //
                idxProtonTracks, idxAntiProtonTracks);

  if (mode == "Lambda1520->Lambda,pi,pi//Standard") {
    /* Lambda -> pi-,P */
    ReconstructV0s_KF(idxPiMinusTracks, idxProtonTracks, -211, 2212, kfLambdas, idxLambdaDaughters);
    /* AntiLambda -> AntiP,pi+ */
    ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
    /* Lambda(1520) -> X,pi-,pi+ */
    ReconstructV0s_KF(idxPiMinusTracks, idxPiPlusTracks, -211, 211, kfPionPairs, idxPionPairDaughters);
  }

  if (mode == "AntiSexaquark,N->AntiLambda,K0S//Standard") {
    /* AntiLambda -> AntiP,pi+ */
    ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
    /* K0S -> pi-,pi+ */
    ReconstructV0s_KF(idxPiMinusTracks, idxPiPlusTracks, -211, 211, kfNeutralShortKaons, idxNeutralShortKaonsDaughters);
  }

  // ProcessOfficialV0s(Indices_MCRec_Signal);

  if (mode == "Lambda1520->Lambda,pi,pi//Standard") {
    /* Lambda(1520) -> Lambda,pi-,pi+ */
    pdgTracks = {-211, 2212, -211, 211};
    Lambda1520Finder(kfLambdas, idxLambdaDaughters,      //
                     kfPionPairs, idxPionPairDaughters,  //
                     pdgTracks, kfLambdas1520);
    /* AntiLambda(1520) -> AntiLambda,pi-,pi+ */
    pdgTracks = {-2212, 211, -211, 211};
    Lambda1520Finder(kfAntiLambdas, idxAntiLambdaDaughters,  //
                     kfPionPairs, idxPionPairDaughters,      //
                     pdgTracks, kfAntiLambdas1520);
  }

  if (mode == "AntiSexaquark,N->AntiLambda,K0S//Standard") {
    /* AntiSexaquark,N -> AntiLambda,K0S */
    pdgTracks = {-2212, 211, -211, 211};
    SexaquarkFinder_ChannelA(kfAntiLambdas, idxAntiLambdaDaughters,               //
                             kfNeutralShortKaons, idxNeutralShortKaonsDaughters,  //
                             pdgTracks, kfAntiSexaquarks);
  }

  // (tree operations) after all variables have been assigned, fill tree
  // fTree->Fill();

  // stream the results the analysis of this event to the output manager
  PostData(1, fOutputListOfTrees);
  PostData(2, fOutputListOfHists);

  /* Clear containers */

  Indices_MCGen_FS.clear();
  Indices_MCGen_FS_Signal.clear();

  idxPiPlusTracks.clear();
  idxPiMinusTracks.clear();
  idxKaonPlusTracks.clear();
  idxKaonMinusTracks.clear();
  idxProtonTracks.clear();
  idxAntiProtonTracks.clear();

  kfLambdas.clear();
  idxLambdaDaughters.clear();

  kfAntiLambdas.clear();
  idxAntiLambdaDaughters.clear();

  kfPionPairs.clear();
  idxPionPairDaughters.clear();

  kfNeutralShortKaons.clear();
  idxNeutralShortKaonsDaughters.clear();

  pdgTracks.clear();

  kfLambdas1520.clear();
  kfAntiLambdas1520.clear();
  kfAntiSexaquarks.clear();
}

/*
 Search for possible contradictions in the input options
 */
void AliAnalysisTaskLambda1520Lpipi::CheckForInputErrors() {
  // (nothing yet)
}

/*
 Search for the relevant MC particles in a single event.
 - input: fMC
 - output: Indices_MCGen_FS, Indices_MCGen_FS_Signal
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessMCGen(std::set<Int_t>& Indices_MCGen_FS, std::set<Int_t>& Indices_MCGen_FS_Signal) {

  AliMCParticle* mcPart;
  AliMCParticle* mcDaughter;
  AliMCParticle* mcGrandDaughter;

  Int_t n_antilambda = 0;
  Int_t n_lambda = 0;
  Int_t n_pion = 0;
  Int_t n_piplus = 0;
  Int_t n_piminus = 0;

  Double_t radius;  // 2d-radius where (anti)lambda(1520) decays

  // (hists)
  if (fMC->GetNumberOfTracks()) fHist_Bookkeeper->Fill(0);

    // loop over MC gen. particles in a single event
#ifdef DEBUG_MODE
  AliInfo("   IDX    PID   NDAU FIRDAU LASDAU MOTHER ORIG");
#endif
  for (Int_t idx_mc = 0; idx_mc < fMC->GetNumberOfTracks(); idx_mc++) {

    mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

    // (cut)
    if (TMath::Abs(mcPart->PdgCode()) != 3124) continue;

    // (hists)
    if (mcPart->PdgCode() == 3124) fHist_Bookkeeper->Fill(1);
    if (mcPart->PdgCode() == -3124) fHist_Bookkeeper->Fill(4);

    // (cut)
    if (mcPart->GetNDaughters() != 3) continue;

#ifdef DEBUG_MODE
    AliInfoF("%6i %6i %6i %6i %6i %6i %8.3f %8.3f %8.3f", idx_mc, mcPart->PdgCode(), mcPart->GetNDaughters(), mcPart->GetDaughterFirst(),
             mcPart->GetDaughterLast(), mcPart->GetMother(), mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
#endif

    // loop over daughters
    for (Int_t idx_dau = mcPart->GetDaughterFirst(); idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {

      mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);

#ifdef DEBUG_MODE
      AliInfoF("%6i %6i %6i %6i %6i %6i %8.3f %8.3f %8.3f", idx_dau, mcDaughter->PdgCode(), mcDaughter->GetNDaughters(),
               mcDaughter->GetDaughterFirst(), mcDaughter->GetDaughterLast(), mcDaughter->GetMother(), mcDaughter->Xv(), mcDaughter->Yv(),
               mcDaughter->Zv());
#endif

      if (mcDaughter->Charge() != 0) Indices_MCGen_FS_Signal.insert(idx_dau);

      if (mcDaughter->PdgCode() == 3122) n_lambda++;
      if (mcDaughter->PdgCode() == -3122) n_antilambda++;
      if (mcDaughter->PdgCode() == 111 || mcDaughter->PdgCode() == 211 || mcDaughter->PdgCode() == -211) n_pion++;
      if (mcDaughter->PdgCode() == 211) n_piplus++;
      if (mcDaughter->PdgCode() == -211) n_piminus++;

      // check lambda's daughters -- loop over lambda(1520)'s granddaughters
      if (TMath::Abs(mcDaughter->PdgCode()) != 3122 || mcDaughter->GetNDaughters() == 0) continue;
      for (Int_t idx_grdau = mcDaughter->GetDaughterFirst(); idx_grdau <= mcDaughter->GetDaughterLast(); idx_grdau++) {

        mcGrandDaughter = (AliMCParticle*)fMC->GetTrack(idx_grdau);

        if (mcGrandDaughter->Charge() != 0) Indices_MCGen_FS_Signal.insert(idx_grdau);

#ifdef DEBUG_MODE
        AliInfoF("%6i %6i %6i %6i %6i %6i %8.3f %8.3f %8.3f", idx_grdau, mcGrandDaughter->PdgCode(), mcGrandDaughter->GetNDaughters(),
                 mcGrandDaughter->GetDaughterFirst(), mcGrandDaughter->GetDaughterLast(), mcGrandDaughter->GetMother(),
                 mcGrandDaughter->Xv(), mcGrandDaughter->Yv(), mcGrandDaughter->Zv());
#endif
      }  // end of loop over granddaughters

    }  // end of loop over daughters

    // use last daughter pointer to calculate the radius
    radius = TMath::Sqrt(mcDaughter->Xv() * mcDaughter->Xv() + mcDaughter->Yv() * mcDaughter->Yv());

    // (hists)
    if (n_lambda == 1 && n_pion == 2) {
      fHist_Bookkeeper->Fill(2);
      if (n_piplus == 1 && n_piminus == 1) {
        fHist_Bookkeeper->Fill(3);
        fHist_TrueLambda_Pt->Fill(mcPart->Pt());
        fHist_TrueLambda_Pz->Fill(mcPart->Pz());
        fHist_TrueLambda_Radius->Fill(radius);
      }
    }
    if (n_antilambda == 1 && n_pion == 2) {
      fHist_Bookkeeper->Fill(5);
      if (n_piplus == 1 && n_piminus == 1) {
        fHist_Bookkeeper->Fill(6);
        fHist_TrueAntiLambda_Pt->Fill(mcPart->Pt());
        fHist_TrueAntiLambda_Pz->Fill(mcPart->Pz());
        fHist_TrueAntiLambda_Radius->Fill(radius);
      }
    }

    // reset counters
    n_antilambda = 0;
    n_lambda = 0;
    n_pion = 0;
    n_piplus = 0;
    n_piminus = 0;
  }  // end of loop over MC particles
}

/*
 Loop over the reconstructed tracks in a single event.
 IMPORTANT: don't trust track->GetID()
 - uses: fESD, fMC
 - input: Indices_MCGen_FS_Signal (pending)
 - output: idxPiPlusTracks,idxPiMinusTracks,idxKaonPlusTracks,idxKaonMinusTracks,idxProtonTracks,idxAntiProtonTracks
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessTracks(std::set<Int_t> Indices_MCGen_FS_Signal,                                    // (pending)
                                                   std::vector<Int_t>& idxPiPlusTracks, std::vector<Int_t>& idxPiMinusTracks,  //
                                                   std::vector<Int_t>& idxKaonPlusTracks, std::vector<Int_t>& idxKaonMinusTracks,  //
                                                   std::vector<Int_t>& idxProtonTracks, std::vector<Int_t>& idxAntiProtonTracks) {

  AliESDtrack* track;

  Int_t idx_true;
  AliMCParticle* truePart;

  Float_t n_sigma_pion;
  Float_t n_sigma_kaon;
  Float_t n_sigma_proton;

  Int_t counter = 0;  // (debug)

  // loop over MC rec. tracks in a single event
#ifdef DEBUG_MODE
  AliInfo("     IDX  trueIDX  truePID nsigmapi  nsigmaP");
#endif
  for (Int_t idx_track = 0; idx_track < fESD->GetNumberOfTracks(); idx_track++) {

    track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));
    idx_true = TMath::Abs(track->GetLabel());

    // (cut) (pending)
    if (!Indices_MCGen_FS_Signal.count(idx_true)) {
      continue;
    }

    counter++;  // (debug)

    n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    // (cut) particle identification
    if (TMath::Abs(n_sigma_pion) > kMaxNSigma_Pion &&  //
        TMath::Abs(n_sigma_kaon) > kMaxNSigma_Kaon &&  //
        TMath::Abs(n_sigma_proton) > kMaxNSigma_Proton) {
      continue;
    }

    if (TMath::Abs(n_sigma_pion) < kMaxNSigma_Pion) {
      if (track->Charge() > 0.)
        idxPiPlusTracks.push_back(idx_track);
      else
        idxPiMinusTracks.push_back(idx_track);
    }

    if (TMath::Abs(n_sigma_kaon) < kMaxNSigma_Kaon) {
      if (track->Charge() > 0.)
        idxKaonPlusTracks.push_back(idx_track);
      else
        idxKaonMinusTracks.push_back(idx_track);
    }

    if (TMath::Abs(n_sigma_proton) < kMaxNSigma_Proton) {
      if (track->Charge() > 0.)
        idxProtonTracks.push_back(idx_track);
      else
        idxAntiProtonTracks.push_back(idx_track);
    }

#ifdef DEBUG_MODE
    truePart = (AliMCParticle*)fMC->GetTrack(idx_true);
    AliInfoF("%8i %8i %8i %8.2f %8.2f", idx_track, idx_true, truePart->PdgCode(), n_sigma_pion, n_sigma_proton);
#endif
  }  // end of loop over tracks

#ifdef DEBUG_MODE
  if (counter >= 4) AliInfo("IT HAPPENED!");
#endif
}

/*
 Loop over all the V0s found by the Official Offline V0 Finder
 - uses: fESD, fMagneticField
 - input: Indices_MCRec_Signal
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessOfficialV0s(std::vector<Int_t> Indices_MCRec_Signal) {

  // define variables & objects
  AliESDv0* V0;

  Int_t idx_pos;
  Int_t idx_neg;

  AliESDtrack* neg_track;
  AliESDtrack* pos_track;
  Float_t nsigma_pos_pion;
  Float_t nsigma_neg_pion;
  Float_t nsigma_antiproton;

  Double_t V0_X, V0_Y, V0_Z;
  Double_t N_Px, N_Py, N_Pz;
  Double_t P_Px, P_Py, P_Pz;

  Double_t pos_energy;
  Double_t neg_energy_asK0;
  Double_t neg_energy_asAL;

  Int_t idx_pos_true;
  Int_t idx_neg_true;
  AliMCParticle* mcPosTrue;
  AliMCParticle* mcNegTrue;

  // get primary vertex of this event
  Double_t PV[3];
  fPrimaryVertex->GetXYZ(PV);

  Double_t b = fMagneticField;
  Double_t xn, xp;
  Float_t dz[2];
  Double_t dca_neg_pv;
  Double_t dca_pos_pv;

  // loop over V0s
  for (Int_t idx_v0 = 0; idx_v0 < fESD->GetNumberOfV0s(); idx_v0++) {

    V0 = fESD->GetV0(idx_v0);

    // (cut) remove on-the-fly V0s
    if (V0->GetOnFlyStatus()) {
      continue;
    }

    idx_pos = V0->GetPindex();
    idx_neg = V0->GetNindex();

    // (cut) if any daughter (neg. or pos.) is not in the indexer
    if (std::find(Indices_MCRec_Signal.begin(), Indices_MCRec_Signal.end(), idx_pos) == Indices_MCRec_Signal.end() ||
        std::find(Indices_MCRec_Signal.begin(), Indices_MCRec_Signal.end(), idx_neg) == Indices_MCRec_Signal.end()) {
      continue;
    }

    neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_pos));
    pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_neg));

    neg_track->GetDZ(PV[0], PV[1], PV[2], b, dz);
    dca_neg_pv = TMath::Sqrt(dz[0] * dz[0] + dz[1] * dz[1]);

    pos_track->GetDZ(PV[0], PV[1], PV[2], b, dz);
    dca_pos_pv = TMath::Sqrt(dz[0] * dz[0] + dz[1] * dz[1]);

    AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;
    // Double_t xn, xp;
    if (Preoptimize(ntp, ptp, &xn, &xp, b)) {
      // Move tracks to a better position if that helps
      nt.PropagateTo(xn, b);
      pt.PropagateTo(xp, b);
    }

#ifdef DEBUG_MODE
    AliInfoF("%i %i dca(neg,pos) = %5.3f, dca(neg,pv) = %5.3f, dca(pos, pv) = %5.3f",  //
             idx_neg,                                                                  //
             idx_pos,                                                                  //
             ntp->GetDCA(ptp, b, xn, xp),                                              //
             dca_neg_pv,                                                               //
             dca_pos_pv);
    //  V0->GetDcaV0Daughters(), dca_neg_pv, dca_pos_pv);
#endif
  }  // end of loop over V0s
}

/*
 Find all V0s.
 Store all lambdas and kaons with no ambiguity.
 - no dependence on PID, nor ambiguity on measurements
 - only dependence: if they were reconstructed or not
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessTrueV0s() {

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

  // get primary vertex of this event
  Double_t PV[3];
  fPrimaryVertex->GetXYZ(PV);

  // loop over MC gen. mothers
  for (Int_t idx_mother = 0; idx_mother < fMC->GetNumberOfTracks(); idx_mother++) {

    mc_mother = (AliMCParticle*)fMC->GetTrack(idx_mother);

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
    // if (fMap_Duplicates.count(mc_mother->GetDaughterFirst()) == 0 || fMap_Duplicates.count(mc_mother->GetDaughterLast()) == 0) {
    // [note: count() in std::map, returns the number of elements matching specific key]
    // continue;
    // }

    // identify daughters as positive/negative
    if (mc_firstdau->PdgCode() > 0 && mc_lastdau->PdgCode() < 0) {
      idx_pos_mc = mc_mother->GetDaughterFirst();
      idx_neg_mc = mc_mother->GetDaughterLast();
    } else {
      idx_pos_mc = mc_mother->GetDaughterLast();
      idx_neg_mc = mc_mother->GetDaughterFirst();
    }

    // do combinations
    /* for (Int_t idx_pos_track : fMap_Duplicates[idx_pos_mc]) {
      for (Int_t idx_neg_track : fMap_Duplicates[idx_neg_mc]) {

        // load MC particle depending if positive/negative
        mc_pos = (AliMCParticle*)fMC->GetTrack(idx_pos_mc);
        mc_neg = (AliMCParticle*)fMC->GetTrack(idx_neg_mc);

        // calculate energy (assuming positive pion)
        pos_energy =
            TMath::Sqrt(mc_pos->Px() * mc_pos->Px() + mc_pos->Py() * mc_pos->Py() + mc_pos->Pz() * mc_pos->Pz() + kMassPion * kMassPion);
        // assuming K0 -> pi+ pi-
        neg_energy_asK0 =
            TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() + kMassPion * kMassPion);
        // assuming anti-lambda -> pi+ anti-proton
        neg_energy_asAL = TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() +
                                      kMassProton * kMassProton);

      }  // end of loop over negative tracks
    }    // end of loop over positive tracks
 */
  }  // end of loop over MC gen. mothers
}

/*
 Find all true V0s which had both of their daughters reconstructed.
 - uses: fESD, fMC, fMagneticField
 - input: idxNegativeTracks, idxPositiveTracks, pdgTrackNeg, pdgTrackPos
 - output: kfV0s, idxDaughters
 */
void AliAnalysisTaskLambda1520Lpipi::ReconstructV0s_KF(std::vector<Int_t> idxNegativeTracks,  //
                                                       std::vector<Int_t> idxPositiveTracks,  //
                                                       Int_t pdgTrackNeg, Int_t pdgTrackPos,  //
                                                       std::vector<KFParticle>& kfV0s,        //
                                                       std::vector<std::vector<Int_t>>& idxDaughters) {
  AliESDtrack* esdTrackNeg;
  AliESDtrack* esdTrackPos;

  AliMCParticle* mcPartDau1;
  AliMCParticle* mcPartDau2;

  AliMCParticle* mcMother1;
  AliMCParticle* mcMother2;

  Double_t xthis, xp;
  Double_t impar[2];
  Bool_t fStatusNegPropToDCA;
  Bool_t fStatusPosPropToDCA;

  TLorentzVector lvTrackNeg;
  TLorentzVector lvTrackPos;
  TLorentzVector lvV0;

  std::vector<Int_t> aux_idx_vector;

#ifdef DEBUG_MODE
  // AliInfo("     mom  mom_pdg   track");
#endif
  for (Int_t& idxTrackNeg : idxNegativeTracks) {
    for (Int_t& idxTrackPos : idxPositiveTracks) {

      // (protection)
      if (idxTrackNeg == idxTrackPos) continue;

      esdTrackNeg = static_cast<AliESDtrack*>(fESD->GetTrack(idxTrackNeg));
      esdTrackPos = static_cast<AliESDtrack*>(fESD->GetTrack(idxTrackPos));

#ifdef DEBUG_MODE
      mcPartDau1 = (AliMCParticle*)fMC->GetTrack(TMath::Abs(esdTrackNeg->GetLabel()));
      mcPartDau2 = (AliMCParticle*)fMC->GetTrack(TMath::Abs(esdTrackPos->GetLabel()));

      if (mcPartDau1->GetMother() != -1) mcMother1 = (AliMCParticle*)fMC->GetTrack(mcPartDau1->GetMother());
      if (mcPartDau2->GetMother() != -1) mcMother2 = (AliMCParticle*)fMC->GetTrack(mcPartDau2->GetMother());

        // AliInfoF("%8i %8i %8i", mcPartDau1->GetMother(), mcMother1->PdgCode(), idxTrackNeg);
        // AliInfoF("%8.3f %8.3f %8.3f", mcPartDau1->Xv(), mcPartDau1->Yv(), mcPartDau1->Zv());

        // AliInfoF("%8i %8i %8i", mcPartDau2->GetMother(), mcMother2->PdgCode(), idxTrackPos);
        // AliInfoF("%8.3f %8.3f %8.3f", mcPartDau2->Xv(), mcPartDau2->Yv(), mcPartDau2->Zv());
#endif

      // (cut) PENDING
      // if (mcPartDau1->GetMother() != mcPartDau2->GetMother()) continue;

      // (cut)
      if (esdTrackNeg->GetDCA(esdTrackPos, fMagneticField, xthis, xp) > kMaxDCA_T1_T2) {
        continue;
      }

      /* Kalman Filter */

      KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
      KFParticle kfDaughterNeg = CreateKFParticle(*esdTrackNeg, pdgTrackNeg, mcPartDau1->Charge());
      KFParticle kfDaughterPos = CreateKFParticle(*esdTrackPos, pdgTrackPos, mcPartDau2->Charge());

      KFParticle kfV0;
      kfV0.SetConstructMethod(2);
      kfV0.AddDaughter(kfDaughterNeg);
      kfV0.AddDaughter(kfDaughterPos);
      kfV0.TransportToDecayVertex();

      // get V0 vertex
      Double_t V0Vertex[3] = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
      Double_t V0VertexErr[3] = {kfV0.GetErrX(), kfV0.GetErrY(), kfV0.GetErrZ()};
      AliESDVertex* esdV0Vertex = new AliESDVertex(V0Vertex, V0VertexErr);

      // clone ESD track object to be able propagate it (~modify it)
      AliESDtrack cloneTrackNeg(*esdTrackNeg);
      AliESDtrack cloneTrackPos(*esdTrackPos);

      // propagate V0 daughters to their DCA w.r.t. V0 vertex
      fStatusNegPropToDCA = cloneTrackNeg.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar) == 1;
      fStatusPosPropToDCA = cloneTrackPos.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar) == 1;

      // reconstruct mother mass
      lvTrackNeg.SetXYZM(cloneTrackNeg.Px(), cloneTrackNeg.Py(), cloneTrackNeg.Pz(), fPDG.GetParticle(pdgTrackNeg)->Mass());
      lvTrackPos.SetXYZM(cloneTrackPos.Px(), cloneTrackPos.Py(), cloneTrackPos.Pz(), fPDG.GetParticle(pdgTrackPos)->Mass());
      lvV0 = lvTrackNeg + lvTrackPos;

      // (pending) add cuts to V0s!

#ifdef DEBUG_MODE
      AliInfoF("%8.3f %8.3f %8.3f %8.3f", kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), lvV0.M());
      AliInfoF("%i %i dca(d1,pv) = %5.3f, dca(d2,pv) = %5.3f, dca(d1,V0) = %5.3f, dca(d2,V0) = %5.3f, dca(V0,pv) = %5.3f",  //
               idxTrackNeg,                                                                                                 //
               idxTrackPos,                                                                                                 //
               kfDaughterNeg.GetDistanceFromVertex(kfPrimaryVertex),                                                        //
               kfDaughterPos.GetDistanceFromVertex(kfPrimaryVertex),                                                        //
               (Float_t)kfDaughterNeg.GetDistanceFromVertex(kfV0),                                                          //
               (Float_t)kfDaughterPos.GetDistanceFromVertex(kfV0),
               Calculate_LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(),        //
                                      kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(),  //
                                      fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ()));
#endif

      // finally, store info
      kfV0s.push_back(kfV0);
      aux_idx_vector.push_back(idxTrackNeg);
      aux_idx_vector.push_back(idxTrackPos);
      idxDaughters.push_back(aux_idx_vector);

      // clear container
      aux_idx_vector.clear();
    }  // end of loop over pos. tracks
  }    // end of loop over neg. tracks
}

/*
 Find all Lambda1520 candidates.
 // (pending) add more description
 - input: kfFirstV0s, idxFirstV0Daughters, kfSecondV0s, idxSecondV0Daughters, pdgDaughters
 - output: kfLambda1520s
 */
void AliAnalysisTaskLambda1520Lpipi::Lambda1520Finder(std::vector<KFParticle> kfFirstV0s,                    //
                                                      std::vector<std::vector<Int_t>> idxFirstV0Daughters,   //
                                                      std::vector<KFParticle> kfSecondV0s,                   //
                                                      std::vector<std::vector<Int_t>> idxSecondV0Daughters,  //
                                                      std::vector<Int_t> pdgDaughters,                       //
                                                      std::vector<KFParticle>& kfLambdas1520) {

  Double_t impar[2];
  Bool_t fStatusT0PropToDCA;
  Bool_t fStatusT1PropToDCA;
  Bool_t fStatusT2PropToDCA;
  Bool_t fStatusT3PropToDCA;

  Int_t idxV0A_NegDau, idxV0A_PosDau;
  Int_t idxV0B_NegDau, idxV0B_PosDau;

  TLorentzVector lvTrack0, lvTrack1, lvTrack2, lvTrack3;
  TLorentzVector lvLambda1520;

  Double_t dca;

#ifdef DEBUG_MODE
  AliInfo("               Mass       Vx       Vy       Vz     Chi2");
#endif
  for (Int_t idxV0A = 0; idxV0A < (Int_t)kfFirstV0s.size(); idxV0A++) {
    for (Int_t idxV0B = 0; idxV0B < (Int_t)kfSecondV0s.size(); idxV0B++) {

      idxV0A_NegDau = idxFirstV0Daughters[idxV0A][0];
      idxV0A_PosDau = idxFirstV0Daughters[idxV0A][1];
      idxV0B_NegDau = idxSecondV0Daughters[idxV0B][0];
      idxV0B_PosDau = idxSecondV0Daughters[idxV0B][1];

      AliInfoF("%i %i %i %i", idxV0A_NegDau, idxV0A_PosDau, idxV0B_NegDau, idxV0B_PosDau);

      // (protection) avoid repetition of daughters
      if (idxV0A_NegDau == idxV0B_NegDau || idxV0A_PosDau == idxV0B_PosDau) continue;

      AliESDtrack* esdTrack0 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_NegDau));  // not used
      AliESDtrack* esdTrack1 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_PosDau));  // not used
      AliESDtrack* esdTrack2 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_NegDau));
      AliESDtrack* esdTrack3 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_PosDau));

      /* Kalman Filter */

      KFParticle kfV0A = kfFirstV0s[idxV0A];
      KFParticle kfV0B = kfSecondV0s[idxV0B];  // not used

      KFParticle kfTrack0 = CreateKFParticle(*esdTrack0, pdgDaughters[0], -1);  // not used
      KFParticle kfTrack1 = CreateKFParticle(*esdTrack1, pdgDaughters[1], 1);   // not used
      KFParticle kfTrack2 = CreateKFParticle(*esdTrack2, pdgDaughters[2], -1);
      KFParticle kfTrack3 = CreateKFParticle(*esdTrack3, pdgDaughters[3], 1);

      KFParticle kfLambda1520;
      kfLambda1520.SetConstructMethod(2);
      kfLambda1520.AddDaughter(kfV0A);
      kfLambda1520.AddDaughter(kfTrack2);
      kfLambda1520.AddDaughter(kfTrack3);
      kfLambda1520.TransportToDecayVertex();

      // clone ESD tracks
      AliESDtrack cloneTrack0(*esdTrack0);
      AliESDtrack cloneTrack1(*esdTrack1);
      AliESDtrack cloneTrack2(*esdTrack2);
      AliESDtrack cloneTrack3(*esdTrack3);

      // get V0 vertex
      Double_t V0Vertex[3] = {kfV0A.GetX(), kfV0A.GetY(), kfV0A.GetZ()};
      Double_t V0VertexErr[3] = {kfV0A.GetErrX(), kfV0A.GetErrY(), kfV0A.GetErrZ()};
      AliESDVertex* esdV0Vertex = new AliESDVertex(V0Vertex, V0VertexErr);
      // propagate V0 daughters to their DCA w.r.t. V0 vertex
      fStatusT0PropToDCA = cloneTrack0.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar) == 1;
      fStatusT1PropToDCA = cloneTrack1.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar) == 1;

      // (question) should I find a way to propagate the V0 to the SV?

      // get Secondary Vertex
      Double_t SecVertex[3] = {kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ()};
      Double_t SecVertexErr[3] = {kfLambda1520.GetErrX(), kfLambda1520.GetErrY(), kfLambda1520.GetErrZ()};
      AliESDVertex* esdSecVertex = new AliESDVertex(SecVertex, SecVertexErr);
      // propagate remaining particles to their DCA w.r.t. Secondary Vertex
      fStatusT2PropToDCA = cloneTrack2.PropagateToDCA(esdSecVertex, fMagneticField, 10, impar) == 1;
      fStatusT3PropToDCA = cloneTrack3.PropagateToDCA(esdSecVertex, fMagneticField, 10, impar) == 1;

      /* Reconstruct Lambda(1520) */

      lvTrack0.SetXYZM(cloneTrack0.Px(), cloneTrack0.Py(), cloneTrack0.Pz(), fPDG.GetParticle(pdgDaughters[0])->Mass());
      lvTrack1.SetXYZM(cloneTrack1.Px(), cloneTrack1.Py(), cloneTrack1.Pz(), fPDG.GetParticle(pdgDaughters[1])->Mass());
      lvTrack2.SetXYZM(cloneTrack2.Px(), cloneTrack2.Py(), cloneTrack2.Pz(), fPDG.GetParticle(pdgDaughters[2])->Mass());
      lvTrack3.SetXYZM(cloneTrack3.Px(), cloneTrack3.Py(), cloneTrack3.Pz(), fPDG.GetParticle(pdgDaughters[3])->Mass());

      lvLambda1520 = lvTrack0 + lvTrack1 + lvTrack2 + lvTrack3;

      dca = Calculate_LinePointDCA(lvLambda1520.Px(), lvLambda1520.Py(), lvLambda1520.Pz(),        //
                                   kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ(),  //
                                   fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

      // (cut) (pending) DCA w.r.t. PV
      if (dca > 0.1) continue;

#ifdef DEBUG_MODE
      AliInfoF("   -     - %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f", lvLambda1520.M(), kfLambda1520.GetX(), kfLambda1520.GetY(),
               kfLambda1520.GetZ(), kfLambda1520.Chi2() / (Double_t)kfLambda1520.GetNDF(), dca);
#endif

      kfLambdas1520.push_back(kfLambda1520);
    }
  }
}

/*
 Reconstruct AntiSexaquark candidates via the Reaction Channel A:
 AntiSexaquark,N -> AntiLambda,K0S -> AntiP,pi+,pi-,pi+
 - input: kfFirstV0s, idxFirstV0Daughters, kfSecondV0s, idxSecondV0Daughters, pdgDaughters
 - output: kfAntiSexaquarks
 */
void AliAnalysisTaskLambda1520Lpipi::SexaquarkFinder_ChannelA(std::vector<KFParticle> kfFirstV0s,                    //
                                                              std::vector<std::vector<Int_t>> idxFirstV0Daughters,   //
                                                              std::vector<KFParticle> kfSecondV0s,                   //
                                                              std::vector<std::vector<Int_t>> idxSecondV0Daughters,  //
                                                              std::vector<Int_t> pdgDaughters,                       //
                                                              std::vector<KFParticle>& kfAntiSexaquarks) {

  // in this case, the interacting nucleon is a neutron
  const Int_t pdgStruckNucleon = 2112;

  Double_t impar[2];
  Bool_t fStatusT0PropToDCA;
  Bool_t fStatusT1PropToDCA;
  Bool_t fStatusT2PropToDCA;
  Bool_t fStatusT3PropToDCA;

  Int_t idxV0A_NegDau, idxV0A_PosDau;
  Int_t idxV0B_NegDau, idxV0B_PosDau;

  TLorentzVector lvNucleon;
  TLorentzVector lvTrack0, lvTrack1, lvTrack2, lvTrack3;
  TLorentzVector lvAntiSexaquark;

#ifdef DEBUG_MODE
  AliInfo("               Mass       Vx       Vy       Vz     Chi2");
#endif
  for (Int_t idxV0A = 0; idxV0A < (Int_t)kfFirstV0s.size(); idxV0A++) {
    for (Int_t idxV0B = 0; idxV0B < (Int_t)kfSecondV0s.size(); idxV0B++) {

      idxV0A_NegDau = idxFirstV0Daughters[idxV0A][0];
      idxV0A_PosDau = idxFirstV0Daughters[idxV0A][1];
      idxV0B_NegDau = idxSecondV0Daughters[idxV0B][0];
      idxV0B_PosDau = idxSecondV0Daughters[idxV0B][1];

      AliInfoF("%i %i %i %i", idxV0A_NegDau, idxV0A_PosDau, idxV0B_NegDau, idxV0B_PosDau);

      // (protection) avoid repetition of daughters
      if (idxV0A_NegDau == idxV0B_NegDau || idxV0A_PosDau == idxV0B_PosDau) continue;

      AliESDtrack* esdTrack0 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_NegDau));  // not used
      AliESDtrack* esdTrack1 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_PosDau));  // not used
      AliESDtrack* esdTrack2 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_NegDau));
      AliESDtrack* esdTrack3 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_PosDau));

      /* Kalman Filter */

      KFParticle kfV0A = kfFirstV0s[idxV0A];
      KFParticle kfV0B = kfSecondV0s[idxV0B];

      KFParticle kfTrack0 = CreateKFParticle(*esdTrack0, pdgDaughters[0], -1);  // not used
      KFParticle kfTrack1 = CreateKFParticle(*esdTrack1, pdgDaughters[1], 1);   // not used
      KFParticle kfTrack2 = CreateKFParticle(*esdTrack2, pdgDaughters[2], -1);
      KFParticle kfTrack3 = CreateKFParticle(*esdTrack3, pdgDaughters[3], 1);

      KFParticle kfAntiSexaquark;
      kfAntiSexaquark.SetConstructMethod(2);
      kfAntiSexaquark.AddDaughter(kfV0A);
      kfAntiSexaquark.AddDaughter(kfV0B);
      kfAntiSexaquark.TransportToDecayVertex();

      // clone ESD tracks
      AliESDtrack cloneTrack0(*esdTrack0);
      AliESDtrack cloneTrack1(*esdTrack1);
      AliESDtrack cloneTrack2(*esdTrack2);
      AliESDtrack cloneTrack3(*esdTrack3);

      // get V0A vertex
      Double_t V0A_Vertex[3] = {kfV0A.GetX(), kfV0A.GetY(), kfV0A.GetZ()};
      Double_t V0A_VertexErr[3] = {kfV0A.GetErrX(), kfV0A.GetErrY(), kfV0A.GetErrZ()};
      AliESDVertex* esdV0A_Vertex = new AliESDVertex(V0A_Vertex, V0A_VertexErr);
      // propagate V0 daughters to their DCA w.r.t. V0 vertex
      fStatusT0PropToDCA = cloneTrack0.PropagateToDCA(esdV0A_Vertex, fMagneticField, 10., impar) == 1;
      fStatusT1PropToDCA = cloneTrack1.PropagateToDCA(esdV0A_Vertex, fMagneticField, 10., impar) == 1;

      // get V0B vertex
      Double_t V0B_Vertex[3] = {kfV0B.GetX(), kfV0B.GetY(), kfV0B.GetZ()};
      Double_t V0B_VertexErr[3] = {kfV0B.GetErrX(), kfV0B.GetErrY(), kfV0B.GetErrZ()};
      AliESDVertex* esdV0B_Vertex = new AliESDVertex(V0B_Vertex, V0B_VertexErr);
      // propagate V0 daughters to their DCA w.r.t. V0 vertex
      fStatusT2PropToDCA = cloneTrack0.PropagateToDCA(esdV0B_Vertex, fMagneticField, 10., impar) == 1;
      fStatusT3PropToDCA = cloneTrack1.PropagateToDCA(esdV0B_Vertex, fMagneticField, 10., impar) == 1;

      // get Secondary Vertex
      Double_t SecVertex[3] = {kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ()};
      Double_t SecVertexErr[3] = {kfAntiSexaquark.GetErrX(), kfAntiSexaquark.GetErrY(), kfAntiSexaquark.GetErrZ()};
      AliESDVertex* esdSecVertex = new AliESDVertex(SecVertex, SecVertexErr);

      // (question) should I find a way to propagate the V0s to the SV?

      /* Reconstruct AntiSexaquark */

      lvNucleon.SetXYZM(0., 0., 0., fPDG.GetParticle(pdgStruckNucleon)->Mass());
      lvTrack0.SetXYZM(cloneTrack0.Px(), cloneTrack0.Py(), cloneTrack0.Pz(), fPDG.GetParticle(pdgDaughters[0])->Mass());
      lvTrack1.SetXYZM(cloneTrack1.Px(), cloneTrack1.Py(), cloneTrack1.Pz(), fPDG.GetParticle(pdgDaughters[1])->Mass());
      lvTrack2.SetXYZM(cloneTrack2.Px(), cloneTrack2.Py(), cloneTrack2.Pz(), fPDG.GetParticle(pdgDaughters[2])->Mass());
      lvTrack3.SetXYZM(cloneTrack3.Px(), cloneTrack3.Py(), cloneTrack3.Pz(), fPDG.GetParticle(pdgDaughters[3])->Mass());

      lvAntiSexaquark = lvTrack0 + lvTrack1 + lvTrack2 + lvTrack3 - lvNucleon;

      // (cut) (pending) DCA w.r.t. PV
      // dca = Calculate_LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(),        //
                                  //  kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),  //
                                  //  fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
      // // if (dca > 0.1) continue;

#ifdef DEBUG_MODE
      AliInfoF("   -     - %8.3f %8.3f %8.3f %8.3f %8.3f", lvAntiSexaquark.M(), kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
               kfAntiSexaquark.GetZ(), kfAntiSexaquark.Chi2() / (Double_t)kfAntiSexaquark.GetNDF());
#endif

      kfAntiSexaquarks.push_back(kfAntiSexaquark);
    }
  }
}

/*** CUSTOM OFFLINE V0 FINDER ***/

/*
 This function pre-optimizes a two-track pair in the XY plane
 and provides two X values for the tracks if successful
 [copied from AliRoot/STEER/ESD/AliV0vertexer.cxx]
 */
Bool_t AliAnalysisTaskLambda1520Lpipi::Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn,
                                                   Double_t* lPreprocessxp, const Double_t b) {
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

/*
 Copied from AliV0ReaderV1::GetHelixCenter
 Get Center of the helix track parametrization
 [based on AliRoot/STEER/ESD/AliV0vertexer.cxx]
 */
void AliAnalysisTaskLambda1520Lpipi::GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b) {

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

/*** UTILITIES ***/

/*
 Calculate transverse impact parameter w.r.t. Primary Vertex
 */
Float_t AliAnalysisTaskLambda1520Lpipi::MCRec_GetImpactParameter(AliESDtrack* track) {

  Double_t PV[3];
  fPrimaryVertex->GetXYZ(PV);

  return (Float_t)track->GetD(PV[0], PV[1], fMagneticField);
}

/*
 Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 This function stores the point of closest approach, and returns its distance to the PV
 [more info at https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line]
 // (question) should I change the units? nope, it's the same
 */
Double_t AliAnalysisTaskLambda1520Lpipi::Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                                Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                                Double_t PV_X, Double_t PV_Y, Double_t PV_Z) {
  TVector3 V0_Momentum(V0_Px, V0_Py, V0_Pz);

  TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
  TVector3 PrimVertex(PV_X, PV_Y, PV_Z);

  TVector3 CrossProduct = (PrimVertex - V0Vertex).Cross(V0_Momentum);

  return CrossProduct.Mag() / V0_Momentum.Mag();
}

/*
 Correct initialization of a KFParticle
 [copied from AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx]
 */
KFParticle AliAnalysisTaskLambda1520Lpipi::CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge) {

  Double_t fP[6];
  track.GetXYZ(fP);
  track.PxPyPz(fP + 3);

  Int_t fQ = track.Charge() * TMath::Abs(charge);
  fP[3] *= TMath::Abs(charge);
  fP[4] *= TMath::Abs(charge);
  fP[5] *= TMath::Abs(charge);

  Double_t pt = 1. / TMath::Abs(track.GetParameter()[4]) * TMath::Abs(charge);
  Double_t cs = TMath::Cos(track.GetAlpha());
  Double_t sn = TMath::Sin(track.GetAlpha());
  Double_t r = TMath::Sqrt((1. - track.GetParameter()[2]) * (1. + track.GetParameter()[2]));

  Double_t m00 = -sn;
  Double_t m10 = cs;
  Double_t m23 = -pt * (sn + track.GetParameter()[2] * cs / r);
  Double_t m43 = -pt * pt * (r * cs - track.GetParameter()[2] * sn);
  Double_t m24 = pt * (cs - track.GetParameter()[2] * sn / r);
  Double_t m44 = -pt * pt * (r * sn + track.GetParameter()[2] * cs);
  Double_t m35 = pt;
  Double_t m45 = -pt * pt * track.GetParameter()[3];

  m43 *= track.GetSign();
  m44 *= track.GetSign();
  m45 *= track.GetSign();

  const Double_t* cTr = track.GetCovariance();
  Double_t fC[21];
  fC[0] = cTr[0] * m00 * m00;
  fC[1] = cTr[0] * m00 * m10;
  fC[2] = cTr[0] * m10 * m10;
  fC[3] = cTr[1] * m00;
  fC[4] = cTr[1] * m10;
  fC[5] = cTr[2];
  fC[6] = m00 * (cTr[3] * m23 + cTr[10] * m43);
  fC[7] = m10 * (cTr[3] * m23 + cTr[10] * m43);
  fC[8] = cTr[4] * m23 + cTr[11] * m43;
  fC[9] = m23 * (cTr[5] * m23 + cTr[12] * m43) + m43 * (cTr[12] * m23 + cTr[14] * m43);
  fC[10] = m00 * (cTr[3] * m24 + cTr[10] * m44);
  fC[11] = m10 * (cTr[3] * m24 + cTr[10] * m44);
  fC[12] = cTr[4] * m24 + cTr[11] * m44;
  fC[13] = m23 * (cTr[5] * m24 + cTr[12] * m44) + m43 * (cTr[12] * m24 + cTr[14] * m44);
  fC[14] = m24 * (cTr[5] * m24 + cTr[12] * m44) + m44 * (cTr[12] * m24 + cTr[14] * m44);
  fC[15] = m00 * (cTr[6] * m35 + cTr[10] * m45);
  fC[16] = m10 * (cTr[6] * m35 + cTr[10] * m45);
  fC[17] = cTr[7] * m35 + cTr[11] * m45;
  fC[18] = m23 * (cTr[8] * m35 + cTr[12] * m45) + m43 * (cTr[13] * m35 + cTr[14] * m45);
  fC[19] = m24 * (cTr[8] * m35 + cTr[12] * m45) + m44 * (cTr[13] * m35 + cTr[14] * m45);
  fC[20] = m35 * (cTr[9] * m35 + cTr[13] * m45) + m45 * (cTr[13] * m35 + cTr[14] * m45);

  KFParticle part;
  part.Create(fP, fC, fQ, mass);

  return part;
}

/*
  Correct initialization of a KFVertex
  [copied from AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx]
*/
KFVertex AliAnalysisTaskLambda1520Lpipi::CreateKFVertex(const AliVVertex& vertex) {

  Double_t param[6];
  vertex.GetXYZ(param);

  Double_t cov[6];
  vertex.GetCovarianceMatrix(cov);

  KFPVertex kfpVtx;
  Float_t paramF[3] = {(Float_t)param[0], (Float_t)param[1], (Float_t)param[2]};
  kfpVtx.SetXYZ(paramF);
  Float_t covF[6] = {(Float_t)cov[0], (Float_t)cov[1], (Float_t)cov[2], (Float_t)cov[3], (Float_t)cov[4], (Float_t)cov[5]};
  kfpVtx.SetCovarianceMatrix(covF);

  KFVertex KFVtx(kfpVtx);

  return KFVtx;
}
