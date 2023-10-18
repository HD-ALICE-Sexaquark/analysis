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

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliAnalysisTaskLambda1520Lpipi.h"

#define N_SIGMA 3.0
#define DEBUG_MODE  // 1 : MC, 2 : Rec., 3 : V0s

class AliAnalysisTaskLambda1520Lpipi;

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
  InitPDGMasses();

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
 Set mass-at-rest of different particles
 */
void AliAnalysisTaskLambda1520Lpipi::InitPDGMasses() {

  kMassNeutron = 0.939565;
  AliInfoF("kMassNeutron = %.6f", kMassNeutron);

  kMassProton = 0.938272;
  AliInfoF("kMassProton = %.6f", kMassProton);

  kMassPion = 0.139570;
  AliInfoF("kMassPion = %.6f", kMassPion);

  kMassKaon = 0.493677;
  AliInfoF("kMassKaon = %.6f", kMassKaon);
}

/*
 Main function, called per each event.
 */
void AliAnalysisTaskLambda1520Lpipi::UserExec(Option_t*) {

  // containers
  std::set<Int_t> Indices_MCGen_FS;
  std::set<Int_t> Indices_MCGen_FS_Signal;  // FS: final-state
  std::vector<Int_t> Indices_MCRec_PID;
  std::vector<Int_t> Indices_MCRec_Signal;

  std::map<Int_t, Int_t> Daughters_IndexToVectorPos;
  std::map<Int_t, Int_t> Daughters_VectorPosToIndex;
  std::vector<AliKFParticle> Daughters;

  std::map<Int_t, std::vector<Int_t>> V0s_VectorPosToDauVectorPos;
  std::vector<AliKFParticle> V0s;

  std::vector<AliKFParticle> Lambda1520s;

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

  // initialize KFparticle
  AliKFParticle::SetField(fESD->GetMagneticField());

  // load primary vertex
  // [note: const_cast removes the const qualifier from the pointer]
  fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());
#ifdef DEBUG_MODE
  AliInfoF("REC. PRIMARY VERTEX = (%8.3f, %8.3f, %8.3f)", fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
#endif

  ProcessMCGen(Indices_MCGen_FS, Indices_MCGen_FS_Signal);

  ProcessMCRec(Indices_MCGen_FS_Signal, Indices_MCRec_Signal);

  ProcessTrueV0s_KF(Indices_MCRec_Signal,         //
                    Daughters_IndexToVectorPos,   //
                    Daughters_VectorPosToIndex,   //
                    Daughters,                    //
                    V0s_VectorPosToDauVectorPos,  //
                    V0s);

  ProcessOfficialV0s(Indices_MCRec_Signal);

  Lambda1520Finder(Daughters,                    //
                   V0s_VectorPosToDauVectorPos,  //
                   V0s,                          //
                   Lambda1520s);

  // (tree operations) after all variables have been assigned, fill tree
  // fTree->Fill();

  // stream the results the analysis of this event to the output manager
  PostData(1, fOutputListOfTrees);
  PostData(2, fOutputListOfHists);

  /* End of event */

  // clean containers
  Indices_MCGen_FS.clear();
  Indices_MCGen_FS_Signal.clear();
  Indices_MCRec_PID.clear();
  Indices_MCRec_Signal.clear();
  Daughters_IndexToVectorPos.clear();
  Daughters_VectorPosToIndex.clear();
  Daughters.clear();
  V0s_VectorPosToDauVectorPos.clear();
  V0s.clear();
  Lambda1520s.clear();
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
 - input: fESD, fMC, Indices_MCGen_FS_Signal
 - output: Indices_MCRec_Signal
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessMCRec(std::set<Int_t> Indices_MCGen_FS_Signal, std::vector<Int_t>& Indices_MCRec_Signal) {

  AliESDtrack* track;
  Int_t idx_true;

  AliMCParticle* truePart;
  Float_t n_sigma_pion;
  Float_t n_sigma_proton;

  Int_t counter = 0;  // (debug)

  // loop over MC rec. tracks in a single event
#ifdef DEBUG_MODE
  AliInfo("     IDX  trueIDX  truePID nsigmapi  nsigmaP");
#endif
  for (Int_t idx_track = 0; idx_track < fESD->GetNumberOfTracks(); idx_track++) {

    track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));
    idx_true = TMath::Abs(track->GetLabel());

    // (cut) PENDING!
    if (!Indices_MCGen_FS_Signal.count(idx_true)) {
      continue;
    }

    counter++;  // (debug)

    truePart = (AliMCParticle*)fMC->GetTrack(idx_true);
    n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

    // (cut)
    if (TMath::Abs(n_sigma_pion) > N_SIGMA && TMath::Abs(n_sigma_proton) > N_SIGMA) continue;

    if (truePart->GetMother() != -1) Indices_MCRec_Signal.push_back(idx_track);

#ifdef DEBUG_MODE
    AliInfoF("%8i %8i %8i %8.2f %8.2f", idx_track, idx_true, truePart->PdgCode(), n_sigma_pion, n_sigma_proton);
#endif
  }  // end of loop over tracks

#ifdef DEBUG_MODE
  if (counter >= 4) AliInfo("IT HAPPENED!");
#endif
}

/*
 Loop over all the V0s found by the Official Offline V0 Finder
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

  Double_t b = fESD->GetMagneticField();
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
 - input: fESD, fMC, Indices_MCRec_Signal
 - output: Daughters_IndexToVectorPos, Daughters_VectorPosToIndex, Daughters, V0s_VectorPosToDauVectorPos, V0s
 */
void AliAnalysisTaskLambda1520Lpipi::ProcessTrueV0s_KF(std::vector<Int_t> Indices_MCRec_Signal,
                                                       std::map<Int_t, Int_t>& Daughters_IndexToVectorPos,
                                                       std::map<Int_t, Int_t>& Daughters_VectorPosToIndex,
                                                       std::vector<AliKFParticle>& Daughters,
                                                       std::map<Int_t, std::vector<Int_t>>& V0s_VectorPosToDauVectorPos,
                                                       std::vector<AliKFParticle>& V0s) {
  Int_t indexDau1;
  Int_t indexDau2;

  AliESDtrack* trackDau1;
  AliESDtrack* trackDau2;

  AliMCParticle* mcPartDau1;
  AliMCParticle* mcPartDau2;

  AliMCParticle* mcMother1;
  AliMCParticle* mcMother2;

  Double_t b = fESD->GetMagneticField();
  Double_t xthis, xp;

#ifdef DEBUG_MODE
  AliInfo("     mom  mom_pdg   track");
#endif
  for (Int_t i = 0; i < (Int_t)Indices_MCRec_Signal.size() - 1; i++) {
    for (Int_t j = i + 1; j < (Int_t)Indices_MCRec_Signal.size(); j++) {

      indexDau1 = Indices_MCRec_Signal[i];
      indexDau2 = Indices_MCRec_Signal[j];

      trackDau1 = static_cast<AliESDtrack*>(fESD->GetTrack(indexDau1));
      trackDau2 = static_cast<AliESDtrack*>(fESD->GetTrack(indexDau2));

      mcPartDau1 = (AliMCParticle*)fMC->GetTrack(TMath::Abs(trackDau1->GetLabel()));
      mcPartDau2 = (AliMCParticle*)fMC->GetTrack(TMath::Abs(trackDau2->GetLabel()));

      if (mcPartDau1->GetMother() != -1) mcMother1 = (AliMCParticle*)fMC->GetTrack(mcPartDau1->GetMother());
      if (mcPartDau2->GetMother() != -1) mcMother2 = (AliMCParticle*)fMC->GetTrack(mcPartDau2->GetMother());

#ifdef DEBUG_MODE
      AliInfoF("%8i %8i %8i", mcPartDau1->GetMother(), mcMother1->PdgCode(), indexDau1);
      AliInfoF("%8.3f %8.3f %8.3f", mcPartDau1->Xv(), mcPartDau1->Yv(), mcPartDau1->Zv());

      AliInfoF("%8i %8i %8i", mcPartDau2->GetMother(), mcMother2->PdgCode(), indexDau2);
      AliInfoF("%8.3f %8.3f %8.3f", mcPartDau2->Xv(), mcPartDau2->Yv(), mcPartDau2->Zv());
#endif

      // (cut) PENDING
      // if (mcPartDau1->GetMother() != mcPartDau2->GetMother()) continue;

      // (cut)
      if (trackDau1->GetDCA(trackDau2, b, xthis, xp) > 2.) continue;

      /* Kalman Filter */

      AliKFVertex kfPrimaryVertex = CreateKFVertex(fPrimaryVertex);

      AliKFParticle kfDaughter1 = CreateKFParticle(*trackDau1, mcPartDau1->PdgCode());
      Daughters_IndexToVectorPos[indexDau1] = (Int_t)Daughters.size();
      Daughters_VectorPosToIndex[(Int_t)Daughters.size()] = indexDau1;
      Daughters.push_back(kfDaughter1);

      AliKFParticle kfDaughter2 = CreateKFParticle(*trackDau2, mcPartDau2->PdgCode());
      Daughters_IndexToVectorPos[indexDau2] = (Int_t)Daughters.size();
      Daughters_VectorPosToIndex[(Int_t)Daughters.size()] = indexDau2;
      Daughters.push_back(kfDaughter2);

      AliKFParticle kfV0(kfDaughter1, kfDaughter2);
      kfV0.TransportToDecayVertex();
      V0s_VectorPosToDauVectorPos[(Int_t)V0s.size()].push_back(Daughters_IndexToVectorPos[indexDau1]);
      V0s_VectorPosToDauVectorPos[(Int_t)V0s.size()].push_back(Daughters_IndexToVectorPos[indexDau2]);
      V0s.push_back(kfV0);

#ifdef DEBUG_MODE
      AliInfoF("%8.3f %8.3f %8.3f %8.3f", kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), kfV0.GetMass());
      AliInfoF("%i %i dca(d1,pv) = %5.3f, dca(d2,pv) = %5.3f, dca(d1,V0) = %5.3f, dca(d2,V0) = %5.3f",  //
               indexDau1,                                                                               //
               indexDau2,                                                                               //
               kfDaughter1.GetDistanceFromVertex(kfPrimaryVertex),                                      //
               kfDaughter2.GetDistanceFromVertex(kfPrimaryVertex),                                      //
               (Float_t)kfDaughter1.GetDistanceFromVertex(kfV0),                                        //
               (Float_t)kfDaughter2.GetDistanceFromVertex(kfV0));
#endif
    }
  }
}

/*
 Find all Lambda1520 candidates. It will pair all available V0s into V0A and V0B, which are split into PRON1 and PRON2 for V0A, and PRON3
 and PRON4 for V0B. Then, two Lambda1520 candidates will be formed: V0A+PRON3+PRON4 and V0B+PRON1+PRON2. The one with a better Chi2 will be
 stored.
 - input: Daughters, V0s_VectorPosToDauVectorPos, V0s
 - output: Lambda1520s
 */
void AliAnalysisTaskLambda1520Lpipi::Lambda1520Finder(std::vector<AliKFParticle> Daughters,
                                                      std::map<Int_t, std::vector<Int_t>> V0s_VectorPosToDauVectorPos,
                                                      std::vector<AliKFParticle> V0s,  //
                                                      std::vector<AliKFParticle>& Lambda1520s) {

#ifdef DEBUG_MODE
  AliInfo("               Mass       Vx       Vy       Vz     Chi2");
#endif
  for (Int_t pos_v0a = 0; pos_v0a < (Int_t)V0s.size() - 1; pos_v0a++) {
    for (Int_t pos_v0b = pos_v0a + 1; pos_v0b < (Int_t)V0s.size(); pos_v0b++) {

      AliKFVertex kfPrimaryVertex_I(*fPrimaryVertex);
      AliKFVertex kfPrimaryVertex_II(*fPrimaryVertex);

      AliKFParticle kfV0A = V0s[pos_v0a];
      AliKFParticle kfProng1 = Daughters[V0s_VectorPosToDauVectorPos[pos_v0a][0]];
      AliKFParticle kfProng2 = Daughters[V0s_VectorPosToDauVectorPos[pos_v0a][1]];

      AliKFParticle kfV0B = V0s[pos_v0b];
      AliKFParticle kfProng3 = Daughters[V0s_VectorPosToDauVectorPos[pos_v0b][0]];
      AliKFParticle kfProng4 = Daughters[V0s_VectorPosToDauVectorPos[pos_v0b][1]];

      AliKFParticle kfLambda1520_I(kfV0A, kfProng3, kfProng4);
      kfPrimaryVertex_I += kfLambda1520_I;
      kfLambda1520_I.SetProductionVertex(kfPrimaryVertex_I);
      kfLambda1520_I.SetNoDecayLength();

      AliKFParticle kfLambda1520_II(kfV0B, kfProng1, kfProng2);
      kfPrimaryVertex_II += kfLambda1520_II;
      kfLambda1520_II.SetProductionVertex(kfPrimaryVertex_II);
      kfLambda1520_II.SetNoDecayLength();

      AliKFParticle kfLambda1520_I_NoVtxSet(kfV0A, kfProng3, kfProng4);
      // kfLambda1520_I_NoVtxSet.SetNoDecayLength();

      AliKFParticle kfLambda1520_II_NoVtxSet(kfV0B, kfProng1, kfProng2);
      // kfLambda1520_II_NoVtxSet.SetNoDecayLength();

#ifdef DEBUG_MODE
      /*
      AliInfoF(" I -     - %8.3f %8.3f %8.3f %8.3f %8.3f", kfLambda1520_I.GetMass(), kfLambda1520_I.GetX(), kfLambda1520_I.GetY(),
               kfLambda1520_I.GetZ(), kfLambda1520_I.GetChi2());
      AliInfoF("II -     - %8.3f %8.3f %8.3f %8.3f %8.3f", kfLambda1520_II.GetMass(), kfLambda1520_II.GetX(), kfLambda1520_II.GetY(),
               kfLambda1520_II.GetZ(), kfLambda1520_II.GetChi2());
      AliInfoF(" I - NVS - %8.3f %8.3f %8.3f %8.3f %8.3f", kfLambda1520_I_NoVtxSet.GetMass(), kfLambda1520_I_NoVtxSet.GetX(),
               kfLambda1520_I_NoVtxSet.GetY(), kfLambda1520_I_NoVtxSet.GetZ(), kfLambda1520_I_NoVtxSet.GetChi2());
      AliInfoF("II - NVS - %8.3f %8.3f %8.3f %8.3f %8.3f", kfLambda1520_II_NoVtxSet.GetMass(), kfLambda1520_II_NoVtxSet.GetX(),
               kfLambda1520_II_NoVtxSet.GetY(), kfLambda1520_II_NoVtxSet.GetZ(), kfLambda1520_II_NoVtxSet.GetChi2());
      */
#endif

      AliKFParticle kfLambda1520 =
          kfLambda1520_I_NoVtxSet.GetChi2() < kfLambda1520_II_NoVtxSet.GetChi2() ? kfLambda1520_I_NoVtxSet : kfLambda1520_II_NoVtxSet;

#ifdef DEBUG_MODE
      AliInfoF("   -     - %8.3f %8.3f %8.3f %8.3f %8.3f", kfLambda1520.GetMass(), kfLambda1520.GetX(), kfLambda1520.GetY(),
               kfLambda1520.GetZ(), kfLambda1520.GetChi2());
#endif

      Lambda1520s.push_back(kfLambda1520);
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

  return (Float_t)track->GetD(PV[0], PV[1], fESD->GetMagneticField());
}

/*
 Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 This function stores the point of closest approach, and returns its distance to the PV
 */
Double_t AliAnalysisTaskLambda1520Lpipi::Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                                Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                                Double_t PV_X, Double_t PV_Y, Double_t PV_Z) {
  TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
  TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
  TVector3 PrimVertex(PV_X, PV_Y, PV_Z);

  TVector3 CrossProduct = (PrimVertex - V0Vertex).Cross(V0Momentum);

  return CrossProduct.Mag() / V0Momentum.Mag();
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
