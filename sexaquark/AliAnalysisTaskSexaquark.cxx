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
#define DEBUG_MODE 1  // 1 : MC, 2 : Rec., 3 : V0s

class AliAnalysisTaskSexaquark;

/*
 Create output objects
 This function is called ONCE at the start of your analysis (RUNTIME)
*/
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {

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

    /*** List of Trees ***/

    fOutputListOfTrees = new TList();
    fOutputListOfTrees->SetOwner(kTRUE);

    fTree = new TTree("Events", "Tree of Events");
    SetBranches();
    fOutputListOfTrees->Add(fTree);

    PostData(1, fOutputListOfTrees);
}

/*
 Set mass-at-rest of different particles
*/
void AliAnalysisTaskSexaquark::InitPDGMasses() {
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
 This function is called per each event.
*/
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

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

    // load primary vertex
    // [note: const_cast removes the const qualifier from the pointer]
    fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());

    // process MC tree
    MCGen_FindRelevant();

    // process reconstructed tracks
    MCRec_FindRelevant();

    // add non-relevant MC particles, that were rec. as relevant
    MCGen_AddNonRelevant();

    // label duplicated tracks
    MCRec_FindDuplicates();

    // V0 finder
    if (fSourceOfV0s == "official") {
        V0_ProcessESD();
    }
    if (fSourceOfV0s == "true") {
        V0_ProcessTrue();
    }
    if (fSourceOfV0s == "custom") {
        V0_ProcessCustom();
    }

#if DEBUG_MODE == 1
    AliInfo("MONTE CARLO TREE");
    AliInfo("");
    AliInfo("         EVT    IDX    PID SIGNAL MOTHER  M_PID GM_PID FIRDAU LASDAU");
#endif
    fEvent.N_MCGen = (Int_t)fVec_Idx_MCGen.size();
    for (Int_t ii = 0; ii < fEvent.N_MCGen; ii++) {
        MCGen_PushBack(ii);
    }

#if DEBUG_MODE > 0
    AliInfo("RECONSTRUCTED TRACKS");
    AliInfo("");
    AliInfo("         EVT    IDX EVTTRU IDXTRU DUPLI?");
#endif
    fEvent.N_MCRec = (Int_t)fVec_Idx_MCRec.size();
    for (Int_t ii = 0; ii < fEvent.N_MCRec; ii++) {
        MCRec_PushBack(ii);
    }

#if DEBUG_MODE == 2
    AliInfo("FOUND V0s");
    AliInfo("");
    AliInfo("   EVTPOS EVTNEG IDXPOS IDXNEG");
#endif
    fEvent.N_V0s = (Int_t)fVec_V0s.size();
    for (Int_t ii = 0; ii < fEvent.N_V0s; ii++) {
        V0_PushBack(ii);
    }

    // (tree operations) after all variables have been assigned, fill tree
    fTree->Fill();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfTrees);

    /*** End of Event: Clear Structures ***/

    fMap_Evt_MCGen.clear();
    fVec_Idx_MCGen.clear();

    fMap_Evt_MCRec.clear();
    fVec_Idx_MCRec.clear();

    fVec_MCRec_IsDuplicate.clear();
    fVec_MCRec_IsSimilar.clear();

    fVec_Idx_MCGen_NonRel.clear();
    fMap_Duplicates.clear();

    fVec_V0s.clear();

    ClearEvent();
}

/*
 Search for possible contradictions in the input options
*/
void AliAnalysisTaskSexaquark::CheckForInputErrors() {
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

/*
 Search for the relevant MC particles in a single event
*/
void AliAnalysisTaskSexaquark::MCGen_FindRelevant() {

    AliMCParticle* mcPart;
    AliMCParticle* mcMother;
    AliMCParticle* mcGrandMother;

    Bool_t isGenPIDRelevant;

    Int_t evt_mother;
    Int_t mother_pid;
    Int_t grandmother_pid;

    // loop over MC gen. particles in a single event
    for (Int_t idx_mc = 0; idx_mc < fMC->GetNumberOfTracks(); idx_mc++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

        isGenPIDRelevant = TMath::Abs(mcPart->PdgCode()) == 211 ||  // pos. and neg. pion
                           mcPart->PdgCode() == 310 ||              // neutral kaon short
                           mcPart->PdgCode() == -3122 ||            // anti-lambda
                           mcPart->PdgCode() == -2212;              // anti-proton

        if (!isGenPIDRelevant) {
            continue;
        }

        fMap_Evt_MCGen[idx_mc] = (Int_t)fVec_Idx_MCGen.size();
        fVec_Idx_MCGen.push_back(idx_mc);
    }  // end of loop over MC particles
}

/*
 (First) Loop over the reconstructed tracks in a single event
 Store tracks into the (duplicated tracks) map
 IMPORTANT: don't trust track->GetID()
*/
void AliAnalysisTaskSexaquark::MCRec_FindRelevant() {

    AliESDtrack* track;

    Float_t n_sigma_pion;
    Float_t n_sigma_proton;
    Bool_t isRecPIDRelevant;

    // (non.rel.->rel.)
    Int_t idx_true;
    AliMCParticle* truePart;
    Bool_t isGenPIDRelevant;

    // loop over MC rec. tracks in a single event
    for (Int_t idx_track = 0; idx_track < fESD->GetNumberOfTracks(); idx_track++) {

        track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

        n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);

        isRecPIDRelevant = (TMath::Abs(n_sigma_pion) < N_SIGMA) ||                         // pos. and neg. pion
                           (TMath::Abs(n_sigma_proton) < N_SIGMA && track->Charge() < 0);  // anti-proton

        idx_true = TMath::Abs(track->GetLabel());
        truePart = (AliMCParticle*)fMC->GetTrack(idx_true);

        isGenPIDRelevant = TMath::Abs(truePart->PdgCode()) == 211 ||  // pos. and neg. pion
                           truePart->PdgCode() == -2212;              // anti-proton

        // consider the cases:
        // - rel.     -> rel.
        // - rel.     -> non-rel. (make sure to remove these by imposing PID again!)
        // - non-rel. -> rel.     (for contamination studies)
        if (isGenPIDRelevant || isRecPIDRelevant) {

            fMap_Evt_MCRec[idx_track] = (Int_t)fVec_Idx_MCRec.size();
            fVec_Idx_MCRec.push_back(idx_track);

            // (non.rel.->rel.)
            // for when the true particle is non-relevant, but reconstructed as relevant
            if (!isGenPIDRelevant) fVec_Idx_MCGen_NonRel.push_back(idx_true);

            // (duplicated tracks) add them to the duplicates map
            fMap_Duplicates[idx_true].push_back(idx_track);
        }
    }  // end of loop over tracks
}

/*
 Calculate transverse impact parameter w.r.t. Primary Vertex
*/
Float_t AliAnalysisTaskSexaquark::MCRec_GetImpactParameter(AliESDtrack* track) {
    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

    return (Float_t)track->GetD(PV[0], PV[1], fESD->GetMagneticField());
}

/*
 Get the evt of the mother particle, the PID of the mother particle, and the PID of the grandmother particle
*/
void AliAnalysisTaskSexaquark::MCGen_GetMotherInfo(AliMCParticle* mcPart,  //
                                                   Int_t& evt_mother, Int_t& mother_pid, Int_t& grandmother_pid) {

    AliMCParticle* mcMother;
    AliMCParticle* mcGrandMother;

    Int_t idx_mother;

    evt_mother = -1;  // it can be:
    //                // -2 : mother is not relevant <-> mother outside the vector
    //                // -1 : particle has no mother <-> particle is primary
    //                //  > 0 : particle has a relevant mother
    mother_pid = 0;
    grandmother_pid = 0;

    // is it a primary?
    if (mcPart->GetMother() != -1) {
        // if not, get mother
        idx_mother = mcPart->GetMother();
        mcMother = (AliMCParticle*)fMC->GetTrack(idx_mother);
        mother_pid = mcMother->PdgCode();
        // but, what if the mother is non-relevant?
        // [note: count() in a std::map, returns the number of elements matching specific key]
        evt_mother = fMap_Evt_MCGen.count(idx_mother) != 0 ? fMap_Evt_MCGen[idx_mother] : -2;
        // is mother a primary?
        if (mcMother->GetMother() != -1) {
            // if not, get grandmother
            mcGrandMother = (AliMCParticle*)fMC->GetTrack(mcMother->GetMother());
            grandmother_pid = mcGrandMother->PdgCode();
        }
    }  // closure
}

/*
 Add the non-relevant MC particles that were reconstructed as relevant
*/
void AliAnalysisTaskSexaquark::MCGen_AddNonRelevant() {

    AliMCParticle* mcPart;

    Int_t evt_mother;
    Int_t mother_pid;
    Int_t grandmother_pid;

    for (Int_t idx_mc : fVec_Idx_MCGen_NonRel) {

        mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

        fMap_Evt_MCGen[idx_mc] = (Int_t)fVec_Idx_MCGen.size();
        fVec_Idx_MCGen.push_back(idx_mc);
        s
    }
}

/*
 Get number of daughters, and the evt of the first and last daughter
*/
void AliAnalysisTaskSexaquark::MCGen_GetDaughtersInfo(AliMCParticle* mcPart,  //
                                                      Int_t& n_daughters, Int_t& evt_first_dau, Int_t& evt_last_dau) {

    n_daughters = mcPart->GetNDaughters();
    evt_first_dau = -1;
    evt_last_dau = -1;  // it can be:
    //                 // -2 : daughter is not relevant
    //                 // -1 : particle has no daughter
    //                 //  > 0 : particle has a relevant daughter

    if (!n_daughters) {
        return;
    }

    Int_t idx_first_dau = mcPart->GetDaughterFirst();
    Int_t idx_last_dau = mcPart->GetDaughterLast();

    evt_first_dau = fMap_Evt_MCGen.count(idx_first_dau) != 0 ? fMap_Evt_MCGen[idx_first_dau] : -2;
    evt_last_dau = fMap_Evt_MCGen.count(idx_last_dau) != 0 ? fMap_Evt_MCGen[idx_last_dau] : -2;
}

/*
 Determine if a certain MC particle comes from a generated sexaquark-nucleon interaction
*/
Bool_t AliAnalysisTaskSexaquark::MCGen_IsSignal(AliMCParticle* mcPart) {

    // check if particle is a primary
    if (mcPart->GetMother() == -1) {
        return mcPart->MCStatusCode() == 6;  // comes from a sexaquark?
    }

    // if it's not primary, it must have a mother
    AliMCParticle* mcMother = (AliMCParticle*)fMC->GetTrack(mcPart->GetMother());

    // check if mother is a primary
    if (mcMother->GetMother() == -1) {
        return mcMother->MCStatusCode() == 6;  // comes from a sexaquark?
    }

    // default values
    return kFALSE;
}

/*
 Loop over the (duplicated tracks) map
 If a MC particle has been reconstructed more than once
 Tag the best particle as "non-duplicate", and the rest as "duplicates"
*/
void AliAnalysisTaskSexaquark::MCRec_FindDuplicates() {

    std::vector<Int_t> vec_indices_tracks;

    Int_t idx_mc;
    AliMCParticle* mcTrue;
    TVector3 p_true;

    AliESDtrack* track;
    TVector3 p_track;

    Int_t idx_best;
    Double_t min_angle;
    Double_t delta_angle;

    /*
    Float_t n_sigma_pion;
    Float_t n_sigma_proton;
    Bool_t isRecPIDRelevant;
    */

    // resize the duplicates vector, and fill with ones
    fVec_MCRec_IsDuplicate.resize(fVec_Idx_MCRec.size(), kTRUE);

    for (auto it = fMap_Duplicates.begin(); it != fMap_Duplicates.end(); it++) {

        idx_mc = it->first;
        vec_indices_tracks = it->second;

        // load true particle
        mcTrue = (AliMCParticle*)fMC->GetTrack(idx_mc);
        p_true.SetXYZ(mcTrue->Px(), mcTrue->Py(), mcTrue->Pz());

        // reset best track and min angle
        idx_best = -1;
        min_angle = TMath::Pi();

        // loop over tracks, look for best candidate
        for (Int_t idx_track : vec_indices_tracks) {

            // load track
            track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));
            p_track.SetXYZ(track->Px(), track->Py(), track->Pz());

            delta_angle = p_true.Angle(p_track);

            if (delta_angle < min_angle) {
                min_angle = delta_angle;
                idx_best = idx_track;
                // PENDING: add PID to this
                /*
                n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
                n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
                isRecPIDRelevant = (TMath::Abs(n_sigma_pion) < N_SIGMA) ||                         // positive and negative pion
                                   (TMath::Abs(n_sigma_proton) < N_SIGMA && track->Charge() < 0);  // anti-proton
                if (isRecPIDRelevant) {
                  idx_best = idx_track;
                }
                */
            }
        }  // end of loop over tracks

        // assign zero to the best track
        fVec_MCRec_IsDuplicate[fMap_Evt_MCRec[idx_best]] = kFALSE;
    }  // end of loop over (duplicated tracks) map
}

/*
 PENDING!
*/
void AliAnalysisTaskSexaquark::MCRec_FindSimilarTracks() {

    Int_t idx_track_a;
    AliESDtrack* track_a;

    Int_t idx_track_b;
    AliESDtrack* track_b;

    Int_t idx_best;

    // resize the similars vector, and fill with ones
    fVec_MCRec_IsSimilar.resize(fVec_Idx_MCRec.size(), kFALSE);

    for (Int_t evt_track_a = 0; evt_track_a < (Int_t)fVec_Idx_MCRec.size() - 1; evt_track_a++) {
        for (Int_t evt_track_b = evt_track_a + 1; evt_track_b < (Int_t)fVec_Idx_MCRec.size(); evt_track_b++) {

            idx_track_a = fVec_Idx_MCRec[evt_track_a];
            idx_track_b = fVec_Idx_MCRec[evt_track_b];

            track_a = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track_a));
            track_b = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track_b));

            if (TMath::Abs(1. - track_a->Px() / track_b->Px()) < 0.05 &&  //
                TMath::Abs(1. - track_a->Py() / track_b->Py()) < 0.05 &&  //
                TMath::Abs(1. - track_a->Pz() / track_b->Pz()) < 0.05) {
                // hola
            }
        }

        // assign zero to the best track
        fVec_MCRec_IsSimilar[fMap_Evt_MCRec[idx_best]] = kTRUE;
    }
}

/*
 Loop over all the V0s found by the Official Offline V0 Finder
*/
void AliAnalysisTaskSexaquark::V0_ProcessESD() {

    // define variables & objects
    AliESDv0* V0;

    V0_tt this_V0;

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
        // => the particle was gen. non-relevant & rec. non-relevant
        if (fMap_Evt_MCRec.count(idx_pos) == 0 || fMap_Evt_MCRec.count(idx_neg) == 0) {
            // [note: count() in std::map, returns the number of elements matching specific key]
            continue;
        }

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_pos));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_neg));

        nsigma_pos_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion));
        nsigma_neg_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kPion));
        nsigma_antiproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kProton));

        V0->GetXYZ(V0_X, V0_Y, V0_Z);
        V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
        V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

        // calculate energy
        // for positive daughter, assuming it's a positive pion
        pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + kMassPion * kMassPion);
        // for negative daughter, assuming it's K0 -> pi+ pi-
        neg_energy_asK0 = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassPion * kMassPion);
        // for negative daughter, assuming it's anti-lambda -> pi+ anti-proton
        neg_energy_asAL = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassProton * kMassProton);

        // get true particles
        idx_pos_true = TMath::Abs(pos_track->GetLabel());
        mcPosTrue = (AliMCParticle*)fMC->GetTrack(idx_pos_true);
        idx_neg_true = TMath::Abs(neg_track->GetLabel());
        mcNegTrue = (AliMCParticle*)fMC->GetTrack(idx_neg_true);

        this_V0.Idx_Pos = fMap_Evt_MCRec[idx_pos];
        this_V0.Idx_Neg = fMap_Evt_MCRec[idx_neg];
        this_V0.Px = V0->Px();
        this_V0.Py = V0->Py();
        this_V0.Pz = V0->Pz();
        this_V0.X = V0_X;
        this_V0.Y = V0_Y;
        this_V0.Z = V0_Z;
        this_V0.Pos_Px = P_Px;
        this_V0.Pos_Py = P_Py;
        this_V0.Pos_Pz = P_Pz;
        this_V0.Neg_Px = N_Px;
        this_V0.Neg_Py = N_Py;
        this_V0.Neg_Pz = N_Pz;
        this_V0.isSignal = MCGen_IsSignal(mcPosTrue) && MCGen_IsSignal(mcNegTrue);
        this_V0.E_asK0 = pos_energy + neg_energy_asK0;
        this_V0.E_asAL = pos_energy + neg_energy_asAL;
        this_V0.couldBeK0 = TMath::Abs(nsigma_pos_pion) < N_SIGMA && TMath::Abs(nsigma_neg_pion) < N_SIGMA;
        this_V0.couldBeAL = TMath::Abs(nsigma_pos_pion) < N_SIGMA && TMath::Abs(nsigma_antiproton) < N_SIGMA;
        this_V0.Chi2 = V0->GetChi2V0();
        this_V0.DCA_Daughters = V0->GetDcaV0Daughters();
        this_V0.IP_wrtPV = V0->GetD(PV[0], PV[1], PV[2]);
        this_V0.CPA_wrtPV = V0->GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
        this_V0.ArmAlpha = V0->AlphaV0();
        this_V0.ArmPt = V0->PtArmV0();
        this_V0.DecayLength = TMath::Sqrt(TMath::Power(V0_X - PV[0], 2) + TMath::Power(V0_Y - PV[1], 2) + TMath::Power(V0_Z - PV[2], 2));
        this_V0.DCA_wrtPV = Calculate_LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]);

        fVec_V0s.push_back(this_V0);
    }  // end of loop over V0s
}

/*
 Find all V0s
 we will store all lambdas and kaons with no ambiguity
 - no dependence on PID, nor ambiguity on measurements
 - only dependence: if they were reconstructed or not
*/
void AliAnalysisTaskSexaquark::V0_ProcessTrue() {

    V0_tt this_V0;

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
        if (fMap_Duplicates.count(mc_mother->GetDaughterFirst()) == 0 || fMap_Duplicates.count(mc_mother->GetDaughterLast()) == 0) {
            // [note: count() in std::map, returns the number of elements matching specific key]
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
        for (Int_t idx_pos_track : fMap_Duplicates[idx_pos_mc]) {
            for (Int_t idx_neg_track : fMap_Duplicates[idx_neg_mc]) {

                // load MC particle depending if positive/negative
                mc_pos = (AliMCParticle*)fMC->GetTrack(idx_pos_mc);
                mc_neg = (AliMCParticle*)fMC->GetTrack(idx_neg_mc);

                // calculate energy (assuming positive pion)
                pos_energy = TMath::Sqrt(mc_pos->Px() * mc_pos->Px() + mc_pos->Py() * mc_pos->Py() + mc_pos->Pz() * mc_pos->Pz() +
                                         kMassPion * kMassPion);
                // assuming K0 -> pi+ pi-
                neg_energy_asK0 = TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() +
                                              kMassPion * kMassPion);
                // assuming anti-lambda -> pi+ anti-proton
                neg_energy_asAL = TMath::Sqrt(mc_neg->Px() * mc_neg->Px() + mc_neg->Py() * mc_neg->Py() + mc_neg->Pz() * mc_neg->Pz() +
                                              kMassProton * kMassProton);

                // (tree operations)
                // fill common format and push_back
                this_V0.Idx_Pos = fMap_Evt_MCRec[idx_pos_track];
                this_V0.Idx_Neg = fMap_Evt_MCRec[idx_neg_track];
                this_V0.Px = mc_mother->Px();
                this_V0.Py = mc_mother->Py();
                this_V0.Pz = mc_mother->Pz();
                this_V0.X = mc_pos->Xv();  // use the origin of the daughters
                this_V0.Y = mc_pos->Yv();
                this_V0.Z = mc_pos->Zv();
                this_V0.Pos_Px = mc_pos->Px();
                this_V0.Pos_Py = mc_pos->Py();
                this_V0.Pos_Pz = mc_pos->Pz();
                this_V0.Neg_Px = mc_neg->Px();
                this_V0.Neg_Py = mc_neg->Py();
                this_V0.Neg_Pz = mc_neg->Pz();
                this_V0.isSignal = MCGen_IsSignal(mc_mother);
                this_V0.E_asK0 = pos_energy + neg_energy_asK0;
                this_V0.E_asAL = pos_energy + neg_energy_asAL;
                this_V0.couldBeK0 = mc_mother->PdgCode() == 310;
                this_V0.couldBeAL = mc_mother->PdgCode() == -3122;
                this_V0.Chi2 = 1.0;          // no KF was used for this
                this_V0.DCA_Daughters = 0.;  // PENDING! too complicated to calculate using true information...
                this_V0.IP_wrtPV = 0.;       // PENDING! still not sure about the definition of this one...
                this_V0.CPA_wrtPV = Calculate_CPA(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),      //
                                                  mc_pos->Xv(), mc_pos->Yv(), mc_pos->Zv(),               //
                                                  PV[0], PV[1], PV[2]);                                   //
                this_V0.ArmAlpha = Calculate_ArmAlpha(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),  //
                                                      mc_pos->Px(), mc_pos->Py(), mc_pos->Pz(),           //
                                                      mc_neg->Px(), mc_neg->Py(), mc_neg->Pz());          //
                this_V0.ArmPt = Calculate_ArmPt(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),        //
                                                mc_neg->Px(), mc_neg->Py(), mc_neg->Pz());                //
                this_V0.DecayLength = TMath::Sqrt(TMath::Power(mc_pos->Xv() - PV[0], 2) + TMath::Power(mc_pos->Yv() - PV[1], 2) +
                                                  TMath::Power(mc_pos->Zv() - PV[2], 2));
                this_V0.DCA_wrtPV = Calculate_LinePointDCA(mc_mother->Px(), mc_mother->Py(), mc_mother->Pz(),  //
                                                           mc_pos->Xv(), mc_pos->Yv(), mc_pos->Zv(),           //
                                                           PV[0], PV[1], PV[2]);

                fVec_V0s.push_back(this_V0);
            }  // end of loop over negative tracks
        }      // end of loop over positive tracks

    }  // end of loop over MC gen. mothers
}

/*
 Custom Offline V0 Finder
 [copied and adapted from AliRoot/STEER/ESD/AliV0vertexer.cxx]
*/
void AliAnalysisTaskSexaquark::V0_ProcessCustom() {

    AliESDtrack* track;
    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    V0_tt this_V0;

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
    Double_t dca_to_pv;

    Int_t idx_pos_true;
    Int_t idx_neg_true;
    AliMCParticle* mcPosTrue;
    AliMCParticle* mcNegTrue;
    Int_t idx_pos_mother;
    Int_t idx_neg_mother;

    // define cuts
    Double_t COV0F_DNmin = 0.1;            // (d=0.1) min imp parameter for the negative daughter (depends on PV!)
    Double_t COV0F_DPmin = 0.1;            // (d=0.1) min imp parameter for the positive daughter (depends on PV!)
    Double_t COV0F_Chi2max = 33.;          // (d=33.) max chi2
    Double_t COV0F_DCAmax = 0.6;           // (d=1.) max DCA between the daughter tracks
    Double_t COV0F_CPAmin = 0.4;           // (d=0.998) max cosine of V0's pointing angle (depends on PV!)
    Double_t COV0F_Rmin = 50.;             // (d=0.9) min radius of the V0 radius
    Double_t COV0F_Rmax = 200.;            // (d=100) max radius of the V0 radius
    Double_t COV0F_Etamax_V0 = 1.0;        // (d=1.) max eta of V0
    Double_t COV0F_DCAmin_V0_wrtPV = 5.;   // min DCA between V0 and PV
    Double_t COV0F_Etamax_NegDau = 1.5;    // max eta of negative daughter
    Double_t COV0F_Etamax_PosDau = 1.5;    // max eta of positive daughter
    Int_t COV0F_NClustersTPCmin_Dau = 50;  // minimun number of clusters in the TPC
    Double_t COV0F_Ptmin = 1.0;            // min Pt of V0

    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

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
    for (Int_t idx_track = 0; idx_track < nentr; idx_track++) {

        track = fESD->GetTrack(idx_track);

        // (cut) if particle is not in the indexer
        // => the particle was gen. non-relevant & rec. non-relevant
        if (fMap_Evt_MCRec.count(idx_track) == 0) {
            // [note: count() returns the number of elements matching a specific key]
            continue;
        }

        // (cut) on status
        ULong64_t status = track->GetStatus();
        if ((status & AliESDtrack::kTPCrefit) == 0) {
            continue;
        }

        // (cut) min. number of TPC clusters
        if ((Int_t)track->GetTPCncls() < COV0F_NClustersTPCmin_Dau) {
            continue;
        }

        // cut on transverse impact parameter
        Double_t d = 0.0;
        if (track->GetSign() < 0.) {
            d = track->GetD(PV[0], PV[1], b);
            if (TMath::Abs(d) < COV0F_DNmin) {
                continue;
            }
            if (TMath::Abs(d) > COV0F_Rmax) {
                continue;
            }
            neg[nneg++] = idx_track;
        } else {
            d = track->GetD(PV[0], PV[1], b);
            if (TMath::Abs(d) < COV0F_DPmin) {
                continue;
            }
            if (TMath::Abs(d) > COV0F_Rmax) {
                continue;
            }
            pos[npos++] = idx_track;
        }
    }  // end of first loop

    // second loop: test all neg-pos pairs
    for (Int_t i = 0; i < nneg; i++) {

        Int_t nidx = neg[i];
        neg_track = fESD->GetTrack(nidx);

        for (Int_t k = 0; k < npos; k++) {

            Int_t pidx = pos[k];
            pos_track = fESD->GetTrack(pidx);

            // (cut) -1. < eta(pi+) < 1.
            if (TMath::Abs(pos_track->Eta()) > COV0F_Etamax_PosDau) {
                continue;
            }

            // (cut) -1.5 < eta(pi-,anti-p) < 1.5
            if (TMath::Abs(neg_track->Eta()) > COV0F_Etamax_NegDau) {
                continue;
            }

            AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;

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
            if (cpa < COV0F_CPAmin) {
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
            if (TMath::Abs(vertex.Eta()) > COV0F_Etamax_V0) {
                continue;
            }

            // vertex.ChangeMassHypothesis(kK0Short);   // borquez: this was here...

            // get momentum of negative and positive components at V0 vertex
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            // (cut) on Pt
            if (TMath::Sqrt(vertex.Px() * vertex.Px() + vertex.Py() * vertex.Py()) < COV0F_Ptmin) {
                continue;
            }

            // get NSigmas for PID
            nsigma_pos_pion = fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion);
            nsigma_neg_pion = fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kPion);
            nsigma_antiproton = fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kProton);

            // (cut) on PID
            // positive daughter must be a positive pion
            if (TMath::Abs(nsigma_pos_pion) > N_SIGMA) {
                continue;
            }
            // negative daughter must be a negative pion or an anti-lambda
            if (TMath::Abs(nsigma_neg_pion) > N_SIGMA && TMath::Abs(nsigma_antiproton) > N_SIGMA) {
                continue;
            }

            // (cut) on sign
            if (pos_track->GetSign() == neg_track->GetSign()) {
                continue;
            }

            // calculate energy
            // for positive daughter, assume it's a positive pion
            pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + kMassPion * kMassPion);
            // for negative daughter, assuming it's K0 -> pi+ pi-
            neg_energy_asK0 = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassPion * kMassPion);
            // for negative daughter, assuming it's anti-lambda -> pi+ anti-proton
            neg_energy_asAL = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + kMassProton * kMassProton);

            // get true particles
            idx_pos_true = TMath::Abs(pos_track->GetLabel());
            mcPosTrue = (AliMCParticle*)fMC->GetTrack(idx_pos_true);
            idx_neg_true = TMath::Abs(neg_track->GetLabel());
            mcNegTrue = (AliMCParticle*)fMC->GetTrack(idx_neg_true);

            // get mother of true particles
            idx_pos_mother = mcPosTrue->GetMother();
            idx_neg_mother = mcNegTrue->GetMother();

            this_V0.Idx_Pos = fMap_Evt_MCRec[vertex.GetPindex()];
            this_V0.Idx_Neg = fMap_Evt_MCRec[vertex.GetNindex()];
            this_V0.Px = vertex.Px();
            this_V0.Py = vertex.Py();
            this_V0.Pz = vertex.Pz();
            this_V0.X = V0_X;
            this_V0.Y = V0_Y;
            this_V0.Z = V0_Z;
            this_V0.Pos_Px = P_Px;
            this_V0.Pos_Py = P_Py;
            this_V0.Pos_Pz = P_Pz;
            this_V0.Neg_Px = N_Px;
            this_V0.Neg_Py = N_Py;
            this_V0.Neg_Pz = N_Pz;
            this_V0.E_asK0 = pos_energy + neg_energy_asK0;
            this_V0.E_asAL = pos_energy + neg_energy_asAL;
            this_V0.couldBeK0 = TMath::Abs(nsigma_pos_pion) < N_SIGMA && TMath::Abs(nsigma_neg_pion) < N_SIGMA;
            this_V0.couldBeAL = TMath::Abs(nsigma_pos_pion) < N_SIGMA && TMath::Abs(nsigma_antiproton) < N_SIGMA;
            this_V0.Chi2 = vertex.GetChi2V0();
            this_V0.DCA_Daughters = vertex.GetDcaV0Daughters();
            this_V0.IP_wrtPV = vertex.GetD(PV[0], PV[1], PV[2]);
            this_V0.CPA_wrtPV = vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
            this_V0.ArmAlpha = vertex.AlphaV0();
            this_V0.ArmPt = vertex.PtArmV0();
            this_V0.DecayLength =
                TMath::Sqrt(TMath::Power(V0_X - PV[0], 2) + TMath::Power(V0_Y - PV[1], 2) + TMath::Power(V0_Z - PV[2], 2));
            this_V0.DCA_wrtPV = dca_to_pv;

            // important: to be a signal V0, the MC particles must come from the same V0 as well
            this_V0.isSignal = MCGen_IsSignal(mcNegTrue) && MCGen_IsSignal(mcPosTrue) && idx_neg_mother == idx_pos_mother;

            fVec_V0s.push_back(this_V0);
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
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

/*** FUNCTIONS FOR THE CUSTOM OFFLINE V0 FINDER ***/

/*
 This function pre-optimizes a two-track pair in the XY plane
 and provides two X values for the tracks if successful
 [copied from AliRoot/STEER/ESD/AliV0vertexer.cxx]
*/
Bool_t AliAnalysisTaskSexaquark::Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn,
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
                lPreprocessDCAxy =
                    TMath::Sqrt(TMath::Power(lCase1NegR[0] - lCase1PosR[0], 2) + TMath::Power(lCase1NegR[1] - lCase1PosR[1], 2) +
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
                lCase2aDCA =
                    TMath::Sqrt(TMath::Power(lCase2aNegR[0] - lCase2aPosR[0], 2) + TMath::Power(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                                TMath::Power(lCase2aNegR[2] - lCase2aPosR[2], 2));
            }
        }

        // Case 2b
        if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX) {
            Double_t lCase2bNegR[3] = {0.}, lCase2bPosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg[1], b, lCase2bNegR) && pt->GetXYZAt(xThisPos[1], b, lCase2bPosR)) {
                lCase2bDCA =
                    TMath::Sqrt(TMath::Power(lCase2bNegR[0] - lCase2bPosR[0], 2) + TMath::Power(lCase2bNegR[1] - lCase2bPosR[1], 2) +
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
 (based on AliRoot/STEER/ESD/AliV0vertexer.cxx)
*/
void AliAnalysisTaskSexaquark::GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b) {

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

/*
 Define fTree's branches, and load them into fEvent (which is an Event_tt object)
*/
void AliAnalysisTaskSexaquark::SetBranches() {
    /* MC particles */
    fTree->Branch("N_MCGen", &fEvent.N_MCGen);
    fTree->Branch("MC_Px", &fEvent.MC_Px);
    fTree->Branch("MC_Py", &fEvent.MC_Py);
    fTree->Branch("MC_Pz", &fEvent.MC_Pz);
    fTree->Branch("MC_X", &fEvent.MC_X);
    fTree->Branch("MC_Y", &fEvent.MC_Y);
    fTree->Branch("MC_Z", &fEvent.MC_Z);
    fTree->Branch("MC_Xf", &fEvent.MC_Xf);
    fTree->Branch("MC_Yf", &fEvent.MC_Yf);
    fTree->Branch("MC_Zf", &fEvent.MC_Zf);
    fTree->Branch("MC_PID", &fEvent.MC_PID);
    fTree->Branch("MC_Mother", &fEvent.MC_Mother);
    fTree->Branch("MC_PID_Mother", &fEvent.MC_PID_Mother);
    fTree->Branch("MC_PID_GrandMother", &fEvent.MC_PID_GrandMother);
    fTree->Branch("MC_NDaughters", &fEvent.MC_NDaughters);
    fTree->Branch("MC_FirstDau", &fEvent.MC_FirstDau);
    fTree->Branch("MC_LastDau", &fEvent.MC_LastDau);
    fTree->Branch("MC_Status", &fEvent.MC_Status);
    fTree->Branch("MC_isSignal", &fEvent.MC_isSignal);
    /* MC Rec. Tracks */
    fTree->Branch("N_MCRec", &fEvent.N_MCRec);
    fTree->Branch("Idx_True", &fEvent.Idx_True);
    fTree->Branch("Rec_Px", &fEvent.Rec_Px);
    fTree->Branch("Rec_Py", &fEvent.Rec_Py);
    fTree->Branch("Rec_Pz", &fEvent.Rec_Pz);
    fTree->Branch("Rec_Charge", &fEvent.Rec_Charge);
    fTree->Branch("Rec_IP_wrtPV", &fEvent.Rec_IP_wrtPV);
    fTree->Branch("Rec_NSigmaPion", &fEvent.Rec_NSigmaPion);
    fTree->Branch("Rec_NSigmaProton", &fEvent.Rec_NSigmaProton);
    fTree->Branch("Rec_NClustersTPC", &fEvent.Rec_NClustersTPC);
    fTree->Branch("Rec_isDuplicate", &fEvent.Rec_isDuplicate);
    fTree->Branch("Rec_isSimilar", &fEvent.Rec_isSimilar);
    fTree->Branch("Rec_isSignal", &fEvent.Rec_isSignal);
    fTree->Branch("Rec_HelixParam0", &fEvent.Rec_HelixParam0);
    fTree->Branch("Rec_HelixParam1", &fEvent.Rec_HelixParam1);
    fTree->Branch("Rec_HelixParam2", &fEvent.Rec_HelixParam2);
    fTree->Branch("Rec_HelixParam3", &fEvent.Rec_HelixParam3);
    fTree->Branch("Rec_HelixParam4", &fEvent.Rec_HelixParam4);
    fTree->Branch("Rec_HelixParam5", &fEvent.Rec_HelixParam5);
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
    fTree->Branch("V0_Chi2", &fEvent.V0_Chi2);
    fTree->Branch("V0_DCA_Daughters", &fEvent.V0_DCA_Daughters);
    fTree->Branch("V0_IP_wrtPV", &fEvent.V0_IP_wrtPV);
    fTree->Branch("V0_CPA_wrtPV", &fEvent.V0_CPA_wrtPV);
    fTree->Branch("V0_ArmAlpha", &fEvent.V0_ArmAlpha);
    fTree->Branch("V0_ArmPt", &fEvent.V0_ArmPt);
    fTree->Branch("V0_DecayLength", &fEvent.V0_DecayLength);
    fTree->Branch("V0_DCA_wrtPV", &fEvent.V0_DCA_wrtPV);
}

/*
 Push one of the MC generated particles into the Event_tt object
*/
void AliAnalysisTaskSexaquark::MCGen_PushBack(Int_t evt_mc) {

    Int_t idx_mc = fVec_Idx_MCGen[evt_mc];
    AliMCParticle* mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

    Int_t evt_mother;
    Int_t mother_pid;
    Int_t grandmother_pid;
    MCGen_GetMotherInfo(mcPart, evt_mother, mother_pid, grandmother_pid);

    Int_t n_daughters;
    Int_t evt_first_dau;
    Int_t evt_last_dau;
    MCGen_GetDaughtersInfo(mcPart, n_daughters, evt_first_dau, evt_last_dau);

    Float_t x_f = -999.;
    Float_t y_f = -999.;
    Float_t z_f = -999.;
    if (n_daughters) {
        AliMCParticle* mcDaughter = (AliMCParticle*)fMC->GetTrack(mcPart->GetDaughterFirst());
        x_f = mcDaughter->Xv();
        y_f = mcDaughter->Yv();
        z_f = mcDaughter->Zv();
    }

    Bool_t is_signal = MCGen_IsSignal(mcPart);

#if DEBUG_MODE == 1
    AliInfoF("%6i %6i %6i %6i %6i %6i %6i %6i %6i", evt_mc, idx_mc, mcPart->PdgCode(), is_signal,  //
             evt_mother, mother_pid, grandmother_pid,                                              //
             evt_first_dau, evt_last_dau);
#endif

    fEvent.MC_Px.push_back(mcPart->Px());
    fEvent.MC_Py.push_back(mcPart->Py());
    fEvent.MC_Pz.push_back(mcPart->Pz());
    fEvent.MC_X.push_back(mcPart->Xv());
    fEvent.MC_Y.push_back(mcPart->Yv());
    fEvent.MC_Z.push_back(mcPart->Zv());
    fEvent.MC_Xf.push_back(x_f);
    fEvent.MC_Yf.push_back(y_f);
    fEvent.MC_Zf.push_back(z_f);
    fEvent.MC_PID.push_back(mcPart->PdgCode());
    fEvent.MC_Mother.push_back(evt_mother);
    fEvent.MC_PID_Mother.push_back(mother_pid);
    fEvent.MC_PID_GrandMother.push_back(grandmother_pid);
    fEvent.MC_NDaughters.push_back(n_daughters);
    fEvent.MC_FirstDau.push_back(evt_first_dau);
    fEvent.MC_LastDau.push_back(evt_last_dau);
    fEvent.MC_isSignal.push_back(is_signal);
    fEvent.MC_Status.push_back(mcPart->MCStatusCode());
}

/*
 Push one of the MC reconstructed tracks into the Event_tt object
*/
void AliAnalysisTaskSexaquark::MCRec_PushBack(Int_t evt_track) {

    Int_t idx_track = fVec_Idx_MCRec[evt_track];
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

    Int_t idx_true = TMath::Abs(track->GetLabel());
    AliMCParticle* mcTrue = (AliMCParticle*)fMC->GetTrack(idx_true);

    Int_t evt_true = fMap_Evt_MCGen[idx_true];
    Bool_t is_signal = MCGen_IsSignal(mcTrue);
    Bool_t is_duplicate = fVec_MCRec_IsDuplicate[evt_track];
    Bool_t is_similar = kFALSE;  // PENDING! fVec_MCRec_IsSimilar[evt_track];
    Float_t impact_param = MCRec_GetImpactParameter(track);
    Double_t helix_param[6];
    track->GetHelixParameters(helix_param, fESD->GetMagneticField());

#if DEBUG_MODE > 0
    AliInfoF("%6i %6i %6i %6i %6i", evt_track, idx_track, evt_true, idx_true, is_duplicate);
#endif

    fEvent.Idx_True.push_back(evt_true);
    fEvent.Rec_Px.push_back(track->Px());
    fEvent.Rec_Py.push_back(track->Py());
    fEvent.Rec_Pz.push_back(track->Pz());
    fEvent.Rec_Charge.push_back(track->Charge());
    fEvent.Rec_IP_wrtPV.push_back(impact_param);
    fEvent.Rec_NSigmaPion.push_back(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    fEvent.Rec_NSigmaProton.push_back(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    fEvent.Rec_NClustersTPC.push_back(track->GetTPCncls());
    fEvent.Rec_isDuplicate.push_back(is_duplicate);
    fEvent.Rec_isSimilar.push_back(is_similar);
    fEvent.Rec_isSignal.push_back(is_signal);
    fEvent.Rec_HelixParam0.push_back(helix_param[0]);
    fEvent.Rec_HelixParam1.push_back(helix_param[1]);
    fEvent.Rec_HelixParam2.push_back(helix_param[2]);
    fEvent.Rec_HelixParam3.push_back(helix_param[3]);
    fEvent.Rec_HelixParam4.push_back(helix_param[4]);
    fEvent.Rec_HelixParam5.push_back(helix_param[5]);
}

/*
 Push one of the found V0s into the Event_tt object
*/
void AliAnalysisTaskSexaquark::V0_PushBack(Int_t evt_v0) {

    V0_tt this_V0 = fVec_V0s[evt_v0];

#if DEBUG_MODE == 2
    Int_t evt_pos = this_V0.Idx_Pos;
    Int_t evt_neg = this_V0.Idx_Neg;
    Int_t idx_pos = fVec_Idx_MCRec[evt_pos];
    Int_t idx_neg = fVec_Idx_MCRec[evt_neg];
    AliInfoF("%6i %6i %6i %6i", evt_pos, evt_neg, idx_pos, idx_neg);
#endif

    fEvent.Idx_Pos.push_back(this_V0.Idx_Pos);
    fEvent.Idx_Neg.push_back(this_V0.Idx_Neg);
    fEvent.V0_Px.push_back(this_V0.Px);
    fEvent.V0_Py.push_back(this_V0.Py);
    fEvent.V0_Pz.push_back(this_V0.Pz);
    fEvent.V0_X.push_back(this_V0.X);
    fEvent.V0_Y.push_back(this_V0.Y);
    fEvent.V0_Z.push_back(this_V0.Z);
    fEvent.Pos_Px.push_back(this_V0.Pos_Px);
    fEvent.Pos_Py.push_back(this_V0.Pos_Py);
    fEvent.Pos_Pz.push_back(this_V0.Pos_Pz);
    fEvent.Neg_Px.push_back(this_V0.Neg_Px);
    fEvent.Neg_Py.push_back(this_V0.Neg_Py);
    fEvent.Neg_Pz.push_back(this_V0.Neg_Pz);
    fEvent.V0_isSignal.push_back(this_V0.isSignal);
    fEvent.V0_E_asK0.push_back(this_V0.E_asK0);
    fEvent.V0_E_asAL.push_back(this_V0.E_asAL);
    fEvent.V0_couldBeK0.push_back(this_V0.couldBeK0);
    fEvent.V0_couldBeAL.push_back(this_V0.couldBeAL);
    fEvent.V0_Chi2.push_back(this_V0.Chi2);
    fEvent.V0_DCA_Daughters.push_back(this_V0.DCA_Daughters);
    fEvent.V0_IP_wrtPV.push_back(this_V0.IP_wrtPV);
    fEvent.V0_CPA_wrtPV.push_back(this_V0.CPA_wrtPV);
    fEvent.V0_ArmAlpha.push_back(this_V0.ArmAlpha);
    fEvent.V0_ArmPt.push_back(this_V0.ArmPt);
    fEvent.V0_DecayLength.push_back(this_V0.DecayLength);
    fEvent.V0_DCA_wrtPV.push_back(this_V0.DCA_wrtPV);
}

/*
 Clear all indices and vectors loaded in the Event_tt object
*/
void AliAnalysisTaskSexaquark::ClearEvent() {
    /* MC particles */
    fEvent.N_MCGen = 0;
    fEvent.MC_Px.clear();
    fEvent.MC_Py.clear();
    fEvent.MC_Pz.clear();
    fEvent.MC_X.clear();
    fEvent.MC_Y.clear();
    fEvent.MC_Z.clear();
    fEvent.MC_Xf.clear();
    fEvent.MC_Yf.clear();
    fEvent.MC_Zf.clear();
    fEvent.MC_PID.clear();
    fEvent.MC_Mother.clear();
    fEvent.MC_PID_Mother.clear();
    fEvent.MC_PID_GrandMother.clear();
    fEvent.MC_NDaughters.clear();
    fEvent.MC_FirstDau.clear();
    fEvent.MC_LastDau.clear();
    fEvent.MC_Status.clear();
    fEvent.MC_isSignal.clear();
    /* MC Rec. Tracks */
    fEvent.N_MCRec = 0;
    fEvent.Idx_True.clear();
    fEvent.Rec_Px.clear();
    fEvent.Rec_Py.clear();
    fEvent.Rec_Pz.clear();
    fEvent.Rec_Charge.clear();
    fEvent.Rec_IP_wrtPV.clear();
    fEvent.Rec_NSigmaPion.clear();
    fEvent.Rec_NSigmaProton.clear();
    fEvent.Rec_NClustersTPC.clear();
    fEvent.Rec_isDuplicate.clear();
    fEvent.Rec_isSimilar.clear();
    fEvent.Rec_isSignal.clear();
    fEvent.Rec_HelixParam0.clear();
    fEvent.Rec_HelixParam1.clear();
    fEvent.Rec_HelixParam2.clear();
    fEvent.Rec_HelixParam3.clear();
    fEvent.Rec_HelixParam4.clear();
    fEvent.Rec_HelixParam5.clear();
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
    fEvent.V0_Chi2.clear();
    fEvent.V0_DCA_Daughters.clear();
    fEvent.V0_IP_wrtPV.clear();
    fEvent.V0_CPA_wrtPV.clear();
    fEvent.V0_ArmAlpha.clear();
    fEvent.V0_ArmPt.clear();
    fEvent.V0_DecayLength.clear();
    fEvent.V0_DCA_wrtPV.clear();
}
