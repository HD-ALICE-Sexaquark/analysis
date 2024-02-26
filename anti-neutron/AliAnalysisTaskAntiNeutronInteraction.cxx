#include "AliAnalysisTaskAntiNeutronInteraction.h"

ClassImp(AliAnalysisTaskAntiNeutronInteraction);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliAnalysisTaskAntiNeutronInteraction::AliAnalysisTaskAntiNeutronInteraction()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fStruckNucleonPDG(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      isMcIdxSignal(),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAwrtPV(0.),
      kMin_AN_Mass(0.),
      kMax_AN_Mass(0.),
      kMin_AN_Pt(0.),
      kMax_AN_Eta(0.),
      kMax_AN_Rapidity(0.),
      kMin_AN_DecayLength(0.),
      kMin_AN_Radius(0.),
      kMin_AN_CPAwrtPV(0.),
      kMax_AN_CPAwrtPV(0.),
      kMin_AN_DCAwrtPV(0.),
      kMax_AN_DCAwrtPV(0.),
      kMax_AN_DCAbtwV0s(0.),
      kMax_AN_DCAv0aSV(0.),
      kMax_AN_DCAv0bSV(0.),
      kMax_AN_DCAv0anegSV(0.),
      kMax_AN_DCAv0aposSV(0.),
      kMax_AN_DCAv0bnegSV(0.),
      kMax_AN_DCAv0bposSV(0.),
      kMax_AN_Chi2ndf(0.) {}

/*
 Constructor, called locally.
*/
AliAnalysisTaskAntiNeutronInteraction::AliAnalysisTaskAntiNeutronInteraction(const char* name, Bool_t IsMC)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fStruckNucleonPDG(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      isMcIdxSignal(),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAwrtPV(0.),
      kMin_AN_Mass(0.),
      kMax_AN_Mass(0.),
      kMin_AN_Pt(0.),
      kMax_AN_Eta(0.),
      kMax_AN_Rapidity(0.),
      kMin_AN_DecayLength(0.),
      kMin_AN_Radius(0.),
      kMin_AN_CPAwrtPV(0.),
      kMax_AN_CPAwrtPV(0.),
      kMin_AN_DCAwrtPV(0.),
      kMax_AN_DCAwrtPV(0.),
      kMax_AN_DCAbtwV0s(0.),
      kMax_AN_DCAv0aSV(0.),
      kMax_AN_DCAv0bSV(0.),
      kMax_AN_DCAv0anegSV(0.),
      kMax_AN_DCAv0aposSV(0.),
      kMax_AN_DCAv0bnegSV(0.),
      kMax_AN_DCAv0bposSV(0.),
      kMax_AN_Chi2ndf(0.) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/*
 Destructor.
*/
AliAnalysisTaskAntiNeutronInteraction::~AliAnalysisTaskAntiNeutronInteraction() {
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskAntiNeutronInteraction::UserCreateOutputObjects() {

    Initialize();

    /** Add mandatory routines **/

    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) {
        AliFatal("ERROR: AliAnalysisManager couldn't be found.");
    }

    AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
    if (!inputHandler) {
        AliFatal("ERROR: AliESDInputHandler couldn't be found.");
    }

    fPIDResponse = inputHandler->GetPIDResponse();

    /** Prepare output histograms **/

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    PrepareTracksHistograms();
    PrepareAntiNeutronHistograms();

    PostData(1, fOutputListOfHists);
}

/*
 Define and add tracks' histograms to the histograms output list.
 - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `fOutputListOfHists`
*/
void AliAnalysisTaskAntiNeutronInteraction::PrepareTracksHistograms() {

    TString tracks_stages[2] = {"MCGen", "Found"};
    TString tracks_sets[5] = {"All", "Selected", "True", "Secondary", "Signal"};
    Int_t tracks_species[4] = {-2212, 321, 211, -211};
    std::map<Int_t, TString> tracks_name = {{-2212, "AntiProton"}, {321, "PosKaon"}, {211, "PiPlus"}, {-211, "PiMinus"}};

    const Int_t N_tracks_props = 11;
    TString tracks_props[N_tracks_props] = {"Pt",           "Pz",           "Eta",             //
                                            "DCAwrtPV",     "NTPCClusters", "Chi2/NClusters",  //
                                            "NSigmaProton", "NSigmaKaon",   "NSigmaPion",
                                            "Status",       "GoldenChi2"};
    Int_t tracks_nbins[N_tracks_props] = {100, 100, 100, 100, 100, 100,  //
                                          100, 100, 100, 16,  100};
    Double_t tracks_min[N_tracks_props] = {0.,   -20., -4.,  0., 0., 0.,  //
                                           -7.5, -7.5, -7.5, 0,  -10};
    Double_t tracks_max[N_tracks_props] = {10., 20., 4.,  200., 1000., 50.,  //
                                           7.5, 7.5, 7.5, 16,   40};

    for (Int_t& species : tracks_species) {
        fHist_Tracks_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", tracks_name[species].Data()), "", 15, 0., 15.);
        fOutputListOfHists->Add(fHist_Tracks_Bookkeep[species]);
    }

    for (TString& stage : tracks_stages) {
        for (TString& set : tracks_sets) {
            for (Int_t& species : tracks_species) {
                for (Int_t prop_idx = 0; prop_idx < N_tracks_props; prop_idx++) {

                    if (set == "Selected") continue;  // set name reserved for 2D histograms

                    if (stage == "MCGen" && set == "True") continue;

                    if (stage == "MCGen" && (tracks_props[prop_idx] == "DCAwrtPV" || tracks_props[prop_idx] == "NTPCClusters" ||
                                             tracks_props[prop_idx] == "Chi2/NClusters" || tracks_props[prop_idx] == "NSigmaProton" ||
                                             tracks_props[prop_idx] == "NSigmaKaon" || tracks_props[prop_idx] == "NSigmaPion" ||
                                             tracks_props[prop_idx] == "Status" || tracks_props[prop_idx] == "GoldenChi2")) {
                        continue;
                    }

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), tracks_name[species].Data(), tracks_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, species, tracks_props[prop_idx]);
                    fHist_Tracks[histKey] = new TH1F(histName, "", tracks_nbins[prop_idx], tracks_min[prop_idx], tracks_max[prop_idx]);
                    fOutputListOfHists->Add(fHist_Tracks[histKey]);
                }
            }
        }
    }

    for (TString& set : tracks_sets) {
        TString histName = Form("Tracks_%s_TPCsignal", set.Data());
        f2DHist_Tracks_TPCsignal[set] = new TH2F(histName, "", 200, 0., 10., 200, 0., 400.);
        fOutputListOfHists->Add(f2DHist_Tracks_TPCsignal[set]);
    }
}

/*
 Define and add anti-neutron histograms to the histograms output list.
*/
void AliAnalysisTaskAntiNeutronInteraction::PrepareAntiNeutronHistograms() {

    fHist_All_AntiNeutron_Pt = new TH1F("All_AntiNeutron_Pt", "", 100, 0., 10.);
    fHist_Signal_AntiNeutron_Pt = new TH1F("Signal_AntiNeutron_Pt", "", 100, 0., 10.);
    fHist_Signal_AntiNeutron_Radius = new TH1F("Signal_AntiNeutron_Radius", "", 100, 0., 200.);
    fHist_True_AntiNeutron_NDaughters = new TH1F("True_AntiNeutron_NDaughters", "", 100, 0., 100.);
    fHist_True_AntiNeutron_Bookkeep = new TH1F("True_AntiNeutron_Bookkeep", "", 10, 0., 10.);

    fOutputListOfHists->Add(fHist_All_AntiNeutron_Pt);
    fOutputListOfHists->Add(fHist_Signal_AntiNeutron_Pt);
    fOutputListOfHists->Add(fHist_Signal_AntiNeutron_Radius);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_NDaughters);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_Bookkeep);
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fMagneticField`, `fPrimaryVertex`, `fSourceOfV0s`, `fReactionChannel`, `fOutputListOfTrees`,
 `fOutputListOfHists`
*/
void AliAnalysisTaskAntiNeutronInteraction::UserExec(Option_t*) {

    // load MC generated event
    fMC = MCEvent();

    if (!fMC) {
        AliFatal("ERROR: AliMCEvent couldn't be found.");
    }

    fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());

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

    if (fIsMC) {
        ProcessMCGen();
        // ProcessSignalInteractions();
    }

    // ProcessTracks();

    // if (fIsMC) ProcessFindableAntiNeutrons();

    // KalmanAntiNeutronFinder();

    ClearContainers();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfHists);
}

/*
 Initialize analysis task.
*/
void AliAnalysisTaskAntiNeutronInteraction::Initialize() {

    DefineTracksCuts("Standard");
    DefineAntiNeutronCuts("Standard");

    AliInfo("Initializing AliAnalysisTaskAntiNeutronInteraction...");
    AliInfo("INPUT OPTIONS:");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);
}

/*
 Define track selection cuts.
 - Input: `cuts_option`
*/
void AliAnalysisTaskAntiNeutronInteraction::DefineTracksCuts(TString cuts_option) {

    kMax_NSigma_Pion = 3.;
    kMax_NSigma_Kaon = 3.;
    kMax_NSigma_Proton = 3.;

    kMax_Track_Eta = 1.;
    kMin_Track_NTPCClusters = 50;
    kMax_Track_Chi2PerNTPCClusters = 7.;
    kTurnedOn_Track_StatusCuts = kTRUE;
    kTurnedOn_Track_RejectKinks = kFALSE;  // PENDING: need to study further
    kMin_Track_DCAwrtPV = 2.;

    kMin_Track_Pt[2212] = 0.4;
    kMin_Track_Pt[321] = 0.0;  // PENDING: to check when analyzing the next reaction channels
    kMin_Track_Pt[211] = 0.2;
}

/*
 Define anti-neutron candidates selection cuts.
 - Uses: `fSourceOfV0s`
 - Input: `cuts_option`
*/
void AliAnalysisTaskAntiNeutronInteraction::DefineAntiNeutronCuts(TString cuts_option) {

    /* Common */

    // none so far

    // kMin_Sexa_Radius = 40;
    // kMin_Sexa_CPAwrtPV = 0.99;
    // kMax_Sexa_CPAwrtPV = 1.;
    // kMax_Sexa_DCAbtwV0s = 2.;
    // kMax_Sexa_DCAv0aSV = 2.;
    // kMax_Sexa_DCAv0bSV = 2.;
    // kMax_Sexa_DCAv0anegSV = 10.;
    // kMax_Sexa_DCAv0aposSV = 10.;
    // kMax_Sexa_DCAv0bnegSV = 10.;
    // kMax_Sexa_DCAv0bposSV = 10.;
    // kMax_Sexa_Chi2ndf = 0.5;
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

/*
 Loop over MC particles in a single event. Store the indices of the signal particles.
 - Uses: `fMC`, `fPDG`, `fMC_PrimaryVertex`
*/
void AliAnalysisTaskAntiNeutronInteraction::ProcessMCGen() {

    AliMCParticle* mcPart;
    Int_t pdg_mc;
    Int_t n_daughters = 0;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    Float_t pt;
    Float_t radius;

    Bool_t is_signal;

    std::map<TString, Int_t> counter_reaction_channels;

    /*** Loop over MC gen. particles in a single event ***/

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

        /* Get PDG Code and store it */

        pdg_mc = mcPart->PdgCode();
        getPdgCode_fromMcIdx[mcIdx] = pdg_mc;

        /* Protection in case of a particle with no momentum */

        if (mcPart->P() < 1E-6) {
            continue;
        }

        /* Focus on anti-neutrons */

        if (pdg_mc != -2112) continue;

        /* Get particle's properties */

        pt = mcPart->Pt();
        fHist_True_AntiNeutron_Bookkeep->Fill(0);
        fHist_True_AntiNeutron_NDaughters->Fill(n_daughters);
        fHist_All_AntiNeutron_Pt->Fill(pt);

        n_daughters = mcPart->GetNDaughters();

        if (!n_daughters) continue;

        fHist_True_AntiNeutron_Bookkeep->Fill(1);
        fHist_Signal_AntiNeutron_Pt->Fill(pt);

        std::vector<Int_t> pdg_antineutron_daughters;

        /* Fill, then, the PDG code of its daughters */

        for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

            mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
            pdg_dau = mcDaughter->PdgCode();

            pdg_antineutron_daughters.push_back(pdg_dau);

            if (mcIdxDaughter == mcPart->GetDaughterFirst()) {
                radius = TMath::Sqrt(TMath::Power(mcDaughter->Xv(), 2) + TMath::Power(mcDaughter->Yv(), 2));
                fHist_Signal_AntiNeutron_Radius->Fill(radius);
            }
        }

        /* Sort PDG codes of daughters, group them into a string and update counter */

        std::sort(pdg_antineutron_daughters.begin(), pdg_antineutron_daughters.end());

        TString sorted_names_antineutron_daughters = "";
        for (Int_t idx_pdg_dau = 0; idx_pdg_dau < (Int_t)pdg_antineutron_daughters.size(); idx_pdg_dau++) {
            if (idx_pdg_dau == 0) {
                sorted_names_antineutron_daughters += Form("%i", pdg_antineutron_daughters[idx_pdg_dau]);
            } else {
                sorted_names_antineutron_daughters += Form(", %i", pdg_antineutron_daughters[idx_pdg_dau]);
            }
        }

        counter_reaction_channels[sorted_names_antineutron_daughters]++;
    }  // end of loop over MC particles

    /* Print result of counter of reaction channels */

    for (auto& x : counter_reaction_channels) {
        AliInfoF("Reaction channel: %s N -> %s, count: %i", fPDG.GetParticle(-2112)->GetName(), x.first.Data(), x.second);
    }
}

/*                */
/**    Tracks    **/
/*** ========== ***/

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisTaskAntiNeutronInteraction::ProcessTracks() {

    AliESDtrack* track;

    Int_t mcIdx;
    Int_t mcPdgCode;
    Int_t mcIdxOfTrueV0;

    Float_t pt, pz, eta;
    Float_t impar_pv[2], dca_wrt_pv;
    Float_t n_tpc_clusters;
    Float_t chi2_over_nclusters;
    Float_t n_sigma_proton;
    Float_t n_sigma_kaon;
    Float_t n_sigma_pion;
    Float_t golden_chi2;

    /* Loop over tracks in a single event */

    for (Int_t esdIdxTrack = 0; esdIdxTrack < fESD->GetNumberOfTracks(); esdIdxTrack++) {

        /* Get track */

        track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxTrack));

        /* Get MC info */

        mcIdx = TMath::Abs(track->GetLabel());
        mcPdgCode = getPdgCode_fromMcIdx[mcIdx];

        /* Fill containers */

        getMcIdx_fromEsdIdx[esdIdxTrack] = mcIdx;

        isEsdIdxSignal[esdIdxTrack] = isMcIdxSignal[mcIdx];

        /* Fill histograms before track selection */

        f2DHist_Tracks_TPCsignal["All"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        /* Apply track selection */

        if (!PassesTrackSelection(track)) continue;

        doesEsdIdxPassTrackSelection[esdIdxTrack] = kTRUE;

        /* Fill containers -- related to findable anti-neutrons */

        /*
                if (isMcIdxSignal[mcIdx]) {
                    getReactionIdx_fromEsdIdx[esdIdxTrack] =
                        isMcIdxDaughterOfSignalV0[mcIdx] ? getReactionIdx_fromMcIdx[getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdx]] :
           getReactionIdx_fromMcIdx[mcIdx]; getEsdIdx_fromReactionIdx[getReactionIdx_fromEsdIdx[esdIdxTrack]].push_back(esdIdxTrack);
                }
        */

        /* Get track properties */

        pt = track->Pt();
        pz = track->Pz();
        eta = track->Eta();

        Float_t xy_impar_wrt_pv, z_impar_wrt_pv;
        track->GetImpactParameters(xy_impar_wrt_pv, z_impar_wrt_pv);  // pre-calculated DCA w.r.t. PV
        dca_wrt_pv = TMath::Sqrt(xy_impar_wrt_pv * xy_impar_wrt_pv + z_impar_wrt_pv * z_impar_wrt_pv);

        n_tpc_clusters = track->GetTPCNcls();  // note: capital N
        chi2_over_nclusters = n_tpc_clusters ? track->GetTPCchi2() / (Double_t)n_tpc_clusters : 999.;
        n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        golden_chi2 = track->GetChi2TPCConstrainedVsGlobal(fPrimaryVertex);

        /* Fill histograms */

        f2DHist_Tracks_TPCsignal["Selected"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        std::vector<Int_t> possiblePID;
        if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton && track->Charge() < 0.) possiblePID.push_back(-2212);
        if (TMath::Abs(n_sigma_kaon) < kMax_NSigma_Kaon && track->Charge() > 0.) possiblePID.push_back(321);
        if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion && track->Charge() < 0.) possiblePID.push_back(-211);
        if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion && track->Charge() > 0.) possiblePID.push_back(211);

        for (Int_t& esdPdgCode : possiblePID) {

            if (esdPdgCode == -2212) esdIndicesOfAntiProtonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 321) esdIndicesOfPosKaonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == -211) esdIndicesOfPiMinusTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 211) esdIndicesOfPiPlusTracks.push_back(esdIdxTrack);

            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pt")]->Fill(pt);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pz")]->Fill(pz);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Eta")]->Fill(eta);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Chi2/NClusters")]->Fill(chi2_over_nclusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "GoldenChi2")]->Fill(golden_chi2);
            PlotStatus(track, "All", esdPdgCode);

            if (mcPdgCode == esdPdgCode) {
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pt")]->Fill(pt);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pz")]->Fill(pz);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Eta")]->Fill(eta);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Chi2/NClusters")]->Fill(chi2_over_nclusters);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
                fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "GoldenChi2")]->Fill(golden_chi2);
                PlotStatus(track, "True", esdPdgCode);
                f2DHist_Tracks_TPCsignal["True"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
            }

            if (mcPdgCode == esdPdgCode && isEsdIdxSignal[esdIdxTrack]) {
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pt")]->Fill(pt);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pz")]->Fill(pz);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Eta")]->Fill(eta);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Chi2/NClusters")]->Fill(chi2_over_nclusters);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
                fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "GoldenChi2")]->Fill(golden_chi2);
                PlotStatus(track, "Signal", esdPdgCode);
                f2DHist_Tracks_TPCsignal["Signal"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
            }
        }

        /* Fill vectors */

        /*
         mcIdxOfTrueV0 = getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdx];
                getMcIdxOfTrueV0_fromEsdIdx[esdIdxTrack] = mcIdxOfTrueV0;

                if (track->Charge() < 0.) {
                    getEsdIdxOfRecNegDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0] = esdIdxTrack;
                } else {
                    getEsdIdxOfRecPosDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0] = esdIdxTrack;
                }
         */
    }  // end of loop over tracks
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskAntiNeutronInteraction::PassesTrackSelection(AliESDtrack* track) {

    Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

    // >> particle identification
    std::vector<Int_t> possiblePID;
    if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton) possiblePID.push_back(2212);
    if (TMath::Abs(n_sigma_kaon) < kMax_NSigma_Kaon) possiblePID.push_back(321);
    if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion) possiblePID.push_back(211);
    if (!possiblePID.size()) return kFALSE;

    // >> eta
    if (TMath::Abs(track->Eta()) > kMax_Track_Eta) return kFALSE;

    // >> TPC clusters
    if (track->GetTPCNcls() < kMin_Track_NTPCClusters) return kFALSE;

    // >> chi2 per TPC cluster
    if (track->GetTPCchi2() / (Double_t)track->GetTPCNcls() > kMax_Track_Chi2PerNTPCClusters) return kFALSE;

    // >> TPC and ITS status
    Bool_t tpc_status =
        (track->GetStatus() & AliESDtrack::kTPCin) && (track->GetStatus() & AliESDtrack::kTPCout) && (track->GetStatus() & AliESDtrack::kTPCrefit);
    Bool_t its_status =
        !(track->GetStatus() & AliESDtrack::kITSin) && !(track->GetStatus() & AliESDtrack::kITSout) && !(track->GetStatus() & AliESDtrack::kITSrefit);
    if (kTurnedOn_Track_StatusCuts && (!tpc_status || !its_status)) return kFALSE;

    // >> reject kinks
    if (kTurnedOn_Track_RejectKinks && track->GetKinkIndex(0) != 0) return kFALSE;

    // >> DCA w.r.t Primary Vertex
    Float_t xy_impar, z_impar;
    track->GetImpactParameters(xy_impar, z_impar);
    Float_t dca_wrt_pv = TMath::Sqrt(xy_impar * xy_impar + z_impar * z_impar);
    if (dca_wrt_pv < kMin_Track_DCAwrtPV) return kFALSE;

    // >> reject low-pT pions
    // -- PENDING: this will also affect tracks with multiple PID hypotheses:
    // -- if a track both id as pion and proton, it will be cut by the highest pT cut
    for (Int_t& esdPdgCode : possiblePID) {
        if (track->Pt() < kMin_Track_Pt[esdPdgCode]) return kFALSE;
    }

    return kTRUE;
}

/*
 Plot the status of the track.
 - Input: `track`, `stage`, `esdPdgCode`
*/
void AliAnalysisTaskAntiNeutronInteraction::PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode) {

    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};

    for (Int_t i = 0; i < 16; i++) {
        if ((track->GetStatus() & StatusCollection[i])) fHist_Tracks[std::make_tuple("Found", set, esdPdgCode, "Status")]->Fill(i);
    }
}

/*                                             */
/**  Anti-Neutrons -- That Require True Info  **/
/*** ======================================= ***/

/*
  Find the anti-neutrons that have all of their *final-state particles* reconstructed.
*/
void AliAnalysisTaskAntiNeutronInteraction::ProcessFindableAntiNeutrons() {

    // PENDING
}

/*                                   */
/**  Anti-Neutron -- Kalman Filter  **/
/*** ============================= ***/

/*
 Using Kalman Filter, reconstruct anti-neutron candidates
 - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`, `fMagneticField`
 - Containers: `kf[AntiLambdas|KaonsZeroShort]`, `getEsdIdxOf[Neg|Pos]Dau_from[AntiLambda|KaonZeroShort]Idx`
*/
void AliAnalysisTaskAntiNeutronInteraction::KalmanAntiNeutronFinder() {

    // PENDING
}

/*
 Apply anti-neutron candidate cuts.
 - Use: `fPrimaryVertex`
 - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskAntiNeutronInteraction::PassesAntiNeutronCuts(KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                                    KFParticle kfAntiLambda, KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg,
                                                                    KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg,
                                                                    KFParticle kfKaonZeroShortPos) {

    return kTRUE;
}

/*                            */
/**  Mathematical Functions  **/
/*** ====================== ***/

/*
 Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 - Input: `lvParticle`, `X`, `Y`, `Z`, `refPointX`, `refPointY`, `refPointZ`
 - Return: `cos(theta)`
 (Based on `AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()`)
*/
Double_t AliAnalysisTaskAntiNeutronInteraction::CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z,  //
                                                                    Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(lvParticle.Angle(posRelativeToRef));
}

/*
 Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
*/
Double_t AliAnalysisTaskAntiNeutronInteraction::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z,  //
                                                                    Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

/*
 Correct initialization of a KFParticle.
 (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
*/
KFParticle AliAnalysisTaskAntiNeutronInteraction::CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge) {

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
 Correct initialization of a KFVertex.
 (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
*/
KFVertex AliAnalysisTaskAntiNeutronInteraction::CreateKFVertex(const AliVVertex& vertex) {

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

/*
 Transport a KFParticle to the point of closest approach w.r.t. another KFParticle.
 - Uses: `fPDG`
 - Input: `kfThis`, `kfOther`, `pdgThis`, `chargeThis`
 - Return: `kfTransported`
*/
KFParticle AliAnalysisTaskAntiNeutronInteraction::TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis) {

    float dsdr[4][6];
    float dS[2];
    kfThis.GetDStoParticle(kfOther, dS, dsdr);

    float mP[8], mC[36];
    kfThis.Transport(dS[0], dsdr[0], mP, mC);

    float mM = fPDG.GetParticle(pdgThis)->Mass();
    float mQ = chargeThis;  // only valid for charged particles with Q = +/- 1

    KFParticle kfTransported;
    kfTransported.Create(mP, mC, mQ, mM);

    return kfTransported;
}

/*               */
/**  Utilities  **/
/*** ========= ***/

/*
 Clear all the containers.
 - Uses:
*/
void AliAnalysisTaskAntiNeutronInteraction::ClearContainers() {

    isMcIdxSignal.clear();
    getPdgCode_fromMcIdx.clear();

    getReactionIdx_fromEsdIdx.clear();
    getEsdIdx_fromReactionIdx.clear();

    doesEsdIdxPassTrackSelection.clear();

    esdIndicesOfAntiProtonTracks.clear();
    esdIndicesOfPosKaonTracks.clear();
    esdIndicesOfPiPlusTracks.clear();
    esdIndicesOfPiMinusTracks.clear();

    getMcIdx_fromEsdIdx.clear();

    isEsdIdxSignal.clear();

    getEsdIdxOfNegDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfPosDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfNegDau_fromKaonZeroShortIdx.clear();
}
