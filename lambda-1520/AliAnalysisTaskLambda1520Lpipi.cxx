#include "AliAnalysisTaskLambda1520Lpipi.h"

ClassImp(AliAnalysisTaskLambda1520Lpipi);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliAnalysisTaskLambda1520Lpipi::AliAnalysisTaskLambda1520Lpipi()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fMagneticField(0.),
      fPDG(),
      fOutputListOfHists(0),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAwrtPV(0.),
      kMin_LStar_Mass(0.),
      kMax_LStar_Mass(0.),
      kMin_LStar_Pt(0.),
      kMax_LStar_Eta(0.),
      kMax_LStar_Rapidity(0.),
      kMin_LStar_DecayLength(0.),
      kMin_LStar_Radius(0.),
      kMin_LStar_CPAwrtPV(0.),
      kMax_LStar_CPAwrtPV(0.),
      kMin_LStar_DCAwrtPV(0.),
      kMax_LStar_DCAwrtPV(0.),
      kMax_LStar_DCAbtwV0s(0.),
      kMax_LStar_DCAv0aSV(0.),
      kMax_LStar_DCAv0bSV(0.),
      kMax_LStar_DCAv0anegSV(0.),
      kMax_LStar_DCAv0aposSV(0.),
      kMax_LStar_DCAv0bnegSV(0.),
      kMax_LStar_DCAv0bposSV(0.),
      kMax_LStar_Chi2ndf(0.) {}

/*
 Constructor, called locally.
*/
AliAnalysisTaskLambda1520Lpipi::AliAnalysisTaskLambda1520Lpipi(const char* name, Bool_t IsMC)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fMagneticField(0.),
      fPDG(),
      fOutputListOfHists(0),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAwrtPV(0.),
      kMin_LStar_Mass(0.),
      kMax_LStar_Mass(0.),
      kMin_LStar_Pt(0.),
      kMax_LStar_Eta(0.),
      kMax_LStar_Rapidity(0.),
      kMin_LStar_DecayLength(0.),
      kMin_LStar_Radius(0.),
      kMin_LStar_CPAwrtPV(0.),
      kMax_LStar_CPAwrtPV(0.),
      kMin_LStar_DCAwrtPV(0.),
      kMax_LStar_DCAwrtPV(0.),
      kMax_LStar_DCAbtwV0s(0.),
      kMax_LStar_DCAv0aSV(0.),
      kMax_LStar_DCAv0bSV(0.),
      kMax_LStar_DCAv0anegSV(0.),
      kMax_LStar_DCAv0aposSV(0.),
      kMax_LStar_DCAv0bnegSV(0.),
      kMax_LStar_DCAv0bposSV(0.),
      kMax_LStar_Chi2ndf(0.) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/*
 Destructor.
*/
AliAnalysisTaskLambda1520Lpipi::~AliAnalysisTaskLambda1520Lpipi() {
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskLambda1520Lpipi::UserCreateOutputObjects() {

    CheckForInputErrors();
    Initialize();

    DefineTracksCuts("Standard");
    DefineV0Cuts("Standard");
    DefineLambda1520Cuts("Standard");

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

    /** Prepare histograms **/

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    PrepareTracksHistograms();
    PrepareV0Histograms();
    PrepareLambda1520Histograms();

    PostData(1, fOutputListOfHists);
}

/*
 Search for possible contradictions in the input options
*/
void AliAnalysisTaskLambda1520Lpipi::CheckForInputErrors() {
    // PENDING
}

/*
 PENDING.
*/
void AliAnalysisTaskLambda1520Lpipi::Initialize() {
    AliInfo("Initializing AliAnalysisTaskLambda1520Lpipi...");
    AliInfo("INPUT OPTIONS:");
    AliInfoF(" >> IsMC = %i", (Int_t)fIsMC);
}

/*
 Define and add tracks' histograms to the histograms output list.
 - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `f2DHist_Tracks_TPCsignal`, `fOutputListOfHists`
*/
void AliAnalysisTaskLambda1520Lpipi::PrepareTracksHistograms() {

    TString tracks_stages[2] = {"MCGen", "Found"};
    TString tracks_sets[4] = {"All", "Selected", "True", "Signal"};
    Int_t tracks_species[4] = {2212, -2212, 211, -211};
    std::map<Int_t, TString> tracks_name = {{2212, "Proton"}, {-2212, "AntiProton"}, {211, "PiPlus"}, {-211, "PiMinus"}};

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

                    if (stage == "MCGen" && set == "True") continue;  // only keep "MCGen All" and "MCGen Signal"

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
 Define and add V0s' histograms to the histograms output list.
 - Uses: `fHist_V0s_Bookkeep`, `fHist_V0s`, `f2DHist_V0s_ArmenterosPodolanski`, `fOutputListOfHists`
*/
void AliAnalysisTaskLambda1520Lpipi::PrepareV0Histograms() {

    TString V0_stages[3] = {"MCGen", "Findable", "Found"};
    TString V0_sets[3] = {"All", "True", "Signal"};
    Int_t V0_species[2] = {3122, -3122};
    std::map<Int_t, TString> V0_name = {{3122, "Lambda"}, {-3122, "AntiLambda"}};

    const Int_t N_V0_props = 16;
    TString V0_props[N_V0_props] = {"Mass", "Radius", "CPAwrtPV", "DCAwrtPV", "DCAbtwDau", "DCAnegV0", "DCAposV0",     "DecayLength",
                                    "Zv",   "Eta",    "Pt",       "Pz",       "ArmQt",     "ArmAlpha", "OpeningAngle", "Chi2ndf"};
    Int_t V0_nbins[N_V0_props] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Double_t V0_min[N_V0_props] = {0., 0., 0., 0., 0., 0., 0., 0., -50., -4., 0., -20., 0., -1., -0.5 * TMath::Pi(), 0.};
    Double_t V0_max[N_V0_props] = {10., 250., 1., 200., 2., 2., 2., 500., 50., 4., 10., 20., 0.3, 1., 0.5 * TMath::Pi(), 0.5};

    for (Int_t& species : V0_species) {
        fHist_V0s_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", V0_name[species].Data()), "", 15, 0., 15.);
        fOutputListOfHists->Add(fHist_V0s_Bookkeep[species]);
    }

    for (TString& stage : V0_stages) {
        for (TString& set : V0_sets) {
            for (Int_t& species : V0_species) {
                for (Int_t prop_idx = 0; prop_idx < N_V0_props; prop_idx++) {

                    if (stage == "Findable" && set == "All") continue;  // only keep "Findable True" and "Findable Signal"
                    if (stage == "MCGen" && set == "True") continue;    // only keep "MCGen All" and "MCGen Signal"

                    if ((stage == "MCGen" || stage == "Findable") && (V0_props[prop_idx] == "DCAbtwDau" || V0_props[prop_idx] == "DCAnegV0" ||
                                                                      V0_props[prop_idx] == "DCAposV0" || V0_props[prop_idx] == "Chi2ndf")) {
                        continue;
                    }

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), V0_name[species].Data(), V0_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, species, V0_props[prop_idx]);
                    Float_t minEdge = V0_props[prop_idx] == "Mass" ? kMin_V0_Mass[species] : V0_min[prop_idx];
                    Float_t maxEdge = V0_props[prop_idx] == "Mass" ? kMax_V0_Mass[species] : V0_max[prop_idx];
                    fHist_V0s[histKey] = new TH1F(histName, "", V0_nbins[prop_idx], minEdge, maxEdge);
                    fOutputListOfHists->Add(fHist_V0s[histKey]);
                }
            }
        }
    }

    for (TString& set : V0_sets) {
        TString histName = Form("Found_%s_V0s_ArmenterosPodolanski", set.Data());
        f2DHist_V0s_ArmenterosPodolanski[set] = new TH2F(histName, "", 200, -1., 1., 200, 0., 0.3);
        fOutputListOfHists->Add(f2DHist_V0s_ArmenterosPodolanski[set]);
    }
}

/*
 Define and add Lambda(1520) histograms to the histograms output list.
 - Uses: `fHist_Lambdas1520_Bookkeep`, `fHist_Lambdas1520`, `fOutputListOfHists`
 PENDING: change properties!!
*/
void AliAnalysisTaskLambda1520Lpipi::PrepareLambda1520Histograms() {

    // "LS" stands for LambdaStar ~ Lambda(1520)
    TString LS_stages[3] = {"MCGen", "Findable", "Found"};
    TString LS_sets[2] = {"All", "Signal"};
    Int_t LS_species[2] = {3124, -3124};
    std::map<Int_t, TString> LS_name = {{3124, "Lambda1520"}, {-3124, "AntiLambda1520"}};

    const Int_t N_LS_props = 18;
    TString LS_props[N_LS_props] = {"Mass",     "Pt",          "Pz",          "Radius",      "Zv",          "Eta",       //
                                    "Rapidity", "DecayLength", "CPAwrtPV",    "DCAwrtPV",    "DCAbtwV0s",   "DCAv0aSV",  //
                                    "DCAv0bSV", "DCAv0anegSV", "DCAv0aposSV", "DCAv0bnegSV", "DCAv0bposSV", "Chi2ndf"};
    Int_t LS_nbins[N_LS_props] = {100, 100, 100, 100, 100, 100,  //
                                  100, 100, 100, 100, 100, 100,  //
                                  100, 100, 100, 100, 100, 100};
    Double_t LS_min[N_LS_props] = {-10., -10., -50., 0., -50., -2.,  //
                                   -5.,  0.,   0.,   0., 0.,   0.,   //
                                   0.,   0.,   0.,   0., 0.,   0.};
    Double_t LS_max[N_LS_props] = {10., 10.,  50., 200., 50., 2.,  //
                                   5.,  500., 1.,  10.,  2.,  2.,  //
                                   2.,  10.,  10., 10.,  10., 1.};

    for (Int_t& species : LS_species) {
        fHist_Lambdas1520_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", LS_name[species].Data()), "", 15, 0., 15.);
        fOutputListOfHists->Add(fHist_Lambdas1520_Bookkeep[species]);
    }

    for (TString& stage : LS_stages) {
        for (Int_t& species : LS_species) {
            for (TString& set : LS_sets) {
                for (Int_t prop_idx = 0; prop_idx < N_LS_props; prop_idx++) {

                    if (stage == "Findable" && set == "Signal") continue;  // only keep "Findable All"
                    if (stage == "MCGen" && set == "Signal") continue;     // only keep "MCGen All"

                    if ((stage == "MCGen" || stage == "Findable") &&
                        (LS_props[prop_idx] == "DCAbtwV0s" || LS_props[prop_idx] == "DCAv0aSV" ||       //
                         LS_props[prop_idx] == "DCAv0bSV" || LS_props[prop_idx] == "DCAv0anegSV" ||     //
                         LS_props[prop_idx] == "DCAv0aposSV" || LS_props[prop_idx] == "DCAv0bnegSV" ||  //
                         LS_props[prop_idx] == "DCAv0bposSV" || LS_props[prop_idx] == "Chi2ndf")) {
                        continue;
                    }

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), LS_name[species].Data(), LS_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, species, LS_props[prop_idx]);
                    fHist_Lambdas1520[histKey] = new TH1F(histName, "", LS_nbins[prop_idx], LS_min[prop_idx], LS_max[prop_idx]);
                    fOutputListOfHists->Add(fHist_Lambdas1520[histKey]);
                }
            }
        }
    }
}

/*
 Define track selection cuts.
 - Input: `cuts_option`
*/
void AliAnalysisTaskLambda1520Lpipi::DefineTracksCuts(TString cuts_option) {

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
    kMin_Track_Pt[211] = 0.2;
}

/*
 Define Lambda1520 candidate cuts.
*/
void AliAnalysisTaskLambda1520Lpipi::DefineV0Cuts(TString cuts_option) {

    // none so far
}

/*
 Define Lambda1520 candidate cuts.
*/
void AliAnalysisTaskLambda1520Lpipi::DefineLambda1520Cuts(TString cuts_option) {

    // none so far
}

/*
 Main function, called per each event.
*/
void AliAnalysisTaskLambda1520Lpipi::UserExec(Option_t*) {

    // load MC generated event
    fMC = MCEvent();

    if (fIsMC && !fMC) {
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

    if (fIsMC) ProcessMCGen();

    // ProcessTracks();

    // if (fIsMC) {
    // ProcessFindableV0s();
    // ProcessFindableLambdas1520();
    // }

    /* Lambda -> pi-,P */

    // KalmanV0Finder(3122, -211, 2212);

    /* AntiLambda -> AntiP,pi+ */

    // KalmanV0Finder(-3122, -2212, 211);

    /* Lambda(1520) -> X,pi-,pi+ */

    // KalmanV0Finder(0, -211, 211);

    /* Lambda(1520) -> Lambda,pi-,pi+ */

    // pdgTracks = {-211, 2212, -211, 211};
    // KalmanLambda1520Finder(kfLambdas, idxLambdaDaughters, kfPionPairs, idxPionPairDaughters, pdgTracks, kfLambdas1520);
    // KalmanLambda1520Finder();

    /* AntiLambda(1520) -> AntiLambda,pi-,pi+ */

    // pdgTracks = {-2212, 211, -211, 211};
    // KalmanLambda1520Finder(kfAntiLambdas, idxAntiLambdaDaughters, kfPionPairs, idxPionPairDaughters, pdgTracks, kfAntiLambdas1520);
    // KalmanLambda1520Finder();

    ClearContainers();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfHists);
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

/*
 Search for the relevant MC particles in a single event.
 - input: fMC
 - output: Indices_MCGen_FS, Indices_MCGen_FS_Signal
*/
void AliAnalysisTaskLambda1520Lpipi::ProcessMCGen() {

    AliMCParticle* mcPart;
    Int_t pdg_mc;

    Int_t n_daughters;
    Int_t mcIdxNeutralDaughter;
    Int_t mcIdxNegDaughter;
    Int_t mcIdxPosDaughter;

    TVector3 originVertex;
    TVector3 decayVertex;

    Float_t pt;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;

    /*** Loop over MC gen. particles in a single event ***/

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);
        pdg_mc = mcPart->PdgCode();
        getPdgCode_fromMcIdx[mcIdx] = pdg_mc;

        n_daughters = mcPart->GetNDaughters();
        mcIdxNeutralDaughter = 0;
        mcIdxNegDaughter = 0;
        mcIdxPosDaughter = 0;

        /** Protection in case of a particle with no momentum **/

        if (mcPart->P() < 1E-6) {
            continue;
        }

        /** Get particle's properties **/

        originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
        decayVertex.SetXYZ(0., 0., 0.);

        pt = mcPart->Pt();
        if (n_daughters) {
            GetDaughtersInfo(mcPart, mcIdxNeutralDaughter, mcIdxNegDaughter, mcIdxPosDaughter, decayVertex);
            cpa_wrt_pv = CosinePointingAngle(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                             fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                      fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        }

        /** Determine if they're signal or not **/

        isMcIdxSignal[mcIdx] = TMath::Abs(pdg_mc) == 3124;
        if (isMcIdxSignal[mcIdx]) getLStarIdx_fromMcIdx[mcIdx] = mcIdx;
        if (mcPart->GetMother() >= 0) {
            isMcIdxSignal[mcIdx] = isMcIdxSignal[mcPart->GetMother()];
            if (isMcIdxSignal[mcIdx]) {
                getLStarIdx_fromMcIdx[mcIdx] = getLStarIdx_fromMcIdx[mcPart->GetMother()];
                getMcIdx_fromLStarIdx[getLStarIdx_fromMcIdx[mcPart->GetMother()]].push_back(mcIdx);
            }
        }

        /** Fill histograms, indices vectors and maps **/

        /* (Anti)Lambda1520 */

        if (TMath::Abs(pdg_mc) == 3124) {

            fHist_Lambdas1520_Bookkeep[pdg_mc]->Fill(0);

            if (n_daughters && mcIdxNeutralDaughter && mcIdxNegDaughter && mcIdxPosDaughter) {

                fHist_Lambdas1520_Bookkeep[pdg_mc]->Fill(1);
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Radius")]->Fill(decayVertex.Perp());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                // fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "Rapidity")]->Fill(mcPart->Rapidity()); // PENDING
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "DecayLength")]->Fill((originVertex - decayVertex).Mag());
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_Lambdas1520[std::make_tuple("MCGen", "All", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
            }
        }  // end Lambda1520 condition

        /* (Anti)Lambda */

        if (TMath::Abs(pdg_mc) == 3122) {

            fHist_V0s_Bookkeep[pdg_mc]->Fill(0);

            if (isMcIdxSignal[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(1);

            if (isMcIdxSignal[mcIdx] && n_daughters && mcIdxNegDaughter && mcIdxPosDaughter) {
                fHist_V0s_Bookkeep[pdg_mc]->Fill(2);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Radius")]->Fill(decayVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayLength")]->Fill((originVertex - decayVertex).Mag());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                // PENDING!
                // fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmQt")]->Fill(  //
                // Calculate_ArmPt(mcPart->Px(), mcPart->Py(), mcPart->Pz(),          //
                // mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz()));
                // fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmAlpha")]->Fill(  //
                // Calculate_ArmAlpha(mcPart->Px(), mcPart->Py(), mcPart->Pz(),          //
                //    mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz(),    //
                //    mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz()));
            }
        }  // end of (Anti)Lambda condition

        /* (Anti)protons and Charged Pions */

        if (TMath::Abs(pdg_mc) == 2212 || TMath::Abs(pdg_mc) == 211) {

            fHist_Tracks_Bookkeep[pdg_mc]->Fill(0);
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());

            if (isMcIdxSignal[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(1);
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }

        }  // end of (anti)protons and charged pions condition

    }  // end of loop over MC particles
}

/*
 Loop over the daughters of an MC particle.
 - Input: `mcPart`
 - Output: `mcIdxNeutralDaughter`, `mcIdxNegDaughter`, `mcIdxPosDaughter`, `decayVertex`
*/
void AliAnalysisTaskLambda1520Lpipi::GetDaughtersInfo(AliMCParticle* mcPart, Int_t& mcIdxNeutralDaughter, Int_t& mcIdxNegDaughter,
                                                      Int_t& mcIdxPosDaughter, TVector3& decayVertex) {

    Int_t pdg_mc = mcPart->PdgCode();

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

        mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
        pdg_dau = mcDaughter->PdgCode();

        if (TMath::Abs(pdg_mc) == 3124) {
            if (TMath::Abs(pdg_dau) == 3122) mcIdxNeutralDaughter = mcIdxDaughter;
            if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
            if (pdg_dau == -211) mcIdxNegDaughter = mcIdxDaughter;
        }

        if (pdg_mc == 3122) {
            if (pdg_dau == -211) mcIdxNegDaughter = mcIdxDaughter;
            if (pdg_dau == 2212) mcIdxPosDaughter = mcIdxDaughter;
        }

        if (pdg_mc == -3122) {
            if (pdg_dau == -2212) mcIdxNegDaughter = mcIdxDaughter;
            if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
        }

        /* In any case, get the decay vertex of the V0 */

        decayVertex.SetXYZ(mcDaughter->Xv(), mcDaughter->Yv(), mcDaughter->Zv());
    }
}

/*            */
/**  Tracks  **/
/*** ====== ***/

/*
 Loop over the reconstructed tracks in a single event.
 IMPORTANT: don't trust track->GetID()
 - uses: fESD, fMC
 - input: Indices_MCGen_FS_Signal (pending)
 - output: idxPiPlusTracks,idxPiMinusTracks,idxKaonPlusTracks,idxKaonMinusTracks,idxProtonTracks,idxAntiProtonTracks
*/
void AliAnalysisTaskLambda1520Lpipi::ProcessTracks() {

    AliESDtrack* track;

    // loop over MC rec. tracks in a single event
    for (Int_t idx_track = 0; idx_track < fESD->GetNumberOfTracks(); idx_track++) {

        track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

        /* Apply track selection cuts */

        if (!PassesTrackSelection(track)) continue;

    }  // end of loop over tracks
}

/*
 Plot the status of the track.
 - Input: `track`, `stage`, `esdPdgCode`
*/
void AliAnalysisTaskLambda1520Lpipi::PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode) {
#ifdef PENDING
    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};

    for (Int_t i = 0; i < 16; i++) {
        if ((track->GetStatus() & StatusCollection[i])) fHist_Tracks[std::make_tuple("Found", set, esdPdgCode, "Status")]->Fill(i);
    }
#endif
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskLambda1520Lpipi::PassesTrackSelection(AliESDtrack* track) {

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

/*                                       */
/**  Findables -- Require MC Gen. Info  **/
/*** ================================= ***/

/*
 Find all V0s.
 Store all lambdas and kaons with no ambiguity.
 - no dependence on PID, nor ambiguity on measurements
 - only dependence: if they were reconstructed or not
*/
void AliAnalysisTaskLambda1520Lpipi::ProcessFindableV0s() {
#ifdef PENDING

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

        antiLambdaNonPi0Channel =
            (mc_firstdau->PdgCode() == -2212 && mc_lastdau->PdgCode() == 211) || (mc_firstdau->PdgCode() == 211 && mc_lastdau->PdgCode() == -2212);

        neutralKaonNonPi0Channel =
            (mc_firstdau->PdgCode() == -211 && mc_lastdau->PdgCode() == 211) || (mc_firstdau->PdgCode() == 211 && mc_lastdau->PdgCode() == -211);

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

    }  // end of loop over MC gen. mothers
#endif
}

/*
 PENDING.
*/
void AliAnalysisTaskLambda1520Lpipi::ProcessFindableLambdas1520() {
    // PENDING
}

/*         */
/**  V0s  **/
/*** === ***/

/*
 Find all true V0s which had both of their daughters reconstructed.
 - uses: fESD, fMC, fMagneticField
 - input: idxNegativeTracks, idxPositiveTracks, pdgTrackNeg, pdgTrackPos
 - output: kfV0s, idxDaughters
*/
void AliAnalysisTaskLambda1520Lpipi::KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {
#ifdef PENDING
    AliESDtrack* esdTrackNeg;
    AliESDtrack* esdTrackPos;

    Bool_t successful_daughters_check;  // PENDING: I still don't understand this...

    Double_t impar_neg[2], impar_pos[2];

    Double_t radius;
    Double_t cpa_wrt_pv;
    Double_t decay_length;
    Double_t arm_qt;
    Double_t arm_alpha;
    Double_t opening_angle;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /*  Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;
    std::vector<KFParticleMother>* kfFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;
    std::vector<Bool_t>* isFoundV0IdxATrueV0;
    std::unordered_map<Int_t, Int_t>* getMcIdx_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        kfFoundV0s = &kfAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
        isFoundV0IdxATrueV0 = &isAntiLambdaIdxATrueV0;
        getMcIdx_fromFoundV0Idx = &getMcIdx_fromAntiLambdaIdx;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        kfFoundV0s = &kfKaonsZeroShort;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
        isFoundV0IdxATrueV0 = &isKaonZeroShortIdxATrueV0;
        getMcIdx_fromFoundV0Idx = &getMcIdx_fromKaonZeroShortIdx;
    }

    /* Loop over all possible pairs of tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            /* Sanity check */

            if (esdIdxNeg == esdIdxPos) continue;

            /* Get tracks */

            esdTrackNeg = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
            esdTrackPos = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

            /* Kalman Filter */

            KFParticle kfDaughterNeg = CreateKFParticle(*esdTrackNeg, fPDG.GetParticle(pdgTrackNeg)->Mass(), (Int_t)esdTrackNeg->Charge());
            KFParticle kfDaughterPos = CreateKFParticle(*esdTrackPos, fPDG.GetParticle(pdgTrackPos)->Mass(), (Int_t)esdTrackPos->Charge());

            KFParticleMother kfV0;
            kfV0.AddDaughter(kfDaughterNeg);
            kfV0.AddDaughter(kfDaughterPos);

            /* Transport V0 and daughters */

            kfV0.TransportToDecayVertex();

            KFParticle kfTransportedNeg = TransportKFParticle(kfDaughterNeg, kfDaughterPos, pdgTrackNeg, (Int_t)esdTrackNeg->Charge());
            KFParticle kfTransportedPos = TransportKFParticle(kfDaughterPos, kfDaughterNeg, pdgTrackPos, (Int_t)esdTrackPos->Charge());

            /* Reconstruct V0 */

            lvTrackNeg.SetXYZM(kfTransportedNeg.Px(), kfTransportedNeg.Py(), kfTransportedNeg.Pz(), fPDG.GetParticle(pdgTrackNeg)->Mass());
            lvTrackPos.SetXYZM(kfTransportedPos.Px(), kfTransportedPos.Py(), kfTransportedPos.Pz(), fPDG.GetParticle(pdgTrackPos)->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            /* Apply cuts */

            if (!PassesV0Cuts(pdgV0, kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;

            /* Calculate variables */

            radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
            cpa_wrt_pv = CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                             fPrimaryVertex->GetZ());
            decay_length = TMath::Sqrt((kfV0.GetX() - fPrimaryVertex->GetX()) * (kfV0.GetX() - fPrimaryVertex->GetX()) +
                                       (kfV0.GetY() - fPrimaryVertex->GetY()) * (kfV0.GetY() - fPrimaryVertex->GetY()) +
                                       (kfV0.GetZ() - fPrimaryVertex->GetZ()) * (kfV0.GetZ() - fPrimaryVertex->GetZ()));
            arm_qt = Calculate_ArmPt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = Calculate_ArmAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                           lvTrackPos.Py(), lvTrackPos.Pz());
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Store V0 */

            kfFoundV0s->push_back(kfV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxPos;

            isFoundV0IdxATrueV0->push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                           getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(9);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex)));
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos)));
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2")]->Fill(kfV0.GetChi2());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2ndf")]->Fill((Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(arm_alpha, arm_qt);

            if (!isFoundV0IdxATrueV0->back()) continue;

            (*getMcIdx_fromFoundV0Idx)[kfFoundV0s->size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

            if (getPdgCode_fromMcIdx[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]] != pdgV0) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(10);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex)));
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos)));
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2")]->Fill(kfV0.GetChi2());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2ndf")]->Fill((Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(12);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex)));
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos)));
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0)));
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2")]->Fill(kfV0.GetChi2());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2ndf")]->Fill((Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(arm_alpha, arm_qt);

        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
#endif
}

/*
 Apply cuts to a V0 candidate.
 - Input: `pdgV0`, `kfV0`, `kfDaughterNeg`, `kfDaughterPos`, `lvV0`, `lvTrackNeg`, `lvTrackPos`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskLambda1520Lpipi::PassesV0Cuts(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfDaughterNeg, KFParticle kfDaughterPos,
                                                    TLorentzVector lvV0, TLorentzVector lvTrackNeg, TLorentzVector lvTrackPos) {
#ifdef PENDING
    Double_t Radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;

    Double_t DecayLength = TMath::Sqrt(TMath::Power(kfV0.GetX() - fPrimaryVertex->GetX(), 2) + TMath::Power(kfV0.GetY() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(kfV0.GetZ() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;

    Double_t CPAwrtPV =
        CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;

    Double_t DCAwrtPV = LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(),
                                     fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;

    Double_t DCAbtwDau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;

    Double_t DCAnegV0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;

    Double_t DCAposV0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;

    Double_t ArmPt = Calculate_ArmPt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
    Double_t ArmAlpha = Calculate_ArmAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                           lvTrackPos.Py(), lvTrackPos.Pz());
    Double_t ArmPtOverAlpha = TMath::Abs(ArmPt / ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;

    Double_t Chi2ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();
    if (kMax_V0_Chi2ndf[pdgV0] && Chi2ndf > kMax_V0_Chi2ndf[pdgV0]) return kFALSE;
#endif
    return kTRUE;
}

/*                */
/**  Lambda1520  **/
/*** ========== ***/

/*
 Find all Lambda1520 candidates.
 // (pending) add more description
 - input: kfFirstV0s, idxFirstV0Daughters, kfSecondV0s, idxSecondV0Daughters, pdgDaughters
 - output: kfLambda1520s
*/
void AliAnalysisTaskLambda1520Lpipi::KalmanLambda1520Finder() {
#ifdef PENDING

    Double_t impar[2];
    Bool_t fStatusT0PropToDCA;
    Bool_t fStatusT1PropToDCA;
    Bool_t fStatusT2PropToDCA;
    Bool_t fStatusT3PropToDCA;

    Int_t idxV0A_NegDau, idxV0A_PosDau;
    Int_t idxV0B_NegDau, idxV0B_PosDau;

    TLorentzVector lvTrack0, lvTrack1, lvTrack2, lvTrack3;
    TLorentzVector lvLambda, lvPionPair;
    TLorentzVector lvLambda1520;

    Bool_t successful_daughters_check;

    Double_t dca_sv_pv, cpa_sv_pv, dca_v0a_sv, cpa_v0a_sv, dca_t2_sv, dca_t3_sv, dca_v0a_v0b;

    for (Int_t idxV0A = 0; idxV0A < (Int_t)kfFirstV0s.size(); idxV0A++) {
        for (Int_t idxV0B = 0; idxV0B < (Int_t)kfSecondV0s.size(); idxV0B++) {

            idxV0A_NegDau = idxFirstV0Daughters[idxV0A][0];
            idxV0A_PosDau = idxFirstV0Daughters[idxV0A][1];
            idxV0B_NegDau = idxSecondV0Daughters[idxV0B][0];
            idxV0B_PosDau = idxSecondV0Daughters[idxV0B][1];

            // (protection) avoid repetition of daughters
            if (idxV0A_NegDau == idxV0B_NegDau || idxV0A_PosDau == idxV0B_PosDau) continue;

            AliESDtrack* esdTrack0 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_NegDau));  // not used
            AliESDtrack* esdTrack1 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_PosDau));  // not used
            AliESDtrack* esdTrack2 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_NegDau));
            AliESDtrack* esdTrack3 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_PosDau));

            /* Kalman Filter */

            KFParticleMother kfV0A = kfFirstV0s[idxV0A];
            KFParticleMother kfV0B = kfSecondV0s[idxV0B];  // not used

            KFParticle kfTrack0 = CreateKFParticle(*esdTrack0,                                           //
                                                   fPDG.GetParticle(pdgDaughters[0])->Mass(),            //
                                                   (Int_t)fPDG.GetParticle(pdgDaughters[0])->Charge());  // not used
            KFParticle kfTrack1 = CreateKFParticle(*esdTrack1,                                           //
                                                   fPDG.GetParticle(pdgDaughters[1])->Mass(),            //
                                                   (Int_t)fPDG.GetParticle(pdgDaughters[1])->Charge());  // not used
            KFParticle kfTrack2 = CreateKFParticle(*esdTrack2,                                           //
                                                   fPDG.GetParticle(pdgDaughters[2])->Mass(),            //
                                                   (Int_t)fPDG.GetParticle(pdgDaughters[2])->Charge());
            KFParticle kfTrack3 = CreateKFParticle(*esdTrack3,                                 //
                                                   fPDG.GetParticle(pdgDaughters[3])->Mass(),  //
                                                   (Int_t)fPDG.GetParticle(pdgDaughters[3])->Charge());

            KFParticleMother kfLambda1520;
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

            lvLambda = lvTrack0 + lvTrack1;
            lvPionPair = lvTrack2 + lvTrack3;
            lvLambda1520 = lvTrack0 + lvTrack1 + lvTrack2 + lvTrack3;

            /* Get properties */

            dca_sv_pv = LinePointDCA(lvLambda1520.Px(), lvLambda1520.Py(), lvLambda1520.Pz(),        //
                                     kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ(),  //
                                     fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            cpa_sv_pv = CosinePointingAngle(lvLambda1520, kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ(),  //
                                            fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_v0a_sv = LinePointDCA(lvLambda.Px(), lvLambda.Py(), lvLambda.Pz(),  //
                                      kfV0A.GetX(), kfV0A.GetY(), kfV0A.GetZ(),     //
                                      kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ());
            cpa_v0a_sv = CosinePointingAngle(lvLambda, kfV0A.GetX(), kfV0A.GetY(), kfV0A.GetZ(),  //
                                             kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ());
            dca_t2_sv = kfTrack2.GetDistanceFromVertexXY(kfLambda1520);
            dca_t3_sv = kfTrack3.GetDistanceFromVertexXY(kfLambda1520);
            dca_v0a_v0b = kfV0A.GetDistanceFromParticleXY(kfV0B);  // (pending)

            // (cut) lambdas1520 selection
            if (!PassesLambda1520Cuts_KF(kfLambda1520, kfV0A, kfTrack2, kfTrack3, lvLambda1520, lvLambda, lvTrack2, lvTrack3)) continue;

            /* Histograms! */
        }
    }
#endif
}

/*
 Apply Lambda(1520) candidate cuts.
*/
Bool_t AliAnalysisTaskLambda1520Lpipi::PassesLambda1520Cuts(KFParticleMother kfLambda1520, TLorentzVector lvLambda1520, KFParticle kfAntiLambda,
                                                            KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg, KFParticle kfAntiLambdaPos,
                                                            KFParticle kfKaonZeroShortNeg, KFParticle kfKaonZeroShortPos) {

    // PENDING

    return kTRUE;
}

/*               */
/**  Utilities  **/
/*** ========= ***/

/*
 Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 This function stores the point of closest approach, and returns its distance to the PV
 [more info at https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line]
*/
Double_t AliAnalysisTaskLambda1520Lpipi::LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                      Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                      Double_t PV_X, Double_t PV_Y, Double_t PV_Z) {
    TVector3 V0_Momentum(V0_Px, V0_Py, V0_Pz);

    TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
    TVector3 PrimVertex(PV_X, PV_Y, PV_Z);

    TVector3 CrossProduct = (PrimVertex - V0Vertex).Cross(V0_Momentum);

    return CrossProduct.Mag() / V0_Momentum.Mag();
}

/*
 This gives the Armenteros-Podolanski alpha:
 asymmetry of longitudinal momentum w.r.t. mother between daughters
 [based on AliRoot/STEER/ESD/AliESDv0::AlphaV0()]
*/
Double_t AliAnalysisTaskLambda1520Lpipi::ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
                                                         Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz,  //
                                                         Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);
    TVector3 momPos(Pos_Px, Pos_Py, Pos_Pz);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    // (protection)
    if (lQlPos + lQlNeg == 0.) {
        return 2;
    }  // closure
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 This gives the Armenteros-Podolanski q_T:
 transverse momentum of either daughter w.r.t. the mother's momentum
 [based on AliRoot/STEER/ESD/AliESDv0::PtArmV0()]
*/
Double_t AliAnalysisTaskLambda1520Lpipi::ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                      Double_t Dau_Px, Double_t Dau_Py, Double_t Dau_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momDau(Dau_Px, Dau_Py, Dau_Pz);

    return momDau.Perp(momTot);
}

/*
 Calculates the cosine of the pointing angle
 of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 [based on AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()]
*/
Double_t AliAnalysisTaskLambda1520Lpipi::CosinePointingAngle(TLorentzVector particles_mom, Double_t X, Double_t Y, Double_t Z,  //
                                                             Double_t refPointX, Double_t refPointY, Double_t refPointZ) {

    // method 1
    // vector between the reference point and the particle's vertex
    Double_t deltaPos[3];
    deltaPos[0] = X - refPointX;
    deltaPos[1] = Y - refPointY;
    deltaPos[2] = Z - refPointZ;

    TVector3 hh(deltaPos[0], deltaPos[1], deltaPos[2]);
    return TMath::Cos(particles_mom.Angle(hh));

    // method 2
    /*
    TVector3 nom_mom = particles_mom.Vect();
    TVector3 nom_pos(X, Y, Z);
    TVector3 ref_pos(refPointX, refPointY, refPointZ);
    TVector3 delta = nom_pos - ref_pos;

    return delta.Dot(nom_mom) / (delta.Mag() * nom_mom.Mag());
    */
}

/*
 Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
*/
Double_t AliAnalysisTaskLambda1520Lpipi::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z,  //
                                                             Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

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

/*
 Transport a KFParticle to the point of closest approach w.r.t. another KFParticle.
 - Uses: `fPDG`
 - Input: `kfThis`, `kfOther`, `pdgThis`, `chargeThis`
 - Return: `kfTransported`
*/
KFParticle AliAnalysisTaskLambda1520Lpipi::TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis) {

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

/*

*/
void AliAnalysisTaskLambda1520Lpipi::ClearContainers() {
    getPdgCode_fromMcIdx.clear();

    isMcIdxSignal.clear();

    getLStarIdx_fromMcIdx.clear();
    getMcIdx_fromLStarIdx.clear();
}
