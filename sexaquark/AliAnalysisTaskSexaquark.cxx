#include "AliAnalysisTaskSexaquark.h"

ClassImp(AliAnalysisTaskSexaquark);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fSourceOfV0s(),
      fSimulationSet(""),
      fReactionID(0),
      fInjectedSexaquarkMass(0),
      fReactionChannel(""),
      fStruckNucleonPDG(0),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      fTree(0),
      fEvent(),
      mcIndicesOfTrueV0s(),
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
      kMin_Sexa_Mass(0.),
      kMax_Sexa_Mass(0.),
      kMin_Sexa_Pt(0.),
      kMax_Sexa_Eta(0.),
      kMax_Sexa_Rapidity(0.),
      kMin_Sexa_DecayLength(0.),
      kMin_Sexa_Radius(0.),
      kMin_Sexa_CPAwrtPV(0.),
      kMax_Sexa_CPAwrtPV(0.),
      kMin_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_Chi2ndf(0.) {}

/*
 Constructor, called locally.
*/
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name, Bool_t IsMC, TString SourceOfV0s, TString SimulationSet)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fSourceOfV0s(SourceOfV0s),
      fSimulationSet(SimulationSet),
      fReactionID(0),
      fInjectedSexaquarkMass(0),
      fReactionChannel(""),
      fStruckNucleonPDG(0),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      fTree(0),
      fEvent(),
      mcIndicesOfTrueV0s(),
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
      kMin_Sexa_Mass(0.),
      kMax_Sexa_Mass(0.),
      kMin_Sexa_Pt(0.),
      kMax_Sexa_Eta(0.),
      kMax_Sexa_Rapidity(0.),
      kMin_Sexa_DecayLength(0.),
      kMin_Sexa_Radius(0.),
      kMin_Sexa_CPAwrtPV(0.),
      kMax_Sexa_CPAwrtPV(0.),
      kMin_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAwrtPV(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_Chi2ndf(0.) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}

/*
 Destructor.
*/
AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() {
    if (fOutputListOfTrees) {
        delete fOutputListOfTrees;
    }
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {

    CheckForInputErrors();
    Initialize();

    DefineTracksCuts("Standard");
    DefineV0Cuts("Standard");
    DefineSexaquarkCuts("Standard");

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

    /** Prepare output  **/

    /* Trees */

    fOutputListOfTrees = new TList();
    fOutputListOfTrees->SetOwner(kTRUE);

    fTree = new TTree("Events", "Tree of Events");
    SetBranches();
    fOutputListOfTrees->Add(fTree);

    PostData(1, fOutputListOfTrees);

    /* Histograms */

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    PrepareTracksHistograms();
    PrepareV0Histograms();
    PrepareAntiSexaquarkHistograms();

    PostData(2, fOutputListOfHists);
}

/*
 Define and add tracks' histograms to the histograms output list.
 - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareTracksHistograms() {

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

                    if (!fIsMC && (stage == "MCGen" || set == "True" || set == "Secondary" || set == "Signal")) continue;

                    if (set == "Selected") continue;  // set name reserved for 2D histograms

                    if (stage == "MCGen" && set == "True") continue;

                    if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S" && species == 321 && set == "Signal") continue;

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
 - Uses: `fSourceOfV0s`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareV0Histograms() {

    TString V0_stages[3] = {"MCGen", "Findable", "Found"};
    TString V0_sets[4] = {"All", "True", "Secondary", "Signal"};
    Int_t V0_species[2] = {-3122, 310};
    std::map<Int_t, TString> V0_name = {{-3122, "AntiLambda"}, {310, "KaonZeroShort"}};

    const Int_t N_V0_props = 18;
    TString V0_props[N_V0_props] = {"Mass", "Radius", "CPAwrtPV", "DCAwrtPV", "DCAbtwDau", "DCAnegV0", "DCAposV0",     "DecayLength",
                                    "Zv",   "Eta",    "Pt",       "Pz",       "ArmQt",     "ArmAlpha", "OpeningAngle", "ImprvDCAbtwDau",
                                    "Chi2", "Chi2ndf"};
    Int_t V0_nbins[N_V0_props] = {100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100};
    Double_t V0_min[N_V0_props] = {0., 0., 0., 0., 0., 0., 0., 0., -50., -4., 0., -20., 0., -1., -0.5 * TMath::Pi(), 0., 0., 0.};
    Double_t V0_max[N_V0_props] = {10., 250., 1., 200., 2., 2., 2., 500., 50., 4., 10., 20., 0.3, 1., 0.5 * TMath::Pi(), 1., 2., 0.5};

    for (Int_t& species : V0_species) {
        fHist_V0s_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", V0_name[species].Data()), "", 50, 0., 50.);
        fOutputListOfHists->Add(fHist_V0s_Bookkeep[species]);
    }

    for (TString& stage : V0_stages) {
        for (TString& set : V0_sets) {
            for (Int_t& species : V0_species) {
                for (Int_t prop_idx = 0; prop_idx < N_V0_props; prop_idx++) {

                    if (!fIsMC && (stage == "MCGen" || stage == "Findable" || set == "True" || set == "Secondary" || set == "Signal")) continue;

                    if (fSourceOfV0s == "kalman" && V0_props[prop_idx] == "ImprvDCAbtwDau") continue;
                    if ((fSourceOfV0s == "offline" || fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "custom") && V0_props[prop_idx] == "Chi2ndf") {
                        continue;
                    }

                    if (stage == "Findable" && set == "All") continue;
                    if (stage == "MCGen" && set == "True") continue;

                    if (stage == "Findable" && (V0_props[prop_idx] == "Mass" || V0_props[prop_idx] == "DCAbtwDau" ||       //
                                                V0_props[prop_idx] == "DCAnegV0" || V0_props[prop_idx] == "DCAposV0" ||    //
                                                V0_props[prop_idx] == "ImprvDCAbtwDau" || V0_props[prop_idx] == "Chi2" ||  //
                                                V0_props[prop_idx] == "Chi2ndf")) {
                        continue;
                    }
                    if (stage == "MCGen" && (V0_props[prop_idx] == "DCAbtwDau" || V0_props[prop_idx] == "DCAnegV0" ||       //
                                             V0_props[prop_idx] == "DCAposV0" || V0_props[prop_idx] == "ImprvDCAbtwDau" ||  //
                                             V0_props[prop_idx] == "Chi2" || V0_props[prop_idx] == "Chi2ndf")) {
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
 Define and add anti-sexaquark histograms to the histograms output list.
 - Uses: `fHist_AntiSexaquarks`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareAntiSexaquarkHistograms() {

    TString AS_stages[3] = {"MCGen", "Findable", "Found"};
    TString AS_sets[2] = {"All", "Signal"};

    const Int_t N_AS_props = 18;
    TString AS_props[N_AS_props] = {"Mass",     "Pt",          "Pz",          "Radius",      "Zv",          "Eta",       //
                                    "Rapidity", "DecayLength", "CPAwrtPV",    "DCAwrtPV",    "DCAbtwV0s",   "DCAv0aSV",  //
                                    "DCAv0bSV", "DCAv0anegSV", "DCAv0aposSV", "DCAv0bnegSV", "DCAv0bposSV", "Chi2ndf"};
    Int_t AS_nbins[N_AS_props] = {100, 100, 100, 100, 100, 100,  //
                                  100, 100, 100, 100, 100, 100,  //
                                  100, 100, 100, 100, 100, 100};
    Double_t AS_min[N_AS_props] = {-10., -10., -50., 0., -50., -2.,  //
                                   -5.,  0.,   0.,   0., 0.,   0.,   //
                                   0.,   0.,   0.,   0., 0.,   0.};
    Double_t AS_max[N_AS_props] = {10., 10.,  50., 200., 50., 2.,  //
                                   5.,  500., 1.,  10.,  2.,  2.,  //
                                   2.,  10.,  10., 10.,  10., 1.};

    fHist_AntiSexaquarks_Bookkeep = new TH1F("AntiSexaquarks_Bookkeep", "", 25, 0., 25.);
    fOutputListOfHists->Add(fHist_AntiSexaquarks_Bookkeep);

    for (TString& stage : AS_stages) {
        for (TString& set : AS_sets) {
            for (Int_t prop_idx = 0; prop_idx < N_AS_props; prop_idx++) {

                if (!fIsMC && (stage == "MCGen" || stage == "Findable" || set == "Signal")) continue;

                if (stage == "Findable" && set == "Signal") continue;
                if (stage == "MCGen" && set == "Signal") continue;

                if (fSourceOfV0s != "kalman" && AS_props[prop_idx] == "Chi2ndf") continue;

                if ((stage == "MCGen" || stage == "Findable") && (AS_props[prop_idx] == "DCAbtwV0s" || AS_props[prop_idx] == "DCAv0aSV" ||       //
                                                                  AS_props[prop_idx] == "DCAv0bSV" || AS_props[prop_idx] == "DCAv0anegSV" ||     //
                                                                  AS_props[prop_idx] == "DCAv0aposSV" || AS_props[prop_idx] == "DCAv0bnegSV" ||  //
                                                                  AS_props[prop_idx] == "DCAv0bposSV" || AS_props[prop_idx] == "Chi2ndf")) {
                    continue;
                }

                TString histName = Form("%s_%s_AntiSexaquark_%s", stage.Data(), set.Data(), AS_props[prop_idx].Data());
                std::tuple<TString, TString, TString> histKey = std::make_tuple(stage, set, AS_props[prop_idx]);
                fHist_AntiSexaquarks[histKey] = new TH1F(histName, "", AS_nbins[prop_idx], AS_min[prop_idx], AS_max[prop_idx]);
                fOutputListOfHists->Add(fHist_AntiSexaquarks[histKey]);
            }
        }
    }
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fMagneticField`, `fPrimaryVertex`, `fSourceOfV0s`, `fReactionChannel`, `fOutputListOfTrees`,
 `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

    /* Load MC Gen. Event and PV */

    if (fIsMC) {
        fMC = MCEvent();
        if (!fMC) AliFatal("ERROR: AliMCEvent couldn't be found.");
        fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());
    }

    /* Load Reconstructed Event, PV and Magnetic Field */

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) AliFatal("ERROR: AliESDEvent couldn't be found.");
    fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());
    fMagneticField = fESD->GetMagneticField();

    // init KFparticle
    KFParticle::SetField(fMagneticField);

    if (fIsMC) {
        ProcessMCGen();
        ProcessSignalInteractions();
    }

    ProcessTracks();

    if (fIsMC) {
        ProcessFindableV0s();
        ProcessFindableSexaquarks();
    }

    if (fSourceOfV0s == "offline") {
        GetV0sFromESD(kFALSE);
    }

    if (fSourceOfV0s == "on-the-fly") {
        GetV0sFromESD(kTRUE);
    }

    if (fSourceOfV0s == "custom") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            CustomV0Finder(-3122);
            CustomV0Finder(310);
        }
    }

    if (fSourceOfV0s == "offline" || fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "custom") {
        if (fReactionID == 'A') {
            GeoSexaquarkFinder_ChannelA();
        }
    }

    if (fSourceOfV0s == "kalman") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            KalmanV0Finder(-3122, -2212, 211);
            KalmanV0Finder(310, -211, 211);
            KalmanSexaquarkFinder_ChannelA();
        }
    }

    /*
        if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+") {
            ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
        }
        if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+,pi-,pi+") {
            ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
            ReconstructV0s_KF(idxPiMinusTracks, idxPosKaonTracks, 0, -211, 321, kfPionKaonPairs, idxPionKaonPairs);
            ReconstructV0s_KF(idxPiMinusTracks, idxPiPlusTracks, 0, -211, 211, kfPionPairs, idxPionPairs);
        }
        if (fReactionChannel == "AntiSexaquark,P->AntiProton,K+,K+,pi0") {
            ReconstructV0s_KF(idxPosKaonTracks, idxPosKaonTracks, 0, 312, 321, kfPosKaonPairs, idxPosKaonPairs);
        }
    */

    // FillTree();

    ClearContainers();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfTrees);
    PostData(2, fOutputListOfHists);
}

/*
 Initialize analysis task.
*/
void AliAnalysisTaskSexaquark::Initialize() {

    /* Parse TString: get reaction ID and injected sexaquark mass from SimulationSet */

    fReactionID = fSimulationSet(0);
    fInjectedSexaquarkMass = ((TString)fSimulationSet(1, fSimulationSet.Length())).Atof();

    /* PENDING: inj. sexaquark mass doesn't make sense when analyzing data... do smth about it... */

    getNegPdgCode_fromV0PdgCode[-3122] = -2212;
    getPosPdgCode_fromV0PdgCode[-3122] = 211;
    getNegPdgCode_fromV0PdgCode[310] = -211;
    getPosPdgCode_fromV0PdgCode[310] = 211;

    /* Assign the rest of properties related to the respective Reaction Channel */

    if (fReactionID == 'A') {
        fReactionChannel = "AntiSexaquark,N->AntiLambda,K0S";
        fProductsPDG = {-3122, 310};
        fFinalStateProductsPDG = {-2212, 211, -211, 211};
        fStruckNucleonPDG = 2112;
    } else if (fReactionID == 'D') {
        fReactionChannel = "AntiSexaquark,P->AntiLambda,K+";
        fProductsPDG = {-3122, 321};
        fFinalStateProductsPDG = {-2212, 211, 321};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'E') {
        fReactionChannel = "AntiSexaquark,P->AntiLambda,K+,pi-,pi+";
        fProductsPDG = {-3122, 321, -211, 211};
        fFinalStateProductsPDG = {-2212, 211, 321, -211, 211};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'H') {
        fReactionChannel = "AntiSexaquark,P->AntiProton,K+,K+,pi0";
        fProductsPDG = {-2212, 321, 321, 111};
        fFinalStateProductsPDG = {-2212, 321, 321, 11, -11, 11, -11, 11, -11, 11, -11};  // not used
        fStruckNucleonPDG = 2212;
    }

    AliInfo("Initializing AliAnalysisTaskSexaquark...");
    AliInfo("INPUT OPTIONS:");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);
    AliInfoF(">> Source of V0s       = %s", fSourceOfV0s.Data());
    AliInfoF(">> Simulation Set      = %s", fSimulationSet.Data());
    AliInfoF(">> Reaction ID         = %c", fReactionID);
    AliInfoF(">> Reaction Channel    = %s", fReactionChannel.Data());
    AliInfoF(">> Inj. Sexaquark Mass = %.2f", fInjectedSexaquarkMass);
}

/*
 Define track selection cuts.
 - Input: `cuts_option`
*/
void AliAnalysisTaskSexaquark::DefineTracksCuts(TString cuts_option) {

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
 Define V0 selection cuts.
 - Uses: `fSourceOfV0s`
 - Input: `cuts_option`
*/
void AliAnalysisTaskSexaquark::DefineV0Cuts(TString cuts_option) {

    /* Common */

    kMin_V0_Mass[-3122] = 1.08;
    kMax_V0_Mass[-3122] = 1.16;
    kMin_V0_Pt[-3122] = 1.0;

    kMin_V0_Mass[310] = 0.4;
    kMax_V0_Mass[310] = 0.6;
    kMin_V0_Pt[310] = 1.0;

    if (fSourceOfV0s == "offline") {

        // none so far
    }

    if (fSourceOfV0s == "on-the-fly") {

        kMax_V0_ImprvDCAbtwDau[-3122] = 0.6;

        kMax_V0_ImprvDCAbtwDau[310] = 0.6;
    }

    if (fSourceOfV0s == "custom") {

        kMax_V0_CPAwrtPV[-3122] = 0.99;
        kMin_V0_DCAwrtPV[-3122] = 4.;
        kMax_V0_Chi2[-3122] = 0.3;

        kMax_V0_Eta[310] = 0.6;
        kMax_V0_CPAwrtPV[310] = 0.99;
        kMax_V0_DCAposV0[310] = 0.15;
        kMin_V0_DCAwrtPV[310] = 4.;
        kMax_V0_Chi2[310] = 0.3;
        kMax_V0_ImprvDCAbtwDau[310] = 0.3;
    }

    if (fSourceOfV0s == "kalman") {

        kMax_V0_Eta[-3122] = 0.9;
        kMin_V0_Radius[-3122] = 20.;
        kMin_V0_DecayLength[-3122] = 40.;
        kMin_V0_CPAwrtPV[-3122] = 0.1;
        kMax_V0_CPAwrtPV[-3122] = 0.99;
        kMin_V0_DCAwrtPV[-3122] = 4.;
        kMax_V0_DCAbtwDau[-3122] = 2.;
        kMax_V0_DCAnegV0[-3122] = 2.;
        kMax_V0_DCAposV0[-3122] = 2.;
        kMax_V0_ArmPtOverAlpha[-3122] = 0.2;
        kMax_V0_Chi2ndf[-3122] = 0.4;

        kMax_V0_Eta[310] = 0.8;
        kMin_V0_Radius[310] = 15.;
        kMin_V0_DecayLength[310] = 30.;
        kMax_V0_DecayLength[310] = 175.;
        kMin_V0_CPAwrtPV[310] = 0.2;
        kMax_V0_CPAwrtPV[310] = 0.99;
        kMin_V0_DCAwrtPV[310] = 4.;
        kMax_V0_DCAbtwDau[310] = 0.25;
        kMax_V0_DCAnegV0[310] = 0.25;
        kMax_V0_DCAposV0[310] = 0.25;
        kMax_V0_Chi2ndf[310] = 0.4;
    }
}

/*
 Define anti-sexaquark candidates selection cuts.
 - Uses: `fSourceOfV0s`
 - Input: `cuts_option`
*/
void AliAnalysisTaskSexaquark::DefineSexaquarkCuts(TString cuts_option) {

    /* Common */

    // none so far

    if (fSourceOfV0s == "offline") {

        // none so far
    }

    if (fSourceOfV0s == "on-the-fly") {

        // none so far
    }

    if (fSourceOfV0s == "custom") {

        kMin_Sexa_Radius = 40;
        kMin_Sexa_CPAwrtPV = 0.99;
        kMax_Sexa_CPAwrtPV = 1.;
        kMax_Sexa_DCAv0aSV = 0.1;
        kMax_Sexa_DCAv0bSV = 0.1;
    }

    if (fSourceOfV0s == "kalman") {

        kMin_Sexa_Radius = 40;
        kMin_Sexa_CPAwrtPV = 0.99;
        kMax_Sexa_CPAwrtPV = 1.;
        kMax_Sexa_DCAbtwV0s = 2.;
        kMax_Sexa_DCAv0aSV = 2.;
        kMax_Sexa_DCAv0bSV = 2.;
        kMax_Sexa_DCAv0anegSV = 10.;
        kMax_Sexa_DCAv0aposSV = 10.;
        kMax_Sexa_DCAv0bnegSV = 10.;
        kMax_Sexa_DCAv0bposSV = 10.;
        kMax_Sexa_Chi2ndf = 0.5;
    }
}

/*
 Search for possible contradictions in the input options
*/
void AliAnalysisTaskSexaquark::CheckForInputErrors() {
    std::array<std::string, 4> validSources = {"offline", "on-the-fly", "custom", "kalman"};
    if (std::find(validSources.begin(), validSources.end(), fSourceOfV0s) == validSources.end()) {
        AliFatal("ERROR: source of V0s must be \"offline\", \"on-the-fly\", \"custom\" or \"kalman\".");
    }
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

/*
 Loop over MC particles in a single event. Store the indices of the signal particles.
 - Uses: `fMC`, `fPDG`, `fMC_PrimaryVertex`
*/
void AliAnalysisTaskSexaquark::ProcessMCGen() {

    AliMCParticle* mcPart;
    Int_t pdg_mc;

    AliMCParticle* mcNegDau;
    AliMCParticle* mcPosDau;

    Int_t n_daughters = 0;
    Int_t mcIdxNegDaughter = 0;
    Int_t mcIdxPosDaughter = 0;

    TVector3 originVertex;
    TVector3 decayVertex;

    Float_t pt;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;
    Float_t armenteros_alpha;
    Float_t armenteros_qt;
    Float_t opening_angle;

    /* Loop over MC gen. particles in a single event */

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);
        pdg_mc = mcPart->PdgCode();
        getPdgCode_fromMcIdx[mcIdx] = pdg_mc;

        n_daughters = mcPart->GetNDaughters();
        mcIdxNegDaughter = 0;
        mcIdxPosDaughter = 0;

        /* Protection in case of a particle with no momentum */

        if (mcPart->P() < 1E-6) {
            continue;
        }

        /* Get particle's properties */

        originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
        decayVertex.SetXYZ(0., 0., 0.);

        pt = mcPart->Pt();
        if (n_daughters) {
            GetDaughtersInfo(mcPart, mcIdxNegDaughter, mcIdxPosDaughter, decayVertex);
            cpa_wrt_pv = CosinePointingAngle(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                             fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                      fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        }

        /* Determine if particle is signal -- by checking if it's a reaction product */
        // NOTE:
        // If they're signal, their MC status number should be among 600 and 619.
        // This is done exclusively for the MC productions to the sexaquark analysis.
        // This was done to recognize each interaction, assigning a different number to each of
        // the 20 antisexaquark-nucleon reactions that were injected in a single event.
        //

        if (!isMcIdxSignal[mcIdx]) {  // prevent overwriting the previously assigned signal status
            isMcIdxSignal[mcIdx] = mcPart->MCStatusCode() >= 600 && mcPart->MCStatusCode() < 620;
            if (isMcIdxSignal[mcIdx]) {
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }
        }

        /* Determine if it's a secondary particle */
        // NOTE:
        // Include signal particles, as well, because despite they were injected as primaries,
        // in reality, they are not
        //

        isMcIdxSecondary[mcIdx] = mcPart->IsSecondaryFromMaterial() || mcPart->IsSecondaryFromWeakDecay() || isMcIdxSignal[mcIdx];

        /** Charged particles: anti-protons, kaons, pions **/

        if (pdg_mc == -2212 || pdg_mc == 321 || pdg_mc == 211 || pdg_mc == -211) {

            /* Fill histograms */

            fHist_Tracks_Bookkeep[pdg_mc]->Fill(0);
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());

            if (isMcIdxSecondary[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(1);
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }

            if (isMcIdxSignal[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(2);
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }
        }

        /** Neutral particles: anti-lambdas and K0S **/

        if (pdg_mc == -3122 || pdg_mc == 310) {

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdg_mc]->Fill(0);

            if (isMcIdxSecondary[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(2);

            if (isMcIdxSignal[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(4);

            /* Check if this V0 can be reconstructed */
            // via its decay channel containing two oppositely charged particles
            // -- anti-lambdas -> anti-proton piplus
            // -- K0S -> piplus piminus
            // (these conditions are defined in `GetDaughtersInfo()`)

            if (!n_daughters || !mcIdxNegDaughter || !mcIdxPosDaughter) continue;

            mcIndicesOfTrueV0s.push_back(mcIdx);

            /* Store mother-daughter info */

            doesMcIdxHaveMother[mcIdxNegDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxNegDaughter] = mcIdx;
            getNegDauMcIdx_fromMcIdx[mcIdx] = mcIdxNegDaughter;
            isMcIdxSignal[mcIdxNegDaughter] = isMcIdxSignal[mcIdx];
            if (isMcIdxSignal[mcIdxNegDaughter]) {
                getReactionIdx_fromMcIdx[mcIdxNegDaughter] = getReactionIdx_fromMcIdx[mcIdx];
                getMcIdx_fromReactionIdx[getReactionIdx_fromMcIdx[mcIdxNegDaughter]].push_back(mcIdxNegDaughter);
            }

            doesMcIdxHaveMother[mcIdxPosDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxPosDaughter] = mcIdx;
            getPosDauMcIdx_fromMcIdx[mcIdx] = mcIdxPosDaughter;
            isMcIdxSignal[mcIdxPosDaughter] = isMcIdxSignal[mcIdx];
            if (isMcIdxSignal[mcIdxPosDaughter]) {
                getReactionIdx_fromMcIdx[mcIdxPosDaughter] = getReactionIdx_fromMcIdx[mcIdx];
                getMcIdx_fromReactionIdx[getReactionIdx_fromMcIdx[mcIdxPosDaughter]].push_back(mcIdxPosDaughter);
            }

            /* Fill histograms */

            mcNegDau = (AliMCParticle*)fMC->GetTrack(mcIdxNegDaughter);
            mcPosDau = (AliMCParticle*)fMC->GetTrack(mcIdxPosDaughter);

            armenteros_qt = ArmenterosQt(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            armenteros_alpha = ArmenterosAlpha(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz(),
                                               mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());
            TVector3 neg(mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            TVector3 pos(mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());
            opening_angle = neg.Angle(pos);

            fHist_V0s_Bookkeep[pdg_mc]->Fill(1);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Mass")]->Fill(mcPart->M());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Radius")]->Fill(decayVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DecayLength")]->Fill((originVertex - decayVertex).Mag());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Zv")]->Fill(decayVertex.Z());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "OpeningAngle")]->Fill(opening_angle);

            if (isMcIdxSecondary[mcIdx]) {
                fHist_V0s_Bookkeep[pdg_mc]->Fill(3);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Radius")]->Fill(decayVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DecayLength")]->Fill((originVertex - decayVertex).Mag());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            }

            if (isMcIdxSignal[mcIdx]) {
                fHist_V0s_Bookkeep[pdg_mc]->Fill(5);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Radius")]->Fill(decayVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayLength")]->Fill((originVertex - decayVertex).Mag());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            }
        }  // end of anti-lambda && K0s condition
    }      // end of loop over MC particles
}

/*
 Loop over the daughters of a MC particle.
 If the particle is a relevant V0 (anti-lambda or K0S) that decays into charged particles,
 store the indices of its daughters.
 - Input: `mcPart`
 - Output: `mcIdxNegDaughter`, `mcIdxPosDaughter`, `decayVertex`
*/
void AliAnalysisTaskSexaquark::GetDaughtersInfo(AliMCParticle* mcPart, Int_t& mcIdxNegDaughter, Int_t& mcIdxPosDaughter, TVector3& decayVertex) {

    Int_t pdg_mc = mcPart->PdgCode();

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

        mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
        pdg_dau = mcDaughter->PdgCode();

        if (pdg_mc == -3122) {
            if (pdg_dau == -2212) mcIdxNegDaughter = mcIdxDaughter;
            if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
        }

        if (pdg_mc == 310) {
            if (pdg_dau == -211) mcIdxNegDaughter = mcIdxDaughter;
            if (pdg_dau == 211) mcIdxPosDaughter = mcIdxDaughter;
        }

        /* In any case, get the decay vertex of the V0, which is the origin vertex of its daughters */

        if (decayVertex.Mag() < 1E-6) decayVertex.SetXYZ(mcDaughter->Xv(), mcDaughter->Yv(), mcDaughter->Zv());
    }
}

/*
 Loop over generated antisexaquark-nucleon reactions.
 - Uses: `fMC`, `fPDG`, `fStruckNucleonPDG`
 - Containers: `getMcIdx_fromReactionIdx`
*/
void AliAnalysisTaskSexaquark::ProcessSignalInteractions() {

    AliMCParticle* mcProduct;

    /* Declare TLorentzVectors */

    TLorentzVector lvProduct;
    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Declare variables */

    TVector3 secondary_vertex(0., 0., 0.);
    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over reactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection in case of general-purpose simulations */

        if (!getMcIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Reset vectors */

        secondary_vertex.SetXYZ(0., 0., 0.);
        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /* Loop over reaction products */

        for (Int_t& mcIdx : getMcIdx_fromReactionIdx[reactionIdx]) {

            mcProduct = (AliMCParticle*)fMC->GetTrack(mcIdx);

            /* Get first generation products, the direct daughters of the reaction */

            if (mcProduct->MCStatusCode() != reactionIdx) continue;

            lvProduct.SetPxPyPzE(mcProduct->Px(), mcProduct->Py(), mcProduct->Pz(), mcProduct->E());

            lvAntiSexaquark = lvAntiSexaquark + lvProduct;

            /* Get sec. vertex */

            if (secondary_vertex.Mag() < 1E-6) secondary_vertex.SetXYZ(mcProduct->Xv(), mcProduct->Yv(), mcProduct->Zv());
        }

        lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest
        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /* Get properties */

        radius = TMath::Sqrt(TMath::Power(mcProduct->Xv(), 2) + TMath::Power(mcProduct->Yv(), 2));
        decay_length =
            TMath::Sqrt(TMath::Power(mcProduct->Xv() - fMC_PrimaryVertex->GetX(), 2) + TMath::Power(mcProduct->Yv() - fMC_PrimaryVertex->GetY(), 2) +
                        TMath::Power(mcProduct->Zv() - fMC_PrimaryVertex->GetZ(), 2));
        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, mcProduct->Xv(), mcProduct->Yv(), mcProduct->Zv(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), mcProduct->Xv(), mcProduct->Yv(), mcProduct->Zv(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill histograms */

        fHist_AntiSexaquarks_Bookkeep->Fill(0);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Radius")]->Fill(radius);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Zv")]->Fill(mcProduct->Zv());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("MCGen", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }  // end of loop over reactions
}

/*                */
/**    Tracks    **/
/*** ========== ***/

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisTaskSexaquark::ProcessTracks() {

    AliESDtrack* track;

    Int_t mcIdx;
    Int_t mcPdgCode;
    Int_t motherMcIdx;

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
        getMcIdx_fromEsdIdx[esdIdxTrack] = mcIdx;
        mcPdgCode = getPdgCode_fromMcIdx[mcIdx];

        /* Fill histograms before track selection */

        f2DHist_Tracks_TPCsignal["All"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        /* Apply track selection */

        doesEsdIdxPassTrackSelection[esdIdxTrack] = PassesTrackSelection(track);

        if (!doesEsdIdxPassTrackSelection[esdIdxTrack]) continue;

        getEsdIdx_fromMcIdx[mcIdx] = esdIdxTrack;  // purposefully done after `PassesTrackSelection()`

        /* Fill containers */

        isEsdIdxSignal[esdIdxTrack] = isMcIdxSignal[mcIdx];
        if (isEsdIdxSignal[esdIdxTrack]) {
            getReactionIdx_fromEsdIdx[esdIdxTrack] = getReactionIdx_fromMcIdx[mcIdx];
            getEsdIdx_fromReactionIdx[getReactionIdx_fromEsdIdx[esdIdxTrack]].push_back(esdIdxTrack);
        }

        doesEsdIdxHaveMother[esdIdxTrack] = doesMcIdxHaveMother[mcIdx];
        if (doesEsdIdxHaveMother[esdIdxTrack]) {
            motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
            getMotherMcIdx_fromEsdIdx[esdIdxTrack] = motherMcIdx;
            if (getNegDauMcIdx_fromMcIdx[motherMcIdx] == mcIdx) getNegDauEsdIdx_fromMcIdx[motherMcIdx] = esdIdxTrack;
            if (getPosDauMcIdx_fromMcIdx[motherMcIdx] == mcIdx) getPosDauEsdIdx_fromMcIdx[motherMcIdx] = esdIdxTrack;
        }

        /* Get track properties */

        pt = track->Pt();
        pz = track->Pz();
        eta = track->Eta();

        Float_t xy_impar_wrt_pv, z_impar_wrt_pv;
        track->GetImpactParameters(xy_impar_wrt_pv, z_impar_wrt_pv);  // pre-calculated DCA w.r.t. PV
        dca_wrt_pv = TMath::Sqrt(xy_impar_wrt_pv * xy_impar_wrt_pv + z_impar_wrt_pv * z_impar_wrt_pv);

        n_tpc_clusters = track->GetTPCNcls();  // NOTE: capital N
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

            if (mcPdgCode != esdPdgCode) continue;

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

            if (!isMcIdxSecondary[mcIdx]) continue;

            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pt")]->Fill(pt);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pz")]->Fill(pz);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Eta")]->Fill(eta);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Chi2/NClusters")]->Fill(chi2_over_nclusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "GoldenChi2")]->Fill(golden_chi2);
            PlotStatus(track, "Secondary", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Secondary"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!isEsdIdxSignal[esdIdxTrack]) continue;

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
        }  // end of loop over PID hypotheses
    }      // end of loop over tracks
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesTrackSelection(AliESDtrack* track) {

    /* Particle ID */

    Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

    std::vector<Int_t> possiblePID;
    if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton) possiblePID.push_back(2212);
    if (TMath::Abs(n_sigma_kaon) < kMax_NSigma_Kaon) possiblePID.push_back(321);
    if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion) possiblePID.push_back(211);
    if ((kMax_NSigma_Proton || kMax_NSigma_Kaon || kMax_NSigma_Pion) && !possiblePID.size()) return kFALSE;

    Double_t Eta = TMath::Abs(track->Eta());
    if (kMax_Track_Eta && Eta > kMax_Track_Eta) return kFALSE;

    Double_t NTPCClusters = track->GetTPCNcls();
    if (kMin_Track_NTPCClusters && NTPCClusters < kMin_Track_NTPCClusters) return kFALSE;

    Double_t Chi2PerNTPCClusters = track->GetTPCchi2() / (Double_t)track->GetTPCNcls();
    if (kMax_Track_Chi2PerNTPCClusters && Chi2PerNTPCClusters > kMax_Track_Chi2PerNTPCClusters) return kFALSE;

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
    Float_t DCAwrtPV = TMath::Sqrt(xy_impar * xy_impar + z_impar * z_impar);
    if (kMin_Track_DCAwrtPV && DCAwrtPV < kMin_Track_DCAwrtPV) return kFALSE;

    // >> reject low-pT pions
    // -- PENDING: this will also affect tracks with multiple PID hypotheses:
    // -- if a track both id as pion and proton, it will be cut by the highest pT cut
    for (Int_t& esdPdgCode : possiblePID) {
        if (kMin_Track_Pt[esdPdgCode] && track->Pt() < kMin_Track_Pt[esdPdgCode]) return kFALSE;
    }

    return kTRUE;
}

/*
 Plot the status of the track.
 - Input: `track`, `stage`, `esdPdgCode`
*/
void AliAnalysisTaskSexaquark::PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode) {

    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};

    for (Int_t i = 0; i < 16; i++) {
        if ((track->GetStatus() & StatusCollection[i])) fHist_Tracks[std::make_tuple("Found", set, esdPdgCode, "Status")]->Fill(i);
    }
}

/*                                                       */
/**  V0s and Anti-Sexaquarks -- That Require True Info  **/
/*** ================================================= ***/

/*
 Find all true V0s that both of their daughters passed track selection.
 Variables are true MC vertexing & kinematics.
 - Uses: `fMC`, `fMC_PrimaryVertex`, `fHist_V0s`, `fHist_V0s_Bookkeep`
*/
void AliAnalysisTaskSexaquark::ProcessFindableV0s() {

    AliMCParticle* mcV0;
    AliMCParticle* mcNegDaughter;
    AliMCParticle* mcPosDaughter;

    Bool_t is_signal;
    Bool_t is_secondary;

    Int_t pdg_code;
    Float_t radius;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;
    Float_t decay_length;
    Float_t arm_qt;
    Float_t arm_alpha;
    Float_t opening_angle;

    /* Loop over true V0s */

    for (Int_t& mcIdxOfTrueV0 : mcIndicesOfTrueV0s) {

        mcV0 = (AliMCParticle*)fMC->GetTrack(mcIdxOfTrueV0);

        /* Check if both daughters were reconstructed */

        if (!getNegDauEsdIdx_fromMcIdx[mcIdxOfTrueV0] || !getPosDauEsdIdx_fromMcIdx[mcIdxOfTrueV0]) {
            continue;
        }

        is_signal = mcV0->MCStatusCode() >= 600 && mcV0->MCStatusCode() < 620;
        is_secondary = mcV0->IsSecondaryFromMaterial() || mcV0->IsSecondaryFromWeakDecay() || is_signal;

        /* Get the V0 decay vertex from one of the true daughters */

        mcNegDaughter = (AliMCParticle*)fMC->GetTrack(getNegDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);
        mcPosDaughter = (AliMCParticle*)fMC->GetTrack(getPosDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);

        /* Get V0 properties */

        pdg_code = mcV0->PdgCode();
        radius = TMath::Sqrt(mcNegDaughter->Xv() * mcNegDaughter->Xv() + mcNegDaughter->Yv() * mcNegDaughter->Yv());
        cpa_wrt_pv = CosinePointingAngle(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                         fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        decay_length = TMath::Sqrt(TMath::Power(mcNegDaughter->Xv() - mcV0->Xv(), 2) + TMath::Power(mcNegDaughter->Yv() - mcV0->Yv(), 2) +
                                   TMath::Power(mcNegDaughter->Zv() - mcV0->Zv(), 2));
        arm_qt = ArmenterosQt(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        arm_alpha = ArmenterosAlpha(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz(),
                                    mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());
        TVector3 neg(mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        TVector3 pos(mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());
        opening_angle = neg.Angle(pos);

        /* Fill bookkeeping histograms */

        if (pdg_code == -3122 || pdg_code == 310) {

            fHist_V0s_Bookkeep[pdg_code]->Fill(6);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "OpeningAngle")]->Fill(opening_angle);

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(7);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "OpeningAngle")]->Fill(opening_angle);

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(8);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "OpeningAngle")]->Fill(opening_angle);
        }
    }  // end of loop over true V0s
}

/*
  Find the anti-sexaquarks that have all of their *final-state particles* reconstructed.
  - For example:
    -- Reaction Channel A: 1 anti-proton, 2 pi+, 1 pi-
    -- Reaction Channel D: 1 anti-proton, 1 pi+, 1 K+
    -- Reaction Channel E: 1 anti-proton, 2 pi+, 1 pi-, 1 K+
*/
void AliAnalysisTaskSexaquark::ProcessFindableSexaquarks() {

    AliMCParticle* mcFSPart;  // FS: final-state
    Int_t mcIdx;
    AliMCParticle* mcFGPart;  // FG: first gen.
    Int_t motherMcIdx;

    /* Declare TLorentzVectors */

    TLorentzVector lvFinalStateParticle;
    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Declare variables */

    TVector3 secondary_vertex(0., 0., 0.);
    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection in case of general-purpose simulations */

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Findable condition */
        // All final state products should have been reconstructed and passed track selection

        std::unordered_multiset<Int_t> pdg_reconstructed_particles;
        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {
            pdg_reconstructed_particles.insert(getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]]);
        }

        std::unordered_multiset<Int_t> expected_pdg_fs_particles(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());

        if (expected_pdg_fs_particles != pdg_reconstructed_particles) continue;

        /* Reset vectors */

        secondary_vertex.SetXYZ(0., 0., 0.);
        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /* Loop, to get the reconstructed final-state particles info */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            mcIdx = getMcIdx_fromEsdIdx[esdIdx];
            mcFSPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

            lvFinalStateParticle.SetPxPyPzE(mcFSPart->Px(), mcFSPart->Py(), mcFSPart->Pz(), mcFSPart->E());

            lvAntiSexaquark = lvAntiSexaquark + lvFinalStateParticle;

            /* Get sec. vertex from the origin of the first-generation products of the reaction */

            if (secondary_vertex.Mag() < 1E-6) {
                if (mcFSPart->MCStatusCode() == reactionIdx) {
                    secondary_vertex.SetXYZ(mcFSPart->Xv(), mcFSPart->Yv(), mcFSPart->Zv());
                } else if (doesMcIdxHaveMother[mcIdx]) {
                    motherMcIdx = getMotherMcIdx_fromMcIdx[mcIdx];
                    mcFGPart = (AliMCParticle*)fMC->GetTrack(motherMcIdx);
                    if (mcFGPart->MCStatusCode() == reactionIdx) secondary_vertex.SetXYZ(mcFGPart->Xv(), mcFGPart->Yv(), mcFGPart->Zv());
                }
            }
        }

        lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest
        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /* Get properties */

        radius = TMath::Sqrt(TMath::Power(secondary_vertex.X(), 2) + TMath::Power(secondary_vertex.Y(), 2));
        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), secondary_vertex.X(), secondary_vertex.Y(),
                                  secondary_vertex.Z(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill histograms  */

        fHist_AntiSexaquarks_Bookkeep->Fill(1);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Radius")]->Fill(radius);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Zv")]->Fill(secondary_vertex.Z());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }
}

/*                             */
/**  V0s -- ALICE V0 Finders  **/
/*** ======================= ***/

/*
 Loop over all the V0s that were found and stored in the ESD file by the Official V0 Finders.
 - Uses: `fESD`, `fMC`, `fPDG`, `fPrimaryVertex`
 - Containers: `esdAntiLambdas`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`, `isAntiLambdaIdxATrueV0`,
 `isEsdIdxDaughterOfTrueV0`, `getMcIdxOfTrueV0_fromEsdIdx`, `isMcIdxSignal`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `isKaonZeroShortIdxATrueV0`
*/
void AliAnalysisTaskSexaquark::GetV0sFromESD(Bool_t onTheFly) {

    AliESDv0* V0;

    Int_t esdIdxPos;
    Int_t esdIdxNeg;

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    Float_t radius;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv, dca_btw_dau;
    Float_t dca_neg_v0, dca_pos_v0;
    Float_t decay_length;
    Float_t opening_angle;
    Float_t arm_alpha, arm_qt;

    AliESDv0 outputV0;

    Double_t xthiss, xpp;
    Float_t impar_neg_v0[2], impar_pos_v0[2];

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Auxiliary containers */

    Bool_t is_true_v0;
    std::vector<AliESDv0>* esdFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    /* Loop over V0s that are stored in the ESD */

    for (Int_t esdIdxV0 = 0; esdIdxV0 < fESD->GetNumberOfV0s(); esdIdxV0++) {

        V0 = fESD->GetV0(esdIdxV0);

        /* Choose between offline and on-the-fly V0s */

        if (!onTheFly && V0->GetOnFlyStatus()) {
            continue;
        }

        if (onTheFly && !V0->GetOnFlyStatus()) {
            continue;
        }

        esdIdxNeg = V0->GetNindex();
        esdIdxPos = V0->GetPindex();

        /* Check if both daughters passed track selection */

        if (!doesEsdIdxPassTrackSelection[esdIdxNeg] || !doesEsdIdxPassTrackSelection[esdIdxPos]) {
            continue;
        }

        /* Get tracks */

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

        /* Apply V0 cuts */

        if (!PassesV0Cuts(-3122, V0, neg_track, pos_track) && !PassesV0Cuts(310, V0, neg_track, pos_track)) {
            continue;
        }

        /* Get V0 properties */

        V0->GetXYZ(V0_X, V0_Y, V0_Z);
        V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
        V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

        neg_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_neg_v0);
        pos_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_pos_v0);

        radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
        cpa_wrt_pv = V0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_wrt_pv =
            LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_btw_dau = TMath::Abs(neg_track->GetDCA(pos_track, fMagneticField, xthiss, xpp));              // DCAxyz
        dca_neg_v0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);  // DCAxyz
        dca_pos_v0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);  // DCAxyz
        decay_length = TMath::Sqrt((V0_X - fPrimaryVertex->GetX()) * (V0_X - fPrimaryVertex->GetX()) +
                                   (V0_Y - fPrimaryVertex->GetY()) * (V0_Y - fPrimaryVertex->GetY()) +
                                   (V0_Z - fPrimaryVertex->GetZ()) * (V0_Z - fPrimaryVertex->GetZ()));

        /* Fill histograms */

        std::vector<Int_t> possiblePID;
        if (PassesV0Cuts(-3122, V0, neg_track, pos_track)) possiblePID.push_back(-3122);
        if (PassesV0Cuts(310, V0, neg_track, pos_track)) possiblePID.push_back(310);

        for (Int_t& v0PdgCode : possiblePID) {

            /* Update TLorentzVectors to get reconstructed mass */

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            arm_qt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Determine containers to use */

            if (v0PdgCode == -3122) {
                esdFoundV0s = &esdAntiLambdas;
                getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
                getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
            } else if (v0PdgCode == 310) {
                esdFoundV0s = &esdKaonsZeroShort;
                getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
                getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
            }

            /* Store V0 */

            outputV0 = *V0;

            esdFoundV0s->push_back(outputV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxPos;

            is_true_v0 = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                         getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                         getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == v0PdgCode;

            /* Fill histograms */

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(30);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Eta")]->Fill(V0->Eta());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pt")]->Fill(V0->Pt());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pz")]->Fill(V0->Pz());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(arm_alpha, arm_qt);

            if (!is_true_v0) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(31);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Eta")]->Fill(V0->Eta());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pt")]->Fill(V0->Pt());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pz")]->Fill(V0->Pz());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(32);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Eta")]->Fill(V0->Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pt")]->Fill(V0->Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pz")]->Fill(V0->Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(33);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Eta")]->Fill(V0->Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pt")]->Fill(V0->Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pz")]->Fill(V0->Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(arm_alpha, arm_qt);
        }
    }  // end of loop over V0s
}

/*
 Custom V0 Finder.
 (Adapted from `AliRoot/STEER/ESD/AliV0vertexer.cxx`)
 - Uses: `fESD`, `fMC`, `fPDG`, `fPrimaryVertex`, `fMagneticField`
 - Input: `pdgV0`
*/
void AliAnalysisTaskSexaquark::CustomV0Finder(Int_t pdgV0) {

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Auxiliary variables */

    Double_t neg_d_wrt_pv, pos_d_wrt_pv;
    Double_t radius;
    Double_t decay_length;
    Double_t cpa_wrt_pv;
    Double_t dca_wrt_pv, dca_btw_dau;
    Double_t dca_neg_v0, dca_pos_v0;
    Double_t arm_qt, arm_alpha;
    Double_t opening_angle;

    Double_t xthiss, xpp;
    Float_t impar_neg_v0[2], impar_pos_v0[2];

    /* Custom V0 Finder's Technical Cuts */

    Double_t COV0F_DCAmax = 0.6;  // (d=1.) max DCA between the daughter tracks
    Double_t COV0F_Rmin = 50.;    // (d=0.9) min radius of the V0 radius
    Double_t COV0F_Rmax = 200.;   // (d=100) max radius of the V0 radius
    Bool_t COV0F_UseImprovedFinding = kTRUE;

    /*  Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    Bool_t is_true_v0;
    std::vector<AliESDv0>* esdFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        esdFoundV0s = &esdAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        esdFoundV0s = &esdKaonsZeroShort;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
    }

    /* Loop over all possible pairs of selected tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            /* Sanity check: just in case */

            if (esdIdxNeg == esdIdxPos) continue;

            /* Get tracks */

            neg_track = fESD->GetTrack(esdIdxNeg);
            pos_track = fESD->GetTrack(esdIdxPos);

            /* --------------------------------------------- */
            /* --- Beginning of Core of Custom V0 Finder --- */

            neg_d_wrt_pv = neg_track->GetD(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fMagneticField);
            if (TMath::Abs(neg_d_wrt_pv) > COV0F_Rmax) continue;
            pos_d_wrt_pv = pos_track->GetD(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fMagneticField);
            if (TMath::Abs(pos_d_wrt_pv) > COV0F_Rmax) continue;

            AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;

            Double_t xn, xp;
            if (COV0F_UseImprovedFinding && Preoptimize(ntp, ptp, &xn, &xp, fMagneticField)) {
                nt.PropagateTo(xn, fMagneticField);
                pt.PropagateTo(xp, fMagneticField);
            }

            Double_t dca_after_preopt = ntp->GetDCA(ptp, fMagneticField, xn, xp);
            if (dca_after_preopt > COV0F_DCAmax) continue;
            if ((xn + xp) > 2 * COV0F_Rmax) continue;
            if ((xn + xp) < 2 * COV0F_Rmin) continue;

            nt.PropagateTo(xn, fMagneticField);
            pt.PropagateTo(xp, fMagneticField);

            /* Note: declaration of AliESDv0 within the loop, because it has important functions on its constructor */
            AliESDv0 vertex(nt, esdIdxNeg, pt, esdIdxPos);
            if (COV0F_UseImprovedFinding) vertex.Refit();

            vertex.SetDcaV0Daughters(dca_after_preopt);

            /* --- End of Core of Custom V0 Finder --- */
            /* --------------------------------------- */

            /* Apply cuts */

            if (!PassesV0Cuts(pdgV0, &vertex, neg_track, pos_track)) continue;

            /* Get V0 properties */

            vertex.GetXYZ(V0_X, V0_Y, V0_Z);
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);

            cpa_wrt_pv = vertex.GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

            dca_wrt_pv = LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                      fPrimaryVertex->GetZ());

            dca_btw_dau = TMath::Abs(neg_track->GetDCA(pos_track, fMagneticField, xthiss, xpp));  // DCAxyz

            neg_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_neg_v0);
            dca_neg_v0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);  // DCAxyz

            pos_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_pos_v0);
            dca_pos_v0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);  // DCAxyz

            decay_length = TMath::Sqrt((V0_X - fPrimaryVertex->GetX()) * (V0_X - fPrimaryVertex->GetX()) +
                                       (V0_Y - fPrimaryVertex->GetY()) * (V0_Y - fPrimaryVertex->GetY()) +
                                       (V0_Z - fPrimaryVertex->GetZ()) * (V0_Z - fPrimaryVertex->GetZ()));

            arm_qt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Store V0 */

            esdFoundV0s->push_back(vertex);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxPos;

            is_true_v0 = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                         getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                         getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(30);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(vertex.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(vertex.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(vertex.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(arm_alpha, arm_qt);

            if (!is_true_v0) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(31);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(vertex.Eta());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(vertex.Pt());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(vertex.Pz());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(32);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(vertex.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(vertex.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(vertex.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(33);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(vertex.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(vertex.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(vertex.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(arm_alpha, arm_qt);
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*
 Apply cuts to a V0 candidate.
 - Input: `pdgV0`, `v0`, `neg_track`, `pos_track`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesV0Cuts(Int_t pdgV0, AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(10);

    Double_t V0_X, V0_Y, V0_Z;
    v0->GetXYZ(V0_X, V0_Y, V0_Z);

    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;
    v0->GetNPxPyPz(N_Px, N_Py, N_Pz);
    v0->GetPPxPyPz(P_Px, P_Py, P_Pz);
    TLorentzVector lvTrackNeg, lvTrackPos;
    lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[pdgV0])->Mass());
    lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[pdgV0])->Mass());
    TLorentzVector lvV0 = lvTrackNeg + lvTrackPos;

    Double_t Radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(11);

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(12);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(13);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(14);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(v0->Xv() - fPrimaryVertex->GetX(), 2) + TMath::Power(v0->Yv() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(v0->Zv() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(15);

    Double_t CPAwrtPV = v0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(16);

    Double_t DCAwrtPV = LinePointDCA(v0->Px(), v0->Py(), v0->Pz(), v0->Xv(), v0->Yv(), v0->Zv(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                     fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(17);

    Float_t impar_neg_v0[2];
    neg_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_neg_v0);
    Double_t DCAnegV0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(18);

    Float_t impar_pos_v0[2];
    pos_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_pos_v0);
    Double_t DCAposV0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(19);

    Double_t ArmPt = ArmenterosQt(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz);
    Double_t ArmAlpha = ArmenterosAlpha(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz, P_Px, P_Py, P_Pz);
    Double_t ArmPtOverAlpha = ArmPt / TMath::Abs(ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(20);

    /* Exclusive for Offline, On-The-Fly, and Custom V0 Finders */

    Double_t Chi2 = v0->GetChi2V0();
    if (kMax_V0_Chi2[pdgV0] && Chi2 > kMax_V0_Chi2[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(21);

    Double_t ImprvDCAbtwDau = v0->GetDcaV0Daughters();
    if (kMax_V0_ImprvDCAbtwDau[pdgV0] && ImprvDCAbtwDau > kMax_V0_ImprvDCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(22);

    return kTRUE;
}

/*
 Apply cuts to a V0 candidate.
 - Input: `pdgV0`, `kfV0`, `kfDaughterNeg`, `kfDaughterPos`, `lvV0`, `lvTrackNeg`, `lvTrackPos`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesV0Cuts(Int_t pdgV0, KFParticleMother kfV0, KFParticle kfDaughterNeg, KFParticle kfDaughterPos,
                                              TLorentzVector lvV0, TLorentzVector lvTrackNeg, TLorentzVector lvTrackPos) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(10);
    Double_t Radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(11);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(12);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(13);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(kfV0.GetX() - fPrimaryVertex->GetX(), 2) + TMath::Power(kfV0.GetY() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(kfV0.GetZ() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(14);

    Double_t CPAwrtPV =
        CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(15);

    Double_t DCAwrtPV = LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(),
                                     fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(16);

    Double_t DCAbtwDau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(17);

    Double_t DCAnegV0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(18);

    Double_t DCAposV0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(19);

    Double_t ArmPt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
    Double_t ArmAlpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
    Double_t ArmPtOverAlpha = TMath::Abs(ArmPt / ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(20);

    Double_t Chi2ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();
    if (kMax_V0_Chi2ndf[pdgV0] && Chi2ndf > kMax_V0_Chi2ndf[pdgV0]) return kFALSE;  // exclusive for KF method
    fHist_V0s_Bookkeep[pdgV0]->Fill(21);

    // there's no improved DCA btw daughters cut
    fHist_V0s_Bookkeep[pdgV0]->Fill(22);

    return kTRUE;
}

/*                                     */
/**  Sexaquark -- Geometrical Finder  **/
/*** =============================== ***/

/*
 Using two AliESDv0 objects, reconstruct the anti-sexaquark candidate.
 - Containers: `esdAntiLambdas`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`
*/
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelA() {

    /* Declare AliRoot objects */

    AliESDv0 AntiLambda;
    AliESDv0 KaonZeroShort;

    AliESDtrack* AntiLambda_NegDau;
    AliESDtrack* AntiLambda_PosDau;
    AliESDtrack* KaonZeroShort_NegDau;
    AliESDtrack* KaonZeroShort_PosDau;

    /* Declare TVector3s  */

    TVector3 AntiLambda_momentum, KaonZeroShort_momentum;
    TVector3 AntiLambda_vertex, KaonZeroShort_vertex;

    TVector3 SecondaryVertex;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Declare auxiliar variables */

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    TVector3 pca1, pca2;
    Float_t radius;
    Float_t decay_length;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;
    Float_t dca_btw_v0s;
    Float_t dca_v0a_sv;
    Float_t dca_v0b_sv;
    Float_t dca_v0a_neg_sv;
    Float_t dca_v0a_pos_sv;
    Float_t dca_v0b_neg_sv;
    Float_t dca_v0b_pos_sv;

    Bool_t is_signal;

    Float_t impar_v0a_neg_sv[2], impar_v0a_pos_sv[2];
    Float_t impar_v0b_neg_sv[2], impar_v0b_pos_sv[2];

    /* Loop over all possible pairs of V0s */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)esdAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)esdKaonsZeroShort.size(); idxKaonZeroShort++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Get AliRoot objects */

            AntiLambda = esdAntiLambdas[idxAntiLambda];
            KaonZeroShort = esdKaonsZeroShort[idxKaonZeroShort];

            AntiLambda_NegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            AntiLambda_PosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));
            KaonZeroShort_NegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_KaonZeroShort));
            KaonZeroShort_PosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_KaonZeroShort));

            /* Transport V0s as if they were straight lines */

            AntiLambda_momentum.SetXYZ(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz());
            KaonZeroShort_momentum.SetXYZ(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz());

            AntiLambda_vertex.SetXYZ(AntiLambda.Xv(), AntiLambda.Yv(), AntiLambda.Zv());
            KaonZeroShort_vertex.SetXYZ(KaonZeroShort.Xv(), KaonZeroShort.Yv(), KaonZeroShort.Zv());

            dca_btw_v0s = Calculate_TwoLinesDCA_v2(AntiLambda_vertex, AntiLambda_momentum, KaonZeroShort_vertex, KaonZeroShort_momentum, pca1, pca2);

            SecondaryVertex = 0.5 * (pca1 + pca2);  // middle point between the two points of closest approach

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda.SetPxPyPzE(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz(), AntiLambda.E());
            lvKaonZeroShort.SetPxPyPzE(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz(), KaonZeroShort.E());

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;

            /* Apply cuts */

            if (!PassesSexaquarkCuts_ChannelA(SecondaryVertex, lvAntiSexaquark, AntiLambda, KaonZeroShort, AntiLambda_NegDau, AntiLambda_PosDau,
                                              KaonZeroShort_NegDau, KaonZeroShort_PosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(SecondaryVertex.X() * SecondaryVertex.X() + SecondaryVertex.Y() * SecondaryVertex.Y());
            decay_length = TMath::Sqrt(TMath::Power(SecondaryVertex.X() - fPrimaryVertex->GetX(), 2) +
                                       TMath::Power(SecondaryVertex.Y() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(SecondaryVertex.Z() - fPrimaryVertex->GetZ(), 2));
            cpa_wrt_pv =
                CosinePointingAngle(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                    SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                      SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_v0a_sv = LinePointDCA(lvAntiLambda.Px(), lvAntiLambda.Py(), lvAntiLambda.Pz(), AntiLambda_vertex.X(), AntiLambda_vertex.Y(),
                                      AntiLambda_vertex.Z(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
            dca_v0b_sv =
                LinePointDCA(lvKaonZeroShort.Px(), lvKaonZeroShort.Py(), lvKaonZeroShort.Pz(), KaonZeroShort_vertex.X(), KaonZeroShort_vertex.Y(),
                             KaonZeroShort_vertex.Z(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());

            AntiLambda_NegDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_neg_sv);
            dca_v0a_neg_sv = TMath::Sqrt(impar_v0a_neg_sv[0] * impar_v0a_neg_sv[0] + impar_v0a_neg_sv[1] * impar_v0a_neg_sv[1]);

            AntiLambda_PosDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_pos_sv);
            dca_v0a_pos_sv = TMath::Sqrt(impar_v0a_pos_sv[0] * impar_v0a_pos_sv[0] + impar_v0a_pos_sv[1] * impar_v0a_pos_sv[1]);

            KaonZeroShort_NegDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_neg_sv);
            dca_v0b_neg_sv = TMath::Sqrt(impar_v0b_neg_sv[0] * impar_v0b_neg_sv[0] + impar_v0b_neg_sv[1] * impar_v0b_neg_sv[1]);

            KaonZeroShort_PosDau->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_pos_sv);
            dca_v0b_pos_sv = TMath::Sqrt(impar_v0b_pos_sv[0] * impar_v0b_pos_sv[0] + impar_v0b_pos_sv[1] * impar_v0b_pos_sv[1]);

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(19);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(SecondaryVertex.Z());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(20);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(SecondaryVertex.Z());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv);

        }  // end of loop over K0S
    }      // end of loop over anti-lambdas

    return;
}

/*
 Apply anti-sexaquark candidate cuts.
 - Use: `fPrimaryVertex`, `fMagneticField`
 - Input: `SecondaryVertex`, `lvAntiSexaquark`, `esdAntiLambda`, `esdKaonZeroShort`, `esdAntiLambdaNeg`, `esdAntiLambdaPos`, `esdKaonZeroShortNeg`,
 `esdKaonZeroShortPos`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark, AliESDv0 esdAntiLambda,
                                                              AliESDv0 esdKaonZeroShort, AliESDtrack* esdAntiLambdaNeg, AliESDtrack* esdAntiLambdaPos,
                                                              AliESDtrack* esdKaonZeroShortNeg, AliESDtrack* esdKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(2);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(3);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(4);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(5);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(6);

    Double_t radius = TMath::Sqrt(SecondaryVertex.X() * SecondaryVertex.X() + SecondaryVertex.Y() * SecondaryVertex.Y());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(7);

    Double_t decay_length = TMath::Sqrt((SecondaryVertex.X() - fPrimaryVertex->GetX()) * (SecondaryVertex.X() - fPrimaryVertex->GetX()) +
                                        (SecondaryVertex.Y() - fPrimaryVertex->GetY()) * (SecondaryVertex.Y() - fPrimaryVertex->GetY()) +
                                        (SecondaryVertex.Z() - fPrimaryVertex->GetZ()) * (SecondaryVertex.Z() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(8);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fPrimaryVertex->GetX(),
                                              fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(9);

    Double_t dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                       SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(10);

    TVector3 AntiLambda_vertex(esdAntiLambda.Xv(), esdAntiLambda.Yv(), esdAntiLambda.Zv());
    TVector3 AntiLambda_momentum(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz());
    TVector3 KaonZeroShort_vertex(esdKaonZeroShort.Xv(), esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv());
    TVector3 KaonZeroShort_momentum(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz());
    TVector3 pca1, pca2;
    Double_t dca_btw_v0s = Calculate_TwoLinesDCA_v2(AntiLambda_vertex, AntiLambda_momentum, KaonZeroShort_vertex, KaonZeroShort_momentum, pca1, pca2);
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(11);

    Double_t dca_v0a_sv = LinePointDCA(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz(), esdAntiLambda.Xv(), esdAntiLambda.Yv(),
                                       esdAntiLambda.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(12);

    Double_t dca_v0b_sv = LinePointDCA(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz(), esdKaonZeroShort.Xv(),
                                       esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(13);

    Float_t impar_v0a_neg_sv[2];
    esdAntiLambdaNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_neg_sv);
    Double_t dca_v0a_neg_sv = TMath::Sqrt(impar_v0a_neg_sv[0] * impar_v0a_neg_sv[0] + impar_v0a_neg_sv[1] * impar_v0a_neg_sv[1]);
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(14);

    Float_t impar_v0a_pos_sv[2];
    esdAntiLambdaPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_pos_sv);
    Double_t dca_v0a_pos_sv = TMath::Sqrt(impar_v0a_pos_sv[0] * impar_v0a_pos_sv[0] + impar_v0a_pos_sv[1] * impar_v0a_pos_sv[1]);
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(15);

    Float_t impar_v0b_neg_sv[2];
    esdKaonZeroShortNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_neg_sv);
    Double_t dca_v0b_neg_sv = TMath::Sqrt(impar_v0b_neg_sv[0] * impar_v0b_neg_sv[0] + impar_v0b_neg_sv[1] * impar_v0b_neg_sv[1]);
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(16);

    Float_t impar_v0b_pos_sv[2];
    esdKaonZeroShortPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_pos_sv);
    Double_t dca_v0b_pos_sv = TMath::Sqrt(impar_v0b_pos_sv[0] * impar_v0b_pos_sv[0] + impar_v0b_pos_sv[1] * impar_v0b_pos_sv[1]);
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(17);

    // since there's no Chi2/ndf cut
    fHist_AntiSexaquarks_Bookkeep->Fill(18);

    return kTRUE;
}

/*                          */
/**  V0s -- Kalman Filter  **/
/*** ==================== ***/

/*
 Find all V0s via Kalman Filter.
 - Uses: `fESD`, `fMC`, `fMagneticField`
 - Input: `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 - Output: `kfV0s`, `idxDaughters`
*/
void AliAnalysisTaskSexaquark::KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {

    AliESDtrack* esdTrackNeg;
    AliESDtrack* esdTrackPos;

    Bool_t successful_daughters_check;  // PENDING: I still don't understand this...

    Double_t impar_neg[2], impar_pos[2];

    Double_t radius;
    Double_t cpa_wrt_pv;
    Double_t decay_length;
    Double_t dca_wrt_pv;
    Double_t dca_btw_dau;
    Double_t dca_neg_v0;
    Double_t dca_pos_v0;
    Double_t arm_qt;
    Double_t arm_alpha;
    Double_t opening_angle;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    Bool_t is_true_v0;
    std::vector<KFParticleMother>* kfFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;

        kfFoundV0s = &kfAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;

        kfFoundV0s = &kfKaonsZeroShort;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
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
            dca_wrt_pv = TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_dau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
            dca_neg_v0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
            dca_pos_v0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));
            decay_length = TMath::Sqrt((kfV0.GetX() - fPrimaryVertex->GetX()) * (kfV0.GetX() - fPrimaryVertex->GetX()) +
                                       (kfV0.GetY() - fPrimaryVertex->GetY()) * (kfV0.GetY() - fPrimaryVertex->GetY()) +
                                       (kfV0.GetZ() - fPrimaryVertex->GetZ()) * (kfV0.GetZ() - fPrimaryVertex->GetZ()));
            arm_qt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Store V0 */

            kfFoundV0s->push_back(kfV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxPos;

            is_true_v0 = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                         getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                         getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(30);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
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

            if (!is_true_v0) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(31);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
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

            if (!isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(32);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2")]->Fill(kfV0.GetChi2());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2ndf")]->Fill((Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(arm_alpha, arm_qt);

            if (!isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]]) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(33);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Radius")]->Fill(radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
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
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel A.
 - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`, `fMagneticField`
 - Containers: `kf[AntiLambdas|KaonsZeroShort]`, `getEsdIdxOf[Neg|Pos]Dau_from[AntiLambda|KaonZeroShort]Idx`
*/
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelA() {

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack *esdKaonZeroShortNegDau, *esdKaonZeroShortPosDau;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;
    TLorentzVector lvKaonZeroShortNegDau, lvKaonZeroShortPosDau;

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;
    Float_t dca_btw_v0s;
    Float_t dca_v0a_sv;
    Float_t dca_v0b_sv;
    Float_t dca_v0a_neg_sv;
    Float_t dca_v0a_pos_sv;
    Float_t dca_v0b_neg_sv;
    Float_t dca_v0b_pos_sv;

    Bool_t is_signal;

    Float_t impar_v0a_neg_sv[2], impar_v0a_pos_sv[2];
    Float_t impar_v0b_neg_sv[2], impar_v0b_pos_sv[2];

    /* Loop over all possible pairs of anti-lambdas and K0S */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)kfKaonsZeroShort.size(); idxKaonZeroShort++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Kalman Filter */

            KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];
            KFParticleMother kfKaonZeroShort = kfKaonsZeroShort[idxKaonZeroShort];

            KFParticleMother kfAntiSexaquark;
            kfAntiSexaquark.AddDaughter(kfAntiLambda);
            kfAntiSexaquark.AddDaughter(kfKaonZeroShort);
            kfAntiSexaquark.TransportToDecayVertex();

            /* Transport V0s to the secondary vertex */

            kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
            kfKaonZeroShort.SetProductionVertex(kfAntiSexaquark);

            kfAntiLambda.TransportToProductionVertex();
            kfKaonZeroShort.TransportToProductionVertex();

            /* Get and transport the final-state AliESDtracks */

            esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));
            esdKaonZeroShortNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_KaonZeroShort));
            esdKaonZeroShortPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_KaonZeroShort));

            KFParticle kfAntiLambdaNegDau =
                CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
            KFParticle kfAntiLambdaPosDau =
                CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());
            KFParticle kfKaonZeroShortNegDau = CreateKFParticle(*esdKaonZeroShortNegDau, fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass(),
                                                                (Int_t)esdKaonZeroShortNegDau->Charge());
            KFParticle kfKaonZeroShortPosDau = CreateKFParticle(*esdKaonZeroShortPosDau, fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass(),
                                                                (Int_t)esdKaonZeroShortPosDau->Charge());

            kfAntiLambdaNegDau = TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                                     (Int_t)esdAntiLambdaNegDau->Charge());
            kfAntiLambdaPosDau = TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                                     (Int_t)esdAntiLambdaPosDau->Charge());
            kfKaonZeroShortNegDau =
                TransportKFParticle(kfKaonZeroShortNegDau, kfKaonZeroShortPosDau, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(),
                                    (Int_t)esdKaonZeroShortNegDau->Charge());
            kfKaonZeroShortPosDau =
                TransportKFParticle(kfKaonZeroShortPosDau, kfKaonZeroShortNegDau, fPDG.GetParticle(fFinalStateProductsPDG[3])->PdgCode(),
                                    (Int_t)esdKaonZeroShortPosDau->Charge());

            lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
            lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
            lvKaonZeroShortNegDau.SetXYZM(kfKaonZeroShortNegDau.Px(), kfKaonZeroShortNegDau.Py(), kfKaonZeroShortNegDau.Pz(),
                                          fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());
            lvKaonZeroShortPosDau.SetXYZM(kfKaonZeroShortPosDau.Px(), kfKaonZeroShortPosDau.Py(), kfKaonZeroShortPosDau.Pz(),
                                          fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass());

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;
            lvKaonZeroShort = lvKaonZeroShortNegDau + lvKaonZeroShortPosDau;

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;

            /* Apply cuts */

            if (!PassesSexaquarkCuts_ChannelA(kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfKaonZeroShort, kfAntiLambdaNegDau, kfAntiLambdaPosDau,
                                              kfKaonZeroShortNegDau, kfKaonZeroShortPosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
            decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                       (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                       (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),  //
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
            dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));

            dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPosDau.GetDistanceFromVertex(kfAntiSexaquark));

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(19);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill((Double_t)kfAntiSexaquark.GetChi2() /
                                                                                   (Double_t)kfAntiSexaquark.GetNDF());  // exclusive

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(20);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0aSV")]->Fill(dca_v0a_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bSV")]->Fill(dca_v0b_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill((Double_t)kfAntiSexaquark.GetChi2() /
                                                                                      (Double_t)kfAntiSexaquark.GetNDF());  // exclusive

        }  // end of loop over K0S
    }      // end of loop over anti-lambdas
}

/*
 Apply anti-sexaquark candidate cuts.
 - Use: `fPrimaryVertex`
 - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg,
                                                              KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg,
                                                              KFParticle kfKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(2);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(3);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(4);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(5);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(6);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(7);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(8);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(9);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(10);

    Double_t dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(11);

    Double_t dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(12);

    Double_t dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(13);

    Double_t dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(14);

    Double_t dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(15);

    Double_t dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(16);

    Double_t dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(17);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(18);

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
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z,  //
                                                       Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(lvParticle.Angle(posRelativeToRef));
}

/*
 Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
*/
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z,  //
                                                       Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*
 This gives the Armenteros-Podolanski alpha.
 (Based on `AliRoot/STEER/ESD/AliESDv0::AlphaV0()`)
*/
Double_t AliAnalysisTaskSexaquark::ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
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
 This gives the Armenteros-Podolanski qT
 (Based on `AliRoot/STEER/ESD/AliESDv0::PtArmV0()`)
*/
Double_t AliAnalysisTaskSexaquark::ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

    return momNeg.Perp(momTot);
}

/*
 Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 This function stores the point of closest approach, and returns its distance to the PV.
 (PENDING: are we sure of that?)
 - Input: `V0_Px`, `V0_Py`, `V0_Pz`, `V0_X`, `V0_Y`, `V0_Z`, `PV_X`, `PV_Y`, `PV_Z`
 - Return: the distance of closest approach to the PV
*/
Double_t AliAnalysisTaskSexaquark::LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                Double_t refPointX, Double_t refPointY, Double_t refPointZ) {

    TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
    TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
    TVector3 RefVertex(refPointX, refPointY, refPointZ);

    TVector3 CrossProduct = (RefVertex - V0Vertex).Cross(V0Momentum);

    return CrossProduct.Mag() / V0Momentum.Mag();

    /*
    // convert input vectors to arrays
    Double_t ref[3] = {refPointX, refPointY, refPointZ};
    Double_t pos0[3] = {V0_X, V0_Y, V0_Z};
    Double_t dir0[3] = {V0_Px, V0_Py, V0_Pz};

    // lambda function
    auto func = [this, &ref, &pos0, &dir0](const Double_t* t) { return SquaredDistancePointToLine(t, ref, pos0, dir0); };
    ROOT::Math::Functor f(func, 1);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", 0., 0.01);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();

    return TMath::Sqrt(SquaredDistancePointToLine(xs, ref, pos0, dir0));
    */
}

/*
 Calculate transverse impact parameter w.r.t. Primary Vertex.
*/
Float_t AliAnalysisTaskSexaquark::MCRec_GetImpactParameter(AliESDtrack* track) {
    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

    return (Float_t)track->GetD(PV[0], PV[1], fESD->GetMagneticField());
}

/*
 Find a common vertex for two neutral particles, given their position and momentum, assuming they propagate in straight lines.
 If a common vertex cannot be found, then return the DCA between the two lines.
 (Based on Marty Cohen's answer in https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d)
*/
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v1(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    if (v3_dir0.Mag() == 0. || v3_dir1.Mag() == 0.) {
        AliWarning("One of the directions is null!");
        return -1.;
    }

    /* Require perpendicularity, that means both lines must intersect eventually */

    if (v3_dir0.Cross(v3_dir1).Mag() <= 0.) {
        AliWarning("TwoLinesDCA :: Not possible to intersect!");
        return -1.;
    }

    TVector3 diff = v3_pos0 - v3_pos1;
    Double_t determinant = v3_dir0.Dot(v3_dir1) * v3_dir0.Dot(v3_dir1) - v3_dir0.Mag2() * v3_dir1.Mag2();
    if (determinant == 0.) {
        AliWarning("TwoLinesDCA :: Lines are parallel!");
        return -1.;
    }
    Double_t t0 = (v3_dir1.Mag2() * v3_dir0.Dot(diff) - v3_dir1.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;
    Double_t t1 = (-1 * v3_dir0.Mag2() * v3_dir1.Dot(diff) + v3_dir0.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;

    TString aux = "";
    if (t0 >= 0. || t1 >= 0.) {
        aux = "*";
    }

    Double_t dca = (diff + v3_dir0 * t0 - v3_dir1 * t1).Mag();
    PCA0 = v3_pos0 + t0 * v3_dir0;
    PCA1 = v3_pos1 + t1 * v3_dir1;
    TVector3 middle_point((PCA0.X() + PCA1.X()) * 0.5, (PCA0.Y() + PCA1.Y()) * 0.5, (PCA0.Z() + PCA1.Z()) * 0.5);

    AliInfoF("dca pca1 pca2 radius = %f (%f, %f, %f) (%f, %f, %f) %f %s", dca, PCA0.X(), PCA0.Y(), PCA0.Z(), PCA1.X(), PCA1.Y(), PCA1.Z(),
             middle_point.Perp(), aux.Data());

    return dca;
}

/*
 Given the parametric form of a line, find the distance between it and a point.
 - Variable: `t`, the parameter of the line.
 - Parameters: `point`, the point to which we want to calculate the distance.
 - Parameters: `linePoint`, `lineDir`, the position and direction of the line.
 - Return: the square of the distance between the point and the line.
*/
Double_t AliAnalysisTaskSexaquark::SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]) {
    return TMath::Power(point[0] - linePoint[0] - t[0] * lineDir[0], 2) +  //
           TMath::Power(point[1] - linePoint[1] - t[0] * lineDir[1], 2) +  //
           TMath::Power(point[2] - linePoint[2] - t[0] * lineDir[2], 2);
}

/*
 Given the parametric form of two lines, find their distance.
 - Variable: `t`, array of size 2, where `t[0]` and `t[1] are the parameters for the first and second line, respectively.
 - Parameters: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 - Return: the square of the distance between the two lines.
*/
Double_t AliAnalysisTaskSexaquark::SquaredDistanceBetweenLines(const Double_t* t, Double_t pos0[], Double_t dir0[], Double_t pos1[],
                                                               Double_t dir1[]) {
    return TMath::Power(dir1[0] * t[1] + pos1[0] - (dir0[0] * t[0] + pos0[0]), 2) +
           TMath::Power(dir1[1] * t[1] + pos1[1] - (dir0[1] * t[0] + pos0[1]), 2) +
           TMath::Power(dir1[2] * t[1] + pos1[2] - (dir0[2] * t[0] + pos0[2]), 2);
}

/*
 Calculate the distance of closest approach between two lines, given their position and momentum, assuming they propagate in straight lines.
 - Input: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 - Output: `CPA0`, `CPA1`, the closest point of approach for the two lines.
 - Return: the distance of closest approach between the two lines.
 */
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    // convert input vectors to arrays
    Double_t pos0[3] = {v3_pos0.X(), v3_pos0.Y(), v3_pos0.Z()};
    Double_t dir0[3] = {v3_dir0.X(), v3_dir0.Y(), v3_dir0.Z()};

    Double_t pos1[3] = {v3_pos1.X(), v3_pos1.Y(), v3_pos1.Z()};
    Double_t dir1[3] = {v3_dir1.X(), v3_dir1.Y(), v3_dir1.Z()};

    // lambda function
    auto func = [this, &pos0, &dir0, &pos1, &dir1](const Double_t* t) { return SquaredDistanceBetweenLines(t, pos0, dir0, pos1, dir1); };
    ROOT::Math::Functor f(func, 2);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetLimitedVariable(0, "t0", 0.0, 0.1, -999., 0.);
    min->SetLimitedVariable(1, "t1", 0.0, 0.1, -999., 0.);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();
    Double_t t0 = xs[0];
    Double_t t1 = xs[1];

    Double_t dca = TMath::Sqrt(SquaredDistanceBetweenLines(xs, pos0, dir0, pos1, dir1));
    PCA0.SetXYZ(dir0[0] * t0 + pos0[0], dir0[1] * t0 + pos0[1], dir0[2] * t0 + pos0[2]);
    PCA1.SetXYZ(dir1[0] * t1 + pos1[0], dir1[1] * t1 + pos1[1], dir1[2] * t1 + pos1[2]);

    return dca;
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

/*
  Calculate position of a point on a track and some derivatives
  (Copied from `AliRoot/STEER/STEERBase/AliExternalTrackParam.cxx`)
 */
void AliAnalysisTaskSexaquark::Evaluate(const Double_t* h, Double_t t, Double_t r[3], Double_t g[3], Double_t gg[3]) {
    Double_t phase = h[4] * t + h[2];
    Double_t sn = TMath::Sin(phase), cs = TMath::Cos(phase);

    r[0] = h[5];
    r[1] = h[0];
    if (TMath::Abs(h[4]) > kAlmost0) {
        r[0] += (sn - h[6]) / h[4];
        r[1] -= (cs - h[7]) / h[4];
    } else {
        r[0] += t * cs;
        r[1] -= -t * sn;
    }
    r[2] = h[1] + h[3] * t;

    g[0] = cs;
    g[1] = sn;
    g[2] = h[3];

    gg[0] = -h[4] * sn;
    gg[1] = h[4] * cs;
    gg[2] = 0.;
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

    Int_t idx_mc = 0;
    AliMCParticle* mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

    Int_t evt_mother;
    Int_t mother_pid;
    Int_t grandmother_pid;
    // MCGen_GetMotherInfo(mcPart, evt_mother, mother_pid, grandmother_pid);

    Int_t n_daughters;
    Int_t evt_first_dau;
    Int_t evt_last_dau;
    // MCGen_GetDaughtersInfo(mcPart, n_daughters, evt_first_dau, evt_last_dau);

    Float_t x_f = -999.;
    Float_t y_f = -999.;
    Float_t z_f = -999.;
    if (n_daughters) {
        AliMCParticle* mcDaughter = (AliMCParticle*)fMC->GetTrack(mcPart->GetDaughterFirst());
        x_f = mcDaughter->Xv();
        y_f = mcDaughter->Yv();
        z_f = mcDaughter->Zv();
    }

    Bool_t is_signal = 0;  // MCGen_IsSignal(mcPart);

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

    Int_t idx_track = 0;  // fVec_Idx_MCRec[evt_track];
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

    Int_t idx_true = TMath::Abs(track->GetLabel());
    AliMCParticle* mcTrue = (AliMCParticle*)fMC->GetTrack(idx_true);

    Int_t evt_true = 0;          // fMap_Evt_MCGen[idx_true];
    Bool_t is_signal = 0;        // MCGen_IsSignal(mcTrue);
    Bool_t is_duplicate = 0;     // fVec_MCRec_IsDuplicate[evt_track];
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

    V0_tt this_V0;  // = fVec_V0s[evt_v0];

#if DEBUG_MODE == 2
    Int_t evt_pos = this_V0.Idx_Pos;
    Int_t evt_neg = this_V0.Idx_Neg;
    Int_t idx_pos = 0;  // fVec_Idx_MCRec[evt_pos];
    Int_t idx_neg = 0;  // fVec_Idx_MCRec[evt_neg];
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

 */
void AliAnalysisTaskSexaquark::FillTree() {
    AliInfo("MONTE CARLO TREE");
    AliInfo("");
    AliInfo("         EVT    IDX    PID SIGNAL MOTHER  M_PID GM_PID FIRDAU LASDAU");
    /*
    fEvent.N_MCGen = (Int_t)fVec_Idx_MCGen.size();
    for (Int_t ii = 0; ii < fEvent.N_MCGen; ii++) {
        MCGen_PushBack(ii);
    }

    AliInfo("RECONSTRUCTED TRACKS");
    AliInfo("");
    AliInfo("         EVT    IDX EVTTRU IDXTRU DUPLI?");
    fEvent.N_MCRec = (Int_t)fVec_Idx_MCRec.size();
    for (Int_t ii = 0; ii < fEvent.N_MCRec; ii++) {
        MCRec_PushBack(ii);
    }

    AliInfo("FOUND V0s");
    AliInfo("");
    AliInfo("   EVTPOS EVTNEG IDXPOS IDXNEG");
    fEvent.N_V0s = (Int_t)fVec_V0s.size();
    for (Int_t ii = 0; ii < fEvent.N_V0s; ii++) {
        V0_PushBack(ii);
    }
 */
    // (tree operations) after all variables have been assigned, fill tree
    fTree->Fill();
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

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

/*
 Correct initialization of a KFParticle.
 (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
*/
KFParticle AliAnalysisTaskSexaquark::CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge) {

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
KFVertex AliAnalysisTaskSexaquark::CreateKFVertex(const AliVVertex& vertex) {

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
KFParticle AliAnalysisTaskSexaquark::TransportKFParticle(KFParticle kfThis, KFParticle kfOther, Int_t pdgThis, Int_t chargeThis) {

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
 Clear all containers.
*/
void AliAnalysisTaskSexaquark::ClearContainers() {

    /*  */
    getPdgCode_fromMcIdx.clear();

    isMcIdxSecondary.clear();
    isMcIdxSignal.clear();

    getReactionIdx_fromMcIdx.clear();
    getMcIdx_fromReactionIdx.clear();

    doesMcIdxHaveMother.clear();
    getMotherMcIdx_fromMcIdx.clear();
    getNegDauMcIdx_fromMcIdx.clear();
    getPosDauMcIdx_fromMcIdx.clear();

    mcIndicesOfTrueV0s.clear();

    /*  */
    getMcIdx_fromEsdIdx.clear();
    doesEsdIdxPassTrackSelection.clear();
    getEsdIdx_fromMcIdx.clear();

    isEsdIdxSignal.clear();
    getReactionIdx_fromEsdIdx.clear();
    getEsdIdx_fromReactionIdx.clear();

    doesEsdIdxHaveMother.clear();
    getMotherMcIdx_fromEsdIdx.clear();
    getNegDauEsdIdx_fromMcIdx.clear();
    getPosDauEsdIdx_fromMcIdx.clear();

    esdIndicesOfAntiProtonTracks.clear();
    esdIndicesOfPosKaonTracks.clear();
    esdIndicesOfPiPlusTracks.clear();
    esdIndicesOfPiMinusTracks.clear();

    /*  */
    getEsdIdxOfNegDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfPosDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfNegDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfPosDau_fromKaonZeroShortIdx.clear();

    esdAntiLambdas.clear();
    esdKaonsZeroShort.clear();

    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
}
