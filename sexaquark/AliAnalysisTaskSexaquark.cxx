#include "AliAnalysisTaskSexaquark.h"

ClassImp(AliAnalysisTaskSexaquark);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      /*  */
      fIsMC(0),
      fSourceOfV0s(),
      fReactionID(0),
      fDoQA(0),
      fStopAfter(),
      fReadSignalLogs(0),
      fReweightPt(0),
      fReweightRadius(0),
      /*  */
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fEventCuts(),
      fMagneticField(0.),
      fRunNumber(0),
      fCentrality(0.),
      /*  */
      fPDG(),
      fList_Trees(0),
      fTree(0),
      fList_QA_Hists(0),
      fList_Tracks_Hists(0),
      fList_V0s_Hists(0),
      fList_Sexaquarks_Hists(0),
      fList_PosKaonPairs_Hists(0),
      /*  */
      fAliEnPath(),
      fLogTree(0),
      fIsFirstEvent(0),
      fPtWeights(0),
      fRadiusWeights(0),
      /*  */
      fProductsPDG(),
      fFinalStateProductsPDG(),
      fStruckNucleonPDG(0),
      getNegPdgCode_fromV0PdgCode(),
      getPosPdgCode_fromV0PdgCode(),
      /*  */
      fHist_QA(),
      f2DHist_QA(),
      fHist_Tracks_Bookkeep(),
      fHist_Tracks(),
      f2DHist_Tracks_TPCsignal(),
      fHist_V0s_Bookkeep(),
      fHist_V0s(),
      f2DHist_V0s_ArmenterosPodolanski(),
      fHist_AntiSexaquarks_Bookkeep(),
      fHist_AntiSexaquarks(),
      fHist_PosKaonPairs_Bookkeep(),
      fHist_PosKaonPairs(),
      /*  */
      getPdgCode_fromMcIdx(),
      isMcIdxSignal(),
      isMcIdxSecondary(),
      getReactionIdx_fromMcIdx(),
      getMcIdx_fromReactionIdx(),
      getPt_fromReactionIdx(),
      getRadius_fromReactionIdx(),
      doesMcIdxHaveMother(),
      getMotherMcIdx_fromMcIdx(),
      getNegDauMcIdx_fromMcIdx(),
      getPosDauMcIdx_fromMcIdx(),
      mcIndicesOfTrueV0s(),
      /*  */
      getMcIdx_fromEsdIdx(),
      esdIdxPassesTrackSelection_AsProton(),
      esdIdxPassesTrackSelection_AsAntiProton(),
      esdIdxPassesTrackSelection_AsPosKaon(),
      esdIdxPassesTrackSelection_AsNegKaon(),
      esdIdxPassesTrackSelection_AsPiPlus(),
      esdIdxPassesTrackSelection_AsPiMinus(),
      getEsdIdx_fromMcIdx(),
      isEsdIdxSignal(),
      getReactionIdx_fromEsdIdx(),
      getEsdIdx_fromReactionIdx(),
      doesEsdIdxHaveMother(),
      getMotherMcIdx_fromEsdIdx(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx(),
      esdIndicesOfAntiProtonTracks(),
      esdIndicesOfPosKaonTracks(),
      esdIndicesOfPiMinusTracks(),
      esdIndicesOfPiPlusTracks(),
      /*  */
      getEsdIdxOfNegDau_fromAntiLambdaIdx(),
      getEsdIdxOfPosDau_fromAntiLambdaIdx(),
      getEsdIdxOfNegDau_fromKaonZeroShortIdx(),
      getEsdIdxOfPosDau_fromKaonZeroShortIdx(),
      getEsdIdxOfNegDau_fromPionPairIdx(),
      getEsdIdxOfPosDau_fromPionPairIdx(),
      esdAntiLambdas(),
      esdKaonsZeroShort(),
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      kfPionPairs(),
      /*  */
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0),
      kTurnedOn_Track_RejectKinks(0),
      kMin_Track_DCA_wrtPV(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.),
      kMin_Track_Pt(),
      /*  */
      kMin_V0_Mass(),
      kMax_V0_Mass(),
      kMax_V0_Eta(),
      kMax_V0_ArmPtOverAlpha(),
      kMin_V0_Pt(),
      kMin_V0_Radius(),
      kMin_V0_DecayLength(),
      kMax_V0_DecayLength(),
      kMin_V0_CPAwrtPV(),
      kMax_V0_CPAwrtPV(),
      kMin_V0_DCAwrtPV(),
      kMax_V0_DCAbtwDau(),
      kMax_V0_DCAnegV0(),
      kMax_V0_DCAposV0(),
      kMax_V0_ImprvDCAbtwDau(),
      kMax_V0_Chi2(),
      kMax_V0_Chi2ndf(),
      /*  */
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
      kMax_Sexa_Chi2ndf(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_DCAbaSV(0.),
      kMax_Sexa_DCAv0ba(0.),
      kMax_Sexa_DCAv0negba(0.),
      kMax_Sexa_DCAv0posba(0.),
      kMax_Sexa_DCAv0SV(0.),
      kMax_Sexa_DCAv0negSV(0.),
      kMax_Sexa_DCAv0posSV(0.),
      kMax_Sexa_DCApkSV(0.),
      kMax_Sexa_DCApmSV(0.),
      kMax_Sexa_DCAppSV(0.),
      kMax_Sexa_DCAv0pk(0.),
      kMax_Sexa_DCAv0pm(0.),
      kMax_Sexa_DCAv0pp(0.),
      kMax_Sexa_DCApkpm(0.),
      kMax_Sexa_DCApkpp(0.),
      /*  */
      kMin_PosKaonPair_Radius(0.),
      kMin_PosKaonPair_Pt(0.),
      kMin_PosKaonPair_DCAwrtPV(0.),
      kMax_PosKaonPair_DCAbtwKaons(0.),
      kMax_PosKaonPair_DCApkaSV(0.),
      kMax_PosKaonPair_DCApkbSV(0.),
      kMax_PosKaonPair_Chi2ndf(0.)
// fEvent()
{}

/*
 Constructor, called locally.
*/
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name)
    : AliAnalysisTaskSE(name),
      /*  */
      fIsMC(0),
      fSourceOfV0s(),
      fReactionID(0),
      fDoQA(0),
      fStopAfter(),
      fReadSignalLogs(0),
      fReweightPt(0),
      fReweightRadius(0),
      /*  */
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPrimaryVertex(0),
      fPIDResponse(0),
      fEventCuts(),
      fMagneticField(0.),
      fRunNumber(0),
      fCentrality(0.),
      /*  */
      fPDG(),
      fList_Trees(0),
      fTree(0),
      fList_QA_Hists(0),
      fList_Tracks_Hists(0),
      fList_V0s_Hists(0),
      fList_Sexaquarks_Hists(0),
      fList_PosKaonPairs_Hists(0),
      /*  */
      fAliEnPath(),
      fLogTree(0),
      fIsFirstEvent(0),
      fPtWeights(0),
      fRadiusWeights(0),
      /*  */
      fProductsPDG(),
      fFinalStateProductsPDG(),
      fStruckNucleonPDG(0),
      getNegPdgCode_fromV0PdgCode(),
      getPosPdgCode_fromV0PdgCode(),
      /*  */
      fHist_QA(),
      f2DHist_QA(),
      fHist_Tracks_Bookkeep(),
      fHist_Tracks(),
      f2DHist_Tracks_TPCsignal(),
      fHist_V0s_Bookkeep(),
      fHist_V0s(),
      f2DHist_V0s_ArmenterosPodolanski(),
      fHist_AntiSexaquarks_Bookkeep(),
      fHist_AntiSexaquarks(),
      fHist_PosKaonPairs_Bookkeep(),
      fHist_PosKaonPairs(),
      /*  */
      getPdgCode_fromMcIdx(),
      isMcIdxSignal(),
      isMcIdxSecondary(),
      getReactionIdx_fromMcIdx(),
      getMcIdx_fromReactionIdx(),
      getPt_fromReactionIdx(),
      getRadius_fromReactionIdx(),
      doesMcIdxHaveMother(),
      getMotherMcIdx_fromMcIdx(),
      getNegDauMcIdx_fromMcIdx(),
      getPosDauMcIdx_fromMcIdx(),
      mcIndicesOfTrueV0s(),
      /*  */
      getMcIdx_fromEsdIdx(),
      esdIdxPassesTrackSelection_AsProton(),
      esdIdxPassesTrackSelection_AsAntiProton(),
      esdIdxPassesTrackSelection_AsPosKaon(),
      esdIdxPassesTrackSelection_AsNegKaon(),
      esdIdxPassesTrackSelection_AsPiPlus(),
      esdIdxPassesTrackSelection_AsPiMinus(),
      getEsdIdx_fromMcIdx(),
      isEsdIdxSignal(),
      getReactionIdx_fromEsdIdx(),
      getEsdIdx_fromReactionIdx(),
      doesEsdIdxHaveMother(),
      getMotherMcIdx_fromEsdIdx(),
      getNegDauEsdIdx_fromMcIdx(),
      getPosDauEsdIdx_fromMcIdx(),
      esdIndicesOfAntiProtonTracks(),
      esdIndicesOfPosKaonTracks(),
      esdIndicesOfPiMinusTracks(),
      esdIndicesOfPiPlusTracks(),
      /*  */
      getEsdIdxOfNegDau_fromAntiLambdaIdx(),
      getEsdIdxOfPosDau_fromAntiLambdaIdx(),
      getEsdIdxOfNegDau_fromKaonZeroShortIdx(),
      getEsdIdxOfPosDau_fromKaonZeroShortIdx(),
      getEsdIdxOfNegDau_fromPionPairIdx(),
      getEsdIdxOfPosDau_fromPionPairIdx(),
      esdAntiLambdas(),
      esdKaonsZeroShort(),
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      kfPionPairs(),
      /*  */
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0),
      kTurnedOn_Track_RejectKinks(0),
      kMin_Track_DCA_wrtPV(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.),
      kMin_Track_Pt(),
      /*  */
      kMin_V0_Mass(),
      kMax_V0_Mass(),
      kMax_V0_Eta(),
      kMax_V0_ArmPtOverAlpha(),
      kMin_V0_Pt(),
      kMin_V0_Radius(),
      kMin_V0_DecayLength(),
      kMax_V0_DecayLength(),
      kMin_V0_CPAwrtPV(),
      kMax_V0_CPAwrtPV(),
      kMin_V0_DCAwrtPV(),
      kMax_V0_DCAbtwDau(),
      kMax_V0_DCAnegV0(),
      kMax_V0_DCAposV0(),
      kMax_V0_ImprvDCAbtwDau(),
      kMax_V0_Chi2(),
      kMax_V0_Chi2ndf(),
      /*  */
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
      kMax_Sexa_Chi2ndf(0.),
      kMax_Sexa_DCAbtwV0s(0.),
      kMax_Sexa_DCAv0aSV(0.),
      kMax_Sexa_DCAv0bSV(0.),
      kMax_Sexa_DCAv0anegSV(0.),
      kMax_Sexa_DCAv0aposSV(0.),
      kMax_Sexa_DCAv0bnegSV(0.),
      kMax_Sexa_DCAv0bposSV(0.),
      kMax_Sexa_DCAbaSV(0.),
      kMax_Sexa_DCAv0ba(0.),
      kMax_Sexa_DCAv0negba(0.),
      kMax_Sexa_DCAv0posba(0.),
      kMax_Sexa_DCAv0SV(0.),
      kMax_Sexa_DCAv0negSV(0.),
      kMax_Sexa_DCAv0posSV(0.),
      kMax_Sexa_DCApkSV(0.),
      kMax_Sexa_DCApmSV(0.),
      kMax_Sexa_DCAppSV(0.),
      kMax_Sexa_DCAv0pk(0.),
      kMax_Sexa_DCAv0pm(0.),
      kMax_Sexa_DCAv0pp(0.),
      kMax_Sexa_DCApkpm(0.),
      kMax_Sexa_DCApkpp(0.),
      /*  */
      kMin_PosKaonPair_Radius(0.),
      kMin_PosKaonPair_Pt(0.),
      kMin_PosKaonPair_DCAwrtPV(0.),
      kMax_PosKaonPair_DCAbtwKaons(0.),
      kMax_PosKaonPair_DCApkaSV(0.),
      kMax_PosKaonPair_DCApkbSV(0.),
      kMax_PosKaonPair_Chi2ndf(0.)
// fEvent()
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());  // fList_Trees
    DefineOutput(2, TList::Class());  // fList_QA_Hists
    DefineOutput(3, TList::Class());  // fList_Tracks_Hists
    DefineOutput(4, TList::Class());  // fList_V0s_Hists
    DefineOutput(5, TList::Class());  // fList_Sexaquarks_Hists
    DefineOutput(6, TList::Class());  // fList_PosKaonPairs_Hists
}

/*
 Destructor.
 - Note: if `TList::SetOwner(kTRUE)` was called, the TList destructor should delete all objects added to it
*/
AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() {
    if (fList_Trees) delete fList_Trees;
    if (fList_QA_Hists) delete fList_QA_Hists;
    if (fList_Tracks_Hists) delete fList_Tracks_Hists;
    if (fList_V0s_Hists) delete fList_V0s_Hists;
    if (fList_Sexaquarks_Hists) delete fList_Sexaquarks_Hists;
    if (fList_PosKaonPairs_Hists) delete fList_PosKaonPairs_Hists;
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {

    AliInfo("!! It begins !!");

    /** Add mandatory routines **/

    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) AliFatal("ERROR: AliAnalysisManager couldn't be found.");

    AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
    if (!inputHandler) AliFatal("ERROR: AliESDInputHandler couldn't be found.");

    fPIDResponse = inputHandler->GetPIDResponse();

    /** Prepare output  **/

    /* Trees */

    fList_Trees = new TList();
    fList_Trees->SetOwner(kTRUE);

    fLogTree = nullptr;
    if (fReadSignalLogs) {
        fLogTree = new TTree("Injected", "Injected");
        fList_Trees->Add(fLogTree);
    }

    fTree = new TTree("Events", "Events");
    // SetBranches();
    fList_Trees->Add(fTree);
    PostData(1, fList_Trees);

    /* Histograms */

    fList_QA_Hists = new TList();
    fList_QA_Hists->SetOwner(kTRUE);
    if (fDoQA) PrepareQAHistograms();
    // ! NOTE:                                                                                                    ! //
    // ! the following lines are commented on purpose, adding these histograms makes sense for debugging purposes ! //
    // ! but not for production purposes, as they would only increase the filesize when merging the output files  ! //
    /*
    if (fReweightPt) {
        for (Int_t cc = 0; cc < (Int_t)fPtWeights.size(); cc++) {
            fList_QA_Hists->Add(fPtWeights[cc]);
        }
    }
    if (fReweightRadius) fList_QA_Hists->Add(fRadiusWeights);
    */
    PostData(2, fList_QA_Hists);

    fList_Tracks_Hists = new TList();
    fList_Tracks_Hists->SetOwner(kTRUE);
    PrepareTracksHistograms();
    PostData(3, fList_Tracks_Hists);

    fList_V0s_Hists = new TList();
    fList_V0s_Hists->SetOwner(kTRUE);
    if (fReactionID != 'H') PrepareV0Histograms();
    PostData(4, fList_V0s_Hists);

    fList_Sexaquarks_Hists = new TList();
    fList_Sexaquarks_Hists->SetOwner(kTRUE);
    if (fReactionID != 'H') PrepareSexaquarkHistograms();
    PostData(5, fList_Sexaquarks_Hists);

    fList_PosKaonPairs_Hists = new TList();
    fList_PosKaonPairs_Hists->SetOwner(kTRUE);
    if (fReactionID == 'H') PreparePosKaonPairHistograms();
    PostData(6, fList_PosKaonPairs_Hists);

    AliInfo("!! It ends !!");
}

/*
 Define and add QA histograms to the histograms output list.
 - Uses: `fHist_QA`, `f2DHist_QA`
*/
void AliAnalysisTaskSexaquark::PrepareQAHistograms() {

    /* 1D Hists */

    // clang-format off
    std::vector<TString> QA_props = {// MC. Gen. Event
                                     "MCGen_VertexX", "MCGen_VertexY", "MCGen_VertexZ",
                                     // Event
                                     "Events_Bookkeep", "VertexX", "VertexY", "VertexZ",
                                     "NTracks", "NSelectedTracks",
                                     // Tracks
                                     "Tracks_DCAxy", "Tracks_DCAz", "Tracks_Pt",
                                     // SPD
                                     "SPD_NTracklets",
                                     // TPC
                                     "TPC_NClusters"};
    std::vector<Int_t> QA_nbins = {100, 100, 200,
                                   5, 100, 100, 200,
                                   2000, 2000,
                                   100, 100, 200,
                                   60,
                                   200};
    std::vector<Float_t> QA_min = {-1, -1, -50,
                                   0., -1, -1, -50,
                                   0., 0.,
                                   -5., -5., 0.,
                                   0.,
                                   0.};
    std::vector<Float_t> QA_max = {1, 1, 50,
                                   5., 1, 1, 50,
                                   20000., 20000.,
                                   5., 5., 100.,
                                   6000,
                                   5000000.};
    // clang-format on

    TString histKey;
    for (Int_t prop_idx = 0; prop_idx < (Int_t)QA_props.size(); prop_idx++) {
        histKey = QA_props[prop_idx];
        fHist_QA[histKey] = new TH1F("QA_" + histKey, "", QA_nbins[prop_idx], QA_min[prop_idx], QA_max[prop_idx]);
        fList_QA_Hists->Add(fHist_QA[histKey]);
    }

    /* 2D Hists */

    // clang-format off
    QA_props = {// Tracks
                "Tracks_PhiVsEta",
                // SPD
                "SPD_MultiplicityVsVertexZ", "SPD_NTrackletsClu0VsNTracklets", "SPD_NTrackletsClu1VsNTracklets", "SPD_PhiVsEta",
                // ITS
                "ITS_LayerNoVsPhi", "ITS_LayerNoVsEta",
                // TPC
                "TPC_NSigmaProtonVsInnerParamP", "TPC_NSigmaKaonVsInnerParamP", "TPC_NSigmaPionVsInnerParamP",
                "TPC_NSigmaProtonVsEta", "TPC_NSigmaKaonVsEta", "TPC_NSigmaPionVsEta"};
    std::vector<Int_t> QA_nbinsx = {100,
                                   10, 200, 200, 400,
                                   200, 40,
                                   200, 200, 200, 200, 200, 200};
    std::vector<Float_t> QA_xlow = {-1.,
                                   -10., 0., 0., -1.,
                                   0., -1.,
                                   0., 0., 0., -1., -1., -1.};
    std::vector<Float_t> QA_xup = {1.,
                                   10., 3000., 3000., 1.,
                                   2 * TMath::Pi(), 1.,
                                   10., 10., 10., 1., 1., 1.};
    std::vector<Int_t> QA_nbinsy = {200,
                                   30, 200, 200, 400,
                                   6, 6,
                                   200, 200, 200, 200, 200, 200};
    std::vector<Float_t> QA_ylow = {0.,
                                    0., 0., 0., 0.,
                                    0., 0.,
                                    -5., -5., -5., -5., -5., -5.};
    std::vector<Float_t> QA_yup = {2 * TMath::Pi(),
                                   6000., 6000., 6000., 2 * TMath::Pi(),
                                   6., 6.,
                                   5., 5., 5., 5., 5., 5.};
    // clang-format on

    for (Int_t prop_idx = 0; prop_idx < (Int_t)QA_props.size(); prop_idx++) {
        histKey = QA_props[prop_idx];
        if (fReactionID == 'A' && histKey.Contains("Kaon")) continue;
        f2DHist_QA[histKey] = new TH2F("QA_" + histKey, "", QA_nbinsx[prop_idx], QA_xlow[prop_idx], QA_xup[prop_idx], QA_nbinsy[prop_idx],
                                       QA_ylow[prop_idx], QA_yup[prop_idx]);
        fList_QA_Hists->Add(f2DHist_QA[histKey]);
    }
}

/*
 Define and add tracks' histograms to the histograms output list.
 - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareTracksHistograms() {

    TString tracks_stages[2] = {"MCGen", "Found"};
    std::vector<TString> tracks_sets = {"All", "True", "Secondary", "Signal"};
    std::set<Int_t> tracks_species(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());
    std::map<Int_t, TString> tracks_name = {{-2212, "AntiProton"}, {321, "PosKaon"}, {211, "PiPlus"}, {-211, "PiMinus"}};

    // clang-format off
    std::vector<TString> tracks_props = {"Px",           "Py",          "Pt",         "Pz",           "Eta",
                                         "DCA_wrtPV",    "DCAxy_wrtPV", "DCAz_wrtPV", "NTPCClusters", "NCrossedRows",
                                         "NFindableCls", "NSharedCls",  "Chi2/NCls",
                                         "NSigmaProton", "NSigmaKaon",  "NSigmaPion", "Status"};
    std::vector<Int_t> tracks_nbins = {100, 100, 100, 100, 100,
                                       100, 100, 100, 160, 160,
                                       100, 100, 140,
                                       100, 100, 100, 16};
    std::vector<Double_t> tracks_min = {-20., -20., 0.,  -20., -1.2,
                                        0.,   0.,   0.,  0.,   0.,
                                        0.,   0.,   0.,
                                        -10., -10,  -10, 0.};
    std::vector<Double_t> tracks_max = {20.,  20.,  10.,  20.,  1.2,
                                        200., 200., 200., 160., 160.,
                                        100., 100., 7.,
                                        10,   10,   10,   16};
    // clang-format on

    for (Int_t species : tracks_species) {
        fHist_Tracks_Bookkeep[species] = new TH1F(Form("%s_Bookkeep", tracks_name[species].Data()), "", 150, 0., 150.);
        fList_Tracks_Hists->Add(fHist_Tracks_Bookkeep[species]);
    }

    for (TString& stage : tracks_stages) {
        for (TString& set : tracks_sets) {
            for (Int_t species : tracks_species) {
                for (Int_t prop_idx = 0; prop_idx < (Int_t)tracks_props.size(); prop_idx++) {

                    if (!fIsMC && (stage == "MCGen" || set == "True" || set == "Secondary" || set == "Signal")) continue;

                    if (stage == "MCGen" && set == "True") continue;
                    if (stage == "MCGen" && !(tracks_props[prop_idx] == "Px" || tracks_props[prop_idx] == "Py" || tracks_props[prop_idx] == "Pt" ||
                                              tracks_props[prop_idx] == "Pz" || tracks_props[prop_idx] == "Eta")) {
                        continue;
                    }

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), tracks_name[species].Data(), tracks_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, species, tracks_props[prop_idx]);
                    fHist_Tracks[histKey] = new TH1F(histName, "", tracks_nbins[prop_idx], tracks_min[prop_idx], tracks_max[prop_idx]);
                    fList_Tracks_Hists->Add(fHist_Tracks[histKey]);
                }
            }
        }
    }

    tracks_sets.insert(tracks_sets.begin() + 1, "Selected");
    for (TString& set : tracks_sets) {
        TString histName = Form("Tracks_%s_TPCsignal", set.Data());
        f2DHist_Tracks_TPCsignal[set] = new TH2F(histName, "", 200, 0., 10., 200, 0., 400.);
        fList_Tracks_Hists->Add(f2DHist_Tracks_TPCsignal[set]);
    }
}

/*
 Define and add V0s' histograms to the histograms output list.
 - Uses: `fSourceOfV0s`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareV0Histograms() {

    TString V0_stages[3] = {"MCGen", "Findable", "Found"};
    TString V0_sets[4] = {"All", "True", "Secondary", "Signal"};

    std::vector<Int_t> V0_species = {-3122, 310, 422};
    std::vector<TString> V0_name = {"AntiLambda", "KaonZeroShort", "PionPair"};

    // clang-format off
    std::vector<TString> V0_props = {"Mass",        "Pt",           "Px",          "Py",       "Pz",        "Phi",
                                     "OriginRadius", "DecayRadius", "Zv",       "Eta",       "Rapidity",
                                     "DecayLength", "OpeningAngle", "ArmQt",    "ArmAlpha",       "CPAwrtPV",    "DCAwrtPV",
                                     "DCAbtwDau",   "DCAnegV0",     "DCAposV0", "ImprvDCAbtwDau", "Chi2",        "Chi2ndf"};
    std::vector<Int_t> V0_nbins = {100, 100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100, 100,
                                   100, 100, 100, 100, 100, 100};
    std::vector<Double_t> V0_min = {0., 0., -20., -20., -20., 0.,
                                    0., 0.,   -50., -4.,  -4.,
                                    0., 0., 0.,   -1.,  0.,   0.,
                                    0., 0., 0.,   0.,   0.,   0.};
    std::vector<Double_t> V0_max = {10.,  10.,         20.,  20., 20., 2. * TMath::Pi(),
                                    250.,        250., 50., 4.,  4.,
                                    500., TMath::Pi(), 0.3,  1.,  1.,  200.,
                                    2.,   2.,          2.,   1.,  2.,  0.5};
    // clang-format on

    for (Int_t sp = 0; sp < (Int_t)V0_species.size(); sp++) {
        if (fReactionID != 'A' && V0_species[sp] == 310) continue;
        if (fReactionID != 'E' && V0_species[sp] == 422) continue;
        fHist_V0s_Bookkeep[V0_species[sp]] = new TH1F(Form("%s_Bookkeep", V0_name[sp].Data()), "", 200, 0., 200.);
        fList_V0s_Hists->Add(fHist_V0s_Bookkeep[V0_species[sp]]);
    }

    for (TString& stage : V0_stages) {
        for (TString& set : V0_sets) {
            for (Int_t sp = 0; sp < (Int_t)V0_species.size(); sp++) {
                for (Int_t prop_idx = 0; prop_idx < (Int_t)V0_props.size(); prop_idx++) {

                    if (!fIsMC && (stage == "MCGen" || stage == "Findable" || set == "True" || set == "Secondary" || set == "Signal")) continue;
                    if (fReactionID != 'A' && V0_species[sp] == 310) continue;
                    if (fReactionID != 'E' && V0_species[sp] == 422) continue;

                    if (fSourceOfV0s == "kalman" && (V0_props[prop_idx] == "ImprvDCAbtwDau" || V0_props[prop_idx] == "Chi2")) continue;
                    if ((fSourceOfV0s == "offline" || fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "custom") && V0_props[prop_idx] == "Chi2ndf") {
                        continue;
                    }

                    if (stage == "Findable" && set == "All") continue;
                    if (stage == "MCGen" && set == "True") continue;

                    if (V0_species[sp] == 422) {
                        if (stage == "MCGen") continue;
                        if (stage == "Findable" || stage == "Found") {
                            if (set == "True" || set == "Secondary") continue;
                        }
                    }

                    if ((stage == "Findable" || stage == "MCGen") &&
                        (V0_props[prop_idx] == "DCAbtwDau" || V0_props[prop_idx] == "DCAnegV0" || V0_props[prop_idx] == "DCAposV0" ||
                         V0_props[prop_idx] == "ImprvDCAbtwDau" || V0_props[prop_idx] == "Chi2" || V0_props[prop_idx] == "Chi2ndf")) {
                        continue;
                    }

                    if (stage == "Found" && V0_props[prop_idx] == "OriginRadius") continue;

                    TString histName = Form("%s_%s_%s_%s", stage.Data(), set.Data(), V0_name[sp].Data(), V0_props[prop_idx].Data());
                    std::tuple<TString, TString, Int_t, TString> histKey = std::make_tuple(stage, set, V0_species[sp], V0_props[prop_idx]);
                    Float_t minEdge = V0_props[prop_idx] == "Mass" ? kMin_V0_Mass[V0_species[sp]] : V0_min[prop_idx];
                    Float_t maxEdge = V0_props[prop_idx] == "Mass" ? kMax_V0_Mass[V0_species[sp]] : V0_max[prop_idx];
                    fHist_V0s[histKey] = new TH1F(histName, "", V0_nbins[prop_idx], minEdge, maxEdge);
                    fList_V0s_Hists->Add(fHist_V0s[histKey]);
                }
            }
        }
    }

    for (TString& set : V0_sets) {
        TString histName = Form("Found_%s_V0s_ArmenterosPodolanski", set.Data());
        f2DHist_V0s_ArmenterosPodolanski[set] = new TH2F(histName, "", 200, -1., 1., 200, 0., 0.3);
        fList_V0s_Hists->Add(f2DHist_V0s_ArmenterosPodolanski[set]);
    }
}

/*
 Define and add anti-sexaquark histograms to the histograms output list.
 - Uses: `fHist_AntiSexaquarks_Bookkeep`, `fList_Sexaquarks_Hists`, `fHist_AntiSexaquarks`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PrepareSexaquarkHistograms() {

    TString SQ_stages[7] = {"MCGen_BI", "MCGen_AI", "Findable", "Found", "wRadiusWeights", "wPtWeights", "wBothWeights"};
    TString SQ_sets[2] = {"All", "Signal"};

    std::vector<TString> SQ_props;
    std::vector<Int_t> SQ_nbins;
    std::vector<Double_t> SQ_min;
    std::vector<Double_t> SQ_max;

    // clang-format off
    if (fReactionID == 'A') {
        SQ_props = {"Mass",     "Pt",       "Px",          "Py",          "Pz",          "Phi",         "Radius",
                    "Zv",       "Eta",      "Rapidity",    "DecayLength", "CPAwrtPV",    "DCAwrtPV",    "DCAbtwV0s",
                    "DCAv0aSV", "DCAv0bSV", "DCAv0anegSV", "DCAv0aposSV", "DCAv0bnegSV", "DCAv0bposSV", "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., -2 * TMath::Pi(), 0.,
                  -50., -2.,  -5.,  0.,   0.,   0.,               0.,
                  0.,   0.,   0.,   0.,   0.,   0.,               0.};
        SQ_max = {10., 10., 10., 10.,  50., 2 * TMath::Pi(), 200.,
                  50., 2.,  5.,  500., 1.,  10.,             2.,
                  2.,  2.,  10., 10.,  10., 10.,             1.};
    } else if (fReactionID == 'D') {
        SQ_props = {"Mass",       "Pt",         "Px",       "Py",          "Pz",         "Phi",        "Radius",
                    "Zv",         "Eta",        "Rapidity", "DecayLength", "CPAwrtPV",   "DCAwrtPV",   "DCAv0ba",
                    "DCAv0negba", "DCAv0posba", "DCAv0SV",  "DCAbaSV",     "DCAv0negSV", "DCAv0posSV", "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0.,
                  -50., -2.,  -5.,  0.,   0.,   0., 0.,
                  0.,   0.,   0.,   0.,   0.,   0., 0.};
        SQ_max = {10., 10., 10., 10.,  50., 2 * TMath::Pi(), 200.,
                  50., 2.,  5.,  500., 1.,  10.,             10.,
                  10., 10., 10., 10.,  10., 10.,             1.};
    } else if (fReactionID == 'E') {
        SQ_props = {"Mass",    "Pt",       "Px",          "Py",       "Pz",       "Phi",     "Radius",     "Zv",
                    "Eta",     "Rapidity", "DecayLength", "CPAwrtPV", "DCAwrtPV", "DCAv0SV", "DCAv0negSV", "DCAv0posSV",
                    "DCApkSV", "DCApmSV",  "DCAppSV",     "DCAv0pk",  "DCAv0pm",  "DCAv0pp", "DCApkpm",    "DCApkpp",
                    "Chi2ndf"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100, 100, 100, 100,
                    100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0., -50.,
                  -2.,  -5.,  0.,   0.,   0.,   0., 0., 0.,
                  0.,   0.,   0.,   0.,   0.,   0., 0., 0.,
                  0.};
        SQ_max = {10., 10., 10.,  10., 50., 2 * TMath::Pi(), 200., 50.,
                  2.,  5.,  500., 1.,  10., 10.,             10.,  10.,
                  10., 10., 10.,  10., 10., 10.,             10.,  1.,
                  30.};
    } else if (fReactionID == 'H') {
        SQ_props = {"Mass", "Pt",       "Px",          "Py",       "Pz",       "Phi",   "Radius",   "Zv",
                    "Eta",  "Rapidity", "DecayLength", "CPAwrtPV", "DCAwrtPV"};
        SQ_nbins = {100, 100, 100, 100, 100, 100, 100, 100,
                    100, 100, 100, 100, 100};
        SQ_min = {-10., -10., -10., -10., -50., 0., 0., -50.,
                  -2.,  -5.,  0.,   0.,   0.};
        SQ_max = {10., 10., 10.,  10., 50., 2 * TMath::Pi(), 200., 50.,
                  2.,  5.,  500., 1.,  10.};
    }
    // clang-format on

    fHist_AntiSexaquarks_Bookkeep = new TH1D("AntiSexaquarks_Bookkeep", "", 100, 0., 100.);
    fList_Sexaquarks_Hists->Add(fHist_AntiSexaquarks_Bookkeep);

    for (TString& stage : SQ_stages) {
        for (TString& set : SQ_sets) {
            for (Int_t prop_idx = 0; prop_idx < (Int_t)SQ_props.size(); prop_idx++) {
                if (!fIsMC && (stage == "MCGen_BI" || stage == "MCGen_AI" || stage == "Findable" || set == "Signal")) continue;
                if (fReactionID == 'H' && (stage == "Findable" || stage == "Found")) continue;

                if (set == "Signal" && (stage == "MCGen_BI" || stage == "MCGen_AI" || stage == "Findable")) continue;

                if (fSourceOfV0s != "kalman" && SQ_props[prop_idx] == "Chi2ndf") continue;

                if (stage == "MCGen_BI" && (SQ_props[prop_idx] == "Zv" || SQ_props[prop_idx] == "Radius" || SQ_props[prop_idx] == "DecayLength"))
                    continue;
                if ((stage == "MCGen_BI" || stage == "MCGen_AI" || stage == "Findable") &&
                    ((SQ_props[prop_idx].Contains("DCA") && SQ_props[prop_idx] != "DCAwrtPV") || SQ_props[prop_idx].Contains("Chi2"))) {
                    continue;
                }

                if ((stage == "wRadiusWeights" || stage == "wPtWeights" || stage == "wBothWeights") && set != "Signal") continue;
                if (stage == "wRadiusWeights" && !fReweightRadius) continue;
                if (stage == "wPtWeights" && !fReweightPt) continue;
                if (stage == "wBothWeights" && (!fReweightPt || !fReweightRadius)) continue;

                TString histName = Form("%s_%s_AntiSexaquark_%s", stage.Data(), set.Data(), SQ_props[prop_idx].Data());
                std::tuple<TString, TString, TString> histKey = std::make_tuple(stage, set, SQ_props[prop_idx]);
                fHist_AntiSexaquarks[histKey] = new TH1D(histName, "", SQ_nbins[prop_idx], SQ_min[prop_idx], SQ_max[prop_idx]);

                if (stage == "wRadiusWeights" || stage == "wPtWeights" || stage == "wBothWeights") fHist_AntiSexaquarks[histKey]->Sumw2();

                fList_Sexaquarks_Hists->Add(fHist_AntiSexaquarks[histKey]);
            }
        }
    }
}

/*
 Define and add histograms for the K+K+ pairs to the histograms output list.
 - Uses: `fSourceOfV0s`, `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::PreparePosKaonPairHistograms() {

    TString stages[3] = {"MCGen", "Findable", "Found"};
    TString sets[2] = {"All", "Signal"};

    // clang-format off
    std::vector<TString> props = {// true + rec. vars.
                                  "Mass", "OriginRadius", "DecayLength", "Xv", "Yv", "Zv", "Eta", "Rapidity", "Px", "Py", "Pz", "Pt", "OpeningAngle",
                                  // rec. vars.
                                  "DCAwrtPV", "DCAbtwKaons", "DCApkaSV", "DCApkbSV", "Chi2ndf"};
    std::vector<Int_t> nbins = {// true + rec. vars.
                                100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                                 // rec. vars.
                                100, 100, 100, 100, 150};
    std::vector<Double_t> min = {// true + rec. vars
                                 0., 0., 0., -50., -300., -300., -4., -1., -10., -10., -10., 0., 0.,
                                 // rec. vars.
                                 0., 0., 0., 0., 0.};
    std::vector<Double_t> max = {// true + rec. vars
                                 10., 250., 500., 50., 300., 300., 4., 1., 10., 10., 10., 15., TMath::Pi(),
                                 // rec. vars.
                                 200., 20., 20., 20., 1.5};
    // clang-format on

    fHist_PosKaonPairs_Bookkeep = new TH1F("PosKaonPairs_Bookkeep", "", 25, 0., 25.);
    fList_PosKaonPairs_Hists->Add(fHist_PosKaonPairs_Bookkeep);

    for (TString& stage : stages) {
        for (TString& set : sets) {
            for (Int_t prop_idx = 0; prop_idx < (Int_t)props.size(); prop_idx++) {

                if (!fIsMC && (stage == "MCGen" || stage == "Findable" || set == "Signal")) continue;

                if (stage == "Findable" && set == "Signal") continue;

                if ((stage == "Findable" || stage == "MCGen") && prop_idx > 12) continue;

                TString histName = Form("%s_%s_PosKaonPair_%s", stage.Data(), set.Data(), props[prop_idx].Data());
                std::tuple<TString, TString, TString> histKey = std::make_tuple(stage, set, props[prop_idx]);
                fHist_PosKaonPairs[histKey] = new TH1F(histName, "", nbins[prop_idx], min[prop_idx], max[prop_idx]);
                fList_PosKaonPairs_Hists->Add(fHist_PosKaonPairs[histKey]);
            }
        }
    }
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fPrimaryVertex`, `fMagneticField`, `fSourceOfV0s`, `fList_Trees`,
 `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

    /* Handle external logs */

    if (fIsFirstEvent) {
        if (fAliEnPath == "") AliInfo("!! No luck finding fAliEnPath !!");
        AliInfoF("!! fAliEnPath: %s !!", fAliEnPath.Data());
        LoadLogsIntoTree();
        // if (LoadLogsIntoTree()) fLogTree->Print();
        fIsFirstEvent = kFALSE;
    }

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
    fRunNumber = fESD->GetRunNumber();

    // print event centrality class
    AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*>(fESD->FindListObject("MultSelection"));
    if (!MultSelection) AliFatal("ERROR: AliMultSelection couldn't be found.");
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    AliInfoF("!! fCentrality = %f !!", fCentrality);

    // init KFparticle
    KFParticle::SetField(fMagneticField);

    if (!PassesEventSelection()) return;

    /* MC */

    if (fIsMC) {
        ProcessMCGen();
        ProcessSignalInteractions();
    }

    if (fStopAfter == "MC") {
        EndOfEvent();
        return;
    }

    /* Tracks */

    ProcessTracks();

    if (fStopAfter == "Tracks") {
        EndOfEvent();
        return;
    }

    /* Findables */

    if (fIsMC) {
        if (fReactionID == 'H') {
            ProcessFindablePosKaonPairs();
        } else {
            ProcessFindableV0s();
            if (fReactionID == 'E') ProcessFindablePionPairs();
            ProcessFindableSexaquarks();
        }
    }

    if (fStopAfter == "Findables") {
        EndOfEvent();
        return;
    }

    /* V0s */

    if (fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "offline") {
        GetV0sFromESD(fSourceOfV0s == "on-the-fly");
        if (fReactionID == 'A') GeoSexaquarkFinder_ChannelA();
        // if (fReactionID == 'D') GeoSexaquarkFinder_ChannelD();
        // if (fReactionID == 'E') GeoSexaquarkFinder_ChannelE();
    }

    if (fSourceOfV0s == "custom") {
        if (fReactionID == 'A') {  // "AntiSexaquark,N->AntiLambda,K0S"
            CustomV0Finder(-3122);
            CustomV0Finder(310);
        }
    }

    if (fSourceOfV0s == "kalman") {
        if (fReactionID == 'A') {  // "AntiSexaquark,N->AntiLambda,K0S"
            KalmanV0Finder(-3122, -2212, 211);
            KalmanV0Finder(310, -211, 211);
        } else if (fReactionID == 'D') {  // "AntiSexaquark,P->AntiLambda,K+"
            KalmanV0Finder(-3122, -2212, 211);
        } else if (fReactionID == 'E') {  // "AntiSexaquark,P->AntiLambda,K+,pi-,pi+"
            KalmanV0Finder(-3122, -2212, 211);
            KalmanV0Finder(422, -211, 211);
        } else if (fReactionID == 'H') {  // "AntiSexaquark,P->AntiProton,K+,K+,pi0"
            KalmanPosKaonPairFinder();
        }
    }

    if (fStopAfter == "V0s") {
        EndOfEvent();
        return;
    }

    /* Sexaquarks */

    // clang-format off
    if (fSourceOfV0s != "kalman") {
        if (fReactionID == 'A') GeoSexaquarkFinder_ChannelA();
        else if (fReactionID == 'D') GeoSexaquarkFinder_ChannelD();
        else if (fReactionID == 'E') GeoSexaquarkFinder_ChannelE();
    } else {
        if (fReactionID == 'A') KalmanSexaquarkFinder_ChannelA();
        else if (fReactionID == 'D') KalmanSexaquarkFinder_ChannelD();
        else if (fReactionID == 'E') KalmanSexaquarkFinder_ChannelE();
    }
    // clang-format on

    EndOfEvent();
}

/*
 End of event.
*/
void AliAnalysisTaskSexaquark::EndOfEvent() {

    // FillTree();

    ClearContainers();

    PostData(1, fList_Trees);
    PostData(2, fList_QA_Hists);
    PostData(3, fList_Tracks_Hists);
    PostData(4, fList_V0s_Hists);
    PostData(5, fList_Sexaquarks_Hists);
    PostData(6, fList_PosKaonPairs_Hists);
}

/*
 User implementation of Notify(). Needed for reading the AliEn path.
 This function is loaded during AliAnalysisManager::Notify().
 It's called after UserCreateOutputObjects(), but before each UserExec().
*/
Bool_t AliAnalysisTaskSexaquark::UserNotify() {
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) AliFatal("!! Analysis Manager not found !!");
    TTree* man_tree = man->GetTree();
    if (!man_tree) AliFatal("!! Analysis Manager Tree not found !!");
    TFile* man_file = man_tree->GetCurrentFile();
    if (!man_file) AliFatal("!! Analysis Manager File not found !!");
    fAliEnPath = man_file->GetName();
    AliInfoF("!! Loaded AliEn Path: %s !!", fAliEnPath.Data());

    fIsFirstEvent = kTRUE;

    return kTRUE;
}

/*
 Initialize analysis task.
*/
void AliAnalysisTaskSexaquark::Initialize() {

    std::array<TString, 4> validSourcesOfV0s = {"offline", "on-the-fly", "custom", "kalman"};

    if (std::find(validSourcesOfV0s.begin(), validSourcesOfV0s.end(), fSourceOfV0s) == validSourcesOfV0s.end()) {
        AliWarning("Source of V0s must be \"offline\", \"on-the-fly\", \"custom\" or \"kalman\".");
        AliWarning("=> Setting it to \"kalman\".");
        fSourceOfV0s = "kalman";
    }

    if (fReactionID == 'H' && (fSourceOfV0s == "on-the-fly" || fSourceOfV0s == "offline")) {
        AliWarning("Source of V0s must be \"custom\" or \"kalman\" for the reaction channel H.");
        AliWarning("=> Setting it to \"kalman\".");
        fSourceOfV0s = "kalman";
    }

    std::array<Char_t, 4> validReactionIDs = {'A', 'D', 'E', 'H'};
    if (std::find(validReactionIDs.begin(), validReactionIDs.end(), fReactionID) == validReactionIDs.end()) {
        AliWarning("Reaction ID must be set.");
        AliWarning("=> Setting it to 'A'.");
        fReactionID = 'A';
    }

    if (fReweightPt && !fReadSignalLogs) AliFatal("Signal logs must be read for reweighting the pT distribution.");

    /* Define cuts */

    DefineTracksCuts("standard");
    // clang-format off
    if (fReactionID != 'H') DefineV0Cuts("standard");
    else DefinePosKaonPairCuts("standard");
    // clang-format on
    DefineSexaquarkCuts("standard");

    getNegPdgCode_fromV0PdgCode[-3122] = -2212;
    getPosPdgCode_fromV0PdgCode[-3122] = 211;
    getNegPdgCode_fromV0PdgCode[310] = -211;
    getPosPdgCode_fromV0PdgCode[310] = 211;

    /* Assign the rest of properties related to the respective Reaction Channel */

    TString ReactionChannel;
    if (fReactionID == 'A') {
        ReactionChannel = "AntiSexaquark,N->AntiLambda,K0S";
        fProductsPDG = {-3122, 310};
        fFinalStateProductsPDG = {-2212, 211, -211, 211};
        fStruckNucleonPDG = 2112;
    } else if (fReactionID == 'D') {
        ReactionChannel = "AntiSexaquark,P->AntiLambda,K+";
        fProductsPDG = {-3122, 321};
        fFinalStateProductsPDG = {-2212, 211, 321};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'E') {
        ReactionChannel = "AntiSexaquark,P->AntiLambda,K+,pi-,pi+";
        fProductsPDG = {-3122, 321, -211, 211};
        fFinalStateProductsPDG = {-2212, 211, 321, -211, 211};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'H') {
        ReactionChannel = "AntiSexaquark,P->AntiProton,K+,K+,pi0";
        fProductsPDG = {-2212, 321, 321, 111};
        fFinalStateProductsPDG = {321, 321};  // focused only on K+K+
        fStruckNucleonPDG = 2212;
    }

    AliInfo("Initializing...");
    AliInfo("Settings:");
    AliInfo("========");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);
    AliInfoF(">> Source of V0s       = \"%s\"", fSourceOfV0s.Data());
    AliInfoF(">> Reaction ID         = \'%c\'", fReactionID);
    AliInfoF(">> Reaction Channel    = \"%s\"", ReactionChannel.Data());
    AliInfoF(">> Stop After          = \"%s\"", fStopAfter.Data());
    AliInfoF(">> DoQA                = %i", (Int_t)fDoQA);
    AliInfoF(">> Read Signal Logs    = %i", (Int_t)fReadSignalLogs);
    AliInfoF(">> ReweightPt          = %i", (Int_t)fReweightPt);
    AliInfoF(">> ReweightRadius      = %i", (Int_t)fReweightRadius);
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
    kMax_Track_Chi2PerNTPCClusters = 2.;
    kTurnedOn_Track_StatusCuts = kTRUE;
    kTurnedOn_Track_RejectKinks = kFALSE;  // PENDING: need to study further
    kMin_Track_DCAxy_wrtPV = 2.;           // PENDING: could be improved for pions with DCAxy > 4 ??

    kMin_Track_Pt[2212] = 0.2;
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

    kMin_V0_Mass[310] = 0.475;
    kMax_V0_Mass[310] = 0.525;
    kMin_V0_Pt[310] = 1.0;

    kMin_V0_Mass[422] = 2 * fPDG.GetParticle(211)->Mass();
    kMax_V0_Mass[422] = 2.;
    kMin_V0_Pt[422] = 0.1;

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
        // kMax_V0_Chi2ndf[-3122] = 0.4;

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
        // kMax_V0_Chi2ndf[310] = 0.4;

        kMax_V0_Eta[422] = 0.9;
        kMin_V0_Radius[422] = 15.;
        kMin_V0_CPAwrtPV[422] = 0.8;
        kMax_V0_CPAwrtPV[422] = 0.98;
        kMin_V0_DCAwrtPV[422] = 4.;
        kMax_V0_DCAbtwDau[422] = 0.3;
        kMax_V0_DCAnegV0[422] = 0.2;
        kMax_V0_DCAposV0[422] = 0.2;
        // kMax_V0_Chi2ndf[422] = 0.1;
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
        if (fReactionID == 'A') {
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
            // kMax_Sexa_Chi2ndf = 0.5;
        } else if (fReactionID == 'D') {
            kMin_Sexa_Mass = 0.;
            kMin_Sexa_Radius = 30;
            kMax_Sexa_Rapidity = 0.8;
            kMin_Sexa_CPAwrtPV = 0.99;
            kMax_Sexa_CPAwrtPV = 1.;
            kMax_Sexa_DCAv0ba = 0.3;
            kMax_Sexa_DCAbaSV = 0.2;
            kMax_Sexa_DCAv0SV = 0.2;
            // kMax_Sexa_Chi2ndf = 1.;
        } else if (fReactionID == 'E') {
            kMin_Sexa_Radius = 30;
            kMin_Sexa_CPAwrtPV = 0.99;
            kMax_Sexa_CPAwrtPV = 1.;
            kMax_Sexa_DCAv0SV = 0.1;
            kMax_Sexa_DCApkSV = 0.1;
            kMax_Sexa_DCApmSV = 0.1;
            kMax_Sexa_DCAppSV = 0.1;
            kMax_Sexa_DCAv0pk = 0.1;
            kMax_Sexa_DCAv0pm = 0.1;
            kMax_Sexa_DCAv0pp = 0.1;
            kMax_Sexa_DCApkpm = 0.1;
            kMax_Sexa_DCApkpp = 0.1;
            // kMax_Sexa_Chi2ndf = 1.;
        }
    }
}

/*
 Define K+K+ pair selection cuts.
 - Uses: `fSourceOfV0s`
 - Input: `cuts_option`
*/
void AliAnalysisTaskSexaquark::DefinePosKaonPairCuts(TString cuts_option) {

    if (fSourceOfV0s == "kalman") {
        kMin_PosKaonPair_Radius = 20;
        kMin_PosKaonPair_Pt = 0.2;
        kMax_PosKaonPair_DCAbtwKaons = 0.2;
        kMax_PosKaonPair_DCApkaSV = 0.2;
        kMax_PosKaonPair_DCApkbSV = 0.2;
        // kMax_PosKaonPair_Chi2ndf = 0.1;
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
    Int_t mother_idx;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    AliMCParticle* mcNegDau;
    AliMCParticle* mcPosDau;

    Int_t n_daughters = 0;
    Int_t mcIdxNegDaughter = 0;
    Int_t mcIdxPosDaughter = 0;

    TVector3 originVertex;
    TVector3 decayVertex;

    Float_t decay_length;
    Float_t opening_angle;
    Float_t dca_wrt_pv, cpa_wrt_pv;
    Float_t armenteros_alpha, armenteros_qt;
    Bool_t reconstructable_v0;

    if (fDoQA) {
        fHist_QA["MCGen_VertexX"]->Fill(fMC_PrimaryVertex->GetX());
        fHist_QA["MCGen_VertexY"]->Fill(fMC_PrimaryVertex->GetY());
        fHist_QA["MCGen_VertexZ"]->Fill(fMC_PrimaryVertex->GetZ());
    }

    /* Loop over MC gen. particles in a single event */

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

        /* Protection against particles with near-zero momentum */

        if (mcPart->P() < 1E-6) continue;

        pdg_mc = mcPart->PdgCode();
        getPdgCode_fromMcIdx[mcIdx] = pdg_mc;

        /* Determine if particle is signal -- by checking if it's a reaction product */
        // NOTE:
        // If they're signal, their MC status number should be between 600 and 619.
        // This is done exclusively for the MC productions to the sexaquark analysis.
        // This was done to recognize each interaction, assigning a different number to each of
        // the 20 antisexaquark-nucleon reactions that were injected in a single event.
        //

        mother_idx = mcPart->GetMother();
        isMcIdxSignal[mcIdx] = mcPart->GetGeneratorIndex() == 2;

        if (isMcIdxSignal[mcIdx]) {
            getReactionIdx_fromMcIdx[mcIdx] = mother_idx < 0 ? mcPart->MCStatusCode() : getReactionIdx_fromMcIdx[mother_idx];
            getMcIdx_fromReactionIdx[getReactionIdx_fromMcIdx[mcIdx]].push_back(mcIdx);
        }

        /* Determine if it's a secondary particle */
        // NOTE:
        // Include signal particles, as well, because despite they were injected as primaries,
        // in reality, they are not
        //

        isMcIdxSecondary[mcIdx] = mcPart->IsSecondaryFromMaterial() || mcPart->IsSecondaryFromWeakDecay() || isMcIdxSignal[mcIdx];

        /** Charged particles: anti-protons, kaons, pions **/

        if (std::find(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end(), pdg_mc) != fFinalStateProductsPDG.end()) {

            /* Fill histograms */

            fHist_Tracks_Bookkeep[pdg_mc]->Fill(0);
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_Tracks[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());

            if (isMcIdxSecondary[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(10);
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }

            if (isMcIdxSignal[mcIdx]) {
                fHist_Tracks_Bookkeep[pdg_mc]->Fill(20);
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_Tracks[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            }
        }  // end of anti-protons, kaons, pions condition

        /** Neutral particles: anti-lambdas and K0S **/

        if (fReactionID != 'H' && (pdg_mc == -3122 || (fReactionID == 'A' && pdg_mc == 310))) {

            /* Fill histograms (reconstructable or not) */

            fHist_V0s_Bookkeep[pdg_mc]->Fill(0);
            if (isMcIdxSecondary[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(10);
            if (isMcIdxSignal[mcIdx]) fHist_V0s_Bookkeep[pdg_mc]->Fill(20);

            /* Check if this V0 can be reconstructed */

            // conditions:
            // (1) has daughters
            // (2) follows the following decay chains:
            //     -- anti-lambdas -> anti-proton pi+
            //     -- K0S -> pi+ pi-

            n_daughters = mcPart->GetNDaughters();
            mcIdxNegDaughter = 0;
            mcIdxPosDaughter = 0;

            // loop over daughters and fill mcIdxNegDaughter and mcIdxPosDaughter in case of relevant decays
            if (n_daughters) {
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
                }
            }  // end of n_daughters condition

            /* In any case, get the decay vertex of the V0, which is the origin vertex of its daughters */

            reconstructable_v0 = n_daughters && mcIdxNegDaughter && mcIdxPosDaughter;
            if (!reconstructable_v0) continue;

            /* Store true V0 and mother-daughter info */

            mcIndicesOfTrueV0s.push_back(mcIdx);

            doesMcIdxHaveMother[mcIdxNegDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxNegDaughter] = mcIdx;
            getNegDauMcIdx_fromMcIdx[mcIdx] = mcIdxNegDaughter;

            doesMcIdxHaveMother[mcIdxPosDaughter] = kTRUE;
            getMotherMcIdx_fromMcIdx[mcIdxPosDaughter] = mcIdx;
            getPosDauMcIdx_fromMcIdx[mcIdx] = mcIdxPosDaughter;

            /* Get particle's properties */

            mcNegDau = (AliMCParticle*)fMC->GetTrack(mcIdxNegDaughter);
            mcPosDau = (AliMCParticle*)fMC->GetTrack(mcIdxPosDaughter);

            armenteros_qt = ArmenterosQt(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            armenteros_alpha = ArmenterosAlpha(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz(),
                                               mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());

            decayVertex.SetXYZ(mcNegDau->Xv(), mcNegDau->Yv(), mcNegDau->Zv());
            cpa_wrt_pv = CosinePointingAngle(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                             fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                      fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

            originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
            decay_length = (originVertex - decayVertex).Mag();

            TVector3 neg(mcNegDau->Px(), mcNegDau->Py(), mcNegDau->Pz());
            TVector3 pos(mcPosDau->Px(), mcPosDau->Py(), mcPosDau->Pz());
            opening_angle = neg.Angle(pos);

            /* Fill histograms (reconstructables only) */

            fHist_V0s_Bookkeep[pdg_mc]->Fill(1);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Mass")]->Fill(mcPart->M());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pt")]->Fill(mcPart->Pt());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Px")]->Fill(mcPart->Px());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Py")]->Fill(mcPart->Py());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Pz")]->Fill(mcPart->Pz());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Phi")]->Fill(mcPart->Phi());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Zv")]->Fill(decayVertex.Z());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Eta")]->Fill(mcPart->Eta());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("MCGen", "All", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (isMcIdxSecondary[mcIdx]) {
                fHist_V0s_Bookkeep[pdg_mc]->Fill(11);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Phi")]->Fill(mcPart->Phi());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Secondary", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
            }

            if (isMcIdxSignal[mcIdx]) {
                fHist_V0s_Bookkeep[pdg_mc]->Fill(21);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Mass")]->Fill(mcPart->M());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pt")]->Fill(mcPart->Pt());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Px")]->Fill(mcPart->Px());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Py")]->Fill(mcPart->Py());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Pz")]->Fill(mcPart->Pz());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Phi")]->Fill(mcPart->Phi());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "OriginRadius")]->Fill(originVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayRadius")]->Fill(decayVertex.Perp());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Zv")]->Fill(decayVertex.Z());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Eta")]->Fill(mcPart->Eta());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "Rapidity")]->Fill(mcPart->Y());
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmQt")]->Fill(armenteros_qt);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "ArmAlpha")]->Fill(armenteros_alpha);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("MCGen", "Signal", pdg_mc, "DCAwrtPV")]->Fill(dca_wrt_pv);
            }
        }  // end of anti-lambda && K0s condition
    }  // end of loop over MC particles
}

/*
 Loop over generated antisexaquark-nucleon reactions.
 - Uses: `fMC`, `fPDG`, `fStruckNucleonPDG`
 - Containers: `getMcIdx_fromReactionIdx`
*/
void AliAnalysisTaskSexaquark::ProcessSignalInteractions() {

    AliMCParticle* mcProduct;

    /** Declare TLorentzVectors **/

    TLorentzVector lvProduct;
    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /** Declare variables **/

    TVector3 secondary_vertex(0., 0., 0.);
    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    Double_t injected_px;
    Double_t injected_py;
    Double_t injected_pz;
    Double_t injected_m;
    if (fReadSignalLogs) {
        fLogTree->SetBranchAddress("Sexa_Px_ini", &injected_px);
        fLogTree->SetBranchAddress("Sexa_Py_ini", &injected_py);
        fLogTree->SetBranchAddress("Sexa_Pz_ini", &injected_pz);
        fLogTree->SetBranchAddress("Sexa_M_ini", &injected_m);
    }
    TLorentzVector lvInjected;

    /* Loop over reactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /** Protection in case of general-purpose simulations **/

        if (!getMcIdx_fromReactionIdx[reactionIdx].size()) continue;

        /** Reset vectors **/

        secondary_vertex.SetXYZ(0., 0., 0.);
        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /** Loop over reaction products **/

        for (Int_t& mcIdx : getMcIdx_fromReactionIdx[reactionIdx]) {

            mcProduct = (AliMCParticle*)fMC->GetTrack(mcIdx);

            /* Get first generation products, the direct daughters of the reaction */

            if (mcProduct->MCStatusCode() != reactionIdx) continue;

            lvProduct.SetPxPyPzE(mcProduct->Px(), mcProduct->Py(), mcProduct->Pz(), mcProduct->E());

            lvAntiSexaquark = lvAntiSexaquark + lvProduct;

            /* Get sec. vertex (if it hasn't been defined yet)*/

            if (secondary_vertex.Mag() < 1E-6) secondary_vertex.SetXYZ(mcProduct->Xv(), mcProduct->Yv(), mcProduct->Zv());
        }

        lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest
        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /** Get properties **/

        radius = secondary_vertex.Perp();
        getRadius_fromReactionIdx[reactionIdx] = radius;

        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), secondary_vertex.X(), secondary_vertex.Y(),
                                  secondary_vertex.Z(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /** Fill histograms **/

        /* With information A.I. -- After Interaction  */

        fHist_AntiSexaquarks_Bookkeep->Fill(0);
        if (fReactionID == 'H') fHist_PosKaonPairs_Bookkeep->Fill(0);
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Px")]->Fill(lvAntiSexaquark.Px());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Py")]->Fill(lvAntiSexaquark.Py());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Radius")]->Fill(radius);  // interaction radius
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Zv")]->Fill(mcProduct->Zv());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("MCGen_AI", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);

        /* With information B.I. -- Before Interaction */

        if (!fReadSignalLogs) continue;

        fLogTree->GetEntry(reactionIdx - 600);

        lvInjected.SetXYZM(injected_px, injected_py, injected_pz, injected_m);
        getPt_fromReactionIdx[reactionIdx] = lvInjected.Pt();

        cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), secondary_vertex.X(), secondary_vertex.Y(),
                                  secondary_vertex.Z(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Mass")]->Fill(lvInjected.M());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Pt")]->Fill(lvInjected.Pt());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Px")]->Fill(lvInjected.Px());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Py")]->Fill(lvInjected.Py());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Pz")]->Fill(lvInjected.Pz());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Phi")]->Fill(lvInjected.Phi());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Eta")]->Fill(lvInjected.Eta());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "Rapidity")]->Fill(lvInjected.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("MCGen_BI", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }  // end of loop over reactions
}

/*                   */
/**  Reconstructed  **/
/*** ============= ***/

/*
 Apply event selection
*/
Bool_t AliAnalysisTaskSexaquark::PassesEventSelection() {

    if (fDoQA) fHist_QA["Events_Bookkeep"]->Fill(0);

    /* MC PV z-vertex range */

    if (fIsMC && TMath::Abs(fMC_PrimaryVertex->GetZ()) > 10.) return kFALSE;
    AliInfo("!! MC Event Passed z-vertex range on PV !!");  // DEBUG

    if (fDoQA) fHist_QA["Events_Bookkeep"]->Fill(1);

    // ! REFERENCE: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGRunList18r ! //
    if (!fIsMC && (fRunNumber == 296749 || fRunNumber == 296750 || fRunNumber == 296849 || fRunNumber == 296890 || fRunNumber == 297029 ||
                   fRunNumber == 297194 || fRunNumber == 297219 || fRunNumber == 297481)) {
        fEventCuts.UseTimeRangeCut();
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny);
        if (!fEventCuts.AcceptEvent(fESD)) return kFALSE;
        AliInfo("!! MC Event Passed z-vertex range on PV !!");  // DEBUG
    }

    /* Trigger Selection */

    Bool_t MB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7);  // kAnyINT
    if (!MB) {
        AliInfoF("!! Event Rejected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);  // DEBUG
        return kFALSE;
    }
    AliInfoF("!! Event Selected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);  // DEBUG

    if (fDoQA) fHist_QA["Events_Bookkeep"]->Fill(2);

    /* rec. PV z-vertex range */

    if (TMath::Abs(fPrimaryVertex->GetZ()) > 10.) return kFALSE;
    AliInfo("!! Event Passed z-vertex range on PV !!");  // DEBUG

    if (fDoQA) {
        fHist_QA["Events_Bookkeep"]->Fill(3);

        /* Vertex */

        fHist_QA["VertexX"]->Fill(fPrimaryVertex->GetX());
        fHist_QA["VertexY"]->Fill(fPrimaryVertex->GetY());
        fHist_QA["VertexZ"]->Fill(fPrimaryVertex->GetZ());

        /* N Tracks */

        Int_t ntracks = fESD->GetNumberOfTracks();
        fHist_QA["NTracks"]->Fill(ntracks);

        /* SPD */

        const AliESDVertex* spd_vertex = fESD->GetPrimaryVertexSPD();
        const AliMultiplicity* mult = fESD->GetMultiplicity();

        fHist_QA["SPD_NTracklets"]->Fill(mult->GetNumberOfTracklets());

        f2DHist_QA["SPD_MultiplicityVsVertexZ"]->Fill(spd_vertex->GetZ(), mult->GetNumberOfTracklets());
        f2DHist_QA["SPD_NTrackletsClu0VsNTracklets"]->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(0));
        f2DHist_QA["SPD_NTrackletsClu1VsNTracklets"]->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(1));

        for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++) {
            Float_t phiTr = mult->GetPhi(iTracklet);
            Float_t etaTr = mult->GetEta(iTracklet);
            f2DHist_QA["SPD_PhiVsEta"]->Fill(etaTr, phiTr);
        }

        /* TPC */

        fHist_QA["TPC_NClusters"]->Fill(fESD->GetNumberOfTPCClusters());
    }

    return kTRUE;
}

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisTaskSexaquark::ProcessTracks() {

    AliESDtrack* track;

    Int_t mcIdx;
    Int_t mcPdgCode;
    Int_t motherMcIdx;

    Float_t dca_wrt_pv, dca_xy_wrt_pv, dca_z_wrt_pv;
    Float_t n_tpc_clusters;
    Float_t n_crossed_rows;
    Float_t chi2_over_nclusters;
    Float_t n_findable_clusters, n_shared_clusters;
    Float_t chi2_over_ncls;

    Float_t n_sigma_proton;
    Float_t n_sigma_kaon;
    Float_t n_sigma_pion;

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Bool_t already_selected;  // as another species?
    UInt_t nSelectedTracks = 0;

    /** Loop over tracks in a single event **/

    for (Int_t esdIdxTrack = 0; esdIdxTrack < fESD->GetNumberOfTracks(); esdIdxTrack++) {

        /* Get track */

        track = fESD->GetTrack(esdIdxTrack);

        /* Get MC info */

        mcIdx = TMath::Abs(track->GetLabel());
        getMcIdx_fromEsdIdx[esdIdxTrack] = mcIdx;
        mcPdgCode = getPdgCode_fromMcIdx[mcIdx];

        /* Fill histograms before track selection */

        f2DHist_Tracks_TPCsignal["All"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        /* QA : Tracks */

        if (fDoQA) {
            const Int_t NITSLayers = 6;
            for (Int_t iLayer = 0; iLayer < NITSLayers; iLayer++) {
                if (track->HasPointOnITSLayer(iLayer)) {
                    f2DHist_QA["ITS_LayerNoVsPhi"]->Fill(track->Phi(), iLayer);
                    f2DHist_QA["ITS_LayerNoVsEta"]->Fill(track->Eta(), iLayer);
                }
            }
            if (TMath::Abs(track->Eta()) < kMax_Track_Eta) {
                track->GetImpactParameters(dca_xy_wrt_pv, dca_z_wrt_pv);  // pre-calculated DCA w.r.t. PV
                fHist_QA["Tracks_DCAxy"]->Fill(dca_xy_wrt_pv);
                fHist_QA["Tracks_DCAz"]->Fill(dca_z_wrt_pv);
                fHist_QA["Tracks_Pt"]->Fill(track->Pt());
                f2DHist_QA["Tracks_PhiVsEta"]->Fill(track->Eta(), track->Phi());
            }
        }  // end of QA

        /* -- DEBUG -- */

        // clang-format off
        if (track->P() > 1E3) {
            Float_t xy_impar, z_impar;
            track->GetImpactParameters(xy_impar, z_impar);
            AliInfoF("idx px py pz true pdg = %d, %f, %f, %f, %d, %d", esdIdxTrack, track->Px(), track->Py(), track->Pz(), mcIdx, mcPdgCode);
            AliInfo("nsigmapro, nsigmak, nsigmapi, eta, nclusters, chi2/nclusters, kink, dcaxy, dcaz, ncrossedrows, sharedcl, sharedcl/cRows, sharedcl/ncl =");
            AliInfoF(" %f, %f, %f, %f, %i, %f, %i, %f, %f, %f, %i, %f, %f",
                     fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton),
                     fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon),
                     fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion),
                     track->Eta(),
                     (Int_t)track->GetTPCNcls(),
                     track->GetTPCchi2() / (Double_t)track->GetTPCNcls(),
                     (Int_t)track->GetKinkIndex(0),
                     xy_impar,
                     z_impar,
                     track->GetTPCCrossedRows(),
                     (Int_t)track->GetTPCnclsS(),
                     track->GetTPCCrossedRows() > 1E-6 ? (Float_t)track->GetTPCnclsS() / track->GetTPCCrossedRows() : 999.,
                     (Float_t)track->GetTPCNcls() > 1E-6 ? (Float_t)track->GetTPCnclsS() / (Float_t)track->GetTPCNcls() : 999.);
            continue;
        }
        // clang-format on

        /* -- END OF DEBUG -- */

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

        n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

        track->GetImpactParameters(dca_xy_wrt_pv, dca_z_wrt_pv);  // pre-calculated DCA w.r.t. PV
        dca_xy_wrt_pv = TMath::Abs(dca_xy_wrt_pv);
        dca_z_wrt_pv = TMath::Abs(dca_z_wrt_pv);
        dca_wrt_pv = TMath::Sqrt(dca_xy_wrt_pv * dca_xy_wrt_pv + dca_z_wrt_pv * dca_z_wrt_pv);
        n_tpc_clusters = track->GetTPCNcls();  // NOTE: capital N
        n_crossed_rows = track->GetTPCCrossedRows();
        n_findable_clusters = (Int_t)track->GetTPCNclsF();
        n_shared_clusters = (Int_t)track->GetTPCnclsS();
        chi2_over_ncls = n_tpc_clusters > 1E-6 ? track->GetTPCchi2() / (Double_t)n_tpc_clusters : 999.;

        /* Loop over PID hypotheses */

        // get only the tracks that have a relevant PID for the chosen channel
        // (note: fFinalStateProductsPDG is copied into a set to remove duplicates)
        std::set<Int_t> possiblePID(fFinalStateProductsPDG.begin(), fFinalStateProductsPDG.end());
        already_selected = kFALSE;

        for (Int_t esdPdgCode : possiblePID) {

            /* Apply track selection */

            is_true = mcPdgCode == esdPdgCode;
            is_secondary = isMcIdxSecondary[mcIdx] && is_true;
            is_signal = isEsdIdxSignal[esdIdxTrack] && is_secondary;

            if (!PassesTrackSelectionAs(esdPdgCode, is_true, is_secondary, is_signal, track)) continue;

            if (esdPdgCode == 2212) esdIdxPassesTrackSelection_AsProton[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -2212) esdIdxPassesTrackSelection_AsAntiProton[esdIdxTrack] = kTRUE;
            if (esdPdgCode == 321) esdIdxPassesTrackSelection_AsPosKaon[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -321) esdIdxPassesTrackSelection_AsNegKaon[esdIdxTrack] = kTRUE;
            if (esdPdgCode == 211) esdIdxPassesTrackSelection_AsPiPlus[esdIdxTrack] = kTRUE;
            if (esdPdgCode == -211) esdIdxPassesTrackSelection_AsPiMinus[esdIdxTrack] = kTRUE;

            if (!already_selected) {
                getEsdIdx_fromMcIdx[mcIdx] = esdIdxTrack;  // purposefully done after selection
                f2DHist_Tracks_TPCsignal["Selected"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
                if (fDoQA) nSelectedTracks++;
                already_selected = kTRUE;
            }

            if (fDoQA) {
                if (esdPdgCode == -2212) {
                    f2DHist_QA["TPC_NSigmaProtonVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_proton);
                    f2DHist_QA["TPC_NSigmaProtonVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_proton);
                }
                if (esdPdgCode == 321) {
                    f2DHist_QA["TPC_NSigmaKaonVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_kaon);
                    f2DHist_QA["TPC_NSigmaKaonVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_kaon);
                }
                if (TMath::Abs(esdPdgCode) == 211) {
                    f2DHist_QA["TPC_NSigmaPionVsInnerParamP"]->Fill(track->GetInnerParam()->P(), n_sigma_pion);
                    f2DHist_QA["TPC_NSigmaPionVsEta"]->Fill(track->GetInnerParam()->Eta(), n_sigma_pion);
                }
            }

            if (esdPdgCode == -2212) esdIndicesOfAntiProtonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 321) esdIndicesOfPosKaonTracks.push_back(esdIdxTrack);
            if (esdPdgCode == -211) esdIndicesOfPiMinusTracks.push_back(esdIdxTrack);
            if (esdPdgCode == 211) esdIndicesOfPiPlusTracks.push_back(esdIdxTrack);

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(50);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "All", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "All", esdPdgCode);

            if (!is_true) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(80);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "True", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "True", esdPdgCode);
            f2DHist_Tracks_TPCsignal["True"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!is_secondary) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(110);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "Secondary", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "Secondary", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Secondary"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!is_signal) continue;

            fHist_Tracks_Bookkeep[esdPdgCode]->Fill(140);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCA_wrtPV")]->Fill(dca_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Found", "Signal", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            PlotStatus(track, "Signal", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Signal"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
        }  // end of loop over PID hypotheses
    }  // end of loop over tracks

    if (fDoQA) fHist_QA["NSelectedTracks"]->Fill(nSelectedTracks);
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesTrackSelectionAs(Int_t HypothesisPID, Bool_t is_true, Bool_t is_secondary, Bool_t is_signal,
                                                        AliESDtrack* track) {

    AliPID::EParticleType ParticleType;
    Double_t kMax_NSigma;
    if (TMath::Abs(HypothesisPID) == 2212) {
        ParticleType = AliPID::kProton;
        kMax_NSigma = kMax_NSigma_Proton;
    }
    if (TMath::Abs(HypothesisPID) == 321) {
        ParticleType = AliPID::kKaon;
        kMax_NSigma = kMax_NSigma_Kaon;
    }
    if (TMath::Abs(HypothesisPID) == 211) {
        ParticleType = AliPID::kPion;
        kMax_NSigma = kMax_NSigma_Pion;
    }

    Double_t n_sigma = fPIDResponse->NumberOfSigmasTPC(track, ParticleType);
    if (TMath::Abs(n_sigma) > kMax_NSigma) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(30);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(60);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(90);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(120);

    Double_t Eta = TMath::Abs(track->Eta());
    if (kMax_Track_Eta && Eta > kMax_Track_Eta) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(31);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(61);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(91);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(121);

    Double_t NTPCClusters = track->GetTPCNcls();
    if (kMin_Track_NTPCClusters && NTPCClusters < kMin_Track_NTPCClusters) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(32);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(62);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(92);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(122);

    Double_t Chi2PerNTPCClusters = NTPCClusters > 1E-6 ? track->GetTPCchi2() / (Double_t)track->GetTPCNcls() : 999;
    if (kMax_Track_Chi2PerNTPCClusters && Chi2PerNTPCClusters > kMax_Track_Chi2PerNTPCClusters) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(33);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(63);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(93);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(123);

    // >> TPC and ITS status
    Bool_t tpc_status =
        (track->GetStatus() & AliESDtrack::kTPCin) && (track->GetStatus() & AliESDtrack::kTPCout) && (track->GetStatus() & AliESDtrack::kTPCrefit);
    if (kTurnedOn_Track_StatusCuts && !tpc_status) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(34);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(64);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(94);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(124);

    Bool_t its_status =
        !(track->GetStatus() & AliESDtrack::kITSin) && !(track->GetStatus() & AliESDtrack::kITSout) && !(track->GetStatus() & AliESDtrack::kITSrefit);
    if (kTurnedOn_Track_StatusCuts && !its_status) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(35);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(65);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(95);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(125);

    // >> reject kinks
    if (kTurnedOn_Track_RejectKinks && track->GetKinkIndex(0) != 0) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(36);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(66);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(96);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(126);

    // >> DCA w.r.t Primary Vertex
    Float_t DCAxy_wrtPV, DCAz_wrtPV;
    track->GetImpactParameters(DCAxy_wrtPV, DCAz_wrtPV);
    Float_t DCA_wrtPV = TMath::Sqrt((Double_t)DCAxy_wrtPV * (Double_t)DCAxy_wrtPV + (Double_t)DCAz_wrtPV * (Double_t)DCAz_wrtPV);
    if (kMin_Track_DCA_wrtPV && DCA_wrtPV < kMin_Track_DCA_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(37);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(67);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(97);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(127);

    if (kMin_Track_DCAxy_wrtPV && TMath::Abs(DCAxy_wrtPV) < kMin_Track_DCAxy_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(38);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(68);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(98);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(128);

    if (kMin_Track_DCAz_wrtPV && TMath::Abs(DCAz_wrtPV) < kMin_Track_DCAz_wrtPV) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(39);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(69);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(99);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(129);

    // >> pT threshold
    if (kMin_Track_Pt[HypothesisPID] && track->Pt() < kMin_Track_Pt[HypothesisPID]) return kFALSE;
    fHist_Tracks_Bookkeep[HypothesisPID]->Fill(40);
    if (is_true) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(70);
    if (is_secondary) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(100);
    if (is_signal) fHist_Tracks_Bookkeep[HypothesisPID]->Fill(130);

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

/*                                                                   */
/**  V0s, Anti-Sexaquarks and K+K+ Pairs -- That Require True Info  **/
/*** ============================================================= ***/

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

    TVector3 neg, pos;

    Int_t pdg_code;
    Float_t origin_radius, decay_radius;
    Float_t decay_length;
    Float_t opening_angle;
    Float_t arm_qt, arm_alpha;
    Float_t cpa_wrt_pv, dca_wrt_pv;

    /* Loop over true V0s */

    for (Int_t& mcIdxOfTrueV0 : mcIndicesOfTrueV0s) {

        mcV0 = (AliMCParticle*)fMC->GetTrack(mcIdxOfTrueV0);

        /* Check if both daughters were reconstructed */

        if (!getNegDauEsdIdx_fromMcIdx[mcIdxOfTrueV0] || !getPosDauEsdIdx_fromMcIdx[mcIdxOfTrueV0]) {
            continue;
        }

        is_signal = mcV0->MCStatusCode() >= 600 && mcV0->MCStatusCode() < 620;
        is_secondary = mcV0->IsSecondaryFromMaterial() || mcV0->IsSecondaryFromWeakDecay() || is_signal;

        /* Get daughters' properties */

        mcNegDaughter = (AliMCParticle*)fMC->GetTrack(getNegDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);
        mcPosDaughter = (AliMCParticle*)fMC->GetTrack(getPosDauMcIdx_fromMcIdx[mcIdxOfTrueV0]);

        neg.SetXYZ(mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        pos.SetXYZ(mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());

        /* Get V0 properties */

        pdg_code = mcV0->PdgCode();
        origin_radius = TMath::Sqrt(mcV0->Xv() * mcV0->Xv() + mcV0->Yv() * mcV0->Yv());
        decay_radius = TMath::Sqrt(mcNegDaughter->Xv() * mcNegDaughter->Xv() + mcNegDaughter->Yv() * mcNegDaughter->Yv());
        decay_length = TMath::Sqrt(TMath::Power(mcNegDaughter->Xv() - mcV0->Xv(), 2) + TMath::Power(mcNegDaughter->Yv() - mcV0->Yv(), 2) +
                                   TMath::Power(mcNegDaughter->Zv() - mcV0->Zv(), 2));
        opening_angle = neg.Angle(pos);
        arm_qt = ArmenterosQt(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz());
        arm_alpha = ArmenterosAlpha(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Px(), mcNegDaughter->Py(), mcNegDaughter->Pz(),
                                    mcPosDaughter->Px(), mcPosDaughter->Py(), mcPosDaughter->Pz());
        cpa_wrt_pv = CosinePointingAngle(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                         fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill bookkeeping histograms */

        if (pdg_code == -3122 || pdg_code == 310) {

            fHist_V0s_Bookkeep[pdg_code]->Fill(30);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "True", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(40);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Secondary", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdg_code]->Fill(50);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Mass")]->Fill(mcV0->M());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pt")]->Fill(mcV0->Pt());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Px")]->Fill(mcV0->Px());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Py")]->Fill(mcV0->Py());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Pz")]->Fill(mcV0->Pz());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Phi")]->Fill(mcV0->Phi());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "OriginRadius")]->Fill(origin_radius);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Zv")]->Fill(mcNegDaughter->Zv());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Eta")]->Fill(mcV0->Eta());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "Rapidity")]->Fill(mcV0->Y());
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Findable", "Signal", pdg_code, "DCAwrtPV")]->Fill(dca_wrt_pv);
        }
    }  // end of loop over true V0s
}

/*
  Find the reconstructed signal pi+ (that don't come from the anti-lambda) and pi-
  (only valid for Reaction Channel E: anti-lambda, K+, pi+, pi-)
*/
void AliAnalysisTaskSexaquark::ProcessFindablePionPairs() {

    Int_t single_pion_pdg_code;
    Int_t single_pion_mother_pdg_code;

    Int_t mcIdxPiPlus = -1;
    Int_t mcIdxPiMinus = -1;
    AliMCParticle* mcPiPlus;
    AliMCParticle* mcPiMinus;

    /* Declare TLorentzVectors and TVector3 */

    TLorentzVector lvPiPlus, lvPiMinus;
    TLorentzVector lvPionPair;
    TVector3 secondary_vertex(0., 0., 0.);

    /* Declare properties */

    Float_t origin_radius;
    Float_t decay_length;
    Float_t opening_angle;
    Float_t arm_qt;
    Float_t arm_alpha;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection against general-purpose simulations */
        // There should be tracks that belong to the interaction //

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Collect true info from rec. signal tracks */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            single_pion_pdg_code = getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]];
            single_pion_mother_pdg_code = 0;
            if (doesEsdIdxHaveMother[esdIdx]) single_pion_mother_pdg_code = getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdx]];

            if (TMath::Abs(single_pion_pdg_code) != 211) continue;
            if (single_pion_mother_pdg_code == -3122) continue;

            // if it hasn't been set yet, set it
            if (mcIdxPiPlus == -1) mcIdxPiPlus = getMcIdx_fromEsdIdx[esdIdx];
            if (mcIdxPiMinus == -1) mcIdxPiMinus = getMcIdx_fromEsdIdx[esdIdx];
        }

        /* Findable condition */
        // At least 1 pi+ and 1 pi- should have been reconstructed and passed track selection

        if (mcIdxPiPlus < 0 || mcIdxPiMinus < 0) continue;

        /* Get MC particles */

        mcPiPlus = (AliMCParticle*)fMC->GetTrack(mcIdxPiPlus);
        mcPiMinus = (AliMCParticle*)fMC->GetTrack(mcIdxPiMinus);

        /* Get true momentum */

        lvPiPlus.SetXYZM(mcPiPlus->Px(), mcPiPlus->Py(), mcPiPlus->Pz(), fPDG.GetParticle(211)->Mass());
        lvPiMinus.SetXYZM(mcPiMinus->Px(), mcPiMinus->Py(), mcPiMinus->Pz(), fPDG.GetParticle(211)->Mass());

        lvPionPair = lvPiPlus + lvPiMinus;

        /* Get sec. vertex */

        secondary_vertex.SetXYZ(mcPiPlus->Xv(), mcPiPlus->Yv(), mcPiPlus->Zv());

        /* Get properties */

        origin_radius = TMath::Sqrt(TMath::Power(secondary_vertex.X(), 2) + TMath::Power(secondary_vertex.Y(), 2));
        decay_length = TMath::Sqrt(TMath::Power(secondary_vertex.X() - fMC_PrimaryVertex->GetX(), 2) +
                                   TMath::Power(secondary_vertex.Y() - fMC_PrimaryVertex->GetY(), 2) +
                                   TMath::Power(secondary_vertex.Z() - fMC_PrimaryVertex->GetZ(), 2));
        opening_angle = lvPionPair.Vect().Angle(lvPiMinus.Vect());
        arm_qt = ArmenterosQt(lvPionPair.Px(), lvPionPair.Py(), lvPionPair.Pz(), lvPiMinus.Px(), lvPiMinus.Py(), lvPiMinus.Pz());
        arm_alpha = ArmenterosAlpha(lvPionPair.Px(), lvPionPair.Py(), lvPionPair.Pz(), lvPiMinus.Px(), lvPiMinus.Py(), lvPiMinus.Pz(), lvPiPlus.Px(),
                                    lvPiPlus.Py(), lvPiPlus.Pz());
        cpa_wrt_pv = CosinePointingAngle(lvPiPlus, secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(), fMC_PrimaryVertex->GetX(),
                                         fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(lvPiPlus.Px(), lvPiPlus.Py(), lvPiPlus.Pz(), secondary_vertex.X(), secondary_vertex.Y(), secondary_vertex.Z(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill histograms  */

        fHist_V0s_Bookkeep[422]->Fill(50);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Mass")]->Fill(lvPionPair.M());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Pt")]->Fill(lvPionPair.Pt());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Px")]->Fill(lvPionPair.Px());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Py")]->Fill(lvPionPair.Py());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Pz")]->Fill(lvPionPair.Pz());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Phi")]->Fill(lvPionPair.Phi());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "OriginRadius")]->Fill(origin_radius);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DecayRadius")]->Fill(origin_radius);  // NOTE: is the same on purpose
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Zv")]->Fill(secondary_vertex.Z());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Eta")]->Fill(lvPionPair.Eta());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "Rapidity")]->Fill(lvPionPair.Rapidity());
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DecayLength")]->Fill(decay_length);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "OpeningAngle")]->Fill(opening_angle);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "ArmQt")]->Fill(arm_qt);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "ArmAlpha")]->Fill(arm_alpha);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_V0s[std::make_tuple("Findable", "Signal", 422, "DCAwrtPV")]->Fill(dca_wrt_pv);
    }
}

/*
  Find the anti-sexaquarks that have all of their *final-state particles* reconstructed. Valid for:
  - Reaction Channel A: 1 anti-proton, 2 pi+, 1 pi-
  - Reaction Channel D: 1 anti-proton, 1 pi+, 1 K+
  - Reaction Channel E: 1 anti-proton, 2 pi+, 1 pi-, 1 K+
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

        fHist_AntiSexaquarks_Bookkeep->Fill(10);
        if (fReactionID == 'H') fHist_PosKaonPairs_Bookkeep->Fill(1);  // PENDING: should be ~ 0, because it's a channel so difficult to rec.
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Mass")]->Fill(lvAntiSexaquark.M());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Px")]->Fill(lvAntiSexaquark.Px());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Py")]->Fill(lvAntiSexaquark.Py());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Radius")]->Fill(radius);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Zv")]->Fill(secondary_vertex.Z());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DecayLength")]->Fill(decay_length);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
        fHist_AntiSexaquarks[std::make_tuple("Findable", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
    }
}

/*
  Find the anti-sexaquarks that have two K+s reconstructed.
  Valid for Reaction Channel H: 1 anti-proton, 2 K+, 1 pi0
*/
void AliAnalysisTaskSexaquark::ProcessFindablePosKaonPairs() {

    Int_t pdgCode;
    Int_t mcIdxPosKaonA = -1;
    Int_t mcIdxPosKaonB = -1;
    AliMCParticle* mcPosKaonA;
    AliMCParticle* mcPosKaonB;
    Int_t motherMcIdx;

    /* Declare TLorentzVectors */

    TLorentzVector lvPosKaonA, lvPosKaonB;
    TLorentzVector lvPosKaonPair;

    /* Declare vertices */

    TVector3 primary_vertex(fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
    TVector3 secondary_vertex(0., 0., 0.);

    /* Loop over interactions */

    for (UInt_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Protection in case of general-purpose simulations */

        if (!getEsdIdx_fromReactionIdx[reactionIdx].size()) continue;

        /* Collect true info from rec. signal tracks */

        for (Int_t& esdIdx : getEsdIdx_fromReactionIdx[reactionIdx]) {

            pdgCode = getPdgCode_fromMcIdx[getMcIdx_fromEsdIdx[esdIdx]];

            if (pdgCode != 321) continue;

            if (mcIdxPosKaonA == -1)
                mcIdxPosKaonA = getMcIdx_fromEsdIdx[esdIdx];
            else if (mcIdxPosKaonB == -1)
                mcIdxPosKaonB = getMcIdx_fromEsdIdx[esdIdx];
        }

        /* Findable condition */
        // Two positive kaons should have been reconstructed and passed track selection

        if (mcIdxPosKaonA < 0 || mcIdxPosKaonB < 0) continue;

        /* Get MC particles */

        mcPosKaonA = (AliMCParticle*)fMC->GetTrack(mcIdxPosKaonA);
        mcPosKaonB = (AliMCParticle*)fMC->GetTrack(mcIdxPosKaonB);

        /* Get true momentum */

        lvPosKaonA.SetXYZM(mcPosKaonA->Px(), mcPosKaonA->Py(), mcPosKaonA->Pz(), fPDG.GetParticle(321)->Mass());
        lvPosKaonB.SetXYZM(mcPosKaonB->Px(), mcPosKaonB->Py(), mcPosKaonB->Pz(), fPDG.GetParticle(321)->Mass());

        lvPosKaonPair = lvPosKaonA + lvPosKaonB;

        /* Get sec. vertex */

        secondary_vertex.SetXYZ(mcPosKaonA->Xv(), mcPosKaonA->Yv(), mcPosKaonA->Zv());
    }

    /* Fill histograms  */

    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Mass")]->Fill(lvPosKaonPair.M());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "OriginRadius")]->Fill(secondary_vertex.Perp());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "DecayLength")]->Fill((secondary_vertex - primary_vertex).Mag());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Xv")]->Fill(secondary_vertex.X());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Yv")]->Fill(secondary_vertex.Y());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Zv")]->Fill(secondary_vertex.Z());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Eta")]->Fill(lvPosKaonPair.Eta());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Px")]->Fill(lvPosKaonPair.Px());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Py")]->Fill(lvPosKaonPair.Py());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Pz")]->Fill(lvPosKaonPair.Pz());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "Pt")]->Fill(lvPosKaonPair.Pt());
    fHist_PosKaonPairs[std::make_tuple("Findable", "All", "OpeningAngle")]->Fill(lvPosKaonA.Vect().Angle(lvPosKaonB.Vect()));
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

    Float_t decay_radius;
    Float_t decay_length;  // decay vertex - primary vertex, could be useful
    Float_t opening_angle;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_btw_dau, dca_neg_v0, dca_pos_v0;

    AliESDv0 outputV0;

    Double_t xthiss, xpp;
    Float_t impar_neg_v0[2], impar_pos_v0[2];

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Auxiliary containers */

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Bool_t passes_cuts_as_antilambda;
    Bool_t passes_cuts_as_kaonzero;

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

        if (!esdIdxPassesTrackSelection_AsAntiProton[esdIdxNeg] && !esdIdxPassesTrackSelection_AsNegKaon[esdIdxNeg] &&
            !esdIdxPassesTrackSelection_AsPiMinus[esdIdxNeg]) {
            continue;
        }

        if (!esdIdxPassesTrackSelection_AsProton[esdIdxPos] && !esdIdxPassesTrackSelection_AsPosKaon[esdIdxPos] &&
            !esdIdxPassesTrackSelection_AsPiPlus[esdIdxPos]) {
            continue;
        }

        /* Get tracks */

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

        /* Apply cuts */

        is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                  getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                  getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == -3122;
        is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
        is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

        passes_cuts_as_antilambda = PassesV0CutsAs(-3122, is_true, is_secondary, is_signal, V0, neg_track, pos_track);

        is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                  getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                  getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == 310;
        is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
        is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

        passes_cuts_as_kaonzero = PassesV0CutsAs(310, is_true, is_secondary, is_signal, V0, neg_track, pos_track);

        if (!passes_cuts_as_antilambda && !passes_cuts_as_kaonzero) {
            continue;
        }

        /* Get properties */

        V0->GetXYZ(V0_X, V0_Y, V0_Z);
        V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
        V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

        neg_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_neg_v0);
        pos_track->GetDZ(V0_X, V0_Y, V0_Z, fMagneticField, impar_pos_v0);

        decay_radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
        decay_length = TMath::Sqrt((V0_X - fPrimaryVertex->GetX()) * (V0_X - fPrimaryVertex->GetX()) +
                                   (V0_Y - fPrimaryVertex->GetY()) * (V0_Y - fPrimaryVertex->GetY()) +
                                   (V0_Z - fPrimaryVertex->GetZ()) * (V0_Z - fPrimaryVertex->GetZ()));

        cpa_wrt_pv = V0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_wrt_pv =
            LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_btw_dau = TMath::Abs(neg_track->GetDCA(pos_track, fMagneticField, xthiss, xpp));              // DCAxyz
        dca_neg_v0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);  // DCAxyz
        dca_pos_v0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);  // DCAxyz

        /* Fill histograms */

        std::vector<Int_t> possiblePID;
        if (passes_cuts_as_antilambda) possiblePID.push_back(-3122);
        if (passes_cuts_as_kaonzero) possiblePID.push_back(310);

        for (Int_t& v0PdgCode : possiblePID) {

            /* Update TLorentzVectors to get reconstructed mass */

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[v0PdgCode])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

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

            /* Fill histograms */

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_true) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(110);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(140);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(V0->AlphaV0(), V0->PtArmV0());

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[v0PdgCode]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmQt")]->Fill(V0->PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ArmAlpha")]->Fill(V0->AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "ImprvDCAbtwDau")]->Fill(V0->GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", v0PdgCode, "Chi2")]->Fill(V0->GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(V0->AlphaV0(), V0->PtArmV0());
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

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Double_t neg_d_wrt_pv, pos_d_wrt_pv;
    Double_t decay_radius;
    Double_t decay_length;  // decay vertex - primary vertex, could be useful
    Double_t opening_angle;
    Double_t cpa_wrt_pv, dca_wrt_pv;
    Double_t dca_btw_dau, dca_neg_v0, dca_pos_v0;

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

            is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                      getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                      getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;
            is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
            is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;

            if (!PassesV0CutsAs(pdgV0, is_true, is_secondary, is_signal, &vertex, neg_track, pos_track)) continue;

            /* Get V0 properties */

            vertex.GetXYZ(V0_X, V0_Y, V0_Z);
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            lvTrackNeg.SetXYZM(N_Px, N_Py, N_Pz, fPDG.GetParticle(getNegPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvTrackPos.SetXYZM(P_Px, P_Py, P_Pz, fPDG.GetParticle(getPosPdgCode_fromV0PdgCode[pdgV0])->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            decay_radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);

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
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            /* Store V0 */

            esdFoundV0s->push_back(vertex);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[esdFoundV0s->size() - 1] = esdIdxPos;

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_true) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(110);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_secondary) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(140);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(V0_Z);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(vertex.PtArmV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(vertex.AlphaV0());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ImprvDCAbtwDau")]->Fill(vertex.GetDcaV0Daughters());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2")]->Fill(vertex.GetChi2V0());
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(vertex.AlphaV0(), vertex.PtArmV0());
        }  // end of loop over pos. tracks
    }  // end of loop over neg. tracks
}

/*
 Apply cuts to a V0 candidate.
 - Input: `pdgV0`, `v0`, `neg_track`, `pos_track`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, AliESDv0* v0, AliESDtrack* neg_track,
                                                AliESDtrack* pos_track) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(60);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(90);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(120);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(150);

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
    fHist_V0s_Bookkeep[pdgV0]->Fill(61);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(91);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(121);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(151);

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(62);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(92);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(122);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(152);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(63);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(93);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(123);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(153);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(64);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(94);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(124);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(154);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(v0->Xv() - fPrimaryVertex->GetX(), 2) + TMath::Power(v0->Yv() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(v0->Zv() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(65);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(95);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(125);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(155);

    Double_t CPAwrtPV = v0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(66);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(96);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(126);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(156);

    Double_t DCAwrtPV = LinePointDCA(v0->Px(), v0->Py(), v0->Pz(), v0->Xv(), v0->Yv(), v0->Zv(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                     fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(67);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(97);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(127);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(157);

    Double_t DCAbtwDau = TMath::Abs(v0->GetDcaV0Daughters());
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(68);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(98);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(128);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(158);

    Float_t impar_neg_v0[2];
    neg_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_neg_v0);
    Double_t DCAnegV0 = TMath::Sqrt(impar_neg_v0[0] * impar_neg_v0[0] + impar_neg_v0[1] * impar_neg_v0[1]);
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(69);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(99);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(129);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(159);

    Float_t impar_pos_v0[2];
    pos_track->GetDZ(v0->Xv(), v0->Yv(), v0->Zv(), fMagneticField, impar_pos_v0);
    Double_t DCAposV0 = TMath::Sqrt(impar_pos_v0[0] * impar_pos_v0[0] + impar_pos_v0[1] * impar_pos_v0[1]);
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(70);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(100);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(130);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(160);

    Double_t ArmPt = ArmenterosQt(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz);
    Double_t ArmAlpha = ArmenterosAlpha(v0->Px(), v0->Py(), v0->Pz(), N_Px, N_Py, N_Pz, P_Px, P_Py, P_Pz);
    Double_t ArmPtOverAlpha = ArmPt / TMath::Abs(ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(71);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(101);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(131);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(161);

    /* Exclusive for Offline, On-The-Fly, and Custom V0 Finders */

    Double_t Chi2 = v0->GetChi2V0();
    if (kMax_V0_Chi2[pdgV0] && Chi2 > kMax_V0_Chi2[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(72);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(102);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(132);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(162);

    Double_t ImprvDCAbtwDau = v0->GetDcaV0Daughters();
    if (kMax_V0_ImprvDCAbtwDau[pdgV0] && ImprvDCAbtwDau > kMax_V0_ImprvDCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(73);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(103);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(133);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(163);

    return kTRUE;
}

/*
 Apply cuts to a V0 candidate.
 - Input: `pdgV0`, `kfV0`, `kfDaughterNeg`, `kfDaughterPos`, `lvV0`, `lvTrackNeg`, `lvTrackPos`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the V0 candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(Int_t pdgV0, Bool_t IsTrue, Bool_t IsSecondary, Bool_t IsSignal, KFParticleMother kfV0,
                                                KFParticle kfDaughterNeg, KFParticle kfDaughterPos, TLorentzVector lvV0, TLorentzVector lvTrackNeg,
                                                TLorentzVector lvTrackPos) {

    fHist_V0s_Bookkeep[pdgV0]->Fill(60);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(90);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(120);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(150);

    Double_t Radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
    if (kMin_V0_Radius[pdgV0] && Radius < kMin_V0_Radius[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(61);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(91);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(121);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(151);

    Double_t Mass = lvV0.M();
    if (kMin_V0_Mass[pdgV0] && Mass < kMin_V0_Mass[pdgV0]) return kFALSE;
    if (kMax_V0_Mass[pdgV0] && Mass > kMax_V0_Mass[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(62);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(92);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(122);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(152);

    Double_t Pt = lvV0.Pt();
    if (kMin_V0_Pt[pdgV0] && Pt < kMin_V0_Pt[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(63);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(93);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(123);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(153);

    Double_t Eta = TMath::Abs(lvV0.Eta());
    if (kMax_V0_Eta[pdgV0] && Eta > kMax_V0_Eta[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(64);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(94);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(124);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(154);

    Double_t DecayLength = TMath::Sqrt(TMath::Power(kfV0.GetX() - fPrimaryVertex->GetX(), 2) + TMath::Power(kfV0.GetY() - fPrimaryVertex->GetY(), 2) +
                                       TMath::Power(kfV0.GetZ() - fPrimaryVertex->GetZ(), 2));
    if (kMin_V0_DecayLength[pdgV0] && DecayLength < kMin_V0_DecayLength[pdgV0]) return kFALSE;
    if (kMax_V0_DecayLength[pdgV0] && DecayLength > kMax_V0_DecayLength[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(65);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(95);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(125);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(155);

    Double_t CPAwrtPV =
        CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_CPAwrtPV[pdgV0] && CPAwrtPV < kMin_V0_CPAwrtPV[pdgV0]) return kFALSE;
    if (kMax_V0_CPAwrtPV[pdgV0] && CPAwrtPV > kMax_V0_CPAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(66);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(96);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(126);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(156);

    Double_t DCAwrtPV = LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(),
                                     fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_V0_DCAwrtPV[pdgV0] && DCAwrtPV < kMin_V0_DCAwrtPV[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(67);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(97);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(127);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(157);

    Double_t DCAbtwDau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
    if (kMax_V0_DCAbtwDau[pdgV0] && DCAbtwDau > kMax_V0_DCAbtwDau[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(68);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(98);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(128);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(158);

    Double_t DCAnegV0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAnegV0[pdgV0] && DCAnegV0 > kMax_V0_DCAnegV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(69);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(99);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(129);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(159);

    Double_t DCAposV0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));
    if (kMax_V0_DCAposV0[pdgV0] && DCAposV0 > kMax_V0_DCAposV0[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(70);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(100);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(130);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(160);

    Double_t ArmPt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
    Double_t ArmAlpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());
    Double_t ArmPtOverAlpha = TMath::Abs(ArmPt / ArmAlpha);
    if (kMax_V0_ArmPtOverAlpha[pdgV0] && ArmPtOverAlpha > kMax_V0_ArmPtOverAlpha[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(71);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(101);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(131);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(161);

    /* Exclusive for KF method */

    Double_t Chi2ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();
    if (kMax_V0_Chi2ndf[pdgV0] && Chi2ndf > kMax_V0_Chi2ndf[pdgV0]) return kFALSE;
    fHist_V0s_Bookkeep[pdgV0]->Fill(72);
    if (IsTrue) fHist_V0s_Bookkeep[pdgV0]->Fill(102);
    if (IsSecondary) fHist_V0s_Bookkeep[pdgV0]->Fill(132);
    if (IsSignal) fHist_V0s_Bookkeep[pdgV0]->Fill(162);

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

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            if (!PassesSexaquarkCuts_ChannelA(is_signal, SecondaryVertex, lvAntiSexaquark, AntiLambda, KaonZeroShort, AntiLambda_NegDau,
                                              AntiLambda_PosDau, KaonZeroShort_NegDau, KaonZeroShort_PosDau)) {
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

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
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

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
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
    }  // end of loop over anti-lambdas

    return;
}

/*
 Apply anti-sexaquark candidate cuts. Exclusive to reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`.
 - Use: `fPrimaryVertex`, `fMagneticField`
 - Input: `SecondaryVertex`, `lvAntiSexaquark`, `esdAntiLambda`, `esdKaonZeroShort`, `esdAntiLambdaNeg`, `esdAntiLambdaPos`, `esdKaonZeroShortNeg`,
 `esdKaonZeroShortPos`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(Bool_t isSignal, TVector3 SecondaryVertex, TLorentzVector lvAntiSexaquark,
                                                              AliESDv0 esdAntiLambda, AliESDv0 esdKaonZeroShort, AliESDtrack* esdAntiLambdaNeg,
                                                              AliESDtrack* esdAntiLambdaPos, AliESDtrack* esdKaonZeroShortNeg,
                                                              AliESDtrack* esdKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(SecondaryVertex.X() * SecondaryVertex.X() + SecondaryVertex.Y() * SecondaryVertex.Y());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((SecondaryVertex.X() - fPrimaryVertex->GetX()) * (SecondaryVertex.X() - fPrimaryVertex->GetX()) +
                                        (SecondaryVertex.Y() - fPrimaryVertex->GetY()) * (SecondaryVertex.Y() - fPrimaryVertex->GetY()) +
                                        (SecondaryVertex.Z() - fPrimaryVertex->GetZ()) * (SecondaryVertex.Z() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fPrimaryVertex->GetX(),
                                              fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    Double_t dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), SecondaryVertex.X(), SecondaryVertex.Y(),
                                       SecondaryVertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    TVector3 AntiLambda_vertex(esdAntiLambda.Xv(), esdAntiLambda.Yv(), esdAntiLambda.Zv());
    TVector3 AntiLambda_momentum(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz());
    TVector3 KaonZeroShort_vertex(esdKaonZeroShort.Xv(), esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv());
    TVector3 KaonZeroShort_momentum(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz());
    TVector3 pca1, pca2;
    Double_t dca_btw_v0s = Calculate_TwoLinesDCA_v2(AntiLambda_vertex, AntiLambda_momentum, KaonZeroShort_vertex, KaonZeroShort_momentum, pca1, pca2);
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0a_sv = LinePointDCA(esdAntiLambda.Px(), esdAntiLambda.Py(), esdAntiLambda.Pz(), esdAntiLambda.Xv(), esdAntiLambda.Yv(),
                                       esdAntiLambda.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_v0b_sv = LinePointDCA(esdKaonZeroShort.Px(), esdKaonZeroShort.Py(), esdKaonZeroShort.Pz(), esdKaonZeroShort.Xv(),
                                       esdKaonZeroShort.Yv(), esdKaonZeroShort.Zv(), SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z());
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Float_t impar_v0a_neg_sv[2];
    esdAntiLambdaNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_neg_sv);
    Double_t dca_v0a_neg_sv = TMath::Sqrt(impar_v0a_neg_sv[0] * impar_v0a_neg_sv[0] + impar_v0a_neg_sv[1] * impar_v0a_neg_sv[1]);
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Float_t impar_v0a_pos_sv[2];
    esdAntiLambdaPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0a_pos_sv);
    Double_t dca_v0a_pos_sv = TMath::Sqrt(impar_v0a_pos_sv[0] * impar_v0a_pos_sv[0] + impar_v0a_pos_sv[1] * impar_v0a_pos_sv[1]);
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Float_t impar_v0b_neg_sv[2];
    esdKaonZeroShortNeg->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_neg_sv);
    Double_t dca_v0b_neg_sv = TMath::Sqrt(impar_v0b_neg_sv[0] * impar_v0b_neg_sv[0] + impar_v0b_neg_sv[1] * impar_v0b_neg_sv[1]);
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Float_t impar_v0b_pos_sv[2];
    esdKaonZeroShortPos->GetDZ(SecondaryVertex.X(), SecondaryVertex.Y(), SecondaryVertex.Z(), fMagneticField, impar_v0b_pos_sv);
    Double_t dca_v0b_pos_sv = TMath::Sqrt(impar_v0b_pos_sv[0] * impar_v0b_pos_sv[0] + impar_v0b_pos_sv[1] * impar_v0b_pos_sv[1]);
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    return kTRUE;
}

/*
 */
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelD() {
    // PENDING
}

/*
 */
void AliAnalysisTaskSexaquark::GeoSexaquarkFinder_ChannelE() {
    // PENDING
}

/*                          */
/**  V0s -- Kalman Filter  **/
/*** ==================== ***/

/*
 Find all V0s via Kalman Filter.
 - Uses: `fESD`, `fMC`
 - Input: `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 - Output: `kfV0s`, `idxDaughters`
*/
void AliAnalysisTaskSexaquark::KalmanV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {

    AliESDtrack* esdTrackNeg;
    AliESDtrack* esdTrackPos;

    Double_t impar_neg[2], impar_pos[2];

    Double_t decay_radius;
    Double_t decay_length;  // decay radius - primary vertex, could be useful
    Double_t opening_angle;
    Double_t arm_qt, arm_alpha;
    Double_t cpa_wrt_pv, dca_wrt_pv;
    Double_t dca_btw_dau, dca_neg_v0, dca_pos_v0;
    Double_t chi2_ndf;

    Bool_t is_true;
    Bool_t is_secondary;
    Bool_t is_signal;

    Int_t mc_idx_neg, mc_idx_pos;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Choose between anti-lambda, K0S or just pion pairs */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    std::vector<KFParticleMother>* kfFoundV0s;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfNegDau_fromFoundV0Idx;
    std::unordered_map<Int_t, Int_t>* getEsdIdxOfPosDau_fromFoundV0Idx;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;

        kfFoundV0s = &kfAntiLambdas;
        getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromAntiLambdaIdx;
        getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromAntiLambdaIdx;
    } else {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
        if (pdgV0 == 310) {
            kfFoundV0s = &kfKaonsZeroShort;
            getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromKaonZeroShortIdx;
            getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromKaonZeroShortIdx;
        } else if (pdgV0 == 422) {
            kfFoundV0s = &kfPionPairs;
            getEsdIdxOfNegDau_fromFoundV0Idx = &getEsdIdxOfNegDau_fromPionPairIdx;
            getEsdIdxOfPosDau_fromFoundV0Idx = &getEsdIdxOfPosDau_fromPionPairIdx;
        }
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

            is_true = doesEsdIdxHaveMother[esdIdxNeg] && doesEsdIdxHaveMother[esdIdxPos] &&
                      getMotherMcIdx_fromEsdIdx[esdIdxNeg] == getMotherMcIdx_fromEsdIdx[esdIdxPos] &&
                      getPdgCode_fromMcIdx[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] == pdgV0;
            is_secondary = isMcIdxSecondary[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_true;
            is_signal = isMcIdxSignal[getMotherMcIdx_fromEsdIdx[esdIdxNeg]] && is_secondary;
            if (pdgV0 == 422) {
                mc_idx_neg = getMcIdx_fromEsdIdx[esdIdxNeg];
                mc_idx_pos = getMcIdx_fromEsdIdx[esdIdxPos];
                is_signal = isMcIdxSignal[mc_idx_neg] && isMcIdxSignal[mc_idx_pos] && getPdgCode_fromMcIdx[mc_idx_neg] == -211 &&
                            getPdgCode_fromMcIdx[mc_idx_pos] == 211 && !doesMcIdxHaveMother[mc_idx_neg] && !doesMcIdxHaveMother[mc_idx_pos] &&
                            getReactionIdx_fromMcIdx[mc_idx_neg] == getReactionIdx_fromMcIdx[mc_idx_pos];
            }

            if (!PassesV0CutsAs(pdgV0, is_true, is_secondary, is_signal, kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;

            /* Calculate V0 properties */

            decay_radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
            decay_length = TMath::Sqrt((kfV0.GetX() - fPrimaryVertex->GetX()) * (kfV0.GetX() - fPrimaryVertex->GetX()) +
                                       (kfV0.GetY() - fPrimaryVertex->GetY()) * (kfV0.GetY() - fPrimaryVertex->GetY()) +
                                       (kfV0.GetZ() - fPrimaryVertex->GetZ()) * (kfV0.GetZ() - fPrimaryVertex->GetZ()));
            opening_angle = lvTrackNeg.Angle(lvTrackPos.Vect());

            arm_qt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz());
            arm_alpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(), lvTrackNeg.Px(), lvTrackNeg.Py(), lvTrackNeg.Pz(), lvTrackPos.Px(),
                                        lvTrackPos.Py(), lvTrackPos.Pz());

            cpa_wrt_pv = CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                             fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfV0.GetDistanceFromVertex(kfPrimaryVertex));

            dca_btw_dau = TMath::Abs(kfDaughterNeg.GetDistanceFromParticle(kfDaughterPos));
            dca_neg_v0 = TMath::Abs(kfDaughterNeg.GetDistanceFromVertex(kfV0));
            dca_pos_v0 = TMath::Abs(kfDaughterPos.GetDistanceFromVertex(kfV0));

            chi2_ndf = (Double_t)kfV0.GetChi2() / (Double_t)kfV0.GetNDF();

            /* Store V0 */

            kfFoundV0s->push_back(kfV0);

            (*getEsdIdxOfNegDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxNeg;
            (*getEsdIdxOfPosDau_fromFoundV0Idx)[kfFoundV0s->size() - 1] = esdIdxPos;

            /* Fill histograms */

            fHist_V0s_Bookkeep[pdgV0]->Fill(80);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "All", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
            f2DHist_V0s_ArmenterosPodolanski["All"]->Fill(arm_alpha, arm_qt);

            if (pdgV0 != 422) {

                if (!is_true) continue;

                fHist_V0s_Bookkeep[pdgV0]->Fill(110);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Mass")]->Fill(lvV0.M());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pt")]->Fill(lvV0.Pt());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Px")]->Fill(lvV0.Px());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Py")]->Fill(lvV0.Py());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Pz")]->Fill(lvV0.Pz());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Phi")]->Fill(lvV0.Phi());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayRadius")]->Fill(decay_radius);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Zv")]->Fill(kfV0.GetZ());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Eta")]->Fill(lvV0.Eta());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmQt")]->Fill(arm_qt);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
                fHist_V0s[std::make_tuple("Found", "True", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
                f2DHist_V0s_ArmenterosPodolanski["True"]->Fill(arm_alpha, arm_qt);

                if (!is_secondary) continue;

                fHist_V0s_Bookkeep[pdgV0]->Fill(140);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Mass")]->Fill(lvV0.M());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pt")]->Fill(lvV0.Pt());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Px")]->Fill(lvV0.Px());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Py")]->Fill(lvV0.Py());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Pz")]->Fill(lvV0.Pz());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Phi")]->Fill(lvV0.Phi());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayRadius")]->Fill(decay_radius);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Zv")]->Fill(kfV0.GetZ());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Eta")]->Fill(lvV0.Eta());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DecayLength")]->Fill(decay_length);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "OpeningAngle")]->Fill(opening_angle);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmQt")]->Fill(arm_qt);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
                fHist_V0s[std::make_tuple("Found", "Secondary", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
                f2DHist_V0s_ArmenterosPodolanski["Secondary"]->Fill(arm_alpha, arm_qt);
            }

            if (!is_signal) continue;

            fHist_V0s_Bookkeep[pdgV0]->Fill(170);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Mass")]->Fill(lvV0.M());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pt")]->Fill(lvV0.Pt());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Px")]->Fill(lvV0.Px());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Py")]->Fill(lvV0.Py());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Pz")]->Fill(lvV0.Pz());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Phi")]->Fill(lvV0.Phi());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayRadius")]->Fill(decay_radius);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Zv")]->Fill(kfV0.GetZ());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Eta")]->Fill(lvV0.Eta());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Rapidity")]->Fill(lvV0.Rapidity());
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DecayLength")]->Fill(decay_length);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "OpeningAngle")]->Fill(opening_angle);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmQt")]->Fill(arm_qt);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "ArmAlpha")]->Fill(arm_alpha);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAbtwDau")]->Fill(dca_btw_dau);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAnegV0")]->Fill(dca_neg_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "DCAposV0")]->Fill(dca_pos_v0);
            fHist_V0s[std::make_tuple("Found", "Signal", pdgV0, "Chi2ndf")]->Fill(chi2_ndf);
            f2DHist_V0s_ArmenterosPodolanski["Signal"]->Fill(arm_alpha, arm_qt);
        }  // end of loop over pos. tracks
    }  // end of loop over neg. tracks
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`
 - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`
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

    Float_t chi2_ndf;

    Bool_t is_signal;

    Int_t reaction_idx;
    Int_t centrality_bin;
    Double_t pt_weight, radius_weight;
    Double_t product_weight;

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

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            if (!PassesSexaquarkCuts_ChannelA(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfKaonZeroShort, kfAntiLambdaNegDau,
                                              kfAntiLambdaPosDau, kfKaonZeroShortNegDau, kfKaonZeroShortPosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
            decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                       (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                       (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                             kfAntiSexaquark.GetZ(),  //
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
            dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));

            dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPosDau.GetDistanceFromVertex(kfAntiSexaquark));

            chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
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
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
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
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

            reaction_idx = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

            if (fReweightRadius) {
                radius_weight = fReweightRadius ? GetRadiusWeight(reaction_idx) : 1.;

                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Radius")]->Fill(radius, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DecayLength")]->Fill(decay_length, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0aSV")]->Fill(dca_v0a_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0bSV")]->Fill(dca_v0b_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, radius_weight);
            }

            if (fReweightPt) {
                centrality_bin = GetCentralityBin(fCentrality);
                pt_weight = GetPtWeight(reaction_idx, centrality_bin);

                // AliInfoF("reaction_idx = %d, fCentrality = %f, centrality_bin = %d, pt_weight = %f", reaction_idx, fCentrality, centrality_bin,
                //  pt_weight); // DEBUG

                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Radius")]->Fill(radius, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DecayLength")]->Fill(decay_length, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0aSV")]->Fill(dca_v0a_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0bSV")]->Fill(dca_v0b_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, pt_weight);
            }

            if (!fReweightRadius || !fReweightPt) continue;

            product_weight = radius_weight * pt_weight;

            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Radius")]->Fill(radius, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DecayLength")]->Fill(decay_length, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAbtwV0s")]->Fill(dca_btw_v0s, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0aSV")]->Fill(dca_v0a_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0bSV")]->Fill(dca_v0b_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0anegSV")]->Fill(dca_v0a_neg_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0aposSV")]->Fill(dca_v0a_pos_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0bnegSV")]->Fill(dca_v0b_neg_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0bposSV")]->Fill(dca_v0b_pos_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, product_weight);
        }  // end of loop over K0S
    }  // end of loop over anti-lambdas
}

/*
 Apply anti-sexaquark candidate cuts. Exclusive to reaction channel A: `AntiSexaquark,N->AntiLambda,K0S`.
 - Use: `fPrimaryVertex`
 - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfKaonZeroShort, KFParticle kfAntiLambdaNeg,
                                                              KFParticle kfAntiLambdaPos, KFParticle kfKaonZeroShortNeg,
                                                              KFParticle kfKaonZeroShortPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_btw_v0s = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfKaonZeroShort));
    if (kMax_Sexa_DCAbtwV0s && dca_btw_v0s > kMax_Sexa_DCAbtwV0s) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0a_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aSV && dca_v0a_sv > kMax_Sexa_DCAv0aSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_v0b_sv = TMath::Abs(kfKaonZeroShort.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bSV && dca_v0b_sv > kMax_Sexa_DCAv0bSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_v0a_neg_sv = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0anegSV && dca_v0a_neg_sv > kMax_Sexa_DCAv0anegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_v0a_pos_sv = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0aposSV && dca_v0a_pos_sv > kMax_Sexa_DCAv0aposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0b_neg_sv = TMath::Abs(kfKaonZeroShortNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bnegSV && dca_v0b_neg_sv > kMax_Sexa_DCAv0bnegSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0b_pos_sv = TMath::Abs(kfKaonZeroShortPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0bposSV && dca_v0b_pos_sv > kMax_Sexa_DCAv0bposSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    return kTRUE;
}

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel D: `AntiSexaquark,P->AntiLambda,K+`
 - Uses: `fStruckNucleonPDG`, `fPDG`, `fPrimaryVertex`, `fESD`
 - Containers: `kfAntiLambdas`, `esdIndicesOfPosKaonTracks`, `getEsdIdxOf[Neg|Pos]Dau_fromAntiLambdaIdx`
*/
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelD() {

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack* esdPosKaon;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;
    TLorentzVector lvPosKaon;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv, dca_wrt_pv;
    Float_t dca_v0_ba;
    Float_t dca_v0_neg_ba, dca_v0_pos_ba;
    Float_t dca_v0_sv, dca_ba_sv;
    Float_t dca_v0_neg_sv, dca_v0_pos_sv;

    Float_t chi2_ndf;

    Bool_t is_signal;

    Int_t reaction_idx;
    Int_t centrality_bin;
    Double_t pt_weight, radius_weight;
    Double_t product_weight;

    Float_t impar_v0a_neg_sv[2], impar_v0a_pos_sv[2];
    Float_t impar_v0b_neg_sv[2], impar_v0b_pos_sv[2];

    /* Loop over all possible pairs of anti-lambdas and positive kaons */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t& esdIdxPosKaon : esdIndicesOfPosKaonTracks) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];

            if (esdIdxNeg_AntiLambda == esdIdxPosKaon || esdIdxPosKaon == esdIdxPos_AntiLambda) {
                continue;
            }

            /* Kalman Filter */

            KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];

            esdPosKaon = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaon));
            KFParticle kfPosKaon = CreateKFParticle(*esdPosKaon, fPDG.GetParticle(fProductsPDG[1])->Mass(), (Int_t)esdPosKaon->Charge());

            KFParticleMother kfAntiSexaquark;
            kfAntiSexaquark.AddDaughter(kfAntiLambda);
            kfAntiSexaquark.AddDaughter(kfPosKaon);
            kfAntiSexaquark.TransportToDecayVertex();

            /* Transport V0s to the secondary vertex */

            kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
            kfPosKaon.SetProductionVertex(kfAntiSexaquark);

            kfAntiLambda.TransportToProductionVertex();
            kfPosKaon.TransportToProductionVertex();

            /* Get and transport AliESDtracks of anti-lambda daughters */

            esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));

            KFParticle kfAntiLambdaNegDau =
                CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
            KFParticle kfAntiLambdaPosDau =
                CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());

            kfAntiLambdaNegDau = TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                                     (Int_t)esdAntiLambdaNegDau->Charge());
            kfAntiLambdaPosDau = TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                                     (Int_t)esdAntiLambdaPosDau->Charge());
            kfPosKaon =
                TransportKFParticle(kfPosKaon, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(), (Int_t)esdPosKaon->Charge());

            lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
            lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                       fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
            lvPosKaon.SetXYZM(kfPosKaon.Px(), kfPosKaon.Py(), kfPosKaon.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;

            lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

            lvAntiSexaquark = lvAntiLambda + lvPosKaon - lvStruckNucleon;

            /* Apply cuts */

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxPosKaon];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPosKaon];

            if (!PassesSexaquarkCuts_ChannelD(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfPosKaon, kfAntiLambdaNegDau,
                                              kfAntiLambdaPosDau)) {
                continue;
            }

            /* Get anti-sexaquark's properties */

            radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
            decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                       (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                       (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                             kfAntiSexaquark.GetZ(),  //
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));

            dca_v0_ba = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
            dca_v0_neg_ba = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfPosKaon));
            dca_v0_pos_ba = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfPosKaon));

            dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
            dca_ba_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
            dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));

            chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();  // exclusive to KF

            /* Fill histograms */

            fHist_AntiSexaquarks_Bookkeep->Fill(50);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0ba")]->Fill(dca_v0_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0negba")]->Fill(dca_v0_neg_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0posba")]->Fill(dca_v0_pos_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0SV")]->Fill(dca_v0_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAbaSV")]->Fill(dca_ba_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0negSV")]->Fill(dca_v0_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0posSV")]->Fill(dca_v0_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

            if (!is_signal) continue;

            fHist_AntiSexaquarks_Bookkeep->Fill(90);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0ba")]->Fill(dca_v0_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0negba")]->Fill(dca_v0_neg_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0posba")]->Fill(dca_v0_pos_ba);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0SV")]->Fill(dca_v0_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAbaSV")]->Fill(dca_ba_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv);
            fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

            reaction_idx = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

            if (fReweightRadius) {
                radius_weight = GetRadiusWeight(reaction_idx);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Radius")]->Fill(radius, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DecayLength")]->Fill(decay_length, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0ba")]->Fill(dca_v0_ba, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0negba")]->Fill(dca_v0_neg_ba, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0posba")]->Fill(dca_v0_pos_ba, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAbaSV")]->Fill(dca_ba_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, radius_weight);
                fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, radius_weight);
            }

            if (fReweightPt) {
                centrality_bin = GetCentralityBin(fCentrality);
                pt_weight = GetPtWeight(reaction_idx, centrality_bin);

                // AliInfoF("reaction_idx = %d, fCentrality = %f, centrality_bin = %d, pt_weight = %f", reaction_idx, fCentrality, centrality_bin,
                //          pt_weight); // DEBUG

                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Radius")]->Fill(radius, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DecayLength")]->Fill(decay_length, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0ba")]->Fill(dca_v0_ba, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0negba")]->Fill(dca_v0_neg_ba, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0posba")]->Fill(dca_v0_pos_ba, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAbaSV")]->Fill(dca_ba_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, pt_weight);
                fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, pt_weight);
            }

            if (!fReweightRadius || !fReweightPt) continue;

            product_weight = radius_weight * pt_weight;

            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Radius")]->Fill(radius, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DecayLength")]->Fill(decay_length, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0ba")]->Fill(dca_v0_ba, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0negba")]->Fill(dca_v0_neg_ba, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0posba")]->Fill(dca_v0_pos_ba, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAbaSV")]->Fill(dca_ba_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, product_weight);
            fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, product_weight);
        }  // end of loop over pos. kaons
    }  // end of loop over anti-lambdas
}

/*
 Apply anti-sexaquark candidate cuts.
 - Use: `fPrimaryVertex`
 - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort|AntiLambdaNeg|AntiLambdaPos|KaonZeroShortNeg|KaonZeroShortPos]`, `lvAntiSexaquark`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelD(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfPosKaon, KFParticle kfAntiLambdaNeg,
                                                              KFParticle kfAntiLambdaPos) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0SV && dca_v0_sv > kMax_Sexa_DCAv0SV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0_ba = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0ba && dca_v0_ba > kMax_Sexa_DCAv0ba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_ba_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAbaSV && dca_ba_sv > kMax_Sexa_DCAbaSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0negSV && dca_v0_neg_sv > kMax_Sexa_DCAv0negSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0posSV && dca_v0_pos_sv > kMax_Sexa_DCAv0posSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0_neg_ba = TMath::Abs(kfAntiLambdaNeg.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0negba && dca_v0_neg_ba > kMax_Sexa_DCAv0negba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0_pos_ba = TMath::Abs(kfAntiLambdaPos.GetDistanceFromVertex(kfPosKaon));
    if (kMax_Sexa_DCAv0posba && dca_v0_pos_ba > kMax_Sexa_DCAv0posba) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    return kTRUE;
}

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel E: `AntiSexaquark,P->AntiLambda,K+,pi-,pi+`
 - PENDING
*/
void AliAnalysisTaskSexaquark::KalmanSexaquarkFinder_ChannelE() {

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxPiPlus, esdIdxPiMinus;
    Int_t esdIdxPosKaon;

    AliESDtrack *esdAntiLambdaNegDau, *esdAntiLambdaPosDau;
    AliESDtrack *esdPiMinus, *esdPiPlus;
    AliESDtrack* esdPosKaon;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvAntiLambdaNegDau, lvAntiLambdaPosDau;

    TLorentzVector lvPiMinus;
    TLorentzVector lvPiPlus;
    TLorentzVector lvPosKaon;

    TLorentzVector lvStruckNucleon;
    TLorentzVector lvAntiSexaquark;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Auxiliary variables */

    Float_t radius;
    Float_t decay_length;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;
    Float_t dca_v0_sv, dca_v0_neg_sv, dca_v0_pos_sv;
    Float_t dca_pk_sv, dca_pm_sv, dca_pp_sv;
    Float_t dca_v0_pk, dca_v0_pm, dca_v0_pp;
    Float_t dca_pk_pm, dca_pk_pp;

    Float_t chi2_ndf;

    Bool_t is_signal;

    Int_t reaction_idx;
    Int_t centrality_bin;
    Double_t pt_weight, radius_weight;
    Double_t product_weight;

    Float_t impar_v0a_neg_sv[2], impar_v0a_pos_sv[2];
    Float_t impar_v0b_neg_sv[2], impar_v0b_pos_sv[2];

    /* Loop over all possible pairs of anti-lambdas and K0S */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxPionPair = 0; idxPionPair < (Int_t)kfPionPairs.size(); idxPionPair++) {
            for (Int_t& esdIdxPosKaon : esdIndicesOfPosKaonTracks) {

                esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
                esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
                esdIdxPiMinus = getEsdIdxOfNegDau_fromPionPairIdx[idxPionPair];
                esdIdxPiPlus = getEsdIdxOfPosDau_fromPionPairIdx[idxPionPair];

                /* Check that no daughter is repeated */

                std::set<Int_t> uniqueIndices = {esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda, esdIdxPosKaon, esdIdxPiMinus, esdIdxPiPlus};
                if (uniqueIndices.size() != 5) continue;

                /* Kalman Filter */

                KFParticleMother kfAntiLambda = kfAntiLambdas[idxAntiLambda];

                esdPosKaon = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaon));
                KFParticle kfPosKaon = CreateKFParticle(*esdPosKaon, fPDG.GetParticle(fProductsPDG[1])->Mass(), (Int_t)esdPosKaon->Charge());

                esdPiMinus = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPiMinus));
                KFParticle kfPiMinus = CreateKFParticle(*esdPiMinus, fPDG.GetParticle(fProductsPDG[2])->Mass(), (Int_t)esdPiMinus->Charge());

                esdPiPlus = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPiPlus));
                KFParticle kfPiPlus = CreateKFParticle(*esdPiPlus, fPDG.GetParticle(fProductsPDG[3])->Mass(), (Int_t)esdPiPlus->Charge());

                KFParticleMother kfAntiSexaquark;
                kfAntiSexaquark.AddDaughter(kfAntiLambda);
                kfAntiSexaquark.AddDaughter(kfPosKaon);
                kfAntiSexaquark.AddDaughter(kfPiMinus);
                kfAntiSexaquark.AddDaughter(kfPiPlus);
                kfAntiSexaquark.TransportToDecayVertex();

                /* Transport V0s to the secondary vertex */

                kfAntiLambda.SetProductionVertex(kfAntiSexaquark);
                kfPosKaon.SetProductionVertex(kfAntiSexaquark);
                kfPiMinus.SetProductionVertex(kfAntiSexaquark);
                kfPiPlus.SetProductionVertex(kfAntiSexaquark);

                kfAntiLambda.TransportToProductionVertex();
                kfPosKaon.TransportToProductionVertex();
                kfPiMinus.TransportToProductionVertex();
                kfPiPlus.TransportToProductionVertex();

                /* Transport KFParticles of final-state particles */

                esdAntiLambdaNegDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
                esdAntiLambdaPosDau = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));

                KFParticle kfAntiLambdaNegDau =
                    CreateKFParticle(*esdAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass(), (Int_t)esdAntiLambdaNegDau->Charge());
                KFParticle kfAntiLambdaPosDau =
                    CreateKFParticle(*esdAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass(), (Int_t)esdAntiLambdaPosDau->Charge());

                kfAntiLambdaNegDau =
                    TransportKFParticle(kfAntiLambdaNegDau, kfAntiLambdaPosDau, fPDG.GetParticle(fFinalStateProductsPDG[0])->PdgCode(),
                                        (Int_t)esdAntiLambdaNegDau->Charge());
                kfAntiLambdaPosDau =
                    TransportKFParticle(kfAntiLambdaPosDau, kfAntiLambdaNegDau, fPDG.GetParticle(fFinalStateProductsPDG[1])->PdgCode(),
                                        (Int_t)esdAntiLambdaPosDau->Charge());
                kfPosKaon =
                    TransportKFParticle(kfPosKaon, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[2])->PdgCode(), (Int_t)esdPosKaon->Charge());
                kfPiMinus =
                    TransportKFParticle(kfPiMinus, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[3])->PdgCode(), (Int_t)esdPiMinus->Charge());
                kfPiPlus =
                    TransportKFParticle(kfPiPlus, kfAntiLambda, fPDG.GetParticle(fFinalStateProductsPDG[4])->PdgCode(), (Int_t)esdPiPlus->Charge());

                lvAntiLambdaNegDau.SetXYZM(kfAntiLambdaNegDau.Px(), kfAntiLambdaNegDau.Py(), kfAntiLambdaNegDau.Pz(),
                                           fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
                lvAntiLambdaPosDau.SetXYZM(kfAntiLambdaPosDau.Px(), kfAntiLambdaPosDau.Py(), kfAntiLambdaPosDau.Pz(),
                                           fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
                lvPosKaon.SetXYZM(kfPosKaon.Px(), kfPosKaon.Py(), kfPosKaon.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());
                lvPiMinus.SetXYZM(kfPiMinus.Px(), kfPiMinus.Py(), kfPiMinus.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass());
                lvPiPlus.SetXYZM(kfPiPlus.Px(), kfPiPlus.Py(), kfPiPlus.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[4])->Mass());

                /* Reconstruct anti-sexaquark candidate */

                lvAntiLambda = lvAntiLambdaNegDau + lvAntiLambdaPosDau;

                lvStruckNucleon.SetPxPyPzE(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

                lvAntiSexaquark = lvAntiLambda + lvPosKaon + lvPiMinus + lvPiPlus - lvStruckNucleon;

                /* Apply cuts */

                std::set<UInt_t> uniqueReactionID = {getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda], getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda],
                                                     getReactionIdx_fromEsdIdx[esdIdxPosKaon], getReactionIdx_fromEsdIdx[esdIdxPiMinus],
                                                     getReactionIdx_fromEsdIdx[esdIdxPiPlus]};
                is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda] && isEsdIdxSignal[esdIdxPosKaon] &&
                            isEsdIdxSignal[esdIdxPiMinus] && isEsdIdxSignal[esdIdxPiPlus] && (Int_t)uniqueReactionID.size() == 1;

                if (!PassesSexaquarkCuts_ChannelE(is_signal, kfAntiSexaquark, lvAntiSexaquark, kfAntiLambda, kfAntiLambdaNegDau, kfAntiLambdaPosDau,
                                                  kfPosKaon, kfPiMinus, kfPiPlus)) {
                    continue;
                }

                /* Get anti-sexaquark's properties */

                radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
                decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                           (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                           (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
                cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                                 kfAntiSexaquark.GetZ(),  //
                                                 fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
                dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));

                dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
                dca_v0_neg_sv = TMath::Abs(kfAntiLambdaNegDau.GetDistanceFromVertex(kfAntiSexaquark));
                dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pk_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pm_sv = TMath::Abs(kfPiMinus.GetDistanceFromVertex(kfAntiSexaquark));
                dca_pp_sv = TMath::Abs(kfPiPlus.GetDistanceFromVertex(kfAntiSexaquark));

                dca_v0_pk = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
                dca_v0_pm = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiMinus));
                dca_v0_pp = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiPlus));

                dca_pk_pm = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiMinus));
                dca_pk_pp = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiPlus));

                chi2_ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();  // exclusive to KF method

                /* Fill histograms */

                fHist_AntiSexaquarks_Bookkeep->Fill(50);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Mass")]->Fill(lvAntiSexaquark.M());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pt")]->Fill(lvAntiSexaquark.Pt());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Px")]->Fill(lvAntiSexaquark.Px());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Py")]->Fill(lvAntiSexaquark.Py());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Pz")]->Fill(lvAntiSexaquark.Pz());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Phi")]->Fill(lvAntiSexaquark.Phi());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Radius")]->Fill(radius);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Zv")]->Fill(kfAntiSexaquark.GetZ());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Eta")]->Fill(lvAntiSexaquark.Eta());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0SV")]->Fill(dca_v0_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0negSV")]->Fill(dca_v0_neg_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0posSV")]->Fill(dca_v0_pos_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApkSV")]->Fill(dca_pk_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApmSV")]->Fill(dca_pm_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAppSV")]->Fill(dca_pp_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0pk")]->Fill(dca_v0_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0pm")]->Fill(dca_v0_pm);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCAv0pp")]->Fill(dca_v0_pp);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApkpm")]->Fill(dca_pk_pm);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "DCApkpp")]->Fill(dca_pk_pp);
                fHist_AntiSexaquarks[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2_ndf);

                if (!is_signal) continue;

                fHist_AntiSexaquarks_Bookkeep->Fill(90);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvAntiSexaquark.M());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Px")]->Fill(lvAntiSexaquark.Px());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Py")]->Fill(lvAntiSexaquark.Py());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Radius")]->Fill(radius);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity());
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0SV")]->Fill(dca_v0_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApkSV")]->Fill(dca_pk_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApmSV")]->Fill(dca_pm_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAppSV")]->Fill(dca_pp_sv);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0pk")]->Fill(dca_v0_pk);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0pm")]->Fill(dca_v0_pm);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCAv0pp")]->Fill(dca_v0_pp);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApkpm")]->Fill(dca_pk_pm);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "DCApkpp")]->Fill(dca_pk_pp);
                fHist_AntiSexaquarks[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2_ndf);

                reaction_idx = getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda];

                if (fReweightRadius) {
                    radius_weight = GetRadiusWeight(reaction_idx);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Radius")]->Fill(radius, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DecayLength")]->Fill(decay_length, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCApkSV")]->Fill(dca_pk_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCApmSV")]->Fill(dca_pm_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAppSV")]->Fill(dca_pp_sv, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0pk")]->Fill(dca_v0_pk, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0pm")]->Fill(dca_v0_pm, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCAv0pp")]->Fill(dca_v0_pp, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCApkpm")]->Fill(dca_pk_pm, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "DCApkpp")]->Fill(dca_pk_pp, radius_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wRadiusWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, radius_weight);
                }

                if (fReweightPt) {
                    centrality_bin = GetCentralityBin(fCentrality);
                    pt_weight = GetPtWeight(reaction_idx, centrality_bin);

                    // AliInfoF("reaction_idx = %d, fCentrality = %f, centrality_bin = %d, pt_weight = %f", reaction_idx, fCentrality, centrality_bin,
                    //          pt_weight); // DEBUG

                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Radius")]->Fill(radius, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DecayLength")]->Fill(decay_length, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCApkSV")]->Fill(dca_pk_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCApmSV")]->Fill(dca_pm_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAppSV")]->Fill(dca_pp_sv, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0pk")]->Fill(dca_v0_pk, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0pm")]->Fill(dca_v0_pm, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCAv0pp")]->Fill(dca_v0_pp, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCApkpm")]->Fill(dca_pk_pm, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "DCApkpp")]->Fill(dca_pk_pp, pt_weight);
                    fHist_AntiSexaquarks[std::make_tuple("wPtWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, pt_weight);
                }

                if (!fReweightRadius || !fReweightPt) continue;

                product_weight = radius_weight * pt_weight;

                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Mass")]->Fill(lvAntiSexaquark.M(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pt")]->Fill(lvAntiSexaquark.Pt(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Px")]->Fill(lvAntiSexaquark.Px(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Py")]->Fill(lvAntiSexaquark.Py(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Pz")]->Fill(lvAntiSexaquark.Pz(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Phi")]->Fill(lvAntiSexaquark.Phi(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Radius")]->Fill(radius, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Zv")]->Fill(kfAntiSexaquark.GetZ(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Eta")]->Fill(lvAntiSexaquark.Eta(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Rapidity")]->Fill(lvAntiSexaquark.Rapidity(), product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DecayLength")]->Fill(decay_length, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "CPAwrtPV")]->Fill(cpa_wrt_pv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0SV")]->Fill(dca_v0_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0negSV")]->Fill(dca_v0_neg_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0posSV")]->Fill(dca_v0_pos_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCApkSV")]->Fill(dca_pk_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCApmSV")]->Fill(dca_pm_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAppSV")]->Fill(dca_pp_sv, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0pk")]->Fill(dca_v0_pk, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0pm")]->Fill(dca_v0_pm, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCAv0pp")]->Fill(dca_v0_pp, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCApkpm")]->Fill(dca_pk_pm, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "DCApkpp")]->Fill(dca_pk_pp, product_weight);
                fHist_AntiSexaquarks[std::make_tuple("wBothWeights", "Signal", "Chi2ndf")]->Fill(chi2_ndf, product_weight);
            }  // end of loop over K+
        }  // end of loop over pi+pi- pairs
    }  // end of loop over anti-lambdas
}

/*
 Apply anti-sexaquark candidate cuts.
 // PENDING
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelE(Bool_t isSignal, KFParticleMother kfAntiSexaquark, TLorentzVector lvAntiSexaquark,
                                                              KFParticle kfAntiLambda, KFParticle kfAntiLambdaNegDau, KFParticle kfAntiLambdaPosDau,
                                                              KFParticle kfPosKaon, KFParticle kfPiMinus, KFParticle kfPiPlus) {

    fHist_AntiSexaquarks_Bookkeep->Fill(20);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(60);

    Double_t mass = lvAntiSexaquark.M();
    if (kMin_Sexa_Mass && mass < kMin_Sexa_Mass) return kFALSE;
    if (kMax_Sexa_Mass && mass > kMax_Sexa_Mass) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(21);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(61);

    Double_t pt = lvAntiSexaquark.Pt();
    if (kMin_Sexa_Pt && pt < kMin_Sexa_Pt) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(22);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(62);

    Double_t eta = TMath::Abs(lvAntiSexaquark.Eta());
    if (kMax_Sexa_Eta && eta > kMax_Sexa_Eta) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(23);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(63);

    Double_t rapidity = TMath::Abs(lvAntiSexaquark.Rapidity());
    if (kMax_Sexa_Rapidity && rapidity > kMax_Sexa_Rapidity) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(24);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(64);

    Double_t radius = TMath::Sqrt(kfAntiSexaquark.GetX() * kfAntiSexaquark.GetX() + kfAntiSexaquark.GetY() * kfAntiSexaquark.GetY());
    if (kMin_Sexa_Radius && radius < kMin_Sexa_Radius) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(25);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(65);

    Double_t decay_length = TMath::Sqrt((kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) * (kfAntiSexaquark.GetX() - fPrimaryVertex->GetX()) +
                                        (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) * (kfAntiSexaquark.GetY() - fPrimaryVertex->GetY()) +
                                        (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()) * (kfAntiSexaquark.GetZ() - fPrimaryVertex->GetZ()));
    if (kMin_Sexa_DecayLength && decay_length < kMin_Sexa_DecayLength) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(26);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(66);

    Double_t cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),
                                              fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_Sexa_CPAwrtPV && cpa_wrt_pv < kMin_Sexa_CPAwrtPV) return kFALSE;
    if (kMax_Sexa_CPAwrtPV && cpa_wrt_pv > kMax_Sexa_CPAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(27);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(67);

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);
    Double_t dca_wrt_pv = TMath::Abs(kfAntiSexaquark.GetDistanceFromVertex(kfPrimaryVertex));
    if (kMin_Sexa_DCAwrtPV && dca_wrt_pv < kMin_Sexa_DCAwrtPV) return kFALSE;
    if (kMax_Sexa_DCAwrtPV && dca_wrt_pv > kMax_Sexa_DCAwrtPV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(28);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(68);

    Double_t dca_v0_sv = TMath::Abs(kfAntiLambda.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0SV && dca_v0_sv > kMax_Sexa_DCAv0SV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(29);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(69);

    Double_t dca_v0_pos_sv = TMath::Abs(kfAntiLambdaPosDau.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAv0posSV && dca_v0_pos_sv > kMax_Sexa_DCAv0posSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(30);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(70);

    Double_t dca_pk_sv = TMath::Abs(kfPosKaon.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCApkSV && dca_pk_sv > kMax_Sexa_DCApkSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(31);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(71);

    Double_t dca_pm_sv = TMath::Abs(kfPiMinus.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCApmSV && dca_pm_sv > kMax_Sexa_DCApmSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(32);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(72);

    Double_t dca_pp_sv = TMath::Abs(kfPiPlus.GetDistanceFromVertex(kfAntiSexaquark));
    if (kMax_Sexa_DCAppSV && dca_pp_sv > kMax_Sexa_DCAppSV) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(33);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(73);

    Double_t dca_v0_pk = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPosKaon));
    if (kMax_Sexa_DCAv0pk && dca_v0_pk > kMax_Sexa_DCAv0pk) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(34);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(74);

    Double_t dca_v0_pm = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiMinus));
    if (kMax_Sexa_DCAv0pm && dca_v0_pm > kMax_Sexa_DCAv0pm) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(35);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(75);

    Double_t dca_v0_pp = TMath::Abs(kfAntiLambda.GetDistanceFromParticle(kfPiPlus));
    if (kMax_Sexa_DCAv0pp && dca_v0_pp > kMax_Sexa_DCAv0pp) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(36);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(76);

    Double_t dca_pk_pm = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiMinus));
    if (kMax_Sexa_DCApkpm && dca_pk_pm > kMax_Sexa_DCApkpm) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(37);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(77);

    Double_t dca_pk_pp = TMath::Abs(kfPosKaon.GetDistanceFromParticle(kfPiPlus));
    if (kMax_Sexa_DCApkpp && dca_pk_pp > kMax_Sexa_DCApkpp) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(38);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(78);

    Double_t chi2ndf = (Double_t)kfAntiSexaquark.GetChi2() / (Double_t)kfAntiSexaquark.GetNDF();
    if (kMax_Sexa_Chi2ndf && chi2ndf > kMax_Sexa_Chi2ndf) return kFALSE;
    fHist_AntiSexaquarks_Bookkeep->Fill(39);
    if (isSignal) fHist_AntiSexaquarks_Bookkeep->Fill(79);

    return kTRUE;
}

/*                                */
/**  Channel H -- Kalman Filter  **/
/*** ========================== ***/

/*
 Find all K+K+ pairs via Kalman Filter.
 - Uses: `fESD`, `esdIndicesOfPosKaonTracks`, `fHist_PosKaonPairs`
*/
void AliAnalysisTaskSexaquark::KalmanPosKaonPairFinder() {

    const Int_t kKaonPDGCode = 321;
    const Float_t kKaonMass = fPDG.GetParticle(kKaonPDGCode)->Mass();

    AliESDtrack* esdPosKaonA;
    AliESDtrack* esdPosKaonB;

    Int_t esdIdxPosKaonA, esdIdxPosKaonB;

    Double_t impar_a[2], impar_b[2];

    Double_t origin_radius;
    Double_t decay_length;  // origin radius - primary vertex, i.e., the decay length of the mother particle
    Double_t dca_wrt_pv;
    Double_t dca_btw_kaons;
    Double_t dca_pka_sv, dca_pkb_sv;
    Double_t opening_angle;
    Double_t chi2ndf;

    /* Define primary vertex as a KFVertex */

    KFVertex kfPrimaryVertex = CreateKFVertex(*fPrimaryVertex);

    /* Declare TLorentzVectors */

    TLorentzVector lvPosKaonA, lvPosKaonB;
    TLorentzVector lvPosKaonPair;

    /* Loop over all possible pairs of tracks */

    for (Int_t i = 0; i < (Int_t)esdIndicesOfPosKaonTracks.size() - 1; i++) {
        for (Int_t j = i + 1; j < (Int_t)esdIndicesOfPosKaonTracks.size(); j++) {

            esdIdxPosKaonA = esdIndicesOfPosKaonTracks[i];
            esdIdxPosKaonB = esdIndicesOfPosKaonTracks[j];

            /* Sanity check */
            // (shouldn't be necessary)

            if (esdIdxPosKaonA == esdIdxPosKaonB) continue;

            /* Get tracks */

            esdPosKaonA = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaonA));
            esdPosKaonB = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPosKaonB));

            /* Kalman Filter */

            KFParticle kfPosKaonA = CreateKFParticle(*esdPosKaonA, kKaonMass, (Int_t)esdPosKaonA->Charge());
            KFParticle kfPosKaonB = CreateKFParticle(*esdPosKaonB, kKaonMass, (Int_t)esdPosKaonB->Charge());

            KFParticleMother kfPosKaonPair;
            kfPosKaonPair.AddDaughter(kfPosKaonA);
            kfPosKaonPair.AddDaughter(kfPosKaonB);

            /* Transport V0 and daughters */

            kfPosKaonPair.TransportToDecayVertex();

            KFParticle kfTransportedPosKaonA = TransportKFParticle(kfPosKaonA, kfPosKaonB, kKaonPDGCode, (Int_t)esdPosKaonA->Charge());
            KFParticle kfTransportedPosKaonB = TransportKFParticle(kfPosKaonB, kfPosKaonA, kKaonPDGCode, (Int_t)esdPosKaonB->Charge());

            /* Reconstruct V0 */

            lvPosKaonA.SetXYZM(kfTransportedPosKaonA.Px(), kfTransportedPosKaonA.Py(), kfTransportedPosKaonA.Pz(), kKaonMass);
            lvPosKaonB.SetXYZM(kfTransportedPosKaonB.Px(), kfTransportedPosKaonB.Py(), kfTransportedPosKaonB.Pz(), kKaonMass);
            lvPosKaonPair = lvPosKaonA + lvPosKaonB;

            /* Apply cuts */

            if (!PassesPosKaonPairCuts(kfPosKaonPair, kfPosKaonA, kfPosKaonB, lvPosKaonPair, lvPosKaonA, lvPosKaonB)) {
                continue;
            }

            /* Calculate variables */

            origin_radius = TMath::Sqrt(kfPosKaonPair.GetX() * kfPosKaonPair.GetX() + kfPosKaonPair.GetY() * kfPosKaonPair.GetY());

            dca_wrt_pv = TMath::Abs(kfPosKaonPair.GetDistanceFromVertex(kfPrimaryVertex));
            dca_btw_kaons = TMath::Abs(kfPosKaonA.GetDistanceFromParticle(kfPosKaonB));
            dca_pka_sv = TMath::Abs(kfPosKaonA.GetDistanceFromVertex(kfPosKaonPair));
            dca_pkb_sv = TMath::Abs(kfPosKaonB.GetDistanceFromVertex(kfPosKaonPair));
            decay_length = TMath::Sqrt((kfPosKaonPair.GetX() - fPrimaryVertex->GetX()) * (kfPosKaonPair.GetX() - fPrimaryVertex->GetX()) +
                                       (kfPosKaonPair.GetY() - fPrimaryVertex->GetY()) * (kfPosKaonPair.GetY() - fPrimaryVertex->GetY()) +
                                       (kfPosKaonPair.GetZ() - fPrimaryVertex->GetZ()) * (kfPosKaonPair.GetZ() - fPrimaryVertex->GetZ()));
            opening_angle = lvPosKaonA.Angle(lvPosKaonB.Vect());

            chi2ndf = (Double_t)kfPosKaonPair.GetChi2() / (Double_t)kfPosKaonPair.GetNDF();  // exclusive for KF method

            /* Fill histograms */

            fHist_PosKaonPairs_Bookkeep->Fill(10);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Mass")]->Fill(lvPosKaonPair.M());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "OriginRadius")]->Fill(origin_radius);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DecayLength")]->Fill(decay_length);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Xv")]->Fill(kfPosKaonPair.GetX());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Yv")]->Fill(kfPosKaonPair.GetY());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Zv")]->Fill(kfPosKaonPair.GetZ());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Eta")]->Fill(lvPosKaonPair.Eta());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Px")]->Fill(lvPosKaonPair.Px());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Py")]->Fill(lvPosKaonPair.Py());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Pz")]->Fill(lvPosKaonPair.Pz());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Pt")]->Fill(lvPosKaonPair.Pt());
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "OpeningAngle")]->Fill(opening_angle);

            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCAbtwKaons")]->Fill(dca_btw_kaons);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCApkaSV")]->Fill(dca_pka_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "DCApkbSV")]->Fill(dca_pkb_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "All", "Chi2ndf")]->Fill(chi2ndf);

            if (!isEsdIdxSignal[esdIdxPosKaonA] || !isEsdIdxSignal[esdIdxPosKaonB] ||
                getReactionIdx_fromEsdIdx[esdIdxPosKaonA] != getReactionIdx_fromEsdIdx[esdIdxPosKaonB]) {
                continue;
            }

            fHist_PosKaonPairs_Bookkeep->Fill(11);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Mass")]->Fill(lvPosKaonPair.M());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "OriginRadius")]->Fill(origin_radius);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DecayLength")]->Fill(decay_length);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Xv")]->Fill(kfPosKaonPair.GetX());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Yv")]->Fill(kfPosKaonPair.GetY());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Zv")]->Fill(kfPosKaonPair.GetZ());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Eta")]->Fill(lvPosKaonPair.Eta());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Rapidity")]->Fill(lvPosKaonPair.Rapidity());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Px")]->Fill(lvPosKaonPair.Px());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Py")]->Fill(lvPosKaonPair.Py());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Pz")]->Fill(lvPosKaonPair.Pz());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Pt")]->Fill(lvPosKaonPair.Pt());
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "OpeningAngle")]->Fill(opening_angle);

            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAwrtPV")]->Fill(dca_wrt_pv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCAbtwKaons")]->Fill(dca_btw_kaons);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCApkaSV")]->Fill(dca_pka_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "DCApkbSV")]->Fill(dca_pkb_sv);
            fHist_PosKaonPairs[std::make_tuple("Found", "Signal", "Chi2ndf")]->Fill(chi2ndf);

        }  // end of loop over pos. tracks
    }  // end of loop over neg. tracks
}

/*
 Apply cuts to a K+K+ pair
 - Input: `kfPosKaonPair`, `kfPosKaonA`, `kfPosKaonB`, `lvPosKaonPair`, `lvPosKaonA`, `lvPosKaonB`
 - Uses: `fPrimaryVertex`
 - Return: `kTRUE` if the K+K+ pair passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesPosKaonPairCuts(KFParticleMother kfPosKaonPair, KFParticle kfPosKaonA, KFParticle kfPosKaonB,
                                                       TLorentzVector lvPosKaonPair, TLorentzVector lvPosKaonA, TLorentzVector lvPosKaonB) {

    fHist_PosKaonPairs_Bookkeep->Fill(2);

    Double_t Radius = TMath::Sqrt(kfPosKaonPair.GetX() * kfPosKaonPair.GetX() + kfPosKaonPair.GetY() * kfPosKaonPair.GetY());
    if (kMin_PosKaonPair_Radius && Radius < kMin_PosKaonPair_Radius) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(3);

    Double_t Pt = lvPosKaonPair.Pt();
    if (kMin_PosKaonPair_Pt && Pt < kMin_PosKaonPair_Pt) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(4);

    Double_t DCAwrtPV = LinePointDCA(lvPosKaonPair.Px(), lvPosKaonPair.Py(), lvPosKaonPair.Pz(), kfPosKaonPair.GetX(), kfPosKaonPair.GetY(),
                                     kfPosKaonPair.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (kMin_PosKaonPair_DCAwrtPV && DCAwrtPV < kMin_PosKaonPair_DCAwrtPV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(5);

    Double_t DCAbtwKaons = TMath::Abs(kfPosKaonA.GetDistanceFromParticle(kfPosKaonB));
    if (kMax_PosKaonPair_DCAbtwKaons && DCAbtwKaons > kMax_PosKaonPair_DCAbtwKaons) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(6);

    Double_t DCApkaSV = TMath::Abs(kfPosKaonA.GetDistanceFromVertex(kfPosKaonPair));
    if (kMax_PosKaonPair_DCApkaSV && DCApkaSV > kMax_PosKaonPair_DCApkaSV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(7);

    Double_t DCApkbSV = TMath::Abs(kfPosKaonB.GetDistanceFromVertex(kfPosKaonPair));
    if (kMax_PosKaonPair_DCApkbSV && DCApkbSV > kMax_PosKaonPair_DCApkbSV) return kFALSE;
    fHist_PosKaonPairs_Bookkeep->Fill(8);

    Double_t Chi2ndf = (Double_t)kfPosKaonPair.GetChi2() / (Double_t)kfPosKaonPair.GetNDF();
    if (kMax_PosKaonPair_Chi2ndf && Chi2ndf > kMax_PosKaonPair_Chi2ndf) return kFALSE;  // exclusive for KF method
    fHist_PosKaonPairs_Bookkeep->Fill(9);

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
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z, Double_t refPointX,
                                                       Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(lvParticle.Angle(posRelativeToRef));
}

/*
 Overload of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
*/
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z, Double_t refPointX,
                                                       Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*
 This gives the Armenteros-Podolanski alpha.
 (Based on `AliRoot/STEER/ESD/AliESDv0::AlphaV0()`)
*/
Double_t AliAnalysisTaskSexaquark::ArmenterosAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz,
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
Double_t AliAnalysisTaskSexaquark::ArmenterosQt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
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
Double_t AliAnalysisTaskSexaquark::LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz, Double_t V0_X, Double_t V0_Y, Double_t V0_Z,
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

/*
 Get the centrality bin from the centrality value. According to:
 - centrality ranges: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90
 - centrality bin: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
*/
Int_t AliAnalysisTaskSexaquark::GetCentralityBin(Float_t CentralityValue) {
    Int_t centrality_bin = std::min(static_cast<Int_t>(0.1 * CentralityValue + 1), 9);
    if (CentralityValue < 5.) centrality_bin = 0;
    return centrality_bin;
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
    /*
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
     */
    /* MC Rec. Tracks */
    /*
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
     */
    /* V0s */
    /*
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
     */
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

    Int_t n_daughters = 0;
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
    AliInfoF("%6i %6i %6i %6i %6i %6i %6i %6i %6i", evt_mc, idx_mc, mcPart->PdgCode(),
             is_signal,  //
             evt_mother, mother_pid,
             grandmother_pid,  //
             evt_first_dau, evt_last_dau);
#endif
    /*
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
         */
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
    /*
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
         */
}

/*
 Push one of the found V0s into the Event_tt object
*/
void AliAnalysisTaskSexaquark::V0_PushBack(Int_t evt_v0) {

    // V0_tt this_V0;  // = fVec_V0s[evt_v0];

#if DEBUG_MODE == 2
    Int_t evt_pos = this_V0.Idx_Pos;
    Int_t evt_neg = this_V0.Idx_Neg;
    Int_t idx_pos = 0;  // fVec_Idx_MCRec[evt_pos];
    Int_t idx_neg = 0;  // fVec_Idx_MCRec[evt_neg];
    AliInfoF("%6i %6i %6i %6i", evt_pos, evt_neg, idx_pos, idx_neg);
#endif
    /*
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
         */
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
    /*
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
     */
    /* MC Rec. Tracks */
    /*
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
     */
    /* V0s */
    /*
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
 */
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

    float dS[2];
    float dsdr[4][6];
    // kfThis.GetDStoParticle(kfOther, dS, dsdr);
    GetDStoParticleBz(fMagneticField, kfThis, kfOther, dS, dsdr);

    float mP[8], mC[36];
    kfThis.Transport(dS[0], dsdr[0], mP, mC);

    float mM = fPDG.GetParticle(pdgThis)->Mass();
    float mQ = chargeThis;  // only valid for charged particles with Q = +/- 1

    KFParticle kfTransported;
    kfTransported.Create(mP, mC, mQ, mM);

    return kfTransported;
}

/*
TEST.
*/
// clang-format off
void AliAnalysisTaskSexaquark::GetDStoParticleBz( float Bz, const KFParticleBase &p, const KFParticleBase &q, float dS[2], float dsdr[4][6], const float* param1, const float* param2 )  const
{
  /** Calculates dS = l/p parameters for two particles, where \n
   ** 1) l - signed distance to the DCA point with the other particle;\n
   ** 2) p - momentum of the particle; \n
   ** under the assumption of the constant homogeneous field Bz. dS[0] is the transport parameter for the current particle,
   ** dS[1] - for the particle "p".
   ** Also calculates partial derivatives dsdr of the parameters dS[0] and dS[1] over the state vectors of the particles:\n
   ** 1) dsdr[0][6] = d(dS[0])/d(param1);\n
   ** 2) dsdr[1][6] = d(dS[0])/d(param2);\n
   ** 3) dsdr[2][6] = d(dS[1])/d(param1);\n
   ** 4) dsdr[3][6] = d(dS[1])/d(param2);\n
   ** where param1 are parameters of the current particle (if the pointer is not provided it is initialised with fP) and
   ** param2 are parameters of the second particle "q" (if the pointer is not provided it is initialised with q.fP). Parameters
   ** param1 and param2 should be either provided both or both set to null pointers.
   ** \param[in] Bz - magnetic field Bz
   ** \param[in] p - second particle
   ** \param[out] dS[2] - transport parameters dS for the current particle (dS[0]) and the second particle "p" (dS[1])
   ** \param[out] dsdr[4][6] - partial derivatives of the parameters dS[0] and dS[1] over the state vectors of the both particles
   ** \param[in] param1 - optional parameter, is used in case if the parameters of the current particles are rotated
   ** to other coordinate system (see GetDStoParticleBy() function), otherwise p.fP are used
   ** \param[in] param2 - optional parameter, is used in case if the parameters of the second particles are rotated
   ** to other coordinate system (see GetDStoParticleBy() function), otherwise q.fP are used
   **/

  float p_fP[8];           ///< Particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}.
  float q_fP[8];           ///< Particle parameters { X, Y, Z, Px, Py, Pz, E, S[=DecayLength/P]}.
    for (Int_t i = 0; i < 8; i++) {
        p_fP[i] = p.GetParameter(i);
        q_fP[i] = q.GetParameter(i);
    }

  if(!param1)
  {
    param1 = p_fP;
    param2 = q_fP;
  }

  //* Get dS to another particle for Bz field
  const float kOvSqr6 = 1.f/sqrt(float(6.f));
  const float kCLight = 0.000299792458f;

  //in XY plane
  //first root
  const float& bq1 = Bz*p.GetQ()*kCLight;
  const float& bq2 = Bz*q.GetQ()*kCLight;

  const bool& isStraight1 = fabs(bq1) < 1.e-8f;
  const bool& isStraight2 = fabs(bq2) < 1.e-8f;

  if( isStraight1 && isStraight2 )
  {
    p.GetDStoParticleLine(q, dS, dsdr);
    return;
  }

  const float& px1 = param1[3];
  const float& py1 = param1[4];
  const float& pz1 = param1[5];

  const float& px2 = param2[3];
  const float& py2 = param2[4];
  const float& pz2 = param2[5];

  const float& pt12 = px1*px1 + py1*py1;
  const float& pt22 = px2*px2 + py2*py2;

  const float& x01 = param1[0];
  const float& y01 = param1[1];
  const float& z01 = param1[2];

  const float& x02 = param2[0];
  const float& y02 = param2[1];
  const float& z02 = param2[2];

  float dS1[2] = {0.f}, dS2[2]={0.f};

  const float& dx0 = (x01 - x02);
  const float& dy0 = (y01 - y02);
  const float& dr02 = dx0*dx0 + dy0*dy0;
  const float& drp1  = dx0*px1 + dy0*py1;
  const float& dxyp1 = dx0*py1 - dy0*px1;
  const float& drp2  = dx0*px2 + dy0*py2;
  const float& dxyp2 = dx0*py2 - dy0*px2;
  const float& p1p2 = px1*px2 + py1*py2;
  const float& dp1p2 = px1*py2 - px2*py1;

  const float& k11 = (bq2*drp1 - dp1p2);
  const float& k21 = (bq1*(bq2*dxyp1 - p1p2) + bq2*pt12);
  const float& k12 = ((bq1*drp2 - dp1p2));
  const float& k22 = (bq2*(bq1*dxyp2 + p1p2) - bq1*pt22);

  const float& kp = (dxyp1*bq2 - dxyp2*bq1 - p1p2);
  const float& kd = dr02/2.f*bq1*bq2 + kp;
  const float& c1 = -(bq1*kd + pt12*bq2);
  const float& c2 = bq2*kd + pt22*bq1;

  float d1 = pt12*pt22 - kd*kd;
  if(d1<0.f)
    d1 = float(0.f);
  d1 = sqrt( d1 );
  float d2 = pt12*pt22 - kd*kd;
  if(d2<0.f)
    d2 = float(0.f);
  d2 = sqrt( d2 );

  // find two points of closest approach in XY plane

  float dS1dR1[2][6];
  float dS2dR2[2][6];

  float dS1dR2[2][6];
  float dS2dR1[2][6];

  float dk11dr1[6] = {bq2*px1, bq2*py1, 0, bq2*dx0 - py2, bq2*dy0 + px2, 0};
  float dk11dr2[6] = {-bq2*px1, -bq2*py1, 0, py1, -px1, 0};
  float dk12dr1[6] = {bq1*px2, bq1*py2, 0, -py2, px2, 0};
  float dk12dr2[6] = {-bq1*px2, -bq1*py2, 0, bq1*dx0 + py1, bq1*dy0 - px1, 0};
  float dk21dr1[6] = {bq1*bq2*py1, -bq1*bq2*px1, 0, 2*bq2*px1 + bq1*(-(bq2*dy0) - px2), 2*bq2*py1 + bq1*(bq2*dx0 - py2), 0};
  float dk21dr2[6] = {-(bq1*bq2*py1), bq1*bq2*px1, 0, -(bq1*px1), -(bq1*py1), 0};
  float dk22dr1[6] = {bq1*bq2*py2, -(bq1*bq2*px2), 0, bq2*px2, bq2*py2, 0};
  float dk22dr2[6] = {-(bq1*bq2*py2), bq1*bq2*px2, 0, bq2*(-(bq1*dy0) + px1) - 2*bq1*px2, bq2*(bq1*dx0 + py1) - 2*bq1*py2, 0};

  float dkddr1[6] = {bq1*bq2*dx0 + bq2*py1 - bq1*py2, bq1*bq2*dy0 - bq2*px1 + bq1*px2, 0, -bq2*dy0 - px2, bq2*dx0 - py2, 0};
  float dkddr2[6] = {-bq1*bq2*dx0 - bq2*py1 + bq1*py2, -bq1*bq2*dy0 + bq2*px1 - bq1*px2, 0, bq1*dy0 - px1, -bq1*dx0 - py1, 0};

  float dc1dr1[6] = {-(bq1*(bq1*bq2*dx0 + bq2*py1 - bq1*py2)), -(bq1*(bq1*bq2*dy0 - bq2*px1 + bq1*px2)), 0, -2*bq2*px1 - bq1*(-(bq2*dy0) - px2), -2*bq2*py1 - bq1*(bq2*dx0 - py2), 0};
  float dc1dr2[6] = {-(bq1*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2)), -(bq1*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2)), 0, -(bq1*(bq1*dy0 - px1)), -(bq1*(-(bq1*dx0) - py1)), 0};

  float dc2dr1[6] = {bq2*(bq1*bq2*dx0 + bq2*py1 - bq1*py2), bq2*(bq1*bq2*dy0 - bq2*px1 + bq1*px2), 0, bq2*(-(bq2*dy0) - px2), bq2*(bq2*dx0 - py2), 0};
  float dc2dr2[6] = {bq2*(-(bq1*bq2*dx0) - bq2*py1 + bq1*py2), bq2*(-(bq1*bq2*dy0) + bq2*px1 - bq1*px2), 0, bq2*(bq1*dy0 - px1) + 2*bq1*px2, bq2*(-(bq1*dx0) - py1) + 2*bq1*py2, 0};

  float dd1dr1[6] = {0,0,0,0,0,0};
  float dd1dr2[6] = {0,0,0,0,0,0};
  if(d1>0)
  {
    for(int i=0; i<6; i++)
    {
      dd1dr1[i] = -kd/d1*dkddr1[i];
      dd1dr2[i] = -kd/d1*dkddr2[i];
    }
    dd1dr1[3] += px1/d1*pt22; dd1dr1[4] += py1/d1*pt22;
    dd1dr2[3] += px2/d1*pt12; dd1dr2[4] += py2/d1*pt12;
  }

  if(!isStraight1)
  {
    dS1[0] = atan2( bq1*(k11*c1 + k21*d1), (bq1*k11*d1*bq1 - k21*c1) )/bq1;
    dS1[1] = atan2( bq1*(k11*c1 - k21*d1), (-bq1*k11*d1*bq1 - k21*c1) )/bq1;

    float a = bq1*(k11*c1 + k21*d1);
    float b = bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        const float dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        const float dbdr1 = bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[0][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = bq1*(k11*c1 - k21*d1);
    b = -bq1*k11*d1*bq1 - k21*c1;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq1*( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - (dk21dr1[iP]*d1 + k21*dd1dr1[iP]) );
        const float dadr2 = bq1*( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - (dk21dr2[iP]*d1 + k21*dd1dr2[iP]) );
        const float dbdr1 = -bq1*bq1*( dk11dr1[iP]*d1 + k11*dd1dr1[iP] ) - ( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -bq1*bq1*( dk11dr2[iP]*d1 + k11*dd1dr2[iP] ) - ( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS1dR2[1][iP] = 1/bq1 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(!isStraight2)
  {
    dS2[0] = atan2( (bq2*k12*c2 + k22*d2*bq2), (bq2*k12*d2*bq2 - k22*c2) )/bq2;
    dS2[1] = atan2( (bq2*k12*c2 - k22*d2*bq2), (-bq2*k12*d2*bq2 - k22*c2) )/bq2;

    float a = bq2*(k12*c2 + k22*d2);
    float b = bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        const float dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        const float dbdr1 = bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        const float dbdr2 = bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[0][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = bq2*(k12*c2 - k22*d2);
    b = -bq2*k12*d2*bq2 - k22*c2;
    for(int iP=0; iP<6; iP++)
    {
      if(( b*b + a*a ) > 0)
      {
        const float dadr1 = bq2*( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - (dk22dr1[iP]*d1 + k22*dd1dr1[iP]) );
        const float dadr2 = bq2*( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - (dk22dr2[iP]*d1 + k22*dd1dr2[iP]) );
        const float dbdr1 = -bq2*bq2*( dk12dr1[iP]*d1 + k12*dd1dr1[iP] ) - (dk22dr1[iP]*c2 + k22*dc2dr1[iP]);
        const float dbdr2 = -bq2*bq2*( dk12dr2[iP]*d1 + k12*dd1dr2[iP] ) - (dk22dr2[iP]*c2 + k22*dc2dr2[iP]);

        dS2dR1[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr1*b - dbdr1*a );
        dS2dR2[1][iP] = 1/bq2 * 1/( b*b + a*a ) * ( dadr2*b - dbdr2*a );
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight1 && (pt12>0.f) )
  {
    dS1[0] = (k11*c1 + k21*d1)/(- k21*c1);
    dS1[1] = (k11*c1 - k21*d1)/(- k21*c1);

    float a = k11*c1 + k21*d1;
    float b = -k21*c1;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] + dk21dr1[iP]*d1 + k21*dd1dr1[iP] );
        const float dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] + dk21dr2[iP]*d1 + k21*dd1dr2[iP] );
        const float dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[0][iP] = 0;
        dS1dR2[0][iP] = 0;
      }
    }

    a = k11*c1 - k21*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk11dr1[iP]*c1 + k11*dc1dr1[iP] - dk21dr1[iP]*d1 - k21*dd1dr1[iP] );
        const float dadr2 = ( dk11dr2[iP]*c1 + k11*dc1dr2[iP] - dk21dr2[iP]*d1 - k21*dd1dr2[iP] );
        const float dbdr1 = -( dk21dr1[iP]*c1 + k21*dc1dr1[iP] );
        const float dbdr2 = -( dk21dr2[iP]*c1 + k21*dc1dr2[iP] );

        dS1dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS1dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS1dR1[1][iP] = 0;
        dS1dR2[1][iP] = 0;
      }
    }
  }
  if(isStraight2 && (pt22>0.f) )
  {
    dS2[0] = (k12*c2 + k22*d2)/(- k22*c2);
    dS2[1] = (k12*c2 - k22*d2)/(- k22*c2);

    float a = k12*c2 + k22*d1;
    float b = -k22*c2;

    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] + dk22dr1[iP]*d1 + k22*dd1dr1[iP] );
        const float dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] + dk22dr2[iP]*d1 + k22*dd1dr2[iP] );
        const float dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        const float dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[0][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[0][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[0][iP] = 0;
        dS2dR2[0][iP] = 0;
      }
    }

    a = k12*c2 - k22*d1;
    for(int iP=0; iP<6; iP++)
    {
      if(b*b > 0)
      {
        const float dadr1 = ( dk12dr1[iP]*c2 + k12*dc2dr1[iP] - dk22dr1[iP]*d1 - k22*dd1dr1[iP] );
        const float dadr2 = ( dk12dr2[iP]*c2 + k12*dc2dr2[iP] - dk22dr2[iP]*d1 - k22*dd1dr2[iP] );
        const float dbdr1 = -( dk22dr1[iP]*c2 + k22*dc2dr1[iP] );
        const float dbdr2 = -( dk22dr2[iP]*c2 + k22*dc2dr2[iP] );

        dS2dR1[1][iP] = dadr1/b - dbdr1*a/(b*b) ;
        dS2dR2[1][iP] = dadr2/b - dbdr2*a/(b*b) ;
      }
      else
      {
        dS2dR1[1][iP] = 0;
        dS2dR2[1][iP] = 0;
      }
    }
  }

  //select a point which is close to the primary vertex (with the smallest r)

  float dr2[2];
  for(int iP = 0; iP<2; iP++)
  {
    const float& bs1 = bq1*dS1[iP];
    const float& bs2 = bq2*dS2[iP];
    float sss = sin(bs1), ccc = cos(bs1);

    const bool& bs1Big = fabs(bs1) > 1.e-8f;
    const bool& bs2Big = fabs(bs2) > 1.e-8f;

    float sB(0.f), cB(0.f);
    if(bs1Big)
    {
      sB = sss/bq1;
      cB = (1.f-ccc)/bq1;
    }
    else
    {
      sB = ((1.f-bs1*kOvSqr6)*(1.f+bs1*kOvSqr6)*dS1[iP]);
      cB = .5f*sB*bs1;
    }

    const float& x1 = param1[0] + sB*px1 + cB*py1;
    const float& y1 = param1[1] - cB*px1 + sB*py1;
    const float& z1 = param1[2] + dS1[iP]*param1[5];

    sss = sin(bs2), ccc = cos(bs2);

    if(bs2Big)
    {
      sB = sss/bq2;
      cB = (1.f-ccc)/bq2;
    }
    else
    {
      sB = ((1.f-bs2*kOvSqr6)*(1.f+bs2*kOvSqr6)*dS2[iP]);
      cB = .5f*sB*bs2;
    }

    const float& x2 = param2[0] + sB*px2 + cB*py2;
    const float& y2 = param2[1] - cB*px2 + sB*py2;
    const float& z2 = param2[2] + dS2[iP]*param2[5];

    float dx = (x1-x2);
    float dy = (y1-y2);
    float dz = (z1-z2);

    dr2[iP] = dx*dx + dy*dy + dz*dz;
  }

  const bool isFirstRoot = dr2[0] < dr2[1];
  if(isFirstRoot)
  {
    dS[0]  = dS1[0];
    dS[1] = dS2[0];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[0][iP];
      dsdr[1][iP] = dS1dR2[0][iP];
      dsdr[2][iP] = dS2dR1[0][iP];
      dsdr[3][iP] = dS2dR2[0][iP];
    }
  }
  else
  {
    dS[0]  = dS1[1];
    dS[1] = dS2[1];

    for(int iP=0; iP<6; iP++)
    {
      dsdr[0][iP] = dS1dR1[1][iP];
      dsdr[1][iP] = dS1dR2[1][iP];
      dsdr[2][iP] = dS2dR1[1][iP];
      dsdr[3][iP] = dS2dR2[1][iP];
    }
  }

  //find correct parts of helices
//   int n1(0);
//   int n2(0);
//   float dzMin = fabs( (z01-z02) + dS[0]*pz1 - dS[1]*pz2 );
//   const float pi2(6.283185307f);

  //TODO optimise for loops for neutral particles
//   const float& i1Float = -bq1/pi2*(z01/pz1+dS[0]);
//   for(int di1=-1; di1<=1; di1++)
//   {
//     int i1(0);
//     if(!isStraight1)
//       i1 = int(i1Float) + di1;
//
//     const float& i2Float = ( ((z01-z02) + (dS[0]+pi2*i1/bq1)*pz1)/pz2 - dS[1]) * bq2/pi2;
//     for(int di2 = -1; di2<=1; di2++)
//     {
//       int i2(0);
//       if(!isStraight2)
//         i2 = int(i2Float) + di2;
//
//       const float& z1 = z01 + (dS[0]+pi2*i1/bq1)*pz1;
//       const float& z2 = z02 + (dS[1]+pi2*i2/bq2)*pz2;
//       const float& dz = fabs( z1-z2 );
//
//       if(dz < dzMin)
//       {
//         n1 = i1;
//         n2 = i2;
//         dzMin = dz;
//       }
//     }
//   }
//
//   if(!isStraight1)
//     dS[0] += float(n1)*pi2/bq1;
//   if(!isStraight2)
//     dS[1] += float(n2)*pi2/bq2;

  //Line correction
  {
    const float& bs1 = bq1*dS[0];
    const float& bs2 = bq2*dS[1];
    float sss = sin(bs1), ccc = cos(bs1);

    const bool& bs1Big = fabs(bs1) > 1.e-8f;
    const bool& bs2Big = fabs(bs2) > 1.e-8f;

    float sB(0.f), cB(0.f);
    if(bs1Big)
    {
      sB = sss/bq1;
      cB = (1.f-ccc)/bq1;
    }
    else
    {
      sB = ((1.f-bs1*kOvSqr6)*(1.f+bs1*kOvSqr6)*dS[0]);
      cB = .5f*sB*bs1;
    }

    const float& x1 = x01 + sB*px1 + cB*py1;
    const float& y1 = y01 - cB*px1 + sB*py1;
    const float& z1 = z01 + dS[0]*pz1;
    const float& ppx1 =  ccc*px1 + sss*py1;
    const float& ppy1 = -sss*px1 + ccc*py1;
    const float& ppz1 = pz1;

    float sss1 = sin(bs2), ccc1 = cos(bs2);

    float sB1(0.f), cB1(0.f);
    if(bs2Big)
    {
      sB1 = sss1/bq2;
      cB1 = (1.f-ccc1)/bq2;
    }
    else
    {
      sB1 = ((1.f-bs2*kOvSqr6)*(1.f+bs2*kOvSqr6)*dS[1]);
      cB1 = .5f*sB1*bs2;
    }

    const float& x2 = x02 + sB1*px2 + cB1*py2;
    const float& y2 = y02 - cB1*px2 + sB1*py2;
    const float& z2 = z02 + dS[1]*pz2;
    const float& ppx2 =  ccc1*px2 + sss1*py2;
    const float& ppy2 = -sss1*px2 + ccc1*py2;
    const float& ppz2 = pz2;

    const float& p12  = ppx1*ppx1 + ppy1*ppy1 + ppz1*ppz1;
    const float& p22  = ppx2*ppx2 + ppy2*ppy2 + ppz2*ppz2;
    const float& lp1p2 = ppx1*ppx2 + ppy1*ppy2 + ppz1*ppz2;

    const float& dx = (x2 - x1);
    const float& dy = (y2 - y1);
    const float& dz = (z2 - z1);

    const float& ldrp1 = ppx1*dx + ppy1*dy + ppz1*dz;
    const float& ldrp2 = ppx2*dx + ppy2*dy + ppz2*dz;

    float detp =  lp1p2*lp1p2 - p12*p22;
    if( fabs(detp)<1.e-4 ) detp = 1; //TODO correct!!!

    //dsdr calculation
    const float a1 = ldrp2*lp1p2 - ldrp1*p22;
    const float a2 = ldrp2*p12 - ldrp1*lp1p2;
    const float lp1p2_ds0 = bq1*( ppx2*ppy1 - ppy2*ppx1);
    const float lp1p2_ds1 = bq2*( ppx1*ppy2 - ppy1*ppx2);
    const float ldrp1_ds0 = -p12 + bq1*(ppy1*dx - ppx1*dy);
    const float ldrp1_ds1 =  lp1p2;
    const float ldrp2_ds0 = -lp1p2;
    const float ldrp2_ds1 =  p22 + bq2*(ppy2*dx - ppx2*dy);
    const float detp_ds0 = 2*lp1p2*lp1p2_ds0;
    const float detp_ds1 = 2*lp1p2*lp1p2_ds1;
    const float a1_ds0 = ldrp2_ds0*lp1p2 + ldrp2*lp1p2_ds0 - ldrp1_ds0*p22;
    const float a1_ds1 = ldrp2_ds1*lp1p2 + ldrp2*lp1p2_ds1 - ldrp1_ds1*p22;
    const float a2_ds0 = ldrp2_ds0*p12 - ldrp1_ds0*lp1p2 - ldrp1*lp1p2_ds0;
    const float a2_ds1 = ldrp2_ds1*p12 - ldrp1_ds1*lp1p2 - ldrp1*lp1p2_ds1;

    // AliInfoF("a1_ds0 a1 detp_ds0 detp a1_ds1 detp_ds1 a2_ds0 a2 a2_ds1 = %f %f %f %f %f %f %f %f %f", a1_ds0, a1, detp_ds0, detp, a1_ds1, detp_ds1, a2_ds0, a2, a2_ds1);

    const float dsl1ds0 = a1_ds0/detp - a1*detp_ds0/(detp*detp);
    const float dsl1ds1 = a1_ds1/detp - a1*detp_ds1/(detp*detp);
    const float dsl2ds0 = a2_ds0/detp - a2*detp_ds0/(detp*detp);
    const float dsl2ds1 = a2_ds1/detp - a2*detp_ds1/(detp*detp);

    float dsldr[4][6];
    for(int iP=0; iP<6; iP++)
    {
      dsldr[0][iP] = dsl1ds0*dsdr[0][iP] + dsl1ds1*dsdr[2][iP];
      dsldr[1][iP] = dsl1ds0*dsdr[1][iP] + dsl1ds1*dsdr[3][iP];
      dsldr[2][iP] = dsl2ds0*dsdr[0][iP] + dsl2ds1*dsdr[2][iP];
      dsldr[3][iP] = dsl2ds0*dsdr[1][iP] + dsl2ds1*dsdr[3][iP];
    }

    for(int iDS=0; iDS<4; iDS++)
      for(int iP=0; iP<6; iP++)
        dsdr[iDS][iP] += dsldr[iDS][iP];

    const float lp1p2_dr0[6] = {0, 0, 0, ccc*ppx2 - ppy2*sss, ccc*ppy2 + ppx2*sss, pz2};
    const float lp1p2_dr1[6] = {0, 0, 0, ccc1*ppx1 - ppy1*sss1, ccc1*ppy1 + ppx1*sss1, pz1};
    const float ldrp1_dr0[6] = {-ppx1, -ppy1, -pz1,  cB*ppy1 - ppx1*sB + ccc*dx - sss*dy, -cB*ppx1-ppy1*sB + sss*dx + ccc*dy, -dS[0]*pz1 + dz};
    const float ldrp1_dr1[6] = { ppx1,  ppy1,  pz1, -cB1*ppy1 + ppx1*sB1, cB1*ppx1 + ppy1*sB1, dS[1]*pz1};
    const float ldrp2_dr0[6] = {-ppx2, -ppy2, -pz2, cB*ppy2 - ppx2*sB, -cB*ppx2-ppy2*sB, -dS[0]*pz2};
    const float ldrp2_dr1[6] = {ppx2, ppy2, pz2, -cB1*ppy2 + ppx2*sB1 + ccc1*dx- sss1*dy, cB1*ppx2 + ppy2*sB1 + sss1*dx + ccc1*dy, dz + dS[1]*pz2};
    const float p12_dr0[6] = {0, 0, 0, 2*px1, 2*py1, 2*pz1};
    const float p22_dr1[6] = {0, 0, 0, 2*px2, 2*py2, 2*pz2};
    float a1_dr0[6], a1_dr1[6], a2_dr0[6], a2_dr1[6], detp_dr0[6], detp_dr1[6];
    for(int iP=0; iP<6; iP++)
    {
      a1_dr0[iP] = ldrp2_dr0[iP]*lp1p2 + ldrp2*lp1p2_dr0[iP] - ldrp1_dr0[iP]*p22;
      a1_dr1[iP] = ldrp2_dr1[iP]*lp1p2 + ldrp2*lp1p2_dr1[iP] - ldrp1_dr1[iP]*p22 - ldrp1*p22_dr1[iP];
      a2_dr0[iP] = ldrp2_dr0[iP]*p12 + ldrp2*p12_dr0[iP] - ldrp1_dr0[iP]*lp1p2 - ldrp1*lp1p2_dr0[iP];
      a2_dr1[iP] = ldrp2_dr1[iP]*p12 - ldrp1_dr1[iP]*lp1p2 - ldrp1*lp1p2_dr1[iP];
      detp_dr0[iP] = 2*lp1p2*lp1p2_dr0[iP] - p12_dr0[iP]*p22;
      detp_dr1[iP] = 2*lp1p2*lp1p2_dr1[iP] - p12*p22_dr1[iP];

      dsdr[0][iP] += a1_dr0[iP]/detp - a1*detp_dr0[iP]/(detp*detp);
      dsdr[1][iP] += a1_dr1[iP]/detp - a1*detp_dr1[iP]/(detp*detp);
      dsdr[2][iP] += a2_dr0[iP]/detp - a2*detp_dr0[iP]/(detp*detp);
      dsdr[3][iP] += a2_dr1[iP]/detp - a2*detp_dr1[iP]/(detp*detp);
    }

    dS[0] += (ldrp2*lp1p2 - ldrp1*p22) /detp;
    dS[1] += (ldrp2*p12 - ldrp1*lp1p2)/detp;
  }
}
// clang-format on

/*
 Clear all containers.
*/
void AliAnalysisTaskSexaquark::ClearContainers() {

    /*  */
    getPdgCode_fromMcIdx.clear();

    isMcIdxSignal.clear();
    isMcIdxSecondary.clear();

    getReactionIdx_fromMcIdx.clear();
    getMcIdx_fromReactionIdx.clear();
    getPt_fromReactionIdx.clear();
    getRadius_fromReactionIdx.clear();

    doesMcIdxHaveMother.clear();
    getMotherMcIdx_fromMcIdx.clear();
    getNegDauMcIdx_fromMcIdx.clear();
    getPosDauMcIdx_fromMcIdx.clear();

    mcIndicesOfTrueV0s.clear();

    /*  */
    getMcIdx_fromEsdIdx.clear();
    esdIdxPassesTrackSelection_AsProton.clear();
    esdIdxPassesTrackSelection_AsAntiProton.clear();
    esdIdxPassesTrackSelection_AsPosKaon.clear();
    esdIdxPassesTrackSelection_AsNegKaon.clear();
    esdIdxPassesTrackSelection_AsPiPlus.clear();
    esdIdxPassesTrackSelection_AsPiMinus.clear();
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
    esdIndicesOfPiMinusTracks.clear();
    esdIndicesOfPiPlusTracks.clear();

    /*  */
    getEsdIdxOfNegDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfPosDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfNegDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfPosDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfNegDau_fromPionPairIdx.clear();
    getEsdIdxOfPosDau_fromPionPairIdx.clear();

    esdAntiLambdas.clear();
    esdKaonsZeroShort.clear();

    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
    kfPionPairs.clear();
}

/*                    */
/**  External Files  **/
/*** ============== ***/

/*
 Hola hola.
 It's loaded after UserCreateOutputObjects() and before UserExec()
 - Uses: `fAliEnPath`, `fLogTree`
 Note: must be executed ONLY when analyzing SIGNAL SIMs, protection PENDING!
*/
Bool_t AliAnalysisTaskSexaquark::LoadLogsIntoTree() {

    TGrid* alien = nullptr;
    if (!gGrid) {
        alien = TGrid::Connect("alien://");
        if (!alien) return kFALSE;
    }

    TString AliEn_Dir = fAliEnPath(0, fAliEnPath.Last('/'));

    TString orig_path = Form("%s/sim.log", AliEn_Dir.Data());
    AliInfoF("!! Copying file %s ... !!", orig_path.Data());

    // ! assuming path ends with format .../LHC23l1a3/A1.73/297595/001/sim.log ! //
    TObjArray* tokens = fAliEnPath.Tokenize("/");
    Int_t AliEn_DirNumber = ((TObjString*)tokens->At(tokens->GetEntries() - 2))->GetString().Atoi();
    Int_t AliEn_RunNumber = ((TObjString*)tokens->At(tokens->GetEntries() - 3))->GetString().Atoi();
    TString AliEn_SimSubSet = ((TObjString*)tokens->At(tokens->GetEntries() - 4))->GetString();

    TString Log_NewBasename = Form("sim_%s_%i_%03i.log", AliEn_SimSubSet.Data(), AliEn_RunNumber, AliEn_DirNumber);

    if (AliEn_Dir.BeginsWith("alien://")) {
        gSystem->Exec(Form("alien.py cp %s file://./%s", orig_path.Data(), Log_NewBasename.Data()));
    } else {
        gSystem->Exec(Form("cp %s ./%s", orig_path.Data(), Log_NewBasename.Data()));
    }

    TString new_path = Form("%s/%s", gSystem->pwd(), Log_NewBasename.Data());
    AliInfoF("!! Reading file %s ... !!", new_path.Data());

    std::ifstream SimLog(new_path);
    if (!SimLog.is_open()) {
        AliInfo("!! Unable to open file !!");
        return kFALSE;
    }

    Char_t ReactionChannelLetter = AliEn_SimSubSet[0];
    Double_t SM = TString(AliEn_SimSubSet(1, 4)).Atof();

    Int_t EventID = -1;
    Int_t ReactionID;
    Int_t NPDGCode;
    Double_t SPx, SPy, SPz;
    Double_t NPx, NPy, NPz;

    if (fLogTree->GetBranch("ReactionChannel")) {
        fLogTree->SetBranchAddress("RunNumber", &AliEn_RunNumber);
        fLogTree->SetBranchAddress("DirNumber", &AliEn_DirNumber);
        fLogTree->SetBranchAddress("ReactionChannel", &ReactionChannelLetter);
        fLogTree->SetBranchAddress("EventID", &EventID);
        fLogTree->SetBranchAddress("ReactionID", &ReactionID);
        fLogTree->SetBranchAddress("Nucl_PDGCode", &NPDGCode);
        fLogTree->SetBranchAddress("Sexa_Px_ini", &SPx);
        fLogTree->SetBranchAddress("Sexa_Py_ini", &SPy);
        fLogTree->SetBranchAddress("Sexa_Pz_ini", &SPz);
        fLogTree->SetBranchAddress("Sexa_M_ini", &SM);
        fLogTree->SetBranchAddress("Fermi_Px", &NPx);
        fLogTree->SetBranchAddress("Fermi_Py", &NPy);
        fLogTree->SetBranchAddress("Fermi_Pz", &NPz);
    } else {
        fLogTree->Branch("RunNumber", &AliEn_RunNumber);
        fLogTree->Branch("DirNumber", &AliEn_DirNumber);
        fLogTree->Branch("ReactionChannel", &ReactionChannelLetter);
        fLogTree->Branch("EventID", &EventID);
        fLogTree->Branch("ReactionID", &ReactionID);
        fLogTree->Branch("Nucl_PDGCode", &NPDGCode);
        fLogTree->Branch("Sexa_Px_ini", &SPx);
        fLogTree->Branch("Sexa_Py_ini", &SPy);
        fLogTree->Branch("Sexa_Pz_ini", &SPz);
        fLogTree->Branch("Sexa_M_ini", &SM);
        fLogTree->Branch("Fermi_Px", &NPx);
        fLogTree->Branch("Fermi_Py", &NPy);
        fLogTree->Branch("Fermi_Pz", &NPz);
    }

    // auxiliary variables
    std::string cstr_line;
    TString tstr_line, csv;
    TObjArray* csv_arr = nullptr;

    /* Read lines */

    while (std::getline(SimLog, cstr_line)) {

        tstr_line = cstr_line;

        // a new event has appeared
        if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) EventID++;

        if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;

        csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
        csv_arr = csv.Tokenize(",");

        ReactionID = dynamic_cast<TObjString*>(csv_arr->At(0))->String().Atoi();
        NPDGCode = ReactionChannelLetter == 'A' ? 2112 : 2212;  // considering only ADEH
        SPx = dynamic_cast<TObjString*>(csv_arr->At(1))->String().Atof();
        SPy = dynamic_cast<TObjString*>(csv_arr->At(2))->String().Atof();
        SPz = dynamic_cast<TObjString*>(csv_arr->At(3))->String().Atof();
        NPx = dynamic_cast<TObjString*>(csv_arr->At(4))->String().Atof();
        NPy = dynamic_cast<TObjString*>(csv_arr->At(5))->String().Atof();
        NPz = dynamic_cast<TObjString*>(csv_arr->At(6))->String().Atof();

        fLogTree->Fill();
    }  // end of loop over lines

    AliInfo("!! Closing file ... !!");
    SimLog.close();

    return kTRUE;
}

/*
  Add pT weights.
*/
void AliAnalysisTaskSexaquark::AddPtWeights(std::vector<TH1D*> ptWeights) {
    for (TH1D* pt_weight_hist : ptWeights) {
        TH1D* clone_pt_weight_hist = dynamic_cast<TH1D*>(pt_weight_hist->Clone());
        clone_pt_weight_hist->Scale(1. / clone_pt_weight_hist->Integral());  // normalize
        // clone_pt_weight_hist->Print(); // DEBUG
        fPtWeights.push_back(clone_pt_weight_hist);
    }
}

/*
 Add radius weights.
*/
void AliAnalysisTaskSexaquark::AddRadiusWeights(TH1F* radiusWeights) {
    fRadiusWeights = dynamic_cast<TH1F*>(radiusWeights->Clone());
    fRadiusWeights->Scale(1. / fRadiusWeights->Integral());  // normalize
    // fRadiusWeights->Print(); // DEBUG
}

/*
  For a signal candidate, given its reaction index, get its true Pt,
  and find its weight on the `fPtWeights` histogram.
*/
Double_t AliAnalysisTaskSexaquark::GetPtWeight(Int_t reactionIdx, Int_t centralityBin) {
    Float_t pt = getPt_fromReactionIdx[reactionIdx];
    return fPtWeights[centralityBin]->GetBinContent(fPtWeights[centralityBin]->FindBin(pt));
}

/*
  For a signal candidate, given its reaction index, get its true interaction radius,
  and find its weight on the `fRadiusWeights` histogram.
*/
Double_t AliAnalysisTaskSexaquark::GetRadiusWeight(Int_t reactionIdx) {
    Float_t radius = getRadius_fromReactionIdx[reactionIdx];
    return fRadiusWeights->GetBinContent(fRadiusWeights->FindBin(radius));
}
