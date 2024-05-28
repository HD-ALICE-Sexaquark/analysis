#include "AliQuickTaskTPCProperties.h"

ClassImp(AliQuickTaskTPCProperties);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliQuickTaskTPCProperties::AliQuickTaskTPCProperties()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.) {}

/*
 Constructor, called locally.
*/
AliQuickTaskTPCProperties::AliQuickTaskTPCProperties(const char* name, Bool_t IsMC)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      fPDG(),
      kMax_NSigma_Pion(0.),
      kMax_NSigma_Kaon(0.),
      kMax_NSigma_Proton(0.),
      kMax_Track_Eta(0.),
      kMin_Track_NTPCClusters(0.),
      kMax_Track_Chi2PerNTPCClusters(0.),
      kTurnedOn_Track_StatusCuts(0.),
      kTurnedOn_Track_RejectKinks(0.),
      kMin_Track_DCAxy_wrtPV(0.),
      kMin_Track_DCAz_wrtPV(0.) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/*
 Destructor.
*/
AliQuickTaskTPCProperties::~AliQuickTaskTPCProperties() {
    if (fOutputListOfHists) delete fOutputListOfHists;
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliQuickTaskTPCProperties::UserCreateOutputObjects() {

    AliInfo("Initializing AliQuickTaskTPCProperties...");
    AliInfo("INPUT OPTIONS:");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);

    DefineTracksCuts("Standard");

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

    /* Histograms */

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    PrepareTracksHistograms();

    PostData(1, fOutputListOfHists);
}

/*
 Define and add tracks' histograms to the histograms output list.
 - Uses: `fHist_Tracks_Bookkeep`, `fHist_Tracks`, `fOutputListOfHists`
*/
void AliQuickTaskTPCProperties::PrepareTracksHistograms() {

    std::vector<TString> tracks_sets = {"Selected", "True", "Secondary"};

    std::vector<Int_t> species_codes = {2212, -2212, 321, 211, -211};
    std::vector<TString> species_names = {"Proton", "AntiProton", "PosKaon", "PiPlus", "PiMinus"};

    // clang-format off
    std::vector<TString> tracks_props = {"Px",          "Py",                "Pt",              "Pz",           "Eta",
                                         "DCAxy_wrtPV", "DCAz_wrtPV",        "NTPCClusters",    "NCrossedRows", "NFindableCls",
                                         "NSharedCls",  "NSharedCls/NXRows", "NSharedCls/NCls", "Chi2/NCls",    "NSigmaProton",
                                         "NSigmaKaon",  "NSigmaPion",        "Status",          "FirstHit",     "LastHit"};
    std::vector<Int_t> tracks_nbins = {100, 100, 100, 100, 100,
                                       100, 100, 160, 160, 100,
                                       100, 100, 100, 100, 100,
                                       100, 100, 16,  160, 160};
    std::vector<Double_t> tracks_min = {-20., -20., 0., -20., -1.2,
                                        0.,   0.,   0., 0.,   0.,
                                        0.,   0.,   0., 0.,   -10.,
                                        -10,  -10,  0., 0.,   0.};
    std::vector<Double_t> tracks_max = {20.,  20.,  10.,  20.,  1.2,
                                        200., 200., 160., 160., 100.,
                                        100., 1.2,  1.2,  7.,   10,
                                        10,   10,   16,   160., 160.};
    // clang-format on

    for (Int_t sp = 0; sp < (Int_t)species_codes.size(); sp++) {
        fHist_Tracks_Bookkeep[species_codes[sp]] = new TH1F(Form("%s_Bookkeep", species_names[sp].Data()), "", 15, 0., 15.);
        fOutputListOfHists->Add(fHist_Tracks_Bookkeep[species_codes[sp]]);
    }

    for (TString& set : tracks_sets) {
        for (Int_t sp = 0; sp < (Int_t)species_codes.size(); sp++) {
            for (Int_t prop_idx = 0; prop_idx < (Int_t)tracks_props.size(); prop_idx++) {

                TString histName = Form("%s_%s_%s", set.Data(), species_names[sp].Data(), tracks_props[prop_idx].Data());
                std::tuple<TString, Int_t, TString> histKey = std::make_tuple(set, species_codes[sp], tracks_props[prop_idx]);
                fHist_Tracks[histKey] = new TH1F(histName, "", tracks_nbins[prop_idx], tracks_min[prop_idx], tracks_max[prop_idx]);
                fOutputListOfHists->Add(fHist_Tracks[histKey]);
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
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fPrimaryVertex`, `fMagneticField`, `fSourceOfV0s`, `fOutputListOfTrees`,
 `fOutputListOfHists`
*/
void AliQuickTaskTPCProperties::UserExec(Option_t*) {

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

    AliInfoF("fMagneticField = %f", fMagneticField);  // DEBUG

    if (!PassesEventSelection()) return;

    ProcessTracks();

    PostData(1, fOutputListOfHists);
}

/*
 Define track selection cuts.
 - Input: `cuts_option`
*/
void AliQuickTaskTPCProperties::DefineTracksCuts(TString cuts_option) {

    // kMax_NSigma_Pion = 3.;
    // kMax_NSigma_Kaon = 3.;
    // kMax_NSigma_Proton = 3.;

    // kMax_Track_Eta = 1.;
    // kMin_Track_NTPCClusters = 50;
    // kMax_Track_Chi2PerNTPCClusters = 7.;
    kTurnedOn_Track_StatusCuts = kTRUE;
    // kTurnedOn_Track_RejectKinks = kFALSE;  // PENDING: need to study further
    // kMin_Track_DCAwrtPV = 2.;

    // kMin_Track_Pt[2212] = 0.4;
    // kMin_Track_Pt[321] = 0.0;  // PENDING: to check when analyzing the next reaction channels
    // kMin_Track_Pt[211] = 0.2;
}

/*
 Apply event selection
*/
Bool_t AliQuickTaskTPCProperties::PassesEventSelection() {

    /* MC PV z-vertex range */

    if (fIsMC) {
        if (TMath::Abs(fMC_PrimaryVertex->GetZ()) > 10.) return kFALSE;
    }
    AliInfo("!! MC Event Passed z-vertex range on PV !!");  // DEBUG

    /* trigger */

    Bool_t MB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7);  // kAnyINT
    if (!MB) {
        AliInfoF("!! Event Rejected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);  // DEBUG
        return kFALSE;
    }
    AliInfoF("!! Event Selected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);  // DEBUG

    /* rec. PV z-vertex range */

    if (TMath::Abs(fPrimaryVertex->GetZ()) > 10.) return kFALSE;
    AliInfo("!! Event Passed z-vertex range on PV !!");  // DEBUG

    return kTRUE;
}

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliQuickTaskTPCProperties::ProcessTracks() {

    AliESDtrack* track;

    AliMCParticle* mcPart;

    Int_t mcIdx;
    Int_t mcPdgCode;
    Bool_t mcIsSecondary;
    Float_t mcOriginRadius;

    Float_t dca_xy_wrt_pv, dca_z_wrt_pv;
    Float_t n_tpc_clusters;
    Float_t n_crossed_rows;
    Float_t chi2_over_nclusters;
    Float_t n_findable_clusters, n_shared_clusters;
    Float_t ratio_sharedcls_nxrows, ratio_sharedcls_ncls;
    Float_t chi2_over_ncls;
    Int_t first_hit, last_hit;

    Float_t n_sigma_proton;
    Float_t n_sigma_kaon;
    Float_t n_sigma_pion;

    /* For debugging */

    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};
    TString StatusName[16] = {"kITSin",  "kITSout",     "kITSrefit",    "kITSpid",   "kTPCin",   "kTPCout",   "kTPCrefit",  "kTPCpid",
                              "kITSupg", "kSkipFriend", "kGlobalMerge", "kMultInV0", "kMultSec", "kEmbedded", "kITSpureSA", "kESDpid"};

    Int_t counter_secfrommaterial_fh_nonzero = 0;
    Int_t counter_printed_tracks = 0;
    Int_t counter_radiusaboveIB = 0;
    Int_t counter_radiusaboveIB_fh_nonzero = 0;
    Int_t counter_radiusbelowIB_fh_nonzero = 0;
    Int_t counter_rs_sec = 0;
    Int_t counter_rs_sec_but_no_mc_sec = 0;

    /** Loop over tracks in a single event **/

    for (Int_t esdIdxTrack = 0; esdIdxTrack < fESD->GetNumberOfTracks(); esdIdxTrack++) {

        /* Get track */

        track = fESD->GetTrack(esdIdxTrack);

        /* Get MC info */

        mcIdx = TMath::Abs(track->GetLabel());
        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);
        mcPdgCode = mcPart->PdgCode();
        mcIsSecondary = mcPart->IsSecondaryFromMaterial();  //  || mcPart->IsSecondaryFromWeakDecay();
        mcOriginRadius = TMath::Sqrt(mcPart->Xv() * mcPart->Xv() + mcPart->Yv() * mcPart->Yv());

        /* Apply track selection */

        if (!PassesTrackSelection(track)) continue;

        track->GetImpactParameters(dca_xy_wrt_pv, dca_z_wrt_pv);  // pre-calculated DCA w.r.t. PV

        if ((TMath::Abs(mcPdgCode) == 2212 || TMath::Abs(mcPdgCode) == 321 || TMath::Abs(mcPdgCode) == 211)) {
            // } &&  mcIsSecondary && mcOriginRadius > 84.8) {
            AliInfoF(">> (Selected) Track %i <<", esdIdxTrack);
            AliInfo("");
            AliInfo("   MC Information");
            AliInfo("   ==============");
            AliInfoF("   MC Index: %i, PDG Code: %i", mcIdx, mcPdgCode);
            // AliInfoF("   Origin (x, y, z): (%.3f, %.3f, %.3f)", mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
            AliInfoF("   Origin Radius: %.3f", mcOriginRadius);
            AliInfoF("   IsSecFromMaterial: %i, IsSecFromWeakDecay: %i", (Int_t)mcPart->IsSecondaryFromMaterial(),
                     (Int_t)mcPart->IsSecondaryFromWeakDecay());
            // AliInfoF("   (Apparent) First Row on TPC: %i", GetFirstTPCRow(mcOriginRadius));
            AliInfo("");
            AliInfo("   Rec. Information");
            AliInfo("   ================");
            AliInfoF("   DCAxy, DCAz = %.3f, %.3f", dca_xy_wrt_pv, dca_z_wrt_pv);
            AliESDVertex* vtx = (AliESDVertex*)fESD->GetPrimaryVertexTracks();
            Bool_t rs_sec = IsRubenSecondary(track, vtx);
            AliInfoF("   Ruben Secondary = %i", (Int_t)rs_sec);
            AliInfo("");
            AliInfo("   ITS Information");
            AliInfo("   ===============");
            AliInfoF("   Has point on layers (0,1,2,3,4,5): (%i, %i, %i, %i, %i, %i)", (Int_t)track->HasPointOnITSLayer(0),
                     (Int_t)track->HasPointOnITSLayer(1), (Int_t)track->HasPointOnITSLayer(2), (Int_t)track->HasPointOnITSLayer(3),
                     (Int_t)track->HasPointOnITSLayer(4), (Int_t)track->HasPointOnITSLayer(5));
            // AliInfoF("   N Clusters: %c", track->GetITSNcls());
            // AliInfoF("   ITS Chi2: %f", track->GetITSchi2());
            // AliInfoF("N Clusters: %c", track->GetITSclusters()); // idk yet
            AliInfo("");
            AliInfo("   TPC Information");
            AliInfo("   ===============");
            AliInfoF("   N Clusters: %i", (Int_t)track->GetTPCNcls());
            // AliInfoF("   N Shared Clusters: %i", (Int_t)track->GetTPCnclsS());
            // AliInfoF("   N Findable Clusters: %i", (Int_t)track->GetTPCNclsF());
            AliInfoF("   N Crossed Rows: %f", track->GetTPCCrossedRows());
            AliInfoF("   TPC Chi2 / N Clusters: %f", track->GetTPCchi2() / (Double_t)track->GetTPCNcls());
            Int_t firstHit, lastHit;
            GetTracksHits(track, firstHit, lastHit);
            AliInfoF("   (Apparent) First Hit: %i", firstHit);
            if (mcOriginRadius > 84.8 && firstHit == 0) AliInfo("   !! Radius above TPC IB, but FH is zero !!");
            if (mcOriginRadius < 84.8 && firstHit > 0) AliInfo("   !! Radius below TPC IB, but FH is non-zero !!");
            AliInfoF("   (Apparent) Last Hit: %i", lastHit);
            AliInfoF("   TPC Points: (%f, %f, %f, %f)", track->GetTPCPoints(0), track->GetTPCPoints(1), track->GetTPCPoints(2),
                     track->GetTPCPoints(3));
            AliInfo("");
            // AliInfo("   Track Status");
            // AliInfo("   ============");
            // for (Int_t i = 0; i < 16; i++) {
            // if ((track->GetStatus() & StatusCollection[i])) AliInfoF("   track->GetStatus() & %s > 0", StatusName[i].Data());
            // }
            // AliInfo("");
            AliInfo("");
            if (mcPart->IsSecondaryFromMaterial() && firstHit > 0) counter_secfrommaterial_fh_nonzero++;
            if (mcOriginRadius > 84.8) counter_radiusaboveIB++;
            if (mcOriginRadius > 84.8 && firstHit > 0) counter_radiusaboveIB_fh_nonzero++;
            if (mcOriginRadius < 84.8 && firstHit > 0) counter_radiusbelowIB_fh_nonzero++;
            if (rs_sec) counter_rs_sec++;
            if (rs_sec && !mcPart->IsSecondaryFromMaterial() && !mcPart->IsSecondaryFromWeakDecay()) counter_rs_sec_but_no_mc_sec++;
            counter_printed_tracks++;
        }

        /*
        track->GetTPCClusterMap().CountBits();
        Double_t FoundOverFindable = (nfindableTPC ? static_cast<Float_t>(nclustersTPC) / static_cast<Float_t>(nfindableTPC) : 0);
        // TPC Quality
        if (!(track->GetTPCNcls() >= 110 && track->GetTPCsignalN() >= 80 && FoundOverFindable > 0.6 &&
              track->GetTPCchi2() / Double_t(nclustersTPC) < 4.))
            return kFALSE;
        // ITS Quality
        if (track->GetITSchi2() / double(track->GetITSNcls()) > 10.) return kFALSE;
        UChar_t ITSShared = track->GetITSSharedClusterMap();
        Int_t nSharedITS = 0;
        for (int i = 0; i < 6; i++)
            if ((ITSShared >> i) & 1) nSharedITS++;
        if (nSharedITS > 3) return kFALSE;

        if (!track1->GetInnerXYZ(tpcEnt1)) continue;
        clu1 = track1->GetTPCClusterMap();
        sha1 = track1->GetTPCSharedMap();
        Int_t label = vTrack->GetTPCLabel();  // Use TPC-only label for TPC-only resolution analysis
        if (label <= 0) return;

        // get the first TPC track reference
        AliTrackReference* ref0 = GetFirstTPCTrackRef(mcParticle);
        if (!ref0) return;
        */

        /* Get track properties */

        n_tpc_clusters = track->GetTPCNcls();  // NOTE: capital N
        n_crossed_rows = track->GetTPCCrossedRows();
        n_findable_clusters = (Int_t)track->GetTPCNclsF();
        n_shared_clusters = (Int_t)track->GetTPCnclsS();
        ratio_sharedcls_nxrows = n_crossed_rows ? (Double_t)n_shared_clusters / n_crossed_rows : 999.;
        ratio_sharedcls_ncls = n_tpc_clusters ? (Double_t)n_shared_clusters / (Double_t)n_tpc_clusters : 999.;
        chi2_over_ncls = n_tpc_clusters ? track->GetTPCchi2() / (Double_t)n_tpc_clusters : 999.;
        GetTracksHits(track, first_hit, last_hit);

        n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

        /* Fill histograms */

        f2DHist_Tracks_TPCsignal["Selected"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

        std::vector<Int_t> possiblePID;
        if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton && track->Charge() > 0.) possiblePID.push_back(2212);
        if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton && track->Charge() < 0.) possiblePID.push_back(-2212);
        if (TMath::Abs(n_sigma_kaon) < kMax_NSigma_Kaon && track->Charge() > 0.) possiblePID.push_back(321);
        if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion && track->Charge() < 0.) possiblePID.push_back(-211);
        if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion && track->Charge() > 0.) possiblePID.push_back(211);

        for (Int_t& esdPdgCode : possiblePID) {

            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSharedCls/NXRows")]->Fill(ratio_sharedcls_nxrows);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSharedCls/NCls")]->Fill(ratio_sharedcls_ncls);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "FirstHit")]->Fill((Float_t)first_hit);
            fHist_Tracks[std::make_tuple("Selected", esdPdgCode, "LastHit")]->Fill((Float_t)last_hit);
            PlotStatus(track, "Selected", esdPdgCode);

            if (mcPdgCode != esdPdgCode) continue;

            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSharedCls/NXRows")]->Fill(ratio_sharedcls_nxrows);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSharedCls/NCls")]->Fill(ratio_sharedcls_ncls);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "FirstHit")]->Fill((Float_t)first_hit);
            fHist_Tracks[std::make_tuple("True", esdPdgCode, "LastHit")]->Fill((Float_t)last_hit);
            PlotStatus(track, "True", esdPdgCode);
            f2DHist_Tracks_TPCsignal["True"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());

            if (!mcIsSecondary) continue;

            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Px")]->Fill(track->Px());
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Py")]->Fill(track->Py());
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Pt")]->Fill(track->Pt());
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Pz")]->Fill(track->Pz());
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Eta")]->Fill(track->Eta());
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "DCAxy_wrtPV")]->Fill(dca_xy_wrt_pv);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "DCAz_wrtPV")]->Fill(dca_z_wrt_pv);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NTPCClusters")]->Fill(n_tpc_clusters);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NCrossedRows")]->Fill(n_crossed_rows);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NFindableCls")]->Fill(n_findable_clusters);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSharedCls")]->Fill(n_shared_clusters);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSharedCls/NXRows")]->Fill(ratio_sharedcls_nxrows);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSharedCls/NCls")]->Fill(ratio_sharedcls_ncls);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "Chi2/NCls")]->Fill(chi2_over_ncls);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSigmaProton")]->Fill(n_sigma_proton);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSigmaKaon")]->Fill(n_sigma_kaon);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "NSigmaPion")]->Fill(n_sigma_pion);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "FirstHit")]->Fill((Float_t)first_hit);
            fHist_Tracks[std::make_tuple("Secondary", esdPdgCode, "LastHit")]->Fill((Float_t)last_hit);
            PlotStatus(track, "Secondary", esdPdgCode);
            f2DHist_Tracks_TPCsignal["Secondary"]->Fill(track->P() * track->GetSign(), track->GetTPCsignal());
        }  // end of loop over PID hypotheses
    }      // end of loop over tracks

    AliInfo("## Summary ##");
    AliInfo("=============");
    AliInfoF(">> N Selected Tracks: %i", counter_printed_tracks);
    AliInfoF(">> N Tracks with SecFromMaterial && FirstHit > 0: %i", counter_secfrommaterial_fh_nonzero);
    AliInfoF(">> N Tracks with OriginRadius > 84.8 cm: %i", counter_radiusaboveIB);
    AliInfoF(">> N Tracks with OriginRadius > 84.8 cm && FirstHit > 0: %i", counter_radiusaboveIB_fh_nonzero);
    AliInfoF(">> N Tracks with OriginRadius < 84.8 cm && FirstHit > 0: %i", counter_radiusbelowIB_fh_nonzero);
    AliInfoF(">> N Tracks with Ruben Secondary: %i", counter_rs_sec);
    AliInfoF(">> N Tracks with Ruben Secondary but no MC Secondary: %i", counter_rs_sec_but_no_mc_sec);
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliQuickTaskTPCProperties::PassesTrackSelection(AliESDtrack* track) {

    /* Particle ID */

    Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

    std::vector<Int_t> possiblePID;
    if (TMath::Abs(n_sigma_proton) < kMax_NSigma_Proton) possiblePID.push_back(2212);
    if (TMath::Abs(n_sigma_kaon) < kMax_NSigma_Kaon) possiblePID.push_back(321);
    if (TMath::Abs(n_sigma_pion) < kMax_NSigma_Pion) possiblePID.push_back(211);
    if ((kMax_NSigma_Proton || kMax_NSigma_Kaon || kMax_NSigma_Pion) && !possiblePID.size()) return kFALSE;

    Float_t Eta = TMath::Abs(track->Eta());
    if (kMax_Track_Eta && Eta > kMax_Track_Eta) return kFALSE;

    Int_t NTPCClusters = (Int_t)track->GetTPCNcls();
    if (kMin_Track_NTPCClusters && NTPCClusters < kMin_Track_NTPCClusters) return kFALSE;

    Float_t Chi2PerNTPCClusters = NTPCClusters > 0 ? track->GetTPCchi2() / (Double_t)NTPCClusters : 999.;
    if (kMax_Track_Chi2PerNTPCClusters && Chi2PerNTPCClusters > kMax_Track_Chi2PerNTPCClusters) return kFALSE;

    // >> TPC and ITS status
    Bool_t tpc_status =
        (track->GetStatus() & AliESDtrack::kTPCin) && (track->GetStatus() & AliESDtrack::kTPCout) && (track->GetStatus() & AliESDtrack::kTPCrefit);
    Bool_t its_status =
        !(track->GetStatus() & AliESDtrack::kITSin) && !(track->GetStatus() & AliESDtrack::kITSout) && !(track->GetStatus() & AliESDtrack::kITSrefit);
    if (kTurnedOn_Track_StatusCuts && (!tpc_status || !its_status)) return kFALSE;

    // No hit on ITS
    Bool_t its_layers = track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1) && track->HasPointOnITSLayer(2) &&
                        track->HasPointOnITSLayer(3) && track->HasPointOnITSLayer(4) && track->HasPointOnITSLayer(5);
    if (kTurnedOn_Track_StatusCuts && its_layers) return kFALSE;

    // >> reject kinks
    if (kTurnedOn_Track_RejectKinks && track->GetKinkIndex(0) != 0) return kFALSE;

    // >> DCA w.r.t Primary Vertex
    Float_t DCAxy_wrtPV, DCAz_wrtPV;
    track->GetImpactParameters(DCAxy_wrtPV, DCAz_wrtPV);
    // Float_t DCAwrtPV = TMath::Sqrt(xy_impar * xy_impar + z_impar * z_impar);
    if (kMin_Track_DCAxy_wrtPV && DCAxy_wrtPV < kMin_Track_DCAxy_wrtPV) return kFALSE;
    if (kMin_Track_DCAz_wrtPV && DCAz_wrtPV < kMin_Track_DCAz_wrtPV) return kFALSE;

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
void AliQuickTaskTPCProperties::PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode) {

    ULong64_t StatusCollection[16] = {AliESDtrack::kITSin,   AliESDtrack::kITSout,     AliESDtrack::kITSrefit,    AliESDtrack::kITSpid,
                                      AliESDtrack::kTPCin,   AliESDtrack::kTPCout,     AliESDtrack::kTPCrefit,    AliESDtrack::kTPCpid,
                                      AliESDtrack::kITSupg,  AliESDtrack::kSkipFriend, AliESDtrack::kGlobalMerge, AliESDtrack::kMultInV0,
                                      AliESDtrack::kMultSec, AliESDtrack::kEmbedded,   AliESDtrack::kITSpureSA,   AliESDtrack::kESDpid};

    for (Int_t i = 0; i < 16; i++) {
        if ((track->GetStatus() & StatusCollection[i])) fHist_Tracks[std::make_tuple(set, esdPdgCode, "Status")]->Fill(i);
    }
}

/*
 (taken from: PWGJE/AliPWG4HighPtTrackQA.cxx)
*/
void AliQuickTaskTPCProperties::GetTracksHits(AliESDtrack* track, Int_t& firstHit, Int_t& lastHit) {

    TBits fTPCClusterMap = track->GetTPCClusterMap();
    firstHit = 0;
    lastHit = 0;

    for (Int_t i = 0; i <= 159; i++) {
        if (fTPCClusterMap[i] > 0) {
            firstHit = i;
            break;
        }
    }
    for (Int_t i = 159; i >= 0; i--) {
        if (fTPCClusterMap[i] > 0) {
            lastHit = i;
            break;
        }
    }
}

/*
  (taken from `PWGGA/GammaConvBase/AliConversionPhotonCuts.cxx`)
*/
Int_t AliQuickTaskTPCProperties::GetFirstTPCRow(Double_t radius) {

    Int_t firstTPCRow = 0;
    Double_t radiusI = 84.8;
    Double_t radiusO = 134.6;
    Double_t radiusOB = 198.;
    Double_t rSizeI = 0.75;
    Double_t rSizeO = 1.;
    Double_t rSizeOB = 1.5;
    Int_t nClsI = 63;
    Int_t nClsIO = 127;

    if (radius > radiusI && radius <= radiusO) return (Int_t)((radius - radiusI) / rSizeI);
    if (radius > radiusO && radius <= radiusOB) return (Int_t)(nClsI + (radius - radiusO) / rSizeO);
    if (radius > radiusOB) return (Int_t)(nClsIO + (radius - radiusOB) / rSizeOB);

    return 0;
}

/*
Taken from AliTrackletAlg.cxx
*/
Bool_t AliQuickTaskTPCProperties::IsRubenSecondary(AliESDtrack* track, const AliVertex* vtx) {
    // RS: check if the track is primary and set the flag
    Float_t fCutPxDrSPDin = 0.1;    // max P*DR for primaries involving at least 1 SPD
    Float_t fCutPxDrSPDout = 0.15;  // max P*DR for primaries not involving any SPD
    Float_t fCutPxDz = 0.2;         // max P*DZ for primaries
    Float_t fCutDCArz = 0.5;        // max DR or DZ for primares

    Double_t cut = (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) ? fCutPxDrSPDin : fCutPxDrSPDout;
    Float_t xz[2];
    track->GetDZ(vtx->GetX(), vtx->GetY(), vtx->GetZ(), fMagneticField, xz);
    if (TMath::Abs(xz[0] * track->P()) > cut || TMath::Abs(xz[1] * track->P()) > fCutPxDz || TMath::Abs(xz[0]) > fCutDCArz ||
        TMath::Abs(xz[1]) > fCutDCArz) {
        // track->SetStatus(AliESDtrack::kMultSec);
        return kTRUE;
    }
    // else
    // track->ResetStatus(AliESDtrack::kMultSec);
    return kFALSE;
}
