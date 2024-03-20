#include "AliAnalysisQuickTask.h"

ClassImp(AliAnalysisQuickTask);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliAnalysisQuickTask::AliAnalysisQuickTask()
    : AliAnalysisTaskSE(),  //
      fIsMC(0),
      fOutputListOfHists(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0) {}

/*
 Constructor, called locally.
*/
AliAnalysisQuickTask::AliAnalysisQuickTask(const char *name, Bool_t IsMC)
    : AliAnalysisTaskSE(name), fIsMC(IsMC), fOutputListOfHists(0), fMC(0), fESD(0), fPIDResponse(0) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/*
 Destructor.
*/
AliAnalysisQuickTask::~AliAnalysisQuickTask() {
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisQuickTask::UserCreateOutputObjects() {

    /** Add mandatory routines **/

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (!man) {
        AliFatal("ERROR: AliAnalysisManager couldn't be found.");
    }

    AliESDInputHandler *inputHandler = (AliESDInputHandler *)(man->GetInputEventHandler());
    if (!inputHandler) {
        AliFatal("ERROR: AliESDInputHandler couldn't be found.");
    }

    fPIDResponse = inputHandler->GetPIDResponse();

    /** Prepare output histograms */

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    fHist_Tracks_Eta = new TH1F("Eta", "", 100, -0.8, 0.8);
    fOutputListOfHists->Add(fHist_Tracks_Eta);

    fHist_Tracks_Status = new TH1F("Status", "", 20, 0., 20);
    fOutputListOfHists->Add(fHist_Tracks_Status);

    /* From SPD */

    hTracklets = new TH1I("hNtracklets", "Tracklet distribution", 800, 0, 8000);
    hTracklets->SetXTitle("# Tracklets");
    fOutputListOfHists->Add(hTracklets);

    hSPDphivsSPDeta = new TH2F("hSPDphivsSPDeta", "Tracklets - #varphi vs #eta", 120, -3., 3, 360, 0., 2 * TMath::Pi());
    hSPDphivsSPDeta->GetXaxis()->SetTitle("Pseudorapidity #eta");
    hSPDphivsSPDeta->GetYaxis()->SetTitle("#varphi [rad]");
    fOutputListOfHists->Add(hSPDphivsSPDeta);

    hVertexZ = new TH1F("hVertexZ", "Vertex Z distribution", 300, -15, 15);
    hVertexZ->SetXTitle("Z Vertex [cm]");
    fOutputListOfHists->Add(hVertexZ);

    hClu0VsTracklet = new TH2F("hClu0VsTracklet", "Clusters Layer 0 Vs Tracklets", 400, 0, 8000, 400, 0, 1.5 * 8000);
    hClu0VsTracklet->SetYTitle("clusters SPD Layer 0");
    hClu0VsTracklet->SetXTitle("tracklets");
    fOutputListOfHists->Add(hClu0VsTracklet);

    hClu1VsTracklet = new TH2F("hClu1VsTracklet", "Clusters Layer 1 Vs Tracklets", 400, 0, 8000, 400, 0, 1.5 * 8000);
    hClu1VsTracklet->SetYTitle("clusters SPD Layer 1");
    hClu1VsTracklet->SetXTitle("tracklets");
    fOutputListOfHists->Add(hClu1VsTracklet);

    hEventsProcessed = new TH1I("hEventsProcessed", "Number of processed events", 1, 0, 1);
    fOutputListOfHists->Add(hEventsProcessed);

    hSPDmultiplicityVsVertexZ = new TH2F("hSPDmultiplicityVsVertexZ", "SPD Trackles vs zvertex", 300, -15, 15, 800, 0, 8000);
    hSPDmultiplicityVsVertexZ->GetYaxis()->SetTitle("SPD multiplicity");
    hSPDmultiplicityVsVertexZ->GetXaxis()->SetTitle("zvertex position");
    fOutputListOfHists->Add(hSPDmultiplicityVsVertexZ);

    /* From CheckESDTracks */

    fHistNEvents = new TH1F("hNEvents", "Number of processed events", 15, -0.5, 14.5);
    fHistNEvents->SetMinimum(0);
    fHistNEvents->GetXaxis()->SetBinLabel(1, "All events");
    fHistNEvents->GetXaxis()->SetBinLabel(2, "PhysSel");
    fHistNEvents->GetXaxis()->SetBinLabel(3, "InCentralityClass");
    fHistNEvents->GetXaxis()->SetBinLabel(4, "Good vertex");
    fHistNEvents->GetXaxis()->SetBinLabel(5, "Pass zSPD-zTrk vert sel");
    fHistNEvents->GetXaxis()->SetBinLabel(6, "|zvert|<10");
    fHistNEvents->GetXaxis()->SetBinLabel(7, "Pileup cut");
    fHistNEvents->GetXaxis()->SetBinLabel(8, "Generated pileup");
    fOutputListOfHists->Add(fHistNEvents);
    /*  */
    hDCAxy = new TH1F("DCAxy", "", 100, -5, 5);
    hDCAxy->SetXTitle("DCAxy (cm)");
    fOutputListOfHists->Add(hDCAxy);

    hDCAz = new TH1F("DCAz", "", 100, -5, 5);
    hDCAz->SetXTitle("DCAz (cm)");
    fOutputListOfHists->Add(hDCAz);
    /*  */
    hNTracks = new TH1F("hNTracks", "Number of tracks in ESD events", 2000, 0., 20000);
    fOutputListOfHists->Add(hNTracks);

    hPt = new TH1F("hPt", "", 2000, 0., 200);
    fOutputListOfHists->Add(hPt);

    hNTracksPerSelectedEvent = new TH1F("hNTracks", "Number of tracks in selected events", 2000, 0., 20000);
    fOutputListOfHists->Add(hNTracksPerSelectedEvent);
    /* */
    hITSLayerVsPhi = new TH2F("hITSLayerVsPhi", "Hit in the ITS Layer vs #phi", 300, 0., 2 * TMath::Pi(), 6, 0., 6.);
    fOutputListOfHists->Add(hITSLayerVsPhi);

    hITSLayerVsEta = new TH2F("hITSLayerVsEta", "Hit in the ITS Layer vs #eta", 40, -1., 1., 6, 0., 6.);
    fOutputListOfHists->Add(hITSLayerVsEta);
    /*  */
    hTPCNSigmaProtonVsInnerParamP =
        new TH2F("hTPCNSigmaProtonVsInnerParamP", "TPC N_{#sigma} Proton vs p^{in} (GeV/c)", 200, 0., 10., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaProtonVsInnerParamP);

    hTPCNSigmaKaonVsInnerParamP = new TH2F("hTPCNSigmaKaonVsInnerParamP", "TPC N_{#sigma} Kaon vs p^{in} (GeV/c)", 200, 0., 10., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaKaonVsInnerParamP);

    hTPCNSigmaPionVsInnerParamP = new TH2F("hTPCNSigmaPionVsInnerParamP", "TPC N_{#sigma} Pion vs p^{in} (GeV/c)", 200, 0., 10., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaPionVsInnerParamP);
    /*  */
    hTPCNSigmaProtonVsEta = new TH2F("hTPCNSigmaProtonVsEta", "TPC N_{#sigma} Proton vs #eta", 40, -1., 1., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaProtonVsEta);

    hTPCNSigmaKaonVsEta = new TH2F("hTPCNSigmaKaonVsEta", "TPC N_{#sigma} Kaon vs #eta", 40, -1., 1., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaKaonVsEta);

    hTPCNSigmaPionVsEta = new TH2F("hTPCNSigmaPionVsEta", "TPC N_{#sigma} Pion vs #eta", 40, -1., 1., 200, -10., 10.);
    fOutputListOfHists->Add(hTPCNSigmaPionVsEta);
    /*  */
    hTPCoutTracks = new TH1F("hTPCoutTracks", "Number of TPCout tracks in ESD events", 2000, 0., 20000);
    fOutputListOfHists->Add(hTPCoutTracks);

    hTPCclusters = new TH1F("hTPCclusters", "Number of TPC clusters in ESD events", 2000, 0., 5000000);
    fOutputListOfHists->Add(hTPCclusters);

    PostData(1, fOutputListOfHists);
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fMagneticField`, `fPrimaryVertex`, `fSourceOfV0s`, `fReactionChannel`, `fOutputListOfTrees`,
 `fOutputListOfHists`
*/
void AliAnalysisQuickTask::UserExec(Option_t *) {

    fESD = dynamic_cast<AliESDEvent *>(InputEvent());

    if (!fESD) {
        AliFatal("ERROR: AliESDEvent couldn't be found.");
    }

    /*  */

    Int_t ntracks = fESD->GetNumberOfTracks();
    hNTracks->Fill(ntracks);

    Bool_t isPhysSel =
        (((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAnyINT);
    if (isPhysSel) hNTracksPerSelectedEvent->Fill(ntracks);  // TEMPORARY

    hTPCoutTracks->Fill(fESD->GetNumberOfTPCTracks());
    hTPCclusters->Fill(fESD->GetNumberOfTPCClusters());

    const AliESDVertex *vertex = fESD->GetPrimaryVertexSPD();
    const AliMultiplicity *mult = fESD->GetMultiplicity();

    if (!vertex) return;
    if (!vertex->GetStatus()) return;
    if (vertex->GetNContributors() < 1) return;

    hEventsProcessed->Fill(0);

    hTracklets->Fill(mult->GetNumberOfTracklets());

    Double_t vtxPos[3] = {999, 999, 999};
    vertex->GetXYZ(vtxPos);
    hVertexZ->Fill(vtxPos[2]);

    hSPDmultiplicityVsVertexZ->Fill(vtxPos[2], mult->GetNumberOfTracklets());
    hClu0VsTracklet->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(0));
    hClu1VsTracklet->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(1));

    for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++) {
        Float_t phiTr = mult->GetPhi(iTracklet);
        Float_t etaTr = mult->GetEta(iTracklet);
        hSPDphivsSPDeta->Fill(etaTr, phiTr);
    }

    ProcessTracks();

    // fMC = MCEvent();
    // if (!fMC) {
    // AliFatal("ERROR: AliMCEvent couldn't be found.");
    // }
    // fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());
    // fMagneticField = fESD->GetMagneticField();
    // fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());
    // DefineTracksCuts("");
    // if (fIsMC) ProcessMCGen();
    // getPdgCode_fromMcIdx.clear();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfHists);
}

/*
 Loop over MC particles in a single event. Store the indices of the signal particles.
 - Uses: `fMC`, `fPDG`, `fMC_PrimaryVertex`
*/
void AliAnalysisQuickTask::ProcessMCGen() {

    AliMCParticle *mcPart;
    Int_t pdg_mc;

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle *)fMC->GetTrack(mcIdx);
        pdg_mc = mcPart->PdgCode();

        if (pdg_mc != 2212) continue;
        if (!mcPart->IsPhysicalPrimary()) continue;
    }
}

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisQuickTask::ProcessTracks() {

    AliESDtrack *track;

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

    for (Int_t esdIdx = 0; esdIdx < fESD->GetNumberOfTracks(); esdIdx++) {

        /* Get track */

        track = fESD->GetTrack(esdIdx);

        hPt->Fill(track->Pt());

        Float_t xy_impar_wrt_pv, z_impar_wrt_pv;
        track->GetImpactParameters(xy_impar_wrt_pv, z_impar_wrt_pv);  // pre-calculated DCA w.r.t. PV
        hDCAxy->Fill(xy_impar_wrt_pv);
        hDCAz->Fill(z_impar_wrt_pv);

        // loop over the ITS hits of a track
        for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
            if (track->HasPointOnITSLayer(iLayer)) {
                hITSLayerVsPhi->Fill(track->Phi(), iLayer);
                hITSLayerVsEta->Fill(track->Eta(), iLayer);
            }
        }

        Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

        if (TMath::Abs(n_sigma_proton) < 3.) {
            hTPCNSigmaProtonVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_proton);
            hTPCNSigmaProtonVsEta->Fill(track->Eta(), n_sigma_proton);
        }
        if (TMath::Abs(n_sigma_kaon) < 3.) {
            hTPCNSigmaKaonVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_kaon);
            hTPCNSigmaKaonVsEta->Fill(track->Eta(), n_sigma_kaon);
        }
        if (TMath::Abs(n_sigma_pion) < 3.) {
            hTPCNSigmaPionVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_pion);
            hTPCNSigmaPionVsEta->Fill(track->Eta(), n_sigma_pion);
        }

        /* Previously */
        /*
        mcIdx = TMath::Abs(track->GetLabel());
        mcPdgCode = getPdgCode_fromMcIdx[mcIdx];
        if (!PassesTrackSelection(track)) continue;
        if (mcPdgCode != 2212) continue;
        fHist_Tracks_Eta->Fill(track->Eta());
        PlotStatus(track);
        */

    }  // end of loop over tracks
}
