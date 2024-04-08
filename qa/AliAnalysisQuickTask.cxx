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
    : AliAnalysisTaskSE(name),  //
      fIsMC(IsMC),
      fOutputListOfHists(0),
      fMC(0),
      fESD(0),
      fPIDResponse(0) {
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
    if (!man) AliFatal("ERROR: AnalysisManager couldn't be found.");

    fInputHandler = (AliESDInputHandler *)(man->GetInputEventHandler());
    if (!fInputHandler) AliFatal("ERROR: InputEventHandler couldn't be found.");

    fPIDResponse = fInputHandler->GetPIDResponse();

    /** Prepare output histograms */

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    PrepareHistograms();

    PostData(1, fOutputListOfHists);
}

/*
 */
void AliAnalysisQuickTask::PrepareHistograms() {

    /* MC Gen. Event */

    hMCGen_VertexZ = new TH1F("hMCGen_VertexZ", "", 200, -15., 15.);
    fOutputListOfHists->Add(hMCGen_VertexZ);

    /* Event */

    hEvents_Bookkeep = new TH1I("hEvents_Bookkeep", "Number of processed events", 3, 0., 3.);
    fOutputListOfHists->Add(hEvents_Bookkeep);

    hVertexZ = new TH1F("hVertexZ", "", 200, -15., 15.);
    fOutputListOfHists->Add(hVertexZ);

    hNTracks = new TH1F("hNTracks", "", 200, 0., 20000.);
    fOutputListOfHists->Add(hNTracks);

    hNSelectedTracks = new TH1F("hNSelectedTracks", "", 2000, 0., 20000.);
    fOutputListOfHists->Add(hNSelectedTracks);

    /* Tracks */

    hTracks_DCAxy = new TH1F("hTracks_DCAxy", "", 100, -5., 5.);
    fOutputListOfHists->Add(hTracks_DCAxy);

    hTracks_DCAz = new TH1F("hTracks_DCAz", "", 100, -5., 5.);
    fOutputListOfHists->Add(hTracks_DCAz);

    hTracks_Pt = new TH1F("hTracks_Pt", "", 200, 0., 100.);
    fOutputListOfHists->Add(hTracks_Pt);

    hTracks_PhiVsEta = new TH2F("hTracks_PhiVsEta", "", 100, -1., 1., 200, 0., 2. * TMath::Pi());
    fOutputListOfHists->Add(hTracks_PhiVsEta);

    /* SPD */

    hSPD_MultiplicityVsVertexZ = new TH2F("hSPD_MultiplicityVsVertexZ", "", 10, -10., 10., 30, 0, 6000.);
    fOutputListOfHists->Add(hSPD_MultiplicityVsVertexZ);

    hSPD_NTracklets = new TH1I("hSPD_NTracklets", "Tracklet distribution", 60, 0, 6000.);
    fOutputListOfHists->Add(hSPD_NTracklets);

    hSPD_NTrackletsClu0VsNTracklets = new TH2F("hSPD_NTrackletsClu0VsNTracklets", "", 200, 0, 3000., 200, 0, 6000.);
    fOutputListOfHists->Add(hSPD_NTrackletsClu0VsNTracklets);

    hSPD_NTrackletsClu1VsNTracklets = new TH2F("hSPD_NTrackletsClu1VsNTracklets", "", 200, 0, 3000., 200, 0, 6000.);
    fOutputListOfHists->Add(hSPD_NTrackletsClu1VsNTracklets);

    hSPD_PhiVsEta = new TH2F("hSPD_PhiVsEta", "", 400, -1., 1., 400, 0., 2 * TMath::Pi());
    fOutputListOfHists->Add(hSPD_PhiVsEta);

    /* ITS */

    hITS_LayerNoVsPhi = new TH2F("hITS_LayerNoVsPhi", "", 200, 0., 2 * TMath::Pi(), 6, 0., 6.);
    fOutputListOfHists->Add(hITS_LayerNoVsPhi);

    hITS_LayerNoVsEta = new TH2F("hITS_LayerNoVsEta", "", 40, -1., 1., 6, 0., 6.);
    fOutputListOfHists->Add(hITS_LayerNoVsEta);

    /* TPC */

    hTPC_NClusters = new TH1F("hTPC_NClusters", "", 200, 0., 5000000.);
    fOutputListOfHists->Add(hTPC_NClusters);

    hTPC_NSigmaProtonVsInnerParamP = new TH2F("hTPC_NSigmaProtonVsInnerParamP", "", 200, 0., 10., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaProtonVsInnerParamP);

    hTPC_NSigmaKaonVsInnerParamP = new TH2F("hTPC_NSigmaKaonVsInnerParamP", "", 200, 0., 10., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaKaonVsInnerParamP);

    hTPC_NSigmaPionVsInnerParamP = new TH2F("hTPC_NSigmaPionVsInnerParamP", "", 200, 0., 10., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaPionVsInnerParamP);

    hTPC_NSigmaProtonVsEta = new TH2F("hTPC_NSigmaProtonVsEta", "", 200, -1., 1., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaProtonVsEta);

    hTPC_NSigmaKaonVsEta = new TH2F("hTPC_NSigmaKaonVsEta", "", 200, -1., 1., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaKaonVsEta);

    hTPC_NSigmaPionVsEta = new TH2F("hTPC_NSigmaPionVsEta", "", 200, -1., 1., 200, -5., 5.);
    fOutputListOfHists->Add(hTPC_NSigmaPionVsEta);
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
*/
void AliAnalysisQuickTask::UserExec(Option_t *) {

    AliInfo("-- new event --");

    /* MC Generated */

    if (fIsMC) {
        fMC = MCEvent();
        if (!fMC) AliFatal("ERROR: AliMCEvent couldn't be found.");
        fMC_PrimaryVertex = const_cast<AliVVertex *>(fMC->GetPrimaryVertex());
    }

    ProcessMCEvent();

    /* ESD */

    fESD = dynamic_cast<AliESDEvent *>(InputEvent());
    if (!fESD) AliFatal("ERROR: AliESDEvent couldn't be found.");
    fPrimaryVertex = const_cast<AliESDVertex *>(fESD->GetPrimaryVertex());

    ProcessEvent();
    ProcessTracks();

    PostData(1, fOutputListOfHists);
}

/*
 */
void AliAnalysisQuickTask::ProcessMCEvent() {
    fMC_PrimaryVertex = const_cast<AliVVertex *>(fMC->GetPrimaryVertex());
    hMCGen_VertexZ->Fill(fMC_PrimaryVertex->GetZ());
}

/*
 */
void AliAnalysisQuickTask::ProcessEvent() {

    hEvents_Bookkeep->Fill(0);

    Bool_t MB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7);
    if (!MB) {
        AliInfoF("!! Event Rejected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);
        return;
    }
    AliInfoF("!! Event Selected -- %u & %u = %u !!", fInputHandler->IsEventSelected(), AliVEvent::kINT7, MB);

    hEvents_Bookkeep->Fill(1);

    if (TMath::Abs(fPrimaryVertex->GetZ()) > 10.) return;

    hEvents_Bookkeep->Fill(2);

    hVertexZ->Fill(fPrimaryVertex->GetZ());

    Int_t ntracks = fESD->GetNumberOfTracks();
    hNTracks->Fill(ntracks);

    /* SPD */

    const AliESDVertex *spd_vertex = fESD->GetPrimaryVertexSPD();
    const AliMultiplicity *mult = fESD->GetMultiplicity();

    hSPD_NTracklets->Fill(mult->GetNumberOfTracklets());

    hSPD_MultiplicityVsVertexZ->Fill(spd_vertex->GetZ(), mult->GetNumberOfTracklets());
    hSPD_NTrackletsClu0VsNTracklets->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(0));
    hSPD_NTrackletsClu1VsNTracklets->Fill(mult->GetNumberOfTracklets(), mult->GetNumberOfITSClusters(1));

    for (Int_t iTracklet = 0; iTracklet < mult->GetNumberOfTracklets(); iTracklet++) {
        Float_t phiTr = mult->GetPhi(iTracklet);
        Float_t etaTr = mult->GetEta(iTracklet);
        hSPD_PhiVsEta->Fill(etaTr, phiTr);
    }

    /* TPC */

    hTPC_NClusters->Fill(fESD->GetNumberOfTPCClusters());
}

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisQuickTask::ProcessTracks() {

    AliESDtrack *track;

    Int_t nSelectedTracks = 0;

    for (Int_t esdIdx = 0; esdIdx < fESD->GetNumberOfTracks(); esdIdx++) {

        track = fESD->GetTrack(esdIdx);

        for (Int_t iLayer = 0; iLayer < 6; iLayer++) {
            if (track->HasPointOnITSLayer(iLayer)) {
                hITS_LayerNoVsPhi->Fill(track->Phi(), iLayer);
                hITS_LayerNoVsEta->Fill(track->Eta(), iLayer);
            }
        }

        if (TMath::Abs(track->Eta()) > 0.9) continue;

        Float_t xy_impar_wrt_pv, z_impar_wrt_pv;
        track->GetImpactParameters(xy_impar_wrt_pv, z_impar_wrt_pv);  // pre-calculated DCA w.r.t. PV
        hTracks_DCAxy->Fill(xy_impar_wrt_pv);
        hTracks_DCAz->Fill(z_impar_wrt_pv);

        hTracks_Pt->Fill(track->Pt());

        hTracks_PhiVsEta->Fill(track->Eta(), track->Phi());

        Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
        Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
        Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);

        if (TMath::Abs(n_sigma_pion) > 3. && TMath::Abs(n_sigma_kaon) > 3. && TMath::Abs(n_sigma_proton) > 3.) continue;

        nSelectedTracks++;

        if (TMath::Abs(n_sigma_pion) < 3.) {
            hTPC_NSigmaPionVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_pion);
            hTPC_NSigmaPionVsEta->Fill(track->GetInnerParam()->Eta(), n_sigma_pion);
        }
        if (TMath::Abs(n_sigma_kaon) < 3.) {
            hTPC_NSigmaKaonVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_kaon);
            hTPC_NSigmaKaonVsEta->Fill(track->GetInnerParam()->Eta(), n_sigma_kaon);
        }
        if (TMath::Abs(n_sigma_proton) < 3.) {
            hTPC_NSigmaProtonVsInnerParamP->Fill(track->GetInnerParam()->P(), n_sigma_proton);
            hTPC_NSigmaProtonVsEta->Fill(track->GetInnerParam()->Eta(), n_sigma_proton);
        }
    }  // end of loop over tracks

    hNSelectedTracks->Fill(nSelectedTracks);
}
