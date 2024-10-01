#include "include/Headers.hxx"
#include "include/Utilities.hxx"

/*
 * Merge all files of a directory as histograms
 */
void Merger_QA(  //
                 // TString input_dir = "/home/ceres/borquez/some/output/A_bkg_15o_latest",  //
                 // TString pattern = "AnalysisResults_24514*.root",                         //
    TString input_dir = "/home/ceres/borquez/some/output/A_bkg_18qr_latest",  //
    TString pattern = "AnalysisResults_2*.root",                              //
    TString output_dir = "output") {

    TDatabasePDG DB_PDG;
    TString DirBaseName = gSystem->BaseName(input_dir);

    TString InputListName;

    /* Prepare output */

    system("mkdir -p " + output_dir);
    TString OutputFilename = output_dir + "/QA_" + DirBaseName + ".root";
    TFile *OutputFile = new TFile(OutputFilename, "RECREATE");

    TList *OutputList_QA_Hists = new TList();

    /***               ***/
    /*** Case 1: Trees ***/
    /***               ***/

    if (!input_dir.Contains("bkg")) {

        InputListName = "Trees";
        TString InputTreeName = "Injected";

        TList *ListOfTrees = new TList();

        TString TreePath = input_dir + "/" + pattern + "/" + InputListName + "/" + InputTreeName;
        std::cout << "TreePath = " << TreePath << std::endl;
        AddObjects(ListOfTrees, TreePath.Data());

        /* Prepare histograms */

        TH1F *hInjected_Pt = new TH1F("hInjected_Pt", ";p_{T} (GeV/c);Counts", 200, -0.5, 5.5);
        TH1F *hInjected_Phi = new TH1F("hInjected_Phi", ";#phi (rad);Counts", 200, -TMath::Pi() - 1., TMath::Pi() + 1.);
        TH1F *hInjected_Rapidity = new TH1F("hInjected_Rapidity", ";Rapidity;Counts", 200, -0.9, 0.9);

        TH1F *hNucleon_P = new TH1F("hNucleon_P", ";p (GeV/c);Counts", 200, 0.2, 0.3);
        TH1F *hNucleon_Phi = new TH1F("hNucleon_Phi", ";#phi (rad);Counts", 200, -TMath::Pi() - 1., TMath::Pi() + 1.);
        TH1F *hNucleon_Theta = new TH1F("hNucleon_Theta", ";#theta (rad);Counts", 200, 0., TMath::Pi());
        TH1F *hNucleon_Px = new TH1F("hNucleon_Px", ";p_{x} (GeV/c);Counts", 200, -0.3, 0.3);
        TH1F *hNucleon_Py = new TH1F("hNucleon_Py", ";p_{y} (GeV/c);Counts", 200, -0.3, 0.3);
        TH1F *hNucleon_Pz = new TH1F("hNucleon_Pz", ";p_{z} (GeV/c);Counts", 200, -0.3, 0.3);

        /* Prepare variables */

        Float_t Injected_Px;
        Float_t Injected_Py;
        Float_t Injected_Pz;
        Float_t Injected_M;

        TLorentzVector Sexaquark;

        Int_t Nucleon_PdgCode;
        Float_t Nucleon_Px;
        Float_t Nucleon_Py;
        Float_t Nucleon_Pz;
        Float_t Nucleon_M;

        TLorentzVector Nucleon;

        /* Loop over trees */

        TIter next(ListOfTrees);
        while (TObject *obj = next()) {

            TTree *CurrentTree = (TTree *)obj;

            /* Set branches */

            CurrentTree->SetBranchAddress("Px", &Injected_Px);
            CurrentTree->SetBranchAddress("Py", &Injected_Py);
            CurrentTree->SetBranchAddress("Pz", &Injected_Pz);
            CurrentTree->SetBranchAddress("M", &Injected_M);

            CurrentTree->SetBranchAddress("PdgCode_Nucleon", &Nucleon_PdgCode);
            CurrentTree->SetBranchAddress("Px_Nucleon", &Nucleon_Px);
            CurrentTree->SetBranchAddress("Py_Nucleon", &Nucleon_Py);
            CurrentTree->SetBranchAddress("Pz_Nucleon", &Nucleon_Pz);

            /* Fill histograms */

            for (Int_t i = 0; i < CurrentTree->GetEntries(); i++) {
                CurrentTree->GetEntry(i);

                Sexaquark.SetXYZM(Injected_Px, Injected_Py, Injected_Pz, Injected_M);

                Nucleon_M = DB_PDG.GetParticle(Nucleon_PdgCode)->Mass();
                Nucleon.SetXYZM(Nucleon_Px, Nucleon_Py, Nucleon_Pz, Nucleon_M);

                hInjected_Pt->Fill(Sexaquark.Pt());
                hInjected_Phi->Fill(Sexaquark.Phi());
                hInjected_Rapidity->Fill(Sexaquark.Rapidity());

                hNucleon_P->Fill(Nucleon.P());
                hNucleon_Phi->Fill(Nucleon.Phi());
                hNucleon_Theta->Fill(Nucleon.Theta());
                hNucleon_Px->Fill(Nucleon.Px());
                hNucleon_Py->Fill(Nucleon.Py());
                hNucleon_Pz->Fill(Nucleon.Pz());
            }
        }

        /* Save histograms */

        OutputList_QA_Hists->Add(hInjected_Pt);
        OutputList_QA_Hists->Add(hInjected_Phi);
        OutputList_QA_Hists->Add(hInjected_Rapidity);

        OutputList_QA_Hists->Add(hNucleon_P);
        OutputList_QA_Hists->Add(hNucleon_Phi);
        OutputList_QA_Hists->Add(hNucleon_Theta);
        OutputList_QA_Hists->Add(hNucleon_Px);
        OutputList_QA_Hists->Add(hNucleon_Py);
        OutputList_QA_Hists->Add(hNucleon_Pz);
    }

    /***                       ***/
    /*** Case 2: QA Hists (1D) ***/
    /***                       ***/

    InputListName = "QA_Hists";
    std::vector<TString> InputHistNames = {"QA_MCGen_VertexZ",   "QA_VertexZ",   "QA_SPD_NTracklets", "QA_NTracks",
                                           "QA_NSelectedTracks", "QA_Tracks_Pt", "QA_Tracks_DCAxy",   "QA_Tracks_DCAz"};
    std::vector<TString> InputHistXAxisTitle = {"MC Gen. PV Zv (cm)", "MC Rec. PV Zv (cm)",   "SPD NTracklets",       "NTracks",
                                                "NSelectedTracks",    "Tracks p_{T} (GeV/c)", "Tracks DCA_{xy} (cm)", "Tracks DCA_{z} (cm)"};

    TString ThisHistName;
    TString HistPath;

    TList *ListOfHists;
    TH1F *ThisHist;

    for (Int_t i = 0; i < (Int_t)InputHistNames.size(); i++) {

        ThisHistName = InputHistNames[i];
        ListOfHists = new TList();

        HistPath = input_dir + "/" + pattern + "/" + InputListName + "/" + ThisHistName;
        AddObjects(ListOfHists, HistPath.Data());

        /* Add up histograms */

        ThisHist = new TH1F("ThisHist", ";;Counts", 1, 0, 1);  // dummy binning, replaced by first hist binning
        ThisHist->SetName(ThisHistName.Data());
        ThisHist->GetXaxis()->SetTitle(InputHistXAxisTitle[i].Data());

        TIter next(ListOfHists);
        while (TObject *obj = next()) {
            TH1F *CurrentHist = (TH1F *)obj;
            if (!ThisHist->GetEntries())
                ThisHist->SetBins(CurrentHist->GetNbinsX(), CurrentHist->GetXaxis()->GetXmin(), CurrentHist->GetXaxis()->GetXmax());
            ThisHist->Add(CurrentHist);
        }

        /* Save histogram */

        OutputList_QA_Hists->Add(ThisHist);

        /* Clean memory */

        delete ListOfHists;
    }

    /***                       ***/
    /*** Case 3: QA Hists (2D) ***/
    /***                       ***/

    InputHistNames = {"QA_SPD_MultiplicityVsVertexZ",
                      "QA_SPD_PhiVsEta",
                      "QA_ITS_LayerNoVsPhi",
                      "QA_ITS_LayerNoVsEta",
                      "QA_TPC_NSigmaProtonVsInnerParamP",
                      "QA_TPC_NSigmaPionVsInnerParamP",
                      "QA_TPC_NSigmaProtonVsEta",
                      "QA_TPC_NSigmaPionVsEta"};
    std::vector<TString> InputHistYAxisTitle = {"SPD # tracklets",  "SPD #phi",           "ITS Layer",      "ITS Layer",
                                                "TPC N #sigma_{p}", "TPC N #sigma_{#pi}", "TPC #sigma_{p}", "TPC #sigma_{#pi}"};
    InputHistXAxisTitle = {"SPD Zv (cm)", "SPD #eta", "#phi", "#eta", "p^{in} (GeV/c)", "p^{in} (GeV/c)", "#eta", "#eta"};

    TH2F *This2DHist;

    for (Int_t i = 0; i < (Int_t)InputHistNames.size(); i++) {

        ThisHistName = InputHistNames[i];
        ListOfHists = new TList();

        HistPath = input_dir + "/" + pattern + "/" + InputListName + "/" + ThisHistName;
        AddObjects(ListOfHists, HistPath.Data());

        /* Add up histograms */

        This2DHist = new TH2F("This2DHist", ";;", 1, 0, 1, 1, 0, 1);  // dummy binning, replaced by first hist binning
        This2DHist->SetName(ThisHistName.Data());
        This2DHist->GetYaxis()->SetTitle(InputHistYAxisTitle[i].Data());
        This2DHist->GetXaxis()->SetTitle(InputHistXAxisTitle[i].Data());

        TIter next(ListOfHists);
        while (TObject *obj = next()) {
            TH2F *CurrentHist = (TH2F *)obj;
            if (!This2DHist->GetEntries()) {
                This2DHist->SetBins(CurrentHist->GetNbinsX(), CurrentHist->GetXaxis()->GetXmin(), CurrentHist->GetXaxis()->GetXmax(),
                                    CurrentHist->GetNbinsY(), CurrentHist->GetYaxis()->GetXmin(), CurrentHist->GetYaxis()->GetXmax());
            }
            This2DHist->Add(CurrentHist);
        }

        /* Save histogram */

        OutputList_QA_Hists->Add(This2DHist);

        /* Clean memory */

        delete ListOfHists;
    }

    /***         ***/
    /*** The end ***/
    /***         ***/

    /* Save lists */

    OutputFile->cd();

    OutputList_QA_Hists->Write("QA_Hists", 1);

    /* Close file */

    OutputFile->Close();
    std::cout << ".. Output file " << OutputFilename << " has been created ." << std::endl;

    /* Clean memory */

    delete OutputFile;
    delete OutputList_QA_Hists;

    delete ThisHist;
    delete This2DHist;
}
