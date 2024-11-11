#include "include/Headers.hxx"
#include "include/Style.cxx"
#include "include/Utilities.hxx"

/*
 *
 */
void QA_InitialChecks(TString input_files = "/home/ceres/borquez/some/output/A1.8_18qr_latest/AnalysisResults_2*.root",  //
                      TString output_dir = "gfx/") {

    /* Fill chain */

    TString InputListName = "Trees";
    TString InputTreeName = "Injected";

    TList *ListOfTrees = new TList();

    TString TreePath = input_files + "/" + InputListName + "/" + InputTreeName;
    AddTreesToList(ListOfTrees, TreePath.Data());

    /* Prepare histograms */

    TH1F *hInjected_Pt = new TH1F("hInjected_Pt", ";p_{T} (GeV/c);Counts", 200, -0.5, 5.5);
    TH1F *hInjected_Phi = new TH1F("hInjected_Phi", ";#phi (rad);Counts", 200, -TMath::Pi() - 1., TMath::Pi() + 1.);
    TH1F *hInjected_Rapidity = new TH1F("hInjected_Rapidity", ";Rapidity;Counts", 200, -0.9, 0.9);

    /* Prepare variables */

    Float_t Injected_Px;
    Float_t Injected_Py;
    Float_t Injected_Pz;
    Float_t Injected_M;

    TLorentzVector Sexaquark;

    /* Loop over trees */

    TIter next(ListOfTrees);
    for (TObject *obj = next(); obj != nullptr; obj = next()) {

        TTree *CurrentTree = (TTree *)obj;

        /* Set branches */

        CurrentTree->SetBranchAddress("Px", &Injected_Px);
        CurrentTree->SetBranchAddress("Py", &Injected_Py);
        CurrentTree->SetBranchAddress("Pz", &Injected_Pz);
        CurrentTree->SetBranchAddress("M", &Injected_M);

        /* Fill histograms */

        for (Int_t i = 0; i < CurrentTree->GetEntries(); i++) {
            CurrentTree->GetEntry(i);

            Sexaquark.SetXYZM(Injected_Px, Injected_Py, Injected_Pz, Injected_M);

            hInjected_Pt->Fill(Sexaquark.Pt());
            hInjected_Phi->Fill(Sexaquark.Phi());
            hInjected_Rapidity->Fill(Sexaquark.Rapidity());
        }
    }

    /* Draw */

    SetMyStyle();
    gStyle->SetOptStat(0);

    TString output_filename;
    TCanvas *c = new TCanvas("c", "c", 3 * 600, 600);

    // split canvas in 2 pads
    c->Divide(3, 1);

    c->cd(1);
    SetMyHistStyle(hInjected_Pt);
    hInjected_Pt->SetFillStyle(0);
    hInjected_Pt->SetLineWidth(3);
    hInjected_Pt->SetLineColor(myBlue);
    hInjected_Pt->Draw();

    c->cd(2);
    SetMyHistStyle(hInjected_Phi);
    hInjected_Phi->SetFillStyle(0);
    hInjected_Phi->SetLineWidth(3);
    hInjected_Phi->SetLineColor(myOrange);
    hInjected_Phi->Draw();

    c->cd(3);
    SetMyHistStyle(hInjected_Rapidity);
    hInjected_Rapidity->SetFillStyle(0);
    hInjected_Rapidity->SetLineWidth(3);
    hInjected_Rapidity->SetLineColor(myGreen);
    hInjected_Rapidity->Draw();

    c->Print(output_dir + "test.png");
}
