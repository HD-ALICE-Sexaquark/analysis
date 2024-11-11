#include "include/Headers.hxx"
#include "include/Style.cxx"

/*
 *
 */
void QA_TwoDimHists(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputSignalFilenames = {"QA_18qr.root", "QA_15o.root"};
    std::vector<TString> InputBkgFilenames = {"QA_bkg_18qr_latest.root", "QA_bkg_15o_latest.root"};
    TString InputListName = "QA_Hists";

    /* Case 1: 18qr signal */

    std::vector<TString> InputHistNames = {"QA_SPD_MultiplicityVsVertexZ",
                                           "QA_SPD_PhiVsEta",
                                           "QA_ITS_LayerNoVsPhi",
                                           "QA_ITS_LayerNoVsEta",
                                           "QA_TPC_NSigmaProtonVsInnerParamP",
                                           "QA_TPC_NSigmaPionVsInnerParamP",
                                           "QA_TPC_NSigmaProtonVsEta",
                                           "QA_TPC_NSigmaPionVsEta"};

    TFile *InputFile;
    TList *ListOfHists;
    TH2F *ThisHist;

    std::vector<TColor> Colors = {myBlue, myOrange, myGreen, myRed};

    TString HistPath;
    TCanvas *c = new TCanvas("c", "c", 2 * 800, 2 * 800);

    for (Int_t hh = 0; hh < (Int_t)InputHistNames.size(); hh++) {

        c->Divide(2, 2);
        SetMy2DStyle();

        for (Int_t ff = 0; ff < (Int_t)InputSignalFilenames.size(); ff++) {

            /* Signal */

            InputFile = TFile::Open(input_dir + "/" + InputSignalFilenames[ff], "READ");
            ListOfHists = (TList *)InputFile->Get(InputListName);
            ThisHist = (TH2F *)ListOfHists->FindObject(InputHistNames[hh]);

            c->cd(2 * ff + 1);

            SetMy2DHistStyle(ThisHist);
            c->SetMargin(0.15, 0.15, 0.15, 0.15);

            ThisHist->Draw("COLZ");

            /* Bkg */

            InputFile = TFile::Open(input_dir + "/" + InputBkgFilenames[ff], "READ");
            ListOfHists = (TList *)InputFile->Get(InputListName);
            ThisHist = (TH2F *)ListOfHists->FindObject(InputHistNames[hh]);

            c->cd(2 * (ff + 1));

            SetMy2DHistStyle(ThisHist);
            c->SetMargin(0.15, 0.15, 0.15, 0.15);

            ThisHist->Draw("COLZ");
        }

        system("mkdir -p " + output_dir);
        c->Print(output_dir + "/" + InputHistNames[hh] + ".png");

        c->Clear();
    }
}
