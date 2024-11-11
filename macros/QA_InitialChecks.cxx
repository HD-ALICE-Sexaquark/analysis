#include "include/Headers.hxx"
#include "include/Style.cxx"

/*
 * QA, Initial checks:
 * - Injected: Pt, Phi, Rapidity
 * - Struck Nucleon: Pt, Phi, Theta, Px, Py, Pz
 */
void QA_InitialChecks(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputFilenames = {"QA_ChannelA.root", "QA_ChannelD.root", "QA_ChannelE.root", "QA_ChannelH.root"};
    TString InputListName = "QA_Hists";

    /* Prepare objects */

    TFile *InputFile;
    TList *ListOfHists;
    TH1F *ThisHist;

    std::vector<TColor> Colors = {myBlue, myOrange, myGreen, myRed};

    /* Case 1: Injected */

    SetMyStyle();
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas("c1", "c1", 3 * 800, 800);
    c1->Divide(3, 1);

    std::vector<TString> InputHistsNames = {"hInjected_Pt", "hInjected_Phi", "hInjected_Rapidity"};

    TLegend *leg = new TLegend(0.25, 0.2, 0.4, 0.35);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    for (Int_t hh = 0; hh < (Int_t)InputHistsNames.size(); hh++) {

        for (Int_t ff = 0; ff < (Int_t)InputFilenames.size(); ff++) {

            InputFile = TFile::Open(input_dir + "/" + InputFilenames[ff], "READ");
            ListOfHists = (TList *)InputFile->Get(InputListName);
            ThisHist = (TH1F *)ListOfHists->FindObject(InputHistsNames[hh]);

            c1->cd(hh + 1);

            SetMyHistStyle(ThisHist);
            ThisHist->SetLineWidth(3);
            ThisHist->SetLineColor(Colors[ff].GetNumber());

            ThisHist->GetYaxis()->SetRangeUser(0, 1.1 * ThisHist->GetMaximum());

            ThisHist->Draw("SAME");

            if (hh == 0) {
                TString channel = InputFilenames[ff](10, 1);
                leg->AddEntry(ThisHist, "Channel " + channel, "l");
                ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ff].GetNumber());
            }
        }

        if (hh == 0) leg->Draw();
    }

    system("mkdir -p " + output_dir);
    c1->Print(output_dir + "/QA_Injected.png");
    delete c1;
    delete leg;

    /* Case 2: Struck Nucleon */

    TCanvas *c2 = new TCanvas("c2", "c2", 3 * 800, 2 * 800);
    c2->Divide(3, 2);

    leg = new TLegend(0.7, 0.7, 0.85, 0.85);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    InputHistsNames = {"hNucleon_Pt", "hNucleon_Phi", "hNucleon_Theta", "hNucleon_Px", "hNucleon_Py", "hNucleon_Pz"};

    for (Int_t hh = 0; hh < (Int_t)InputHistsNames.size(); hh++) {

        for (Int_t ff = 0; ff < (Int_t)InputFilenames.size(); ff++) {

            InputFile = TFile::Open(input_dir + "/" + InputFilenames[ff], "READ");
            ListOfHists = (TList *)InputFile->Get(InputListName);
            ThisHist = (TH1F *)ListOfHists->FindObject(InputHistsNames[hh]);

            c2->cd(hh + 1);

            SetMyHistStyle(ThisHist);
            ThisHist->SetLineWidth(3);
            ThisHist->SetLineColor(Colors[ff].GetNumber());
            ThisHist->Draw("SAME");

            if (hh == 0) {
                TString channel = InputFilenames[ff](10, 1);
                leg->AddEntry(ThisHist, "Channel " + channel, "l");
                ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ff].GetNumber());
            }
        }

        if (hh == 0) leg->Draw();
    }

    c2->Print(output_dir + "/QA_Nucleon.png");
    delete c2;
    delete leg;
}
