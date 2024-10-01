#include "include/Headers.hxx"
#include "include/Style.hxx"

/*
 *
 */
void QA_ChannelD(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputFilenames = {"ChannelD1.73.root", "ChannelD1.8.root", "ChannelD1.87.root", "ChannelD1.94.root", "ChannelD2.01.root"};
    TString ListName = "V0s_Hists";

    /* Case 1: AntiLambda */

    TString V0sName = "AntiLambda";
    std::vector<TString> V0sSet = {"All", "Signal"};
    std::vector<TColor> Colors = {myOrange, myPurple};

    TLegend *leg;
    TFile *InputFile;
    TList *ListOfHists;

    TString HistName;
    TH1F *ThisHist;

    TCanvas *c1 = new TCanvas("c1", "c1", 3 * 800, 2 * 800);

    c1->Divide(3, 2);
    SetMyStyle();
    gStyle->SetOptStat(0);
    c1->SetMargin(0.15, 0.15, 0.15, 0.15);

    for (Int_t ff = 0; ff < (Int_t)InputFilenames.size(); ff++) {

        leg = new TLegend(0.65, 0.7, 0.8, 0.85);
        TString channel = InputFilenames[ff](7, InputFilenames[ff].Length() - 5 - 7);
        leg->SetTextFont(62);
        leg->SetTextSize(0.03);
        leg->SetHeader(channel, "C");
        leg->SetBorderSize(0);

        for (Int_t ss = 0; ss < (Int_t)V0sSet.size(); ss++) {

            InputFile = TFile::Open(input_dir + "/" + InputFilenames[ff], "READ");
            ListOfHists = (TList *)InputFile->Get(ListName);

            HistName = "Found_" + V0sSet[ss] + "_" + V0sName + "_Mass";
            ThisHist = (TH1F *)ListOfHists->FindObject(HistName);

            c1->cd(ff + 1);

            SetMyHistStyle(ThisHist);
            ThisHist->SetLineWidth(3);
            ThisHist->SetLineColor(Colors[ss].GetNumber());
            ThisHist->GetXaxis()->SetLabelSize(0.03);

            ThisHist->Draw("SAME HIST");

            leg->AddEntry(ThisHist, V0sSet[ss] == "All" ? "Signal+Bkg" : V0sSet[ss], "l");
            ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ss].GetNumber());
        }

        leg->Draw();
    }

    system("mkdir -p " + output_dir);
    c1->Print(output_dir + "/ChannelD_" + V0sName + "_Mass.png");

    /* Case 2: Sexaquarks Bookkeep */

    ListName = "Sexaquarks_Hists";
    HistName = "AntiSexaquarks_Bookkeep";
    Colors = {myBlue, myOrange, myGreen, myRed, myPurple};

    TCanvas *c2 = new TCanvas("c2", "c2", 3 * 800, 2 * 800);
    SetMyStyle();
    gStyle->SetOptStat(0);
    c2->SetMargin(0.15, 0.15, 0.15, 0.15);

    leg = new TLegend(0.7, 0.67, 0.8, 0.82);
    leg->SetTextFont(62);
    leg->SetTextSize(0.02);
    leg->SetHeader("Channel D", "C");
    leg->SetBorderSize(0);

    for (Int_t ff = 0; ff < (Int_t)InputFilenames.size(); ff++) {

        InputFile = TFile::Open(input_dir + "/" + InputFilenames[ff], "READ");
        ListOfHists = (TList *)InputFile->Get(ListName);

        ThisHist = (TH1F *)ListOfHists->FindObject(HistName);

        gPad->SetLogy();

        SetMyHistStyle(ThisHist);
        ThisHist->SetLineWidth(3);
        ThisHist->SetLineColor(Colors[ff].GetNumber());
        ThisHist->SetMinimum(1);
        ThisHist->GetXaxis()->SetLabelSize(0);

        ThisHist->Draw("SAME HIST");

        TString mass = InputFilenames[ff](8, InputFilenames[ff].Length() - 5 - 8);
        leg->AddEntry(ThisHist, mass, "l");
        ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ff].GetNumber());
    }

    leg->Draw();

    system("mkdir -p " + output_dir);
    c2->Print(output_dir + "/ChannelD_" + HistName + ".png");

    delete leg;
    c2->Clear();
    delete c2;
}
