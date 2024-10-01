#include "include/Headers.hxx"
#include "include/Style.hxx"

/*
 *
 */
void QA_ChannelH(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputFilenames = {"ChannelH1.73.root", "ChannelH1.8.root", "ChannelH1.87.root", "ChannelH1.94.root", "ChannelH2.01.root"};
    std::vector<TColor> Colors = {myBlue, myOrange, myGreen, myRed, myPurple};

    TString HistName = "PosKaonPairs_Bookkeep";
    TString ListName = "PosKaonPairs_Hists";

    TLegend *leg;
    TFile *InputFile;
    TList *ListOfHists;

    TH1F *ThisHist;

    TCanvas *c1 = new TCanvas("c1", "c1", 3 * 800, 2 * 800);

    SetMyStyle();
    gStyle->SetOptStat(0);
    c1->SetMargin(0.15, 0.15, 0.15, 0.15);

    leg = new TLegend(0.7, 0.67, 0.8, 0.82);
    leg->SetTextFont(62);
    leg->SetTextSize(0.02);
    leg->SetHeader("Channel H", "C");
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
    c1->Print(output_dir + "/ChannelH_" + HistName + ".png");

    delete leg;
    c1->Clear();
    delete c1;
}
