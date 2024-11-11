#include "include/Headers.hxx"
#include "include/Style.cxx"

/*
 * QA, Interaction Radius
 */
void QA_InteractionRadius(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputFilenames = {"ChannelA.root", "ChannelD.root", "ChannelE.root", "ChannelH.root"};
    TString InputListName = "Sexaquarks_Hists";
    TString InputHistName = "MCGen_All_AntiSexaquark_Radius";

    /* Prepare objects */

    TFile *InputFile;
    TList *ListOfHists;
    TH1D *ThisHist;

    std::vector<TColor> Colors = {myBlue, myOrange, myGreen, myRed};

    /* Case: Interaction Radius */

    SetMyStyle();
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 800, 800);

    TLegend *leg = new TLegend(0.2, 0.2, 0.4, 0.35);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    for (Int_t ff = 0; ff < (Int_t)InputFilenames.size(); ff++) {

        InputFile = new TFile(input_dir + "/" + InputFilenames[ff], "READ");
        ListOfHists = (TList *)InputFile->Get(InputListName);
        ThisHist = (TH1D *)ListOfHists->FindObject(InputHistName);

        SetMyHistStyle(ThisHist);
        ThisHist->SetLineWidth(3);
        ThisHist->SetLineColor(Colors[ff].GetNumber());
        ThisHist->Draw("SAME");

        TString channel = InputFilenames[ff](7, 1);
        leg->AddEntry(ThisHist, "Channel " + channel, "l");
        ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ff].GetNumber());

        InputFile->Close();
        delete InputFile;
    }

    leg->Draw();

    system("mkdir -p " + output_dir);
    c->Print(output_dir + "/QA_InteractionRadius.png");

    delete c;
}
