#include "include/Headers.hxx"
#include "include/Style.cxx"

using namespace ROOT;

/*
 *
 */
void SimplePlot_Centrality() {

    /** Input **/

    // TString InputFilename = "../output/local_signalMC_15o_full_kalmanA1.94/AnalysisResults_merged.root";
    TString InputFilename = "../output/local_signalMC_18qr_full_kalmanA1.94/AnalysisResults_merged.root";
    TFile *InputFile = TFile::Open(InputFilename);

    /* Tree */

    RDataFrame RDF_Events("Events", InputFilename);

    auto Hist_Centrality = RDF_Events.Histo1D({"Centrality", ";Centrality (%);N Events", 22, 0., 110.}, "Centrality");

    auto Hist_IsMB = RDF_Events.Filter("IsMB").Histo1D({"Centrality_IsMB", ";;", 22, 0., 110.}, "Centrality");
    auto Hist_IsHighMultV0 = RDF_Events.Filter("IsHighMultV0").Histo1D({"Centrality_IsHighMultV0", ";;", 22, 0., 110.}, "Centrality");
    auto Hist_IsHighMultSPD = RDF_Events.Filter("IsHighMultSPD").Histo1D({"Centrality_IsHighMultSPD", ";;", 22, 0., 110.}, "Centrality");
    auto Hist_IsCentral = RDF_Events.Filter("IsCentral").Histo1D({"Centrality_IsCentral", ";;", 22, 0., 110.}, "Centrality");
    auto Hist_IsSemiCentral = RDF_Events.Filter("IsSemiCentral").Histo1D({"Centrality_IsSemiCentral", ";;", 22, 0., 110.}, "Centrality");

    /* QA Histograms */

    TList *List_QA_Hists = (TList *)InputFile->Get("QA_Hists");
    if (!List_QA_Hists) return;

    TH1F *Hist_QA_Centrality = (TH1F *)List_QA_Hists->FindObject("QA_Centrality");
    TH1F *Hist_QA_Centrality_AfterCuts = (TH1F *)List_QA_Hists->FindObject("QA_Centrality_AfterCuts");

    /** Draw **/

    SetMyStyle();
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 1080, 1080);

    /* Tree */

    Hist_Centrality->GetYaxis()->SetTitleSize(0.04);
    Hist_Centrality->GetYaxis()->SetTitleOffset(1.4);
    Hist_Centrality->GetYaxis()->SetMaxDigits(3);
    Hist_Centrality->GetYaxis()->SetTitleFont(62);

    Hist_Centrality->GetXaxis()->SetTitleSize(0.04);
    Hist_Centrality->GetXaxis()->SetTitleOffset(1.2);
    Hist_Centrality->GetXaxis()->SetTitleFont(62);
    Hist_Centrality->GetXaxis()->SetNdivisions(11);

    Hist_Centrality->SetTitle("");
    Hist_Centrality->SetLineWidth(4);
    Hist_Centrality->SetLineColor(myTransparentBlue.GetNumber());
    Hist_Centrality->SetFillStyle(0);

    Hist_IsMB->SetTitle("");
    Hist_IsMB->SetLineWidth(4);
    Hist_IsMB->SetLineColor(myTransparentOrange.GetNumber());
    Hist_IsMB->SetFillStyle(0);

    Hist_IsHighMultV0->SetTitle("");
    Hist_IsHighMultV0->SetLineWidth(4);
    Hist_IsHighMultV0->SetLineColor(myTransparentGreen.GetNumber());
    Hist_IsHighMultV0->SetFillStyle(0);

    Hist_IsHighMultSPD->SetTitle("");
    Hist_IsHighMultSPD->SetLineWidth(4);
    Hist_IsHighMultSPD->SetLineColor(myTransparentRed.GetNumber());
    Hist_IsHighMultSPD->SetFillStyle(0);

    Hist_IsCentral->SetTitle("");
    Hist_IsCentral->SetLineWidth(4);
    Hist_IsCentral->SetLineColor(myTransparentPurple.GetNumber());
    Hist_IsCentral->SetFillStyle(0);

    Hist_IsSemiCentral->SetTitle("");
    Hist_IsSemiCentral->SetLineWidth(4);
    Hist_IsSemiCentral->SetLineColor(myTransparentPink.GetNumber());
    Hist_IsSemiCentral->SetFillStyle(0);

    Hist_Centrality->Draw();
    Hist_IsMB->Draw("SAME");
    Hist_IsHighMultV0->Draw("SAME");
    Hist_IsHighMultSPD->Draw("SAME");
    Hist_IsCentral->Draw("SAME");
    Hist_IsSemiCentral->Draw("SAME");

    c->Print("gfx/Centrality_fromEventsTree.png");
    c->Clear();

    /* QA Histograms */

    Hist_QA_Centrality->SetLineWidth(4);
    Hist_QA_Centrality->SetLineColor(myTransparentOlive.GetNumber());

    Hist_QA_Centrality_AfterCuts->SetLineWidth(4);
    Hist_QA_Centrality_AfterCuts->SetLineColor(myTransparentCyan.GetNumber());

    Hist_QA_Centrality->Draw();
    Hist_QA_Centrality_AfterCuts->Draw("SAME");

    c->Print("gfx/Centrality_fromQAHists.png");
}
