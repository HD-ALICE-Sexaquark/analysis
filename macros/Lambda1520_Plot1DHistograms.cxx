#include "include/HistUtilities.hxx"
#include "include/Style.hxx"

void Lambda1520_Plot1DHistograms() {

  /* Parameters */

  const Int_t N_HISTS = 6;
  TString InputFiles = "/misc/alidata121/alice_u/borquez/analysis/output/lambda-1520/*/AnalysisResults_*.root";  // uses wildcard

  TString GfxDir = "/misc/alidata121/alice_u/borquez/analysis/gfx/lambda-1520";
  system("mkdir -p " + GfxDir);

  /*** MERGE ***/

  // histogram lists
  //   TList *hList_Bookkeeper = new TList();
  //   AddHistToList(hList_Bookkeeper, InputFiles + "/Bookkeeper");
  TList *hList_TrueLambda_Pt = new TList();
  AddHistsToList(hList_TrueLambda_Pt, InputFiles + "/Hists/TrueLambda_Pt");
  TList *hList_TrueLambda_Pz = new TList();
  AddHistsToList(hList_TrueLambda_Pz, InputFiles + "/Hists/TrueLambda_Pz");
  TList *hList_TrueLambda_Radius = new TList();
  AddHistsToList(hList_TrueLambda_Radius, InputFiles + "/Hists/TrueLambda_Radius");
  TList *hList_TrueAntiLambda_Pt = new TList();
  AddHistsToList(hList_TrueAntiLambda_Pt, InputFiles + "/Hists/TrueAntiLambda_Pt");
  TList *hList_TrueAntiLambda_Pz = new TList();
  AddHistsToList(hList_TrueAntiLambda_Pz, InputFiles + "/Hists/TrueAntiLambda_Pz");
  TList *hList_TrueAntiLambda_Radius = new TList();
  AddHistsToList(hList_TrueAntiLambda_Radius, InputFiles + "/Hists/TrueAntiLambda_Radius");

  // merged histograms
  //   TH1I *mHist_Bookkeeper = new TH1I("Bookkeeper", "Bookkeeper", 7, 0., 7.);
  TH1F *mHist_TrueLambda_Pt = new TH1F("TrueLambda_Pt", "TrueLambda_Pt", 100, 0., 10.);
  TH1F *mHist_TrueLambda_Pz = new TH1F("TrueLambda_Pz", "TrueLambda_Pz", 100, 0., 20.);
  TH1F *mHist_TrueLambda_Radius = new TH1F("TrueLambda_Radius", "TrueLambda_Radius", 100, 0., 10.);
  TH1F *mHist_TrueAntiLambda_Pt = new TH1F("TrueAntiLambda_Pt", "TrueAntiLambda_Pt", 100, 0., 10.);
  TH1F *mHist_TrueAntiLambda_Pz = new TH1F("TrueAntiLambda_Pz", "TrueAntiLambda_Pz", 100, 0., 20.);
  TH1F *mHist_TrueAntiLambda_Radius = new TH1F("TrueAntiLambda_Radius", "TrueAntiLambda_Radius", 100, 0., 10.);

  // group them
  TList *CurrentHistList[N_HISTS] = {hList_TrueLambda_Pt,     hList_TrueLambda_Pz,     hList_TrueLambda_Radius,
                                     hList_TrueAntiLambda_Pt, hList_TrueAntiLambda_Pz, hList_TrueAntiLambda_Radius};
  TH1F *CurrentMergedHist[N_HISTS] = {mHist_TrueLambda_Pt,     mHist_TrueLambda_Pz,     mHist_TrueLambda_Radius,
                                      mHist_TrueAntiLambda_Pt, mHist_TrueAntiLambda_Pz, mHist_TrueAntiLambda_Radius};

  // declare iterator
  TListIter *ListIter;

  for (Int_t i = 0; i < N_HISTS; i++) {
    // then iterate over each list, and add their histograms to the respective final hist
    ListIter = new TListIter(CurrentHistList[i]);

    while (TH1F *hist = (TH1F *)ListIter->Next()) {
      CurrentMergedHist[i]->Add(hist);
    }

    CurrentMergedHist[i]->GetYaxis()->SetRangeUser(0, 1.5 * CurrentMergedHist[i]->GetMaximum());
    CurrentMergedHist[i]->GetYaxis()->SetTitle("Counts");
    CurrentMergedHist[i]->GetYaxis()->SetTitleFont(62);
    // CurrentMergedHist[i]->GetXaxis()->SetTitle("CPA (Gen. #bar{#Lambda})");
    // CurrentMergedHist[i]->GetXaxis()->SetTitleFont(62);
  }

  /*** DRAW ***/

  SetMyStyle();

  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  for (TH1F *CurrentHistogram : CurrentMergedHist) {
    c->cd();

    SetMyHistStyle(CurrentHistogram);
    CurrentHistogram->SetLineColor(myBlue);
    CurrentHistogram->SetLineWidth(2);
    CurrentHistogram->Draw("HIST");

    // set statistics box
    gStyle->SetOptStat("e");
    /*
    // not working... weird...
    TPaveStats *st = (TPaveStats *)CurrentHistogram->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.85);
    st->SetY2NDC(0.9);
    */

    c->Print(GfxDir + "/" + CurrentHistogram->GetName() + ".png");
  }
}
