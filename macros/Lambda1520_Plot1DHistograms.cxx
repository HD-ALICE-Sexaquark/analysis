void Lambda1520_Plot1DHistograms() {

  // define locations
  TString InputDir = "/misc/alidata121/alice_u/borquez/analysis/output/";
  TString InputFiles = "AnalysisResults_" + SimulationSet + "_run*.root";  // uses wildcard

  TString OutputDir = "/misc/alidata121/alice_u/borquez/analysis/output_merged/";
  TString OutputFilename = "AnalysisResults_" + SimulationSet + "_merged.root";

  TString GfxDir = "/misc/alidata121/alice_u/borquez/analysis/gfx/";

  /*** MERGE ***/

  // sim. rec.
  TList *HistList_Mass = new TList();
  AddHistToList(HistList_Mass, InputDir + InputFiles + "/fHist_Sexaquark_M");

  TList *HistList_Momentum = new TList();
  AddHistToList(HistList_Momentum, InputDir + InputFiles + "/fHist_Sexaquark_P");

  TList *HistList_Pt = new TList();
  AddHistToList(HistList_Pt, InputDir + InputFiles + "/fHist_Sexaquark_Pt");

  // sim. gen.
  TList *HistList_MC_Mass = new TList();
  AddHistToList(HistList_MC_Mass, InputDir + InputFiles + "/fHistMC_Sexaquark_M");

  TList *HistList_MC_Momentum = new TList();
  AddHistToList(HistList_MC_Momentum, InputDir + InputFiles + "/fHistMC_Sexaquark_P");

  TList *HistList_MC_Pt = new TList();
  AddHistToList(HistList_MC_Pt, InputDir + InputFiles + "/fHistMC_Sexaquark_Pt");

  // experimental hists
  TList *HistList_MC_AntiL_CPA = new TList();
  AddHistToList(HistList_MC_AntiL_CPA, InputDir + InputFiles + "/fHistMC_AntiL_CPA");

  TList *HistList_MC_K0_CPA = new TList();
  AddHistToList(HistList_MC_K0_CPA, InputDir + InputFiles + "/fHistMC_K0_CPA");

  // some options
  Double_t PDG_M = 1.8;
  Double_t MIN_CPA = 0.9;
  const Int_t N_HISTS = 8;

  // define merged/final histograms
  // sim. rec.
  TH1F *MergedHist_Mass = new TH1F("AntiS_M", "AntiS_M", 100, 0., 5.);
  TH1F *MergedHist_Momentum = new TH1F("AntiS_P", "AntiS_P", 8, 0., 8.);
  TH1F *MergedHist_Pt = new TH1F("AntiS_Pt", "AntiS_Pt", 5, 0., 5.);
  // sim. gen.
  TH1F *MergedHist_MC_Mass = new TH1F("MC_AntiS_M", "MC_AntiS_M", 100, 0., 5.);
  TH1F *MergedHist_MC_Momentum = new TH1F("MC_AntiS_P", "MC_AntiS_P", 8, 0., 8.);
  TH1F *MergedHist_MC_Pt = new TH1F("MC_AntiS_Pt", "MC_AntiS_Pt", 5, 0., 5.);
  // experimental hists
  TH1F *MergedHist_MC_AntiL_CPA = new TH1F("MC_AntiL_CPA", "MC_AntiL_CPA", 200, -1., 1.);
  TH1F *MergedHist_MC_K0_CPA = new TH1F("MC_K0_CPA", "MC_K0_CPA", 200, -1., 1.);

  // group them
  TList *CurrentHistList[N_HISTS] = {HistList_Mass,        HistList_Momentum, HistList_Pt,           HistList_MC_Mass,
                                     HistList_MC_Momentum, HistList_MC_Pt,    HistList_MC_AntiL_CPA, HistList_MC_K0_CPA};
  TH1F *CurrentMergedHist[N_HISTS] = {MergedHist_Mass,        MergedHist_Momentum, MergedHist_Pt,           MergedHist_MC_Mass,
                                      MergedHist_MC_Momentum, MergedHist_MC_Pt,    MergedHist_MC_AntiL_CPA, MergedHist_MC_K0_CPA};

  // declare iterator
  TListIter *ListIter;

  // create output file
  TFile *OutputFile = new TFile(OutputDir + OutputFilename, "RECREATE");

  for (Int_t i = 0; i < N_HISTS; i++) {
    // then iterate over each list, and add their histograms to the respective final hist
    ListIter = new TListIter(CurrentHistList[i]);

    while (TH1F *hist = (TH1F *)ListIter->Next()) {
      CurrentMergedHist[i]->Add(hist);
    }

    CurrentMergedHist[i]->GetYaxis()->SetRangeUser(0, 1.5 * CurrentMergedHist[i]->GetMaximum());
    CurrentMergedHist[i]->GetYaxis()->SetTitle("Counts");
    CurrentMergedHist[i]->GetYaxis()->SetTitleFont(62);
    if (((TString)CurrentMergedHist[i]->GetName()).Contains("_M")) {
      CurrentMergedHist[i]->GetXaxis()->SetTitle("#it{M}_{" + SimulationLatex + "} (GeV/c^{2})");
    } else if (((TString)CurrentMergedHist[i]->GetName()).Contains("_Pt")) {
      CurrentMergedHist[i]->GetXaxis()->SetTitle("#it{p}_{T} (#bar{S}) (GeV/c)");
    } else if (((TString)CurrentMergedHist[i]->GetName()).Contains("_P")) {
      CurrentMergedHist[i]->GetXaxis()->SetTitle("#it{p} (#bar{S}) (GeV/c)");
    } else if (((TString)CurrentMergedHist[i]->GetName()).Contains("_K0_CPA")) {
      CurrentMergedHist[i]->GetXaxis()->SetTitle("CPA (Gen. K^{0}_{S})");
    } else if (((TString)CurrentMergedHist[i]->GetName()).Contains("_AntiL_CPA")) {
      CurrentMergedHist[i]->GetXaxis()->SetTitle("CPA (Gen. #bar{#Lambda})");
    }
    CurrentMergedHist[i]->GetXaxis()->SetTitleFont(62);

    // save onto output file
    CurrentMergedHist[i]->Write();
  }

  /*** DRAW ***/

  SetMyStyle();

  TString CanvasName;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  for (TH1F *CurrentHistogram : CurrentMergedHist) {
    c->cd();

    CanvasName = SimulationSet + "_" + (TString)CurrentHistogram->GetName();

    SetMyHistStyle(CurrentHistogram);
    CurrentHistogram->SetLineColor(myBlue);
    CurrentHistogram->SetLineWidth(2);
    CurrentHistogram->Draw("HIST");

    // set statistics box
    gStyle->SetOptStat("e");
    TPaveStats *st = (TPaveStats *)CurrentHistogram->FindObject("stats");
    st->SetX1NDC(0.7);
    st->SetX2NDC(0.9);
    st->SetY1NDC(0.85);
    st->SetY2NDC(0.9);

    // draw vertical lines
    if (((TString)CurrentHistogram->GetName()).Contains("_M")) {
      DrawVerticalLine(PDG_M, myRed, kDashed, 2);
    } else if (((TString)CurrentHistogram->GetName()).Contains("_CPA")) {
      DrawVerticalLine(MIN_CPA, myRed, kDashed, 2);
    }

    c->Print(GfxDir + CanvasName + ".png");
  }

  // close output file
  OutputFile->Close();

  // print output file path
  std::cout << "The following file has been created: " << OutputFilename << std::endl;
}
