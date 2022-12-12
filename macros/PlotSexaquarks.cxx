#include "include/Headers.hxx"
#include "include/Style.hxx"

// this macro is a mess...
// -- A. Bórquez

//_____________________________________________________________________________
void PlotSexaquarks(  //
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/SexaquarkResults_TrueV0s_*.root",
    TString output_dir = "gfx/signal+bkg/true_v0s/", Bool_t useLogScale = kFALSE) {
  // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/SexaquarkResults_OfficialV0s_000.root",
  // TString output_dir = "gfx/signal+bkg/official_v0s/", Bool_t useLogScale = kTRUE) {

  // single file
  // TFile *input_file = new TFile(input_filename);
  // TTree *input_tree = (TTree *)input_file->Get("Sexaquarks");

  // wildcard!
  TChain *input_chain = new TChain("hola");
  input_chain->Add(input_filename + "/Sexaquarks");

  // define input sexaquark variables
  Bool_t Sexa_isSignal;
  Float_t Sexa_M;
  input_chain->SetBranchAddress("isSignal", &Sexa_isSignal);
  input_chain->SetBranchAddress("M", &Sexa_M);

  // in case the output dir doesn't exist, create it
  system("mkdir -p " + output_dir);

  // define histograms
  TH1F *inv_mass_sexa = new TH1F("inv_mass_sexa", "inv_mass_sexa", 150, 0, 3);
  TH1F *inv_mass_sexa_signal = new TH1F("inv_mass_sexa_signal", "inv_mass_sexa_signal", 150, 0, 3);

  // loop over candidates
  for (Int_t candidate = 0; candidate < input_chain->GetEntries(); candidate++) {

    input_chain->GetEntry(candidate);

    // fill hists
    inv_mass_sexa->Fill(Sexa_M);
    if (Sexa_isSignal) {
      std::cout << "PlotSexaquarks :: A signal anti-sexaquark was found!!" << std::endl;
      inv_mass_sexa_signal->Fill(Sexa_M);
    }
  }  // end of loop over candidates

  /* Draw */

  // set style to everything
  SetMyStyle();
  gStyle->SetOptStat(0);

  // don't forget to add more histograms when necessary
  std::vector<TH1F *> list_hists = {inv_mass_sexa, inv_mass_sexa_signal};

  for (TH1F *hist : list_hists) {
    SetMyHistStyle(hist);
    hist->SetLineWidth(0);
    hist->SetFillStyle(kSolid);
  }

  // the only hists to override this are:
  inv_mass_sexa->SetFillStyle(0);
  inv_mass_sexa->SetLineWidth(2);
  inv_mass_sexa->SetLineColor(kBlack);
  inv_mass_sexa->SetMinimum(1);
  inv_mass_sexa->GetYaxis()->SetTitle("Counts");
  inv_mass_sexa->GetXaxis()->SetTitle("Reconstructed Inv. Mass [GeV/c^{2}]");

  TString output_filename;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  // inv mass sexaquarks
  c->SetLogy((Int_t)useLogScale);
  inv_mass_sexa->Draw();
  inv_mass_sexa_signal->SetFillColor(myRed);
  inv_mass_sexa_signal->Draw("SAME");
  inv_mass_sexa_signal->Draw("AXIS SAME");
  inv_mass_sexa_signal->SetMinimum((Int_t)useLogScale);

  TLegend *l = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  // l->SetHeader("#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected ", "C");
  l->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l->AddEntry(inv_mass_sexa, "#bar{S} candidates", "l");
  l->AddEntry(inv_mass_sexa_signal, "#bar{S} signal", "f");
  l->SetTextSize(0.025);
  l->SetBorderSize(0);
  l->SetFillStyle(1001);
  l->Draw();

  output_filename = inv_mass_sexa->GetName();
  c->Print(output_dir + output_filename + ".png");
}
