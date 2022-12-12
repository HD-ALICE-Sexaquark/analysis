#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/Utilities.hxx"

// this macro is a mess...
// -- A. Bórquez

//_____________________________________________________________________________
void PlotTrueV0s(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_TrueV0s_000.root",
    TString output_dir = "gfx/signal+bkg/true_v0s/") {

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  // in case the output dir doesn't exist, create it
  system("mkdir -p " + output_dir);

  // [based on "Event_tt" from "AliAnalysisTaskSexaquark.h"]
  /* define MC variables */                  //
  Int_t N_MCGen;                             // number of MC particles
  std::vector<Float_t> *MC_Px = 0;           // x-component of true momentum
  std::vector<Float_t> *MC_Py = 0;           // y-component of true momentum
  std::vector<Float_t> *MC_Pz = 0;           // z-component of true momentum
  std::vector<Float_t> *MC_X = 0;            // x-coordinate of generation vertex
  std::vector<Float_t> *MC_Y = 0;            // y-coordinate of generation vertex
  std::vector<Float_t> *MC_Z = 0;            // z-coordinate of generation vertex
  std::vector<Int_t> *MC_PID = 0;            // PDG code
  std::vector<Int_t> *MC_Mother = 0;         // index of mother
  std::vector<Int_t> *MC_FirstDau = 0;       // index of first daughter
  std::vector<Int_t> *MC_LastDau = 0;        // index of last daughter
  std::vector<Bool_t> *MC_isSignal = 0;      // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
  /* define track variables  */              //
  std::vector<Int_t> *Idx_True = 0;          // index of true MC particle
  std::vector<Bool_t> *Rec_isDuplicate = 0;  // kTRUE if track is a duplicate, kFALSE if not
  /* define V0 variables */                  //
  Int_t N_V0s;                               // number of formed V0s
  std::vector<Int_t> *Idx_Pos = 0;           // index of positive daughter
  std::vector<Int_t> *Idx_Neg = 0;           // index of negative daughter
  std::vector<Float_t> *V0_Px = 0;           // x-component of V0 momentum
  std::vector<Float_t> *V0_Py = 0;           // y-component of V0 momentum
  std::vector<Float_t> *V0_Pz = 0;           // z-component of V0 momentum
  std::vector<Float_t> *V0_X = 0;            // x-coordinate of V0
  std::vector<Float_t> *V0_Y = 0;            // y-coordinate of V0
  std::vector<Float_t> *V0_Z = 0;            // z-coordinate of V0
  std::vector<Float_t> *Pos_Px = 0;          // x-component of positive track momentum at V0 position
  std::vector<Float_t> *Pos_Py = 0;          // y-component of positive track momentum at V0 position
  std::vector<Float_t> *Pos_Pz = 0;          // z-component of positive track momentum at V0 position
  std::vector<Float_t> *Neg_Px = 0;          // x-component of negative track momentum at V0 position
  std::vector<Float_t> *Neg_Py = 0;          // y-component of negative track momentum at V0 position
  std::vector<Float_t> *Neg_Pz = 0;          // z-component of negative track momentum at V0 position
  std::vector<Bool_t> *V0_isSignal = 0;      // kTRUE if signal, kFALSE if background
  std::vector<Float_t> *V0_E_asK0 = 0;       // energy of V0, assuming it's a K0
  std::vector<Float_t> *V0_E_asAL = 0;       // energy of V0, assuming it's an anti-lambda
  std::vector<Bool_t> *V0_couldBeK0 = 0;     // kTRUE if K0 candidate, kFALSE if not
  std::vector<Bool_t> *V0_couldBeAL = 0;     // kTRUE if anti-lambda candidate, kFALSE if not

  // define histograms
  TH1F *inv_mass_V0 = new TH1F("inv_mass_V0", "inv_mass_V0", 90, 0, 3);
  TH1F *inv_mass_V0_sexa_signal = new TH1F("inv_mass_V0_sexa_signal", "inv_mass_V0_sexa_signal", 90, 0, 3);
  TH1F *radius_K0 = new TH1F("radius_K0", "radius_K0", 90, 0, 180);
  TH1F *radius_K0_sexa_signal = new TH1F("radius_K0_sexa_signal", "radius_K0_sexa_signal", 90, 0, 180);
  TH1F *pt_K0 = new TH1F("pt_K0", "pt_K0", 100, 0., 5.);
  TH1F *pt_K0_sexa_signal = new TH1F("pt_K0_sexa_signal", "pt_K0_sexa_signal", 100, 0., 5.);
  TH1F *radius_AL = new TH1F("radius_AL", "radius_AL", 90, 0, 180);
  TH1F *radius_AL_sexa_signal = new TH1F("radius_AL_sexa_signal", "radius_AL_sexa_signal", 90, 0, 180);
  TH1F *pt_AL = new TH1F("pt_AL", "pt_AL", 100, 0., 5.);
  TH1F *pt_AL_sexa_signal = new TH1F("pt_AL_sexa_signal", "pt_AL_sexa_signal", 100, 0., 5.);

  // (01) inv mass K0 & anti-lambdas without signal cut
  // TH1F *inv_mass_V0_one_dupli = new TH1F("inv_mass_V0_one_dupli", "inv_mass_V0_one_dupli", 150, 0, 3);
  // TH1F *inv_mass_V0_two_dupli = new TH1F("inv_mass_V0_two_dupli", "inv_mass_V0_two_dupli", 150, 0, 3);
  // (02) inv mass K0 & anti-lambdas with signal cut
  // TH1F *inv_mass_V0_sexa_signal_one_dupli = new TH1F("inv_mass_V0_sexa_signal_one_dupli", "inv_mass_V0_sexa_signal_one_dupli", 150, 0,
  // 3); TH1F *inv_mass_V0_sexa_signal_two_dupli = new TH1F("inv_mass_V0_sexa_signal_two_dupli", "inv_mass_V0_sexa_signal_two_dupli", 150,
  // 0, 3); (03) eff. vs R of K0 without signal cut TH1F *n_gen_K0_vs_R = new TH1F("n_gen_K0_vs_R", "n_gen_K0_vs_R", 35, 5, 180); TH1F
  // *eff_K0_vs_R = new TH1F("eff_K0_vs_R", "eff_K0_vs_R", 35, 5, 180); (04) eff. vs R of K0 with signal cut TH1F *n_gen_K0_vs_R_sexa_signal
  // = new TH1F("n_gen_K0_vs_R_sexa_signal", "n_gen_K0_vs_R_sexa_signal", 35, 5, 180); TH1F *eff_K0_vs_R_sexa_signal = new
  // TH1F("eff_K0_vs_R_sexa_signal", "eff_K0_vs_R_sexa_signal", 35, 5, 180); (05) eff. vs Pt of K0 without signal cut TH1F *n_gen_K0_vs_Pt =
  // new TH1F("n_gen_K0_vs_Pt", "n_gen_K0_vs_Pt", 50, 0., 5.); TH1F *eff_K0_vs_Pt = new TH1F("eff_K0_vs_Pt", "eff_K0_vs_Pt", 50, 0., 5.);
  // (06) eff. vs Pt of K0 with signal cut
  // TH1F *n_gen_K0_vs_Pt_sexa_signal = new TH1F("n_gen_K0_vs_Pt_sexa_signal", "n_gen_K0_vs_Pt_sexa_signal", 50, 0., 5.);
  // TH1F *eff_K0_vs_Pt_sexa_signal = new TH1F("eff_K0_vs_Pt_sexa_signal", "eff_K0_vs_Pt_sexa_signal", 50, 0., 5.);
  // (07) eff. vs R of AL without signal cut
  // TH1F *n_gen_AL_vs_R = new TH1F("n_gen_AL_vs_R", "n_gen_AL_vs_R", 35, 5, 180);
  // TH1F *eff_AL_vs_R = new TH1F("eff_AL_vs_R", "eff_AL_vs_R", 35, 5, 180);
  // (08) eff. vs R of AL with signal cut
  // TH1F *n_gen_AL_vs_R_sexa_signal = new TH1F("n_gen_AL_vs_R_sexa_signal", "n_gen_AL_vs_R_sexa_signal", 35, 5, 180);
  // TH1F *eff_AL_vs_R_sexa_signal = new TH1F("eff_AL_vs_R_sexa_signal", "eff_AL_vs_R_sexa_signal", 35, 5, 180);
  // (09) eff. vs Pt of AL without signal cut
  // TH1F *n_gen_AL_vs_Pt = new TH1F("n_gen_AL_vs_Pt", "n_gen_AL_vs_Pt", 50, 0., 5.);
  // TH1F *eff_AL_vs_Pt = new TH1F("eff_AL_vs_Pt", "eff_AL_vs_Pt", 50, 0., 5.);
  // (10) eff. vs Pt of AL with signal cut
  // TH1F *n_gen_AL_vs_Pt_sexa_signal = new TH1F("n_gen_AL_vs_Pt_sexa_signal", "n_gen_AL_vs_Pt_sexa_signal", 50, 0., 5.);
  // TH1F *eff_AL_vs_Pt_sexa_signal = new TH1F("eff_AL_vs_Pt_sexa_signal", "eff_AL_vs_Pt_sexa_signal", 50, 0., 5.);

  // define iterator
  TListIter *list_of_trees_it = new TListIter(list_of_trees);

  // loop over collected trees
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

    // load MC particles branches
    this_tree->SetBranchAddress("N_MCGen", &N_MCGen);
    this_tree->SetBranchAddress("MC_Px", &MC_Px);
    this_tree->SetBranchAddress("MC_Py", &MC_Py);
    this_tree->SetBranchAddress("MC_Pz", &MC_Pz);
    this_tree->SetBranchAddress("MC_X", &MC_X);
    this_tree->SetBranchAddress("MC_Y", &MC_Y);
    this_tree->SetBranchAddress("MC_Z", &MC_Z);
    this_tree->SetBranchAddress("MC_PID", &MC_PID);
    this_tree->SetBranchAddress("MC_Mother", &MC_Mother);
    this_tree->SetBranchAddress("MC_FirstDau", &MC_FirstDau);
    this_tree->SetBranchAddress("MC_LastDau", &MC_LastDau);
    this_tree->SetBranchAddress("MC_isSignal", &MC_isSignal);
    // load track branches
    this_tree->SetBranchAddress("Idx_True", &Idx_True);
    this_tree->SetBranchAddress("Rec_isDuplicate", &Rec_isDuplicate);
    // load V0 branches
    this_tree->SetBranchAddress("N_V0s", &N_V0s);
    this_tree->SetBranchAddress("Idx_Pos", &Idx_Pos);
    this_tree->SetBranchAddress("Idx_Neg", &Idx_Neg);
    this_tree->SetBranchAddress("V0_Px", &V0_Px);
    this_tree->SetBranchAddress("V0_Py", &V0_Py);
    this_tree->SetBranchAddress("V0_Pz", &V0_Pz);
    this_tree->SetBranchAddress("V0_X", &V0_X);
    this_tree->SetBranchAddress("V0_Y", &V0_Y);
    this_tree->SetBranchAddress("V0_Z", &V0_Z);
    this_tree->SetBranchAddress("Pos_Px", &Pos_Px);
    this_tree->SetBranchAddress("Pos_Py", &Pos_Py);
    this_tree->SetBranchAddress("Pos_Pz", &Pos_Pz);
    this_tree->SetBranchAddress("Neg_Px", &Neg_Px);
    this_tree->SetBranchAddress("Neg_Py", &Neg_Py);
    this_tree->SetBranchAddress("Neg_Pz", &Neg_Pz);
    this_tree->SetBranchAddress("V0_isSignal", &V0_isSignal);
    this_tree->SetBranchAddress("V0_E_asK0", &V0_E_asK0);
    this_tree->SetBranchAddress("V0_E_asAL", &V0_E_asAL);
    this_tree->SetBranchAddress("V0_couldBeK0", &V0_couldBeK0);
    this_tree->SetBranchAddress("V0_couldBeAL", &V0_couldBeAL);

    // loop over events
    for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      /*** Count True V0s ***/

      Int_t PID_pos_dau;
      Int_t PID_neg_dau;

      Bool_t exactlyOneTrackIsDupli;  // is exactly one track duplicated?
      Bool_t bothTracksAreDupli;      // are both tracks duplicated?

      Float_t true_mass_asK0;
      Float_t true_mass_asAL;
      Float_t true_R;
      Float_t true_Pt;

      // loop over V0s
      for (Int_t v0 = 0; v0 < N_V0s; v0++) {

        // load variables
        PID_pos_dau = (*MC_PID)[(*Idx_True)[(*Idx_Pos)[v0]]];
        PID_neg_dau = (*MC_PID)[(*Idx_True)[(*Idx_Neg)[v0]]];

        exactlyOneTrackIsDupli = ((*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && !(*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]]) ||
                                 (!(*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]]);
        bothTracksAreDupli = (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]];

        true_mass_asK0 = TMath::Sqrt((*V0_E_asK0)[v0] * (*V0_E_asK0)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                     (*V0_Pz)[v0] * (*V0_Pz)[v0]);
        true_mass_asAL = TMath::Sqrt((*V0_E_asAL)[v0] * (*V0_E_asAL)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                     (*V0_Pz)[v0] * (*V0_Pz)[v0]);
        true_R = TMath::Sqrt((*V0_X)[v0] * (*V0_X)[v0] + (*V0_Y)[v0] * (*V0_Y)[v0]);
        true_Pt = TMath::Sqrt((*V0_Px)[v0] * (*V0_Px)[v0] + (*V0_Py)[v0] * (*V0_Py)[v0]);

        // (cut) on K0 that have relevant daughters, i.e., the K0 -> pi+ pi- channel
        if ((*V0_couldBeK0)[v0] && PID_pos_dau == 211 && PID_neg_dau == -211) {
          inv_mass_V0->Fill(true_mass_asK0);
          pt_K0->Fill(true_Pt);
          radius_K0->Fill(true_R);
          if ((*V0_isSignal)[v0]) {
            inv_mass_V0_sexa_signal->Fill(true_mass_asK0);
            pt_K0_sexa_signal->Fill(true_Pt);
            radius_K0_sexa_signal->Fill(true_R);
          }
        }

        // (cut) on anti-lambdas that have relevant daughters, i.e., the anti-lambda -> anti-proton pi+ channel
        if ((*V0_couldBeAL)[v0] && PID_pos_dau == 211 && PID_neg_dau == -2212) {
          inv_mass_V0->Fill(true_mass_asAL);
          pt_AL->Fill(true_Pt);
          radius_AL->Fill(true_R);
          if ((*V0_isSignal)[v0]) {
            inv_mass_V0_sexa_signal->Fill(true_mass_asAL);
            pt_AL_sexa_signal->Fill(true_Pt);
            radius_AL_sexa_signal->Fill(true_R);
          }
        }

      }  // end of loop over V0s

      /*** Count Generated V0s (K0 & anti-lambdas) ***/
      /*
          Float_t dau_X;
          Float_t dau_Y;

          // loop over MC particles
          for (Int_t mc = 0; mc < N_MCGen; mc++) {

            if ((*MC_PID)[mc] == 310) {
              // (cut) on K0 that have relevant daughters
              // i.e., the K0 -> pi+ pi- channel
              if ((*MC_FirstDau)[mc] > 0 && (*MC_LastDau)[mc] > 0) {
                // (03) eff. vs R of K0 without signal cut
                dau_X = (*MC_X)[(*MC_FirstDau)[mc]];
                dau_Y = (*MC_Y)[(*MC_FirstDau)[mc]];
                n_gen_K0_vs_R->Fill(TMath::Sqrt(dau_X * dau_X + dau_Y * dau_Y));
                // (05) eff. vs Pt of K0 without signal cut
                n_gen_K0_vs_Pt->Fill(TMath::Sqrt((*MC_Px)[mc] * (*MC_Px)[mc] + (*MC_Py)[mc] * (*MC_Py)[mc]));
                // (04) eff. vs R of K0 with signal cut
                // (06) eff. vs Pt of K0 with signal cut
                if ((*MC_isSignal)[mc]) {
                  n_gen_K0_vs_R_sexa_signal->Fill(TMath::Sqrt(dau_X * dau_X + dau_Y * dau_Y));
                  n_gen_K0_vs_Pt_sexa_signal->Fill(TMath::Sqrt((*MC_Px)[mc] * (*MC_Px)[mc] + (*MC_Py)[mc] * (*MC_Py)[mc]));
                }
              }
            }  // end of K0 condition

            if ((*MC_PID)[mc] == -3122) {
              // (cut) on anti-lambdas that have relevant daughters
              // i.e., the anti-lambda -> anti-proton pi+ channel
              if ((*MC_FirstDau)[mc] > 0 && (*MC_LastDau)[mc] > 0) {
                // (07) eff. vs R of AL without signal cut
                dau_X = (*MC_X)[(*MC_FirstDau)[mc]];
                dau_Y = (*MC_Y)[(*MC_FirstDau)[mc]];
                n_gen_AL_vs_R->Fill(TMath::Sqrt(dau_X * dau_X + dau_Y * dau_Y));
                // (09) eff. vs Pt of anti-lambda without signal cut
                n_gen_AL_vs_Pt->Fill(TMath::Sqrt((*MC_Px)[mc] * (*MC_Px)[mc] + (*MC_Py)[mc] * (*MC_Py)[mc]));
                // (08) eff. vs R of AL with signal cut
                // (10) eff. vs Pt of anti-lambda with signal cut
                if ((*MC_isSignal)[mc]) {
                  n_gen_AL_vs_R_sexa_signal->Fill(TMath::Sqrt(dau_X * dau_X + dau_Y * dau_Y));
                  n_gen_AL_vs_Pt_sexa_signal->Fill(TMath::Sqrt((*MC_Px)[mc] * (*MC_Px)[mc] + (*MC_Py)[mc] * (*MC_Py)[mc]));
                }
              }
            }  // end of anti-lambda condition

          }  // end of loop over MC particles
       */
    }  // end of loop over events

  }  // end of loop over trees

  /* Draw */

  // set style to everything
  SetMyStyle();
  gStyle->SetOptStat(0);

  // don't forget to add more histograms when necessary
  std::vector<TH1F *> list_hists = {inv_mass_V0, inv_mass_V0_sexa_signal, pt_K0,     pt_K0_sexa_signal,    pt_AL, pt_AL_sexa_signal,
                                    radius_K0,   radius_K0_sexa_signal,   radius_AL, radius_AL_sexa_signal};

  for (TH1F *hist : list_hists) {
    SetMyHistStyle(hist);
    hist->SetLineWidth(0);
    hist->SetFillStyle(kSolid);
  }

  // the only hists to override this are:
  inv_mass_V0->SetFillStyle(0);
  inv_mass_V0->SetLineWidth(2);
  inv_mass_V0->SetLineColor(kBlack);
  inv_mass_V0->SetMinimum(1);
  inv_mass_V0->GetYaxis()->SetTitle("Counts");
  inv_mass_V0->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-} when K_{S}^{0}, #pi^{+}#bar{p} when #bar{#Lambda}) [GeV/c^{2}]");
  pt_K0->SetFillStyle(0);
  pt_K0->SetLineWidth(2);
  pt_K0->SetLineColor(kBlack);
  pt_K0->SetMinimum(1);
  pt_K0->GetYaxis()->SetTitle("Counts");
  pt_K0->GetXaxis()->SetTitle("p_{T}(K^{0}_{S}) [GeV/c]");
  pt_AL->SetFillStyle(0);
  pt_AL->SetLineWidth(2);
  pt_AL->SetLineColor(kBlack);
  pt_AL->SetMinimum(1);
  pt_AL->GetYaxis()->SetTitle("Counts");
  pt_AL->GetXaxis()->SetTitle("p_{T}(#bar{#Lambda}) [GeV/c]");
  radius_K0->SetFillStyle(0);
  radius_K0->SetLineWidth(2);
  radius_K0->SetLineColor(kBlack);
  radius_K0->SetMinimum(1);
  radius_K0->GetYaxis()->SetTitle("Counts");
  radius_K0->GetXaxis()->SetTitle("Fiducial Radius (K^{0}_{S}) [cm]");
  radius_AL->SetFillStyle(0);
  radius_AL->SetLineWidth(2);
  radius_AL->SetLineColor(kBlack);
  radius_AL->SetMinimum(1);
  radius_AL->GetYaxis()->SetTitle("Counts");
  radius_AL->GetXaxis()->SetTitle("Fiducial Radius (#bar{#Lambda}) [cm]");

  TString output_filename;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  c->SetLogy(1);

  // inv mass V0s
  inv_mass_V0->Draw();
  inv_mass_V0_sexa_signal->SetFillColor(myRed);
  inv_mass_V0_sexa_signal->Draw("SAME");
  inv_mass_V0_sexa_signal->Draw("AXIS SAME");
  TLegend *l = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  // l->SetHeader("#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected ", "C");
  l->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l->AddEntry(inv_mass_V0, "True V0s", "l");
  l->AddEntry(inv_mass_V0_sexa_signal, "#bar{S} signal", "f");
  l->SetTextSize(0.02);
  l->SetBorderSize(0);
  l->SetFillStyle(1001);
  l->Draw();
  output_filename = inv_mass_V0->GetName();
  c->Print(output_dir + output_filename + ".png");

  // pt k0s
  pt_K0->Draw();
  pt_K0_sexa_signal->SetFillColor(myRed);
  pt_K0_sexa_signal->Draw("SAME");
  pt_K0_sexa_signal->Draw("AXIS SAME");
  TLegend *l1 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l1->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l1->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l1->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l1->AddEntry(pt_K0, "True K^{0}_{S}", "l");
  l1->AddEntry(pt_K0_sexa_signal, "#bar{S} signal", "f");
  l1->SetTextSize(0.02);
  l1->SetBorderSize(0);
  l1->SetFillStyle(1001);
  l1->Draw();
  output_filename = pt_K0->GetName();
  c->Print(output_dir + output_filename + ".png");

  // radius k0s
  radius_K0->Draw();
  radius_K0_sexa_signal->SetFillColor(myRed);
  radius_K0_sexa_signal->Draw("SAME");
  radius_K0_sexa_signal->Draw("AXIS SAME");
  TLegend *l2 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l2->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l2->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l2->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l2->AddEntry(radius_K0, "True K^{0}_{S}", "l");
  l2->AddEntry(radius_K0_sexa_signal, "#bar{S} signal", "f");
  l2->SetTextSize(0.02);
  l2->SetBorderSize(0);
  l2->SetFillStyle(1001);
  l2->Draw();
  output_filename = radius_K0->GetName();
  c->Print(output_dir + output_filename + ".png");

  // pt anti-lambda
  pt_AL->Draw();
  pt_AL_sexa_signal->SetFillColor(myRed);
  pt_AL_sexa_signal->Draw("SAME");
  pt_AL_sexa_signal->Draw("AXIS SAME");
  TLegend *l3 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l3->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l3->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l3->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l3->AddEntry(pt_AL, "True #bar{#Lambda}", "l");
  l3->AddEntry(pt_AL_sexa_signal, "#bar{S} signal", "f");
  l3->SetTextSize(0.02);
  l3->SetBorderSize(0);
  l3->SetFillStyle(1001);
  l3->Draw();
  output_filename = pt_AL->GetName();
  c->Print(output_dir + output_filename + ".png");

  // radius anti-lambda
  radius_AL->Draw();
  radius_AL_sexa_signal->SetFillColor(myRed);
  radius_AL_sexa_signal->Draw("SAME");
  radius_AL_sexa_signal->Draw("AXIS SAME");
  TLegend *l4 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l4->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l4->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l4->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l4->AddEntry(radius_AL, "True #bar{#Lambda}", "l");
  l4->AddEntry(radius_AL_sexa_signal, "#bar{S} signal", "f");
  l4->SetTextSize(0.02);
  l4->SetBorderSize(0);
  l4->SetFillStyle(1001);
  l4->Draw();
  output_filename = radius_AL->GetName();
  c->Print(output_dir + output_filename + ".png");

  /*
    // (checking for duplicated tracks...)
    inv_mass_V0->Draw();
    inv_mass_V0_one_dupli->SetFillColor(myBlue);
    inv_mass_V0_one_dupli->Draw("SAME");
    inv_mass_V0_two_dupli->SetFillColor(myGreen);
    inv_mass_V0_two_dupli->Draw("SAME");
    output_filename = (TString)inv_mass_V0->GetName() + "_duplis";
    c->Print(output_dir + output_filename + ".png");

    // (checking for duplicated tracks in signal plots...)
    inv_mass_V0_sexa_signal->SetFillStyle(0);
    inv_mass_V0_sexa_signal->Draw();
    inv_mass_V0_sexa_signal_one_dupli->SetFillColor(myBlue);
    inv_mass_V0_sexa_signal_one_dupli->Draw("SAME");
    inv_mass_V0_sexa_signal_two_dupli->SetFillColor(myGreen);
    inv_mass_V0_sexa_signal_two_dupli->Draw("SAME");
    output_filename = (TString)inv_mass_V0_sexa_signal->GetName() + "_duplis";
    c->Print(output_dir + output_filename + ".png");

    // (03) eff. vs R of K0 without signal cut
    eff_K0_vs_R->Divide(n_rec_K0_vs_R, n_gen_K0_vs_R, 1, 1, "B");
    eff_K0_vs_R->Draw();
    c->Print(output_dir + eff_K0_vs_R->GetName() + ".png");

    // (04) eff. vs R of K0 with signal cut
    eff_K0_vs_R_sexa_signal->Divide(n_rec_K0_vs_R_sexa_signal, n_gen_K0_vs_R_sexa_signal, 1, 1, "B");
    eff_K0_vs_R_sexa_signal->Draw();
    c->Print(output_dir + eff_K0_vs_R_sexa_signal->GetName() + ".png");

    // (05) eff. vs Pt of K0 without signal cut
    eff_K0_vs_Pt->Divide(n_rec_K0_vs_Pt, n_gen_K0_vs_Pt, 1, 1, "B");
    eff_K0_vs_Pt->Draw();
    c->Print(output_dir + eff_K0_vs_Pt->GetName() + ".png");

    // (06) eff. vs Pt of K0 with signal cut
    eff_K0_vs_Pt_sexa_signal->Divide(n_rec_K0_vs_Pt_sexa_signal, n_gen_K0_vs_Pt_sexa_signal, 1, 1, "B");
    eff_K0_vs_Pt_sexa_signal->Draw();
    c->Print(output_dir + eff_K0_vs_Pt_sexa_signal->GetName() + ".png");

    // (07) eff. vs R of AL without signal cut
    eff_AL_vs_R->Divide(n_rec_AL_vs_R, n_gen_AL_vs_R, 1, 1, "B");
    eff_AL_vs_R->Draw();
    c->Print(output_dir + eff_AL_vs_R->GetName() + ".png");

    // (08) eff. vs R of AL with signal cut
    eff_AL_vs_R_sexa_signal->Divide(n_rec_AL_vs_R_sexa_signal, n_gen_AL_vs_R_sexa_signal, 1, 1, "B");
    eff_AL_vs_R_sexa_signal->Draw();
    c->Print(output_dir + eff_AL_vs_R_sexa_signal->GetName() + ".png");

    // (09) eff. vs Pt of AL without signal cut
    eff_AL_vs_Pt->Divide(n_rec_AL_vs_Pt, n_gen_AL_vs_Pt, 1, 1, "B");
    eff_AL_vs_Pt->Draw();
    c->Print(output_dir + eff_AL_vs_Pt->GetName() + ".png");

    // (10) eff. vs Pt of AL with signal cut
    eff_AL_vs_Pt_sexa_signal->Divide(n_rec_AL_vs_Pt_sexa_signal, n_gen_AL_vs_Pt_sexa_signal, 1, 1, "B");
    eff_AL_vs_Pt_sexa_signal->Draw();
    c->Print(output_dir + eff_AL_vs_Pt_sexa_signal->GetName() + ".png");
  */
}
