#include "include/Headers.hxx"
#include "include/Math.hxx"
#include "include/Style.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

// An attempt to do some Quality Assurance directly over True/Official/Custom V0 Finders
// with no additional cut implemented.
// -- A. Bórquez

//_____________________________________________________________________________
void V0s_Plot1DHistograms(
    // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_OfficialV0s_000.root",
    // TString output_dir = "gfx/signal+bkg/official_v0s/") {
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297590/AnalysisResults_CustomV0s_*.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

  /*** Process Input ***/

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  Event_tt this_event;

  /*** Prepare Output ***/

  // define v0 variables
  const Int_t n_v0_vars = 10;
  TString v0_var_name[n_v0_vars] = {"dca(pi+,pi-)", "dca(pi+,anti-proton)",  //
                                    "m(K0)",        "m(AL)",                 //
                                    "y(K0)",        "y(AL)",                 //
                                    "radius(V0)",   "pt(V0)",                //
                                    "dca(V0,PV)",   "cpa(V0,PV)"};           // PENDING: DCA(daughters, V0), armenteros vars
  TString v0_axis_name[n_v0_vars] = {"dca(pi+,pi-)", "dca(pi+,anti-proton)",  //
                                    "m(K0)",        "m(AL)",                 //
                                    "y(K0)",        "y(AL)",                 //
                                    "radius(V0)",   "pt(V0)",                //
                                    "dca(V0,PV)",   "CPA(V0 w.r.t. PV)"};           // PENDING!
  Int_t v0_var_n_bins[n_v0_vars] = {75,  75,                                 //
                                    100, 100,                                //
                                    100, 100,                                //
                                    100, 100,                                //
                                    100, 100};
  Float_t v0_var_min_x[n_v0_vars] = {0.,  0.,   //
                                     0.4, 1.0,  //
                                     -5., -5.,  //
                                     0.,  0.,   //
                                     0.,  0.4};
  Float_t v0_var_max_x[n_v0_vars] = {1.5,  1.5,  //
                                     0.6,  1.2,  //
                                     5.,   5.,   //
                                     200., 10.,  //
                                     100., 1.};

  // create v0 histograms
  TH1F *hist_v0[n_v0_vars];
  TH1F *hist_v0_AL[n_v0_vars];
  TH1F *hist_v0_K0[n_v0_vars];
  TH1F *hist_v0_Sexa[n_v0_vars];
  for (Int_t vv = 0; vv < n_v0_vars; vv++) {
    hist_v0[vv] = new TH1F(v0_var_name[vv], v0_var_name[vv],  //
                           v0_var_n_bins[vv], v0_var_min_x[vv], v0_var_max_x[vv]);
    hist_v0_K0[vv] = new TH1F(v0_var_name[vv] + "_K0", v0_var_name[vv] + "_K0",  //
                              v0_var_n_bins[vv], v0_var_min_x[vv], v0_var_max_x[vv]);
    hist_v0_AL[vv] = new TH1F(v0_var_name[vv] + "_AL", v0_var_name[vv] + "_AL",  //
                              v0_var_n_bins[vv], v0_var_min_x[vv], v0_var_max_x[vv]);
    hist_v0_Sexa[vv] = new TH1F(v0_var_name[vv] + "_Sexa", v0_var_name[vv] + "_Sexa",  //
                                v0_var_n_bins[vv], v0_var_min_x[vv], v0_var_max_x[vv]);
  }

  // define daughter variables
  const Int_t n_dau_vars = 15;
  TString dau_var_name[n_dau_vars] = {"pt(pi+)",        "pt(pi-)",        "pt(anti-p)",         //
                                      "n_sigma(pi+)",   "n_sigma(pi-)",   "n_sigma(anti-p)",    //
                                      "n_cls_tpc(pi+)", "n_cls_tpc(pi-)", "n_cls_tpc(anti-p)",  //
                                      "ip_wrt_pv(pi+)", "ip_wrt_pv(pi-)", "ip_wrt_pv(anti-p)",  //
                                      "eta(pi+)",       "eta(pi-)",       "eta(anti-p)"};       //
  Int_t dau_var_n_bins[n_dau_vars] = {100, 100, 100,                                            //
                                      100, 100, 100,                                            //
                                      100, 100, 100,                                            //
                                      100, 100, 100,                                            //
                                      100, 100, 100};
  Float_t dau_var_min_x[n_dau_vars] = {0.,    0.,    0.,     //
                                       -5.,   -5.,   -5.,    //
                                       0.,    0.,    0.,     //
                                       -200., -200., -200.,  //
                                       -5.,   -5.,   -5.};
  Float_t dau_var_max_x[n_dau_vars] = {10.,  10.,  10.,   //
                                       5.,   5.,   5.,    //
                                       200., 200., 200.,  //
                                       200., 200., 200.,  //
                                       5.,   5.,   5.};

  // create daughter histograms
  TH1F *hist_dau[n_dau_vars];
  TH1F *hist_dau_AL[n_dau_vars];
  TH1F *hist_dau_K0[n_dau_vars];
  TH1F *hist_dau_Sexa[n_dau_vars];
  for (Int_t vv = 0; vv < n_dau_vars; vv++) {
    hist_dau[vv] = new TH1F(dau_var_name[vv], dau_var_name[vv],  //
                            dau_var_n_bins[vv], dau_var_min_x[vv], dau_var_max_x[vv]);
    hist_dau_K0[vv] = new TH1F(dau_var_name[vv] + "_K0", dau_var_name[vv] + "_K0",  //
                               dau_var_n_bins[vv], dau_var_min_x[vv], dau_var_max_x[vv]);
    hist_dau_AL[vv] = new TH1F(dau_var_name[vv] + "_AL", dau_var_name[vv] + "_AL",  //
                               dau_var_n_bins[vv], dau_var_min_x[vv], dau_var_max_x[vv]);
    hist_dau_Sexa[vv] = new TH1F(dau_var_name[vv] + "_Sexa", dau_var_name[vv] + "_Sexa",  //
                                 dau_var_n_bins[vv], dau_var_min_x[vv], dau_var_max_x[vv]);
  }

  // in case the output dir doesn't exist, create it
  system("mkdir -p " + output_dir);

  // loop over collected trees
  TListIter *list_of_trees_it = new TListIter(list_of_trees);
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

    // load branches
    LoadBranches(this_tree, this_event);

    // loop over events
    for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      // loop over V0s
      for (Int_t v0 = 0; v0 < this_event.N_V0s; v0++) {

        // get values
        Bool_t is_K0_signal = kFALSE;
        Bool_t is_AL_signal = kFALSE;
        if ((*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] > 0 &&  //
            (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]] > 0) {
          is_K0_signal = (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]]] == 310 &&  //
                         (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]]] == 310 &&  //
                         (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] ==
                             (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]];
          is_AL_signal = (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]]] == -3122 &&  //
                         (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]]] == -3122 &&  //
                         (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] ==
                             (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]];
        }
        Bool_t is_Sexa_signal = (*this_event.V0_isSignal)[v0];
        Float_t v0_dca_daus = (*this_event.V0_DCA_Daughters)[v0];
        Float_t v0_y_ask0 = 0.5 * TMath::Log(((*this_event.V0_E_asK0)[v0] + (*this_event.V0_Pz)[v0]) /  //
                                             ((*this_event.V0_E_asK0)[v0] - (*this_event.V0_Pz)[v0]));
        Float_t v0_y_asal = 0.5 * TMath::Log(((*this_event.V0_E_asAL)[v0] + (*this_event.V0_Pz)[v0]) /  //
                                             ((*this_event.V0_E_asAL)[v0] - (*this_event.V0_Pz)[v0]));
        Float_t v0_momentum = TMath::Sqrt((*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] +  //
                                          (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0] +  //
                                          (*this_event.V0_Pz)[v0] * (*this_event.V0_Pz)[v0]);
        Float_t v0_mass_ask0 = TMath::Sqrt((*this_event.V0_E_asK0)[v0] * (*this_event.V0_E_asK0)[v0] - v0_momentum * v0_momentum);
        Float_t v0_mass_asal = TMath::Sqrt((*this_event.V0_E_asAL)[v0] * (*this_event.V0_E_asAL)[v0] - v0_momentum * v0_momentum);
        Float_t v0_radius = TMath::Sqrt((*this_event.V0_X)[v0] * (*this_event.V0_X)[v0] + (*this_event.V0_Y)[v0] * (*this_event.V0_Y)[v0]);
        Float_t v0_pt = TMath::Sqrt((*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] + (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0]);
        Float_t v0_dca_v0pv = (*this_event.V0_DCA_wrtPV)[v0];
        Float_t v0_cpa_v0pv = (*this_event.V0_CPA_wrtPV)[v0];

        // "dca(pi+,pi-)" "dca(pi+,anti-proton)"
        // "m(K0)"        "m(AL)"
        // "y(K0)"        "y(AL)"
        // "radius(V0)"   "pt(V0)"
        // "dca(V0,PV)"   "cpa(V0,PV)"

        // (fill hist) BOTH V0s
        if ((*this_event.V0_couldBeK0)[v0] || (*this_event.V0_couldBeAL)[v0]) {
          hist_v0[6]->Fill(v0_radius);
          hist_v0[7]->Fill(v0_pt);
          hist_v0[8]->Fill(v0_dca_v0pv);
          hist_v0[9]->Fill(v0_cpa_v0pv);
          if (is_K0_signal) {
            hist_v0_K0[6]->Fill(v0_radius);
            hist_v0_K0[7]->Fill(v0_pt);
            hist_v0_K0[8]->Fill(v0_dca_v0pv);
            // hist_v0_K0[9]->Fill(v0_cpa_v0pv); // PENDING!
          }
          if (is_AL_signal) {
            hist_v0_AL[6]->Fill(v0_radius);
            hist_v0_AL[7]->Fill(v0_pt);
            hist_v0_AL[8]->Fill(v0_dca_v0pv);
            // hist_v0_AL[9]->Fill(v0_cpa_v0pv); // PENDING!
          }
          if (is_Sexa_signal) {
            hist_v0_Sexa[6]->Fill(v0_radius);
            hist_v0_Sexa[7]->Fill(v0_pt);
            hist_v0_Sexa[8]->Fill(v0_dca_v0pv);
            hist_v0_Sexa[9]->Fill(v0_cpa_v0pv);
          }
        }

        // (fill hist) only K0s
        if ((*this_event.V0_couldBeK0)[v0]) {
          hist_v0[0]->Fill(v0_dca_daus);
          hist_v0[2]->Fill(v0_mass_ask0);
          hist_v0[4]->Fill(v0_y_ask0);
          if (is_K0_signal) {
            hist_v0_K0[0]->Fill(v0_dca_daus);
            hist_v0_K0[2]->Fill(v0_mass_ask0);
            hist_v0_K0[4]->Fill(v0_y_ask0);
          }
          if (is_Sexa_signal) {
            hist_v0_Sexa[0]->Fill(v0_dca_daus);
            hist_v0_Sexa[2]->Fill(v0_mass_ask0);
            hist_v0_Sexa[4]->Fill(v0_y_ask0);
          }
        }

        // (fill hist) only AL
        if ((*this_event.V0_couldBeAL)[v0]) {
          hist_v0[1]->Fill(v0_dca_daus);
          hist_v0[3]->Fill(v0_mass_asal);
          hist_v0[5]->Fill(v0_y_asal);
          if (is_AL_signal) {
            hist_v0_AL[1]->Fill(v0_dca_daus);
            hist_v0_AL[3]->Fill(v0_mass_asal);
            hist_v0_AL[5]->Fill(v0_y_asal);
          }
          if (is_Sexa_signal) {
            hist_v0_Sexa[1]->Fill(v0_dca_daus);
            hist_v0_Sexa[3]->Fill(v0_mass_asal);
            hist_v0_Sexa[5]->Fill(v0_y_asal);
          }
        }

        // loop over daughters of V0
        for (Int_t dau : {(*this_event.Idx_Neg)[v0], (*this_event.Idx_Pos)[v0]}) {

          // get values
          is_K0_signal = kFALSE;
          is_AL_signal = kFALSE;
          if ((*this_event.MC_Mother)[(*this_event.Idx_True)[dau]] > 0) {
            is_K0_signal = (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[dau]]] == 310;
            is_AL_signal = (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[dau]]] == -3122;
          }
          is_Sexa_signal = (*this_event.Rec_isSignal)[dau];
          Float_t dau_pt =
              TMath::Sqrt((*this_event.Rec_Px)[dau] * (*this_event.Rec_Px)[dau] + (*this_event.Rec_Py)[dau] * (*this_event.Rec_Py)[dau]);
          Float_t dau_charge = (*this_event.Rec_Charge)[dau];
          Float_t dau_nsigma_pion = (*this_event.Rec_NSigmaPion)[dau];
          Float_t dau_nsigma_proton = (*this_event.Rec_NSigmaProton)[dau];
          Float_t dau_ip_wrt_pv = (*this_event.Rec_IP_wrtPV)[dau];
          Float_t dau_n_cls_tpc = (*this_event.Rec_NClustersTPC)[dau];
          Float_t dau_theta = 0.5 * TMath::Pi() - TMath::ATan((*this_event.Rec_Pz)[dau]);
          Float_t dau_eta = -TMath::Log(TMath::Tan(0.5 * dau_theta));

          // (fill hist) PI PLUS
          if (TMath::Abs(dau_nsigma_pion) < 3.0 && dau_charge > 0.) {
            hist_dau[0]->Fill(dau_pt);
            hist_dau[3]->Fill(dau_nsigma_pion);
            hist_dau[6]->Fill(dau_n_cls_tpc);
            hist_dau[9]->Fill(dau_ip_wrt_pv);
            hist_dau[12]->Fill(dau_eta);
            if (is_K0_signal) {
              hist_dau_K0[0]->Fill(dau_pt);
              hist_dau_K0[3]->Fill(dau_nsigma_pion);
              hist_dau_K0[6]->Fill(dau_n_cls_tpc);
              hist_dau_K0[9]->Fill(dau_ip_wrt_pv);
              hist_dau_K0[12]->Fill(dau_eta);
            }
            if (is_AL_signal) {
              hist_dau_AL[0]->Fill(dau_pt);
              hist_dau_AL[3]->Fill(dau_nsigma_pion);
              hist_dau_AL[6]->Fill(dau_n_cls_tpc);
              hist_dau_AL[9]->Fill(dau_ip_wrt_pv);
              hist_dau_AL[12]->Fill(dau_eta);
            }
            if (is_Sexa_signal) {
              hist_dau_Sexa[0]->Fill(dau_pt);
              hist_dau_Sexa[3]->Fill(dau_nsigma_pion);
              hist_dau_Sexa[6]->Fill(dau_n_cls_tpc);
              hist_dau_Sexa[9]->Fill(dau_ip_wrt_pv);
              hist_dau_Sexa[12]->Fill(dau_eta);
            }
          }

          // (fill hist) PI MINUS
          if (TMath::Abs(dau_nsigma_pion) < 3.0 && dau_charge < 0.) {
            hist_dau[1]->Fill(dau_pt);
            hist_dau[4]->Fill(dau_nsigma_pion);
            hist_dau[7]->Fill(dau_n_cls_tpc);
            hist_dau[10]->Fill(dau_ip_wrt_pv);
            hist_dau[13]->Fill(dau_eta);
            if (is_K0_signal) {
              hist_dau_K0[1]->Fill(dau_pt);
              hist_dau_K0[4]->Fill(dau_nsigma_pion);
              hist_dau_K0[7]->Fill(dau_n_cls_tpc);
              hist_dau_K0[10]->Fill(dau_ip_wrt_pv);
              hist_dau_K0[13]->Fill(dau_eta);
            }
            if (is_Sexa_signal) {
              hist_dau_Sexa[1]->Fill(dau_pt);
              hist_dau_Sexa[4]->Fill(dau_nsigma_pion);
              hist_dau_Sexa[7]->Fill(dau_n_cls_tpc);
              hist_dau_Sexa[10]->Fill(dau_ip_wrt_pv);
              hist_dau_Sexa[13]->Fill(dau_eta);
            }
          }

          // (fill hist) ANTI-PROTON
          if (TMath::Abs(dau_nsigma_proton) < 3.0 && dau_charge < 0.) {
            hist_dau[2]->Fill(dau_pt);
            hist_dau[5]->Fill(dau_nsigma_proton);
            hist_dau[8]->Fill(dau_n_cls_tpc);
            hist_dau[11]->Fill(dau_ip_wrt_pv);
            hist_dau[14]->Fill(dau_eta);
            if (is_AL_signal) {
              hist_dau_AL[2]->Fill(dau_pt);
              hist_dau_AL[5]->Fill(dau_nsigma_proton);
              hist_dau_AL[8]->Fill(dau_n_cls_tpc);
              hist_dau_AL[11]->Fill(dau_ip_wrt_pv);
              hist_dau_AL[14]->Fill(dau_eta);
            }
            if (is_Sexa_signal) {
              hist_dau_Sexa[2]->Fill(dau_pt);
              hist_dau_Sexa[5]->Fill(dau_nsigma_proton);
              hist_dau_Sexa[8]->Fill(dau_n_cls_tpc);
              hist_dau_Sexa[11]->Fill(dau_ip_wrt_pv);
              hist_dau_Sexa[14]->Fill(dau_eta);
            }
          }

        }  // end of loop over daughters

      }  // end of loop over V0s

    }  // end of loop over events

  }  // end of loop over trees

  /* Draw */

  // set style to everything
  SetMyStyle();
  gStyle->SetOptStat(0);

  TString output_filename;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  c->SetLogy(1);

  // loop over daughter variables
  for (Int_t vv = 0; vv < n_dau_vars; vv++) {

    // "pt(pi+)"      "pt(pi-)"      "pt(anti-p)"
    // "n_sigma(pi+)" "n_sigma(pi-)" "n_sigma(anti-p)"
    // "eta(pi+)"     "eta(pi-)"     "eta(anti-p)"

    SetMyHistStyle(hist_dau[vv]);
    hist_dau[vv]->SetFillStyle(0);
    hist_dau[vv]->SetLineWidth(3);
    hist_dau[vv]->SetLineColor(myBlack);
    hist_dau[vv]->SetMinimum(1);
    hist_dau[vv]->SetMaximum(1000000);  // adjusted for one file
    hist_dau[vv]->GetYaxis()->SetTitle("Counts");
    hist_dau[vv]->GetXaxis()->SetTitle(dau_var_name[vv]);
    hist_dau[vv]->Draw();

    if (vv != 2 && vv != 5 && vv != 8) {
      hist_dau_K0[vv]->SetLineWidth(0);
      hist_dau_K0[vv]->SetFillStyle(kSolid);
      hist_dau_K0[vv]->SetFillColorAlpha(myViolet, 0.5);
      hist_dau_K0[vv]->SetLineWidth(0);
      hist_dau_K0[vv]->Draw("SAME");
    }

    if (vv != 1 && vv != 4 && vv != 7) {
      hist_dau_AL[vv]->SetLineWidth(0);
      hist_dau_AL[vv]->SetFillStyle(kSolid);
      hist_dau_AL[vv]->SetFillColorAlpha(myOrange, 0.5);
      hist_dau_AL[vv]->SetLineWidth(0);
      hist_dau_AL[vv]->Draw("SAME");
    }

    hist_dau_Sexa[vv]->SetLineWidth(0);
    hist_dau_Sexa[vv]->SetFillStyle(kSolid);
    hist_dau_Sexa[vv]->SetFillColorAlpha(myGreen, 0.5);
    hist_dau_Sexa[vv]->SetLineWidth(0);
    hist_dau_Sexa[vv]->Draw("SAME");
    hist_dau_Sexa[vv]->Draw("AXIS SAME");

    // DrawText(Form("N(Offline V0s) = %i", (Int_t)hist_V0[vv]->GetEntries()), myBlack, 0.2, 0.84, kHAlignLeft);
    // DrawText(Form("N(Offline V0s, p > 1.5) = %i", (Int_t)hist_V0_V0Signal[vv]->GetEntries()), myViolet, 0.2, 0.79, kHAlignLeft);
    // DrawText(Form("N(True %s) = %i", particle_opt_latex.Data(), (Int_t)hist_V0_V0Signal[vv]->GetEntries()), myViolet, 0.2, 0.79,
    //          kHAlignLeft);
    // DrawText(Form("N(#bar{S} Signal) = %i", (Int_t)hist_V0_SexaSignal[vv]->GetEntries()), myGreen, 0.2, 0.74, kHAlignLeft);

    output_filename = hist_dau[vv]->GetName();
    c->Print(output_dir + output_filename + ".png");
  }  // end of loop over daughter variables

  // loop over v0 variables
  for (Int_t vv = 0; vv < n_v0_vars; vv++) {

    // "dca(pi+,pi-)" "dca(pi+,anti-proton)"
    // "m(K0)"        "m(AL)"
    // "y(K0)"        "y(AL)"
    // "radius(V0)"   "pt(V0)"
    // "dca(V0,PV)"   "cpa(V0,PV)"

    SetMyHistStyle(hist_v0[vv]);
    hist_v0[vv]->SetFillStyle(0);
    hist_v0[vv]->SetLineWidth(3);
    hist_v0[vv]->SetLineColor(myBlack);
    hist_v0[vv]->SetMinimum(1);
    hist_v0[vv]->SetMaximum(1000000);  // adjusted for one file
    hist_v0[vv]->GetYaxis()->SetTitle("Counts");
    hist_v0[vv]->GetXaxis()->SetTitle(v0_axis_name[vv]);
    hist_v0[vv]->Draw();

    if (vv != 1 && vv != 3 && vv != 5) {
      hist_v0_K0[vv]->SetLineWidth(0);
      hist_v0_K0[vv]->SetFillStyle(kSolid);
      hist_v0_K0[vv]->SetFillColorAlpha(myViolet, 0.5);
      hist_v0_K0[vv]->SetLineWidth(0);
      hist_v0_K0[vv]->Draw("SAME");
    }

    if (vv != 0 && vv != 2 && vv != 4) {
      hist_v0_AL[vv]->SetLineWidth(0);
      hist_v0_AL[vv]->SetFillStyle(kSolid);
      hist_v0_AL[vv]->SetFillColorAlpha(myOrange, 0.5);
      hist_v0_AL[vv]->SetLineWidth(0);
      hist_v0_AL[vv]->Draw("SAME");
    }

    hist_v0_Sexa[vv]->SetLineWidth(0);
    hist_v0_Sexa[vv]->SetFillStyle(kSolid);
    hist_v0_Sexa[vv]->SetFillColorAlpha(myGreen, 0.5);
    hist_v0_Sexa[vv]->SetLineWidth(0);
    hist_v0_Sexa[vv]->Draw("SAME");
    hist_v0_Sexa[vv]->Draw("AXIS SAME");

    // PENDING!
    if (vv == 9) {

      TLegend *l =  new TLegend(0.2, 0.75, 0.6, 0.9, "", "NDC");  // x1,y1,x2,y2
      l->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
      l->AddEntry((TObject *)0, "V0 Candidates", "");
      l->AddEntry(hist_v0[vv], "Signal + Bkg", "l");
      l->AddEntry(hist_v0_Sexa[vv], "Signal", "f");
      l->SetTextSize(0.03);
      l->SetBorderSize(0);
      l->SetFillStyle(1001);
      l->Draw();
    }

    output_filename = hist_v0[vv]->GetName();
    c->Print(output_dir + output_filename + ".png");
  }  // end of loop over v0 variables
}
