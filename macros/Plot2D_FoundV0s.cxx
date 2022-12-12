#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/Utilities.hxx"

// -- A. Bórquez

//_____________________________________________________________________________
void Plot2D_FoundV0s(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/AnalysisResults_OfficialV0s_*.root",
    TString output_dir = "gfx/signal+bkg/official_v0s/") {

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  // in case the output dir doesn't exist, create it
  system("mkdir -p " + output_dir);

  // [based on "Event_tt" from "AliAnalysisTaskSexaquark.h"]
  /* define V0 variables */                    //
  Int_t N_V0s;                                 // number of formed V0s
  std::vector<Int_t> *Idx_Pos = 0;             // index of positive daughter
  std::vector<Int_t> *Idx_Neg = 0;             // index of negative daughter
  std::vector<Float_t> *V0_Px = 0;             // x-component of V0 momentum
  std::vector<Float_t> *V0_Py = 0;             // y-component of V0 momentum
  std::vector<Float_t> *V0_Pz = 0;             // z-component of V0 momentum
  std::vector<Float_t> *V0_X = 0;              // x-coordinate of V0
  std::vector<Float_t> *V0_Y = 0;              // y-coordinate of V0
  std::vector<Float_t> *V0_Z = 0;              // z-coordinate of V0
  std::vector<Float_t> *Pos_Px = 0;            // x-component of positive track momentum at V0 position
  std::vector<Float_t> *Pos_Py = 0;            // y-component of positive track momentum at V0 position
  std::vector<Float_t> *Pos_Pz = 0;            // z-component of positive track momentum at V0 position
  std::vector<Float_t> *Neg_Px = 0;            // x-component of negative track momentum at V0 position
  std::vector<Float_t> *Neg_Py = 0;            // y-component of negative track momentum at V0 position
  std::vector<Float_t> *Neg_Pz = 0;            // z-component of negative track momentum at V0 position
  std::vector<Bool_t> *V0_isSignal = 0;        // kTRUE if signal, kFALSE if background
  std::vector<Float_t> *V0_E_asK0 = 0;         // energy of V0, assuming it's a K0
  std::vector<Float_t> *V0_E_asAL = 0;         // energy of V0, assuming it's an anti-lambda
  std::vector<Bool_t> *V0_couldBeK0 = 0;       // kTRUE if K0 candidate, kFALSE if not
  std::vector<Bool_t> *V0_couldBeAL = 0;       // kTRUE if anti-lambda candidate, kFALSE if not
  std::vector<Bool_t> *V0_onFlyStatus = 0;     // kTRUE if comes from the on-the-fly finder, kFALSE if not
  std::vector<Float_t> *V0_Chi2 = 0;           // chi2 value from the Kalman Filter (?)
  std::vector<Float_t> *V0_DCA_Daughters = 0;  // distance of closest approach between daughters
  std::vector<Float_t> *V0_IP_wrtPV = 0;       // impact parameter w.r.t. Primary Vertex
  std::vector<Float_t> *V0_CPA_wrtPV = 0;      // cosine of pointing angle w.r.t. Primary Vertex
  std::vector<Float_t> *V0_ArmAlpha = 0;       // Armenteros-Podolanski variable alpha
  std::vector<Float_t> *V0_ArmPt = 0;          // Armenteros-Podolanski variable Pt
  std::vector<Float_t> *V0_DecayLength = 0;    // distance between PV and V0

  // define histograms
  // (01) inv mass K0
  TH2F *inv_mass_K0_vs_Chi2 = new TH2F("inv_mass_K0_vs_Chi2", "inv_mass_K0_vs_Chi2", 132, 0, 33, 150, 0, 3);
  TH2F *inv_mass_K0_vs_DCA_Daughters = new TH2F("inv_mass_K0_vs_DCA_Daughters", "inv_mass_K0_vs_DCA_Daughters", 150, 0, 1.5, 150, 0, 3);
  TH2F *inv_mass_K0_vs_IP_wrtPV = new TH2F("inv_mass_K0_vs_IP_wrtPV", "inv_mass_K0_vs_IP_wrtPV", 150, 0., 1.5, 150, 0, 3);
  TH2F *inv_mass_K0_vs_CPA_wrtPV = new TH2F("inv_mass_K0_vs_CPA_wrtPV", "inv_mass_K0_vs_CPA_wrtPV", 200, 0, 2, 150, 0, 3);
  TH2F *inv_mass_K0_vs_DecayLength = new TH2F("inv_mass_K0_vs_DecayLength", "inv_mass_K0_vs_DecayLength", 150, 0, 50, 150, 0, 3);
  TH2F *inv_mass_K0_vs_Pt = new TH2F("inv_mass_K0_vs_Pt", "inv_mass_K0_vs_Pt", 150, 0, 5, 150, 0, 3);

  // (02) inv mass anti-lambdas
  TH2F *inv_mass_AL_vs_Chi2 = new TH2F("inv_mass_AL_vs_Chi2", "inv_mass_AL_vs_Chi2", 132, 0, 33, 150, 0, 3);
  TH2F *inv_mass_AL_vs_DCA_Daughters = new TH2F("inv_mass_AL_vs_DCA_Daughters", "inv_mass_AL_vs_DCA_Daughters", 150, 0, 1.5, 150, 0, 3);
  TH2F *inv_mass_AL_vs_IP_wrtPV = new TH2F("inv_mass_AL_vs_IP_wrtPV", "inv_mass_AL_vs_IP_wrtPV", 150, 0., 1.5, 150, 0, 3);
  TH2F *inv_mass_AL_vs_CPA_wrtPV = new TH2F("inv_mass_AL_vs_CPA_wrtPV", "inv_mass_AL_vs_CPA_wrtPV", 200, 0, 2, 150, 0, 3);
  TH2F *inv_mass_AL_vs_DecayLength = new TH2F("inv_mass_AL_vs_DecayLength", "inv_mass_AL_vs_DecayLength", 150, 0, 50, 150, 0, 3);
  TH2F *inv_mass_AL_vs_Pt = new TH2F("inv_mass_AL_vs_Pt", "inv_mass_AL_vs_Pt", 150, 0, 5, 150, 0, 3);

  // (armenteros-podolanski plot)
  TH2F *ArmPt_vs_ArmAlpha = new TH2F("armenteros_podolanski", "armenteros_podolanski", 200, -1., 1., 150, 0, 0.3);
  TH2F *ArmPt_vs_ArmAlpha_K0 = new TH2F("armenteros_podolanski_K0", "armenteros_podolanski_K0", 200, -1., 1., 150, 0, 0.3);
  TH2F *ArmPt_vs_ArmAlpha_AL = new TH2F("armenteros_podolanski_AL", "armenteros_podolanski_AL", 200, -1., 1., 150, 0, 0.3);

  // define iterator
  TListIter *list_of_trees_it = new TListIter(list_of_trees);

  // loop over collected trees
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

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
    this_tree->SetBranchAddress("V0_onFlyStatus", &V0_onFlyStatus);
    this_tree->SetBranchAddress("V0_Chi2", &V0_Chi2);
    this_tree->SetBranchAddress("V0_DCA_Daughters", &V0_DCA_Daughters);
    this_tree->SetBranchAddress("V0_IP_wrtPV", &V0_IP_wrtPV);
    this_tree->SetBranchAddress("V0_CPA_wrtPV", &V0_CPA_wrtPV);
    this_tree->SetBranchAddress("V0_ArmAlpha", &V0_ArmAlpha);
    this_tree->SetBranchAddress("V0_ArmPt", &V0_ArmPt);
    this_tree->SetBranchAddress("V0_DecayLength", &V0_DecayLength);

    // loop over events
    for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      /*** Count Found V0s ***/

      Float_t V0_Pt;

      Float_t mass_asK0;
      Float_t mass_asAL;

      // loop over V0s
      for (Int_t v0 = 0; v0 < N_V0s; v0++) {

        // (cut)
        // reject V0s found online/on-the-fly
        if ((*V0_onFlyStatus)[v0]) {
          continue;
        }

        V0_Pt = TMath::Sqrt((*V0_Px)[v0] * (*V0_Px)[v0] + (*V0_Py)[v0] * (*V0_Py)[v0]);
        mass_asK0 = TMath::Sqrt((*V0_E_asK0)[v0] * (*V0_E_asK0)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                (*V0_Pz)[v0] * (*V0_Pz)[v0]);
        mass_asAL = TMath::Sqrt((*V0_E_asAL)[v0] * (*V0_E_asAL)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                (*V0_Pz)[v0] * (*V0_Pz)[v0]);

        if ((*V0_couldBeK0)[v0]) {
          inv_mass_K0_vs_Chi2->Fill((*V0_Chi2)[v0], mass_asK0);
          inv_mass_K0_vs_DCA_Daughters->Fill((*V0_DCA_Daughters)[v0], mass_asK0);
          inv_mass_K0_vs_IP_wrtPV->Fill((*V0_IP_wrtPV)[v0], mass_asK0);
          inv_mass_K0_vs_CPA_wrtPV->Fill((*V0_CPA_wrtPV)[v0], mass_asK0);
          inv_mass_K0_vs_DecayLength->Fill((*V0_DecayLength)[v0], mass_asK0);
          inv_mass_K0_vs_Pt->Fill(V0_Pt, mass_asK0);
          ArmPt_vs_ArmAlpha_K0->Fill((*V0_ArmAlpha)[v0], (*V0_ArmPt)[v0]);
        }

        if ((*V0_couldBeAL)[v0]) {
          inv_mass_AL_vs_Chi2->Fill((*V0_Chi2)[v0], mass_asAL);
          inv_mass_AL_vs_DCA_Daughters->Fill((*V0_DCA_Daughters)[v0], mass_asAL);
          inv_mass_AL_vs_IP_wrtPV->Fill((*V0_IP_wrtPV)[v0], mass_asAL);
          inv_mass_AL_vs_CPA_wrtPV->Fill((*V0_CPA_wrtPV)[v0], mass_asAL);
          inv_mass_AL_vs_DecayLength->Fill((*V0_DecayLength)[v0], mass_asAL);
          inv_mass_AL_vs_Pt->Fill(V0_Pt, mass_asAL);
          ArmPt_vs_ArmAlpha_AL->Fill((*V0_ArmAlpha)[v0], (*V0_ArmPt)[v0]);
        }

        // (armenteros-podolanski plot)
        ArmPt_vs_ArmAlpha->Fill((*V0_ArmAlpha)[v0], (*V0_ArmPt)[v0]);
      }  // end of loop over V0s

    }  // end of loop over events

  }  // end of loop over trees

  /* Draw */

  // set style to everything
  SetMy2DStyle();

  // don't forget to add more histograms when necessary
  std::vector<TH2F *> list_hists = {inv_mass_K0_vs_Chi2,      inv_mass_K0_vs_DCA_Daughters, inv_mass_K0_vs_IP_wrtPV,
                                    inv_mass_K0_vs_CPA_wrtPV, inv_mass_K0_vs_DecayLength,   inv_mass_K0_vs_Pt,
                                    inv_mass_AL_vs_Chi2,      inv_mass_AL_vs_DCA_Daughters, inv_mass_AL_vs_IP_wrtPV,
                                    inv_mass_AL_vs_CPA_wrtPV, inv_mass_AL_vs_DecayLength,   inv_mass_AL_vs_Pt,
                                    ArmPt_vs_ArmAlpha,        ArmPt_vs_ArmAlpha_K0,         ArmPt_vs_ArmAlpha_AL};

  for (TH2F *hist : list_hists) {
    SetMy2DHistStyle(hist);
  }

  TString output_filename;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);
  c->SetFrameLineWidth(2);

  for (TH2F *hist : list_hists) {
    if (((TString)hist->GetName()).CompareTo("inv_mass_K0_vs_DecayLength") == 0 ||
        ((TString)hist->GetName()).CompareTo("inv_mass_AL_vs_DecayLength") == 0) {
      c->SetLogz(1);
    }
    hist->Draw("COLZ");
    output_filename = hist->GetName();
    c->Print(output_dir + output_filename + ".png");
    if (((TString)hist->GetName()).CompareTo("inv_mass_K0_vs_DecayLength") == 0 ||
        ((TString)hist->GetName()).CompareTo("inv_mass_AL_vs_DecayLength") == 0) {
      c->SetLogz(0);
    }
  }
}
