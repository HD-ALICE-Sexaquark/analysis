#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/Utilities.hxx"

// this macro is also a mess...
// -- A. Bórquez

//_____________________________________________________________________________
void PlotFoundV0s(  //
                    // TString input_filename =
                    // "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_CustomV0s_000.root",
                    //                 TString output_dir = "gfx/signal+bkg/custom_v0s/") {
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595_v1/AnalysisResults_CustomV0s_000.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  // in case the output dir doesn't exist, create it
  system("mkdir -p " + output_dir);

  // [based on "Event_tt" from "AliAnalysisTaskSexaquark.h"]
  /* define MC variables */                  //
  std::vector<Int_t> *MC_PID = 0;            // PDG code
  std::vector<Int_t> *MC_Mother = 0;         // index of mother
  std::vector<Bool_t> *MC_isSignal = 0;      // index of mother
  /* define track variables  */              //
  std::vector<Int_t> *Idx_True = 0;          // index of true MC particle
  std::vector<Float_t> *Rec_Px = 0;          //
  std::vector<Float_t> *Rec_Py = 0;          //
  std::vector<Float_t> *Rec_Pz = 0;          //
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
  // (01) inv mass K0
  TH1F *inv_mass_K0 = new TH1F("inv_mass_K0", "inv_mass_K0", 150, 0, 3);
  // (duplicated tracks)
  TH1F *inv_mass_K0_one_dupli = new TH1F("inv_mass_K0_one_dupli", "inv_mass_K0_one_dupli", 150, 0, 3);
  TH1F *inv_mass_K0_two_dupli = new TH1F("inv_mass_K0_two_dupli", "inv_mass_K0_two_dupli", 150, 0, 3);
  TH1F *inv_mass_K0_no_dupli = new TH1F("inv_mass_K0_no_dupli", "inv_mass_K0_no_dupli", 150, 0, 3);
  // (ambiguous tracks)
  TH1F *inv_mass_K0_also_AL = new TH1F("inv_mass_K0_also_AL", "inv_mass_K0_also_AL", 150, 0, 3);
  TH1F *inv_mass_K0_no_ambi = new TH1F("inv_mass_K0_no_ambi", "inv_mass_K0_no_ambi", 150, 0, 3);
  // (both)
  TH1F *inv_mass_K0_no_dupli_no_ambi = new TH1F("inv_mass_K0_no_dupli_no_ambi", "inv_mass_K0_no_dupli_no_ambi", 150, 0, 3);
  // (K0 signal)
  TH1F *inv_mass_K0_signal = new TH1F("inv_mass_K0_signal", "inv_mass_K0_signal", 150, 0, 3);
  TH1F *alt_mass_K0_signal = new TH1F("alt_mass_K0_signal", "alt_mass_K0_signal", 150, 0, 3);
  // (02) inv mass K0 with sexaquark signal cut
  TH1F *inv_mass_K0_sexa_signal = new TH1F("inv_mass_K0_sexa_signal", "inv_mass_K0_sexa_signal", 150, 0, 3);
  TH1F *alt_mass_K0_sexa_signal = new TH1F("alt_mass_K0_sexa_signal", "alt_mass_K0_sexa_signal", 150, 0, 3);

  // (03) inv mass anti-lambdas
  TH1F *inv_mass_AL = new TH1F("inv_mass_AL", "inv_mass_AL", 150, 0, 5);
  // (duplicated tracks)
  TH1F *inv_mass_AL_one_dupli = new TH1F("inv_mass_AL_one_dupli", "inv_mass_AL_one_dupli", 150, 0, 3);
  TH1F *inv_mass_AL_two_dupli = new TH1F("inv_mass_AL_two_dupli", "inv_mass_AL_two_dupli", 150, 0, 3);
  TH1F *inv_mass_AL_no_dupli = new TH1F("inv_mass_AL_no_dupli", "inv_mass_AL_no_dupli", 150, 0, 3);
  // (ambiguous tracks)
  TH1F *inv_mass_AL_also_K0 = new TH1F("inv_mass_AL_also_K0", "inv_mass_AL_also_K0", 150, 0, 3);
  TH1F *inv_mass_AL_no_ambi = new TH1F("inv_mass_AL_no_ambi", "inv_mass_AL_no_ambi", 150, 0, 3);
  // (both)
  TH1F *inv_mass_AL_no_dupli_no_ambi = new TH1F("inv_mass_AL_no_dupli_no_ambi", "inv_mass_AL_no_dupli_no_ambi", 150, 0, 3);
  // (anti-lambda signal)
  TH1F *inv_mass_AL_signal = new TH1F("inv_mass_AL_signal", "inv_mass_AL_signal", 150, 0, 3);
  TH1F *alt_mass_AL_signal = new TH1F("alt_mass_AL_signal", "alt_mass_AL_signal", 150, 0, 3);
  // (04) inv mass anti-lambdas with sexaquark signal cut
  TH1F *inv_mass_AL_sexa_signal = new TH1F("inv_mass_AL_sexa_signal", "inv_mass_AL_sexa_signal", 150, 0, 5);
  TH1F *alt_mass_AL_sexa_signal = new TH1F("alt_mass_AL_sexa_signal", "alt_mass_AL_sexa_signal", 150, 0, 5);

  // Pt of (K0,anti-lambda) (with/without) sexaquark signal cut
  TH1F *pt_K0 = new TH1F("pt_K0", "pt_K0", 150, 0, 5.);
  TH1F *pt_K0_sexa_signal = new TH1F("pt_K0_sexa_signal", "pt_K0_sexa_signal", 150, 0, 5.);
  TH1F *pt_AL = new TH1F("pt_AL", "pt_AL", 150, 0, 5.);
  TH1F *pt_AL_sexa_signal = new TH1F("pt_AL_sexa_signal", "pt_AL_sexa_signal", 150, 0, 5.);

  // radius of (K0,anti-lambda) (with/without) sexaquark signal cut
  TH1F *radius_K0 = new TH1F("radius_K0", "radius_K0", 180, 0, 180);
  TH1F *radius_K0_sexa_signal = new TH1F("radius_K0_sexa_signal", "radius_K0_sexa_signal", 180, 0, 180);
  TH1F *radius_AL = new TH1F("radius_AL", "radius_AL", 180, 0, 180);
  TH1F *radius_AL_sexa_signal = new TH1F("radius_AL_sexa_signal", "radius_AL_sexa_signal", 180, 0, 180);

  // define iterator
  TListIter *list_of_trees_it = new TListIter(list_of_trees);

  // loop over collected trees
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

    // load MC particles branches
    this_tree->SetBranchAddress("MC_PID", &MC_PID);
    this_tree->SetBranchAddress("MC_Mother", &MC_Mother);
    this_tree->SetBranchAddress("MC_isSignal", &MC_isSignal);
    // load track branches
    this_tree->SetBranchAddress("Idx_True", &Idx_True);
    this_tree->SetBranchAddress("Rec_Px", &Rec_Px);
    this_tree->SetBranchAddress("Rec_Py", &Rec_Py);
    this_tree->SetBranchAddress("Rec_Pz", &Rec_Pz);
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

    printf("Mom_PID Pos_PID Neg_PID Pos_isSignal Neg_isSignal  V0_isSignal V0_couldBeK0 V0_couldBeAL       M_asK0       M_asAL\n");

    // loop over events
    for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      /*** Count Found V0s ***/

      Bool_t exactlyOneTrackIsDupli;  // is exactly one track duplicated?
      Bool_t bothTracksAreDupli;      // are both tracks duplicated?

      Float_t V0_Pt;
      Float_t V0_Radius;

      Float_t mass_asK0;
      Float_t mass_asAL;

      Float_t isK0Signal;
      Float_t isALSignal;

      // (test)
      Float_t kMassProton = 0.938272;
      Float_t kMassPion = 0.139570;

      Float_t alt_pos_energy;

      Float_t alt_neg_energy_asK0;
      Float_t alt_energy_asK0;
      Float_t alt_Px;
      Float_t alt_Py;
      Float_t alt_Pz;
      Float_t alt_mass_asK0;

      Float_t alt_neg_energy_asAL;
      Float_t alt_energy_asAL;
      Float_t alt_mass_asAL;

      Int_t mother_pid;
      Int_t pos_pid;
      Int_t neg_pid;
      Int_t pos_issignal;
      Int_t neg_issignal;

      // loop over V0s
      for (Int_t v0 = 0; v0 < N_V0s; v0++) {

        // load variables
        exactlyOneTrackIsDupli = ((*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && !(*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]]) ||
                                 (!(*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]]);
        bothTracksAreDupli = (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Pos)[v0]]] && (*Rec_isDuplicate)[(*Idx_True)[(*Idx_Neg)[v0]]];

        V0_Pt = TMath::Sqrt((*V0_Px)[v0] * (*V0_Px)[v0] + (*V0_Py)[v0] * (*V0_Py)[v0]);
        V0_Radius = TMath::Sqrt((*V0_X)[v0] * (*V0_X)[v0] + (*V0_Y)[v0] * (*V0_Y)[v0]);

        mass_asK0 = TMath::Sqrt((*V0_E_asK0)[v0] * (*V0_E_asK0)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                (*V0_Pz)[v0] * (*V0_Pz)[v0]);
        mass_asAL = TMath::Sqrt((*V0_E_asAL)[v0] * (*V0_E_asAL)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                (*V0_Pz)[v0] * (*V0_Pz)[v0]);

        alt_pos_energy =
            TMath::Sqrt(kMassPion * kMassPion + (*Rec_Px)[(*Idx_Pos)[v0]] * (*Rec_Px)[(*Idx_Pos)[v0]] +
                        (*Rec_Py)[(*Idx_Pos)[v0]] * (*Rec_Py)[(*Idx_Pos)[v0]] + (*Rec_Pz)[(*Idx_Pos)[v0]] * (*Rec_Pz)[(*Idx_Pos)[v0]]);
        alt_neg_energy_asK0 =
            TMath::Sqrt(kMassPion * kMassPion + (*Rec_Px)[(*Idx_Neg)[v0]] * (*Rec_Px)[(*Idx_Neg)[v0]] +
                        (*Rec_Py)[(*Idx_Neg)[v0]] * (*Rec_Py)[(*Idx_Neg)[v0]] + (*Rec_Pz)[(*Idx_Neg)[v0]] * (*Rec_Pz)[(*Idx_Neg)[v0]]);
        alt_neg_energy_asAL =
            TMath::Sqrt(kMassProton * kMassProton + (*Rec_Px)[(*Idx_Neg)[v0]] * (*Rec_Px)[(*Idx_Neg)[v0]] +
                        (*Rec_Py)[(*Idx_Neg)[v0]] * (*Rec_Py)[(*Idx_Neg)[v0]] + (*Rec_Pz)[(*Idx_Neg)[v0]] * (*Rec_Pz)[(*Idx_Neg)[v0]]);

        alt_energy_asK0 = alt_pos_energy + alt_neg_energy_asK0;
        alt_energy_asAL = alt_pos_energy + alt_neg_energy_asAL;
        alt_Px = (*Rec_Px)[(*Idx_Pos)[v0]] + (*Rec_Px)[(*Idx_Neg)[v0]];
        alt_Py = (*Rec_Py)[(*Idx_Pos)[v0]] + (*Rec_Py)[(*Idx_Neg)[v0]];
        alt_Pz = (*Rec_Pz)[(*Idx_Pos)[v0]] + (*Rec_Pz)[(*Idx_Neg)[v0]];
        alt_mass_asK0 = TMath::Sqrt(alt_energy_asK0 * alt_energy_asK0 - alt_Px * alt_Px - alt_Py * alt_Py - alt_Pz * alt_Pz);
        alt_mass_asAL = TMath::Sqrt(alt_energy_asAL * alt_energy_asAL - alt_Px * alt_Px - alt_Py * alt_Py - alt_Pz * alt_Pz);

        pos_pid = (*MC_PID)[(*Idx_True)[(*Idx_Pos)[v0]]];
        neg_pid = (*MC_PID)[(*Idx_True)[(*Idx_Neg)[v0]]];
        pos_issignal = (Int_t)(*MC_isSignal)[(*Idx_True)[(*Idx_Pos)[v0]]];
        neg_issignal = (Int_t)(*MC_isSignal)[(*Idx_True)[(*Idx_Neg)[v0]]];

        if ((*MC_Mother)[(*Idx_True)[(*Idx_Pos)[v0]]] > 0) {
          mother_pid = (*MC_PID)[(*MC_Mother)[(*Idx_True)[(*Idx_Pos)[v0]]]];
          isK0Signal = mother_pid == 310;
          isALSignal = mother_pid == -3122;
        } else {
          isK0Signal = kFALSE;
          isALSignal = kFALSE;
        }

        if (isK0Signal) {
          inv_mass_K0_signal->Fill(mass_asK0);
          alt_mass_K0_signal->Fill(alt_mass_asK0);
        }

        if (isALSignal) {
          inv_mass_AL_signal->Fill(mass_asAL);
          alt_mass_AL_signal->Fill(alt_mass_asAL);
        }

        if ((*V0_couldBeK0)[v0]) {
          inv_mass_K0->Fill(mass_asK0);
          pt_K0->Fill(V0_Pt);
          radius_K0->Fill(V0_Radius);
          if (exactlyOneTrackIsDupli) {
            inv_mass_K0_one_dupli->Fill(mass_asK0);
          } else if (bothTracksAreDupli) {
            inv_mass_K0_two_dupli->Fill(mass_asK0);
          } else {
            inv_mass_K0_no_dupli->Fill(mass_asK0);
            if (!(*V0_couldBeAL)[v0]) {
              inv_mass_K0_no_dupli_no_ambi->Fill(mass_asK0);
            }
          }
          if ((*V0_couldBeAL)[v0]) {
            inv_mass_K0_also_AL->Fill(mass_asK0);
          } else {
            inv_mass_K0_no_ambi->Fill(mass_asK0);
          }
          if ((*V0_isSignal)[v0]) {
            inv_mass_K0_sexa_signal->Fill(mass_asK0);
            alt_mass_K0_sexa_signal->Fill(alt_mass_asK0);
            pt_K0_sexa_signal->Fill(V0_Pt);
            radius_K0_sexa_signal->Fill(V0_Radius);
          }
        }

        if ((*V0_couldBeAL)[v0]) {
          inv_mass_AL->Fill(mass_asAL);
          pt_AL->Fill(V0_Pt);
          radius_AL->Fill(V0_Radius);
          if (exactlyOneTrackIsDupli) {
            inv_mass_AL_one_dupli->Fill(mass_asAL);
          } else if (bothTracksAreDupli) {
            inv_mass_AL_two_dupli->Fill(mass_asAL);
          } else {
            inv_mass_AL_no_dupli->Fill(mass_asAL);
            if (!(*V0_couldBeK0)[v0]) {
              inv_mass_AL_no_dupli_no_ambi->Fill(mass_asAL);
            }
          }
          if ((*V0_couldBeK0)[v0]) {
            inv_mass_AL_also_K0->Fill(mass_asAL);
          } else {
            inv_mass_AL_no_ambi->Fill(mass_asAL);
          }
          if ((*V0_isSignal)[v0]) {
            printf("%7i %7i %7i %12i %12i %12i %12i %12i %12f %12f\n", mother_pid, pos_pid, neg_pid, pos_issignal, neg_issignal,
                   (Int_t)(*V0_isSignal)[v0], (Int_t)(*V0_couldBeK0)[v0], (Int_t)(*V0_couldBeAL)[v0], mass_asK0, mass_asAL);
            inv_mass_AL_sexa_signal->Fill(mass_asAL);
            alt_mass_AL_sexa_signal->Fill(alt_mass_asAL);
            pt_AL_sexa_signal->Fill(V0_Pt);
            radius_AL_sexa_signal->Fill(V0_Radius);
          }
        }

      }  // end of loop over V0s

    }  // end of loop over events

  }  // end of loop over trees

  /* Draw */

  // set style to everything
  SetMyStyle();
  gStyle->SetOptStat(0);

  // don't forget to add more histograms when necessary
  std::vector<TH1F *> list_hists = {inv_mass_K0,
                                    inv_mass_K0_sexa_signal,
                                    inv_mass_K0_one_dupli,
                                    inv_mass_K0_two_dupli,
                                    inv_mass_K0_no_dupli,
                                    inv_mass_K0_also_AL,
                                    inv_mass_K0_no_ambi,
                                    inv_mass_K0_no_dupli_no_ambi,
                                    inv_mass_K0_signal,
                                    alt_mass_K0_signal,
                                    inv_mass_AL,
                                    inv_mass_AL_sexa_signal,
                                    inv_mass_AL_one_dupli,
                                    inv_mass_AL_two_dupli,
                                    inv_mass_AL_no_dupli,
                                    inv_mass_AL_also_K0,
                                    inv_mass_AL_no_ambi,
                                    inv_mass_AL_no_dupli_no_ambi,
                                    inv_mass_AL_signal,
                                    alt_mass_AL_signal,
                                    pt_K0,
                                    pt_K0_sexa_signal,
                                    pt_AL,
                                    pt_AL_sexa_signal,
                                    radius_K0,
                                    radius_K0_sexa_signal,
                                    radius_AL,
                                    radius_AL_sexa_signal};

  for (TH1F *hist : list_hists) {
    SetMyHistStyle(hist);
    hist->SetLineWidth(0);
    hist->SetFillStyle(kSolid);
  }

  // the only hists to override this are:
  inv_mass_K0->SetFillStyle(0);
  inv_mass_K0->SetLineWidth(2);
  inv_mass_K0->SetLineColor(kBlack);
  inv_mass_K0->SetMinimum(1);
  inv_mass_K0->GetYaxis()->SetTitle("Counts");
  inv_mass_K0->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  inv_mass_AL->SetFillStyle(0);
  inv_mass_AL->SetLineWidth(2);
  inv_mass_AL->SetLineColor(kBlack);
  inv_mass_AL->SetMinimum(1);
  inv_mass_AL->GetYaxis()->SetTitle("Counts");
  inv_mass_AL->GetXaxis()->SetTitle("m(#pi^{+}#bar{p}) [GeV/c^{2}]");
  alt_mass_K0_signal->SetFillStyle(0);
  alt_mass_K0_signal->SetLineWidth(2);
  alt_mass_K0_signal->SetLineColor(myBlue);
  alt_mass_K0_signal->GetYaxis()->SetTitle("Counts");
  alt_mass_K0_signal->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV/c^{2}]");
  alt_mass_AL_signal->SetFillStyle(0);
  alt_mass_AL_signal->SetLineWidth(2);
  alt_mass_AL_signal->SetLineColor(myBlue);
  alt_mass_AL_signal->GetYaxis()->SetTitle("Counts");
  alt_mass_AL_signal->GetXaxis()->SetTitle("m(#pi^{+}#bar{p}) [GeV/c^{2}]");
  pt_K0->SetFillStyle(0);
  pt_K0->SetLineWidth(2);
  pt_K0->SetLineColor(kBlack);
  pt_K0->SetMinimum(1);
  pt_K0->GetYaxis()->SetTitle("Counts");
  pt_K0->GetXaxis()->SetTitle("p_{T}(#pi^{+}#pi^{-}) [GeV/c]");
  pt_AL->SetFillStyle(0);
  pt_AL->SetLineWidth(2);
  pt_AL->SetLineColor(kBlack);
  pt_AL->SetMinimum(1);
  pt_AL->GetYaxis()->SetTitle("Counts");
  pt_AL->GetXaxis()->SetTitle("p_{T}(#pi^{+}#bar{p}) [GeV/c]");
  radius_K0->SetFillStyle(0);
  radius_K0->SetLineWidth(2);
  radius_K0->SetLineColor(kBlack);
  radius_K0->SetMinimum(1);
  radius_K0->GetYaxis()->SetTitle("Counts");
  radius_K0->GetXaxis()->SetTitle("Fiducial Radius(#pi^{+}#pi^{-}) [cm]");
  radius_AL->SetFillStyle(0);
  radius_AL->SetLineWidth(2);
  radius_AL->SetLineColor(kBlack);
  radius_AL->SetMinimum(1);
  radius_AL->GetYaxis()->SetTitle("Counts");
  radius_AL->GetXaxis()->SetTitle("Fiducial Radius(#pi^{+}#bar{p}) [cm]");

  TString output_filename;
  TCanvas *c = new TCanvas("c", "c", 1080, 1080);

  c->SetLogy(1);
  // (01) inv mass K0 without signal cut
  // (02) inv mass K0 with signal cut
  inv_mass_K0->Draw();
  inv_mass_K0_sexa_signal->SetFillColor(myRed);
  inv_mass_K0_sexa_signal->Draw("SAME");
  inv_mass_K0_sexa_signal->Draw("AXIS SAME");
  TLegend *l = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  // l->SetHeader("#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected ", "C");
  l->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l->AddEntry(inv_mass_K0, "K^{0}_{S} candidates", "l");
  l->AddEntry(inv_mass_K0_sexa_signal, "#bar{S} signal", "f");
  l->SetTextSize(0.02);
  l->SetBorderSize(0);
  l->SetFillStyle(1001);
  l->Draw();
  output_filename = inv_mass_K0->GetName();
  c->Print(output_dir + output_filename + ".png");

  inv_mass_K0->Draw();
  alt_mass_K0_sexa_signal->SetFillColor(myRed);
  alt_mass_K0_sexa_signal->Draw("SAME");
  alt_mass_K0_sexa_signal->Draw("AXIS SAME");
  TLegend *l7 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  // l->SetHeader("#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected ", "C");
  l7->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l7->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l7->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l7->AddEntry(inv_mass_K0, "K^{0}_{S} candidates", "l");
  l7->AddEntry(alt_mass_K0_sexa_signal, "#bar{S} signal", "f");
  l7->SetTextSize(0.02);
  l7->SetBorderSize(0);
  l7->SetFillStyle(1001);
  l7->Draw();
  output_filename = alt_mass_K0_sexa_signal->GetName();
  c->Print(output_dir + output_filename + ".png");
  c->SetLogy(0);

  // (duplicated tracks)
  inv_mass_K0->Draw();
  inv_mass_K0_one_dupli->SetFillColor(myBlue);
  inv_mass_K0_one_dupli->Draw("SAME");
  inv_mass_K0_two_dupli->SetFillColor(myGreen);
  inv_mass_K0_two_dupli->Draw("SAME");
  inv_mass_K0_two_dupli->Draw("AXIS SAME");
  DrawText("one daughter is duplicated", myBlue, 0.8);
  DrawText("both daughters are duplicated", myGreen, 0.75);
  output_filename = (TString)inv_mass_K0->GetName() + "_dupli";
  c->Print(output_dir + output_filename + ".png");

  inv_mass_K0->Draw();
  inv_mass_K0_no_dupli->SetFillColor(myYellow);
  inv_mass_K0_no_dupli->Draw("SAME");
  inv_mass_K0_no_dupli->Draw("AXIS SAME");
  DrawText("none of the daughters is duplicated", myYellow, 0.8);
  output_filename = inv_mass_K0_no_dupli->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (ambiguous tracks)
  inv_mass_K0->Draw();
  inv_mass_K0_also_AL->SetFillColor(myOrange);
  inv_mass_K0_also_AL->Draw("SAME");
  inv_mass_K0_also_AL->Draw("AXIS SAME");
  DrawText("could also be anti-lambda", myOrange, 0.8);
  output_filename = (TString)inv_mass_K0->GetName() + "_ambi";
  c->Print(output_dir + output_filename + ".png");

  inv_mass_K0->Draw();
  inv_mass_K0_no_ambi->SetFillColor(myViolet);
  inv_mass_K0_no_ambi->Draw("SAME");
  inv_mass_K0_no_ambi->Draw("AXIS SAME");
  DrawText("exclusively neutral kaon short", myViolet, 0.8);
  output_filename = inv_mass_K0_no_ambi->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (both)
  inv_mass_K0->Draw();
  inv_mass_K0_no_dupli_no_ambi->SetFillColor(myMagenta);
  inv_mass_K0_no_dupli_no_ambi->Draw("SAME");
  inv_mass_K0_no_dupli_no_ambi->Draw("AXIS SAME");
  DrawText("no dupli no ambi", myMagenta, 0.8);
  output_filename = inv_mass_K0_no_dupli_no_ambi->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (K0 signal)
  inv_mass_K0->Draw();
  inv_mass_K0_signal->SetFillColor(myBlack);
  inv_mass_K0_signal->Draw("SAME");
  inv_mass_K0_signal->Draw("AXIS SAME");
  DrawText("K0S signal", myBlack, 0.8);
  output_filename = inv_mass_K0_signal->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (testing another mass definition)
  alt_mass_K0_signal->Draw();
  alt_mass_K0_signal->Draw("AXIS SAME");
  DrawText("(taking momentum from TPC)", myBlack, 0.8);
  output_filename = alt_mass_K0_signal->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (03) inv mass anti-lambdas without signal cut
  // (04) inv mass anti-lambdas with signal cut
  c->SetLogy(1);
  inv_mass_AL->Draw();
  inv_mass_AL_sexa_signal->SetFillColor(myRed);
  inv_mass_AL_sexa_signal->Draw("SAME");
  inv_mass_AL_sexa_signal->Draw("AXIS SAME");
  TLegend *l2 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l2->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l2->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l2->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l2->AddEntry(inv_mass_AL, "#bar{#Lambda} candidates", "l");
  l2->AddEntry(inv_mass_AL_sexa_signal, "#bar{S} signal", "f");
  l2->SetTextSize(0.02);
  l2->SetBorderSize(0);
  l2->SetFillStyle(1001);
  l2->Draw();
  output_filename = inv_mass_AL->GetName();
  c->Print(output_dir + output_filename + ".png");

  inv_mass_AL->Draw();
  alt_mass_AL_sexa_signal->SetFillColor(myRed);
  alt_mass_AL_sexa_signal->Draw("SAME");
  alt_mass_AL_sexa_signal->Draw("AXIS SAME");
  TLegend *l8 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l8->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l8->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l8->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l8->AddEntry(inv_mass_AL, "#bar{#Lambda} candidates", "l");
  l8->AddEntry(alt_mass_AL_sexa_signal, "#bar{S} signal", "f");
  l8->SetTextSize(0.02);
  l8->SetBorderSize(0);
  l8->SetFillStyle(1001);
  l8->Draw();
  output_filename = alt_mass_AL_sexa_signal->GetName();
  c->Print(output_dir + output_filename + ".png");
  c->SetLogy(0);

  // (duplicated tracks)
  inv_mass_AL->Draw();
  inv_mass_AL_one_dupli->SetFillColor(myBlue);
  inv_mass_AL_one_dupli->Draw("SAME");
  inv_mass_AL_two_dupli->SetFillColor(myGreen);
  inv_mass_AL_two_dupli->Draw("SAME");
  inv_mass_AL_two_dupli->Draw("AXIS SAME");
  DrawText("one daughter is duplicated", myBlue, 0.8);
  DrawText("both daughters are duplicated", myGreen, 0.75);
  output_filename = (TString)inv_mass_AL->GetName() + "_dupli";
  c->Print(output_dir + output_filename + ".png");

  inv_mass_AL->Draw();
  inv_mass_AL_no_dupli->SetFillColor(myYellow);
  inv_mass_AL_no_dupli->Draw("SAME");
  inv_mass_AL_no_dupli->Draw("AXIS SAME");
  DrawText("none of the daughters is duplicated", myYellow, 0.8);
  output_filename = inv_mass_AL_no_dupli->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (ambiguous tracks)
  inv_mass_AL->Draw();
  inv_mass_AL_also_K0->SetFillColor(myOrange);
  inv_mass_AL_also_K0->Draw("SAME");
  inv_mass_AL_also_K0->Draw("AXIS SAME");
  DrawText("could also be K0", myOrange, 0.8);
  output_filename = (TString)inv_mass_AL->GetName() + "_ambi";
  c->Print(output_dir + output_filename + ".png");

  inv_mass_AL->Draw();
  inv_mass_AL_no_ambi->SetFillColor(myViolet);
  inv_mass_AL_no_ambi->Draw("SAME");
  inv_mass_AL_no_ambi->Draw("AXIS SAME");
  DrawText("exclusively anti-lambda", myViolet, 0.8);
  output_filename = inv_mass_AL_no_ambi->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (both)
  inv_mass_AL->Draw();
  inv_mass_AL_no_dupli_no_ambi->SetFillColor(myMagenta);
  inv_mass_AL_no_dupli_no_ambi->Draw("SAME");
  inv_mass_AL_no_dupli_no_ambi->Draw("AXIS SAME");
  DrawText("no dupli no ambi", myMagenta, 0.8);
  output_filename = inv_mass_AL_no_dupli_no_ambi->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (both)
  inv_mass_AL->Draw();
  inv_mass_AL_signal->SetFillColor(myBlack);
  inv_mass_AL_signal->Draw("SAME");
  inv_mass_AL_signal->Draw("AXIS SAME");
  DrawText("anti-lambda signal", myBlack, 0.8);
  output_filename = inv_mass_AL_signal->GetName();
  c->Print(output_dir + output_filename + ".png");

  // (testing another mass definition)
  alt_mass_AL_signal->Draw();
  alt_mass_AL_signal->Draw("AXIS SAME");
  DrawText("(taking momentum from TPC)", myBlack, 0.8);
  output_filename = alt_mass_AL_signal->GetName();
  c->Print(output_dir + output_filename + ".png");

  // pt
  c->SetLogy(1);
  pt_K0->Draw();
  pt_K0_sexa_signal->SetFillColor(myRed);
  pt_K0_sexa_signal->Draw("SAME");
  pt_K0_sexa_signal->Draw("AXIS SAME");
  output_filename = pt_K0->GetName();
  TLegend *l3 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l3->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l3->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l3->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l3->AddEntry(pt_K0, "K^{0}_{S} candidates", "l");
  l3->AddEntry(pt_K0_sexa_signal, "#bar{S} signal", "f");
  l3->SetTextSize(0.02);
  l3->SetBorderSize(0);
  l3->SetFillStyle(1001);
  l3->Draw();
  c->Print(output_dir + output_filename + ".png");

  pt_AL->Draw();
  pt_AL_sexa_signal->SetFillColor(myRed);
  pt_AL_sexa_signal->Draw("SAME");
  pt_AL_sexa_signal->Draw("AXIS SAME");
  output_filename = pt_AL->GetName();
  TLegend *l4 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l4->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l4->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l4->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l4->AddEntry(pt_AL, "#bar{#Lambda} candidates", "l");
  l4->AddEntry(pt_AL_sexa_signal, "#bar{S} signal", "f");
  l4->SetTextSize(0.02);
  l4->SetBorderSize(0);
  l4->SetFillStyle(0);
  l4->Draw();
  c->Print(output_dir + output_filename + ".png");

  // pt
  radius_K0->Draw();
  radius_K0_sexa_signal->SetFillColor(myRed);
  radius_K0_sexa_signal->Draw("SAME");
  radius_K0_sexa_signal->Draw("AXIS SAME");
  output_filename = radius_K0->GetName();
  TLegend *l5 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l5->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l5->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l5->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l5->AddEntry(radius_K0, "K^{0}_{S} candidates", "l");
  l5->AddEntry(radius_K0_sexa_signal, "#bar{S} signal", "f");
  l5->SetTextSize(0.02);
  l5->SetBorderSize(0);
  l5->SetFillStyle(1001);
  l5->Draw();
  c->Print(output_dir + output_filename + ".png");

  radius_AL->Draw();
  radius_AL_sexa_signal->SetFillColor(myRed);
  radius_AL_sexa_signal->Draw("SAME");
  radius_AL_sexa_signal->Draw("AXIS SAME");
  output_filename = radius_AL->GetName();
  TLegend *l6 = new TLegend(0.5, 0.75, 0.9, 0.9, "", "NDC");  // x1,y1,x2,y2
  l6->AddEntry((TObject *)0, "ALICE Run 2 Simulation", "");
  l6->AddEntry((TObject *)0, "#bar{S} + n #rightarrow #bar{#Lambda} + K^{0}_{S} injected over HIJING", "");
  l6->AddEntry((TObject *)0, "Pb-Pb, #sqrt{s}=5 TeV", "");
  l6->AddEntry(radius_AL, "#bar{#Lambda} candidates", "l");
  l6->AddEntry(radius_AL_sexa_signal, "#bar{S} signal", "f");
  l6->SetTextSize(0.02);
  l6->SetBorderSize(0);
  l6->SetFillStyle(1001);
  l6->Draw();
  c->Print(output_dir + output_filename + ".png");
  c->SetLogy(0);
}
