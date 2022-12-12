#include "macros/include/Headers.hxx"
#include "macros/include/TreeFunctions.hxx"
#include "macros/include/Utilities.hxx"

/*** Declaration of functions ***/

Double_t TwoLinesDCA(TVector3 pos1, TVector3 dir1, TVector3 pos2, TVector3 dir2, TVector3 &P1, TVector3 &P2);

/*** Main ***/

//________________________________________________________________________
void SexaquarkFinder(  //
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595_v1/AnalysisResults_TrueV0s_000.root",
    TString output_filename = "output/signal+bkg/297595_v1/SexaquarkResults_TrueV0s_000.root") {
  // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/AnalysisResults_OfficialV0s_000.root",
  // TString output_filename = "output/signal+bkg/SexaquarkResults_OfficialV0s_000.root") {
  //
  // Just another macro, but more important than the others...
  //

  /*** Process Input ***/

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  Event_tt input_event;

  /*** Prepare Output ***/

  // define output file and tree
  TFile *output_file = new TFile(output_filename, "RECREATE");
  TTree *output_tree = new TTree("Sexaquarks", "Sexaquark Candidates");

  // output variables (these are not vectors, we're converting into a flat format!)
  /* Anti-Sexaquark Candidate */  //
  Int_t event;                    // event index
  Int_t Idx_V0A;                  // index of first daughter
  Int_t Idx_V0B;                  // index of second daughter
  Float_t Sexa_E;                 // energy of anti-sexaquark candidate
  Float_t Sexa_Px;                // x-component of anti-sexaquark candidate
  Float_t Sexa_Py;                // y-component of anti-sexaquark candidate
  Float_t Sexa_Pz;                // z-component of anti-sexaquark candidate
  Float_t Sexa_M;                 // inv. mass of anti-sexaquark candidate
  Float_t Sexa_X;                 // x-coordinate of anti-sexaquark candidate (secondary vertex)
  Float_t Sexa_Y;                 // y-coordinate of anti-sexaquark candidate (secondary vertex)
  Float_t Sexa_Z;                 // z-coordinate of anti-sexaquark candidate (secondary vertex)
  Float_t Sexa_DCA;               // distance of closest approach after extrapolation of V0s
  Bool_t Sexa_isSignal;           // kTRUE if signal, kFALSE if background
  /* First V0: "V0A" */           //
  Int_t V0A_Idx_Pos;
  Int_t V0A_Idx_Neg;
  Float_t V0A_Px;
  Float_t V0A_Py;
  Float_t V0A_Pz;
  Float_t V0A_X;
  Float_t V0A_Y;
  Float_t V0A_Z;
  Float_t V0A_Pos_Px;
  Float_t V0A_Pos_Py;
  Float_t V0A_Pos_Pz;
  Float_t V0A_Neg_Px;
  Float_t V0A_Neg_Py;
  Float_t V0A_Neg_Pz;
  Bool_t V0A_isSignal;
  Float_t V0A_E_asK0;
  Float_t V0A_E_asAL;
  Bool_t V0A_couldBeK0;
  Bool_t V0A_couldBeAL;
  Bool_t V0A_onFlyStatus;
  Float_t V0A_Chi2;
  Float_t V0A_DCA_Daughters;
  Float_t V0A_IP_wrtPV;
  Float_t V0A_CPA_wrtPV;
  Float_t V0A_ArmAlpha;
  Float_t V0A_ArmPt;
  Float_t V0A_DecayLength;
  Bool_t V0A_Pos_isDuplicate;
  Bool_t V0A_Neg_isDuplicate;
  Float_t V0A_Pos_Rec_Px;
  Float_t V0A_Pos_Rec_Py;
  Float_t V0A_Pos_Rec_Pz;
  Float_t V0A_Neg_Rec_Px;
  Float_t V0A_Neg_Rec_Py;
  Float_t V0A_Neg_Rec_Pz;
  /* Second V0: "V0B" */  //
  Int_t V0B_Idx_Pos;
  Int_t V0B_Idx_Neg;
  Float_t V0B_Px;
  Float_t V0B_Py;
  Float_t V0B_Pz;
  Float_t V0B_X;
  Float_t V0B_Y;
  Float_t V0B_Z;
  Float_t V0B_Pos_Px;
  Float_t V0B_Pos_Py;
  Float_t V0B_Pos_Pz;
  Float_t V0B_Neg_Px;
  Float_t V0B_Neg_Py;
  Float_t V0B_Neg_Pz;
  Bool_t V0B_isSignal;
  Float_t V0B_E_asK0;
  Float_t V0B_E_asAL;
  Bool_t V0B_couldBeK0;
  Bool_t V0B_couldBeAL;
  Bool_t V0B_onFlyStatus;
  Float_t V0B_Chi2;
  Float_t V0B_DCA_Daughters;
  Float_t V0B_IP_wrtPV;
  Float_t V0B_CPA_wrtPV;
  Float_t V0B_ArmAlpha;
  Float_t V0B_ArmPt;
  Float_t V0B_DecayLength;
  Bool_t V0B_Pos_isDuplicate;
  Bool_t V0B_Neg_isDuplicate;
  Float_t V0B_Pos_Rec_Px;
  Float_t V0B_Pos_Rec_Py;
  Float_t V0B_Pos_Rec_Pz;
  Float_t V0B_Neg_Rec_Px;
  Float_t V0B_Neg_Rec_Py;
  Float_t V0B_Neg_Rec_Pz;

  // link output branches

  /* Anti-Sexaquark Candidate */

  output_tree->Branch("event", &event);
  output_tree->Branch("Idx_V0A", &Idx_V0A);
  output_tree->Branch("Idx_V0B", &Idx_V0B);
  output_tree->Branch("E", &Sexa_E);
  output_tree->Branch("Px", &Sexa_Px);
  output_tree->Branch("Py", &Sexa_Py);
  output_tree->Branch("Pz", &Sexa_Pz);
  output_tree->Branch("M", &Sexa_M);
  output_tree->Branch("X", &Sexa_X);
  output_tree->Branch("Y", &Sexa_Y);
  output_tree->Branch("Z", &Sexa_Z);
  output_tree->Branch("DCA", &Sexa_DCA);
  output_tree->Branch("isSignal", &Sexa_isSignal);

  /* First V0: "V0A" */

  output_tree->Branch("V0A_Idx_Pos", &V0A_Idx_Pos);
  output_tree->Branch("V0A_Idx_Neg", &V0A_Idx_Neg);
  output_tree->Branch("V0A_Px", &V0A_Px);
  output_tree->Branch("V0A_Py", &V0A_Py);
  output_tree->Branch("V0A_Pz", &V0A_Pz);
  output_tree->Branch("V0A_X", &V0A_X);
  output_tree->Branch("V0A_Y", &V0A_Y);
  output_tree->Branch("V0A_Z", &V0A_Z);
  output_tree->Branch("V0A_Pos_Px", &V0A_Pos_Px);
  output_tree->Branch("V0A_Pos_Py", &V0A_Pos_Py);
  output_tree->Branch("V0A_Pos_Pz", &V0A_Pos_Pz);
  output_tree->Branch("V0A_Neg_Px", &V0A_Neg_Px);
  output_tree->Branch("V0A_Neg_Py", &V0A_Neg_Py);
  output_tree->Branch("V0A_Neg_Pz", &V0A_Neg_Pz);
  output_tree->Branch("V0A_isSignal", &V0A_isSignal);
  output_tree->Branch("V0A_E_asK0", &V0A_E_asK0);
  output_tree->Branch("V0A_E_asAL", &V0A_E_asAL);
  output_tree->Branch("V0A_couldBeK0", &V0A_couldBeK0);
  output_tree->Branch("V0A_couldBeAL", &V0A_couldBeAL);
  output_tree->Branch("V0A_onFlyStatus", &V0A_onFlyStatus);
  output_tree->Branch("V0A_Chi2", &V0A_Chi2);
  output_tree->Branch("V0A_DCA_Daughters", &V0A_DCA_Daughters);
  output_tree->Branch("V0A_IP_wrtPV", &V0A_IP_wrtPV);
  output_tree->Branch("V0A_CPA_wrtPV", &V0A_CPA_wrtPV);
  output_tree->Branch("V0A_ArmAlpha", &V0A_ArmAlpha);
  output_tree->Branch("V0A_ArmPt", &V0A_ArmPt);
  output_tree->Branch("V0A_DecayLength", &V0A_DecayLength);

  output_tree->Branch("V0A_Pos_isDuplicate", &V0A_Pos_isDuplicate);
  output_tree->Branch("V0A_Neg_isDuplicate", &V0A_Neg_isDuplicate);
  output_tree->Branch("V0A_Pos_Rec_Px", &V0A_Pos_Rec_Px);
  output_tree->Branch("V0A_Pos_Rec_Py", &V0A_Pos_Rec_Py);
  output_tree->Branch("V0A_Pos_Rec_Pz", &V0A_Pos_Rec_Pz);
  output_tree->Branch("V0A_Neg_Rec_Px", &V0A_Neg_Rec_Px);
  output_tree->Branch("V0A_Neg_Rec_Py", &V0A_Neg_Rec_Py);
  output_tree->Branch("V0A_Neg_Rec_Pz", &V0A_Neg_Rec_Pz);

  /* Second V0: "V0B" */

  output_tree->Branch("V0B_Idx_Pos", &V0B_Idx_Pos);
  output_tree->Branch("V0B_Idx_Neg", &V0B_Idx_Neg);
  output_tree->Branch("V0B_Px", &V0B_Px);
  output_tree->Branch("V0B_Py", &V0B_Py);
  output_tree->Branch("V0B_Pz", &V0B_Pz);
  output_tree->Branch("V0B_X", &V0B_X);
  output_tree->Branch("V0B_Y", &V0B_Y);
  output_tree->Branch("V0B_Z", &V0B_Z);
  output_tree->Branch("V0B_Pos_Px", &V0B_Pos_Px);
  output_tree->Branch("V0B_Pos_Py", &V0B_Pos_Py);
  output_tree->Branch("V0B_Pos_Pz", &V0B_Pos_Pz);
  output_tree->Branch("V0B_Neg_Px", &V0B_Neg_Px);
  output_tree->Branch("V0B_Neg_Py", &V0B_Neg_Py);
  output_tree->Branch("V0B_Neg_Pz", &V0B_Neg_Pz);
  output_tree->Branch("V0B_isSignal", &V0B_isSignal);
  output_tree->Branch("V0B_E_asK0", &V0B_E_asK0);
  output_tree->Branch("V0B_E_asAL", &V0B_E_asAL);
  output_tree->Branch("V0B_couldBeK0", &V0B_couldBeK0);
  output_tree->Branch("V0B_couldBeAL", &V0B_couldBeAL);
  output_tree->Branch("V0B_onFlyStatus", &V0B_onFlyStatus);
  output_tree->Branch("V0B_Chi2", &V0B_Chi2);
  output_tree->Branch("V0B_DCA_Daughters", &V0B_DCA_Daughters);
  output_tree->Branch("V0B_IP_wrtPV", &V0B_IP_wrtPV);
  output_tree->Branch("V0B_CPA_wrtPV", &V0B_CPA_wrtPV);
  output_tree->Branch("V0B_ArmAlpha", &V0B_ArmAlpha);
  output_tree->Branch("V0B_ArmPt", &V0B_ArmPt);
  output_tree->Branch("V0B_DecayLength", &V0B_DecayLength);

  output_tree->Branch("V0B_Pos_isDuplicate", &V0B_Pos_isDuplicate);
  output_tree->Branch("V0B_Neg_isDuplicate", &V0B_Neg_isDuplicate);
  output_tree->Branch("V0B_Pos_Rec_Px", &V0B_Pos_Rec_Px);
  output_tree->Branch("V0B_Pos_Rec_Py", &V0B_Pos_Rec_Py);
  output_tree->Branch("V0B_Pos_Rec_Pz", &V0B_Pos_Rec_Pz);
  output_tree->Branch("V0B_Neg_Rec_Px", &V0B_Neg_Rec_Px);
  output_tree->Branch("V0B_Neg_Rec_Py", &V0B_Neg_Rec_Py);
  output_tree->Branch("V0B_Neg_Rec_Pz", &V0B_Neg_Rec_Pz);

  /*** Start ***/

  // define iterator
  TListIter *list_of_trees_it = new TListIter(list_of_trees);

  // loop over collected trees
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

    // load branches
    LoadBranches(this_tree, input_event);

    // define variables & objects
    TVector3 pos_A;
    TVector3 dir_A;
    TVector3 pos_B;
    TVector3 dir_B;

    TVector3 pca_A;  // point of closest approach
    TVector3 pca_B;
    TVector3 secondary_vtx;
    Double_t half_distance;

    // loop over events
    for (event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      // (debug)
      printf("SexaquarkFinder :: Processing event #%i\n", event);
      /*
      printf("V0s FROM ANALYSIS RESULTS\n");
      printf("\n");
      printf("Event Index  V0_Px  V0_Py  V0_Pz   V0_X   V0_Y   V0_Z P_Idx   P_Px   P_Py   P_Pz N_Idx   N_Px   N_Py   N_Pz\n");
      */

      // (debug)
      // printf("SexaquarkFinder :: ANTI-SEXAQUARK CANDIDATES\n");
      // printf("SexaquarkFinder :: \n");
      // printf("SexaquarkFinder :: Idx_V0A Idx_V0B    SV_X    SV_Y    SV_Z     DCA       M\n");

      // loop over V0s in an event
      for (Int_t aa = 0; aa < input_event.N_V0s - 1; aa++) {
        for (Int_t bb = aa + 1; bb < input_event.N_V0s; bb++) {

          /*
          // (debug)
          printf("%5i %5i %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %5i %6.2f %6.2f %6.2f %5i %6.2f %6.2f %6.2f\n", event, aa, (*V0_Px)[aa],
                 (*V0_Py)[aa], (*V0_Pz)[aa], (*V0_X)[aa], (*V0_Y)[aa], (*V0_Z)[aa], (*Idx_Pos)[aa], (*Pos_Px)[aa], (*Pos_Py)[aa],
                 (*Pos_Pz)[aa], (*Idx_Neg)[aa], (*Neg_Px)[aa], (*Neg_Py)[aa], (*Neg_Pz)[aa]);
          */

          // (cut) first, check that all four granddaughters are different
          if ((*input_event.Idx_Pos)[aa] == (*input_event.Idx_Pos)[bb] || (*input_event.Idx_Neg)[aa] == (*input_event.Idx_Neg)[bb]) {
            continue;
          }

          // (cut) we're not allowing ambiguity!
          // both particles must be different V0s!
          if (((*input_event.V0_couldBeK0)[aa] && (*input_event.V0_couldBeAL)[bb]) ||
              ((*input_event.V0_couldBeK0)[bb] && (*input_event.V0_couldBeAL)[aa])) {
            // nothing
          } else {
            continue;
          }

          pos_A.SetXYZ((*input_event.V0_X)[aa], (*input_event.V0_Y)[aa], (*input_event.V0_Z)[aa]);
          dir_A.SetXYZ((*input_event.V0_Px)[aa], (*input_event.V0_Py)[aa], (*input_event.V0_Pz)[aa]);

          pos_B.SetXYZ((*input_event.V0_X)[bb], (*input_event.V0_Y)[bb], (*input_event.V0_Z)[bb]);
          dir_B.SetXYZ((*input_event.V0_Px)[bb], (*input_event.V0_Py)[bb], (*input_event.V0_Pz)[bb]);

          // (assign) from both momentum and coordinates of V0s,
          // extrapolate them to a common vertex,
          // if there's none, calculate DCA and extract PCAs
          Sexa_DCA = TwoLinesDCA(pos_A, dir_A, pos_B, dir_B, pca_A, pca_B);

          // calculate the secondary_vtx point between PCAs
          half_distance = (pca_B - pca_A).Mag() / 2.;
          secondary_vtx = pca_A + half_distance * (pca_B - pca_A).Unit();

          // (assign) both V0s must be signal for the candidate to be signal
          Sexa_isSignal = (*input_event.V0_isSignal)[aa] && (*input_event.V0_isSignal)[bb];

          // (assign) rest of variables
          Idx_V0A = aa;
          Idx_V0B = bb;
          if ((*input_event.V0_couldBeK0)[aa]) {
            Sexa_E = (*input_event.V0_E_asK0)[aa] + (*input_event.V0_E_asAL)[bb] - 0.939565;  // assuming neutron
          }
          if ((*input_event.V0_couldBeAL)[aa]) {
            Sexa_E = (*input_event.V0_E_asAL)[aa] + (*input_event.V0_E_asK0)[bb] - 0.939565;  // assuming neutron
          }
          Sexa_Px = (dir_A + dir_B).X();
          Sexa_Py = (dir_A + dir_B).Y();
          Sexa_Pz = (dir_A + dir_B).Z();
          Float_t pre_mass = Sexa_E * Sexa_E - Sexa_Px * Sexa_Px - Sexa_Py * Sexa_Py - Sexa_Pz * Sexa_Pz;
          Sexa_M = pre_mass >= 0 ? TMath::Sqrt(pre_mass) : -1;
          Sexa_X = secondary_vtx.X();
          Sexa_Y = secondary_vtx.Y();
          Sexa_Z = secondary_vtx.Z();

          // (assign) first V0: "V0A"
          V0A_Idx_Pos = (*input_event.Idx_Pos)[aa];
          V0A_Idx_Neg = (*input_event.Idx_Neg)[aa];
          V0A_Px = (*input_event.V0_Px)[aa];
          V0A_Py = (*input_event.V0_Py)[aa];
          V0A_Pz = (*input_event.V0_Pz)[aa];
          V0A_X = (*input_event.V0_X)[aa];
          V0A_Y = (*input_event.V0_Y)[aa];
          V0A_Z = (*input_event.V0_Z)[aa];
          V0A_Pos_Px = (*input_event.Pos_Px)[aa];
          V0A_Pos_Py = (*input_event.Pos_Py)[aa];
          V0A_Pos_Pz = (*input_event.Pos_Pz)[aa];
          V0A_Neg_Px = (*input_event.Neg_Px)[aa];
          V0A_Neg_Py = (*input_event.Neg_Py)[aa];
          V0A_Neg_Pz = (*input_event.Neg_Pz)[aa];
          V0A_isSignal = (*input_event.V0_isSignal)[aa];
          V0A_E_asK0 = (*input_event.V0_E_asK0)[aa];
          V0A_E_asAL = (*input_event.V0_E_asAL)[aa];
          V0A_couldBeK0 = (*input_event.V0_couldBeK0)[aa];
          V0A_couldBeAL = (*input_event.V0_couldBeAL)[aa];
          V0A_onFlyStatus = (*input_event.V0_onFlyStatus)[aa];
          V0A_Chi2 = (*input_event.V0_Chi2)[aa];
          V0A_DCA_Daughters = (*input_event.V0_DCA_Daughters)[aa];
          V0A_IP_wrtPV = (*input_event.V0_IP_wrtPV)[aa];
          V0A_CPA_wrtPV = (*input_event.V0_CPA_wrtPV)[aa];
          V0A_ArmAlpha = (*input_event.V0_ArmAlpha)[aa];
          V0A_ArmPt = (*input_event.V0_ArmPt)[aa];
          V0A_DecayLength = (*input_event.V0_DecayLength)[aa];

          V0A_Pos_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Pos)[aa]];
          V0A_Neg_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Neg)[aa]];
          V0A_Pos_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Pos)[aa]];
          V0A_Pos_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Pos)[aa]];
          V0A_Pos_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Pos)[aa]];
          V0A_Neg_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Neg)[aa]];
          V0A_Neg_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Neg)[aa]];
          V0A_Neg_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Neg)[aa]];

          // (assign) second V0: "V0B"
          V0B_Idx_Pos = (*input_event.Idx_Pos)[bb];
          V0B_Idx_Neg = (*input_event.Idx_Neg)[bb];
          V0B_Px = (*input_event.V0_Px)[bb];
          V0B_Py = (*input_event.V0_Py)[bb];
          V0B_Pz = (*input_event.V0_Pz)[bb];
          V0B_X = (*input_event.V0_X)[bb];
          V0B_Y = (*input_event.V0_Y)[bb];
          V0B_Z = (*input_event.V0_Z)[bb];
          V0B_Pos_Px = (*input_event.Pos_Px)[bb];
          V0B_Pos_Py = (*input_event.Pos_Py)[bb];
          V0B_Pos_Pz = (*input_event.Pos_Pz)[bb];
          V0B_Neg_Px = (*input_event.Neg_Px)[bb];
          V0B_Neg_Py = (*input_event.Neg_Py)[bb];
          V0B_Neg_Pz = (*input_event.Neg_Pz)[bb];
          V0B_isSignal = (*input_event.V0_isSignal)[bb];
          V0B_E_asK0 = (*input_event.V0_E_asK0)[bb];
          V0B_E_asAL = (*input_event.V0_E_asAL)[bb];
          V0B_couldBeK0 = (*input_event.V0_couldBeK0)[bb];
          V0B_couldBeAL = (*input_event.V0_couldBeAL)[bb];
          V0B_onFlyStatus = (*input_event.V0_onFlyStatus)[bb];
          V0B_Chi2 = (*input_event.V0_Chi2)[bb];
          V0B_DCA_Daughters = (*input_event.V0_DCA_Daughters)[bb];
          V0B_IP_wrtPV = (*input_event.V0_IP_wrtPV)[bb];
          V0B_CPA_wrtPV = (*input_event.V0_CPA_wrtPV)[bb];
          V0B_ArmAlpha = (*input_event.V0_ArmAlpha)[bb];
          V0B_ArmPt = (*input_event.V0_ArmPt)[bb];
          V0B_DecayLength = (*input_event.V0_DecayLength)[bb];

          V0B_Pos_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Pos)[bb]];
          V0B_Neg_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Neg)[bb]];
          V0B_Pos_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Pos)[bb]];
          V0B_Pos_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Pos)[bb]];
          V0B_Pos_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Pos)[bb]];
          V0B_Neg_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Neg)[bb]];
          V0B_Neg_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Neg)[bb]];
          V0B_Neg_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Neg)[bb]];

          // (debug)
          // printf("SexaquarkFinder :: %7i %7i %7.3f %7.3f %7.3f %7.3f %7.3f\n", aa, bb, secondary_vtx.X(), secondary_vtx.Y(),
          //  secondary_vtx.Z(), Sexa_DCA, Sexa_M);

          // fill tree
          output_tree->Fill();
        }
      }

      // (debug)
      printf("SexaquarkFinder ::\n");
    }  // end of loop over events

  }  // end of loop over trees

  output_tree->Write();
  output_file->Save();

  output_file->Close();
}

/*** Functions ***/

//________________________________________________________________________
Double_t TwoLinesDCA(TVector3 pos1, TVector3 dir1, TVector3 pos2, TVector3 dir2, TVector3 &P1, TVector3 &P2) {
  //
  // find a common vertex for two neutral particles, given their
  // this function stores the origin vertices for each particle, and returns the DCA between them
  //

  // require perpendicularity <-> lines must intersect
  if (dir1.Cross(dir2).Mag() > 0.) {
    // based on https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d
    // (Marty Cohen's answer)

    TVector3 diff = pos1 - pos2;
    Double_t determinant = dir1.Dot(dir2) * dir1.Dot(dir2) - dir1.Mag2() * dir2.Mag2();
    Double_t sol1 = (dir2.Mag2() * dir1.Dot(diff) - dir2.Dot(diff) * dir2.Dot(dir1)) / determinant;
    Double_t sol2 = (-1 * dir1.Mag2() * dir2.Dot(diff) + dir1.Dot(diff) * dir2.Dot(dir1)) / determinant;

    // endpoints of the closest line
    P1 = pos1 + sol1 * dir1;
    P2 = pos2 + sol2 * dir2;

    // the distance of the closest line
    return (diff + dir1 * sol1 - dir2 * sol2).Mag();
  }
  return -1.;
}
