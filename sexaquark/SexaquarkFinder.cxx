#include "macros/include/Headers.hxx"
#include "macros/include/Math.hxx"
#include "macros/include/TreeFunctions.hxx"
#include "macros/include/Utilities.hxx"

#include "SexaquarkFinder.h"

//________________________________________________________________________
void SexaquarkFinder(  //
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_CustomV0s_*.root",
    TString output_filename = "output/signal+bkg/297595/SexaquarkResults_CustomV0s.root") {
  // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_OfficialV0s_*.root",
  // TString output_filename = "output/signal+bkg/297595/SexaquarkResults_OfficialV0s.root") {
  //
  // Just another macro, but more important than the others...
  //

  /*** Process Input ***/

  TList *list_of_trees = new TList();
  AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

  Event_tt input_event;

  /*** Prepare Output ***/

  TFile *output_file = new TFile(output_filename, "RECREATE");
  TTree *output_tree = new TTree("Sexaquarks", "Sexaquark Candidates");

  Sexa_tt this_sexa;
  SetSexaBranches(output_tree, this_sexa);

  /*** Start ***/

  Int_t run_number = ((TString)input_filename(input_filename.Index("/2") + 1, 6)).Atoi();
  Int_t dir_number = ((TString)input_filename(input_filename.Index("_00") + 1, 3)).Atoi();

  // (debug)
  printf("SexaquarkFinder :: Processing run number %i\n", run_number);

  // loop over collected trees
  TListIter *list_of_trees_it = new TListIter(list_of_trees);
  while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

    LoadBranches(this_tree, input_event);

    // (debug)
    printf("SexaquarkFinder :: Processing dir %i\n", dir_number);

    // loop over events
    for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

      this_tree->GetEntry(event);

      // (debug)
      printf("SexaquarkFinder :: Processing event #%i\n", event);

      // define variables & objects
      TVector3 tv3_pos_K0;
      TVector3 tv3_mom_K0;
      TVector3 tv3_pos_AL;
      TVector3 tv3_mom_AL;

      TVector3 pca_A;  // point of closest approach
      TVector3 pca_B;
      TVector3 secondary_vtx;
      Double_t half_distance;

      Double_t energy_K0;
      Double_t energy_AL;
      Double_t mass_K0;
      Double_t mass_AL;

      // (debug)
      printf("SexaquarkFinder :: ANTI-SEXAQUARK CANDIDATES\n");
      printf("SexaquarkFinder :: \n");
      printf("SexaquarkFinder :: Idx_V0A Idx_V0B SIG_V0A SIG_V0B SIGSEXA    SV_X    SV_Y    SV_Z     DCA       M\n");

      // combinatorics
      for (Int_t aa = 0; aa < input_event.N_V0s - 1; aa++) {
        for (Int_t bb = aa + 1; bb < input_event.N_V0s; bb++) {

          // (assign) identifiers
          this_sexa.RunNumber = run_number;
          this_sexa.DirNumber = dir_number;
          this_sexa.Event = event;

          // (assign) both V0s must be signal for the candidate to be signal
          // PENDING: is this right? check comment below
          this_sexa.Sexa_isSignal = (*input_event.V0_isSignal)[aa] && (*input_event.V0_isSignal)[bb];

          // (cut) check that all four granddaughters are different particles
          if ((*input_event.Idx_Pos)[aa] == (*input_event.Idx_Pos)[bb] ||  //
              (*input_event.Idx_Pos)[aa] == (*input_event.Idx_Neg)[bb] ||  //
              (*input_event.Idx_Neg)[aa] == (*input_event.Idx_Pos)[bb] ||  //
              (*input_event.Idx_Neg)[aa] == (*input_event.Idx_Neg)[bb]) {
            // PENDING: in reality, didn't lose any signal it...
            //          the fact that the 4 particles are signal, doesn't mean it's a proper signal candidate
            //          because we also need to consider if they come from the same/different V0 when corresponds
            printf("SexaquarkFinder :: At least one of the 4 daughters was repeated!\n");
            continue;
          }

          // (cut) check that one of the V0s is indeed a K0 candidate and the other an anti-lambda candidate
          // (borquez comment: can be extended to allow ambiguity)
          if ((*input_event.V0_couldBeK0)[aa] && (*input_event.V0_couldBeAL)[bb]) {
            // assign V0_A to the K0
            tv3_pos_K0.SetXYZ((*input_event.V0_X)[aa], (*input_event.V0_Y)[aa], (*input_event.V0_Z)[aa]);
            tv3_mom_K0.SetXYZ((*input_event.V0_Px)[aa], (*input_event.V0_Py)[aa], (*input_event.V0_Pz)[aa]);
            energy_K0 = (*input_event.V0_E_asK0)[aa];
            // assign V0_B to the AL
            tv3_pos_AL.SetXYZ((*input_event.V0_X)[bb], (*input_event.V0_Y)[bb], (*input_event.V0_Z)[bb]);
            tv3_mom_AL.SetXYZ((*input_event.V0_Px)[bb], (*input_event.V0_Py)[bb], (*input_event.V0_Pz)[bb]);
            energy_AL = (*input_event.V0_E_asAL)[bb];
          } else if ((*input_event.V0_couldBeAL)[aa] && (*input_event.V0_couldBeK0)[bb]) {
            // assign V0_A to the K0
            tv3_pos_K0.SetXYZ((*input_event.V0_X)[bb], (*input_event.V0_Y)[bb], (*input_event.V0_Z)[bb]);
            tv3_mom_K0.SetXYZ((*input_event.V0_Px)[bb], (*input_event.V0_Py)[bb], (*input_event.V0_Pz)[bb]);
            energy_K0 = (*input_event.V0_E_asK0)[bb];
            // assign V0_B to the AL
            tv3_pos_AL.SetXYZ((*input_event.V0_X)[aa], (*input_event.V0_Y)[aa], (*input_event.V0_Z)[aa]);
            tv3_mom_AL.SetXYZ((*input_event.V0_Px)[aa], (*input_event.V0_Py)[aa], (*input_event.V0_Pz)[aa]);
            energy_AL = (*input_event.V0_E_asAL)[aa];
          } else {
            if (this_sexa.Sexa_isSignal) printf("SexaquarkFinder :: Lost a signal candidate because because of rec. PID\n");
            continue;
          }

          // (assign) from both momentum and coordinates of V0s,
          // extrapolate them to a common vertex,
          // if there's none, calculate DCA and extract PCAs
          this_sexa.Sexa_DCA = TwoLinesDCA(tv3_pos_K0, tv3_mom_K0, tv3_pos_AL, tv3_mom_AL, pca_A, pca_B);

          // (cut) when the two lines cannot intersect in the correct direction
          if (this_sexa.Sexa_DCA < 0.) {
            continue;
          }

          // calculate the secondary_vtx point between PCAs
          half_distance = (pca_B - pca_A).Mag() / 2.;
          secondary_vtx = pca_A + half_distance * (pca_B - pca_A).Unit();

          // (assign) rest of variables
          this_sexa.Idx_V0A = aa;
          this_sexa.Idx_V0B = bb;
          this_sexa.Sexa_E = energy_K0 + energy_AL - 0.939565;  // assuming struck neutron at rest
          this_sexa.Sexa_Px = (tv3_mom_K0 + tv3_mom_AL).X();
          this_sexa.Sexa_Py = (tv3_mom_K0 + tv3_mom_AL).Y();
          this_sexa.Sexa_Pz = (tv3_mom_K0 + tv3_mom_AL).Z();

          // (cut) remove off-shell sexaquarks
          Float_t pre_mass = this_sexa.Sexa_E * this_sexa.Sexa_E - this_sexa.Sexa_Px * this_sexa.Sexa_Px -
                             this_sexa.Sexa_Py * this_sexa.Sexa_Py - this_sexa.Sexa_Pz * this_sexa.Sexa_Pz;
          if (pre_mass < 0.) {
            if (this_sexa.Sexa_isSignal) printf("SexaquarkFinder :: Lost a signal candidate because of reconstructed off-shell mass\n");
            continue;
          }

          this_sexa.Sexa_M = TMath::Sqrt(pre_mass);
          this_sexa.Sexa_X = secondary_vtx.X();
          this_sexa.Sexa_Y = secondary_vtx.Y();
          this_sexa.Sexa_Z = secondary_vtx.Z();

          // (assign) first V0: "V0A"
          this_sexa.V0A_Idx_Pos = (*input_event.Idx_Pos)[aa];
          this_sexa.V0A_Idx_Neg = (*input_event.Idx_Neg)[aa];
          this_sexa.V0A_Px = (*input_event.V0_Px)[aa];
          this_sexa.V0A_Py = (*input_event.V0_Py)[aa];
          this_sexa.V0A_Pz = (*input_event.V0_Pz)[aa];
          this_sexa.V0A_X = (*input_event.V0_X)[aa];
          this_sexa.V0A_Y = (*input_event.V0_Y)[aa];
          this_sexa.V0A_Z = (*input_event.V0_Z)[aa];
          this_sexa.V0A_Pos_Px = (*input_event.Pos_Px)[aa];
          this_sexa.V0A_Pos_Py = (*input_event.Pos_Py)[aa];
          this_sexa.V0A_Pos_Pz = (*input_event.Pos_Pz)[aa];
          this_sexa.V0A_Neg_Px = (*input_event.Neg_Px)[aa];
          this_sexa.V0A_Neg_Py = (*input_event.Neg_Py)[aa];
          this_sexa.V0A_Neg_Pz = (*input_event.Neg_Pz)[aa];
          this_sexa.V0A_isSignal = (*input_event.V0_isSignal)[aa];
          this_sexa.V0A_E_asK0 = (*input_event.V0_E_asK0)[aa];
          this_sexa.V0A_E_asAL = (*input_event.V0_E_asAL)[aa];
          this_sexa.V0A_couldBeK0 = (*input_event.V0_couldBeK0)[aa];
          this_sexa.V0A_couldBeAL = (*input_event.V0_couldBeAL)[aa];
          this_sexa.V0A_Chi2 = (*input_event.V0_Chi2)[aa];
          this_sexa.V0A_DCA_Daughters = (*input_event.V0_DCA_Daughters)[aa];
          this_sexa.V0A_IP_wrtPV = (*input_event.V0_IP_wrtPV)[aa];
          this_sexa.V0A_CPA_wrtPV = (*input_event.V0_CPA_wrtPV)[aa];
          this_sexa.V0A_ArmAlpha = (*input_event.V0_ArmAlpha)[aa];
          this_sexa.V0A_ArmPt = (*input_event.V0_ArmPt)[aa];
          this_sexa.V0A_DecayLength = (*input_event.V0_DecayLength)[aa];

          this_sexa.V0A_Pos_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Neg_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Neg)[aa]];
          this_sexa.V0A_Pos_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Pos_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Pos_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Neg_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Neg)[aa]];
          this_sexa.V0A_Neg_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Neg)[aa]];
          this_sexa.V0A_Neg_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Neg)[aa]];
          this_sexa.V0A_Pos_NSigmaPion = (*input_event.Rec_NSigmaPion)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Pos_NSigmaProton = (*input_event.Rec_NSigmaProton)[(*input_event.Idx_Pos)[aa]];
          this_sexa.V0A_Neg_NSigmaPion = (*input_event.Rec_NSigmaPion)[(*input_event.Idx_Neg)[aa]];
          this_sexa.V0A_Neg_NSigmaProton = (*input_event.Rec_NSigmaProton)[(*input_event.Idx_Neg)[aa]];

          this_sexa.V0A_Mother_PID =
              (*input_event.MC_PID_GrandMother)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[aa]]];  // pos. or neg. is the same
          this_sexa.V0A_PID =
              (*input_event.MC_PID_Mother)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[aa]]];       // pos. or neg. is the same
          this_sexa.V0A_Pos_PID = (*input_event.MC_PID)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[aa]]];
          this_sexa.V0A_Neg_PID = (*input_event.MC_PID)[(*input_event.Idx_True)[(*input_event.Idx_Neg)[aa]]];

          // (assign) second V0: "V0B"
          this_sexa.V0B_Idx_Pos = (*input_event.Idx_Pos)[bb];
          this_sexa.V0B_Idx_Neg = (*input_event.Idx_Neg)[bb];
          this_sexa.V0B_Px = (*input_event.V0_Px)[bb];
          this_sexa.V0B_Py = (*input_event.V0_Py)[bb];
          this_sexa.V0B_Pz = (*input_event.V0_Pz)[bb];
          this_sexa.V0B_X = (*input_event.V0_X)[bb];
          this_sexa.V0B_Y = (*input_event.V0_Y)[bb];
          this_sexa.V0B_Z = (*input_event.V0_Z)[bb];
          this_sexa.V0B_Pos_Px = (*input_event.Pos_Px)[bb];
          this_sexa.V0B_Pos_Py = (*input_event.Pos_Py)[bb];
          this_sexa.V0B_Pos_Pz = (*input_event.Pos_Pz)[bb];
          this_sexa.V0B_Neg_Px = (*input_event.Neg_Px)[bb];
          this_sexa.V0B_Neg_Py = (*input_event.Neg_Py)[bb];
          this_sexa.V0B_Neg_Pz = (*input_event.Neg_Pz)[bb];
          this_sexa.V0B_isSignal = (*input_event.V0_isSignal)[bb];
          this_sexa.V0B_E_asK0 = (*input_event.V0_E_asK0)[bb];
          this_sexa.V0B_E_asAL = (*input_event.V0_E_asAL)[bb];
          this_sexa.V0B_couldBeK0 = (*input_event.V0_couldBeK0)[bb];
          this_sexa.V0B_couldBeAL = (*input_event.V0_couldBeAL)[bb];
          this_sexa.V0B_Chi2 = (*input_event.V0_Chi2)[bb];
          this_sexa.V0B_DCA_Daughters = (*input_event.V0_DCA_Daughters)[bb];
          this_sexa.V0B_IP_wrtPV = (*input_event.V0_IP_wrtPV)[bb];
          this_sexa.V0B_CPA_wrtPV = (*input_event.V0_CPA_wrtPV)[bb];
          this_sexa.V0B_ArmAlpha = (*input_event.V0_ArmAlpha)[bb];
          this_sexa.V0B_ArmPt = (*input_event.V0_ArmPt)[bb];
          this_sexa.V0B_DecayLength = (*input_event.V0_DecayLength)[bb];

          this_sexa.V0B_Pos_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Neg_isDuplicate = (*input_event.Rec_isDuplicate)[(*input_event.Idx_Neg)[bb]];
          this_sexa.V0B_Pos_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Pos_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Pos_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Neg_Rec_Px = (*input_event.Rec_Px)[(*input_event.Idx_Neg)[bb]];
          this_sexa.V0B_Neg_Rec_Py = (*input_event.Rec_Py)[(*input_event.Idx_Neg)[bb]];
          this_sexa.V0B_Neg_Rec_Pz = (*input_event.Rec_Pz)[(*input_event.Idx_Neg)[bb]];
          this_sexa.V0B_Pos_NSigmaPion = (*input_event.Rec_NSigmaPion)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Pos_NSigmaProton = (*input_event.Rec_NSigmaProton)[(*input_event.Idx_Pos)[bb]];
          this_sexa.V0B_Neg_NSigmaPion = (*input_event.Rec_NSigmaPion)[(*input_event.Idx_Neg)[bb]];
          this_sexa.V0B_Neg_NSigmaProton = (*input_event.Rec_NSigmaProton)[(*input_event.Idx_Neg)[bb]];

          this_sexa.V0B_Mother_PID =
              (*input_event.MC_PID_GrandMother)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[bb]]];  // pos. or neg. is the same
          this_sexa.V0B_PID =
              (*input_event.MC_PID_Mother)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[bb]]];       // pos. or neg. is the same
          this_sexa.V0B_Pos_PID = (*input_event.MC_PID)[(*input_event.Idx_True)[(*input_event.Idx_Pos)[bb]]];
          this_sexa.V0B_Neg_PID = (*input_event.MC_PID)[(*input_event.Idx_True)[(*input_event.Idx_Neg)[bb]]];

          // (debug)
          printf("SexaquarkFinder :: %7i %7i %7i %7i %7i %7.3f %7.3f %7.3f %7.3f %7.3f\n", aa, bb, this_sexa.V0A_isSignal,
                 this_sexa.V0B_isSignal, this_sexa.Sexa_isSignal, secondary_vtx.X(), secondary_vtx.Y(), secondary_vtx.Z(),
                 this_sexa.Sexa_DCA, this_sexa.Sexa_M);

          // fill tree
          output_tree->Fill();
        }
      }

      // (debug)
      printf("SexaquarkFinder ::\n");
    }  // end of loop over events

    // increase dir number
    dir_number++;
  }  // end of loop over trees

  output_tree->Write();
  output_file->Save();

  output_file->Close();
}
