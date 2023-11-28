#include "include/Headers.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

// Calculate efficiency of Custom V0 Finder
// -- A. Bórquez
// -- 03.Jul.2023

// ```
// V0s Calculate Custom V0 Finder Efficiency
// =========================================
// >> n_mc_signal_nice_reco_kaons       = 117939
// >> n_mc_signal_nice_reco_antilambdas = 109036
// >> n_rec_signal_v0s                  = 57398
// => n_rec_signal_v0s / (n_mc_signal_nice_reco_kaons + n_mc_signal_nice_reco_antilambdas) = 0.25288248
// ```

//_____________________________________________________________________________
void V0s_CalcCustomV0FinderEfficiency(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/*/AnalysisResults_CustomV0s_*.root") {
    // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/246052/AnalysisResults_CustomV0s_000.root") {

    /*** Process Input ***/

    TList *list_of_trees = new TList();
    AddTreesToList(list_of_trees, input_filename + "/Trees/Events");  // PENDING: consumes too much memory!!

    Event_tt this_event;

    // declare counters
    Int_t n_mc_signal_nice_reco_kaons = 0;
    Int_t n_mc_signal_nice_reco_antilambdas = 0;
    Int_t n_rec_signal_v0s = 0;

    // aux variables
    Int_t mc_index_dau1, mc_index_dau2;
    Int_t mc_pdg_dau1, mc_pdg_dau2;

    Bool_t cond_mc_nice_channel_kaon;
    Bool_t cond_mc_nice_channel_antilambda;

    std::set<Int_t> mc_indices_reco;  // indices of MC particles that were reconstructed
    Bool_t cond_mc_were_reco;

    // loop over collected trees
    TListIter *list_of_trees_it = new TListIter(list_of_trees);
    while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

        // load branches
        LoadBranches(this_tree, this_event);

        // loop over events
        for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

            this_tree->GetEntry(event);

            /* Reconstructed Particles */

            // collect mc gen. particles that were reconstructed
            // loop over mc rec. particles
            for (Int_t rec = 0; rec < this_event.N_MCRec; rec++) {
                mc_indices_reco.insert((*this_event.Idx_True)[rec]);
            }  // end of loop over rec.

            /* MC Generated Particles */

            // loop over mc generated particles
            for (Int_t mc = 0; mc < this_event.N_MCGen; mc++) {

                // (cut) particle must be signal
                if (!(*this_event.MC_isSignal)[mc]) {
                    continue;
                }

                // (cut) must have decayed into relevant particles
                mc_index_dau1 = (*this_event.MC_FirstDau)[mc];
                mc_index_dau2 = (*this_event.MC_LastDau)[mc];
                if (mc_index_dau1 <= 0 || mc_index_dau2 <= 0) {
                    continue;
                }

                // (cut) both daughters (charged particles) must have been reconstructed
                cond_mc_were_reco = mc_indices_reco.count(mc_index_dau1) && mc_indices_reco.count(mc_index_dau2);
                if (!cond_mc_were_reco) {
                    continue;
                }

                // (cut) nice channels
                // . K0S -> pi+ pi-
                // . AL  -> anti-p pi+
                mc_pdg_dau1 = (*this_event.MC_PID)[mc_index_dau1];
                mc_pdg_dau2 = (*this_event.MC_PID)[mc_index_dau2];

                cond_mc_nice_channel_kaon = (mc_pdg_dau1 == 211 && mc_pdg_dau2 == -211) || (mc_pdg_dau1 == -211 && mc_pdg_dau2 == 211);
                cond_mc_nice_channel_antilambda =
                    (mc_pdg_dau1 == -2212 && mc_pdg_dau2 == 211) || (mc_pdg_dau1 == 211 && mc_pdg_dau2 == -2212);

                if ((*this_event.MC_PID)[mc] == 310 && cond_mc_nice_channel_kaon) n_mc_signal_nice_reco_kaons++;
                if ((*this_event.MC_PID)[mc] == -3122 && cond_mc_nice_channel_antilambda) n_mc_signal_nice_reco_antilambdas++;
            }  // end of loop over mc

            /* Found V0s */

            // loop over V0s
            for (Int_t v0 = 0; v0 < this_event.N_V0s; v0++) {

                // (cut) V0 must be signal
                if (!(*this_event.V0_isSignal)[v0]) {
                    continue;
                }

                if ((*this_event.V0_couldBeAL)[v0] || (*this_event.V0_couldBeK0)[v0]) n_rec_signal_v0s++;
            }  // end of loop over V0s

            // reset variables
        }  // end of loop over events

    }  // end of loop over trees

    // print results
    printf("V0s Calculate Custom V0 Finder Efficiency\n");
    printf("=========================================\n");
    printf(">> n_mc_signal_nice_reco_kaons       = %i\n", n_mc_signal_nice_reco_kaons);
    printf(">> n_mc_signal_nice_reco_antilambdas = %i\n", n_mc_signal_nice_reco_antilambdas);
    printf(">> n_rec_signal_v0s                  = %i\n", n_rec_signal_v0s);
}
