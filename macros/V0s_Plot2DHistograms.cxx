#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

// -- A. Bórquez

//_____________________________________________________________________________
void V0s_Plot2DHistograms(
    TString input_filename = "/home/ceres/borquez/work/analysis/output/signal+bkg/297595/AnalysisResults_CustomV0s_*.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

    /*** Process Input ***/

    TList *list_of_trees = new TList();
    AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

    Event_tt this_event;

    /*** Prepare Output ***/

    // define v0 variables
    const Int_t n_v0_hists = 1;

    TString v0_var_name_x[n_v0_hists] = {"arm_alpha(V0)"};
    Int_t v0_var_n_bins_x[n_v0_hists] = {100};
    Float_t v0_var_min_x[n_v0_hists] = {-1.};
    Float_t v0_var_max_x[n_v0_hists] = {1.};

    TString v0_var_name_y[n_v0_hists] = {"arm_pt(V0)"};
    Int_t v0_var_n_bins_y[n_v0_hists] = {100};
    Float_t v0_var_min_y[n_v0_hists] = {0.};
    Float_t v0_var_max_y[n_v0_hists] = {1.};

    // create v0 histograms
    TH2F *hist_v0[n_v0_hists];
    TH2F *hist_v0_AL[n_v0_hists];
    TH2F *hist_v0_K0[n_v0_hists];
    TH2F *hist_v0_Sexa[n_v0_hists];
    for (Int_t vv = 0; vv < n_v0_hists; vv++) {
        hist_v0[vv] = new TH2F(v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv],           //
                               v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv],           //
                               v0_var_n_bins_x[vv], v0_var_min_x[vv], v0_var_max_x[vv],  //
                               v0_var_n_bins_y[vv], v0_var_min_y[vv], v0_var_max_y[vv]);
        hist_v0_AL[vv] = new TH2F(v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_AL",   //
                                  v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_AL",   //
                                  v0_var_n_bins_x[vv], v0_var_min_x[vv], v0_var_max_x[vv],  //
                                  v0_var_n_bins_y[vv], v0_var_min_y[vv], v0_var_max_y[vv]);
        hist_v0_K0[vv] = new TH2F(v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_K0",   //
                                  v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_K0",   //
                                  v0_var_n_bins_x[vv], v0_var_min_x[vv], v0_var_max_x[vv],  //
                                  v0_var_n_bins_y[vv], v0_var_min_y[vv], v0_var_max_y[vv]);
        hist_v0_Sexa[vv] = new TH2F(v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_Sexa",  //
                                    v0_var_name_y[vv] + "_vs_" + v0_var_name_x[vv] + "_Sexa",  //
                                    v0_var_n_bins_x[vv], v0_var_min_x[vv], v0_var_max_x[vv],   //
                                    v0_var_n_bins_y[vv], v0_var_min_y[vv], v0_var_max_y[vv]);
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
                    is_K0_signal =
                        (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]]] == 310 &&  //
                        (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]]] == 310 &&  //
                        (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] ==
                            (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]];
                    is_AL_signal =
                        (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]]] == -3122 &&  //
                        (*this_event.MC_PID)[(*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]]] == -3122 &&  //
                        (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] ==
                            (*this_event.MC_Mother)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]];
                }
                Bool_t is_Sexa_signal = (*this_event.V0_isSignal)[v0];
                Float_t V0_Pt =
                    TMath::Sqrt((*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] + (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0]);
                Float_t mass_asK0 = TMath::Sqrt(
                    (*this_event.V0_E_asK0)[v0] * (*this_event.V0_E_asK0)[v0] - (*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] -
                    (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0] - (*this_event.V0_Pz)[v0] * (*this_event.V0_Pz)[v0]);
                Float_t mass_asAL = TMath::Sqrt(
                    (*this_event.V0_E_asAL)[v0] * (*this_event.V0_E_asAL)[v0] - (*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] -
                    (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0] - (*this_event.V0_Pz)[v0] * (*this_event.V0_Pz)[v0]);

                Bool_t armenteros_condition = (*this_event.V0_ArmPt)[v0] < 0.2 * TMath::Abs((*this_event.V0_ArmAlpha)[v0]);

                // (fill hist) both V0s
                // if (((*this_event.V0_couldBeK0)[v0] && armenteros_condition) || (*this_event.V0_couldBeAL)[v0]) {
                if ((*this_event.V0_couldBeK0)[v0] || (*this_event.V0_couldBeAL)[v0]) {
                    hist_v0[0]->Fill((*this_event.V0_ArmAlpha)[v0], (*this_event.V0_ArmPt)[v0]);
                    if (is_AL_signal) {
                        hist_v0_AL[0]->Fill((*this_event.V0_ArmAlpha)[v0], (*this_event.V0_ArmPt)[v0]);
                    }
                    if (is_K0_signal) {
                        hist_v0_K0[0]->Fill((*this_event.V0_ArmAlpha)[v0], (*this_event.V0_ArmPt)[v0]);
                    }
                    if (is_Sexa_signal) {
                        hist_v0_Sexa[0]->Fill((*this_event.V0_ArmAlpha)[v0], (*this_event.V0_ArmPt)[v0]);
                    }
                }

            }  // end of loop over V0s

        }  // end of loop over events

    }  // end of loop over trees

    /* Draw */

    // set style to everything
    SetMy2DStyle();

    gStyle->SetOptStat(10);

    TString output_filename;
    TCanvas *c = new TCanvas("c", "c", 1080, 1080);
    c->SetFrameLineWidth(2);
    c->Divide(2, 2);

    for (Int_t vv = 0; vv < n_v0_hists; vv++) {
        c->cd(1);
        SetMy2DHistStyle(hist_v0[vv]);
        hist_v0[vv]->Draw("COLZ");

        c->cd(2);
        SetMy2DHistStyle(hist_v0_AL[vv]);
        hist_v0_AL[vv]->Draw("COLZ");

        c->cd(3);
        SetMy2DHistStyle(hist_v0_K0[vv]);
        hist_v0_K0[vv]->Draw("COLZ");

        c->cd(4);
        SetMy2DHistStyle(hist_v0_Sexa[vv]);
        hist_v0_Sexa[vv]->Draw("COLZ");

        output_filename = hist_v0[vv]->GetName();
        c->Print(output_dir + output_filename + ".png");
    }
}
