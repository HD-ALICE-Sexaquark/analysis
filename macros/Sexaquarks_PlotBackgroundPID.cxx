#include "include/Headers.hxx"
#include "include/Math.hxx"
#include "include/Style.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

//_____________________________________________________________________________
void Sexaquarks_PlotBackgroundPID(
    // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/*/SexaquarkResults_CustomV0s_*.root",
    TString input_filename = "/home/ceres/borquez/work/analysis/output/Sexa_CustomV0.list",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

    /*** Process Input ***/

    TString file_extension = input_filename(input_filename.First("."), input_filename.Length() - 1);

    TChain *input_chain = new TChain("");
    if (file_extension == ".list") {
        std::ifstream text_file(input_filename);
        TString root_file;
        while (text_file >> root_file) {
            input_chain->Add(root_file + "/Sexaquarks");
        }
        text_file.close();
    } else if (file_extension == ".root") {
        input_chain->Add(input_filename + "/Sexaquarks");
    } else {
        std::cerr << "ERROR: the input file should be a .list file or a .root file" << std::endl;
        return;
    }

    // define input sexaquark variables
    Bool_t is_signal;
    Float_t dca;
    Float_t mass;
    Float_t px;
    Float_t py;
    Float_t pt;
    Float_t sv_x;
    Float_t sv_y;
    Float_t radius;
    Int_t run_number;
    Int_t dir_number;
    Int_t event;
    Int_t idx_v0a;
    Int_t idx_v0b;
    Bool_t aa_is_signal;
    Bool_t bb_is_signal;
    Int_t v0a_mother_pid;
    Int_t v0a_pid;
    Int_t v0a_pos_pid;
    Int_t v0a_neg_pid;
    Int_t v0b_mother_pid;
    Int_t v0b_pid;
    Int_t v0b_pos_pid;
    Int_t v0b_neg_pid;

    input_chain->SetBranchAddress("isSignal", &is_signal);
    input_chain->SetBranchAddress("DCA", &dca);
    input_chain->SetBranchAddress("M", &mass);
    input_chain->SetBranchAddress("Px", &px);
    input_chain->SetBranchAddress("Py", &py);
    input_chain->SetBranchAddress("X", &sv_x);
    input_chain->SetBranchAddress("Y", &sv_y);

    input_chain->SetBranchAddress("RunNumber", &run_number);
    input_chain->SetBranchAddress("DirNumber", &dir_number);
    input_chain->SetBranchAddress("Event", &event);
    input_chain->SetBranchAddress("Idx_V0A", &idx_v0a);
    input_chain->SetBranchAddress("Idx_V0B", &idx_v0b);
    input_chain->SetBranchAddress("V0A_isSignal", &aa_is_signal);
    input_chain->SetBranchAddress("V0B_isSignal", &bb_is_signal);

    input_chain->SetBranchAddress("V0A_Mother_PID", &v0a_mother_pid);
    input_chain->SetBranchAddress("V0A_PID", &v0a_pid);
    input_chain->SetBranchAddress("V0A_Pos_PID", &v0a_pos_pid);
    input_chain->SetBranchAddress("V0A_Neg_PID", &v0a_neg_pid);
    input_chain->SetBranchAddress("V0B_Mother_PID", &v0b_mother_pid);
    input_chain->SetBranchAddress("V0B_PID", &v0b_pid);
    input_chain->SetBranchAddress("V0B_Pos_PID", &v0b_pos_pid);
    input_chain->SetBranchAddress("V0B_Neg_PID", &v0b_neg_pid);

    // create sets
    std::set<Int_t> possible_v0_mothers_pid;
    std::set<Int_t> possible_v0s_pid;
    std::set<Int_t> possible_daughters_pid;

    /*** Prepare Output ***/

    // define sexaquark variables
    const Int_t n_sexa_vars = 4;
    TString sexa_var_name[n_sexa_vars] = {"m(AL,K0)", "dca(AL,K0)", "radius(S)", "pt(S)"};
    TString sexa_axis_name[n_sexa_vars] = {"Invariant Mass #bar{#Lambda}K^{0}_{S} (GeV/c^{2})", "dca(AL,K0)", "radius(S)", "pt(S)"};
    Int_t sexa_var_n_bins[n_sexa_vars] = {200, 50, 200, 200};
    Float_t sexa_var_min_x[n_sexa_vars] = {0., 0., 0., 0.};
    Float_t sexa_var_max_x[n_sexa_vars] = {10., 10., 200., 10.};

    // create sexa histograms
    TH1F *hist_sexa_true_bkg[n_sexa_vars];
    for (Int_t ss = 0; ss < n_sexa_vars; ss++) {
        hist_sexa_true_bkg[ss] = new TH1F(sexa_var_name[ss] + "_true-bkg_PID", sexa_var_name[ss] + "_true-bkg_PID",  //
                                          sexa_var_n_bins[ss], sexa_var_min_x[ss], sexa_var_max_x[ss]);
    }

    // in case the output dir doesn't exist, create it
    system("mkdir -p " + output_dir);

    // (debug)
    printf("Sexaquarks_PlotBackgroundPID :: RN     D E   V0A V0B PID(MomV0A V0A    V0APos V0ANeg MomV0B V0B    V0BPos V0BNeg)\n");

    // (counters)
    Int_t n_v0s_same_mother = 0;
    TString same_mother_str;

    // loop over candidates
    for (Int_t candidate = 0; candidate < input_chain->GetEntries(); candidate++) {

        input_chain->GetEntry(candidate);

        pt = TMath::Sqrt(px * px + py * py);
        radius = TMath::Sqrt(sv_x * sv_x + sv_y * sv_y);

        // (cut) to consider only true bkg
        if (!aa_is_signal && !bb_is_signal) {
            // (plot) before cut
            hist_sexa_true_bkg[1]->Fill(dca);
            // (cut) on dca
            if (0. < dca && dca < 0.5) {
                // (plot) before cut
                hist_sexa_true_bkg[2]->Fill(radius);
                // (cut) on radius
                if (radius > 50.) {
                    // (plot) before cut
                    hist_sexa_true_bkg[3]->Fill(pt);
                    // (cut) on pt
                    if (pt > 1.) {
                        hist_sexa_true_bkg[0]->Fill(mass);

                        // count
                        if (v0a_mother_pid == v0b_mother_pid && v0b_mother_pid != 0) {
                            n_v0s_same_mother++;
                            same_mother_str = "-> HERE!!";
                        } else {
                            same_mother_str = "";
                        }

                        // (debug)
                        printf("Sexaquarks_PlotBackgroundPID :: %6i %1i %3i %3i %3i     %6i %6i %6i %6i %6i %6i %6i %6i %s\n",  //
                               run_number, dir_number, event, idx_v0a, idx_v0b,                                                 //
                               v0a_mother_pid, v0a_pid, v0a_pos_pid, v0a_neg_pid,                                               //
                               v0b_mother_pid, v0b_pid, v0b_pos_pid, v0b_neg_pid, same_mother_str.Data());

                        // insert to V0 mothers sets
                        possible_v0_mothers_pid.insert(v0a_mother_pid);
                        possible_v0_mothers_pid.insert(v0b_mother_pid);
                        // insert to V0 sets
                        possible_v0s_pid.insert(v0a_pid);
                        possible_v0s_pid.insert(v0b_pid);
                        // insert to daughters sets
                        possible_daughters_pid.insert(v0a_pos_pid);
                        possible_daughters_pid.insert(v0a_neg_pid);
                        possible_daughters_pid.insert(v0b_pos_pid);
                        possible_daughters_pid.insert(v0b_neg_pid);
                    }
                }
            }
        }
    }  // end of loop over candidates

    /* Print info */

    printf("Sexaquarks_PlotBackgroundPID ::\n");
    printf("Sexaquarks_PlotBackgroundPID :: SUMMARY\n");
    printf("Sexaquarks_PlotBackgroundPID :: =======\n");

    printf("Sexaquarks_PlotBackgroundPID :: Number of true-bkg candidates two V0s with the same mother: %i\n", n_v0s_same_mother);

    printf("Sexaquarks_PlotBackgroundPID :: V0 Mothers' PID: ");
    for (auto it = possible_v0_mothers_pid.begin(); it != possible_v0_mothers_pid.end(); ++it) {
        Int_t element = *it;
        printf("%i, ", element);
    }
    printf("\n");

    printf("Sexaquarks_PlotBackgroundPID :: V0s' PID: ");
    for (auto it = possible_v0s_pid.begin(); it != possible_v0s_pid.end(); ++it) {
        Int_t element = *it;
        printf("%i, ", element);
    }
    printf("\n");

    printf("Sexaquarks_PlotBackgroundPID :: Daughters' PID: ");
    for (auto it = possible_daughters_pid.begin(); it != possible_daughters_pid.end(); ++it) {
        Int_t element = *it;
        printf("%i, ", element);
    }
    printf("\n");

    /* Draw */

    // set style to everything
    SetMyStyle();
    gStyle->SetOptStat(0);

    TString output_filename;
    TCanvas *c = new TCanvas("c", "c", 1080, 1080);

    c->SetLogy(1);

    // loop over sexaquark variables
    for (Int_t ss = 0; ss < n_sexa_vars; ss++) {

        if (ss == 0) c->SetLogy(0);

        SetMyHistStyle(hist_sexa_true_bkg[ss]);
        hist_sexa_true_bkg[ss]->SetFillStyle(0);
        hist_sexa_true_bkg[ss]->SetLineWidth(3);
        hist_sexa_true_bkg[ss]->SetLineColor(myBlack);
        if (ss) {
            hist_sexa_true_bkg[ss]->SetMinimum(1);
            hist_sexa_true_bkg[ss]->SetMaximum(10 * hist_sexa_true_bkg[ss]->GetMaximum());  // adjusted for one file
        }
        hist_sexa_true_bkg[ss]->GetYaxis()->SetTitle("Counts");
        hist_sexa_true_bkg[ss]->GetXaxis()->SetTitle(sexa_axis_name[ss]);
        hist_sexa_true_bkg[ss]->Draw();

        /*
            hist_sexa_true_bkg[ss]->SetLineWidth(0);
            hist_sexa_true_bkg[ss]->SetFillStyle(kSolid);
            hist_sexa_true_bkg[ss]->SetFillColorAlpha(highlight_color, 0.5);
            hist_sexa_true_bkg[ss]->SetLineWidth(0);
            hist_sexa_true_bkg[ss]->Draw("SAME");
            hist_sexa_true_bkg[ss]->Draw("AXIS SAME");
         */

        DrawText(Form("N(true-bkg) = %i", (Int_t)hist_sexa_true_bkg[ss]->GetEntries()), myBlack, 0.55, 0.8);
        // DrawText(Form("N(%s) = %i", highlight_opt.Data(), (Int_t)hist_sexa_true_bkg[ss]->GetEntries()), highlight_color, 0.55, 0.75);

        output_filename = hist_sexa_true_bkg[ss]->GetName();
        c->Print(output_dir + output_filename + ".png");

        if (ss == 0) c->SetLogy(1);
    }  // end of loop over sexaquark variables

    printf("Sexaquarks_PlotBackgroundPID ::\n");
}
