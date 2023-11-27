#include "include/Headers.hxx"
#include "include/Math.hxx"
#include "include/Style.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

//_____________________________________________________________________________
void Sexaquarks_Plot1DHistograms(
    // TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/*/SexaquarkResults_CustomV0s_*.root",
    TString input_filename = "/home/ceres/borquez/work/analysis/output/Sexa_CustomV0.list",
    TString highlight_opt = "partial-bkg",  // signal | bkg | partial-bkg | true-bkg
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
    Bool_t aa_is_signal;
    Bool_t bb_is_signal;
    Float_t dca;
    Float_t mass;
    Float_t px;
    Float_t py;
    Float_t pt;
    Float_t sv_x;
    Float_t sv_y;
    Float_t radius;
    input_chain->SetBranchAddress("V0A_isSignal", &aa_is_signal);
    input_chain->SetBranchAddress("V0B_isSignal", &bb_is_signal);
    input_chain->SetBranchAddress("isSignal", &is_signal);
    input_chain->SetBranchAddress("DCA", &dca);
    input_chain->SetBranchAddress("M", &mass);
    input_chain->SetBranchAddress("Px", &px);
    input_chain->SetBranchAddress("Py", &py);
    input_chain->SetBranchAddress("X", &sv_x);
    input_chain->SetBranchAddress("Y", &sv_y);

    // define highlight condition
    Bool_t highlight_condition;
    Bool_t only_signal = highlight_opt == "signal";
    Bool_t only_bkg = highlight_opt == "bkg";
    Bool_t only_true_bkg = highlight_opt == "true-bkg";
    Bool_t only_partial_bkg = highlight_opt == "partial-bkg";

    Color_t highlight_color;
    if (highlight_opt == "signal") {
        highlight_color = myGreen;
    } else if (highlight_opt == "bkg") {
        highlight_color = myBlue;
    } else if (highlight_opt == "true-bkg") {
        highlight_color = myViolet;
    } else if (highlight_opt == "partial-bkg") {
        highlight_color = myOrange;
    }

    /*** Prepare Output ***/

    // define sexaquark variables
    const Int_t n_sexa_vars = 4;
    TString sexa_var_name[n_sexa_vars] = {"m(AL,K0)", "dca(AL,K0)", "radius(S)", "pt(S)"};
    TString sexa_axis_name[n_sexa_vars] = {"Invariant Mass #bar{#Lambda}K^{0}_{S} (GeV/c^{2})", "dca(AL,K0)", "radius(S)", "pt(S)"};
    Int_t sexa_var_n_bins[n_sexa_vars] = {200, 50, 200, 200};
    Float_t sexa_var_min_x[n_sexa_vars] = {0., 0., 0., 0.};
    Float_t sexa_var_max_x[n_sexa_vars] = {10., 10., 200., 10.};

    // create sexa histograms
    TH1F *hist_sexa[n_sexa_vars];
    TH1F *hist_sexa_highlight[n_sexa_vars];
    for (Int_t ss = 0; ss < n_sexa_vars; ss++) {
        hist_sexa[ss] = new TH1F(sexa_var_name[ss], sexa_var_name[ss],  //
                                 sexa_var_n_bins[ss], sexa_var_min_x[ss], sexa_var_max_x[ss]);
        hist_sexa_highlight[ss] = new TH1F(sexa_var_name[ss] + "_" + highlight_opt, sexa_var_name[ss] + "_" + highlight_opt,  //
                                           sexa_var_n_bins[ss], sexa_var_min_x[ss], sexa_var_max_x[ss]);
    }

    // in case the output dir doesn't exist, create it
    system("mkdir -p " + output_dir);

    // loop over candidates
    for (Int_t candidate = 0; candidate < input_chain->GetEntries(); candidate++) {

        input_chain->GetEntry(candidate);

        pt = TMath::Sqrt(px * px + py * py);
        radius = TMath::Sqrt(sv_x * sv_x + sv_y * sv_y);

        highlight_condition = (only_signal && is_signal) ||                         //
                              (only_bkg && !is_signal) ||                           //
                              (only_true_bkg && !aa_is_signal && !bb_is_signal) ||  //
                              (only_partial_bkg && ((!aa_is_signal && bb_is_signal) || (aa_is_signal && !bb_is_signal)));

        hist_sexa[1]->Fill(dca);
        if (highlight_condition) hist_sexa_highlight[1]->Fill(dca);

        // (cut) on dca
        if (dca < 0.5) {

            hist_sexa[2]->Fill(radius);
            if (highlight_condition) hist_sexa_highlight[2]->Fill(radius);

            // (cut) on radius
            if (radius > 50.) {

                hist_sexa[3]->Fill(pt);
                if (highlight_condition) hist_sexa_highlight[3]->Fill(pt);

                // (cut) on pt
                if (pt > 1.) {
                    hist_sexa[0]->Fill(mass);
                    if (highlight_condition) hist_sexa_highlight[0]->Fill(mass);
                }
            }
        }
    }  // end of loop over candidates

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

        SetMyHistStyle(hist_sexa[ss]);
        hist_sexa[ss]->SetFillStyle(0);
        hist_sexa[ss]->SetLineWidth(3);
        hist_sexa[ss]->SetLineColor(myBlack);
        if (ss) {
            hist_sexa[ss]->SetMinimum(1);
            hist_sexa[ss]->SetMaximum(10 * hist_sexa[ss]->GetMaximum());  // adjusted for one file
        }
        hist_sexa[ss]->GetYaxis()->SetTitle("Counts");
        hist_sexa[ss]->GetXaxis()->SetTitle(sexa_axis_name[ss]);
        hist_sexa[ss]->Draw();

        hist_sexa_highlight[ss]->SetLineWidth(0);
        hist_sexa_highlight[ss]->SetFillStyle(kSolid);
        hist_sexa_highlight[ss]->SetFillColorAlpha(highlight_color, 0.5);
        hist_sexa_highlight[ss]->SetLineWidth(0);
        hist_sexa_highlight[ss]->Draw("SAME");
        hist_sexa_highlight[ss]->Draw("AXIS SAME");

        DrawText(Form("N(signal+bkg) = %i", (Int_t)hist_sexa[ss]->GetEntries()), myBlack, 0.55, 0.8);
        DrawText(Form("N(%s) = %i", highlight_opt.Data(), (Int_t)hist_sexa_highlight[ss]->GetEntries()), highlight_color, 0.55, 0.75);

        output_filename = hist_sexa_highlight[ss]->GetName();
        c->Print(output_dir + output_filename + ".png");

        if (ss == 0) c->SetLogy(1);
    }  // end of loop over sexaquark variables
}
