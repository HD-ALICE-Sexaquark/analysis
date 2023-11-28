#include "include/Headers.hxx"
#include "include/Math.hxx"
#include "include/Style.hxx"
#include "include/Utilities.hxx"

// I don't know why do I keep on trying...
// -- A. Bórquez

//_____________________________________________________________________________
void Sexaquarks_CombineAndPlot(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_CustomV0s_000.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

    // input
    TList *list_of_trees = new TList();
    AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

    // in case the output dir doesn't exist, create it
    system("mkdir -p " + output_dir);

    // [based on "Event_tt" from "AliAnalysisTaskSexaquark.h"]
    /* define MC variables */                    //
    std::vector<Int_t> *MC_PID = 0;              // PDG code
    std::vector<Int_t> *MC_Mother = 0;           // index of mother
    std::vector<Bool_t> *MC_isSignal = 0;        // index of mother
    /* define track variables  */                //
    Int_t N_MCRec;                               //
    std::vector<Int_t> *Idx_True = 0;            // index of true MC particle
    std::vector<Float_t> *Rec_Px = 0;            //
    std::vector<Float_t> *Rec_Py = 0;            //
    std::vector<Float_t> *Rec_Pz = 0;            //
    std::vector<Short_t> *Rec_Charge = 0;        //
    std::vector<Float_t> *Rec_NSigmaPion = 0;    //
    std::vector<Float_t> *Rec_NSigmaProton = 0;  //
    std::vector<Bool_t> *Rec_isDuplicate = 0;    // kTRUE if track is a duplicate, kFALSE if not
    std::vector<Bool_t> *Rec_isSignal = 0;       //
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
    std::vector<Float_t> *V0_DCA_wrtPV = 0;      // DCA wrt PV

    // define sexa variables
    const Int_t n_sexa_vars = 10;
    TString sexa_var_name[n_sexa_vars] = {"dca(K0,AL)",  "dca(K0,S_v)",   //
                                          "dca(AL,S_v)", "dca(S,P_v)",    //
                                          "cpa(K0,S_v)", "cpa(AL,S_v)",   //
                                          "cpa(S,P_v)",  "pt(S)",         //
                                          "m(S)",        "radius(S_v)"};  // PENDING: "dist(K0_v,S_v)" "dist(AL_v,S_v)"
    Int_t sexa_var_n_bins[n_sexa_vars] = {100, 100,                       //
                                          100, 100,                       //
                                          100, 100,                       //
                                          100, 100,                       //
                                          100, 100};
    Float_t sexa_var_min_x[n_sexa_vars] = {0.,  0.,   //
                                           0.,  0.,   //
                                           -1., -1.,  //
                                           -1., 0.,   //
                                           0.,  0.};
    Float_t sexa_var_max_x[n_sexa_vars] = {10., 10.,  //
                                           10., 10.,  //
                                           1.,  1.,   //
                                           1.,  10.,  //
                                           5.,  200.};

    // create sexa histograms
    TH1F *hist_sexa[n_sexa_vars];
    TH1F *hist_sexa_signal[n_sexa_vars];
    for (Int_t vv = 0; vv < n_sexa_vars; vv++) {
        hist_sexa[vv] = new TH1F(sexa_var_name[vv], sexa_var_name[vv],  //
                                 sexa_var_n_bins[vv], sexa_var_min_x[vv], sexa_var_max_x[vv]);
        hist_sexa_signal[vv] = new TH1F(sexa_var_name[vv] + "_signal", sexa_var_name[vv] + "_signal",  //
                                        sexa_var_n_bins[vv], sexa_var_min_x[vv], sexa_var_max_x[vv]);
    }

    // define iterator
    TListIter *list_of_trees_it = new TListIter(list_of_trees);

    // loop over collected trees
    while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

        // load MC particles branches
        this_tree->SetBranchAddress("MC_PID", &MC_PID);
        this_tree->SetBranchAddress("MC_Mother", &MC_Mother);
        this_tree->SetBranchAddress("MC_isSignal", &MC_isSignal);
        // load track branches
        this_tree->SetBranchAddress("N_MCRec", &N_MCRec);
        this_tree->SetBranchAddress("Idx_True", &Idx_True);
        this_tree->SetBranchAddress("Rec_Px", &Rec_Px);
        this_tree->SetBranchAddress("Rec_Py", &Rec_Py);
        this_tree->SetBranchAddress("Rec_Pz", &Rec_Pz);
        this_tree->SetBranchAddress("Rec_Charge", &Rec_Charge);
        this_tree->SetBranchAddress("Rec_NSigmaPion", &Rec_NSigmaPion);
        this_tree->SetBranchAddress("Rec_NSigmaProton", &Rec_NSigmaProton);
        this_tree->SetBranchAddress("Rec_isDuplicate", &Rec_isDuplicate);
        this_tree->SetBranchAddress("Rec_isSignal", &Rec_isSignal);
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
        this_tree->SetBranchAddress("V0_DCA_wrtPV", &V0_DCA_wrtPV);

        // loop over events
        for (Int_t event = 0; event < this_tree->GetEntries(); event++) {
            // for (Int_t event = 0; event < 5; event++) {

            this_tree->GetEntry(event);

            /* V0s */

            // combine V0s
            for (Int_t v0_A = 0; v0_A < N_V0s - 1; v0_A++) {
                for (Int_t v0_B = v0_A + 1; v0_B < N_V0s; v0_B++) {

                    // (cut) none of them must on-the-fly (PENDING: to delete later)
                    if ((*V0_onFlyStatus)[v0_A] || (*V0_onFlyStatus)[v0_B]) {
                        continue;
                    }

                    // (cut) first, check that all four granddaughters are different
                    if ((*Idx_Pos)[v0_A] == (*Idx_Pos)[v0_B] || (*Idx_Neg)[v0_A] == (*Idx_Neg)[v0_B]) {
                        continue;
                    }

                    // define vectors
                    TVector3 tv3_k0_mom;
                    TVector3 tv3_k0_pos;
                    TVector3 tv3_k0_pca;  // tmp
                    TVector3 tv3_al_mom;
                    TVector3 tv3_al_pos;
                    TVector3 tv3_al_pca;  // tmp
                    TVector3 tv3_sexaquark_mom;
                    TVector3 tv3_sexaquark_pos;       // secondary vertex
                    TVector3 tv3_origin(0., 0., 0.);  // PENDING: it should be the PV, not origin

                    // declare variables
                    Float_t sexa_dca_daus;
                    Float_t sexa_dca_k0sv;
                    Float_t sexa_dca_alsv;
                    Float_t sexa_dca_spv;
                    Float_t sexa_cpa_k0sv;
                    Float_t sexa_cpa_alsv;
                    Float_t sexa_cpa_spv;
                    Float_t sexa_momentum;
                    Float_t sexa_pt;
                    Float_t sexa_energy;
                    Float_t sexa_premass;  // tmp
                    Float_t sexa_mass;
                    Float_t sexa_radius;
                    Bool_t is_signal;

                    // (cut) we're not allowing ambiguity!
                    // both particles must be different V0s!
                    if ((*V0_couldBeK0)[v0_A] && (*V0_couldBeAL)[v0_B] && !(*V0_couldBeAL)[v0_A] && !(*V0_couldBeK0)[v0_B]) {
                        tv3_k0_mom.SetXYZ((*V0_Px)[v0_A], (*V0_Py)[v0_A], (*V0_Pz)[v0_A]);
                        tv3_k0_pos.SetXYZ((*V0_X)[v0_A], (*V0_Y)[v0_A], (*V0_Z)[v0_A]);
                        tv3_al_mom.SetXYZ((*V0_Px)[v0_B], (*V0_Py)[v0_B], (*V0_Pz)[v0_B]);
                        tv3_al_pos.SetXYZ((*V0_X)[v0_B], (*V0_Y)[v0_B], (*V0_Z)[v0_B]);
                        sexa_energy = (*V0_E_asK0)[v0_A] + (*V0_E_asAL)[v0_B] - 0.939565;  // assuming struck neutron at rest
                    } else if ((*V0_couldBeK0)[v0_B] && (*V0_couldBeAL)[v0_A] && !(*V0_couldBeK0)[v0_A] && !(*V0_couldBeAL)[v0_B]) {
                        tv3_al_mom.SetXYZ((*V0_Px)[v0_A], (*V0_Py)[v0_A], (*V0_Pz)[v0_A]);
                        tv3_al_pos.SetXYZ((*V0_X)[v0_A], (*V0_Y)[v0_A], (*V0_Z)[v0_A]);
                        tv3_k0_mom.SetXYZ((*V0_Px)[v0_B], (*V0_Py)[v0_B], (*V0_Pz)[v0_B]);
                        tv3_k0_pos.SetXYZ((*V0_X)[v0_B], (*V0_Y)[v0_B], (*V0_Z)[v0_B]);
                        sexa_energy = (*V0_E_asAL)[v0_A] + (*V0_E_asK0)[v0_B] - 0.939565;  // assuming struck neutron at rest
                    } else {
                        continue;
                    }

                    // get variables
                    sexa_dca_daus = TwoLinesDCA(tv3_k0_pos, tv3_k0_mom, tv3_al_pos, tv3_al_mom, tv3_k0_pca, tv3_al_pca);
                    Double_t half_distance = (tv3_al_pca - tv3_k0_pca).Mag() / 2.;
                    tv3_sexaquark_pos = tv3_k0_pca + half_distance * (tv3_al_pca - tv3_k0_pca).Unit();
                    sexa_radius = TMath::Sqrt(tv3_sexaquark_pos.X() * tv3_sexaquark_pos.X() +  //
                                              tv3_sexaquark_pos.Y() * tv3_sexaquark_pos.Y());

                    tv3_sexaquark_mom = tv3_k0_mom + tv3_al_mom;
                    sexa_momentum = tv3_sexaquark_mom.Mag();
                    sexa_premass = sexa_energy * sexa_energy - sexa_momentum * sexa_momentum;  // (protection)
                    sexa_mass = sexa_premass >= 0 ? TMath::Sqrt(sexa_premass) : -1;
                    sexa_pt = TMath::Sqrt(tv3_sexaquark_mom.X() * tv3_sexaquark_mom.X() +  //
                                          tv3_sexaquark_mom.Y() * tv3_sexaquark_mom.Y());

                    sexa_dca_k0sv = LinePointDCA(tv3_k0_mom, tv3_k0_pos, tv3_sexaquark_pos);
                    sexa_dca_alsv = LinePointDCA(tv3_al_mom, tv3_al_pos, tv3_sexaquark_pos);
                    sexa_dca_spv = LinePointDCA(tv3_sexaquark_mom, tv3_sexaquark_pos, tv3_origin);
                    sexa_cpa_k0sv = CPA(tv3_k0_mom, tv3_k0_pos, tv3_sexaquark_pos);
                    sexa_cpa_alsv = CPA(tv3_al_mom, tv3_al_pos, tv3_sexaquark_pos);
                    sexa_cpa_spv = CPA(tv3_sexaquark_mom, tv3_sexaquark_pos, tv3_origin);
                    is_signal = (*V0_isSignal)[v0_A] && (*V0_isSignal)[v0_B];

                    // "dca(K0,AL)" "dca(K0,S_v)" "dca(AL,S_v)" "dca(S,P_v)"
                    // "cpa(K0,S_v)" "cpa(AL,S_v)" "cpa(S,P_v)"
                    // "pt(S)" "m(S)" "radius(S_v)"

                    // (fill hists)
                    hist_sexa[0]->Fill(sexa_dca_daus);
                    hist_sexa[1]->Fill(sexa_dca_k0sv);
                    hist_sexa[2]->Fill(sexa_dca_alsv);
                    hist_sexa[3]->Fill(sexa_dca_spv);
                    hist_sexa[4]->Fill(sexa_cpa_k0sv);
                    hist_sexa[5]->Fill(sexa_cpa_alsv);
                    hist_sexa[6]->Fill(sexa_cpa_spv);
                    hist_sexa[7]->Fill(sexa_pt);
                    hist_sexa[8]->Fill(sexa_mass);
                    hist_sexa[9]->Fill(sexa_radius);
                    if (is_signal) {
                        hist_sexa_signal[0]->Fill(sexa_dca_daus);
                        hist_sexa_signal[1]->Fill(sexa_dca_k0sv);
                        hist_sexa_signal[2]->Fill(sexa_dca_alsv);
                        hist_sexa_signal[3]->Fill(sexa_dca_spv);
                        hist_sexa_signal[4]->Fill(sexa_cpa_k0sv);
                        hist_sexa_signal[5]->Fill(sexa_cpa_alsv);
                        hist_sexa_signal[6]->Fill(sexa_cpa_spv);
                        hist_sexa_signal[7]->Fill(sexa_pt);
                        hist_sexa_signal[8]->Fill(sexa_mass);
                        hist_sexa_signal[9]->Fill(sexa_radius);
                    }

                }  // end of loop over V0_B
            }      // end of loop over V0_A

        }  // end of loop over events

    }  // end of loop over trees

    /* Draw */

    // set style to everything
    SetMyStyle();
    gStyle->SetOptStat(0);

    TString output_filename;
    TCanvas *c = new TCanvas("c", "c", 1080, 1080);

    c->SetLogy(1);

    // loop over sexaquark variables
    for (Int_t vv = 0; vv < n_sexa_vars; vv++) {

        // "dca(K0,AL)" "dca(K0,S_v)" "dca(AL,S_v)" "dca(S,P_v)"
        // "cpa(K0,S_v)" "cpa(AL,S_v)" "cpa(S,P_v)"
        // "pt(S)" "m(S)" "radius(S_v)"

        SetMyHistStyle(hist_sexa[vv]);
        hist_sexa[vv]->SetFillStyle(0);
        hist_sexa[vv]->SetLineWidth(3);
        hist_sexa[vv]->SetLineColor(myBlack);
        hist_sexa[vv]->SetMinimum(1);
        hist_sexa[vv]->SetMaximum(1000000);  // adjusted for one file
        hist_sexa[vv]->GetYaxis()->SetTitle("Counts");
        hist_sexa[vv]->GetXaxis()->SetTitle(sexa_var_name[vv]);
        hist_sexa[vv]->Draw();

        hist_sexa_signal[vv]->SetLineWidth(0);
        hist_sexa_signal[vv]->SetFillStyle(kSolid);
        hist_sexa_signal[vv]->SetFillColorAlpha(myGreen, 0.5);
        hist_sexa_signal[vv]->SetLineWidth(0);
        hist_sexa_signal[vv]->Draw("SAME");
        hist_sexa_signal[vv]->Draw("AXIS SAME");

        output_filename = hist_sexa[vv]->GetName();
        c->Print(output_dir + output_filename + ".png");
    }  // end of loop over sexaquark variables
}
