#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/Utilities.hxx"

// this macro is also a mess...
// -- A. Bórquez

//_____________________________________________________________________________
void V0s_InvariantMass(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595_v1/AnalysisResults_CustomV0s_000.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

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
    std::vector<Int_t> *Idx_True = 0;            // index of true MC particle
    std::vector<Float_t> *Rec_Px = 0;            //
    std::vector<Float_t> *Rec_Py = 0;            //
    std::vector<Float_t> *Rec_Pz = 0;            //
    std::vector<Bool_t> *Rec_isDuplicate = 0;    // kTRUE if track is a duplicate, kFALSE if not
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
    // inv mass K0
    TH1F *inv_mass_K0 = new TH1F("inv_mass_K0", "inv_mass_K0", 80, 0.4, 0.6);
    TH1F *inv_mass_K0_K0Signal = new TH1F("inv_mass_K0_K0Signal", "inv_mass_K0_K0Signal", 80, 0.4, 0.6);
    TH1F *inv_mass_K0_SexaSignal = new TH1F("inv_mass_K0_SexaSignal", "inv_mass_K0_SexaSignal", 80, 0.4, 0.6);
    // inv mass anti-lambdas
    TH1F *inv_mass_AL = new TH1F("inv_mass_AL", "inv_mass_AL", 80, 1.095, 1.135);
    TH1F *inv_mass_AL_ALSignal = new TH1F("inv_mass_AL_ALSignal", "inv_mass_AL_ALSignal", 80, 1.095, 1.135);
    TH1F *inv_mass_AL_SexaSignal = new TH1F("inv_mass_AL_SexaSignal", "inv_mass_AL_SexaSignal", 80, 1.095, 1.135);

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

            // declare variables
            Bool_t same_mc_mother;
            Bool_t mc_is_K0;
            Bool_t mc_is_AL;
            Bool_t is_sexa_signal;
            Float_t mass_asK0;
            Float_t mass_asAL;

            // declare possible cuts
            Bool_t cut_radius;
            Bool_t cut_armenteros;
            Bool_t cut_dca_daughters;
            Bool_t cut_no_duplicate;
            Bool_t cut_cpa_wrt_pv;

            Bool_t K0_Selection;
            Bool_t AL_Selection;

            // loop over V0s
            for (Int_t v0 = 0; v0 < N_V0s; v0++) {

                mass_asK0 = TMath::Sqrt((*V0_E_asK0)[v0] * (*V0_E_asK0)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                        (*V0_Pz)[v0] * (*V0_Pz)[v0]);
                mass_asAL = TMath::Sqrt((*V0_E_asAL)[v0] * (*V0_E_asAL)[v0] - (*V0_Px)[v0] * (*V0_Px)[v0] - (*V0_Py)[v0] * (*V0_Py)[v0] -
                                        (*V0_Pz)[v0] * (*V0_Pz)[v0]);

                same_mc_mother = (*MC_Mother)[(*Idx_True)[(*Idx_Pos)[v0]]] == (*MC_Mother)[(*Idx_True)[(*Idx_Neg)[v0]]];

                mc_is_K0 = (*MC_PID)[(*MC_Mother)[(*Idx_True)[(*Idx_Pos)[v0]]]] == 310 &&
                           (*MC_PID)[(*MC_Mother)[(*Idx_True)[(*Idx_Neg)[v0]]]] == 310;
                mc_is_AL = (*MC_PID)[(*MC_Mother)[(*Idx_True)[(*Idx_Pos)[v0]]]] == -3122 &&
                           (*MC_PID)[(*MC_Mother)[(*Idx_True)[(*Idx_Neg)[v0]]]] == -3122;

                is_sexa_signal = (*V0_isSignal)[v0];

                // define cuts
                cut_no_duplicate = !(*Rec_isDuplicate)[(*Idx_Pos)[v0]] && !(*Rec_isDuplicate)[(*Idx_Neg)[v0]];

                cut_radius = TMath::Sqrt((*V0_X)[v0] * (*V0_X)[v0] + (*V0_Y)[v0] * (*V0_Y)[v0]) > 5 &&
                             TMath::Sqrt((*V0_X)[v0] * (*V0_X)[v0] + (*V0_Y)[v0] * (*V0_Y)[v0]) < 180;

                cut_armenteros = 0.2 * TMath::Abs((*V0_ArmAlpha)[v0]) > (*V0_ArmPt)[v0];

                cut_dca_daughters = (*V0_DCA_Daughters)[v0] < .25;  // cm

                cut_cpa_wrt_pv = (*V0_CPA_wrtPV)[v0];

                // define selections: CHOOSE cuts
                K0_Selection = (*V0_couldBeK0)[v0] && cut_no_duplicate;
                AL_Selection = (*V0_couldBeAL)[v0] && cut_no_duplicate;

                if (K0_Selection) {
                    inv_mass_K0->Fill(mass_asK0);
                    if (same_mc_mother && mc_is_K0) {
                        inv_mass_K0_K0Signal->Fill(mass_asK0);
                        if (is_sexa_signal) {
                            inv_mass_K0_SexaSignal->Fill(mass_asK0);
                        }
                    }
                }

                if (AL_Selection) {
                    inv_mass_AL->Fill(mass_asAL);
                    if (same_mc_mother && mc_is_AL) {
                        inv_mass_AL_ALSignal->Fill(mass_asAL);
                        if (is_sexa_signal) {
                            inv_mass_AL_SexaSignal->Fill(mass_asAL);
                        }
                    }
                }

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

    /* HIST 1 */
    // inv mass neutral kaon short

    SetMyHistStyle(inv_mass_K0);
    inv_mass_K0->SetFillStyle(0);
    inv_mass_K0->SetLineWidth(3);
    inv_mass_K0->SetLineColor(myBlack);
    inv_mass_K0->SetMinimum(1);
    inv_mass_K0->SetMaximum(1000000);  // adjusted for one file
    inv_mass_K0->GetYaxis()->SetTitle("Counts");
    inv_mass_K0->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    inv_mass_K0->Draw();

    inv_mass_K0_K0Signal->SetLineWidth(0);
    inv_mass_K0_K0Signal->SetFillStyle(kSolid);
    inv_mass_K0_K0Signal->SetFillColor(myViolet);
    inv_mass_K0_K0Signal->SetLineWidth(0);
    inv_mass_K0_K0Signal->Draw("SAME");

    inv_mass_K0_SexaSignal->SetLineWidth(0);
    inv_mass_K0_SexaSignal->SetFillStyle(kSolid);
    inv_mass_K0_SexaSignal->SetFillColor(myGreen);
    inv_mass_K0_SexaSignal->SetLineWidth(0);
    inv_mass_K0_SexaSignal->Draw("SAME");
    inv_mass_K0_SexaSignal->Draw("AXIS SAME");

    DrawText(Form("N(Rec. K^{0}_{S} Selection) = %i", (Int_t)inv_mass_K0->GetEntries()), myBlack, 0.55, 0.84);
    DrawText(Form("N(True K^{0}_{S}) = %i", (Int_t)inv_mass_K0_K0Signal->GetEntries()), myViolet, 0.55, 0.79);
    DrawText(Form("N(#bar{S} Signal) = %i", (Int_t)inv_mass_K0_SexaSignal->GetEntries()), myGreen, 0.55, 0.74);

    output_filename = inv_mass_K0->GetName();
    c->Print(output_dir + output_filename + ".png");

    /* HIST 2 */
    // inv mass anti-lambdas

    SetMyHistStyle(inv_mass_AL);
    inv_mass_AL->SetFillStyle(0);
    inv_mass_AL->SetLineWidth(3);
    inv_mass_AL->SetLineColor(myBlack);
    inv_mass_AL->SetMinimum(1);
    inv_mass_AL->SetMaximum(10000);  // adjusted for one file
    inv_mass_AL->GetYaxis()->SetTitle("Counts");
    inv_mass_AL->GetXaxis()->SetTitle("m(#pi^{+}#bar{p}) [GeV/c^{2}]");
    inv_mass_AL->Draw();

    inv_mass_AL_ALSignal->SetLineWidth(0);
    inv_mass_AL_ALSignal->SetFillStyle(kSolid);
    inv_mass_AL_ALSignal->SetFillColor(myViolet);
    inv_mass_AL_ALSignal->SetLineWidth(0);
    inv_mass_AL_ALSignal->Draw("SAME");

    inv_mass_AL_SexaSignal->SetLineWidth(0);
    inv_mass_AL_SexaSignal->SetFillStyle(kSolid);
    inv_mass_AL_SexaSignal->SetFillColor(myGreen);
    inv_mass_AL_SexaSignal->SetLineWidth(0);
    inv_mass_AL_SexaSignal->Draw("SAME");
    inv_mass_AL_SexaSignal->Draw("AXIS SAME");

    DrawText(Form("N(Rec. #bar{#Lambda} Selection) = %i", (Int_t)inv_mass_AL->GetEntries()), myBlack, 0.55, 0.84);
    DrawText(Form("N(True #bar{#Lambda}) = %i", (Int_t)inv_mass_AL_ALSignal->GetEntries()), myViolet, 0.55, 0.79);
    DrawText(Form("N(#bar{S} Signal) = %i", (Int_t)inv_mass_AL_SexaSignal->GetEntries()), myGreen, 0.55, 0.74);

    output_filename = inv_mass_AL->GetName();
    c->Print(output_dir + output_filename + ".png");
}
