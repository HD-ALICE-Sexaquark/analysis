#include "include/Headers.hxx"
#include "include/Style.hxx"
#include "include/TreeFunctions.hxx"
#include "include/Utilities.hxx"

// V0 Analysis, Background Studies, Particle Identification
// -- A. Bórquez

//_____________________________________________________________________________
void V0s_PlotPIDContamination(
    TString input_filename = "/misc/alidata121/alice_u/borquez/analysis/output/signal+bkg/297595/AnalysisResults_CustomV0s_000.root",
    TString output_dir = "gfx/signal+bkg/custom_v0s/") {

    /*** Process Input ***/

    TList *list_of_trees = new TList();
    AddTreesToList(list_of_trees, input_filename + "/Trees/Events");

    Event_tt this_event;

    TString particle_opt = "K0";
    Int_t particle_opt_pdg_code = 310;
    TString particle_opt_latex = "K^{0}_{S}";
    Float_t particle_opt_min_mass = 0.4;
    Float_t particle_opt_max_mass = 0.6;

    const Int_t n_bins = 80;

    /*** Process Output ***/

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

    Int_t pos_particles_pdg[5] = {211, 2212, 321, -11, 0};
    TString pos_particles_name[5] = {"Pip", "P", "Kp", "Posi", "Else"};
    Int_t neg_particles_pdg[5] = {-211, -2212, -321, 11, 0};
    TString neg_particles_name[5] = {"Pim", "AntiP", "Km", "Elec", "Else"};

    std::map<Int_t, TString> fs_particles_latex;
    fs_particles_latex[211] = "#pi^{+}";
    fs_particles_latex[-211] = "#pi^{-}";
    fs_particles_latex[2212] = "p";
    fs_particles_latex[-2212] = "#bar{p}";
    fs_particles_latex[321] = "K^{+}";
    fs_particles_latex[-321] = "K^{-}";
    fs_particles_latex[11] = "e^{-}";
    fs_particles_latex[-11] = "e^{+}";
    fs_particles_latex[0] = "else";

    // create histograms
    TH1F *nominal_hist = new TH1F(particle_opt + "_BS_PID",  //
                                  particle_opt + "_BS_PID",  //
                                  n_bins,                    //
                                  particle_opt_min_mass,     //
                                  particle_opt_max_mass);
    std::map<std::pair<Int_t, Int_t>, TH1F *> test_hist;

    for (Int_t pp = 0; pp < 5; pp++) {
        for (Int_t nn = 0; nn < 5; nn++) {
            test_hist[{pos_particles_pdg[pp], neg_particles_pdg[nn]}] =
                new TH1F(particle_opt + "_" + pos_particles_name[pp] + "_" + neg_particles_name[nn],  //
                         particle_opt + "_" + pos_particles_name[pp] + "_" + neg_particles_name[nn],  //
                         n_bins,                                                                      //
                         particle_opt_min_mass,                                                       //
                         particle_opt_max_mass);
        }
    }

    // loop over collected trees
    TListIter *list_of_trees_it = new TListIter(list_of_trees);
    while (TTree *this_tree = (TTree *)list_of_trees_it->Next()) {

        // load branches
        LoadBranches(this_tree, this_event);

        // loop over events
        for (Int_t event = 0; event < this_tree->GetEntries(); event++) {

            this_tree->GetEntry(event);

            // define variables
            Float_t v0_mass;
            Float_t v0_energy;
            Bool_t v0_selection;

            Bool_t pos_particle_known;
            Bool_t neg_particle_known;

            // loop over V0s
            for (Int_t v0 = 0; v0 < this_event.N_V0s; v0++) {

                // get mass
                v0_energy = (particle_opt == "K0") * (*this_event.V0_E_asK0)[v0] + (particle_opt == "AL") * (*this_event.V0_E_asAL)[v0];
                v0_mass = TMath::Sqrt(v0_energy * v0_energy -                              //
                                      (*this_event.V0_Px)[v0] * (*this_event.V0_Px)[v0] -  //
                                      (*this_event.V0_Py)[v0] * (*this_event.V0_Py)[v0] -  //
                                      (*this_event.V0_Pz)[v0] * (*this_event.V0_Pz)[v0]);

                // define selections
                v0_selection = (*this_event.V0_couldBeK0)[v0] && !(*this_event.V0_onFlyStatus)[v0] &&  //
                               v0_mass > particle_opt_min_mass && v0_mass < particle_opt_max_mass;

                pos_particle_known = (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] == 211 ||   //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] == 2212 ||  //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] == 321 ||   //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Pos)[v0]]] == -11;

                neg_particle_known = (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]] == -211 ||   //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]] == -2212 ||  //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]] == -321 ||   //
                                     (*this_event.MC_PID)[(*this_event.Idx_True)[(*this_event.Idx_Neg)[v0]]] == 11;

                // fill histograms
                // apply selections
                if (v0_selection) {
                    nominal_hist->Fill(v0_mass);
                    // fill each corr. hist
                    if (pos_particle_known && neg_particle_known) {
                        test_hist[{(*MC_PID)[(*Idx_True)[(*Idx_Pos)[v0]]], (*MC_PID)[(*Idx_True)[(*Idx_Neg)[v0]]]}]->Fill(v0_mass);
                    } else {
                        test_hist[{0, 0}]->Fill(v0_mass);
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

    /*** PI PLUS ***/

    SetMyHistStyle(nominal_hist);
    nominal_hist->SetFillStyle(0);
    nominal_hist->SetLineWidth(3);
    nominal_hist->SetLineColor(myBlack);
    nominal_hist->SetMinimum(1);
    nominal_hist->SetMaximum(10000000);  // adjusted for one file
    nominal_hist->GetYaxis()->SetTitle("Counts");
    nominal_hist->GetXaxis()->SetTitle("m(#pi^{+}#pi^{-}) [GeV/c^{2}]");
    nominal_hist->Draw();

    // test 1: pi+ pi-
    test_hist[{211, -211}]->SetFillStyle(0);
    test_hist[{211, -211}]->SetLineWidth(3);
    test_hist[{211, -211}]->SetLineColor(myBlue);
    test_hist[{211, -211}]->Draw("SAME");

    // test 3: pi+ k-
    test_hist[{211, -321}]->SetFillStyle(0);
    test_hist[{211, -321}]->SetLineWidth(3);
    test_hist[{211, -321}]->SetLineColor(myGreen);
    test_hist[{211, -321}]->Draw("SAME");

    // test 4: pi+ elec
    test_hist[{211, 11}]->SetFillStyle(0);
    test_hist[{211, 11}]->SetLineWidth(3);
    test_hist[{211, 11}]->SetLineColor(myOrange);
    test_hist[{211, 11}]->Draw("SAME");

    // test 2: pi+ anti-p
    test_hist[{211, -2212}]->SetFillStyle(0);
    test_hist[{211, -2212}]->SetLineWidth(3);
    test_hist[{211, -2212}]->SetLineColor(myRed);
    test_hist[{211, -2212}]->Draw("SAME");

    // test 5: else
    test_hist[{0, 0}]->SetFillStyle(0);
    test_hist[{0, 0}]->SetLineWidth(3);
    test_hist[{0, 0}]->SetLineColor(myViolet);
    test_hist[{0, 0}]->Draw("SAME");

    TLegend *leg = new TLegend(0.2, 0.73, 0.5, 0.92, "", "NDC");      // x1,y1,x2,y2
    leg->AddEntry(nominal_hist, Form("rec %s%s",                      //
                                     fs_particles_latex[211].Data(),  //
                                     fs_particles_latex[-211].Data()));
    leg->AddEntry(test_hist[{211, -211}], Form("rec %s%s, true %s%s",           //
                                               fs_particles_latex[211].Data(),  //
                                               fs_particles_latex[-211].Data(),
                                               fs_particles_latex[211].Data(),     //
                                               fs_particles_latex[-211].Data()));  //
    leg->AddEntry(test_hist[{211, -321}], Form("rec %s%s, true %s%s",              //
                                               fs_particles_latex[211].Data(),     //
                                               fs_particles_latex[-211].Data(),
                                               fs_particles_latex[211].Data(),     //
                                               fs_particles_latex[-321].Data()));  //
    leg->AddEntry(test_hist[{211, 11}], Form("rec %s%s, true %s%s",                //
                                             fs_particles_latex[211].Data(),       //
                                             fs_particles_latex[-211].Data(),
                                             fs_particles_latex[211].Data(),     //
                                             fs_particles_latex[11].Data()));    //
    leg->AddEntry(test_hist[{211, -2212}], Form("rec %s%s, true %s%s",           //
                                                fs_particles_latex[211].Data(),  //
                                                fs_particles_latex[-211].Data(),
                                                fs_particles_latex[211].Data(),      //
                                                fs_particles_latex[-2212].Data()));  //
    leg->AddEntry(test_hist[{0, 0}], Form("rec %s%s, else*",                         //
                                          fs_particles_latex[211].Data(),            //
                                          fs_particles_latex[-211].Data()));
    leg->SetTextSize(0.025);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    output_filename = nominal_hist->GetName();
    c->Print(output_dir + output_filename + "_piplus.png");

    /*** PROTON ***/

    SetMyHistStyle(nominal_hist);
    nominal_hist->SetFillStyle(0);
    nominal_hist->SetLineWidth(3);
    nominal_hist->SetLineColor(myBlack);
    nominal_hist->SetMinimum(1);
    nominal_hist->SetMaximum(10000000);  // adjusted for one file
    nominal_hist->GetYaxis()->SetTitle("Counts");
    // nominal_hist->GetXaxis()->SetTitle("Decay Length(#pi^{+}#pi^{-}) [cm]"); // PENDING
    nominal_hist->Draw();

    // test 1: p pi-
    test_hist[{2212, -211}]->SetFillStyle(0);
    test_hist[{2212, -211}]->SetLineWidth(3);
    test_hist[{2212, -211}]->SetLineColor(myBlue);
    test_hist[{2212, -211}]->Draw("SAME");

    // test 3: p k-
    test_hist[{2212, -321}]->SetFillStyle(0);
    test_hist[{2212, -321}]->SetLineWidth(3);
    test_hist[{2212, -321}]->SetLineColor(myGreen);
    test_hist[{2212, -321}]->Draw("SAME");

    // test 4: p elec
    test_hist[{2212, 11}]->SetFillStyle(0);
    test_hist[{2212, 11}]->SetLineWidth(3);
    test_hist[{2212, 11}]->SetLineColor(myOrange);
    test_hist[{2212, 11}]->Draw("SAME");

    // test 2: p anti-p
    test_hist[{2212, -2212}]->SetFillStyle(0);
    test_hist[{2212, -2212}]->SetLineWidth(3);
    test_hist[{2212, -2212}]->SetLineColor(myRed);
    test_hist[{2212, -2212}]->Draw("SAME");

    // test 5: else
    test_hist[{0, 0}]->SetFillStyle(0);
    test_hist[{0, 0}]->SetLineWidth(3);
    test_hist[{0, 0}]->SetLineColor(myViolet);
    test_hist[{0, 0}]->Draw("SAME");

    leg->Clear();
    leg->AddEntry(nominal_hist, Form("rec %s%s",                      //
                                     fs_particles_latex[211].Data(),  //
                                     fs_particles_latex[-211].Data()));
    leg->AddEntry(test_hist[{2212, -211}], Form("rec %s%s, true %s%s",           //
                                                fs_particles_latex[211].Data(),  //
                                                fs_particles_latex[-211].Data(),
                                                fs_particles_latex[2212].Data(),    //
                                                fs_particles_latex[-211].Data()));  //
    leg->AddEntry(test_hist[{2212, -321}], Form("rec %s%s, true %s%s",              //
                                                fs_particles_latex[211].Data(),     //
                                                fs_particles_latex[-211].Data(),
                                                fs_particles_latex[2212].Data(),    //
                                                fs_particles_latex[-321].Data()));  //
    leg->AddEntry(test_hist[{2212, 11}], Form("rec %s%s, true %s%s",                //
                                              fs_particles_latex[211].Data(),       //
                                              fs_particles_latex[-211].Data(),
                                              fs_particles_latex[2212].Data(),    //
                                              fs_particles_latex[11].Data()));    //
    leg->AddEntry(test_hist[{2212, -2212}], Form("rec %s%s, true %s%s",           //
                                                 fs_particles_latex[211].Data(),  //
                                                 fs_particles_latex[-211].Data(),
                                                 fs_particles_latex[2212].Data(),     //
                                                 fs_particles_latex[-2212].Data()));  //
    leg->AddEntry(test_hist[{0, 0}], Form("rec %s%s, else*",                          //
                                          fs_particles_latex[211].Data(),             //
                                          fs_particles_latex[-211].Data()));
    leg->SetTextSize(0.025);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    output_filename = nominal_hist->GetName();
    c->Print(output_dir + output_filename + "_proton.png");

    /*** K PLUS ***/

    SetMyHistStyle(nominal_hist);
    nominal_hist->SetFillStyle(0);
    nominal_hist->SetLineWidth(3);
    nominal_hist->SetLineColor(myBlack);
    nominal_hist->SetMinimum(1);
    nominal_hist->SetMaximum(10000000);  // adjusted for one file
    nominal_hist->GetYaxis()->SetTitle("Counts");
    // nominal_hist->GetXaxis()->SetTitle("Decay Length(#pi^{+}#pi^{-}) [cm]"); // PENDING
    nominal_hist->Draw();

    // test 1: p pi-
    test_hist[{321, -211}]->SetFillStyle(0);
    test_hist[{321, -211}]->SetLineWidth(3);
    test_hist[{321, -211}]->SetLineColor(myBlue);
    test_hist[{321, -211}]->Draw("SAME");

    // test 3: p k-
    test_hist[{321, -321}]->SetFillStyle(0);
    test_hist[{321, -321}]->SetLineWidth(3);
    test_hist[{321, -321}]->SetLineColor(myGreen);
    test_hist[{321, -321}]->Draw("SAME");

    // test 4: p elec
    test_hist[{321, 11}]->SetFillStyle(0);
    test_hist[{321, 11}]->SetLineWidth(3);
    test_hist[{321, 11}]->SetLineColor(myOrange);
    test_hist[{321, 11}]->Draw("SAME");

    // test 2: p anti-p
    test_hist[{321, -2212}]->SetFillStyle(0);
    test_hist[{321, -2212}]->SetLineWidth(3);
    test_hist[{321, -2212}]->SetLineColor(myRed);
    test_hist[{321, -2212}]->Draw("SAME");

    // test 5: else
    test_hist[{0, 0}]->SetFillStyle(0);
    test_hist[{0, 0}]->SetLineWidth(3);
    test_hist[{0, 0}]->SetLineColor(myViolet);
    test_hist[{0, 0}]->Draw("SAME");

    leg->Clear();
    leg->AddEntry(nominal_hist, Form("rec %s%s",                      //
                                     fs_particles_latex[211].Data(),  //
                                     fs_particles_latex[-211].Data()));
    leg->AddEntry(test_hist[{321, -211}], Form("rec %s%s, true %s%s",                //
                                               fs_particles_latex[211].Data(),       //
                                               fs_particles_latex[-211].Data(),      //
                                               fs_particles_latex[321].Data(),       //
                                               fs_particles_latex[-211].Data()));    //
    leg->AddEntry(test_hist[{321, -321}], Form("rec %s%s, true %s%s",                //
                                               fs_particles_latex[211].Data(),       //
                                               fs_particles_latex[-211].Data(),      //
                                               fs_particles_latex[321].Data(),       //
                                               fs_particles_latex[-321].Data()));    //
    leg->AddEntry(test_hist[{321, 11}], Form("rec %s%s, true %s%s",                  //
                                             fs_particles_latex[211].Data(),         //
                                             fs_particles_latex[-211].Data(),        //
                                             fs_particles_latex[321].Data(),         //
                                             fs_particles_latex[11].Data()));        //
    leg->AddEntry(test_hist[{321, -2212}], Form("rec %s%s, true %s%s",               //
                                                fs_particles_latex[211].Data(),      //
                                                fs_particles_latex[-211].Data(),     //
                                                fs_particles_latex[321].Data(),      //
                                                fs_particles_latex[-2212].Data()));  //
    leg->AddEntry(test_hist[{0, 0}], Form("rec %s%s, else*",                         //
                                          fs_particles_latex[211].Data(),            //
                                          fs_particles_latex[-211].Data()));
    leg->SetTextSize(0.025);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    output_filename = nominal_hist->GetName();
    c->Print(output_dir + output_filename + "_kplus.png");

    /*** POSITRON ***/

    SetMyHistStyle(nominal_hist);
    nominal_hist->SetFillStyle(0);
    nominal_hist->SetLineWidth(3);
    nominal_hist->SetLineColor(myBlack);
    nominal_hist->SetMinimum(1);
    nominal_hist->SetMaximum(10000000);  // adjusted for one file
    nominal_hist->GetYaxis()->SetTitle("Counts");
    // nominal_hist->GetXaxis()->SetTitle("Decay Length(#pi^{+}#pi^{-}) [cm]"); // PENDING
    nominal_hist->Draw();

    // test 1: p pi-
    test_hist[{-11, -211}]->SetFillStyle(0);
    test_hist[{-11, -211}]->SetLineWidth(3);
    test_hist[{-11, -211}]->SetLineColor(myBlue);
    test_hist[{-11, -211}]->Draw("SAME");

    // test 3: p k-
    test_hist[{-11, -321}]->SetFillStyle(0);
    test_hist[{-11, -321}]->SetLineWidth(3);
    test_hist[{-11, -321}]->SetLineColor(myGreen);
    test_hist[{-11, -321}]->Draw("SAME");

    // test 4: p elec
    test_hist[{-11, 11}]->SetFillStyle(0);
    test_hist[{-11, 11}]->SetLineWidth(3);
    test_hist[{-11, 11}]->SetLineColor(myOrange);
    test_hist[{-11, 11}]->Draw("SAME");

    // test 2: p anti-p
    test_hist[{-11, -2212}]->SetFillStyle(0);
    test_hist[{-11, -2212}]->SetLineWidth(3);
    test_hist[{-11, -2212}]->SetLineColor(myRed);
    test_hist[{-11, -2212}]->Draw("SAME");

    // test 5: else
    test_hist[{0, 0}]->SetFillStyle(0);
    test_hist[{0, 0}]->SetLineWidth(3);
    test_hist[{0, 0}]->SetLineColor(myViolet);
    test_hist[{0, 0}]->Draw("SAME");

    leg->Clear();
    leg->AddEntry(nominal_hist, Form("rec %s%s",                      //
                                     fs_particles_latex[211].Data(),  //
                                     fs_particles_latex[-211].Data()));
    leg->AddEntry(test_hist[{-11, -211}], Form("rec %s%s, true %s%s",                //
                                               fs_particles_latex[211].Data(),       //
                                               fs_particles_latex[-211].Data(),      //
                                               fs_particles_latex[-11].Data(),       //
                                               fs_particles_latex[-211].Data()));    //
    leg->AddEntry(test_hist[{-11, -321}], Form("rec %s%s, true %s%s",                //
                                               fs_particles_latex[211].Data(),       //
                                               fs_particles_latex[-211].Data(),      //
                                               fs_particles_latex[-11].Data(),       //
                                               fs_particles_latex[-321].Data()));    //
    leg->AddEntry(test_hist[{-11, 11}], Form("rec %s%s, true %s%s",                  //
                                             fs_particles_latex[211].Data(),         //
                                             fs_particles_latex[-211].Data(),        //
                                             fs_particles_latex[-11].Data(),         //
                                             fs_particles_latex[11].Data()));        //
    leg->AddEntry(test_hist[{-11, -2212}], Form("rec %s%s, true %s%s",               //
                                                fs_particles_latex[211].Data(),      //
                                                fs_particles_latex[-211].Data(),     //
                                                fs_particles_latex[-11].Data(),      //
                                                fs_particles_latex[-2212].Data()));  //
    leg->AddEntry(test_hist[{0, 0}], Form("rec %s%s, else*",                         //
                                          fs_particles_latex[211].Data(),            //
                                          fs_particles_latex[-211].Data()));
    leg->SetTextSize(0.025);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    output_filename = nominal_hist->GetName();
    c->Print(output_dir + output_filename + "_positron.png");
}
