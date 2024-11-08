#include "include/Headers.hxx"
#include "include/Style.cxx"

using namespace ROOT;

/*
 *
 */
void Trees_RDF_ResponseMatrix(TString InputFileName = "../output/local_signalMC_15o_full_kalmanA1.8/AnalysisResults_merged.root",  //
                              TString OutputFilename = "output_rdf/ResponseMatrix_local_signalMC_15o_full_kalmanA1.8.root") {

    std::cout << "!! Trees_RDF_ResponseMatrix !! Starting !!" << std::endl;

    /*** Input ***/

    /* Input: Get events tree */

    RDataFrame RDF_Events("Events", InputFileName);

    std::vector<Int_t> Events_RunNumber = RDF_Events.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Events_DirNumber = RDF_Events.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Events_EventNumber = RDF_Events.Take<Int_t>("EventNumber").GetValue();
    std::vector<Float_t> Events_Centrality = RDF_Events.Take<Float_t>("Centrality").GetValue();
    const Int_t N_Events = (Int_t)Events_RunNumber.size();

    std::map<std::tuple<Int_t, Int_t, Int_t>, Float_t> Map_Centrality;
    for (Int_t i = 0; i < N_Events; i++) {
        Map_Centrality[std::make_tuple(Events_RunNumber[i], Events_DirNumber[i], Events_EventNumber[i])] = Events_Centrality[i];
    }

    auto Lambda_GetCentrality = [&Map_Centrality](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber) {
        return Map_Centrality[std::make_tuple(RunNumber, DirNumber, EventNumber)];
    };

    /* Input: Get sexaquark trees */

    RDataFrame RDF_Injected("Injected", InputFileName);
    RDataFrame RDF_Sexaquarks("Sexaquarks", InputFileName);

    auto fRDF_Injected = RDF_Injected  //
                             .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})
                             .Filter("Centrality < 90.");
    auto fRDF_Sexaquarks = RDF_Sexaquarks  //
                               .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})
                               .Filter("Centrality < 90.")
                               .Filter("IsSignal");

    /*** Preparation ***/

    /* Preparation: same-size vectors */

    std::vector<Int_t> Injected_RunNumber = fRDF_Injected.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Injected_DirNumber = fRDF_Injected.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Injected_EventNumber = fRDF_Injected.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> Injected_ReactionID = fRDF_Injected.Take<Int_t>("ReactionID").GetValue();
    std::vector<Float_t> Injected_Px = fRDF_Injected.Take<Float_t>("Px").GetValue();
    std::vector<Float_t> Injected_Py = fRDF_Injected.Take<Float_t>("Py").GetValue();
    const Int_t N_Injected = (Int_t)Injected_ReactionID.size();

    std::vector<Int_t> Sexaquark_RunNumber = fRDF_Sexaquarks.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Sexaquark_DirNumber = fRDF_Sexaquarks.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Sexaquark_EventNumber = fRDF_Sexaquarks.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> Sexaquark_ReactionID = fRDF_Sexaquarks.Take<Int_t>("ReactionID").GetValue();
    std::vector<Float_t> Sexaquark_Px = fRDF_Sexaquarks.Take<Float_t>("Px").GetValue();
    std::vector<Float_t> Sexaquark_Py = fRDF_Sexaquarks.Take<Float_t>("Py").GetValue();
    const Int_t N_Sexaquarks = (Int_t)Sexaquark_ReactionID.size();

    /* Preparation: True Pt maps */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_PtTrue;
    for (Int_t i = 0; i < N_Injected; i++) {
        Map_PtTrue[std::make_tuple(Injected_RunNumber[i], Injected_DirNumber[i], Injected_EventNumber[i], Injected_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)Injected_Px[i] * (Double_t)Injected_Px[i] + (Double_t)Injected_Py[i] * (Double_t)Injected_Py[i]);
    }

    auto Lambda_GetPtTrue = [&Map_PtTrue](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) -> Float_t {
        return Map_PtTrue[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Preparation: Extend trees */

    auto efRDF_Injected = fRDF_Injected  //
                              .Define("Pt", "static_cast<Float_t>(TMath::Sqrt((Double_t)Px * (Double_t)Px + (Double_t)Py * (Double_t)Py))");
    auto efRDF_Sexaquarks = fRDF_Sexaquarks
                                .Define("PtReco", "static_cast<Float_t>(TMath::Sqrt((Double_t)Px * (Double_t)Px + (Double_t)Py * (Double_t)Py))")  //
                                .Define("PtTrue", Lambda_GetPtTrue, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"});

    /*** Response Matrix ***/

    Int_t n_bins = 24;
    Float_t pt_min = 0.;
    Float_t pt_max = 6.;
    Float_t pt_step = (pt_max - pt_min) / (Float_t)n_bins;

    TH2F* Hist_ResponseMatrix =
        new TH2F("ResponseMatrix", Form("ResponseMatrix;p_{T}(true) / %.2f GeV/c;p_{T}(reco) / %.2f GeV/c", pt_step, pt_step),  //
                 n_bins, pt_min, pt_max, n_bins, pt_min, pt_max);

    Double_t N_ij;
    Double_t N_j;
    Float_t A_ij;

    Float_t truth_bin_min, truth_bin_max;
    Float_t reco_bin_min, reco_bin_max;

    /*  */
    for (Int_t truth_bin_j = 0; truth_bin_j < n_bins; truth_bin_j++) {
        truth_bin_min = pt_min + truth_bin_j * pt_step;
        truth_bin_max = truth_bin_min + pt_step;

        auto auxRDF_MCGen = efRDF_Injected  //
                                .Filter(Form("%f < Pt && Pt < %f", truth_bin_min, truth_bin_max));

        N_j = auxRDF_MCGen.Count().GetValue();
        // std::cout << "!! Trees_RDF_ResponseMatrix !! N_" << truth_bin_j << " (" << truth_bin_min << " < PtTrue < " << truth_bin_max << ") = " <<
        // N_j << std::endl;

        /*  */
        for (Int_t reco_bin_i = 0; reco_bin_i < n_bins; reco_bin_i++) {
            reco_bin_min = pt_min + reco_bin_i * pt_step;
            reco_bin_max = reco_bin_min + pt_step;

            auto auxRDF_MCRec = efRDF_Sexaquarks                                                              //
                                    .Filter(Form("%f < PtTrue && PtTrue < %f && %f < PtReco && PtReco < %f",  //
                                                 truth_bin_min, truth_bin_max, reco_bin_min, reco_bin_max));

            N_ij = auxRDF_MCRec.Count().GetValue();
            // std::cout << "!! Trees_RDF_ResponseMatrix !! >> N_" << reco_bin_i << "_" << truth_bin_j << " (" << reco_bin_min << " < PtReco < " <<
            // reco_bin_max << ", " << truth_bin_min << " < TruePt < " << truth_bin_max << ") = " << N_ij << std::endl;

            A_ij = N_j ? static_cast<Float_t>(N_ij / N_j) : 0.;
            Hist_ResponseMatrix->SetBinContent(truth_bin_j + 1, reco_bin_i + 1, A_ij);
            // std::cout << "!! Trees_RDF_ResponseMatrix !! >> >> A_ij = " << A_ij << std::endl;
        }
    }

    /*** Draw ***/

    /* Draw: Response Matrix */

    SetMy2DStyle();
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    TCanvas* c = new TCanvas("c", "c", 1080, 2 * 1080);
    c->Update();

    c->Divide(1, 2);

    c->cd(1);
    gPad->SetPad(0.0, 0.5, 1.0, 1.0);
    gPad->SetTopMargin(0.10);
    gPad->SetBottomMargin(0.10);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);

    SetMy2DHistStyle(Hist_ResponseMatrix);
    Hist_ResponseMatrix->SetNdivisions(-406, "xy");
    Hist_ResponseMatrix->Draw("COLZ");
    DrawHorizontalLine(5., myRed.GetNumber());
    DrawVerticalLine(5., myRed.GetNumber());

    /* Draw: Efficiency vs p_T(true) */

    c->cd(2);
    gPad->SetPad(0.0, 0.0, 1.0, 0.5);
    gPad->SetTopMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);

    TH1D* Hist_Efficiency = Hist_ResponseMatrix->ProjectionX("Efficiency", 1, n_bins);
    Hist_Efficiency->Scale(100.);
    Hist_Efficiency->SetMinimum(0.);
    Hist_Efficiency->SetLineWidth(2);
    Hist_Efficiency->SetLineColor(myBlue.GetNumber());
    Hist_Efficiency->SetNdivisions(-406, "x");
    Hist_Efficiency->GetYaxis()->SetTitle("Efficiency / %");

    SetMyHistStyle(Hist_Efficiency);

    Hist_Efficiency->Draw("E");
    DrawVerticalLine(5., myRed.GetNumber());

    TString PictureFilename = OutputFilename;
    PictureFilename.ReplaceAll(".root", ".png");
    c->Print(PictureFilename);

    /*** Output ***/

    TFile* OutputFile = TFile::Open(OutputFilename, "RECREATE");

    Hist_ResponseMatrix->Write();
    Hist_Efficiency->Write();

    OutputFile->Close();

    std::cout << "!! Trees_RDF_ResponseMatrix !! TFile " << OutputFilename << " has been created !!" << std::endl;

    /* Print how many event-loops were executed */

    std::cout << "!! Trees_RDF_ResponseMatrix !! NRuns !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! ===== !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! RDF_Events       = " << RDF_Events.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! RDF_Injected     = " << RDF_Injected.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! fRDF_Injected    = " << fRDF_Injected.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! efRDF_Injected   = " << efRDF_Injected.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! RDF_Sexaquarks   = " << RDF_Sexaquarks.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! fRDF_Sexaquarks  = " << fRDF_Sexaquarks.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! efRDF_Sexaquarks = " << efRDF_Sexaquarks.GetNRuns() << " !!" << std::endl;

    std::cout << "!! Trees_RDF_ResponseMatrix !! Finished !!" << std::endl;
}
