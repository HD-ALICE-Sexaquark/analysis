#include "include/Headers.hxx"
#include "include/Style.hxx"

using namespace ROOT;

/*
 *
 */
void Trees_RDF_ResponseMatrix(TString InputFileName = "../output/local_signalMC_15o_full_kalmanA1.8/AnalysisResults_merged.root",  //
                              TString OutputFileName = "output_rdf/EfficiencyResults_15o_test_kalmanA1.94_INNER.root") {

    std::cout << "!! Trees_RDF_ResponseMatrix !! Starting !!" << std::endl;

    /*** Input ***/

    /* Input: Get trees */

    RDataFrame RDF_Injected("Injected", InputFileName);
    RDataFrame RDF_Sexaquarks("Sexaquarks", InputFileName);

    /* Input: Filter tree */

    auto fRDF_Sexaquarks = RDF_Sexaquarks.Filter("IsSignal");

    /*** Preparation ***/

    /* Preparation: same-size vectors */

    std::vector<Int_t> Injected_RunNumber = RDF_Injected.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Injected_DirNumber = RDF_Injected.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Injected_EventNumber = RDF_Injected.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> Injected_ReactionID = RDF_Injected.Take<Int_t>("ReactionID").GetValue();
    std::vector<Float_t> Injected_Px = RDF_Injected.Take<Float_t>("Px").GetValue();
    std::vector<Float_t> Injected_Py = RDF_Injected.Take<Float_t>("Py").GetValue();
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

    auto eRDF_Injected = RDF_Injected  //
                             .Define("Pt", "static_cast<Float_t>(TMath::Sqrt((Double_t)Px * (Double_t)Px + (Double_t)Py * (Double_t)Py))");
    auto feRDF_Sexaquarks = fRDF_Sexaquarks
                                .Define("PtReco", "static_cast<Float_t>(TMath::Sqrt((Double_t)Px * (Double_t)Px + (Double_t)Py * (Double_t)Py))")  //
                                .Define("PtTrue", Lambda_GetPtTrue, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"});

    /*** Response Matrix ***/

    Int_t n_bins = 4;
    Float_t pt_min = 0.;
    Float_t pt_max = 8.;
    Float_t pt_step = (pt_max - pt_min) / (Float_t)n_bins;

    Float_t truth_bin_min, truth_bin_max;
    Float_t reco_bin_min, reco_bin_max;

    for (Int_t truth_bin_j = 0; truth_bin_j < n_bins; truth_bin_j++) {
        truth_bin_min = pt_min + truth_bin_j * pt_step;
        truth_bin_max = truth_bin_min + pt_step;

        auto auxRDF_MCGen = eRDF_Injected  //
                                .Filter(Form("%f < Pt && Pt < %f", truth_bin_min, truth_bin_max));

        std::cout << "!! Trees_RDF_ResponseMatrix !! N_" << truth_bin_j << "(" << truth_bin_min << " < PtTrue < " << truth_bin_max
                  << ") = " << auxRDF_MCGen.Count().GetValue() << std::endl;

        for (Int_t reco_bin_i = 0; reco_bin_i < n_bins; reco_bin_i++) {
            reco_bin_min = pt_min + reco_bin_i * pt_step;
            reco_bin_max = reco_bin_min + pt_step;

            auto auxRDF_MCRec = feRDF_Sexaquarks                                                              //
                                    .Filter(Form("%f < PtTrue && PtTrue < %f && %f < PtReco && PtReco < %f",  //
                                                 truth_bin_min, truth_bin_max, reco_bin_min, reco_bin_max));

            std::cout << "!! Trees_RDF_ResponseMatrix !! >> N_" << reco_bin_i << "_" << truth_bin_j << "(" << reco_bin_min << " < PtReco < "
                      << reco_bin_max << ", " << truth_bin_min << " < TruePt < " << truth_bin_max << ") = " << auxRDF_MCRec.Count().GetValue()
                      << std::endl;
        }
    }

    /*** Output ***/

    /*
        TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");

        HistFound_PtReco->Write();
        HistInj_PtReco->Write();
        HistEfficiency_Pt->Write();
    */

    /* Print how many event-loops were executed */

    std::cout << "!! Trees_RDF_ResponseMatrix !! NRuns !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! ===== !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! RDF_Injected     = " << RDF_Injected.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! eRDF_Injected    = " << eRDF_Injected.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! RDF_Sexaquarks   = " << RDF_Sexaquarks.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! fRDF_Sexaquarks  = " << fRDF_Sexaquarks.GetNRuns() << " !!" << std::endl;
    std::cout << "!! Trees_RDF_ResponseMatrix !! feRDF_Sexaquarks = " << feRDF_Sexaquarks.GetNRuns() << " !!" << std::endl;

    std::cout << "!! Trees_RDF_ResponseMatrix !! Finished !!" << std::endl;
}
