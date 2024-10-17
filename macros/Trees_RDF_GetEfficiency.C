#include "include/Headers.hxx"
#include "include/Style.hxx"

using namespace ROOT;

const Int_t N_CORES = 48;

void Trees_RDF_GetEfficiency(TString InputFileName = "../output/A1.8_18qr_local/AnalysisResults_merged.root",  //
                             TString OutputFileName = "output_rdf/EfficiencyResults.root",                     //
                             TString PhotonConvRadiusPath = "../radius-weights/output_bkg/merged.root",        //
                             TString BlastWavePath = "blastwave_bt.root") {

    /* Input */

    EnableImplicitMT(N_CORES);

    RDataFrame RDF_Events("Events", InputFileName);
    // RDataFrame RDF_MCParticles("MCParticles", InputFileName);
    RDataFrame RDF_Injected("Injected", InputFileName);
    RDataFrame RDF_Sexaquarks("Sexaquarks", InputFileName);

    std::vector<Int_t> Events_RunNumber = RDF_Events.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Events_DirNumber = RDF_Events.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Events_EventNumber = RDF_Events.Take<Int_t>("EventNumber").GetValue();
    std::vector<Float_t> Events_Centrality = RDF_Events.Take<Float_t>("Centrality").GetValue();

    // std::vector<Int_t> MCParticles_RunNumber = RDF_MCParticles.Take<Int_t>("RunNumber").GetValue();
    // std::vector<Int_t> MCParticles_DirNumber = RDF_MCParticles.Take<Int_t>("DirNumber").GetValue();
    // std::vector<Int_t> MCParticles_EventNumber = RDF_MCParticles.Take<Int_t>("EventNumber").GetValue();
    // std::vector<Int_t> MCParticles_ReactionID = RDF_MCParticles.Take<Int_t>("ReactionID").GetValue();
    // std::vector<Float_t> MCParticles_Xv_i = RDF_MCParticles.Take<Float_t>("Xv_i").GetValue();
    // std::vector<Float_t> MCParticles_Yv_i = RDF_MCParticles.Take<Float_t>("Yv_i").GetValue();

    /* Make hash tables for Centrality, that use RunNumber,DirNumber,EventNumber as keys */

    std::map<std::tuple<Int_t, Int_t, Int_t>, Float_t> Hash_Centrality;
    for (Int_t i = 0; i < (Int_t)Events_RunNumber.size(); i++) {
        Hash_Centrality[std::make_tuple(Events_RunNumber[i], Events_DirNumber[i], Events_EventNumber[i])] = Events_Centrality[i];
    }

    auto GetCentrality = [&Hash_Centrality](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber) {
        return Hash_Centrality[std::make_tuple(RunNumber, DirNumber, EventNumber)];
    };

    /* Extend trees */

    auto eRDF_Injected = RDF_Injected
                             .Define("Pt", "TMath::Sqrt(Px*Px + Py*Py)")  //
                             .Define("Centrality", GetCentrality, {"RunNumber", "DirNumber", "EventNumber"});
    auto eRDF_Sexaquarks = RDF_Sexaquarks                                   //
                               .Define("Pt", "TMath::Sqrt(Px*Px + Py*Py)")  //
                               .Define("Radius", "TMath::Sqrt(SV_Xv*SV_Xv + SV_Yv*SV_Yv)")
                               .Define("Centrality", GetCentrality, {"RunNumber", "DirNumber", "EventNumber"});

    /* Extend trees */

    auto Hist_TruePt = eRDF_Injected.Filter("Centrality < 90.").Histo1D({"TruePt", ";;", 100, 0., 10.}, "Pt");
    auto Hist_RecPt = eRDF_Sexaquarks.Filter("Centrality < 90.").Filter("IsSignal").Histo1D({"RecPt", ";;", 100, 0., 10.}, "Pt");

    TH1D *Hist_EfficiencyPt = new TH1D("EfficiencyPt", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    Hist_EfficiencyPt->Divide(Hist_RecPt.GetPtr(), Hist_TruePt.GetPtr(), 1., 1., "B");

    auto Hist_Centrality = eRDF_Sexaquarks.Filter("Centrality < 90.").Filter("IsSignal").Histo1D({"Centrality", ";;", 10, 0., 100.}, "Centrality");

    auto Hist_RecRadius = eRDF_Sexaquarks.Filter("Centrality < 90.").Filter("IsSignal").Histo1D({"RecRadius", ";;", 100, 0., 200.}, "Radius");

    /* Output */

    TFile *OutputFile = TFile::Open(OutputFileName, "RECREATE");

    Hist_TruePt->Write();
    Hist_RecPt->Write();
    Hist_EfficiencyPt->Write();
    Hist_Centrality->Write();
    Hist_RecRadius->Write();

    /* Print how many event-loops were executed */

    std::cout << "NRuns" << std::endl;
    std::cout << "=====" << std::endl;
    std::cout << "RDF_Events      = " << RDF_Events.GetNRuns() << std::endl;
    std::cout << "eRDF_Injected   = " << eRDF_Injected.GetNRuns() << std::endl;
    std::cout << "eRDF_Sexaquarks = " << eRDF_Sexaquarks.GetNRuns() << std::endl;
}
