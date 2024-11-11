#include "include/Headers.hxx"
#include "include/Style.cxx"

using namespace ROOT;

/*
 *
 */
void Trees_TPCMaps() {

    TString InputFileName = "../sexaquark/attempts/local_signalMC_18qr_test_kalmanA1.94/AnalysisResults.root";

    RDataFrame RDF_Tracks("Tracks", InputFileName);
    RDataFrame RDF_MCParticles("MCParticles", InputFileName);

    std::cout << "N_RDF_Tracks: " << RDF_Tracks.Count().GetValue() << std::endl;

    /*  */

    auto bkgRDF_Tracks = RDF_Tracks.Filter("!IsSignal");
    std::cout << "N_bkgRDF_Tracks: " << bkgRDF_Tracks.Count().GetValue() << std::endl;

    auto signalRDF_Tracks = RDF_Tracks.Filter("IsSignal");
    std::cout << "N_signalRDF_Tracks: " << signalRDF_Tracks.Count().GetValue() << std::endl;

    /*  */

    std::vector<Int_t> MCParticles_RunNumber = RDF_MCParticles.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> MCParticles_DirNumber = RDF_MCParticles.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> MCParticles_EventNumber = RDF_MCParticles.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> MCParticles_Idx = RDF_MCParticles.Take<Int_t>("Idx").GetValue();
    std::vector<Float_t> MCParticles_Xv_i = RDF_MCParticles.Take<Float_t>("Xv_i").GetValue();
    std::vector<Float_t> MCParticles_Yv_i = RDF_MCParticles.Take<Float_t>("Yv_i").GetValue();
    const Int_t N_MCParticles = (Int_t)MCParticles_RunNumber.size();

    /*  */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_RadiusTrue;
    for (Int_t i = 0; i < N_MCParticles; i++) {
        Map_RadiusTrue[std::make_tuple(MCParticles_RunNumber[i], MCParticles_DirNumber[i], MCParticles_EventNumber[i], MCParticles_Idx[i])] =
            (Float_t)TMath::Sqrt((Double_t)MCParticles_Xv_i[i] * (Double_t)MCParticles_Xv_i[i] +
                                 (Double_t)MCParticles_Yv_i[i] * (Double_t)MCParticles_Yv_i[i]);
    }

    auto Lambda_GetRadiusTrue = [&Map_RadiusTrue](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t IdxMC) -> Float_t {
        return Map_RadiusTrue[std::make_tuple(RunNumber, DirNumber, EventNumber, IdxMC)];
    };

    /*  */

    auto feRDF_Tracks = bkgRDF_Tracks.Define("RadiusTrue", Lambda_GetRadiusTrue, {"RunNumber", "DirNumber", "EventNumber", "Idx_True"});
    auto fefRDF_Tracks = feRDF_Tracks.Filter("RadiusTrue > 85.");
    std::cout << "N_fefRDF_Tracks: " << fefRDF_Tracks.Count().GetValue() << std::endl;

    std::vector<Float_t> Tracks_RadiusTrue = fefRDF_Tracks.Take<Float_t>("RadiusTrue").GetValue();
    std::vector<TBits> Tracks_TPCFitMap = fefRDF_Tracks.Take<TBits>("TPCFitMap").GetValue();
    std::vector<TBits> Tracks_TPCClusterMap = fefRDF_Tracks.Take<TBits>("TPCClusterMap").GetValue();
    const Int_t N_Tracks = (Int_t)Tracks_RadiusTrue.size();

    /*  */

    TH2F *h_test1 = new TH2F("Radius_vs_FH_TPCFitMap", "Radius_vs_FH_TPCFitMap", 159, 0, 159, 200, 0, 200);
    TH2F *h_test2 = new TH2F("Radius_vs_FH_TPCClusterMap", "Radius_vs_FH_TPCClusterMap", 159, 0, 159, 200, 0, 200);

    for (Int_t i = 0; i < N_Tracks; i++) {
        h_test1->Fill(Tracks_TPCFitMap[i].FirstSetBit(), Tracks_RadiusTrue[i]);
        h_test2->Fill(Tracks_TPCClusterMap[i].FirstSetBit(), Tracks_RadiusTrue[i]);
    }

    /*  */

    TCanvas *c = new TCanvas("c", "c", 1080, 1080);

    SetMy2DStyle();

    h_test1->Draw("COLZ");
    c->Print(Form("gfx/%s.png", h_test1->GetName()));
    c->Clear();

    h_test2->Draw("COLZ");
    c->Print(Form("gfx/%s.png", h_test2->GetName()));
}
