#include "include/Headers.hxx"
#include "include/Style.hxx"

using namespace ROOT;

void Trees_RDF_GetEfficiency(TString InputFileName = "../output/AnalysisResults_A1.8_local.root",         //
                             TString OutputFileName = "output_rdf/EfficiencyResults_A1.8.root",           //
                             TString Path_PhotonConvRadius = "../radius-weights/output_bkg/merged.root",  //
                             TString Path_BlastWave = "blastwave.root") {

    /*** Input ***/

    /* Input: Get trees */

    RDataFrame RDF_Events("Events", InputFileName);
    RDataFrame RDF_MCParticles("MCParticles", InputFileName);
    RDataFrame RDF_Injected("Injected", InputFileName);
    RDataFrame RDF_Sexaquarks("Sexaquarks", InputFileName);

    const Float_t SexaquarkMass = RDF_Injected.Take<Float_t>("Mass").GetValue()[0];

    /* Filter trees (1) */

    auto fRDF_MCParticles = RDF_MCParticles.Filter("IsSignal").Filter("Idx_Mother == -1");
    auto fRDF_Sexaquarks = RDF_Sexaquarks.Filter("IsSignal");

    /* Input: Get Radius weights */

    TFile* File_PhotonConvRadius = TFile::Open(Path_PhotonConvRadius, "READ");
    if (!File_PhotonConvRadius || File_PhotonConvRadius->IsZombie()) {
        std::cout << "!! ERROR !! Trees_RDF_GetEfficiency !! Couldn't open TFile " << Path_PhotonConvRadius << " !!" << std::endl;
        return;
    }

    TString Name_Hist_PhotonConv = "MCGen_PhotonConversions_Radius";
    TH1F* Hist_PhotonConv = dynamic_cast<TH1F*>(File_PhotonConvRadius->Get("Hists")->FindObject(Name_Hist_PhotonConv));
    if (!Hist_PhotonConv) {
        std::cout << "!! ERROR !! Trees_RDF_GetEfficiency !! Couldn't find TH1F " << Name_Hist_PhotonConv << " in TFile " << Path_PhotonConvRadius
                  << " !!" << std::endl;
        return;
    }

    // clone the histogram, so I can close the file without any problem
    TH1F* Hist_RadiusWeights = dynamic_cast<TH1F*>(Hist_PhotonConv->Clone());
    Hist_RadiusWeights->SetDirectory(0);                             // detach from current directory
    Hist_RadiusWeights->Scale(1. / Hist_RadiusWeights->Integral());  // normalize it

    File_PhotonConvRadius->Close();

    auto Lambda_GetRadiusWeight = [&Hist_RadiusWeights](Float_t RadiusTrue) {  //
        return Hist_RadiusWeights->GetBinContent(Hist_RadiusWeights->FindBin(RadiusTrue));
    };

    /* Input: Get (Pt, Centrality) weights */

    TFile* File_BlastWave = TFile::Open(Path_BlastWave, "READ");
    const Int_t NCentralityClasses = 10;
    TString CentralityClasses[NCentralityClasses] = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"};
    std::vector<TH1D*> HistsVec_PtWeights;
    TString Name_Hist_BlastWave;

    for (Int_t cc = 0; cc < NCentralityClasses; cc++) {
        Name_Hist_BlastWave = Form("Hist_BlastWave_%.2f_%s", SexaquarkMass, CentralityClasses[cc].Data());
        TH1D* Hist_BlastWave = dynamic_cast<TH1D*>(File_BlastWave->Get(Name_Hist_BlastWave));
        if (!Hist_BlastWave) {
            std::cout << "!! ERROR !! Trees_RDF_GetEfficiency !! Couldn't find TH1D " << Name_Hist_BlastWave << " in TFile " << Path_BlastWave
                      << " !!" << std::endl;
            return;
        }
        // clone the histograms, so I can close the file without any problem
        TH1D* AuxHist_PtWeights = dynamic_cast<TH1D*>(Hist_BlastWave->Clone());
        AuxHist_PtWeights->SetDirectory(0);                            // detach from current directory
        AuxHist_PtWeights->Scale(1. / AuxHist_PtWeights->Integral());  // normalize
        HistsVec_PtWeights.push_back(AuxHist_PtWeights);
    }

    File_BlastWave->Close();

    auto Lambda_GetCentralityBin = [](Float_t Centrality) {  //
        return Centrality >= 5. ? static_cast<Int_t>(TMath::Floor(0.1 * Centrality + 1)) : 0;
    };

    auto Lambda_GetPtWeight = [&HistsVec_PtWeights](Float_t PtTrue, Int_t CentralityBin) {
        return HistsVec_PtWeights[CentralityBin]->GetBinContent(HistsVec_PtWeights[CentralityBin]->FindBin(PtTrue));
    };

    /*** Preparation ***/

    /* Preparation: same-size vectors */

    std::vector<Int_t> Events_RunNumber = RDF_Events.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> Events_DirNumber = RDF_Events.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> Events_EventNumber = RDF_Events.Take<Int_t>("EventNumber").GetValue();
    std::vector<Float_t> Events_Centrality = RDF_Events.Take<Float_t>("Centrality").GetValue();
    const Int_t N_Events = (Int_t)Events_RunNumber.size();

    std::vector<Int_t> MCParticles_RunNumber = fRDF_MCParticles.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> MCParticles_DirNumber = fRDF_MCParticles.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> MCParticles_EventNumber = fRDF_MCParticles.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> MCParticles_ReactionID = fRDF_MCParticles.Take<Int_t>("ReactionID").GetValue();
    std::vector<Float_t> MCParticles_Xv_i = fRDF_MCParticles.Take<Float_t>("Xv_i").GetValue();
    std::vector<Float_t> MCParticles_Yv_i = fRDF_MCParticles.Take<Float_t>("Yv_i").GetValue();
    const Int_t N_MCParticles = (Int_t)MCParticles_RunNumber.size();

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
    std::vector<Float_t> Sexaquark_SV_Xv = fRDF_Sexaquarks.Take<Float_t>("SV_Xv").GetValue();
    std::vector<Float_t> Sexaquark_SV_Yv = fRDF_Sexaquarks.Take<Float_t>("SV_Yv").GetValue();
    const Int_t N_Sexaquarks = (Int_t)Sexaquark_ReactionID.size();

    /* Preparation: Centrality maps */

    std::map<std::tuple<Int_t, Int_t, Int_t>, Float_t> Map_Centrality;
    for (Int_t i = 0; i < N_Events; i++) {
        Map_Centrality[std::make_tuple(Events_RunNumber[i], Events_DirNumber[i], Events_EventNumber[i])] = Events_Centrality[i];
    }

    auto Lambda_GetCentrality = [&Map_Centrality](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber) {
        return Map_Centrality[std::make_tuple(RunNumber, DirNumber, EventNumber)];
    };

    /* Preparation: True Radius maps */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_RadiusTrue;
    for (Int_t i = 0; i < N_MCParticles; i++) {
        Map_RadiusTrue[std::make_tuple(MCParticles_RunNumber[i], MCParticles_DirNumber[i], MCParticles_EventNumber[i], MCParticles_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)MCParticles_Xv_i[i] * (Double_t)MCParticles_Xv_i[i] +
                                 (Double_t)MCParticles_Yv_i[i] * (Double_t)MCParticles_Yv_i[i]);
    }

    auto Lambda_GetRadiusTrue = [&Map_RadiusTrue](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Map_RadiusTrue[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Preparation: True Pt maps */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_PtTrue;
    for (Int_t i = 0; i < N_Injected; i++) {
        Map_PtTrue[std::make_tuple(Injected_RunNumber[i], Injected_DirNumber[i], Injected_EventNumber[i], Injected_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)Injected_Px[i] * (Double_t)Injected_Px[i] + (Double_t)Injected_Py[i] * (Double_t)Injected_Py[i]);
    }

    auto Lambda_GetPtTrue = [&Map_PtTrue](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Map_PtTrue[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Preparation: Reconstructed Pt and Radius mapss */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_PtReco;
    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_RadiusReco;
    for (Int_t i = 0; i < N_Sexaquarks; i++) {
        Map_PtReco[std::make_tuple(Sexaquark_RunNumber[i], Sexaquark_DirNumber[i], Sexaquark_EventNumber[i], Sexaquark_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)Sexaquark_Px[i] * (Double_t)Sexaquark_Px[i] + (Double_t)Sexaquark_Py[i] * (Double_t)Sexaquark_Py[i]);
        Map_RadiusReco[std::make_tuple(Sexaquark_RunNumber[i], Sexaquark_DirNumber[i], Sexaquark_EventNumber[i], Sexaquark_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)Sexaquark_SV_Xv[i] * (Double_t)Sexaquark_SV_Xv[i] +
                                 (Double_t)Sexaquark_SV_Yv[i] * (Double_t)Sexaquark_SV_Yv[i]);
    }

    auto Lambda_GetPtReco = [&Map_PtReco](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Map_PtReco[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };
    auto Lambda_GetRadiusReco = [&Map_RadiusReco](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Map_RadiusReco[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Extend trees */

    auto eRDF_Injected = RDF_Injected
                             .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})                //
                             .Define("CentralityBin", Lambda_GetCentralityBin, {"Centrality"})                                     //
                             .Define("PtTrue", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))")                             //
                             .Define("PtReco", Lambda_GetPtReco, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})          //
                             .Define("PtWeights", Lambda_GetPtWeight, {"PtTrue", "CentralityBin"})                                 //
                             .Define("RadiusTrue", Lambda_GetRadiusTrue, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})  //
                             .Define("RadiusReco", Lambda_GetRadiusReco, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})  //
                             .Define("RadiusWeights", Lambda_GetRadiusWeight, {"RadiusTrue"})                                      //
                             .Define("BothWeights", "PtWeights * RadiusWeights");

    auto feRDF_Sexaquarks = fRDF_Sexaquarks
                                .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})                //
                                .Define("CentralityBin", Lambda_GetCentralityBin, {"Centrality"})                                     //
                                .Define("PtReco", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))")                             //
                                .Define("PtTrue", Lambda_GetPtTrue, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})          //
                                .Define("PtWeights", Lambda_GetPtWeight, {"PtTrue", "CentralityBin"})                                 //
                                .Define("RadiusReco", "TMath::Sqrt(SV_Xv * SV_Xv + SV_Yv * SV_Yv)")                                   //
                                .Define("RadiusTrue", Lambda_GetRadiusTrue, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})  //
                                .Define("RadiusWeights", Lambda_GetRadiusWeight, {"RadiusTrue"})                                      //
                                .Define("BothWeights", "PtWeights * RadiusWeights");

    /* Filter trees (2) */

    auto fRDF_Injected = eRDF_Injected.Filter("Centrality < 90.");
    auto fefRDF_Sexaquarks = feRDF_Sexaquarks.Filter("Centrality < 90.").Filter("IsSignal");

    /*** Histograms ***/

    /* Histograms: Distributions (with and without weights) */

    auto HistFound_PtReco = fefRDF_Sexaquarks.Histo1D({"Found_PtReco", ";;", 100, 0., 10.}, "PtReco");
    auto HistFound_PtReco_wPtWeights = fefRDF_Sexaquarks.Histo1D({"Found_PtReco_wPtWeights", ";;", 100, 0., 10.}, "PtReco", "PtWeights");
    auto HistFound_PtReco_wRadiusWeights = fefRDF_Sexaquarks.Histo1D({"Found_PtReco_wRadiusWeights", ";;", 100, 0., 10.}, "PtReco", "RadiusWeights");
    auto HistFound_PtReco_wBothWeights = fefRDF_Sexaquarks.Histo1D({"Found_PtReco_wBothWeights", ";;", 100, 0., 10.}, "PtReco", "BothWeights");

    auto HistInj_PtReco = fRDF_Injected.Histo1D({"Injected_PtReco", ";;", 100, 0., 10.}, "PtReco");
    auto HistInj_PtReco_wPtWeights = fRDF_Injected.Histo1D({"Injected_PtReco_wPtWeights", ";;", 100, 0., 10.}, "PtReco", "PtWeights");
    auto HistInj_PtReco_wRadiusWeights = fRDF_Injected.Histo1D({"Injected_PtReco_wRadiusWeights", ";;", 100, 0., 10.}, "PtReco", "RadiusWeights");
    auto HistInj_PtReco_wBothWeights = fRDF_Injected.Histo1D({"Injected_PtReco_wBothWeights", ";;", 100, 0., 10.}, "PtReco", "BothWeights");

    auto HistFound_RadiusReco = fefRDF_Sexaquarks.Histo1D({"Found_RadiusReco", ";;", 100, 0., 200.}, "RadiusReco");
    auto HistFound_RadiusReco_wPtWeights = fefRDF_Sexaquarks.Histo1D({"Found_RadiusReco_wPtWeights", ";;", 100, 0., 200.}, "RadiusReco", "PtWeights");
    auto HistFound_RadiusReco_wRadiusWeights =
        fefRDF_Sexaquarks.Histo1D({"Found_RadiusReco_wRadiusWeights", ";;", 100, 0., 200.}, "RadiusReco", "RadiusWeights");
    auto HistFound_RadiusReco_wBothWeights =
        fefRDF_Sexaquarks.Histo1D({"Found_RadiusReco_wBothWeights", ";;", 100, 0., 200.}, "RadiusReco", "BothWeights");

    auto HistInj_RadiusReco = fRDF_Injected.Histo1D({"Injected_RadiusReco", ";;", 100, 0., 200.}, "RadiusReco");
    auto HistInj_RadiusReco_wPtWeights = fRDF_Injected.Histo1D({"Injected_RadiusReco_wPtWeights", ";;", 100, 0., 200.}, "RadiusReco", "PtWeights");
    auto HistInj_RadiusReco_wRadiusWeights =
        fRDF_Injected.Histo1D({"Injected_RadiusReco_wRadiusWeights", ";;", 100, 0., 200.}, "RadiusReco", "RadiusWeights");
    auto HistInj_RadiusReco_wBothWeights =
        fRDF_Injected.Histo1D({"Injected_RadiusReco_wBothWeights", ";;", 100, 0., 200.}, "RadiusReco", "BothWeights");

    /* Histograms: Efficiency (with and without weights) */

    TH1D* HistEfficiency_Pt = new TH1D("Efficiency_Pt", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    HistEfficiency_Pt->Divide(HistFound_PtReco.GetPtr(), HistInj_PtReco.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Pt_wPtWeights = new TH1D("Efficiency_Pt_wPtWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    HistEfficiency_Pt_wPtWeights->Divide(HistFound_PtReco_wPtWeights.GetPtr(), HistInj_PtReco_wPtWeights.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Pt_wRadiusWeights = new TH1D("Efficiency_Pt_wRadiusWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    HistEfficiency_Pt_wRadiusWeights->Divide(HistFound_PtReco_wRadiusWeights.GetPtr(), HistInj_PtReco_wRadiusWeights.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Pt_wBothWeights = new TH1D("Efficiency_Pt_wBothWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    HistEfficiency_Pt_wBothWeights->Divide(HistFound_PtReco_wBothWeights.GetPtr(), HistInj_PtReco_wBothWeights.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Radius = new TH1D("Efficiency_Radius", ";Radius (cm);Efficiency", 100, 0., 200.);
    HistEfficiency_Radius->Divide(HistFound_RadiusReco.GetPtr(), HistInj_RadiusReco.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Radius_wPtWeights = new TH1D("Efficiency_Radius_wPtWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    HistEfficiency_Radius_wPtWeights->Divide(HistFound_RadiusReco_wPtWeights.GetPtr(), HistInj_RadiusReco_wPtWeights.GetPtr(), 1., 1., "B");

    TH1D* HistEfficiency_Radius_wRadiusWeights = new TH1D("Efficiency_Radius_wRadiusWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    HistEfficiency_Radius_wRadiusWeights->Divide(HistFound_RadiusReco_wRadiusWeights.GetPtr(), HistInj_RadiusReco_wRadiusWeights.GetPtr(), 1., 1.,
                                                 "B");

    TH1D* HistEfficiency_Radius_wBothWeights = new TH1D("Efficiency_Radius_wBothWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    HistEfficiency_Radius_wBothWeights->Divide(HistFound_RadiusReco_wBothWeights.GetPtr(), HistInj_RadiusReco_wBothWeights.GetPtr(), 1., 1., "B");

    /* Output */

    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");

    HistFound_PtReco->Write();
    HistFound_PtReco_wPtWeights->Write();
    HistFound_PtReco_wRadiusWeights->Write();
    HistFound_PtReco_wBothWeights->Write();

    HistInj_PtReco->Write();
    HistInj_PtReco_wPtWeights->Write();
    HistInj_PtReco_wRadiusWeights->Write();
    HistInj_PtReco_wBothWeights->Write();

    HistFound_RadiusReco->Write();
    HistFound_RadiusReco_wPtWeights->Write();
    HistFound_RadiusReco_wRadiusWeights->Write();
    HistFound_RadiusReco_wBothWeights->Write();

    HistInj_RadiusReco->Write();
    HistInj_RadiusReco_wPtWeights->Write();
    HistInj_RadiusReco_wRadiusWeights->Write();
    HistInj_RadiusReco_wBothWeights->Write();

    HistEfficiency_Pt->Write();
    HistEfficiency_Pt_wPtWeights->Write();
    HistEfficiency_Pt_wRadiusWeights->Write();
    HistEfficiency_Pt_wBothWeights->Write();

    HistEfficiency_Radius->Write();
    HistEfficiency_Radius_wPtWeights->Write();
    HistEfficiency_Radius_wRadiusWeights->Write();
    HistEfficiency_Radius_wBothWeights->Write();

    /* Print how many event-loops were executed */

    std::cout << "NRuns" << std::endl;
    std::cout << "=====" << std::endl;
    std::cout << "RDF_Events        = " << RDF_Events.GetNRuns() << std::endl;
    std::cout << "RDF_MCParticles   = " << RDF_MCParticles.GetNRuns() << std::endl;
    std::cout << "fRDF_MCParticles  = " << fRDF_MCParticles.GetNRuns() << std::endl;
    std::cout << "RDF_Injected      = " << RDF_Injected.GetNRuns() << std::endl;
    std::cout << "eRDF_Injected     = " << eRDF_Injected.GetNRuns() << std::endl;
    std::cout << "RDF_Sexaquarks    = " << RDF_Sexaquarks.GetNRuns() << std::endl;
    std::cout << "fRDF_Sexaquarks   = " << fRDF_Sexaquarks.GetNRuns() << std::endl;
    std::cout << "feRDF_Sexaquarks  = " << feRDF_Sexaquarks.GetNRuns() << std::endl;
    std::cout << "fefRDF_Sexaquarks = " << fefRDF_Sexaquarks.GetNRuns() << std::endl;
}
