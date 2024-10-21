#include "include/Headers.hxx"
#include "include/Style.hxx"

using namespace ROOT;

void Trees_RDF_GetEfficiency(TString InputFileName = "../output/AnalysisResults_A1.73_local.root",        //
                             TString OutputFileName = "output_rdf/EfficiencyResults_A1.73.root",          //
                             TString Path_PhotonConvRadius = "../radius-weights/output_bkg/merged.root",  //
                             TString Path_BlastWave = "blastwave_bt.root") {

    /*** Input ***/

    /* Input: Get trees */

    RDataFrame RDF_Events("Events", InputFileName);
    RDataFrame RDF_MCParticles("MCParticles", InputFileName);
    RDataFrame RDF_Injected("Injected", InputFileName);
    RDataFrame RDF_Sexaquarks("Sexaquarks", InputFileName);

    auto fRDF_MCParticles = RDF_MCParticles.Filter("IsSignal").Filter("Idx_Mother == -1");
    const Float_t SexaquarkMass = RDF_Injected.Take<Float_t>("Mass").GetValue()[0];

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

    auto Lambda_GetRadiusWeight = [&Hist_RadiusWeights](Float_t TrueRadius) {  //
        return Hist_RadiusWeights->GetBinContent(Hist_RadiusWeights->FindBin(TrueRadius));
    };

    /* Input: Get (Pt, Centrality) weights */

    TFile* File_BlastWave = TFile::Open(Path_BlastWave, "READ");
    std::vector<TString> CentralityBins = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"};
    std::vector<TH1D*> HistsVec_PtWeights;
    TString Name_Hist_BlastWave;

    for (Int_t cc = 0; cc < (Int_t)CentralityBins.size(); cc++) {
        Name_Hist_BlastWave = Form("BlastWaveHist_%.2f_%s", SexaquarkMass, CentralityBins[cc].Data());
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

    auto Lambda_GetPtWeight = [&HistsVec_PtWeights](Float_t TruePt, Int_t CentralityBin) {
        return HistsVec_PtWeights[CentralityBin]->GetBinContent(HistsVec_PtWeights[CentralityBin]->FindBin(TruePt));
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

    /* Preparation: Centrality hash table */

    std::map<std::tuple<Int_t, Int_t, Int_t>, Float_t> Hash_Centrality;
    for (Int_t i = 0; i < N_Events; i++) {
        Hash_Centrality[std::make_tuple(Events_RunNumber[i], Events_DirNumber[i], Events_EventNumber[i])] = Events_Centrality[i];
    }

    auto Lambda_GetCentrality = [&Hash_Centrality](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber) {
        return Hash_Centrality[std::make_tuple(RunNumber, DirNumber, EventNumber)];
    };

    /* Preparation: Radius hash table */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Hash_Radius;
    for (Int_t i = 0; i < N_MCParticles; i++) {
        Hash_Radius[std::make_tuple(MCParticles_RunNumber[i], MCParticles_DirNumber[i], MCParticles_EventNumber[i], MCParticles_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)MCParticles_Xv_i[i] * (Double_t)MCParticles_Xv_i[i] +
                                 (Double_t)MCParticles_Yv_i[i] * (Double_t)MCParticles_Yv_i[i]);
    }

    auto Lambda_GetRadius = [&Hash_Radius](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Hash_Radius[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Preparation: true Pt has table */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Hash_TruePt;
    for (Int_t i = 0; i < N_Injected; i++) {
        Hash_TruePt[std::make_tuple(Injected_RunNumber[i], Injected_DirNumber[i], Injected_EventNumber[i], Injected_ReactionID[i])] =
            (Float_t)TMath::Sqrt((Double_t)Injected_Px[i] * (Double_t)Injected_Px[i] + (Double_t)Injected_Py[i] * (Double_t)Injected_Py[i]);
    }

    auto Lambda_GetTruePt = [&Hash_TruePt](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t ReactionID) {
        return Hash_TruePt[std::make_tuple(RunNumber, DirNumber, EventNumber, ReactionID)];
    };

    /* Extend trees */

    auto eRDF_Injected = RDF_Injected
                             .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})        //
                             .Define("CentralityBin", Lambda_GetCentralityBin, {"Centrality"})                             //
                             .Define("Pt", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))")                         //
                             .Define("PtWeights", Lambda_GetPtWeight, {"Pt", "CentralityBin"})                             //
                             .Define("Radius", Lambda_GetRadius, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})  //
                             .Define("RadiusWeights", Lambda_GetRadiusWeight, {"Radius"})                                  //
                             .Define("BothWeights", "PtWeights * RadiusWeights");

    auto eRDF_Sexaquarks = RDF_Sexaquarks
                               .Define("Centrality", Lambda_GetCentrality, {"RunNumber", "DirNumber", "EventNumber"})            //
                               .Define("CentralityBin", Lambda_GetCentralityBin, {"Centrality"})                                 //
                               .Define("Pt", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))")                             //
                               .Define("TruePt", Lambda_GetTruePt, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})      //
                               .Define("PtWeights", Lambda_GetPtWeight, {"TruePt", "CentralityBin"})                             //
                               .Define("Radius", "TMath::Sqrt(SV_Xv * SV_Xv + SV_Yv * SV_Yv)")                                   //
                               .Define("TrueRadius", Lambda_GetRadius, {"RunNumber", "DirNumber", "EventNumber", "ReactionID"})  //
                               .Define("RadiusWeights", Lambda_GetRadiusWeight, {"TrueRadius"})                                  //
                               .Define("BothWeights", "PtWeights * RadiusWeights");

    /* Filter trees */

    auto fRDF_Injected = eRDF_Injected.Filter("Centrality < 90.");
    auto fRDF_Sexaquarks = eRDF_Sexaquarks.Filter("Centrality < 90.").Filter("IsSignal");

    /*** Histograms ***/

    /* Histograms: Distributions (with and without weights) */

    auto Hist_TruePt = fRDF_Injected.Histo1D({"TruePt", ";;", 100, 0., 10.}, "Pt");
    auto Hist_TruePt_wPtWeights = fRDF_Injected.Histo1D({"TruePt_wPtWeights", ";;", 100, 0., 10.}, "Pt", "PtWeights");
    auto Hist_TruePt_wRadiusWeights = fRDF_Injected.Histo1D({"TruePt_wRadiusWeights", ";;", 100, 0., 10.}, "Pt", "RadiusWeights");
    auto Hist_TruePt_wBothWeights = fRDF_Injected.Histo1D({"TruePt_wBothWeights", ";;", 100, 0., 10.}, "Pt", "BothWeights");

    auto Hist_RecPt = fRDF_Sexaquarks.Histo1D({"RecPt", ";;", 100, 0., 10.}, "Pt");
    auto Hist_RecPt_wPtWeights = fRDF_Sexaquarks.Histo1D({"RecPt_wPtWeights", ";;", 100, 0., 10.}, "Pt", "PtWeights");
    auto Hist_RecPt_wRadiusWeights = fRDF_Sexaquarks.Histo1D({"RecPt_wRadiusWeights", ";;", 100, 0., 10.}, "Pt", "RadiusWeights");
    auto Hist_RecPt_wBothWeights = fRDF_Sexaquarks.Histo1D({"RecPt_wBothWeights", ";;", 100, 0., 10.}, "Pt", "BothWeights");

    auto Hist_TrueRadius = fRDF_Injected.Histo1D({"TrueRadius", ";;", 100, 0., 200.}, "Radius");
    auto Hist_TrueRadius_wPtWeights = fRDF_Injected.Histo1D({"TrueRadius_wPtWeights", ";;", 100, 0., 200.}, "Radius", "PtWeights");
    auto Hist_TrueRadius_wRadiusWeights = fRDF_Injected.Histo1D({"TrueRadius_wRadiusWeights", ";;", 100, 0., 200.}, "Radius", "RadiusWeights");
    auto Hist_TrueRadius_wBothWeights = fRDF_Injected.Histo1D({"TrueRadius_wBothWeights", ";;", 100, 0., 200.}, "Radius", "BothWeights");

    auto Hist_RecRadius = fRDF_Sexaquarks.Histo1D({"RecRadius", ";;", 100, 0., 200.}, "Radius");
    auto Hist_RecRadius_wPtWeights = fRDF_Sexaquarks.Histo1D({"RecRadius_wPtWeights", ";;", 100, 0., 200.}, "Radius", "PtWeights");
    auto Hist_RecRadius_wRadiusWeights = fRDF_Sexaquarks.Histo1D({"RecRadius_wRadiusWeights", ";;", 100, 0., 200.}, "Radius", "RadiusWeights");
    auto Hist_RecRadius_wBothWeights = fRDF_Sexaquarks.Histo1D({"RecRadius_wBothWeights", ";;", 100, 0., 200.}, "Radius", "BothWeights");

    /* Histograms: Efficiency (with and without weights) */

    TH1D* Hist_EfficiencyPt = new TH1D("EfficiencyPt", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    Hist_EfficiencyPt->Divide(Hist_RecPt.GetPtr(), Hist_TruePt.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyPt_wPtWeights = new TH1D("EfficiencyPt_wPtWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    Hist_EfficiencyPt_wPtWeights->Divide(Hist_RecPt_wPtWeights.GetPtr(), Hist_TruePt_wPtWeights.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyPt_wRadiusWeights = new TH1D("EfficiencyPt_wRadiusWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    Hist_EfficiencyPt_wRadiusWeights->Divide(Hist_RecPt_wRadiusWeights.GetPtr(), Hist_TruePt_wRadiusWeights.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyPt_wBothWeights = new TH1D("EfficiencyPt_wBothWeights", ";p_{T} (GeV/c);Efficiency", 100, 0., 10.);
    Hist_EfficiencyPt_wBothWeights->Divide(Hist_RecPt_wBothWeights.GetPtr(), Hist_TruePt_wBothWeights.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyRadius = new TH1D("EfficiencyRadius", ";Radius (cm);Efficiency", 100, 0., 200.);
    Hist_EfficiencyRadius->Divide(Hist_RecRadius.GetPtr(), Hist_TrueRadius.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyRadius_wPtWeights = new TH1D("EfficiencyRadius_wPtWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    Hist_EfficiencyRadius_wPtWeights->Divide(Hist_RecRadius_wPtWeights.GetPtr(), Hist_TrueRadius_wPtWeights.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyRadius_wRadiusWeights = new TH1D("EfficiencyRadius_wRadiusWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    Hist_EfficiencyRadius_wRadiusWeights->Divide(Hist_RecRadius_wRadiusWeights.GetPtr(), Hist_TrueRadius_wRadiusWeights.GetPtr(), 1., 1., "B");

    TH1D* Hist_EfficiencyRadius_wBothWeights = new TH1D("EfficiencyRadius_wBothWeights", ";Radius (cm);Efficiency", 100, 0., 200.);
    Hist_EfficiencyRadius_wBothWeights->Divide(Hist_RecRadius_wBothWeights.GetPtr(), Hist_TrueRadius_wBothWeights.GetPtr(), 1., 1., "B");

    /* Output */

    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");

    Hist_TruePt->Write();
    Hist_TruePt_wPtWeights->Write();
    Hist_TruePt_wRadiusWeights->Write();
    Hist_TruePt_wBothWeights->Write();

    Hist_RecPt->Write();
    Hist_RecPt_wPtWeights->Write();
    Hist_RecPt_wRadiusWeights->Write();
    Hist_RecPt_wBothWeights->Write();

    Hist_TrueRadius->Write();
    Hist_TrueRadius_wPtWeights->Write();
    Hist_TrueRadius_wRadiusWeights->Write();
    Hist_TrueRadius_wBothWeights->Write();

    Hist_RecRadius->Write();
    Hist_RecRadius_wPtWeights->Write();
    Hist_RecRadius_wRadiusWeights->Write();
    Hist_RecRadius_wBothWeights->Write();

    Hist_EfficiencyPt->Write();
    Hist_EfficiencyPt_wPtWeights->Write();
    Hist_EfficiencyPt_wRadiusWeights->Write();
    Hist_EfficiencyPt_wBothWeights->Write();

    Hist_EfficiencyRadius->Write();
    Hist_EfficiencyRadius_wPtWeights->Write();
    Hist_EfficiencyRadius_wRadiusWeights->Write();
    Hist_EfficiencyRadius_wBothWeights->Write();

    /* Print how many event-loops were executed */

    std::cout << "NRuns" << std::endl;
    std::cout << "=====" << std::endl;
    std::cout << "RDF_Events        = " << RDF_Events.GetNRuns() << std::endl;
    std::cout << "fRDF_MCParticles  = " << fRDF_MCParticles.GetNRuns() << std::endl;
    std::cout << "eRDF_Injected     = " << eRDF_Injected.GetNRuns() << std::endl;
    std::cout << "eRDF_Sexaquarks   = " << eRDF_Sexaquarks.GetNRuns() << std::endl;
}
