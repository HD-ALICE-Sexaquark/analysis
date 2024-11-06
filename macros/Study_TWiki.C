#include "include/Headers.hxx"
#include "include/Style.hxx"

using namespace ROOT;

/*
 * Trying to make sense of all the information that appears on the AliDPG TWiki (https://twiki.cern.ch/twiki/bin/view/ALICE/AliceDPG).
 */
void Study_TWiki() {
    std::cout << "!! Study_TWiki !! Started !!" << std::endl;

    /** Input **/

    TString InputFilename = "../output/local_signalMC_15o_full_kalmanA1.8/AnalysisResults_merged.root";
    TFile *InputFile = TFile::Open(InputFilename);

    /***** Trees *****/

    /*** Events ***/

    RDataFrame RDF_Events("Events", InputFilename);

    auto Hist_Centrality = RDF_Events.Histo1D("Centrality");
    auto Hist_IsGenPileup = RDF_Events.Histo1D({"IsGenPileup", ";;", 2, 0., 1.}, "IsGenPileup");
    auto Hist_IsSBCPileup = RDF_Events.Histo1D({"IsSBCPileup", ";;", 2, 0., 1.}, "IsSBCPileup");

    auto Hist_PV_NContributors = RDF_Events.Histo1D("PV_NContributors");
    auto Hist_PV_ZvErr_FromSPD = RDF_Events.Histo1D("PV_ZvErr_FromSPD");
    auto Hist_PV_ZvErr_FromTracks = RDF_Events.Histo1D("PV_ZvErr_FromTracks");
    auto Hist_PV_Zv_FromSPD = RDF_Events.Histo1D("PV_Zv_FromSPD");
    auto Hist_PV_Zv_FromTracks = RDF_Events.Histo1D("PV_Zv_FromTracks");
    auto Hist_PV_Dispersion = RDF_Events.Histo1D("PV_Dispersion");
    auto Hist_NTPCClusters = RDF_Events.Histo1D("NTPCClusters");

    /* Events: NContributors */

    auto f0_RDF_Events = RDF_Events.Filter("PV_NContributors == 0");
    auto f1_RDF_Events = RDF_Events.Filter("PV_NContributors == 1");

    std::cout << "!! Study_TWiki !! Events with No PV Rec. Contributors    = " << f0_RDF_Events.Count().GetValue() << std::endl;
    std::cout << "!! Study_TWiki !! Events with Just 1 PV Rec. Contributor = " << f1_RDF_Events.Count().GetValue() << std::endl;

    /*** MC Particles ***/

    RDataFrame RDF_MCParticles("MCParticles", InputFilename);
    auto eRDF_MCParticles = RDF_MCParticles.Define("Pt", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))");

    std::vector<Int_t> MCParticles_RunNumber = eRDF_MCParticles.Take<Int_t>("RunNumber").GetValue();
    std::vector<Int_t> MCParticles_DirNumber = eRDF_MCParticles.Take<Int_t>("DirNumber").GetValue();
    std::vector<Int_t> MCParticles_EventNumber = eRDF_MCParticles.Take<Int_t>("EventNumber").GetValue();
    std::vector<Int_t> MCParticles_Idx = eRDF_MCParticles.Take<Int_t>("Idx").GetValue();
    std::vector<Int_t> MCParticles_PdgCode = eRDF_MCParticles.Take<Int_t>("PdgCode").GetValue();
    std::vector<Float_t> MCParticles_Pt = eRDF_MCParticles.Take<Float_t>("Pt").GetValue();
    const Int_t N_MCParticles = (Int_t)MCParticles_RunNumber.size();

    /* MC Particles: maps */

    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Int_t> Map_PdgCode;
    std::map<std::tuple<Int_t, Int_t, Int_t, Int_t>, Float_t> Map_TruePt;
    for (Int_t i = 0; i < N_MCParticles; i++) {
        Map_PdgCode[std::make_tuple(MCParticles_RunNumber[i], MCParticles_DirNumber[i], MCParticles_EventNumber[i], MCParticles_Idx[i])] =
            MCParticles_PdgCode[i];
        Map_TruePt[std::make_tuple(MCParticles_RunNumber[i], MCParticles_DirNumber[i], MCParticles_EventNumber[i], MCParticles_Idx[i])] =
            MCParticles_Pt[i];
    }

    auto Lambda_GetPdgCode = [&Map_PdgCode](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t Idx) -> Int_t {
        return Map_PdgCode[std::make_tuple(RunNumber, DirNumber, EventNumber, Idx)];
    };
    auto Lambda_GetTruePt = [&Map_TruePt](Int_t RunNumber, Int_t DirNumber, Int_t EventNumber, Int_t Idx) -> Float_t {
        return Map_TruePt[std::make_tuple(RunNumber, DirNumber, EventNumber, Idx)];
    };

    /* MC Particles: hists */

    auto Hist_MCParticles_IsOOBPileup = RDF_MCParticles.Histo1D({"IsOOBPileup", ";;Counts", 2, 0., 1.}, "IsOOBPileup");

    /*** Tracks ***/

    RDataFrame RDF_Tracks("Tracks", InputFilename);
    auto fRDF_Tracks = RDF_Tracks.Filter("IsSignal");
    // auto fRDF_Tracks = RDF_Tracks  //
    //    .Define("TruePdgCode", Lambda_GetPdgCode, {"RunNumber", "DirNumber", "EventNumber", "Idx_True"})
    //    .Filter("!IsSignal && !IsSecondary && IsTrue && TruePdgCode == 211");
    auto feRDF_Tracks =
        fRDF_Tracks  //
            .Define("Pt_True", Lambda_GetTruePt, {"RunNumber", "DirNumber", "EventNumber", "Idx_True"})
            .Define("Pt_Main", "static_cast<Float_t>(TMath::Sqrt(Px * Px + Py * Py))")
            .Define("Pt_Inner", "static_cast<Float_t>(TMath::Sqrt(Px_Inner * Px_Inner + Py_Inner * Py_Inner))")
            .Define("Pt_TPCInner", "static_cast<Float_t>(TMath::Sqrt(Px_TPCInner * Px_TPCInner + Py_TPCInner * Py_TPCInner))")
            .Define("Pt_Outer", "static_cast<Float_t>(TMath::Sqrt(Px_Outer * Px_Outer + Py_Outer * Py_Outer))")
            .Define("Pt_Constrained", "static_cast<Float_t>(TMath::Sqrt(Px_Constrained * Px_Constrained + Py_Constrained * Py_Constrained))")
            .Define("Res_Pt_Main", "(Pt_Main - Pt_True) / Pt_True")
            .Define("Res_Pt_Inner", "(Pt_Inner - Pt_True) / Pt_True")
            .Define("Res_Pt_TPCInner", "(Pt_TPCInner - Pt_True) / Pt_True")
            .Define("Res_Pt_Outer", "(Pt_Outer - Pt_True) / Pt_True")
            .Define("Res_Pt_Constrained", "(Pt_Constrained - Pt_True) / Pt_True");

    /* Tracks: hists */

    auto Hist_Res_Pt_Main = feRDF_Tracks.Histo1D({"Res_Pt_Main", ";;Counts", 200, -10., 10.}, "Res_Pt_Main");
    auto Hist_Res_Pt_Inner = feRDF_Tracks.Histo1D({"Res_Pt_Inner", ";;Counts", 200, -10., 10.}, "Res_Pt_Inner");
    auto Hist_Res_Pt_TPCInner = feRDF_Tracks.Histo1D({"Res_Pt_TPCInner", ";;Counts", 200, -10., 10.}, "Res_Pt_TPCInner");
    auto Hist_Res_Pt_Outer = feRDF_Tracks.Histo1D({"Res_Pt_Outer", ";;Counts", 200, -10., 10.}, "Res_Pt_Outer");
    auto Hist_Res_Pt_Constrained = feRDF_Tracks.Histo1D({"Res_Pt_Constrained", ";;Counts", 200, -10., 10.}, "Res_Pt_Constrained");

    /***** Draw *****/

    SetMyStyle();
    // gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "c", 1080, 1080);

    TString OutputFilename;
    TH1D *AllHists[26] = {
        Hist_Centrality.GetPtr(),       Hist_IsGenPileup.GetPtr(),      Hist_IsSBCPileup.GetPtr(),         Hist_MCParticles_IsOOBPileup.GetPtr(),
        Hist_PV_NContributors.GetPtr(), Hist_PV_ZvErr_FromSPD.GetPtr(), Hist_PV_ZvErr_FromTracks.GetPtr(), Hist_PV_Zv_FromSPD.GetPtr(),
        Hist_PV_Zv_FromTracks.GetPtr(), Hist_PV_Dispersion.GetPtr(),    Hist_NTPCClusters.GetPtr(),        Hist_Res_Pt_Main.GetPtr(),
        Hist_Res_Pt_Inner.GetPtr(),     Hist_Res_Pt_TPCInner.GetPtr(),  Hist_Res_Pt_Outer.GetPtr(),        Hist_Res_Pt_Constrained.GetPtr()};

    for (auto hist : AllHists) {
        hist->SetTitle("");
        hist->SetLineWidth(4);
        hist->SetLineColor(myBlue.GetNumber());
        hist->SetFillStyle(0);

        hist->Draw();

        OutputFilename = TString::Format("gfx/TWiki_%s.png", hist->GetName());
        c->Print(OutputFilename);
        c->Clear();
    }

    std::cout << "!! Study_TWiki !! Finished !!" << std::endl;
}
