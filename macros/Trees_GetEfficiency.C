#include "include/Headers.hxx"

Int_t GetCentralityBin(Float_t CentralityValue);

// Indices
//     MajorIndex
//         Events, Injected, Sexaquarks : "RunNumber * 1000 + DirNumber"
//         MCParticles                  : "RunNumber"
//                                        (^^ custom index ^^)
//     MinorIndex
//         Events               : "EventNumber"
//         Injected, Sexaquarks : "EventNumber * 1000 + ReactionID"
//         MCParticles          : "DirNumber * 100 * 1000 + EventNumber * 1000 + Status"
//                                (^^ custom index ^^)

/*
 * Process an indexed `AnalysisResults.root` file.
 * Produce (Pt, Radius) histograms for (true, reconstructed) anti-sexaquarks,
 * as well as their (unweighted, weighted) efficiency distributions.
 */
void Trees_GetEfficiency(TString InputFileName = "AnalysisResults_indexed.root",                                                    //
                         TString PhotonConvRadiusPath = "/home/ceres/borquez/work/analysis/radius-weights/output_bkg/merged.root",  //
                         TString BlastWavePath = "/home/ceres/borquez/work/analysis/macros/blastwave_bt.root") {

    /** Input **/

    TFile* InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << InputFileName << " !!" << std::endl;
        return;
    }

    std::map<TString, TTree*> Tree;
    std::vector<TString> TreeNames = {"Injected", "Events", "MCParticles", "Sexaquarks"};

    for (const auto& treeName : TreeNames) {
        Tree[treeName] = (TTree*)InputFile->Get(treeName);
        if (!Tree[treeName]) {
            std::cout << "!! ERROR !! Couldn't find TTree " << treeName << " in " << InputFileName << " !!" << std::endl;
            InputFile->Close();
            return;
        }
    }

    /** Set Branches **/

    /* Injected */

    Int_t Injected_RunNumber, Injected_DirNumber, Injected_EventNumber, Injected_ReactionID;
    Float_t Injected_Px, Injected_Py, Injected_Pz, Injected_M;

    Tree["Injected"]->SetBranchAddress("RunNumber", &Injected_RunNumber);
    Tree["Injected"]->SetBranchAddress("DirNumber", &Injected_DirNumber);
    Tree["Injected"]->SetBranchAddress("EventNumber", &Injected_EventNumber);
    Tree["Injected"]->SetBranchAddress("ReactionID", &Injected_ReactionID);
    Tree["Injected"]->SetBranchAddress("Px", &Injected_Px);
    Tree["Injected"]->SetBranchAddress("Py", &Injected_Py);
    Tree["Injected"]->SetBranchAddress("Pz", &Injected_Pz);
    Tree["Injected"]->SetBranchAddress("M", &Injected_M);

    /* MC Particles */

    // added custom index for easy access when looping over injected anti-sexaquarks
    // "Status" is preferred over "ReactionID" to identify the first-gen products of the interaction,
    // because it contains the secondary vertex info
    TString MCParticles_MajorIndex = "RunNumber";
    TString MCParticles_MinorIndex = "DirNumber * 100 * 1000 + EventNumber * 1000 + Status";
    TTreeIndex* MCParticles_Index = new TTreeIndex(Tree["MCParticles"], MCParticles_MajorIndex, MCParticles_MinorIndex);
    Tree["MCParticles"]->SetTreeIndex(MCParticles_Index);

    Float_t MC_Xv_i, MC_Yv_i;

    Tree["MCParticles"]->SetBranchAddress("Xv_i", &MC_Xv_i);
    Tree["MCParticles"]->SetBranchAddress("Yv_i", &MC_Yv_i);

    /* Events */

    Float_t Event_Centrality;

    Tree["Events"]->SetBranchAddress("Centrality", &Event_Centrality);

    /* Sexaquarks */

    Int_t Sexaquark_RunNumber, Sexaquark_DirNumber, Sexaquark_EventNumber, Sexaquark_ReactionID;
    Float_t Sexaquark_Xv_SV, Sexaquark_Yv_SV;
    Float_t Sexaquark_Px, Sexaquark_Py;

    Tree["Sexaquarks"]->SetBranchAddress("RunNumber", &Sexaquark_RunNumber);
    Tree["Sexaquarks"]->SetBranchAddress("DirNumber", &Sexaquark_DirNumber);
    Tree["Sexaquarks"]->SetBranchAddress("EventNumber", &Sexaquark_EventNumber);
    Tree["Sexaquarks"]->SetBranchAddress("ReactionID", &Sexaquark_ReactionID);
    Tree["Sexaquarks"]->SetBranchAddress("Px", &Sexaquark_Px);
    Tree["Sexaquarks"]->SetBranchAddress("Py", &Sexaquark_Py);
    Tree["Sexaquarks"]->SetBranchAddress("Xv_SV", &Sexaquark_Xv_SV);
    Tree["Sexaquarks"]->SetBranchAddress("Yv_SV", &Sexaquark_Yv_SV);

    /** Auxiliar Variables **/

    Float_t Pt;
    Float_t TruePt;
    Float_t Radius;
    Float_t TrueRadius;

    Int_t RunNumber;
    Int_t DirNumber;
    Int_t EventNumber;
    Int_t ReactionID;
    Int_t mcIdx;

    Float_t SexaquarkMass = 0.;

    /** Part 1a: True **/

    std::map<TString, TH1F*> TrueDistr;
    TrueDistr["Pt"] = new TH1F("True_Pt", "True_Pt", 100, 0., 5.);
    TrueDistr["Radius"] = new TH1F("True_Radius", "True_Radius", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Injected"]->GetEntries(); i++) {
        Tree["Injected"]->GetEntry(i);
        /*  */
        if (SexaquarkMass < 1E-4) SexaquarkMass = Injected_M;
        /*  */
        TruePt = TMath::Sqrt(Injected_Px * Injected_Px + Injected_Py * Injected_Py);
        /*  */
        RunNumber = Injected_RunNumber;
        DirNumber = Injected_DirNumber;
        EventNumber = Injected_EventNumber;
        ReactionID = Injected_ReactionID;
        mcIdx = Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber, DirNumber * 100 * 1000 + EventNumber * 1000 + ReactionID);
        if (mcIdx >= 0) {
            Tree["MCParticles"]->GetEntry(mcIdx);
            TrueRadius = TMath::Sqrt(MC_Xv_i * MC_Xv_i + MC_Yv_i * MC_Yv_i);
        }
        /*  */
        TrueDistr["Pt"]->Fill(TruePt);
        TrueDistr["Radius"]->Fill(TrueRadius);
    }

    /** Part 1b: Reconstructed **/

    std::map<TString, TH1F*> RecDistr;
    RecDistr["Pt"] = new TH1F("Rec_Pt", "Rec_Pt", 100, 0., 5.);
    RecDistr["Radius"] = new TH1F("Rec_Radius", "Rec_Radius", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Sexaquarks"]->GetEntries(); i++) {
        Tree["Sexaquarks"]->GetEntry(i);
        /*  */
        Pt = TMath::Sqrt(Sexaquark_Px * Sexaquark_Px + Sexaquark_Py * Sexaquark_Py);
        /*  */
        Radius = TMath::Sqrt(Sexaquark_Xv_SV * Sexaquark_Xv_SV + Sexaquark_Yv_SV * Sexaquark_Yv_SV);
        /*  */
        RecDistr["Pt"]->Fill(Pt);
        RecDistr["Radius"]->Fill(Radius);
    }

    /** Part 1c: Efficiency **/

    TH1F* EfficiencyPt = new TH1F("EfficiencyPt", "EfficiencyPt", 100, 0., 5.);
    TH1F* EfficiencyRadius = new TH1F("EfficiencyRadius", "EfficiencyRadius", 100, 0., 200.);

    EfficiencyPt->Divide(RecDistr["Pt"], TrueDistr["Pt"], 1., 1., "B");
    EfficiencyRadius->Divide(RecDistr["Radius"], TrueDistr["Radius"], 1., 1., "B");

    /** Part 2a: Get Radius Weights **/

    TFile* PhotonConvRadiusFile = TFile::Open(PhotonConvRadiusPath, "READ");
    if (!PhotonConvRadiusFile || PhotonConvRadiusFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << PhotonConvRadiusPath << " !!" << std::endl;
        return;
    }

    TString PhotonConvHistName = "MCGen_PhotonConversions_Radius";
    TH1F* PhotonConvHist = dynamic_cast<TH1F*>(PhotonConvRadiusFile->Get("Hists")->FindObject(PhotonConvHistName));
    if (!PhotonConvHist) {
        std::cout << "!! ERROR !! Couldn't find histogram " << PhotonConvHistName << " in " << PhotonConvRadiusPath << " !!" << std::endl;
        return;
    }

    // clone correctly the histogram, so I can close the file without any problem
    TH1F* RadiusWeights = dynamic_cast<TH1F*>(PhotonConvHist->Clone());
    RadiusWeights->SetDirectory(0);  // detach from current directory
    RadiusWeights->Scale(1. / RadiusWeights->Integral());

    PhotonConvRadiusFile->Close();

    /** Part 2b: Get (Pt, Centrality) Weights **/

    TFile* BlastWaveFile = TFile::Open(BlastWavePath, "READ");
    std::vector<TString> CentralityBins = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"};
    std::vector<TH1D*> PtWeights;

    for (Int_t cc = 0; cc < (Int_t)CentralityBins.size(); cc++) {
        TString BlastWaveHistName = Form("BlastWaveHist_%.2f_%s", SexaquarkMass, CentralityBins[cc].Data());
        TH1D* BlastWaveHist = dynamic_cast<TH1D*>(BlastWaveFile->Get(BlastWaveHistName));
        /*  */
        TH1D* AuxPtWeights = dynamic_cast<TH1D*>(BlastWaveHist->Clone());
        AuxPtWeights->SetDirectory(0);  // detach from current directory
        AuxPtWeights->Scale(1. / AuxPtWeights->Integral());
        PtWeights.push_back(AuxPtWeights);
    }

    BlastWaveFile->Close();

    /** Part 3: Reweighted Efficiencies **/

    Int_t eventIdx;
    Int_t centralityBin;
    Int_t mcInjected;

    Float_t weightRadius;
    Float_t weightPt;

    // True //

    std::map<TString, TH1F*> TrueDistr_wRadiusWeights;
    TrueDistr_wRadiusWeights["Pt"] = new TH1F("True_Pt_wRadiusWeights", "True_Pt_wRadiusWeights", 100, 0., 5.);
    TrueDistr_wRadiusWeights["Radius"] = new TH1F("True_Radius_wRadiusWeights", "True_Radius_wRadiusWeights", 100, 0., 200.);

    std::map<TString, TH1F*> TrueDistr_wPtWeights;
    TrueDistr_wPtWeights["Pt"] = new TH1F("True_Pt_wPtWeights", "True_Pt_wPtWeights", 100, 0., 5.);
    TrueDistr_wPtWeights["Radius"] = new TH1F("True_Radius_wPtWeights", "True_Radius_wPtWeights", 100, 0., 200.);

    std::map<TString, TH1F*> TrueDistr_wBothWeights;
    TrueDistr_wBothWeights["Pt"] = new TH1F("True_Pt_wBothWeights", "True_Pt_wBothWeights", 100, 0., 5.);
    TrueDistr_wBothWeights["Radius"] = new TH1F("True_Radius_wBothWeights", "True_Radius_wBothWeights", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Injected"]->GetEntries(); i++) {
        Tree["Injected"]->GetEntry(i);
        /*  */
        TruePt = TMath::Sqrt(Injected_Px * Injected_Px + Injected_Py * Injected_Py);
        /*  */
        RunNumber = Injected_RunNumber;
        DirNumber = Injected_DirNumber;
        EventNumber = Injected_EventNumber;
        ReactionID = Injected_ReactionID;
        /*  */
        eventIdx = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);
        if (eventIdx >= 0) {
            Tree["Events"]->GetEntry(eventIdx);
            centralityBin = GetCentralityBin(Event_Centrality);
            weightPt = PtWeights[centralityBin]->GetBinContent(PtWeights[centralityBin]->FindBin(TruePt));
        }
        /*  */
        mcIdx = Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber, DirNumber * 100 * 1000 + EventNumber * 1000 + ReactionID);
        if (mcIdx >= 0) {
            Tree["MCParticles"]->GetEntry(mcIdx);
            TrueRadius = TMath::Sqrt(MC_Xv_i * MC_Xv_i + MC_Yv_i * MC_Yv_i);
            weightRadius = RadiusWeights->GetBinContent(RadiusWeights->FindBin(TrueRadius));
        }
        /*  */
        TrueDistr_wRadiusWeights["Pt"]->Fill(TruePt, weightRadius);
        TrueDistr_wRadiusWeights["Radius"]->Fill(TrueRadius, weightRadius);
        TrueDistr_wPtWeights["Pt"]->Fill(TruePt, weightPt);
        TrueDistr_wPtWeights["Radius"]->Fill(TrueRadius, weightPt);
        TrueDistr_wBothWeights["Pt"]->Fill(TruePt, weightRadius * weightPt);
        TrueDistr_wBothWeights["Radius"]->Fill(TrueRadius, weightRadius * weightPt);
    }

    // Reconstructed //

    std::map<TString, TH1F*> RecDistr_wRadiusWeights;
    RecDistr_wRadiusWeights["Pt"] = new TH1F("Rec_Pt_wRadiusWeights", "Rec_Pt_wRadiusWeights", 100, 0., 5.);
    RecDistr_wRadiusWeights["Radius"] = new TH1F("Rec_Radius_wRadiusWeights", "Rec_Radius_wRadiusWeights", 100, 0., 200.);

    std::map<TString, TH1F*> RecDistr_wPtWeights;
    RecDistr_wPtWeights["Pt"] = new TH1F("Rec_Pt_wPtWeights", "Rec_Pt_wPtWeights", 100, 0., 5.);
    RecDistr_wPtWeights["Radius"] = new TH1F("Rec_Radius_wPtWeights", "Rec_Radius_wPtWeights", 100, 0., 200.);

    std::map<TString, TH1F*> RecDistr_wBothWeights;
    RecDistr_wBothWeights["Pt"] = new TH1F("Rec_Pt_wBothWeights", "Rec_Pt_wBothWeights", 100, 0., 5.);
    RecDistr_wBothWeights["Radius"] = new TH1F("Rec_Radius_wBothWeights", "Rec_Radius_wBothWeights", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Sexaquarks"]->GetEntries(); i++) {
        Tree["Sexaquarks"]->GetEntry(i);
        /*  */
        Pt = TMath::Sqrt(Sexaquark_Px * Sexaquark_Px + Sexaquark_Py * Sexaquark_Py);
        /*  */
        Radius = TMath::Sqrt(Sexaquark_Xv_SV * Sexaquark_Xv_SV + Sexaquark_Yv_SV * Sexaquark_Yv_SV);
        /*  */
        RunNumber = Sexaquark_RunNumber;
        DirNumber = Sexaquark_DirNumber;
        EventNumber = Sexaquark_EventNumber;
        ReactionID = Sexaquark_ReactionID;
        /*  */
        eventIdx = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);
        if (eventIdx >= 0) {
            Tree["Events"]->GetEntry(eventIdx);
            centralityBin = GetCentralityBin(Event_Centrality);
        }
        /*  */
        mcInjected = Tree["Injected"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 1000 + ReactionID);
        if (mcInjected >= 0) {
            Tree["Injected"]->GetEntry(mcInjected);
            TruePt = TMath::Sqrt(Injected_Px * Injected_Px + Injected_Py * Injected_Py);
            /*  */
            mcIdx = Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber, DirNumber * 100 * 1000 + EventNumber * 1000 + ReactionID);
            if (mcIdx >= 0) {
                Tree["MCParticles"]->GetEntry(mcIdx);
                TrueRadius = TMath::Sqrt(MC_Xv_i * MC_Xv_i + MC_Yv_i * MC_Yv_i);
            }
        }
        /*  */
        weightRadius = RadiusWeights->GetBinContent(RadiusWeights->FindBin(TrueRadius));
        weightPt = PtWeights[centralityBin]->GetBinContent(PtWeights[centralityBin]->FindBin(TruePt));
        /*  */
        RecDistr_wRadiusWeights["Pt"]->Fill(Pt, weightRadius);
        RecDistr_wRadiusWeights["Radius"]->Fill(Radius, weightRadius);
        RecDistr_wPtWeights["Pt"]->Fill(Pt, weightPt);
        RecDistr_wPtWeights["Radius"]->Fill(Radius, weightPt);
        RecDistr_wBothWeights["Pt"]->Fill(Pt, weightRadius * weightPt);
        RecDistr_wBothWeights["Radius"]->Fill(Radius, weightRadius * weightPt);
    }

    // Efficiency (with Radius Weights) //

    TH1F* EfficiencyPt_wRadiusWeights = new TH1F("EfficiencyPt_wRadiusWeights", "EfficiencyPt_wRadiusWeights", 100, 0., 5.);
    TH1F* EfficiencyRadius_wRadiusWeights = new TH1F("EfficiencyRadius_wRadiusWeights", "EfficiencyRadius_wRadiusWeights", 100, 0., 200.);

    EfficiencyPt_wRadiusWeights->Divide(RecDistr_wRadiusWeights["Pt"], TrueDistr_wRadiusWeights["Pt"], 1., 1., "B");
    EfficiencyRadius_wRadiusWeights->Divide(RecDistr_wRadiusWeights["Radius"], TrueDistr_wRadiusWeights["Radius"], 1., 1., "B");

    // Efficiency (with Pt Weights) //

    TH1F* EfficiencyPt_wPtWeights = new TH1F("EfficiencyPt_wPtWeights", "EfficiencyPt_wPtWeights", 100, 0., 5.);
    TH1F* EfficiencyRadius_wPtWeights = new TH1F("EfficiencyRadius_wPtWeights", "EfficiencyRadius_wPtWeights", 100, 0., 200.);

    EfficiencyPt_wPtWeights->Divide(RecDistr_wPtWeights["Pt"], TrueDistr_wPtWeights["Pt"], 1., 1., "B");
    EfficiencyRadius_wPtWeights->Divide(RecDistr_wPtWeights["Radius"], TrueDistr_wPtWeights["Radius"], 1., 1., "B");

    // Efficiency (with Both Weights) //

    TH1F* EfficiencyPt_wBothWeights = new TH1F("EfficiencyPt_wBothWeights", "EfficiencyPt_wBothWeights", 100, 0., 5.);
    TH1F* EfficiencyRadius_wBothWeights = new TH1F("EfficiencyRadius_wBothWeights", "EfficiencyRadius_wBothWeights", 100, 0., 200.);

    EfficiencyPt_wBothWeights->Divide(RecDistr_wBothWeights["Pt"], TrueDistr_wBothWeights["Pt"], 1., 1., "B");
    EfficiencyRadius_wBothWeights->Divide(RecDistr_wBothWeights["Radius"], TrueDistr_wBothWeights["Radius"], 1., 1., "B");

    /** Part 4: Cross-Check Average Efficiency **/

    Double_t AverageEfficiencyPt = 0.;
    Double_t AverageEfficiencyPt_wRadiusWeights = 0.;
    Double_t AverageEfficiencyPt_wPtWeights = 0.;
    Double_t AverageEfficiencyPt_wBothWeights = 0.;
    Int_t NonEmptyBinsPt = 0.;
    for (Int_t xbin = 1; xbin <= EfficiencyPt->GetNbinsX(); xbin++) {
        // if (EfficiencyPt->GetBinContent(xbin) > 0.) { // OPTIONAL
        AverageEfficiencyPt += EfficiencyPt->GetBinContent(xbin);
        AverageEfficiencyPt_wRadiusWeights += EfficiencyPt_wRadiusWeights->GetBinContent(xbin);
        AverageEfficiencyPt_wPtWeights += EfficiencyPt_wPtWeights->GetBinContent(xbin);
        AverageEfficiencyPt_wBothWeights += EfficiencyPt_wBothWeights->GetBinContent(xbin);
        NonEmptyBinsPt++;
        // }
    }
    AverageEfficiencyPt /= (Double_t)NonEmptyBinsPt;
    AverageEfficiencyPt_wRadiusWeights /= (Double_t)NonEmptyBinsPt;
    AverageEfficiencyPt_wPtWeights /= (Double_t)NonEmptyBinsPt;
    AverageEfficiencyPt_wBothWeights /= (Double_t)NonEmptyBinsPt;

    Double_t AverageEfficiencyRadius = 0.;
    Double_t AverageEfficiencyRadius_wRadiusWeights = 0.;
    Double_t AverageEfficiencyRadius_wPtWeights = 0.;
    Double_t AverageEfficiencyRadius_wBothWeights = 0.;
    Int_t NonEmptyBinsRadius = 0.;
    for (Int_t xbin = 1; xbin <= EfficiencyRadius->GetNbinsX(); xbin++) {
        // if (EfficiencyRadius->GetBinContent(xbin) > 0.) { // OPTIONAL
        AverageEfficiencyRadius += EfficiencyRadius->GetBinContent(xbin);
        AverageEfficiencyRadius_wRadiusWeights += EfficiencyRadius_wRadiusWeights->GetBinContent(xbin);
        AverageEfficiencyRadius_wPtWeights += EfficiencyRadius_wPtWeights->GetBinContent(xbin);
        AverageEfficiencyRadius_wBothWeights += EfficiencyRadius_wBothWeights->GetBinContent(xbin);
        NonEmptyBinsRadius++;
        // }
    }
    AverageEfficiencyRadius /= (Double_t)NonEmptyBinsRadius;
    AverageEfficiencyRadius_wRadiusWeights /= (Double_t)NonEmptyBinsRadius;
    AverageEfficiencyRadius_wPtWeights /= (Double_t)NonEmptyBinsRadius;
    AverageEfficiencyRadius_wBothWeights /= (Double_t)NonEmptyBinsRadius;

    std::cout << "Average Efficiency (Pt)                         : " << AverageEfficiencyPt << std::endl;
    std::cout << "Average Efficiency (Radius)                     : " << AverageEfficiencyRadius << std::endl;

    std::cout << "Average Efficiency (Pt) with Radius Weights     : " << AverageEfficiencyPt_wRadiusWeights << std::endl;
    std::cout << "Average Efficiency (Radius) with Radius Weights : " << AverageEfficiencyRadius_wRadiusWeights << std::endl;

    std::cout << "Average Efficiency (Pt) with Pt Weights         : " << AverageEfficiencyPt_wPtWeights << std::endl;
    std::cout << "Average Efficiency (Radius) with Pt Weights     : " << AverageEfficiencyRadius_wPtWeights << std::endl;

    std::cout << "Average Efficiency (Pt) with Both Weights       : " << AverageEfficiencyPt_wBothWeights << std::endl;
    std::cout << "Average Efficiency (Radius) with Both Weights   : " << AverageEfficiencyRadius_wBothWeights << std::endl;

    /** Output File **/

    TString OutputFileName = InputFileName.ReplaceAll("Analysis", "Efficiency").ReplaceAll("_indexed", "");
    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");

    TrueDistr["Pt"]->Write();
    RecDistr["Pt"]->Write();
    EfficiencyPt->Write();

    TrueDistr["Radius"]->Write();
    RecDistr["Radius"]->Write();
    EfficiencyRadius->Write();

    RadiusWeights->Write();
    for (const auto& SinglePtWeights : PtWeights) SinglePtWeights->Write();

    RecDistr_wRadiusWeights["Pt"]->Write();
    TrueDistr_wRadiusWeights["Pt"]->Write();
    EfficiencyPt_wRadiusWeights->Write();

    RecDistr_wRadiusWeights["Radius"]->Write();
    TrueDistr_wRadiusWeights["Radius"]->Write();
    EfficiencyRadius_wRadiusWeights->Write();

    RecDistr_wPtWeights["Pt"]->Write();
    TrueDistr_wPtWeights["Pt"]->Write();
    EfficiencyPt_wPtWeights->Write();

    RecDistr_wPtWeights["Radius"]->Write();
    TrueDistr_wPtWeights["Radius"]->Write();
    EfficiencyRadius_wPtWeights->Write();

    RecDistr_wBothWeights["Pt"]->Write();
    TrueDistr_wBothWeights["Pt"]->Write();
    EfficiencyPt_wBothWeights->Write();

    RecDistr_wBothWeights["Radius"]->Write();
    TrueDistr_wBothWeights["Radius"]->Write();
    EfficiencyRadius_wBothWeights->Write();

    OutputFile->Close();
    InputFile->Close();
}

/*
 * Get the centrality bin from the centrality value. According to:
 * - centrality ranges: 0-5, 5-10, 10-20, 20-30, 30-40, 40-50, 50-60, 60-70, 70-80, 80-90
 * - centrality bin: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
 */
Int_t GetCentralityBin(Float_t CentralityValue) {
    Int_t centrality_bin = std::min(static_cast<Int_t>(0.1 * CentralityValue + 1), 9);
    if (CentralityValue < 5.) centrality_bin = 0;
    return centrality_bin;
}
