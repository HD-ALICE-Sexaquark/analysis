#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "../sexaquark/AliAnalysisTaskSexaquark_Structs.h"

// clang-format off
// Indices
//     MajorIndex
//         Injected:    "Injected.RunNumber * 1000 + Injected.DirNumber"
//         Events:      "Event.RunNumber * 1000 + Event.DirNumber"
//         MCParticles: "MCParticle.RunNumber * 1000 + MCParticle.DirNumber"
//         Tracks:      "Track.RunNumber * 1000 + Track.DirNumber"
//         V0s:         "V0.RunNumber * 1000 + V0.DirNumber"
//         Sexaquarks:  "Sexaquark.RunNumber * 1000 + Sexaquark.DirNumber"
//     MinorIndex
//         Injected:    "Injected.EventNumber * 1000 + Injected.ReactionID"
//         Events:      "Event.Number"
//         MCParticles: "MCParticle.EventNumber * 100000 + MCParticle.Idx" OR "MCParticle.EventNumber * 1000 + MCParticle.ReactionID"
//         Tracks:      "Track.EventNumber * 100000 + Track.Idx"
//         V0s:         "V0.EventNumber * 100000 + V0.Idx"
//         Sexaquarks:  "Sexaquark.EventNumber * 1000 + Sexaquark.ReactionID"
// clang-format on

void Trees_GetEfficiency(TString InputFileName = "AnalysisResults_indexed.root") {

    /** Input **/

    TFile* InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << InputFileName << " !!" << std::endl;
        return;
    }

    std::map<TString, TTree*> Tree;
    std::vector<TString> TreeNames = {"Injected",    "Events",      "MCParticles",    "PiPluses",  "PiMinuses",
                                      "AntiProtons", "AntiLambdas", "KaonsZeroShort", "Sexaquarks"};

    for (const auto& treeName : TreeNames) {
        Tree[treeName] = (TTree*)InputFile->Get(treeName);
        if (!Tree[treeName]) {
            std::cout << "!! ERROR !! Couldn't find TTree " << treeName << " in " << InputFileName << " !!" << std::endl;
            InputFile->Close();
            return;
        }
    }

    /** Set Branches **/

    Injected_tt* thisInjected = nullptr;
    Tree["Injected"]->SetBranchAddress("Injected", &thisInjected);

    Event_tt* thisEvent = nullptr;
    Tree["Events"]->SetBranchAddress("Event", &thisEvent);

    MC_tt* thisMC = nullptr;
    Tree["MCParticles"]->SetBranchAddress("MCParticle", &thisMC);

    V0_tt* thisAL = nullptr;
    Tree["AntiLambdas"]->SetBranchAddress("AntiLambda", &thisAL);

    V0_tt* thisK0S = nullptr;
    Tree["KaonsZeroShort"]->SetBranchAddress("KaonZeroShort", &thisK0S);

    RecSexaquark_aa* thisSexaquark = nullptr;
    Tree["Sexaquarks"]->SetBranchAddress("Sexaquark", &thisSexaquark);

    /** Auxiliar Variables **/

    Float_t Pt;
    Float_t Radius;

    Float_t RunNumber;
    Float_t DirNumber;
    Float_t EventNumber;
    Float_t ReactionID;
    Int_t mcIdx;

    /** Part 1: Reconstructed **/

    std::map<TString, TH1F*> RecDistr;
    RecDistr["Pt"] = new TH1F("Rec_Pt", "Rec_Pt", 100, 0., 5.);
    RecDistr["Radius"] = new TH1F("Rec_Radius", "Rec_Radius", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Sexaquarks"]->GetEntries(); i++) {
        Tree["Sexaquarks"]->GetEntry(i);
        /*  */
        Pt = TMath::Sqrt(thisSexaquark->Px * thisSexaquark->Px + thisSexaquark->Py * thisSexaquark->Py);
        RecDistr["Pt"]->Fill(Pt);
        /*  */
        Radius = TMath::Sqrt(thisSexaquark->Xv_SV * thisSexaquark->Xv_SV + thisSexaquark->Yv_SV * thisSexaquark->Yv_SV);
        RecDistr["Radius"]->Fill(Radius);
    }

    /** Part 2: True **/

    std::map<TString, TH1F*> TrueDistr;
    TrueDistr["Pt"] = new TH1F("True_Pt", "True_Pt", 100, 0., 5.);
    TrueDistr["Radius"] = new TH1F("True_Radius", "True_Radius", 100, 0., 200.);

    for (Int_t i = 0; i < Tree["Injected"]->GetEntries(); i++) {
        Tree["Injected"]->GetEntry(i);
        /*  */
        Pt = TMath::Sqrt(thisInjected->Px * thisInjected->Px + thisInjected->Py * thisInjected->Py);
        TrueDistr["Pt"]->Fill(Pt);
        /*  */
        RunNumber = thisInjected->RunNumber;
        DirNumber = thisInjected->DirNumber;
        EventNumber = thisInjected->EventNumber;
        ReactionID = thisInjected->ReactionID;
        mcIdx = Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 1000 + ReactionID);
        if (mcIdx < 0) continue;
        Tree["MCParticles"]->GetEntry(mcIdx);
        Radius = TMath::Sqrt(thisMC->Xv_i * thisMC->Xv_i + thisMC->Yv_i * thisMC->Yv_i);
        std::cout << thisMC->PdgCode << " " << Radius << std::endl;
        TrueDistr["Radius"]->Fill(Radius);
    }

    /** Output File **/

    TString OutputFileName = InputFileName.ReplaceAll("Analysis", "Efficiency").ReplaceAll("_indexed", "");
    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");
    for (const auto& histName : {"Pt", "Radius"}) {
        RecDistr[histName]->Write();
        TrueDistr[histName]->Write();
    }
    OutputFile->Close();

    InputFile->Close();
}
