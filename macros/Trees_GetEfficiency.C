#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TTreeIndex.h"

#include "../sexaquark/AliAnalysisTaskSexaquark_Structs.h"

// clang-format off
/* 
 * Indices
 *     MajorIndex
 *         Injected:    "Injected.RunNumber * 1000 + Injected.DirNumber"
 *         Events:      "Event.RunNumber * 1000 + Event.DirNumber"
 *         MCParticles: "MCParticle.RunNumber"
 *                      (^^ custom index ^^)
 *         Sexaquarks:  "Sexaquark.RunNumber * 1000 + Sexaquark.DirNumber"
 *     MinorIndex
 *         Injected:    "Injected.EventNumber * 1000 + Injected.ReactionID"
 *         Events:      "Event.Number"
 *         MCParticles: "MCParticle.DirNumber * 100 * 1000 + MCParticle.EventNumber * 1000 + MCParticle.Status"
 *                      (^^ custom index ^^)
 *         Sexaquarks:  "Sexaquark.EventNumber * 1000 + Sexaquark.ReactionID"
 */
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
        mcIdx = Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber, DirNumber * 100 * 1000 + EventNumber * 1000 + ReactionID);
        if (mcIdx >= 0) {
            Tree["MCParticles"]->GetEntry(mcIdx);
            Radius = TMath::Sqrt(thisMC->Xv_i * thisMC->Xv_i + thisMC->Yv_i * thisMC->Yv_i);
            TrueDistr["Radius"]->Fill(Radius);
        }
    }

    /** Part 3: Efficiency **/

    TH1F* Eff_Pt = new TH1F("Eff_Pt", "Eff_Pt", 100, 0., 5.);
    TH1F* Eff_Radius = new TH1F("Eff_Radius", "Eff_Radius", 100, 0., 200.);

    Eff_Pt->Divide(RecDistr["Pt"], TrueDistr["Pt"], 1., 1., "B");
    Eff_Radius->Divide(RecDistr["Radius"], TrueDistr["Radius"], 1., 1., "B");

    /** Output File **/

    TString OutputFileName = InputFileName.ReplaceAll("Analysis", "Efficiency").ReplaceAll("_indexed", "");
    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");

    TrueDistr["Pt"]->Write();
    RecDistr["Pt"]->Write();
    Eff_Pt->Write();

    TrueDistr["Radius"]->Write();
    RecDistr["Radius"]->Write();
    Eff_Radius->Write();

    OutputFile->Close();

    InputFile->Close();
}
