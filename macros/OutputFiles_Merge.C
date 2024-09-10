#include "include/Headers.hxx"

#include "../sexaquark/AliAnalysisTaskSexaquark_Structs.h"

void OutputFiles_Merge(TString InputDir, TString OutputFilename = "AnalysisResults_merged2.root") {

    TSystemDirectory dir(InputDir, InputDir);
    TList *files = dir.GetListOfFiles();

    std::cout << "OutputFiles_Merge :: Merging " << files->GetEntries() << " files into " << OutputFilename << "..." << std::endl;

    /* Define General Properties */

    TIter next(files);
    TString StrList_Trees = "Trees";
    std::vector<TString> Str_Chains = {"Injected",    "Events",      "MCParticles",    "PiPluses",  "PiMinuses",
                                       "AntiProtons", "AntiLambdas", "KaonsZeroShort", "Sexaquarks"};
    std::vector<TString> StrList_Hists = {"QA_Hists", "Tracks_Hists", "V0s_Hists", "Sexaquarks_Hists", "PosKaonPairs_Hists"};

    /* Create Output File */

    TFile *OutputFile = new TFile(OutputFilename, "RECREATE");

    TList *OutputList_Trees = new TList();
    std::vector<TChain *> Output_Chains;
    for (Int_t tt = 0; tt < (Int_t)Str_Chains.size(); tt++) Output_Chains.push_back(new TChain(Str_Chains[tt]));  // TO CONTINUE

    std::vector<TList *> OutputList_Hists;
    for (Int_t hh = 0; hh < (Int_t)StrList_Hists.size(); hh++) OutputList_Hists.push_back(new TList());

    /* Define Properties of Current File */

    TSystemFile *file;
    TString filename, filepath;
    TFile *InputFile;

    TList *InputList_Trees;
    TList *InputList_QAHists;

    // loop over files
    while ((file = (TSystemFile *)next())) {

        // get index of file within the TList
        Int_t index = files->IndexOf(file);
        if (index > 2) break;
        std::cout << "OutputFiles_Merge :: Processing file " << index << " of " << files->GetEntries() << "..." << std::endl;

        filename = file->GetName();
        filepath = InputDir + "/" + filename;

        if (!filename.Contains("AnalysisResults_")) continue;

        InputFile = new TFile(filepath, "READ");

        /* Trees */

        InputList_Trees = (TList *)InputFile->Get(StrList_Trees);

        // loop over trees
        TIter nextTree(InputList_Trees);
        TTree *auxTree;
        while ((auxTree = (TTree *)nextTree())) {
            TString treeName = auxTree->GetName();
            OutputFile->cd();
            if (!OutputList_Trees->FindObject(treeName)) {
                // if doesn't exist, clone it
                TTree *newTree = auxTree->CloneTree();
                OutputList_Trees->Add(newTree);
            } else {
                // if already exists, merge tree
                TTree *oldTree = (TTree *)OutputList_Trees->FindObject(treeName);
                TList *auxList = new TList();
                auxList->Add(auxTree);
                oldTree->Merge(auxList);
                delete auxList;
            }
        }

        /* Histograms */

        // loop over histograms lists
        TList *OutputList;
        for (Int_t hh = 0; hh < (Int_t)StrList_Hists.size(); hh++) {

            OutputList = OutputList_Hists[hh];
            InputList_QAHists = (TList *)InputFile->Get(StrList_Hists[hh]);

            TIter nextHist(InputList_QAHists);
            TObject *auxObj;

            // loop over hists
            while ((auxObj = (TObject *)nextHist())) {
                OutputFile->cd();
                if (!OutputList->FindObject(auxObj->GetName())) {
                    // if doesn't exist, clone it
                    TObject *newHist = auxObj->Clone();
                    OutputList->Add(newHist);
                } else {
                    // if already exists, add hist
                    TH1 *oldHist = (TH1 *)OutputList->FindObject(auxObj->GetName());
                    TH1 *newHist = (TH1 *)auxObj;
                    oldHist->Add(newHist);
                }
            }
        }

        InputFile->Close();
    }

    OutputFile->cd();
    OutputList_Trees->Write("Trees", TObject::kSingleKey);
    for (Int_t hh = 0; hh < (Int_t)StrList_Hists.size(); hh++) OutputList_Hists[hh]->Write(StrList_Hists[hh], TObject::kSingleKey);

    OutputFile->Close();
}
