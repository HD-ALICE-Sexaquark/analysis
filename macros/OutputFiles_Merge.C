#include "include/Headers.hxx"

void OutputFiles_Merge(TString InputDir, TString OutputFilename = "AnalysisResults_merged.root") {

    TSystemDirectory dir(InputDir, InputDir);
    TList *files = dir.GetListOfFiles();

    /* Define General Properties */

    std::vector<TString> StrInputFiles;

    TIter next(files);
    TString StrList_Trees = "Trees";
    // TString StrList_QAHists = "QAHists";
    // TString StrList_TrackHists = "Track_Hists";
    // TString StrList_V0sHists = "V0s_Hists";
    // TString StrList_SexaquarksHists = "Sexaquarks_Hists";
    // TString StrList_PosKaonPairsHists = "PosKaonPairs_Hists";

    /* Create Output File */

    TFile *OutputFile = new TFile(OutputFilename, "RECREATE");

    TList *OutputList_Trees = new TList();
    // TList *List_QAHists = new TList();
    // TList *List_TrackHists = new TList();
    // TList *List_V0sHists = new TList();
    // TList *List_SexaquarksHists = new TList();
    // TList *List_PosKaonPairsHists = new TList();

    /* Define Properties of Current File */

    TSystemFile *file;
    TString filename, filepath;
    TFile *InputFile;

    TList *InputList_Trees;

    Int_t FileCounter = 0;
    while ((file = (TSystemFile *)next())) {
        if (FileCounter > 1) break;

        filename = file->GetName();
        filepath = InputDir + "/" + filename;

        if (!filename.Contains("AnalysisResults_")) continue;

        std::cout << "file: " << filename << std::endl;
        StrInputFiles.push_back(filename);

        InputFile = new TFile(filepath, "READ");

        /* Trees */

        InputList_Trees = (TList *)InputFile->Get(StrList_Trees);

        TIter nextTree(InputList_Trees);
        TTree *obj;

        while ((obj = (TTree *)nextTree())) {
            OutputFile->cd();
            if (FileCounter == 0) {
                // add it
                TTree *newTree = obj->CloneTree();
                OutputList_Trees->Add(newTree);
            } else {
                // merge tree of the same name
                TTree *oldTree = (TTree *)OutputList_Trees->FindObject(obj->GetName());
                TList *auxList = new TList();
                auxList->Add(obj);
                oldTree->Merge(auxList);
                delete auxList;
            }
        }

        FileCounter++;
        InputFile->Close();
    }

    std::cout << "Merging " << StrInputFiles.size() << " files!" << std::endl;

    OutputFile->cd();
    OutputList_Trees->Write("Trees", TObject::kSingleKey);

    OutputFile->Close();
}
