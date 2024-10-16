#include "include/Headers.hxx"

/*
 * Process an `AnalysisResults.root` file. Create indexed versions of each tree to improve data access efficiency.
 * Build primary and secondary indices for each tree.
 * Save the indexed trees to a new output file with `_indexed` appended to the original filename.
 */
void Trees_BuildIndices(TString InputFileName = "../output/A1.8_18qr_local/AnalysisResults_merged.root") {

    TFile* InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Trees_BuildIndices !! Couldn't open TFile " << InputFileName << " !!" << std::endl;
        return;
    }

    /* Open Output File */

    TString OutputFileName = InputFileName;
    OutputFileName.ReplaceAll("_merged.root", "_indexed.root");
    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");
    if (!OutputFile) {
        std::cout << "!! ERROR !! Trees_BuildIndices !! Couldn't create TFile " << OutputFileName << " !!" << std::endl;
        InputFile->Close();
        return;
    }

    /* Define Major and Minor Indices */

    TString MajorIndex = "RunNumber * 1000 + DirNumber";

    TString MinorIndex;
    TString MinorIndex_General = "EventNumber * 10000000 + Idx";
    TString MinorIndex_Sexaquark = "EventNumber * 1000 + ReactionID";
    TString MinorIndex_Event = "EventNumber";

    /* Loop over trees */

    TTree* InputEventsTree;
    TTree* NewEventsTree;

    std::vector<TString> TreeNames = {"Injected", "Events", "MCParticles", "Tracks", "V0s", "Sexaquarks"};
    // std::vector<TString> TreeNames = {"Tracks", "V0s", "Sexaquarks"};
    // std::vector<TString> TreeNames = {"MCParticles"};

    for (Int_t i = 0; i < (Int_t)TreeNames.size(); i++) {
        InputEventsTree = InputFile->Get<TTree>(TreeNames[i]);
        if (!InputEventsTree) {
            std::cout << "!! ERROR !! Trees_BuildIndices !! Couldn't find TTree " << TreeNames[i] << " in TFile " << InputFileName << " !!"
                      << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        std::cout << "!! INFO  !! Trees_BuildIndices !! Input TTree " << TreeNames[i] << " has " << InputEventsTree->GetEntries() << " entries !!"
                  << std::endl;  // DEBUG

        NewEventsTree = InputEventsTree->CloneTree();
        if (!NewEventsTree) {
            std::cout << "!! ERROR !! Trees_BuildIndices !! Failed to clone TTree " << TreeNames[i] << " !!" << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        NewEventsTree->SetDirectory(0);

        std::cout << "!! INFO  !! Trees_BuildIndices !! Cloned TTree " << TreeNames[i] << " has " << NewEventsTree->GetEntries() << " entries !!"
                  << std::endl;  // DEBUG

        MinorIndex = MinorIndex_General;
        if (TreeNames[i] == "Sexaquarks" || TreeNames[i] == "Injected") MinorIndex = MinorIndex_Sexaquark;
        if (TreeNames[i] == "Events") MinorIndex = MinorIndex_Event;

        NewEventsTree->BuildIndex(MajorIndex, MinorIndex);
        std::cout << "!! INFO  !! Trees_BuildIndices !! Indexing complete for TTree " << NewEventsTree->GetName() << " !!" << std::endl;

        NewEventsTree->Write();
        std::cout << "!! INFO  !! Trees_BuildIndices !! Indexed TTree " << NewEventsTree->GetName() << " stored in TFile " << OutputFileName << " !!"
                  << std::endl;  // DEBUG
    }

    OutputFile->Close();
    InputFile->Close();

    std::cout << "!! INFO  !! Trees_BuildIndices !! All new trees have been stored in TFile " << OutputFileName << " !!" << std::endl;

    return;
}
