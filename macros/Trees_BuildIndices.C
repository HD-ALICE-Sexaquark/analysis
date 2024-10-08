#include "include/Headers.hxx"

/*
 * Process an `AnalysisResults.root` file. Create indexed versions of each tree to improve data access efficiency.
 * Build primary and secondary indices for each tree.
 * Save the indexed trees to a new output file with `_indexed` appended to the original filename.
 */
void Trees_BuildIndices(TString InputFileName = "AnalysisResults.root", TString ReactionLetter = "A") {

    TFile* InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << InputFileName << " !!" << std::endl;
        return;
    }

    /* Get Trees List */

    TString TreesListName = "Trees";
    TList* TreesList = (TList*)InputFile->Get(TreesListName);
    if (!TreesList) {
        std::cout << "!! ERROR !! Couldn't find TList " << TreesListName << " in " << InputFileName << " !!" << std::endl;
        InputFile->Close();
        return;
    }

    /* Open Output File */

    TString OutputFileName = InputFileName.ReplaceAll(".root", "_indexed.root");
    TFile* OutputFile = TFile::Open(OutputFileName, "RECREATE");
    if (!OutputFile) {
        std::cout << "!! ERROR !! Couldn't create file " << OutputFileName << " !!" << std::endl;
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

    std::vector<TString> TreeNames = {"Injected", "Events", "MCParticles", "Sexaquarks", "AntiLambdas", "AntiProtons", "PiPluses"};
    if (ReactionLetter == "A") {
        TreeNames.push_back("KaonsZeroShort");
        TreeNames.push_back("PiMinuses");
    } else if (ReactionLetter == "D") {
        TreeNames.push_back("PosKaons");
    } else if (ReactionLetter == "E") {
        TreeNames.push_back("PionPairs");
        TreeNames.push_back("PiMinuses");
        TreeNames.push_back("PosKaons");
    }

    for (Int_t i = 0; i < (Int_t)TreeNames.size(); i++) {
        InputEventsTree = (TTree*)TreesList->FindObject(TreeNames[i]);
        if (!InputEventsTree) {
            std::cout << "!! ERROR !! Couldn't find TTree " << TreeNames[i] << " in " << TreesListName << " !!" << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        // std::cout << "!! INFO !! Input TTree " << TreeNames[i] << " has " << InputEventsTree->GetEntries() << " entries !!" << std::endl; // DEBUG

        OutputFile->cd();

        NewEventsTree = (TTree*)InputEventsTree->CloneTree();
        if (!NewEventsTree) {
            std::cout << "!! ERROR !! Failed to clone TTree " << TreeNames[i] << " !!" << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        // std::cout << "!! INFO !! Cloned TTree " << TreeNames[i] << " has " << NewEventsTree->GetEntries() << " entries !!" << std::endl; // DEBUG

        MinorIndex = MinorIndex_General;
        if (TreeNames[i] == "Sexaquarks" || TreeNames[i] == "Injected") MinorIndex = MinorIndex_Sexaquark;
        if (TreeNames[i] == "Events") MinorIndex = MinorIndex_Event;

        NewEventsTree->BuildIndex(MajorIndex, MinorIndex);
        std::cout << "!! INFO !! Indexing complete for TTree " << NewEventsTree->GetName() << " !!" << std::endl;

        NewEventsTree->Write();
        // std::cout << "!! INFO !! Indexed TTree " << NewEventsTree->GetName() << " stored in " << OutputFileName << " !!" << std::endl; // DEBUG
    }

    OutputFile->Close();
    InputFile->Close();

    std::cout << "!! INFO !! All new trees have been stored in " << OutputFileName << " !!" << std::endl;

    return;
}
