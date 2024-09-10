#include "include/Headers.hxx"

void Trees_BuildIndices(TString InputFileName = "AnalysisResults.root") {

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
    TString MinorIndex_Event = "Number";

    /* Loop over trees */

    TTree* InputEventsTree;
    TTree* NewEventsTree;

    std::vector<TString> TreeNames = {"Injected",    "Events",      "MCParticles",    "PiPluses",  "PiMinuses",
                                      "AntiProtons", "AntiLambdas", "KaonsZeroShort", "Sexaquarks"};

    for (Int_t i = 0; i < (Int_t)TreeNames.size(); i++) {
        InputEventsTree = (TTree*)TreesList->FindObject(TreeNames[i]);
        if (!InputEventsTree) {
            std::cout << "!! ERROR !! Couldn't find TTree " << TreeNames[i] << " in " << TreesListName << " !!" << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        std::cout << "!! INFO !! Input TTree " << TreeNames[i] << " has " << InputEventsTree->GetEntries() << " entries !!" << std::endl;

        OutputFile->cd();

        NewEventsTree = (TTree*)InputEventsTree->CloneTree();
        if (!NewEventsTree) {
            std::cout << "!! ERROR !! Failed to clone TTree " << TreeNames[i] << " !!" << std::endl;
            InputFile->Close();
            OutputFile->Close();
            return;
        }
        std::cout << "!! INFO !! Cloned TTree " << TreeNames[i] << " has " << NewEventsTree->GetEntries() << " entries !!" << std::endl;

        MinorIndex = MinorIndex_General;
        if (TreeNames[i] == "Sexaquarks" || TreeNames[i] == "Injected") MinorIndex = MinorIndex_Sexaquark;
        if (TreeNames[i] == "Events") MinorIndex = MinorIndex_Event;

        NewEventsTree->BuildIndex(MajorIndex, MinorIndex);
        std::cout << "!! INFO !! Indexing complete for TTree " << NewEventsTree->GetName() << " !!" << std::endl;

        if (TreeNames[i] == "MCParticles") {
            // Note:
            // added custom index for easy access when looping over injected anti-sexaquarks
            // "Status" preferred over "ReactionID" to identify the first-gen products of the interaction,
            // which are the ones that contain the secondary vertex info
            TString CustomMajorIndex = "RunNumber";
            TString CustomMinorIndex = "DirNumber * 100 * 1000 + EventNumber * 1000 + Status";
            NewEventsTree->BuildIndex(CustomMajorIndex, CustomMinorIndex);
            std::cout << "!! INFO !! Custom index built for TTree " << NewEventsTree->GetName() << " !!" << std::endl;
        }

        NewEventsTree->Write();
        std::cout << "!! INFO !! Indexed TTree " << NewEventsTree->GetName() << " stored in " << OutputFileName << " !!" << std::endl;
    }

    OutputFile->Close();
    InputFile->Close();

    return;
}
