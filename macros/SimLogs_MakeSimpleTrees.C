#include "include/Headers.hxx"

/*
 Read the `sim.log` files, extract initial anti-sexaquark and nucleon kinematics, prior to their reaction, and put it into a simple tree
 */
void SimLogs_MakeSimpleTrees() {

    /* Hard-coded parameters */

    TString pathToSimDir = "/home/ceres/borquez/work/simulations/samples/";
    TSystemDirectory SimDir("SimDir", pathToSimDir);

    /* Define output file */

    TFile *SimpleFile = new TFile("test.root", "RECREATE");

    /* Define tree structure */

    TTree *SimpleTree = new TTree("IC", "IC");

    Int_t eventID;
    Int_t runNumber;
    Char_t reactionChannel;
    std::vector<Int_t> reactionID;
    std::vector<Double_t> SM;
    std::vector<Double_t> SPx, SPy, SPz;
    std::vector<Int_t> NPDGCode;
    std::vector<Double_t> NPx, NPy, NPz;

    SimpleTree->Branch("eventID", &eventID);
    SimpleTree->Branch("runNumber", &runNumber);
    SimpleTree->Branch("reactionChannel", &reactionChannel);
    SimpleTree->Branch("reactionID", &reactionID);
    // SimpleTree->Branch("SM", &SM);
    SimpleTree->Branch("SPx", &SPx);
    SimpleTree->Branch("SPy", &SPy);
    SimpleTree->Branch("SPz", &SPz);
    SimpleTree->Branch("NPDGCode", &NPDGCode);
    SimpleTree->Branch("NPx", &NPx);
    SimpleTree->Branch("NPy", &NPy);
    SimpleTree->Branch("NPz", &NPz);

    /* Loop over directories and files */

    // RCMS = <ReactionChannel><Mass>_<Server>

    TList *listOfRcmsDir = SimDir.GetListOfFiles();
    if (!listOfRcmsDir) return;

    TIter nextRcmsDir(listOfRcmsDir);
    TObject *currentRcsmDir = nullptr;

    // auxiliary variables
    std::string str_line;
    Int_t current_eventID = -1;
    TString csv;
    TObjArray *csv_arr = nullptr;

    while ((currentRcsmDir = nextRcmsDir())) {

        // protection, TString::CompareTo() returns zero when two strings are identical
        if (!static_cast<TString>(currentRcsmDir->GetName()).CompareTo(".") ||
            !static_cast<TString>(currentRcsmDir->GetName()).CompareTo("..")) {
            continue;
        }

        TString pathToRcmsDir = pathToSimDir + currentRcsmDir->GetName() + "/";
        TSystemDirectory RcmsDir("RcmsDir", pathToRcmsDir);

        TList *listOfRunDir = RcmsDir.GetListOfFiles();
        if (!listOfRunDir) return;

        TIter nextRunDir(listOfRunDir);
        TObject *currentRunDir = nullptr;

        while ((currentRunDir = nextRunDir())) {

            if (!static_cast<TString>(currentRunDir->GetName()).CompareTo(".") ||
                !static_cast<TString>(currentRunDir->GetName()).CompareTo("..")) {
                continue;
            }

            TString pathToRunDir = pathToRcmsDir + currentRunDir->GetName();
            TString pathToSimLog = pathToRunDir + "/sim.log";

            std::ifstream SimLog(pathToSimLog);
            if (!SimLog.is_open()) {
                std::cerr << "Unable to open file: " << pathToSimLog << std::endl;
                continue;
            }

            while (std::getline(SimLog, str_line)) {

                TString tstr_line = str_line;

                // a new event has appeared
                if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) {

                    // fill variables
                    eventID = current_eventID;
                    runNumber = static_cast<TString>(currentRunDir->GetName()).Atoi();
                    reactionChannel = static_cast<TString>(currentRcsmDir->GetName())(0);
                    SimpleTree->Fill();

                    // clear vectors
                    reactionID.clear();
                    // SM.clear();
                    SPx.clear();
                    SPy.clear();
                    SPz.clear();
                    NPDGCode.clear();
                    NPx.clear();
                    NPy.clear();
                    NPz.clear();

                    // start new event
                    current_eventID++;
                }

                if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;

                std::cout << tstr_line << std::endl;  // debug

                csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
                csv_arr = csv.Tokenize(",");

                reactionID.push_back(dynamic_cast<TObjString *>(csv_arr->At(0))->String().Atoi());
                // SM.push_back();
                SPx.push_back(dynamic_cast<TObjString *>(csv_arr->At(1))->String().Atof());
                SPy.push_back(dynamic_cast<TObjString *>(csv_arr->At(2))->String().Atof());
                SPz.push_back(dynamic_cast<TObjString *>(csv_arr->At(3))->String().Atof());
                NPDGCode.push_back(reactionChannel == 'A' ? 2112
                                                          : 2212);  // considering only ADEH // BUG it's using the prev. reac. channel
                NPx.push_back(dynamic_cast<TObjString *>(csv_arr->At(4))->String().Atof());
                NPy.push_back(dynamic_cast<TObjString *>(csv_arr->At(5))->String().Atof());
                NPz.push_back(dynamic_cast<TObjString *>(csv_arr->At(6))->String().Atof());

            }  // end of loop over lines

            SimLog.close();
        }  // end of loop over run dirs

    }  // end of loop over RCMS dirs

    SimpleTree->Write();

    SimpleFile->Save();
    SimpleFile->Close();

    std::cout << " " << std::endl;
}