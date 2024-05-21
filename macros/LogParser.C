#include "include/Headers.hxx"

/*
 Read the `sim.log` files, to extract:
 - the injected anti-sexaquark properties
 - and the struck nucleon kinematics
 prior to their reaction, and store them into a simple tree
*/
void LogParser(TString DataSet = "LHC23l1a3",  // could also be "LHC23l1b"
               TString SimulationSet = "A1.8") {

    /* Hard-coded parameters */

    TString pathToSimDir = "/home/ceres/borquez/some/sims/" + DataSet + "/" + SimulationSet;
    TSystemDirectory SimDir("SimDir", pathToSimDir);

    /* Define output file */

    TFile *SimpleFile = new TFile(DataSet + ".root", "RECREATE");

    /* Define tree structure */

    TTree *SimpleTree = new TTree("IC", "IC");

    Char_t reactionChannel;
    Double_t SM;
    Int_t runNumber;
    Int_t subDirNumber;
    Int_t eventID = -1;
    Int_t reactionID;
    Int_t NPDGCode;
    Double_t SPx, SPy, SPz;
    Double_t NPx, NPy, NPz;

    SimpleTree->Branch("reactionChannel", &reactionChannel);
    SimpleTree->Branch("SM", &SM);
    SimpleTree->Branch("runNumber", &runNumber);
    SimpleTree->Branch("subDirNumber", &subDirNumber);
    SimpleTree->Branch("eventID", &eventID);
    SimpleTree->Branch("reactionID", &reactionID);
    SimpleTree->Branch("NPDGCode", &NPDGCode);
    SimpleTree->Branch("SPx", &SPx);
    SimpleTree->Branch("SPy", &SPy);
    SimpleTree->Branch("SPz", &SPz);
    SimpleTree->Branch("NPx", &NPx);
    SimpleTree->Branch("NPy", &NPy);
    SimpleTree->Branch("NPz", &NPz);

    // auxiliary variables
    std::string str_line;
    TString csv;
    TObjArray *csv_arr = nullptr;

    /* Loop over run numbers */

    TList *listOfRunDirs = SimDir.GetListOfFiles();
    if (!listOfRunDirs) return;

    TIter next(listOfRunDirs);
    while (TSystemDirectory *currentRunDir = (TSystemDirectory *)next()) {

        TString runDir = currentRunDir->GetName();
        if (!runDir.CompareTo(".") || !runDir.CompareTo("..")) continue;

        TSystemDirectory RunDir("RunDir", currentRunDir->GetTitle());

        TList *listOfSubDirs = RunDir.GetListOfFiles();
        if (!listOfSubDirs) return;

        /* Loop over sub directories */

        TIter nextt(listOfSubDirs);
        while (TSystemDirectory *currentSubDir = (TSystemDirectory *)nextt()) {

            TString subDir = currentSubDir->GetName();
            if (!subDir.CompareTo(".") || !subDir.CompareTo("..")) continue;

            TString subDirPath = currentSubDir->GetTitle();

            /* Open file */

            std::ifstream SimLog(subDirPath + "/sim.log");
            if (!SimLog.is_open()) {
                std::cerr << "Unable to open file: " << subDirPath + "/sim.log" << std::endl;
                continue;
            }

            /* Updated file and path-related variables */

            reactionChannel = SimulationSet(0);
            SM = ((TString)SimulationSet(1, SimulationSet.Length())).Atof();
            runNumber = runDir.Atoi();
            subDirNumber = subDir.Atoi();

            /* Read lines */

            while (std::getline(SimLog, str_line)) {

                TString tstr_line = str_line;

                // a new event has appeared
                if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) eventID++;

                if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;

                csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
                csv_arr = csv.Tokenize(",");

                reactionID = dynamic_cast<TObjString *>(csv_arr->At(0))->String().Atoi();
                NPDGCode = reactionChannel == 'A' ? 2112 : 2212;  // considering only ADEH
                SPx = dynamic_cast<TObjString *>(csv_arr->At(1))->String().Atof();
                SPy = dynamic_cast<TObjString *>(csv_arr->At(2))->String().Atof();
                SPz = dynamic_cast<TObjString *>(csv_arr->At(3))->String().Atof();
                NPx = dynamic_cast<TObjString *>(csv_arr->At(4))->String().Atof();
                NPy = dynamic_cast<TObjString *>(csv_arr->At(5))->String().Atof();
                NPz = dynamic_cast<TObjString *>(csv_arr->At(6))->String().Atof();

                std::cout << reactionChannel << " " << SM << " " << runNumber << " " << subDirNumber << " " << eventID << " " << reactionID << " "
                          << NPDGCode << " " << SPx << " " << SPy << " " << SPz << " " << NPx << " " << NPy << " " << NPz << std::endl;  // debug

                SimpleTree->Fill();
            }  // end of loop over lines

            SimLog.close();
        }  // end of loop over run dirs

    }  // end of loop over run dirs

    /* Write tree onto file, then save, then close */

    SimpleTree->Write();

    SimpleFile->Save();
    SimpleFile->Close();
}
