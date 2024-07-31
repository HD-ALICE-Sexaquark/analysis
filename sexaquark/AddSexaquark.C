/* --- */

TTree *ReadLogs(TString filename);
TH1D *GetPtWeights(TString filename, Float_t SexaquarkMass);
TH1F *GetRadiusWeights(TString filename);

/* --- */

AliAnalysisTaskSexaquark *AddSexaquark(Bool_t IsMC = kTRUE, TString SourceOfV0s = "offline", Char_t ReactionID = 'A', Float_t SexaquarkMass = 1.8,
                                       Bool_t ReadSignalLogs = kTRUE, Bool_t DoQA = kFALSE, Bool_t ReweightPt = kFALSE,
                                       Bool_t ReweightRadius = kFALSE) {
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return 0x0;

    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) return 0x0;

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    // fileName += ":Sexaquark";  // create a subfolder in the file

    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(IsMC);

    // now we create an instance of your task
    AliAnalysisTaskSexaquark *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (!task) return 0x0;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);

    task->IsMC(IsMC);
    task->SetSourceOfV0s(SourceOfV0s);
    task->SetReactionID(ReactionID);
    task->ReadSignalLogs(ReadSignalLogs);
    task->DoQA(DoQA);
    task->ReweightPt(ReweightPt);
    task->ReweightRadius(ReweightRadius);
    task->Initialize();

    TTree *CsvTree;
    if (ReadSignalLogs) {
        CsvTree = ReadLogs("sim.log");
        task->AddLogTree(CsvTree);
    }

    TH1D *PtWeightsHist = GetPtWeights("blastwave.root", SexaquarkMass);
    task->AddPtWeights(PtWeightsHist);

    TH1F *RadiusWeightsHist = GetRadiusWeights("radius-weights.root");
    task->AddRadiusWeights(RadiusWeightsHist);

    /* --- */

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Trees", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("QA_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 3, mgr->CreateContainer("Tracks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 4, mgr->CreateContainer("V0s_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 5, mgr->CreateContainer("Sexaquarks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task, 6, mgr->CreateContainer("PosKaonPairs_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}

/*
  Read a log file, extract initial conditions data as in CSV format and populate TTree.
*/
TTree *ReadLogs(TString filename) {

    TTree *t = new TTree("Injected", "Injected");  // i.e., before the antisexaquark-nucleon interaction

    Char_t ReactionChannel;
    Int_t RunNumber;

    Int_t EventID = -1;
    Int_t ReactionID;
    Int_t NPDGCode;
    Double_t SPx, SPy, SPz, SM;
    Double_t NPx, NPy, NPz;

    t->Branch("ReactionChannel", &ReactionChannel);
    t->Branch("RunNumber", &RunNumber);
    t->Branch("EventID", &EventID);
    t->Branch("ReactionID", &ReactionID);
    t->Branch("Nucl_PDGCode", &NPDGCode);
    t->Branch("Sexa_Px_ini", &SPx);
    t->Branch("Sexa_Py_ini", &SPy);
    t->Branch("Sexa_Pz_ini", &SPz);
    t->Branch("Sexa_M_ini", &SM);
    t->Branch("Fermi_Px", &NPx);
    t->Branch("Fermi_Py", &NPy);
    t->Branch("Fermi_Pz", &NPz);

    // auxiliary variables
    std::string cstr_line;
    TString tstr_line, csv;
    TObjArray *csv_arr = nullptr;

    std::ifstream SimLog(filename);
    if (!SimLog.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return nullptr;
    }

    /* Read lines */

    while (std::getline(SimLog, cstr_line)) {

        tstr_line = cstr_line;

        if (!ReactionChannel) {
            if (tstr_line.Contains("I-AliGenSexaquarkReaction::PrintParameters:   Chosen Reaction Channel : ")) {
                ReactionChannel = tstr_line(71);
            }
        }

        if (!SM) {
            if (tstr_line.Contains("I-AliGenSexaquarkReaction::PrintParameters:   * M = ")) {
                SM = TString(tstr_line(51, 59)).Atof();  // assuming M with length 8
            }
        }

        if (!RunNumber) {
            if (tstr_line.Contains("I-AliSimulation::ProcessEnvironmentVars: Run number = ")) {
                RunNumber = TString(tstr_line(54, 59)).Atoi();  // assuming RN with length 6
                printf("Gotten the following RN : %i", RunNumber);
            }
        }

        // a new event has appeared
        if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) EventID++;

        if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;

        csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
        csv_arr = csv.Tokenize(",");

        ReactionID = dynamic_cast<TObjString *>(csv_arr->At(0))->String().Atoi();
        NPDGCode = ReactionChannel == 'A' ? 2112 : 2212;  // considering only ADEH
        SPx = dynamic_cast<TObjString *>(csv_arr->At(1))->String().Atof();
        SPy = dynamic_cast<TObjString *>(csv_arr->At(2))->String().Atof();
        SPz = dynamic_cast<TObjString *>(csv_arr->At(3))->String().Atof();
        NPx = dynamic_cast<TObjString *>(csv_arr->At(4))->String().Atof();
        NPy = dynamic_cast<TObjString *>(csv_arr->At(5))->String().Atof();
        NPz = dynamic_cast<TObjString *>(csv_arr->At(6))->String().Atof();

        t->Fill();
    }  // end of loop over lines

    SimLog.close();

    return t;
}

/*
  Get a histogram with the pt weights for the blast wave.
*/
TH1D *GetPtWeights(TString filename, Float_t SexaquarkMass) {

    TFile *f = TFile::Open(filename);
    if (!f) {
        std::cerr << "GetPtWeights :: Unable to open file: " << filename << std::endl;
        return nullptr;
    }

    TString hist_name = Form("BlastWaveHist_%.2f", SexaquarkMass);
    TH1D *h = dynamic_cast<TH1D *>(f->Get(hist_name));
    if (!h) {
        std::cerr << "GetPtWeights :: Unable to find histogram " << hist_name << " in file " << filename << std::endl;
        return nullptr;
    }

    return h;
}

/*
  Get a histogram with the radius weights of photon conversions.
*/
TH1F *GetRadiusWeights(TString filename) {

    TFile *f = TFile::Open(filename);
    if (!f) {
        std::cerr << "GetRadiusWeights :: Unable to open file " << filename << std::endl;
        return nullptr;
    }

    TString hist_name = "MCGen_PhotonConversions_Radius";
    TH1F *h = dynamic_cast<TH1F *>(f->Get("Hists")->FindObject(hist_name));
    if (!h) {
        std::cerr << "GetRadiusWeights :: Unable to find histogram " << hist_name << " in file " << filename << std::endl;
        return nullptr;
    }

    return h;
}
