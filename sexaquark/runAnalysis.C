#include <fstream>

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"

#include "AliAnalysisTaskPIDResponse.h"
#include "AliMultSelectionTask.h"
#include "AliPhysicsSelectionTask.h"

#include "AliAnalysisTaskSexaquark.h"

void runAnalysis(TString Mode,            // "local", "grid"
                 Bool_t GridTestMode,     // (only valid when Mode == "grid")
                 Bool_t IsMC,             //
                 TString ProductionName,  // for data: "LHC15o", "LHC18q", "LHC18r"
                                          // for signal MC: "LHC23l1a3", "LHC23l1b3"
                                          // for gen. purp. MC: "LHC20e3a", "LHC20j6a"
                 TString RunNumbersList,  // path to file with run numbers
                 TString SourceOfV0s,     // "kalman", "custom", "on-the-fly", "offline"
                 TString SimulationSet,   // format: "<A,D,E,H><1.73,1.8,1.87,1.94,2.01>" e.g. "A1.73"
                 Bool_t DoQA,             //
                 Bool_t ReadSignalLogs,   // (only valid when analyzing signal MC)
                 Bool_t ReweightPt,       // (only valid when analyzing signal MC and ReadSignalLogs enabled)
                 Bool_t ReweightRadius,   // (only valid when analyzing signal MC)
                 Int_t ChooseNEvents = 0) {

    /* Determine Further Options */

    Int_t PassNumber = 3;  // default for 18qr and anchored sims
    if (ProductionName == "LHC15o" || ProductionName == "LHC23l1b3" || ProductionName == "LHC20j6a") PassNumber = 2;

    Char_t ReactionID = SimulationSet[0];
    Float_t SexaquarkMass = ((TString)SimulationSet(1, SimulationSet.Length() - 1)).Atof();

    TString ProductionYear = "20" + ProductionName(3, 2);

    TString GridDataDir = Form("/alice/sim/%s/%s", ProductionYear.Data(), ProductionName.Data());
    if (ProductionName.Contains("23l1")) GridDataDir += Form("/%s", SimulationSet.Data());  // signal MC
    TString GridDataPattern = "/*/AliESDs.root";

    TString GridWorkingDir = "work";
    TString GridOutputDir = "sexaquark";
    if (!IsMC)
        GridOutputDir += "_data";
    else if (ProductionName.Contains("23l1"))
        GridOutputDir += Form("_signal_mc_%s", SimulationSet.Data());
    else
        GridOutputDir += "_bkg_mc";

    std::vector<Int_t> GridRunNumbers;
    std::ifstream RunNumbersFile(RunNumbersList);
    if (!RunNumbersFile.is_open()) {
        std::cerr << "!! runAnalysis.C !! Error !! Unable to open file " << RunNumbersList << " !!" << std::endl;
        return;
    }
    Int_t SingleRN;
    while (RunNumbersFile >> SingleRN) GridRunNumbers.push_back(SingleRN);
    RunNumbersFile.close();

    std::cout << "!! runAnalysis.C !! Passed further options !!" << std::endl;

    /* Check for Input Errors */

    if (!IsMC && (DoQA || ReadSignalLogs || ReweightPt || ReweightRadius)) {
        std::cerr << "!! runAnalysis.C !! Error !! DoQA, ReadSignalLogs, ReweightPt, and ReweightRadius are only valid for MC !!" << std::endl;
        return;
    }

    if (!ProductionName.Contains("LHC23l1") && (ReadSignalLogs || ReweightPt || ReweightRadius)) {
        std::cerr << "!! runAnalysis.C !! Error !! ReadSignalLogs, ReweightPt, and ReweightRadius are only valid for signal MC !!" << std::endl;
        return;
    }

    if (!IsMC && (ProductionName != "LHC15o" && ProductionName != "LHC18q" && ProductionName != "LHC18r")) {
        std::cerr << "!! runAnalysis.C !! Error !! Make sure to put a valid ProductionName for data !!" << std::endl;
        return;
    }

    if (ReweightPt && !ReadSignalLogs) {
        std::cerr << "!! runAnalysis.C !! Error !! ReweightPt is only valid when ReadSignalLogs is enabled !!" << std::endl;
        return;
    }

    std::cout << "!! runAnalysis.C !! Passed input checks !!" << std::endl;

    /* Start */

    gInterpreter->ProcessLine(".include .");
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    gInterpreter->ProcessLine(".include $KFPARTICLE_ROOT/include");

    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager_Sexaquark");

    std::cout << "!! runAnalysis.C !! Passed start !!" << std::endl;

    /* Grid Connection */

    AliAnalysisAlien *alienHandler = nullptr;

    if (Mode == "grid") {
        alienHandler = new AliAnalysisAlien();
        alienHandler->SetCheckCopy(kFALSE);
        alienHandler->AddIncludePath(
            "-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$KFPARTICLE_ROOT/include");
        alienHandler->SetAdditionalLibs("AliAnalysisTaskSexaquark.cxx AliAnalysisTaskSexaquark.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskSexaquark.cxx");
        alienHandler->SetAliPhysicsVersion("vAN-20240807_O2-1");
        alienHandler->SetExecutableCommand("aliroot -l -q -b");
        alienHandler->SetGridDataDir(GridDataDir);
        if (!IsMC) alienHandler->SetRunPrefix("000");
        for (Int_t &RN : GridRunNumbers) alienHandler->AddRunNumber(RN);
        alienHandler->SetDataPattern(GridDataPattern);
        alienHandler->SetTTL(3600);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMergeViaJDL(kFALSE);
        // alienHandler->SetMaxMergeStages(1);
        alienHandler->SetGridWorkingDir(GridWorkingDir);
        alienHandler->SetGridOutputDir(GridOutputDir);
        alienHandler->SetAnalysisMacro("TaskSexaquark.C");
        alienHandler->SetJDLName("TaskSexaquark.jdl");
        alienHandler->SetExecutable("TaskSexaquark.sh");
        mgr->SetGridHandler(alienHandler);
    }

    std::cout << "!! runAnalysis.C !! Passed grid connection !!" << std::endl;

    /* Input Handlers */

    AliESDInputHandler *esdH = new AliESDInputHandler();
    esdH->SetNeedField();  // necessary to get GoldenChi2
    mgr->SetInputEventHandler(esdH);

    AliMCEventHandler *mcH = nullptr;
    if (IsMC) {
        mcH = new AliMCEventHandler();
        mcH->SetReadTR(kFALSE);
        mgr->SetMCtruthEventHandler(mcH);
    }

    std::cout << "!! runAnalysis.C !! Passed creation of input handlers !!" << std::endl;

    /* Add Helper Tasks */

    TString TaskPhysicsSelection_Options = Form("(%i)", (Int_t)IsMC);
    AliPhysicsSelectionTask *TaskPhysicsSelection = reinterpret_cast<AliPhysicsSelectionTask *>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C" + TaskPhysicsSelection_Options));
    if (!TaskPhysicsSelection) return;

    TString TaskPIDResponse_Options = Form("(%i, 1, 1, \"%i\")", (Int_t)IsMC, PassNumber);
    AliAnalysisTaskPIDResponse *TaskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
        gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C" + TaskPIDResponse_Options));
    if (!TaskPIDResponse) return;

    std::cout << "!! runAnalysis.C !! Passed adding of helper tasks !!" << std::endl;

    /* Add Sexaquark Task */

    gInterpreter->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");

    TString TaskSexaquark_Options = Form("(%i, \"%s\", \'%c\', %.2f, %i, %i, %i, %i)", (Int_t)IsMC, SourceOfV0s.Data(), ReactionID, SexaquarkMass,
                                         (Int_t)DoQA, (Int_t)ReadSignalLogs, (Int_t)ReweightPt, (Int_t)ReweightRadius);
    AliAnalysisTaskSexaquark *TaskSexaquark =
        reinterpret_cast<AliAnalysisTaskSexaquark *>(gInterpreter->ExecuteMacro("AddSexaquark.C" + TaskSexaquark_Options));
    if (!TaskSexaquark) return;

    std::cout << "!! runAnalysis.C !! Passed adding of main task !!" << std::endl;

    /* Init Analysis Manager */

    mgr->SetDebugLevel(3);
    if (!mgr->InitAnalysis()) return;

    std::cout << "!! runAnalysis.C !! Passed InitAnalysis !!" << std::endl;

    /* Start Analysis */

    if (Mode == "grid") {
        if (GridTestMode) {
            alienHandler->SetNtestFiles(5);
            alienHandler->SetRunMode("test");
        } else {
            alienHandler->SetSplitMaxInputFileNumber(55);
            alienHandler->SetRunMode("full");
        }
        mgr->StartAnalysis("grid");
    } else if (Mode == "local") {
        /* WORK IN PROGRESS */
        TChain *chain = new TChain("esdTree");
        if (!ChooseNEvents) {
            std::vector<TString> localFiles = {"/home/ceres/borquez/some/sims/LHC20e3a/297595/001/AliESDs.root"};
            for (Int_t ifile = 0; ifile < (Int_t)localFiles.size(); ifile++) {
                chain->AddFile(localFiles[ifile]);
            }
            mgr->StartAnalysis("local", chain);  // read all events
        } else {
            mgr->StartAnalysis("local", chain, ChooseNEvents);  // read first NEvents
        }
    } else {
        std::cerr << "!! runAnalysis.C !! Error !! Invalid mode " << Mode << " !!" << std::endl;
        return;
    }

    std::cout << "!! runAnalysis.C !! Passed StartAnalysis !!" << std::endl;
}
