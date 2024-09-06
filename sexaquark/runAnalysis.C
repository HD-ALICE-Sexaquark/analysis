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
                 TString LocalInputPath,  // (only valid when Mode == "local") dir that contains ProductionName dirs
                 Bool_t GridTestMode,     // (only valid when Mode == "grid")
                 Bool_t IsMC,             //
                 TString ProductionName,  // for data: "LHC15o", "LHC18q", "LHC18r"
                                          // for signal MC: "LHC23l1a3", "LHC23l1b3"
                                          // for gen. purp. MC: "LHC20e3a", "LHC20j6a"
                 TString RunNumbersList,  // path to file with run numbers
                 TString SourceOfV0s,     // "kalman", "custom", "on-the-fly", "offline"
                 TString SimulationSet,   // format: "<A,D,E,H><1.73,1.8,1.87,1.94,2.01>" e.g. "A1.73"
                 TString StopAfter,       // "MC", "Tracks", "Findables", "V0s"
                 Bool_t DoQA,             //
                 Bool_t ReadSignalLogs,   // (only valid when analyzing signal MC)
                 Int_t ChooseNEvents = 0) {

    /* Check for Input Errors */

    if (Mode != "local" && Mode != "grid") {
        std::cerr << "!! runAnalysis.C !! Error !! Invalid mode " << Mode << " !!" << std::endl;
        return;
    }

    if (Mode == "local") {
        if (LocalInputPath == "") {
            std::cerr << "!! runAnalysis.C !! Error !! LocalInputPath cannot be empty on local mode !!" << std::endl;
            return;
        }
        if (GridTestMode) {
            std::cerr << "!! runAnalysis.C !! Error !! GridTestMode is only valid on local mode !!" << std::endl;
            return;
        }
    }

    if (Mode == "grid") {
        if (LocalInputPath != "") {
            std::cerr << "!! runAnalysis.C !! Error !! LocalInputPath is only valid on local mode !!" << std::endl;
            return;
        }
        if (ChooseNEvents) {
            std::cerr << "!! runAnalysis.C !! Error !! ChooseNEvents is only valid on local mode !!" << std::endl;
            return;
        }
    }

    if (IsMC && (ProductionName == "LHC15o" || ProductionName == "LHC18q" || ProductionName == "LHC18r")) {
        std::cerr << "!! runAnalysis.C !! Error !! Make sure to put a valid ProductionName for simulations !!" << std::endl;
        return;
    }

    if (!IsMC) {
        if (DoQA || ReadSignalLogs) {
            std::cerr << "!! runAnalysis.C !! Error !! DoQA and ReadSignalLogs are only valid for MC !!" << std::endl;
            return;
        }
        if (ProductionName != "LHC15o" && ProductionName != "LHC18q" && ProductionName != "LHC18r") {
            std::cerr << "!! runAnalysis.C !! Error !! Make sure to put a valid ProductionName for data !!" << std::endl;
            return;
        }
    }

    if (!ProductionName.Contains("LHC23l1") && ReadSignalLogs) {
        std::cerr << "!! runAnalysis.C !! Error !! ReadSignalLogs is only valid for signal MC !!" << std::endl;
        return;
    }

    std::cout << "!! runAnalysis.C !! Passed initial input checks !!" << std::endl;

    /* Determine Further Options */

    Int_t PassNumber = 3;  // default for 18qr and anchored sims
    if (ProductionName == "LHC15o" || ProductionName == "LHC23l1b3" || ProductionName == "LHC20j6a") PassNumber = 2;

    Char_t ReactionID = SimulationSet[0];
    Float_t SexaquarkMass = ((TString)SimulationSet(1, SimulationSet.Length() - 1)).Atof();

    TString ProductionYear = "20" + ProductionName(3, 2);

    TString LocalDataDir = Form("%s/%s", LocalInputPath.Data(), ProductionName.Data());
    if (ProductionName.Contains("23l1")) LocalDataDir += Form("/%s", SimulationSet.Data());  // signal MC

    TString GridDataDir = Form("/alice/sim/%s/%s", ProductionYear.Data(), ProductionName.Data());
    if (ProductionName.Contains("23l1")) GridDataDir += Form("/%s", SimulationSet.Data());  // signal MC
    TString GridDataPattern = "/*/AliESDs.root";

    TString GridWorkingDir = "work/Sexaquark";
    if (!IsMC)
        GridWorkingDir += Form("_Data_Channel%c", ReactionID);
    else if (ProductionName.Contains("23l1"))
        GridWorkingDir += Form("_MC_%s_%s", ProductionName.Data(), SimulationSet.Data());
    else
        GridWorkingDir += Form("_MC_%s_Channel%c", ProductionName.Data(), ReactionID);
    TString GridOutputDir = "output";

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

    /* Start */

    gInterpreter->ProcessLine(".include .");
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    gInterpreter->ProcessLine(".include $KFPARTICLE_ROOT/include");

    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisManager_Sexaquark");

    std::cout << "!! runAnalysis.C !! Created AliAnalysisManager !!" << std::endl;

    /* Grid Connection */

    AliAnalysisAlien *alienHandler = nullptr;

    if (Mode == "grid") {
        alienHandler = new AliAnalysisAlien();
        alienHandler->SetCheckCopy(kFALSE);
        alienHandler->AddIncludePath(
            "-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$KFPARTICLE_ROOT/include");
        alienHandler->SetAdditionalLibs("AliAnalysisTaskSexaquark.cxx AliAnalysisTaskSexaquark.h AliAnalysisTaskSexaquark_Structs.h");
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
        alienHandler->SetJDLName("TaskSexaquark.jdl");
        alienHandler->SetExecutable("TaskSexaquark.sh");

        mgr->SetGridHandler(alienHandler);

        std::cout << "!! runAnalysis.C !! Passed grid connection !!" << std::endl;
    }

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

    TString TaskCentrality_Options = "";  // nothing
    AliMultSelectionTask *TaskCentrality = reinterpret_cast<AliMultSelectionTask *>(
        gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C" + TaskCentrality_Options));
    if (!TaskCentrality) return;

    TString TaskPIDResponse_Options = Form("(%i, 1, 1, \"%i\")", (Int_t)IsMC, PassNumber);
    AliAnalysisTaskPIDResponse *TaskPIDResponse = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
        gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C" + TaskPIDResponse_Options));
    if (!TaskPIDResponse) return;

    std::cout << "!! runAnalysis.C !! Passed addition of helper tasks !!" << std::endl;

    /* Add Sexaquark Task */

    gInterpreter->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");

    TString TaskSexaquark_Options = Form("(%i, \"%s\", \'%c\', %.2f, \"%s\", %i, %i)", (Int_t)IsMC, SourceOfV0s.Data(), ReactionID, SexaquarkMass,
                                         StopAfter.Data(), (Int_t)DoQA, (Int_t)ReadSignalLogs);
    AliAnalysisTaskSexaquark *TaskSexaquark =
        reinterpret_cast<AliAnalysisTaskSexaquark *>(gInterpreter->ExecuteMacro("AddSexaquark.C" + TaskSexaquark_Options));
    if (!TaskSexaquark) return;

    std::cout << "!! runAnalysis.C !! Passed addition of main task !!" << std::endl;

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
            // alienHandler->SetNtestFiles(5);               // TEST
            // alienHandler->SetSplitMaxInputFileNumber(5);  // TEST
            alienHandler->SetSplitMaxInputFileNumber(35);
            alienHandler->SetRunMode("full");
        }
        mgr->StartAnalysis("grid");
    } else {  // local mode
        TChain *chain = new TChain("esdTree");
        TString FilePath;
        for (Int_t &RN : GridRunNumbers) {
            for (Int_t DN = 0; DN < 6; DN++) {  // local directories (hardcoded)
                FilePath = Form("%s/%i/%03i/AliESDs.root", LocalDataDir.Data(), RN, DN);
                chain->AddFile(FilePath);
            }
        }
        if (!ChooseNEvents)
            mgr->StartAnalysis("local", chain);  // read all events
        else
            mgr->StartAnalysis("local", chain, ChooseNEvents);  // read first NEvents
    }

    std::cout << "!! runAnalysis.C !! Passed StartAnalysis !!" << std::endl;
}
