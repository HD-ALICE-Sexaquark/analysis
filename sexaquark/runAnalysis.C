#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskPIDResponse.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"

#include "AliAnalysisTaskSexaquark.h"

void runAnalysis(Bool_t IsMC, TString RunPeriod, TString SourceOfV0s, TString SimulationSet, Bool_t ReadSignalLogs, Bool_t DoQA,
                 Int_t ChooseNEvents = 0) {

    // tell root where to look for headers
    gInterpreter->ProcessLine(".include ${ROOTSYS}/include");
    gInterpreter->ProcessLine(".include ${ALICE_ROOT}/include");
    gInterpreter->ProcessLine(".include ${ALICE_PHYSICS}/include");
    gInterpreter->ProcessLine(".include ${KFPARTICLE_ROOT}/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");

    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    esdH->SetNeedField();  // necessary to get GoldenChi2

    AliMCEventHandler *mcHand = nullptr;
    if (IsMC) {
        mcHand = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcHand);
    }

    // choose options depending on data period
    Int_t n_pass;
    if (RunPeriod == "LHC15o") {
        n_pass = 2;
    } else if (RunPeriod == "LHC18q" || RunPeriod == "LHC18r") {
        n_pass = 3;
    } else {
        printf("Data period not recognized");
        return;
    }

    // CDB CONNECT
    /*
    // - to use this, one needs to clone https://gitlab.cern.ch/alisw/AliRootOCDB
    // - then, one needs to set the environment variable ALIROOT_OCDB_ROOT to the directory where the repo is cloned
    AliTaskCDBconnect *taskCDB = reinterpret_cast<AliTaskCDBconnect *>(
        gInterpreter->ExecuteMacro("${ALICE_PHYSICS}/PWGPP/PilotTrain/AddTaskCDBconnect.C(\"local://${ALIROOT_OCDB_ROOT}/OCDB\")"));
    taskCDB->SetFallBackToRaw(kTRUE);
    */

    // load event selection task
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");

    // load PID task
    TString pid_response_path = gSystem->ExpandPathName("${ALICE_ROOT}/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
        gInterpreter->ExecuteMacro(Form("%s(%i, 1, 1, \"%i\")", pid_response_path.Data(), (Int_t)IsMC, n_pass)));

    // load sexaquark task
    gInterpreter->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");

    Char_t ReactionID = SimulationSet(0);
    Float_t SexaquarkMass = ((TString)SimulationSet(1, SimulationSet.Length())).Atof();
    AliAnalysisTaskSexaquark *task = reinterpret_cast<AliAnalysisTaskSexaquark *>(
        gInterpreter->ExecuteMacro(Form("AddSexaquark.C(%i, \"%s\", \'%c\', %.2f, %i, %i)", (Int_t)IsMC, SourceOfV0s.Data(), ReactionID,
                                        SexaquarkMass, (Int_t)ReadSignalLogs, (Int_t)DoQA)));

    if (!mgr->InitAnalysis()) {
        return;
    }

    // print info on screen
    mgr->SetDebugLevel(3);
    mgr->PrintStatus();

    // if you want to run locally, we need to define some input
    TChain *chain = new TChain("esdTree");
    // add a few files to the chain (change this so that your local files are added)
    chain->Add("AliESDs.root");

    if (!ChooseNEvents)
        mgr->StartAnalysis("local", chain);  // read all events
    else
        mgr->StartAnalysis("local", chain, ChooseNEvents);  // read the first NEvents

    // clean memory
    delete mgr;
    delete esdH;
    delete mcHand;
    delete chain;
}
