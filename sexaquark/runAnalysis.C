#include "AliAnalysisTaskSexaquark.h"

void runAnalysis(Bool_t IsMC, TString SourceOfV0s, TString ReactionChannel, Int_t NEvents = 0) {

    // tell root where to look for headers
    gInterpreter->ProcessLine(".include ${ROOTSYS}/include");
    gInterpreter->ProcessLine(".include ${ALICE_ROOT}/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");

    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);

    AliMCEventHandler *mcHand = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHand);

    // load PID task
    TString pid_response_path = gSystem->ExpandPathName("${ALICE_ROOT}/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(
        gInterpreter->ExecuteMacro(Form("%s(%i, 1, 1, \"3\")", pid_response_path.Data(), (Int_t)IsMC))); // (pending) pass number?

    // load sexaquark task
    gInterpreter->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");
    AliAnalysisTaskSexaquark *task = reinterpret_cast<AliAnalysisTaskSexaquark *>(gInterpreter->ExecuteMacro(Form(
        "AddSexaquark.C(\"name\", %i, \"%s\", \"%s\")", (Int_t)IsMC, SourceOfV0s.Data(), ReactionChannel.Data())));

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

    if (NEvents) {
        // read the first NEvents
        mgr->StartAnalysis("local", chain, NEvents);
    } else {
        // read all events
        mgr->StartAnalysis("local", chain);
    }
}
