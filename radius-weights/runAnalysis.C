#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"

#include "AliQuickTaskRadiusWeights.h"

void runAnalysis(Int_t ChooseNEvents = 0) {

    // tell root where to look for headers
    gInterpreter->ProcessLine(".include ${ROOTSYS}/include");
    gInterpreter->ProcessLine(".include ${ALICE_ROOT}/include");
    gInterpreter->ProcessLine(".include ${ALICE_PHYSICS}/include");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");

    AliESDInputHandler *esdH = new AliESDInputHandler();
    mgr->SetInputEventHandler(esdH);
    esdH->SetNeedField();  // necessary to get GoldenChi2

    AliMCEventHandler *mcHand = nullptr;
    mcHand = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHand);

    // load event selection task
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");

    // load task
    gInterpreter->LoadMacro("AliQuickTaskRadiusWeights.cxx++g");

    AliQuickTaskRadiusWeights *task = reinterpret_cast<AliQuickTaskRadiusWeights *>(gInterpreter->ExecuteMacro("AddRadiusWeights.C"));

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
