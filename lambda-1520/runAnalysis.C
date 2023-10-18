#include "AliAnalysisTaskLambda1520Lpipi.h"

void runAnalysis(Bool_t IsMC, Int_t NEvents = 0) {

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
  // AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE, Bool_t tuneOnData=kTRUE, TString recoPass="2", etc...)
  TString pid_response_path = gSystem->ExpandPathName("${ALICE_ROOT}/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(gInterpreter->ExecuteMacro(
      Form("%s(%i, 1, 1, \"3\")", pid_response_path.Data(), (Int_t)IsMC)));  // PENDING: check what happens with the pass!!

  // load task
  gInterpreter->LoadMacro("AliAnalysisTaskLambda1520Lpipi.cxx++g");
  AliAnalysisTaskLambda1520Lpipi *task = reinterpret_cast<AliAnalysisTaskLambda1520Lpipi *>(
      gInterpreter->ExecuteMacro(Form("AddTaskLambda1520Lpipi.C(\"name\", %i)", (Int_t)IsMC)));

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

  // start the analysis locally, reading the events from the TChain
  if (NEvents) {
    mgr->StartAnalysis("local", chain, NEvents);
  } else {
    mgr->StartAnalysis("local", chain);
  }
}
