#include "AliAnalysisTaskSexaquark.h"

void runAnalysis(Bool_t IsMC, Bool_t HasSexaquark, TString SourceOfV0s, TString ReactionChannel) {

  // since we will compile a class, tell root where to look for headers
#if !defined(__CINT__) || defined(__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

  // create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");

  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

  AliMCEventHandler *mcHand = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHand);

  // compile the class and load the add task macro
  // here we have to differentiate between using the just-in-time compiler
  // from root6, or the interpreter of root5
#if !defined(__CINT__) || defined(__CLING__)
  // load PID task
  TMacro PIDadd(gSystem->ExpandPathName("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
  AliAnalysisTaskPIDResponse *PIDresponseTask = reinterpret_cast<AliAnalysisTaskPIDResponse *>(PIDadd.Exec());
  // load sexaquark task
  gInterpreter->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");
  AliAnalysisTaskSexaquark *task = reinterpret_cast<AliAnalysisTaskSexaquark *>(gInterpreter->ExecuteMacro(Form(
      "AddSexaquark.C(\"name\", %i, %i, \"%s\", \"%s\")", (Int_t)IsMC, (Int_t)HasSexaquark, SourceOfV0s.Data(), ReactionChannel.Data())));
#else
  // load PID task
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AddTaskPIDResponse(IsMC);
  // load sexaquark task
  gROOT->LoadMacro("AliAnalysisTaskSexaquark.cxx++g");
  gROOT->LoadMacro("AddSexaquark.C");
  AliAnalysisTaskSexaquark *task = AddSexaquark();
#endif

  // vital!
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
  mgr->StartAnalysis("local", chain);
}
