#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSexaquark.h"

AliAnalysisTaskSexaquark *AddSexaquark(Bool_t IsMC = kTRUE, TString SourceOfV0s = "kalman", Bool_t DoQA = kFALSE, Bool_t ReadSignalLogs = kTRUE) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return nullptr;

    AliAnalysisTaskSexaquark *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (!task) return nullptr;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);  // from `AliPhysicsSelectionTask`

    task->IsMC(IsMC);
    task->SetSourceOfV0s(SourceOfV0s);
    task->DoQA(DoQA);
    task->ReadSignalLogs(ReadSignalLogs);
    task->Initialize();

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    TString filename = AliAnalysisManager::GetCommonFileName();
    AliAnalysisManager::EAliAnalysisContType output_container = AliAnalysisManager::kOutputContainer;

    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Events", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("Injected", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 3, mgr->CreateContainer("MC", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 4, mgr->CreateContainer("Tracks", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 5, mgr->CreateContainer("V0s", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 6, mgr->CreateContainer("KaonPairs", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 7, mgr->CreateContainer("Sexaquarks_ALPK", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 8, mgr->CreateContainer("Sexaquarks_ALK0", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 9, mgr->CreateContainer("Sexaquarks_ALKKX", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 10, mgr->CreateContainer("Sexaquarks_ALPKPP", TTree::Class(), output_container, filename.Data()));

    mgr->ConnectOutput(task, 11, mgr->CreateContainer("QA_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 12, mgr->CreateContainer("Tracks_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 13, mgr->CreateContainer("V0s_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 14, mgr->CreateContainer("KaonPairs_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 15, mgr->CreateContainer("Sexaquarks_ALPK_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 16, mgr->CreateContainer("Sexaquarks_ALK0_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 17, mgr->CreateContainer("Sexaquarks_ALKKX_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 18, mgr->CreateContainer("Sexaquarks_ALPKPP_Hists", TList::Class(), output_container, filename.Data()));

    return task;
}
