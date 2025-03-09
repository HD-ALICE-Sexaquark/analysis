#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSexaquark.h"

AliAnalysisTaskSexaquark *AddTaskSexaquark(Bool_t is_mc, Bool_t is_signal_mc) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (mgr == nullptr) return nullptr;

    auto *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (task == nullptr) return nullptr;

    task->SelectCollisionCandidates(AliVEvent::kINT7);  // from `AliPhysicsSelectionTask`
    task->Initialize(is_mc, is_signal_mc);

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Events", TTree::Class(), mgr->kOutputContainer, mgr->GetCommonFileName()));

    return task;
}
