#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSexaquark.h"

AliAnalysisTaskSexaquark *AddSexaquark(Bool_t IsMC = kTRUE, Bool_t IsSignalMC = kTRUE, TString SourceOfV0s = "kalman", Bool_t DoQA = kFALSE) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return nullptr;

    AliAnalysisTaskSexaquark *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (!task) return nullptr;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);  // from `AliPhysicsSelectionTask`

    task->IsMC(IsMC, IsSignalMC);
    task->SetSourceOfV0s(SourceOfV0s);
    task->DoQA(DoQA);
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
    mgr->ConnectOutput(task, 6, mgr->CreateContainer("Sexaquarks_ALK0", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 7, mgr->CreateContainer("Sexaquarks_ALPK", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 8, mgr->CreateContainer("Sexaquarks_ALPKPP", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 9, mgr->CreateContainer("Sexaquarks_PKPKX", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 10, mgr->CreateContainer("Sexaquarks_LK0", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 11, mgr->CreateContainer("Sexaquarks_LNK", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 12, mgr->CreateContainer("Sexaquarks_LNKPP", TTree::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 13, mgr->CreateContainer("Sexaquarks_NKNKX", TTree::Class(), output_container, filename.Data()));

    mgr->ConnectOutput(task, 14, mgr->CreateContainer("QA_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 15, mgr->CreateContainer("AntiProton_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 16, mgr->CreateContainer("Proton_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 17, mgr->CreateContainer("NegKaon_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 18, mgr->CreateContainer("PosKaon_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 19, mgr->CreateContainer("PiMinus_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 20, mgr->CreateContainer("PiPlus_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 21, mgr->CreateContainer("AntiLambda_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 22, mgr->CreateContainer("Lambda_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 23, mgr->CreateContainer("KaonZeroShort_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 24, mgr->CreateContainer("PionPair_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 25, mgr->CreateContainer("Sexaquark_ALK0_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 26, mgr->CreateContainer("Sexaquark_ALPK_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 27, mgr->CreateContainer("Sexaquark_ALPKPP_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 28, mgr->CreateContainer("Sexaquark_PKPKX_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 29, mgr->CreateContainer("Sexaquark_LK0_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 30, mgr->CreateContainer("Sexaquark_LNK_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 31, mgr->CreateContainer("Sexaquark_LNKPP_Hists", TList::Class(), output_container, filename.Data()));
    mgr->ConnectOutput(task, 32, mgr->CreateContainer("Sexaquark_NKNKX_Hists", TList::Class(), output_container, filename.Data()));

    return task;
}
