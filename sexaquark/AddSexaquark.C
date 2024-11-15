AliAnalysisTaskSexaquark *AddSexaquark(Bool_t IsMC = kTRUE, TString SourceOfV0s = "kalman", Char_t ReactionChannel = 'A', Float_t SexaquarkMass = 1.8,
                                       Bool_t DoQA = kFALSE, Bool_t ReadSignalLogs = kTRUE) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return nullptr;

    AliAnalysisTaskSexaquark *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (!task) return nullptr;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);  // from `AliPhysicsSelectionTask`

    task->IsMC(IsMC);
    task->SetSourceOfV0s(SourceOfV0s);
    task->SetReactionChannel(ReactionChannel);
    task->DoQA(DoQA);
    task->ReadSignalLogs(ReadSignalLogs);
    task->Initialize();

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    TString filename = AliAnalysisManager::GetCommonFileName();

    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Tree_Events", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("Tree_Injected", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 3, mgr->CreateContainer("Tree_MC", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 4, mgr->CreateContainer("Tree_Tracks", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 5, mgr->CreateContainer("Tree_V0s", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 6, mgr->CreateContainer("Tree_Sexaquarks", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 7, mgr->CreateContainer("Tree_KaonPairs", TTree::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));

    mgr->ConnectOutput(task, 8, mgr->CreateContainer("QA_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 9, mgr->CreateContainer("Tracks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 10, mgr->CreateContainer("V0s_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 11, mgr->CreateContainer("Sexaquarks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 12, mgr->CreateContainer("PosKaonPairs_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));

    return task;
}
