AliAnalysisTaskSexaquark *AddSexaquark(Bool_t IsMC = kTRUE, TString SourceOfV0s = "kalman", Char_t ReactionID = 'A', Float_t SexaquarkMass = 1.8,
                                       TString StopAfter = "", Bool_t DoQA = kFALSE, Bool_t ReadSignalLogs = kTRUE) {

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) return nullptr;

    AliAnalysisTaskSexaquark *task = new AliAnalysisTaskSexaquark("AliAnalysisTaskSexaquark");
    if (!task) return nullptr;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);  // from `AliPhysicsSelectionTask`

    task->IsMC(IsMC);
    task->SetSourceOfV0s(SourceOfV0s);
    task->SetReactionID(ReactionID);
    task->DoQA(DoQA);
    task->StopAfter(StopAfter);
    task->ReadSignalLogs(ReadSignalLogs);
    task->Initialize();

    mgr->AddTask(task);

    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

    TString filename = AliAnalysisManager::GetCommonFileName();
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Trees", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 2, mgr->CreateContainer("QA_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 3, mgr->CreateContainer("Tracks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 4, mgr->CreateContainer("V0s_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 5, mgr->CreateContainer("Sexaquarks_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));
    mgr->ConnectOutput(task, 6, mgr->CreateContainer("PosKaonPairs_Hists", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));

    return task;
}
