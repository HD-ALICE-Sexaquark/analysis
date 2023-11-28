AliAnalysisTaskSexaquark* AddSexaquark(TString name = "name", Bool_t IsMC = kTRUE, Bool_t HasSexaquark = kTRUE,
                                       TString SourceOfV0s = "official", TString ReactionChannel = "A") {
    // get the manager via the static access member. since it's static, you don't need
    // to create an instance of the class here to call the function
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }

    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    // fileName += ":Sexaquark";  // create a subfolder in the file

    // now we create an instance of your task
    AliAnalysisTaskSexaquark* task =
        new AliAnalysisTaskSexaquark(name.Data(), IsMC, HasSexaquark, (TString)SourceOfV0s, ReactionChannel.Data()[0]);

    if (!task) {
        return 0x0;
    }

    // trigger configuration, commented on purpose
    // task->SelectCollisionCandidates(AliVEvent::kAnyINT);

    // add your task to the manager
    mgr->AddTask(task);
    // connect task's to the manager
    mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task, 1, mgr->CreateContainer("Trees", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
