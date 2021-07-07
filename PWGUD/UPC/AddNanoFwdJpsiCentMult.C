AliAnalysisTaskNanoFwdJpsiCentMult* AddNanoFwdJpsiCentMult(const char* suffix = "")
{
TString name = "NanoFwdJpsiCentMult";
name += suffix;     
// get the manager via the static access member. since it's static, you don't need
// to create an instance of the class here to call the function
AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
if (!mgr) {
    return 0x0;
}
// get the input event handler, again via a static method. 
// this handler is part of the managing system and feeds events
// to your task
if (!mgr->GetInputEventHandler()) {
    return 0x0;
}

// by default, a file is open for writing. here, we get the filename
TString fileName = AliAnalysisManager::GetCommonFileName();
fileName += ":NanoFwdJpsiCentMult";     // create a subfolder in the file
fileName += suffix;      	 // specific subfolder name for each subwagon
// now we create an instance of your task
AliAnalysisTaskNanoFwdJpsiCentMult* task = new AliAnalysisTaskNanoFwdJpsiCentMult(name.Data());
if(!task) return 0x0;

// add your task to the manager
mgr->AddTask(task);
// your task needs input: here we connect the manager to your task
mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
// same for the output
TString fAllTracksTreeName = "fAllTracksTree";
fAllTracksTreeName += suffix;
TString fOutputListName = "fOutputList";
fOutputListName += suffix;
TString fTrgTreeName = "fTrgTree";
fTrgTreeName += suffix;

mgr->ConnectOutput(task,1,mgr->CreateContainer(fAllTracksTreeName.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));
mgr->ConnectOutput(task,2,mgr->CreateContainer(fOutputListName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
mgr->ConnectOutput(task,3,mgr->CreateContainer(fTrgTreeName.Data(), TTree::Class(), AliAnalysisManager::kOutputContainer,fileName.Data()));

// in the end, this macro returns a pointer to your task. this will be convenient later on
// when you will run your analysis in an analysis train on grid
return task;
}
