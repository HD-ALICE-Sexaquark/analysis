#ifndef AliAnalysisQuickTask_H
#define AliAnalysisQuickTask_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include "TArray.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TH1.h"
#include "TH1F.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMultiplicity.h"
#include "AliPIDResponse.h"

#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "AliVEvent.h"
#include "AliVVertex.h"

class AliPIDResponse;

class AliAnalysisQuickTask : public AliAnalysisTaskSE {
   public:
    AliAnalysisQuickTask();
    AliAnalysisQuickTask(const char* name, Bool_t IsMC);
    virtual ~AliAnalysisQuickTask();

    virtual void UserCreateOutputObjects();
    void PrepareHistograms();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option) { return; }

    /* MC Generated */
    void ProcessMCEvent();

    /* Tracks */
    void ProcessEvent();
    void ProcessTracks();

   private:
    /* Input options */
    Bool_t fIsMC;  //

    /* AliRoot objects */
    AliMCEvent* fMC;                // MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent* fESD;              // reconstructed event
    AliESDInputHandler* fInputHandler;
    AliPIDResponse* fPIDResponse;  // pid response object
    AliESDVertex* fPrimaryVertex;  // primary vertex

    /* ROOT objects */
    TList* fOutputListOfHists;

    /** Histograms **/
    /* MC Gen. Event */
    TH1F* hMCGen_VertexZ;
    /* Event */
    TH1I* hEvents_Bookkeep;
    TH1F* hVertexZ;
    TH1F* hNTracks;
    TH1F* hNSelectedTracks;
    /* Tracks */
    TH1F* hTracks_DCAxy;
    TH1F* hTracks_DCAz;
    TH1F* hTracks_Pt;
    TH2F* hTracks_PhiVsEta;
    /* SPD */
    TH2F* hSPD_MultiplicityVsVertexZ;
    TH1I* hSPD_NTracklets;
    TH2F* hSPD_NTrackletsClu0VsNTracklets;
    TH2F* hSPD_NTrackletsClu1VsNTracklets;
    TH2F* hSPD_PhiVsEta;
    /* ITS */
    TH2F* hITS_LayerNoVsPhi;
    TH2F* hITS_LayerNoVsEta;
    /* TPC */
    TH1F* hTPC_NClusters;
    TH2F* hTPC_NSigmaProtonVsInnerParamP;
    TH2F* hTPC_NSigmaKaonVsInnerParamP;
    TH2F* hTPC_NSigmaPionVsInnerParamP;
    TH2F* hTPC_NSigmaProtonVsEta;
    TH2F* hTPC_NSigmaKaonVsEta;
    TH2F* hTPC_NSigmaPionVsEta;

    AliAnalysisQuickTask(const AliAnalysisQuickTask&);             // not implemented
    AliAnalysisQuickTask& operator=(const AliAnalysisQuickTask&);  // not implemented

    ClassDef(AliAnalysisQuickTask, 1);
};

#endif
