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
    AliAnalysisQuickTask(const char *name, Bool_t IsMC);
    virtual ~AliAnalysisQuickTask();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option) { return; }

    /* MC Generated */
    void ProcessMCGen();

    /* Tracks */
    void ProcessTracks();

   private:
    /* Input options */
    Bool_t fIsMC;  //

    /* AliRoot objects */
    AliMCEvent *fMC;                // MC event
    AliVVertex *fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent *fESD;              // reconstructed event
    // AliESDInputHandler *fInputHandler;
    AliInputEventHandler *fInputHandler;
    AliPIDResponse *fPIDResponse;  // pid response object
    AliESDVertex *fPrimaryVertex;  // primary vertex

    /* ROOT objects */
    TList *fOutputListOfHists;

    /* Histograms */

    TH1I *hTracklets;
    TH2F *hSPDphivsSPDeta;
    TH1F *hVertexZ;
    TH2F *hClu0VsTracklet;
    TH2F *hClu1VsTracklet;
    TH1I *hEventsProcessed;
    TH2F *hSPDmultiplicityVsVertexZ;

    TH1F *fHistNEvents;

    TH1F *hDCAxy;
    TH1F *hDCAz;

    TH1F *hNTracks;
    TH1F *hPt;
    TH1F *hNTracksPerSelectedEvent;

    TH2F *hITSLayerVsPhi;
    TH2F *hITSLayerVsEta;

    TH2F *hTPCNSigmaProtonVsInnerParamP;
    TH2F *hTPCNSigmaKaonVsInnerParamP;
    TH2F *hTPCNSigmaPionVsInnerParamP;
    TH2F *hTPCNSigmaProtonVsEta;
    TH2F *hTPCNSigmaKaonVsEta;
    TH2F *hTPCNSigmaPionVsEta;

    TH1F *hTPCoutTracks;
    TH1F *hTPCclusters;

    AliAnalysisQuickTask(const AliAnalysisQuickTask &);             // not implemented
    AliAnalysisQuickTask &operator=(const AliAnalysisQuickTask &);  // not implemented

    ClassDef(AliAnalysisQuickTask, 1);
};

#endif
