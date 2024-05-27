#ifndef AliQuickTaskTPCProperties_H
#define AliQuickTaskTPCProperties_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <set>
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
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliExternalTrackParam.h"
#include "AliHelix.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
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
#include "AliVVertex.h"

class AliPIDResponse;

class AliQuickTaskTPCProperties : public AliAnalysisTaskSE {
   public:
    AliQuickTaskTPCProperties();
    AliQuickTaskTPCProperties(const char* name, Bool_t IsMC);
    virtual ~AliQuickTaskTPCProperties();
    virtual void Terminate(Option_t* option) { return; }

    /* Initialization */
    virtual void UserCreateOutputObjects();
    void PrepareTracksHistograms();

    /* Main */
    virtual void UserExec(Option_t* option);

    /* Cuts */
    void DefineTracksCuts(TString cuts_option);

    /* MC Generated */
    void ProcessMCGen();

    /* Tracks */
    Bool_t PassesEventSelection();
    void ProcessTracks();
    Bool_t PassesTrackSelection(AliESDtrack* track);
    void PlotStatus(AliESDtrack* track, TString set, Int_t esdPdgCode);
    void GetTracksHits(AliESDtrack* track, Int_t& firstHit, Int_t& lastHit);

   private:
    /* Input options */
    Bool_t fIsMC;  //

    /* AliRoot objects */
    AliMCEvent* fMC;                // MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex
    AliESDEvent* fESD;              // reconstructed event
    AliESDVertex* fPrimaryVertex;   // primary vertex
    AliPIDResponse* fPIDResponse;   // pid response object
    Double_t fMagneticField;        // magnetic field

    /* ROOT objects */
    TDatabasePDG fPDG;
    TList* fOutputListOfHists;

    /** Tracks Histograms **/

    /*
     key: `pdg code`
     - `0` : MC tracks
     - `1` : secondary MC tracks
     (PENDING)
     -- PENDING: in a sim.set of reactionID='A', signal K+ appear where they shouldn't exist, because of signal particles that were mis-id
     -- Note: should I only, besides if they're signal or not, use their true id? I think so, because it's a bookkeeping of TRUE particles
    */
    std::unordered_map<Int_t, TH1F*> fHist_Tracks_Bookkeep;
    // key: `set, pdg, property`, value: histogram
    std::map<std::tuple<TString, Int_t, TString>, TH1F*> fHist_Tracks;
    // key: `set`, value: histogram
    std::map<TString, TH2F*> f2DHist_Tracks_TPCsignal;

    /*** Cuts ***/

    /* Track selection */
    Float_t kMax_NSigma_Pion;
    Float_t kMax_NSigma_Kaon;
    Float_t kMax_NSigma_Proton;
    Float_t kMax_Track_Eta;
    Float_t kMin_Track_NTPCClusters;
    Float_t kMax_Track_Chi2PerNTPCClusters;
    Bool_t kTurnedOn_Track_StatusCuts;
    Bool_t kTurnedOn_Track_RejectKinks;
    Float_t kMin_Track_DCAxy_wrtPV;
    Float_t kMin_Track_DCAz_wrtPV;
    std::unordered_map<Int_t, Float_t> kMin_Track_Pt;

    AliQuickTaskTPCProperties(const AliQuickTaskTPCProperties&);             // not implemented
    AliQuickTaskTPCProperties& operator=(const AliQuickTaskTPCProperties&);  // not implemented

    ClassDef(AliQuickTaskTPCProperties, 1);
};

#endif
