#ifndef AliQuickTaskRadiusWeights_H
#define AliQuickTaskRadiusWeights_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include "TChain.h"
#include "TH1.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"

#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliVVertex.h"

class AliQuickTaskRadiusWeights : public AliAnalysisTaskSE {
   public:
    AliQuickTaskRadiusWeights();
    AliQuickTaskRadiusWeights(const char* name);
    virtual ~AliQuickTaskRadiusWeights();
    virtual void Terminate(Option_t* option) { return; }

    /* Main ~ executed at runtime */
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);

    /* Histograms */
    void PrepareHistograms();

    /* MC Generated */
    void ProcessMCGen();

   private:
    /* Input options ~ persistent */

    /* AliRoot objects */
    AliMCEvent* fMC;                // MC event
    AliVVertex* fMC_PrimaryVertex;  // MC gen. (or true) primary vertex

    /* ROOT objects */
    TList* fOutputListOfHists;  //!

    /** Topography Histograms **/

    TH1F* fHist_MCGen_SFM_Radius;  //!
    TH2F* fHist_MCGen_SFM_YvsX;    //!

    AliQuickTaskRadiusWeights(const AliQuickTaskRadiusWeights&);             // not implemented
    AliQuickTaskRadiusWeights& operator=(const AliQuickTaskRadiusWeights&);  // not implemented

    ClassDef(AliQuickTaskRadiusWeights, 1);
};

#endif
