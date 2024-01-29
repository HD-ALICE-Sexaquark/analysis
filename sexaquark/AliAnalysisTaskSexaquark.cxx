#include <algorithm>
#include <array>
#include <iostream>
#include <map>
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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

#define HomogeneousField  // homogeneous field in z direction, required by KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFVertex.h"

#include "AliAnalysisTaskSexaquark.h"

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      fIsMC(0),
      fSourceOfV0s(),
      fSimulationSet(""),
      fReactionID(0),
      fInjectedSexaquarkMass(0),
      fReactionChannel(""),
      fStruckNucleonPDG(0),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      kMaxNSigma_Pion(0.),
      kMaxNSigma_Kaon(0.),
      kMaxNSigma_Proton(0.),
      fPDG(),
      fTree(0),
      fEvent(),
      mcIndicesOfTrueV0s(),
      isMcIdxSignal(),
      isMcIdxDaughterOfSecondaryV0(),
      isMcIdxDaughterOfTrueV0(),
      isMcIdxDaughterOfSignalV0(),
      getMcIdxOfTrueV0_fromMcIdxOfDau(),
      getMcIdxOfNegDau_fromMcIdxOfTrueV0(),
      getMcIdxOfPosDau_fromMcIdxOfTrueV0() {}

/*
 Constructor, called locally.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name, Bool_t IsMC, TString SourceOfV0s, TString SimulationSet)
    : AliAnalysisTaskSE(name),
      fIsMC(IsMC),
      fSourceOfV0s(SourceOfV0s),
      fSimulationSet(SimulationSet),
      fReactionID(0),
      fInjectedSexaquarkMass(0),
      fReactionChannel(""),
      fStruckNucleonPDG(0),
      fOutputListOfTrees(0),
      fOutputListOfHists(0),
      fMC(0),
      fMC_PrimaryVertex(0),
      fESD(0),
      fPIDResponse(0),
      fPrimaryVertex(0),
      fMagneticField(0.),
      kMaxNSigma_Pion(0.),
      kMaxNSigma_Kaon(0.),
      kMaxNSigma_Proton(0.),
      fPDG(),
      fTree(0),
      fEvent(),
      mcIndicesOfTrueV0s(),
      isMcIdxSignal(),
      isMcIdxDaughterOfSecondaryV0(),
      isMcIdxDaughterOfTrueV0(),
      isMcIdxDaughterOfSignalV0(),
      getMcIdxOfTrueV0_fromMcIdxOfDau(),
      getMcIdxOfNegDau_fromMcIdxOfTrueV0(),
      getMcIdxOfPosDau_fromMcIdxOfTrueV0() {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    CheckForInputErrors();
    Initialize();
}

/*
 Destructor.
 */
AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() {
    if (fOutputListOfTrees) {
        delete fOutputListOfTrees;
    }
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {

    // add the PID routine
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) {
        AliFatal("ERROR: AliAnalysisManager couldn't be found.");
    }

    AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
    if (!inputHandler) {
        AliFatal("ERROR: AliESDInputHandler couldn't be found.");
    }

    fPIDResponse = inputHandler->GetPIDResponse();

    /*** Trees ***/

    fOutputListOfTrees = new TList();
    fOutputListOfTrees->SetOwner(kTRUE);

    fTree = new TTree("Events", "Tree of Events");
    SetBranches();
    fOutputListOfTrees->Add(fTree);

    PostData(1, fOutputListOfTrees);

    /*** Histograms ***/

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE);

    /* MC Gen. histograms */

    fHist_True_Signal_SecVertex_Radius = new TH1F("True_Signal_SecVertex_Radius", "", 100, 0., 200.);
    fHist_True_InjBkg_SecVertex_Radius = new TH1F("True_InjBkg_SecVertex_Radius", "", 100, 0., 200.);

    fOutputListOfHists->Add(fHist_True_Signal_SecVertex_Radius);
    fOutputListOfHists->Add(fHist_True_InjBkg_SecVertex_Radius);

    fHist_True_FermiSexaquark_Pt = new TH1F("True_FermiSexaquark_Pt", "", 140, -7., 7.);
    fHist_True_FermiSexaquark_Mass = new TH1F("True_FermiSexaquark_Mass", "", 100, -5., 5.);

    fOutputListOfHists->Add(fHist_True_FermiSexaquark_Pt);
    fOutputListOfHists->Add(fHist_True_FermiSexaquark_Mass);

    fHist_True_Signal_AntiLambda_Pt = new TH1F("True_Signal_AntiLambda_Pt", "", 100, 0., 10.);
    fHist_True_Signal_KaonZeroShort_Pt = new TH1F("True_Signal_KaonZeroShort_Pt", "", 100, 0., 10.);
    fHist_True_Signal_AntiLambda_DCAwrtPV = new TH1F("True_Signal_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_True_Signal_KaonZeroShort_DCAwrtPV = new TH1F("True_Signal_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);
    fHist_True_Signal_AntiLambda_CPAwrtPV = new TH1F("True_Signal_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_True_Signal_KaonZeroShort_CPAwrtPV = new TH1F("True_Signal_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);

    fHist_True_AntiLambda_Pt = new TH1F("True_AntiLambda_Pt", "", 100, 0., 10.);
    fHist_True_KaonZeroShort_Pt = new TH1F("True_KaonZeroShort_Pt", "", 100, 0., 10.);
    fHist_True_AntiLambda_DCAwrtPV = new TH1F("True_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_True_KaonZeroShort_DCAwrtPV = new TH1F("True_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);
    fHist_True_AntiLambda_CPAwrtPV = new TH1F("True_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_True_KaonZeroShort_CPAwrtPV = new TH1F("True_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);

    fOutputListOfHists->Add(fHist_True_AntiLambda_Pt);
    fOutputListOfHists->Add(fHist_True_KaonZeroShort_Pt);
    fOutputListOfHists->Add(fHist_True_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_True_KaonZeroShort_DCAwrtPV);
    fOutputListOfHists->Add(fHist_True_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_True_KaonZeroShort_CPAwrtPV);

    fOutputListOfHists->Add(fHist_True_Signal_AntiLambda_Pt);
    fOutputListOfHists->Add(fHist_True_Signal_KaonZeroShort_Pt);
    fOutputListOfHists->Add(fHist_True_Signal_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_True_Signal_KaonZeroShort_DCAwrtPV);
    fOutputListOfHists->Add(fHist_True_Signal_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_True_Signal_KaonZeroShort_CPAwrtPV);

    fHist_True_AntiLambda_Bookkeep = new TH1F("True_AntiLambda_Bookkeep", "", 10, 0., 10.);
    fHist_True_KaonZeroShort_Bookkeep = new TH1F("True_KaonZeroShort_Bookkeep", "", 10, 0., 10.);
    fHist_True_AntiProton_Bookkeep = new TH1F("True_AntiProton_Bookkeep", "", 10, 0., 10.);
    fHist_True_PosKaon_Bookkeep = new TH1F("True_PosKaon_Bookkeep", "", 10, 0., 10.);
    fHist_True_PiPlus_Bookkeep = new TH1F("True_PiPlus_Bookkeep", "", 10, 0., 10.);
    fHist_True_PiMinus_Bookkeep = new TH1F("True_PiMinus_Bookkeep", "", 10, 0., 10.);

    fOutputListOfHists->Add(fHist_True_AntiLambda_Bookkeep);
    fOutputListOfHists->Add(fHist_True_KaonZeroShort_Bookkeep);
    fOutputListOfHists->Add(fHist_True_AntiProton_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PosKaon_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PiPlus_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PiMinus_Bookkeep);

    fHist_True_AntiNeutron_Pt = new TH1F("True_AntiNeutron_Pt", "", 100, 0., 10.);
    fHist_True_AntiNeutronThatInteracted_Pt = new TH1F("True_AntiNeutronThatInteracted_Pt", "", 100, 0., 10.);
    fHist_True_AntiNeutron_NDaughters = new TH1F("True_AntiNeutron_NDaughters", "", 100, 0., 100.);
    fHist_True_AntiNeutron_PDGDaughters = new TH1F("True_AntiNeutron_PDGDaughters", "", 7000, -3500., 3500.);
    fHist_True_AntiNeutron_Bookkeep = new TH1F("True_AntiNeutron_Bookkeep", "", 10, 0., 10.);

    fOutputListOfHists->Add(fHist_True_AntiNeutron_Pt);
    fOutputListOfHists->Add(fHist_True_AntiNeutronThatInteracted_Pt);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_NDaughters);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_PDGDaughters);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_Bookkeep);

    fHist_True_InjBkgProducts_Bookkeep = new TH1F("True_InjBkgProducts_Bookkeep", "", 50, 0., 50.);
    fHist_True_AllSecondaries_MothersPDG = new TH1F("True_AllSecondaries_MothersPDG", "", 7000, -3500, 3500.);

    fOutputListOfHists->Add(fHist_True_InjBkgProducts_Bookkeep);
    fOutputListOfHists->Add(fHist_True_AllSecondaries_MothersPDG);

    /* rec. tracks histograms (that require MC Gen. information) */

    fHist_Rec_Signal_AntiProton_Pt = new TH1F("Rec_Signal_AntiProton_Pt", "", 100, 0., 10.);
    fHist_Rec_Signal_PosKaon_Pt = new TH1F("Rec_Signal_PosKaon_Pt", "", 100, 0., 10.);
    fHist_Rec_Signal_PiPlus_Pt = new TH1F("Rec_Signal_PiPlus_Pt", "", 100, 0., 10.);
    fHist_Rec_Signal_PiMinus_Pt = new TH1F("Rec_Signal_PiMinus_Pt", "", 100, 0., 10.);

    fOutputListOfHists->Add(fHist_Rec_Signal_AntiProton_Pt);
    fOutputListOfHists->Add(fHist_Rec_Signal_PosKaon_Pt);
    fOutputListOfHists->Add(fHist_Rec_Signal_PiPlus_Pt);
    fOutputListOfHists->Add(fHist_Rec_Signal_PiMinus_Pt);

    /* rec. tracks histograms (that don't require MC Gen. information) */

    fHist_Rec_AntiProton_Pt = new TH1F("Rec_AntiProton_Pt", "", 100, 0., 10.);
    fHist_Rec_PosKaon_Pt = new TH1F("Rec_PosKaon_Pt", "", 100, 0., 10.);
    fHist_Rec_PiPlus_Pt = new TH1F("Rec_PiPlus_Pt", "", 100, 0., 10.);
    fHist_Rec_PiMinus_Pt = new TH1F("Rec_PiMinus_Pt", "", 100, 0., 10.);

    fOutputListOfHists->Add(fHist_Rec_AntiProton_Pt);
    fOutputListOfHists->Add(fHist_Rec_PosKaon_Pt);
    fOutputListOfHists->Add(fHist_Rec_PiPlus_Pt);
    fOutputListOfHists->Add(fHist_Rec_PiMinus_Pt);

    f2DHist_TPC_Signal = new TH2F("TPC_Signal", "", 500, 0., 10., 500, 0., 1000.);

    fOutputListOfHists->Add(f2DHist_TPC_Signal);

    /* V0s histograms (that require MC Gen. information) */

    fHist_Findable_True_AntiLambda_Mass = new TH1F("Findable_True_AntiLambda_Mass", "", 100, 0., 10.);
    fHist_Findable_True_AntiLambda_Radius = new TH1F("Findable_True_AntiLambda_Radius", "", 100, 0., 200.);
    fHist_Findable_True_AntiLambda_CPAwrtPV = new TH1F("Findable_True_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_Findable_True_AntiLambda_DCAwrtPV = new TH1F("Findable_True_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_Findable_True_KaonZeroShort_Mass = new TH1F("Findable_True_KaonZeroShort_Mass", "", 100, 0., 10.);
    fHist_Findable_True_KaonZeroShort_Radius = new TH1F("Findable_True_KaonZeroShort_Radius", "", 100, 0., 200.);
    fHist_Findable_True_KaonZeroShort_CPAwrtPV = new TH1F("Findable_True_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);
    fHist_Findable_True_KaonZeroShort_DCAwrtPV = new TH1F("Findable_True_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);

    fHist_Findable_Signal_AntiLambda_Mass = new TH1F("Findable_Signal_AntiLambda_Mass", "", 100, 0., 10.);
    fHist_Findable_Signal_AntiLambda_Radius = new TH1F("Findable_Signal_AntiLambda_Radius", "", 100, 0., 200.);
    fHist_Findable_Signal_AntiLambda_CPAwrtPV = new TH1F("Findable_Signal_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_Findable_Signal_AntiLambda_DCAwrtPV = new TH1F("Findable_Signal_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_Findable_Signal_KaonZeroShort_Mass = new TH1F("Findable_Signal_KaonZeroShort_Mass", "", 100, 0., 10.);
    fHist_Findable_Signal_KaonZeroShort_Radius = new TH1F("Findable_Signal_KaonZeroShort_Radius", "", 100, 0., 200.);
    fHist_Findable_Signal_KaonZeroShort_CPAwrtPV = new TH1F("Findable_Signal_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);
    fHist_Findable_Signal_KaonZeroShort_DCAwrtPV = new TH1F("Findable_Signal_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);

    fHist_Found_True_AntiLambda_Mass = new TH1F("Found_True_AntiLambda_Mass", "", 100, 0., 10.);
    fHist_Found_True_AntiLambda_Radius = new TH1F("Found_True_AntiLambda_Radius", "", 100, 0., 200.);
    fHist_Found_True_AntiLambda_CPAwrtPV = new TH1F("Found_True_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_True_AntiLambda_DCAwrtPV = new TH1F("Found_True_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_Found_True_KaonZeroShort_Mass = new TH1F("Found_True_KaonZeroShort_Mass", "", 100, 0., 10.);
    fHist_Found_True_KaonZeroShort_Radius = new TH1F("Found_True_KaonZeroShort_Radius", "", 100, 0., 200.);
    fHist_Found_True_KaonZeroShort_CPAwrtPV = new TH1F("Found_True_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_True_KaonZeroShort_DCAwrtPV = new TH1F("Found_True_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);

    fHist_Found_Signal_AntiLambda_Mass = new TH1F("Found_Signal_AntiLambda_Mass", "", 100, 0., 10.);
    fHist_Found_Signal_AntiLambda_Radius = new TH1F("Found_Signal_AntiLambda_Radius", "", 100, 0., 200.);
    fHist_Found_Signal_AntiLambda_CPAwrtPV = new TH1F("Found_Signal_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_Signal_AntiLambda_DCAwrtPV = new TH1F("Found_Signal_AntiLambda_DCAwrtPV", "", 100, 0., 10.);
    fHist_Found_Signal_KaonZeroShort_Mass = new TH1F("Found_Signal_KaonZeroShort_Mass", "", 100, 0., 10.);
    fHist_Found_Signal_KaonZeroShort_Radius = new TH1F("Found_Signal_KaonZeroShort_Radius", "", 100, 0., 200.);
    fHist_Found_Signal_KaonZeroShort_CPAwrtPV = new TH1F("Found_Signal_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_Signal_KaonZeroShort_DCAwrtPV = new TH1F("Found_Signal_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);

    fOutputListOfHists->Add(fHist_Findable_True_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_Findable_True_AntiLambda_Radius);
    fOutputListOfHists->Add(fHist_Findable_True_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_True_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_True_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_Findable_True_KaonZeroShort_Radius);
    fOutputListOfHists->Add(fHist_Findable_True_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_True_KaonZeroShort_DCAwrtPV);

    fOutputListOfHists->Add(fHist_Findable_Signal_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_Findable_Signal_AntiLambda_Radius);
    fOutputListOfHists->Add(fHist_Findable_Signal_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_Signal_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_Signal_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_Findable_Signal_KaonZeroShort_Radius);
    fOutputListOfHists->Add(fHist_Findable_Signal_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_Signal_KaonZeroShort_DCAwrtPV);

    fOutputListOfHists->Add(fHist_Found_True_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_Found_True_AntiLambda_Radius);
    fOutputListOfHists->Add(fHist_Found_True_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_True_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Found_True_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_Found_True_KaonZeroShort_Radius);
    fOutputListOfHists->Add(fHist_Found_True_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_True_KaonZeroShort_DCAwrtPV);

    fOutputListOfHists->Add(fHist_Found_Signal_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiLambda_Radius);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Found_Signal_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_Found_Signal_KaonZeroShort_Radius);
    fOutputListOfHists->Add(fHist_Found_Signal_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_Signal_KaonZeroShort_DCAwrtPV);

    /* V0s histograms (that don't require MC Gen. information) */

    fHist_Found_AntiLambda_Mass = new TH1F("Found_AntiLambda_Mass", "", 100, 0., 10.);
    fHist_Found_AntiLambda_Radius = new TH1F("Found_AntiLambda_Radius", "", 100, 0., 200.);
    fHist_Found_AntiLambda_CPAwrtPV = new TH1F("Found_AntiLambda_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_AntiLambda_DCAwrtPV = new TH1F("Found_AntiLambda_DCAwrtPV", "", 100, 0., 10.);

    fHist_Found_KaonZeroShort_Mass = new TH1F("Found_KaonZeroShort_Mass", "", 100, 0., 10.);
    fHist_Found_KaonZeroShort_Radius = new TH1F("Found_KaonZeroShort_Radius", "", 100, 0., 200.);
    fHist_Found_KaonZeroShort_CPAwrtPV = new TH1F("Found_KaonZeroShort_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_KaonZeroShort_DCAwrtPV = new TH1F("Found_KaonZeroShort_DCAwrtPV", "", 100, 0., 10.);

    fOutputListOfHists->Add(fHist_Found_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_Found_AntiLambda_Radius);
    fOutputListOfHists->Add(fHist_Found_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_AntiLambda_DCAwrtPV);

    fOutputListOfHists->Add(fHist_Found_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_Found_KaonZeroShort_Radius);
    fOutputListOfHists->Add(fHist_Found_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_KaonZeroShort_DCAwrtPV);

    /* anti-sexaquark histograms (that require MC Gen. information) */

    fHist_Findable_AntiSexaquark_Mass = new TH1F("Findable_AntiSexaquark_Mass", "", 150, -5., 10.);
    fHist_Findable_AntiSexaquark_Pt = new TH1F("Findable_AntiSexaquark_Pt", "", 100, -2.5, 7.5);
    fHist_Findable_AntiSexaquark_DCAwrtPV = new TH1F("Findable_AntiSexaquark_DCAwrtPV", "", 100, 0., 10.);
    fHist_Findable_AntiSexaquark_CPAwrtPV = new TH1F("Findable_AntiSexaquark_CPAwrtPV", "", 100, -1., 1.);
    fHist_Findable_AntiSexaquark_Radius = new TH1F("Findable_AntiSexaquark_Radius", "", 100, 0., 200.);

    fHist_Found_Signal_AntiSexaquark_Mass = new TH1F("Found_Signal_AntiSexaquark_Mass", "", 150, -5., 10.);
    fHist_Found_Signal_AntiSexaquark_Pt = new TH1F("Found_Signal_AntiSexaquark_Pt", "", 100, -2.5, 7.5);
    fHist_Found_Signal_AntiSexaquark_DCAbtwDau = new TH1F("Found_Signal_AntiSexaquark_DCAbtwDau", "", 100, 0., 10.);
    fHist_Found_Signal_AntiSexaquark_DCAwrtPV = new TH1F("Found_Signal_AntiSexaquark_DCAwrtPV", "", 100, 0., 10.);
    fHist_Found_Signal_AntiSexaquark_CPAwrtPV = new TH1F("Found_Signal_AntiSexaquark_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_Signal_AntiSexaquark_Radius = new TH1F("Found_Signal_AntiSexaquark_Radius", "", 100, 0., 200.);

    fOutputListOfHists->Add(fHist_Findable_AntiSexaquark_Mass);
    fOutputListOfHists->Add(fHist_Findable_AntiSexaquark_Pt);
    fOutputListOfHists->Add(fHist_Findable_AntiSexaquark_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_AntiSexaquark_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Findable_AntiSexaquark_Radius);

    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_Mass);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_Pt);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_DCAbtwDau);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_Signal_AntiSexaquark_Radius);

    /* anti-sexaquark histograms (that don't require MC Gen. information) */

    fHist_Found_AntiSexaquark_Mass = new TH1F("Found_AntiSexaquark_Mass", "", 150, -5., 10.);
    fHist_Found_AntiSexaquark_Pt = new TH1F("Found_AntiSexaquark_Pt", "", 100, -2.5, 7.5);
    fHist_Found_AntiSexaquark_DCAbtwDau = new TH1F("Found_AntiSexaquark_DCAbtwDau", "", 100, 0., 10.);
    fHist_Found_AntiSexaquark_DCAwrtPV = new TH1F("Found_AntiSexaquark_DCAwrtPV", "", 100, 0., 10.);
    fHist_Found_AntiSexaquark_CPAwrtPV = new TH1F("Found_AntiSexaquark_CPAwrtPV", "", 100, -1., 1.);
    fHist_Found_AntiSexaquark_Radius = new TH1F("Found_AntiSexaquark_Radius", "", 100, 0., 200.);

    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_Mass);
    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_Pt);
    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_DCAbtwDau);
    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_DCAwrtPV);
    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_CPAwrtPV);
    fOutputListOfHists->Add(fHist_Found_AntiSexaquark_Radius);

    PostData(2, fOutputListOfHists);
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
 - Uses: `fIsMC`, `fMC_PrimaryVertex`, `fESD`, `fMagneticField`, `fPrimaryVertex`, `fSourceOfV0s`, `fReactionChannel`, `fOutputListOfTrees`,
 `fOutputListOfHists`
*/
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

    // load MC generated event
    fMC = MCEvent();

    if (!fMC) {
        AliFatal("ERROR: AliMCEvent couldn't be found.");
    }

    fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());

    // load reconstructed event
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());

    if (!fESD) {
        AliFatal("ERROR: AliESDEvent couldn't be found.");
    }

    // load magnetic field
    fMagneticField = fESD->GetMagneticField();

    // initialize KFparticle
    KFParticle::SetField(fMagneticField);

    // load primary vertex
    // [note: const_cast removes the const qualifier from the pointer]
    fPrimaryVertex = const_cast<AliESDVertex*>(fESD->GetPrimaryVertex());

    if (fIsMC) {
        ProcessMCGen();
        ProcessSignalInteractions();
    }

    ProcessTracks();

    if (fIsMC) {
        ProcessFindableV0s();
        ProcessFindableSexaquarks();
    }

    if (fSourceOfV0s == "offline") {
        OfficialV0Finder(kFALSE);
    }

    if (fSourceOfV0s == "online") {
        OfficialV0Finder(kTRUE);
    }

    if (fSourceOfV0s == "custom") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            CustomV0Finder(-3122, -2212, 211);
            CustomV0Finder(310, -211, 211);
        }
    }

    if (fSourceOfV0s == "offline" || fSourceOfV0s == "online" || fSourceOfV0s == "custom") {
        if (fReactionID == 'A') {
            SexaquarkFinder_ChannelA_Geo();
        }
    }

    if (fSourceOfV0s == "kalman") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            ReconstructV0s_KF(-3122, -2212, 211);
            ReconstructV0s_KF(310, -211, 211);
            SexaquarkFinder_ChannelA_KF();
        }
    }

    /*
        if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+") {
            ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
        }
        if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+,pi-,pi+") {
            ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
            ReconstructV0s_KF(idxPiMinusTracks, idxPosKaonTracks, 0, -211, 321, kfPionKaonPairs, idxPionKaonPairs);
            ReconstructV0s_KF(idxPiMinusTracks, idxPiPlusTracks, 0, -211, 211, kfPionPairs, idxPionPairs);
        }
        if (fReactionChannel == "AntiSexaquark,P->AntiProton,K+,K+,pi0") {
            ReconstructV0s_KF(idxPosKaonTracks, idxPosKaonTracks, 0, 312, 321, kfPosKaonPairs, idxPosKaonPairs);
        }
    */

    // FillTree();

    ClearContainers();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfTrees);
    PostData(2, fOutputListOfHists);
}

/*
 Initialize analysis task.
*/
void AliAnalysisTaskSexaquark::Initialize() {

    /* Parse TString: get reaction ID and injected sexaquark mass from SimulationSet */

    fReactionID = fSimulationSet(0);
    fInjectedSexaquarkMass = ((TString)fSimulationSet(1, fSimulationSet.Length())).Atof();

    /* PENDING: inj. sexaquark mass doesn't make sense when analyzing data... do smth about it... */

    /* Assign the rest of properties related to the respective Reaction Channel */

    if (fReactionID == 'A') {
        fReactionChannel = "AntiSexaquark,N->AntiLambda,K0S";
        fProductsPDG = {-3122, 310};
        fFinalStateProductsPDG = {-2212, 211, -211, 211};
        fStruckNucleonPDG = 2112;
    } else if (fReactionID == 'D') {
        fReactionChannel = "AntiSexaquark,P->AntiLambda,K+";
        fProductsPDG = {-3122, 321};
        fFinalStateProductsPDG = {-2212, 211, 321};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'E') {
        fReactionChannel = "AntiSexaquark,P->AntiLambda,K+,pi-,pi+";
        fProductsPDG = {-3122, 321, -211, 211};
        fFinalStateProductsPDG = {-2212, 211, 321, -211, 211};
        fStruckNucleonPDG = 2212;
    } else if (fReactionID == 'H') {
        fReactionChannel = "AntiSexaquark,P->AntiProton,K+,K+,pi0";
        fProductsPDG = {-2212, 321, 321, 111};
        fFinalStateProductsPDG = {-2212, 321, 321, 11, -11, 11, -11, 11, -11, 11, -11};  // not used
        fStruckNucleonPDG = 2212;
    }

    DefineCuts("Standard");

    AliInfo("Initializing AliAnalysisTaskSexaquark...");
    AliInfo("INPUT OPTIONS:");
    AliInfoF(">> IsMC                = %i", (Int_t)fIsMC);
    AliInfoF(">> Source of V0s       = %s", fSourceOfV0s.Data());
    AliInfoF(">> Simulation Set      = %s", fSimulationSet.Data());
    AliInfoF(">> Reaction ID         = %c", fReactionID);
    AliInfoF(">> Reaction Channel    = %s", fReactionChannel.Data());
    AliInfoF(">> Inj. Sexaquark Mass = %.2f", fInjectedSexaquarkMass);
}

/*
 - Uses: `fReactionChannel`
*/
void AliAnalysisTaskSexaquark::DefineCuts(TString cuts_option) {
    kMaxNSigma_Pion = 3.;
    kMaxNSigma_Kaon = 3.;
    kMaxNSigma_Proton = 3.;
    if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S" && cuts_option == "Standard") {
    }
    if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+" && cuts_option == "Standard") {
    }
    if (fReactionChannel == "AntiSexaquark,P->AntiLambda,K+,pi-,pi+" && cuts_option == "Standard") {
    }
    if (fReactionChannel == "AntiSexaquark,P->AntiProton,K+,K+,pi0" && cuts_option == "Standard") {
    }
}

/*
 Search for possible contradictions in the input options
*/
void AliAnalysisTaskSexaquark::CheckForInputErrors() {
    std::array<std::string, 4> validSources = {"offline", "online", "custom", "kalman"};
    if (std::find(validSources.begin(), validSources.end(), fSourceOfV0s) == validSources.end()) {
        AliFatal("ERROR: source of V0s must be \"offline\", \"online\", \"custom\" or \"kalman\".");
    }
    if (!fIsMC) {
        AliFatal("ERROR: data cannot have access to true MC information.");
    }
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

/*
 Loop over MC particles in a single event. Store the indices of the signal particles.
 - Uses: `fMC`, `fPDG`, `fMC_PrimaryVertex`
*/
void AliAnalysisTaskSexaquark::ProcessMCGen() {

    AliMCParticle* mcPart;
    Int_t pdg_mc;

    AliMCParticle* mcMother;
    Int_t pdg_mother;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    Float_t antineutron_radius;
    Int_t n_daughters = 0;

    Int_t mcIdxNegDaughter = 0;
    Int_t mcIdxPosDaughter = 0;

    TVector3 originVertex;
    TVector3 decayVertex;

    Float_t pt;
    Float_t radius;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;
    Bool_t is_secondary;
    Bool_t is_signal;

    std::vector<Int_t> relevantPdg = {-3122, 310, -2212, 321, 211, -211, -2112, 111};

    /*** Loop over MC gen. particles in a single event ***/

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);
        pdg_mc = mcPart->PdgCode();
        n_daughters = mcPart->GetNDaughters();

        mcIdxNegDaughter = 0;
        mcIdxPosDaughter = 0;
        decayVertex.SetXYZ(0., 0., 0.);

        /* Check if the particle is relevant for the analysis */

        if (!std::any_of(relevantPdg.begin(), relevantPdg.end(), [pdg_mc](Int_t pdg_rel) { return pdg_rel == pdg_mc; })) {
            continue;
        }

        if (mcPart->P() < 1E-6) {
            continue;
        }

        /* Get particle's properties */

        pt = mcPart->Pt();
        originVertex.SetXYZ(mcPart->Xv(), mcPart->Yv(), mcPart->Zv());
        radius = originVertex.Perp();

        if (n_daughters) {
            GetDaughtersInfo(mcPart, mcIdxNegDaughter, mcIdxPosDaughter, decayVertex);
            cpa_wrt_pv = CosinePointingAngle(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                             fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), decayVertex.X(), decayVertex.Y(), decayVertex.Z(),
                                      fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        }

        // If they're signal, their MC status number should be among 600 and 619.
        // This is done exclusively for the sexaquark MC productions.
        // This was done to recognize each interaction, assigning a different number to each of
        // the 20 antisexaquark-nucleon reactions that were injected in a single event .
        is_signal = mcPart->MCStatusCode() >= 600 && mcPart->MCStatusCode() < 620;

        // Include signal particles, as well, because they were injected as primaries, and we intend that they are not
        is_secondary = mcPart->IsSecondaryFromMaterial() || mcPart->IsSecondaryFromWeakDecay() || is_signal;

        /** Fill histograms, indices vector and maps **/

        /* anti-lambdas */

        if (pdg_mc == -3122) {

            fHist_True_AntiLambda_Bookkeep->Fill(0);
            fHist_True_AntiLambda_Pt->Fill(pt);

            if (n_daughters) {
                fHist_True_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_True_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
            }

            if (mcIdxNegDaughter && mcIdxPosDaughter) {
                mcIndicesOfTrueV0s.push_back(mcIdx);
                isMcIdxDaughterOfTrueV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfTrueV0[mcIdxPosDaughter] = kTRUE;
                getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdxNegDaughter] = mcIdx;
                getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdxPosDaughter] = mcIdx;
                getMcIdxOfNegDau_fromMcIdxOfTrueV0[mcIdx] = mcIdxNegDaughter;
                getMcIdxOfPosDau_fromMcIdxOfTrueV0[mcIdx] = mcIdxPosDaughter;
            }

            if (is_secondary) {
                fHist_True_AntiLambda_Bookkeep->Fill(1);
            }

            // if (is_secondary && n_daughters) {
            //     // PENDING
            // }

            if (is_secondary && mcIdxNegDaughter && mcIdxPosDaughter) {
                isMcIdxDaughterOfSecondaryV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfSecondaryV0[mcIdxPosDaughter] = kTRUE;
            }

            if (is_signal) {
                fHist_True_AntiLambda_Bookkeep->Fill(2);
                fHist_True_Signal_AntiLambda_Pt->Fill(pt);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }

            if (is_signal && n_daughters) {
                fHist_True_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_True_Signal_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
            }

            if (is_signal && mcIdxNegDaughter && mcIdxPosDaughter) {
                isMcIdxSignal[mcIdxNegDaughter] = kTRUE;
                isMcIdxSignal[mcIdxPosDaughter] = kTRUE;
                isMcIdxDaughterOfSignalV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfSignalV0[mcIdxPosDaughter] = kTRUE;
            }
        }

        /* k0s */

        if (pdg_mc == 310) {

            fHist_True_KaonZeroShort_Bookkeep->Fill(0);
            fHist_True_KaonZeroShort_Pt->Fill(pt);

            if (n_daughters) {
                fHist_True_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_True_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
            }

            if (mcIdxNegDaughter && mcIdxPosDaughter) {
                mcIndicesOfTrueV0s.push_back(mcIdx);
                isMcIdxDaughterOfTrueV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfTrueV0[mcIdxPosDaughter] = kTRUE;
                getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdxNegDaughter] = mcIdx;
                getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdxPosDaughter] = mcIdx;
                getMcIdxOfNegDau_fromMcIdxOfTrueV0[mcIdx] = mcIdxNegDaughter;
                getMcIdxOfPosDau_fromMcIdxOfTrueV0[mcIdx] = mcIdxPosDaughter;
            }

            if (is_secondary) {
                fHist_True_KaonZeroShort_Bookkeep->Fill(1);
            }

            // if (is_secondary && n_daughters) {
            //     // PENDING
            // }

            if (is_secondary && mcIdxNegDaughter && mcIdxPosDaughter) {
                isMcIdxDaughterOfSecondaryV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfSecondaryV0[mcIdxPosDaughter] = kTRUE;
            }

            if (is_signal) {
                fHist_True_KaonZeroShort_Bookkeep->Fill(2);
                fHist_True_Signal_KaonZeroShort_Pt->Fill(pt);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }

            if (is_signal && n_daughters) {
                fHist_True_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_True_Signal_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
            }

            if (is_signal && mcIdxNegDaughter && mcIdxPosDaughter) {
                isMcIdxSignal[mcIdxNegDaughter] = kTRUE;
                isMcIdxSignal[mcIdxPosDaughter] = kTRUE;
                isMcIdxDaughterOfSignalV0[mcIdxNegDaughter] = kTRUE;
                isMcIdxDaughterOfSignalV0[mcIdxPosDaughter] = kTRUE;
            }
        }

        /* anti-protons */

        if (pdg_mc == -2212) {
            fHist_True_AntiProton_Bookkeep->Fill(0);
            if (is_secondary) fHist_True_AntiProton_Bookkeep->Fill(1);
            if (is_signal) {
                fHist_True_AntiProton_Bookkeep->Fill(2);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }
        }

        /* positive kaons */

        if (pdg_mc == 321) {
            fHist_True_PosKaon_Bookkeep->Fill(0);
            if (is_secondary) fHist_True_PosKaon_Bookkeep->Fill(1);
            if (is_signal) {
                fHist_True_PosKaon_Bookkeep->Fill(2);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }
        }

        /* positive pions */

        if (pdg_mc == 211) {
            fHist_True_PiPlus_Bookkeep->Fill(0);
            if (is_secondary) fHist_True_PiPlus_Bookkeep->Fill(1);
            if (is_signal) {
                fHist_True_PiPlus_Bookkeep->Fill(2);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }
        }

        /* negative pions */

        if (pdg_mc == -211) {
            fHist_True_PiMinus_Bookkeep->Fill(0);
            if (is_secondary) fHist_True_PiMinus_Bookkeep->Fill(1);
            if (is_signal) {
                fHist_True_PiMinus_Bookkeep->Fill(2);
                isMcIdxSignal[mcIdx] = kTRUE;
                getReactionIdx_fromMcIdx[mcIdx] = mcPart->MCStatusCode();
                getMcIdx_fromReactionIdx[mcPart->MCStatusCode()].push_back(mcIdx);
            }
        }

        /* anti-neutrons */
        // PENDING

        if (pdg_mc == -2112) {

            fHist_True_AntiNeutron_Bookkeep->Fill(0);
            fHist_True_AntiNeutron_Pt->Fill(pt);

            if (n_daughters) {

                fHist_True_AntiNeutron_Bookkeep->Fill(1);
                fHist_True_AntiNeutron_NDaughters->Fill(mcPart->GetNDaughters());
                fHist_True_AntiNeutronThatInteracted_Pt->Fill(pt);

                /* Fill, then, the PDG code of its daughters */

                for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

                    mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
                    pdg_dau = mcDaughter->PdgCode();

                    if (TMath::Abs(pdg_dau) < 3500) {
                        fHist_True_AntiNeutron_PDGDaughters->Fill(pdg_dau);
                    } else {
                        fHist_True_AntiNeutron_PDGDaughters->Fill(0);
                    }

                    if (mcIdxDaughter == mcPart->GetDaughterFirst()) {
                        antineutron_radius = TMath::Sqrt(TMath::Power(mcDaughter->Xv(), 2) + TMath::Power(mcDaughter->Yv(), 2));
                        fHist_True_InjBkg_SecVertex_Radius->Fill(antineutron_radius);
                    }
                }
            }
        }
    }  // end of loop over MC particles
}

/*
 Loop over the daughters of an MC particle.
 - Input: `mcPart`
 - Output: `mcIdxNegDaughter`, `mcIdxPosDaughter`, `decayVertex`
*/
void AliAnalysisTaskSexaquark::GetDaughtersInfo(AliMCParticle* mcPart, Int_t& mcIdxNegDaughter, Int_t& mcIdxPosDaughter, TVector3& decayVertex) {

    Int_t pdg_mc = mcPart->PdgCode();

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {

        mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
        pdg_dau = mcDaughter->PdgCode();

        if ((pdg_mc == -3122 && pdg_dau == -2212) || (pdg_mc == 310 && pdg_dau == -211)) {
            mcIdxNegDaughter = mcIdxDaughter;
        } else if ((pdg_mc == -3122 || pdg_mc == 310) && pdg_dau == 211) {
            mcIdxPosDaughter = mcIdxDaughter;
        }

        /* In any case, get the decay vertex of the V0 */

        if (decayVertex.Mag() == 0) decayVertex.SetXYZ(mcDaughter->Xv(), mcDaughter->Yv(), mcDaughter->Zv());
    }
}

/*
 Loop over antisexaquark-nucleon interactions.
 - Uses: `fMC`, `fPDG`, `fStruckNucleonPDG`
 - Containers: `getMcIdx_fromReactionIdx`
*/
void AliAnalysisTaskSexaquark::ProcessSignalInteractions() {

    AliMCParticle* mcProduct;
    Int_t mcIdxProduct;

    /* Declare TLorentzVectors */

    TLorentzVector lvProduct;

    TLorentzVector lvStruckNucleon(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

    TLorentzVector lvAntiSexaquark;

    /* Loop over interactions */

    for (Int_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Reset anti-sexaquark TLorentzVector */

        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        /* Loop over reaction products */

        for (Int_t productIdx = 0; productIdx < (Int_t)getMcIdx_fromReactionIdx[reactionIdx].size(); productIdx++) {

            mcIdxProduct = getMcIdx_fromReactionIdx[reactionIdx][productIdx];
            mcProduct = (AliMCParticle*)fMC->GetTrack(mcIdxProduct);

            lvProduct.SetPxPyPzE(mcProduct->Px(), mcProduct->Py(), mcProduct->Pz(), mcProduct->E());

            lvAntiSexaquark = lvAntiSexaquark + lvProduct;
        }

        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /* Fill histograms */

        fHist_True_FermiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
        fHist_True_FermiSexaquark_Mass->Fill(lvAntiSexaquark.M());
        fHist_True_Signal_SecVertex_Radius->Fill(TMath::Sqrt(TMath::Power(mcProduct->Xv(), 2) + TMath::Power(mcProduct->Yv(), 2)));
    }
}

/*                */
/**    Tracks    **/
/*** ========== ***/

/*
 Loop over the reconstructed tracks in a single event.
*/
void AliAnalysisTaskSexaquark::ProcessTracks() {

    AliESDtrack* track;

    Int_t mcIdx;
    Int_t mcIdxOfTrueV0;

    Float_t n_sigma_pion;
    Float_t n_sigma_kaon;
    Float_t n_sigma_proton;

    Float_t pt;

    /* Loop over tracks in a single event */

    for (Int_t esdIdxTrack = 0; esdIdxTrack < fESD->GetNumberOfTracks(); esdIdxTrack++) {

        /* Get track */

        track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxTrack));

        /* Get MC info */

        mcIdx = TMath::Abs(track->GetLabel());

        /* Fill containers */

        getMcIdx_fromEsdIdx[esdIdxTrack] = mcIdx;

        isEsdIdxSignal[esdIdxTrack] = isMcIdxSignal[mcIdx];
        isEsdIdxDaughterOfSecondaryV0[esdIdxTrack] = isMcIdxDaughterOfSecondaryV0[mcIdx];
        isEsdIdxDaughterOfTrueV0[esdIdxTrack] = isMcIdxDaughterOfTrueV0[mcIdx];

        /* Apply track selection */

        if (!PassesTrackSelection(track, n_sigma_proton, n_sigma_kaon, n_sigma_pion)) continue;

        /* Fill containers -- related to findable anti-sexaquarks */

        if (isMcIdxSignal[mcIdx]) {
            getReactionIdx_fromEsdIdx[esdIdxTrack] =
                isMcIdxDaughterOfSignalV0[mcIdx] ? getReactionIdx_fromMcIdx[getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdx]] : getReactionIdx_fromMcIdx[mcIdx];
            getEsdIdx_fromReactionIdx[getReactionIdx_fromEsdIdx[esdIdxTrack]].push_back(esdIdxTrack);
        }

        /* Fill histograms */

        pt = track->Pt();

        f2DHist_TPC_Signal->Fill(track->P(), track->GetTPCsignal());

        if (n_sigma_proton < kMaxNSigma_Proton && track->Charge() < 0.) {
            esdIndicesOfAntiProtonTracks.push_back(esdIdxTrack);
            fHist_Rec_AntiProton_Pt->Fill(pt);
            if (isEsdIdxSignal[esdIdxTrack]) fHist_Rec_Signal_AntiProton_Pt->Fill(pt);
        }

        if (n_sigma_kaon < kMaxNSigma_Kaon && track->Charge() > 0.) {
            esdIndicesOfPosKaonTracks.push_back(esdIdxTrack);
            fHist_Rec_PosKaon_Pt->Fill(pt);
            if (isEsdIdxSignal[esdIdxTrack]) fHist_Rec_Signal_PosKaon_Pt->Fill(pt);  // PENDING
        }

        if (n_sigma_pion < kMaxNSigma_Pion && track->Charge() > 0.) {
            esdIndicesOfPiPlusTracks.push_back(esdIdxTrack);
            fHist_Rec_PiPlus_Pt->Fill(pt);
            if (isEsdIdxSignal[esdIdxTrack]) fHist_Rec_Signal_PiPlus_Pt->Fill(pt);
        }

        if (n_sigma_pion < kMaxNSigma_Pion && track->Charge() < 0.) {
            esdIndicesOfPiMinusTracks.push_back(esdIdxTrack);
            fHist_Rec_PiMinus_Pt->Fill(pt);
            if (isEsdIdxSignal[esdIdxTrack]) fHist_Rec_Signal_PiMinus_Pt->Fill(pt);
        }

        /* Fill vectors */

        mcIdxOfTrueV0 = getMcIdxOfTrueV0_fromMcIdxOfDau[mcIdx];
        getMcIdxOfTrueV0_fromEsdIdx[esdIdxTrack] = mcIdxOfTrueV0;

        if (track->Charge() < 0.) {
            getEsdIdxOfRecNegDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0] = esdIdxTrack;
        } else {
            getEsdIdxOfRecPosDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0] = esdIdxTrack;
        }

    }  // end of loop over tracks
}

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`
 - Output: `n_sigma_pion`, `n_sigma_kaon`, `n_sigma_proton`
 */
Bool_t AliAnalysisTaskSexaquark::PassesTrackSelection(AliESDtrack* track, Float_t& n_sigma_proton, Float_t& n_sigma_kaon, Float_t& n_sigma_pion) {

    n_sigma_proton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    n_sigma_kaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
    n_sigma_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));

    // >> particle identification
    if (n_sigma_pion > kMaxNSigma_Pion && n_sigma_kaon > kMaxNSigma_Kaon && n_sigma_proton > kMaxNSigma_Proton) {
        return kFALSE;
    }

    /* Fill bookkeeping histograms */
    /* -- PENDING: in a sim.set of reactionID='A', signal K+ appear where they shouldn't exist, because of signal particles that were mis-id */
    /* -- Note: should I only, besides if they're signal or not, use their true id? I think so, because it's a bookkeeping of TRUE particles*/

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(3);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(3);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(3);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(3);

    // >> eta
    if (TMath::Abs(track->Eta()) > 0.9) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(4);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(4);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(4);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(4);

    // >> TPC clusters
    if (track->GetTPCNcls() < 50) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(5);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(5);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(5);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(5);

    // >> chi2 per TPC cluster
    if (track->GetTPCchi2() / (Float_t)track->GetTPCNcls() > 7.0) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(6);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(6);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(6);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(6);

    // >> require TPC refit
    if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(7);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(7);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(7);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(7);

    // >> reject kinks
    if (track->GetKinkIndex(0) != 0) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(8);
    if (n_sigma_kaon < kMaxNSigma_Kaon && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(8);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(8);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(8);

    // >> reject low-pT pions
    if (track->Pt() < 0.2 && n_sigma_pion < kMaxNSigma_Pion) {
        return kFALSE;
    }

    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(9);
    if (n_sigma_pion < kMaxNSigma_Pion && isEsdIdxSignal[track->GetID()] && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(9);

    return kTRUE;
}

/*                                                       */
/**  V0s and Anti-Sexaquarks -- That Require True Info  **/
/*** ================================================= ***/

/*
 Find all true V0s that had both of their daughters reconstructed. Useful for bookkeeping.
 - Uses: `fTrueV0s_McIdx`, `fMC`, `fPDG`, `fESD`, `fMC_PrimaryVertex`
 - Notes: to get any rec. value from the true V0s, I need to do vertexing first. That means:
   (a) true vertexing & true values => trivial and already done by `ProcessMCGen()`
   (b) true vertexing & rec. values => feels artificial...
   (c) standard vertexing (a-la official V0 finders)
   (d) KF vertexing
   Considering (c) and (d), then this function should be used only as bookkeeping, and for calling the real "reconstructers"
   I'm using (a) for now. (PENDING)
*/
void AliAnalysisTaskSexaquark::ProcessFindableV0s() {

    AliMCParticle* mcV0;
    AliMCParticle* mcNegDaughter;

    Bool_t is_secondary;  // PENDING
    Bool_t is_signal;

    Float_t mass;
    Float_t radius;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    /* Loop over true V0s */

    for (Int_t& mcIdxOfTrueV0 : mcIndicesOfTrueV0s) {

        mcV0 = (AliMCParticle*)fMC->GetTrack(mcIdxOfTrueV0);

        /* Check if both daughters were reconstructed */

        if (!getEsdIdxOfRecNegDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0] || !getEsdIdxOfRecPosDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0]) {
            continue;
        }

        is_signal = mcV0->MCStatusCode() >= 600 && mcV0->MCStatusCode() < 620;

        /* Get the V0 decay vertex from one of the true daughters */

        mcNegDaughter = (AliMCParticle*)fMC->GetTrack(getMcIdxOfNegDau_fromMcIdxOfTrueV0[mcIdxOfTrueV0]);

        /* Get V0 properties */

        mass = fPDG.GetParticle(mcV0->PdgCode())->Mass();  // PENDING
        radius = TMath::Sqrt(mcNegDaughter->Xv() * mcNegDaughter->Xv() + mcNegDaughter->Yv() * mcNegDaughter->Yv());
        cpa_wrt_pv = CosinePointingAngle(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                         fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());
        dca_wrt_pv = LinePointDCA(mcV0->Px(), mcV0->Py(), mcV0->Pz(), mcNegDaughter->Xv(), mcNegDaughter->Yv(), mcNegDaughter->Zv(),
                                  fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(), fMC_PrimaryVertex->GetZ());

        /* Fill bookkeeping histograms */

        if (mcV0->PdgCode() == -3122) {
            fHist_True_AntiLambda_Bookkeep->Fill(3);
            fHist_Findable_True_AntiLambda_Mass->Fill(mass);
            fHist_Findable_True_AntiLambda_Radius->Fill(radius);
            fHist_Findable_True_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Findable_True_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

            if (!is_signal) continue;

            fHist_True_AntiLambda_Bookkeep->Fill(4);
            fHist_Findable_Signal_AntiLambda_Mass->Fill(mass);
            fHist_Findable_Signal_AntiLambda_Radius->Fill(radius);
            fHist_Findable_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Findable_Signal_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
        }

        if (mcV0->PdgCode() == 310) {
            fHist_True_KaonZeroShort_Bookkeep->Fill(3);
            fHist_Findable_True_KaonZeroShort_Mass->Fill(mass);
            fHist_Findable_True_KaonZeroShort_Radius->Fill(radius);
            fHist_Findable_True_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Findable_True_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

            if (!is_signal) continue;

            fHist_True_KaonZeroShort_Bookkeep->Fill(4);
            fHist_Findable_Signal_KaonZeroShort_Mass->Fill(mass);
            fHist_Findable_Signal_KaonZeroShort_Radius->Fill(radius);
            fHist_Findable_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Findable_Signal_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
        }
    }
}

/*
  Find the anti-sexaquarks that have all of their *final-state particles* reconstructed.
  - For example:
    -- Reaction Channel A: 1 anti-proton, 2 pi+, 1 pi-
    -- Reaction Channel D: 1 anti-proton, 1 pi+, 1 K+
    -- Reaction Channel E: 1 anti-proton, 2 pi+, 1 pi-, 1 K+
    -- Reaction Channel H: PENDING
*/
void AliAnalysisTaskSexaquark::ProcessFindableSexaquarks() {

    Int_t esdIdxFinalStateParticle;

    AliMCParticle* mcFinalStateParticle;
    Int_t mcIdxFinalStateParticle;

    AliESDtrack* finalStateTrack;
    Float_t n_sigma_proton;
    Float_t n_sigma_kaon;
    Float_t n_sigma_pion;

    // DEBUG
    /*
    std::map<Int_t, Int_t> counter_mc_idx;
    */

    /* Declare TLorentzVectors */

    TLorentzVector lvFinalStateParticle;

    TLorentzVector lvStruckNucleon(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

    TLorentzVector lvAntiSexaquark;

    /* Loop over interactions */

    for (Int_t reactionIdx = 600; reactionIdx < 620; reactionIdx++) {

        /* Reset anti-sexaquark TLorentzVector */

        lvAntiSexaquark.SetPxPyPzE(0., 0., 0., 0.);

        if (getEsdIdx_fromReactionIdx[reactionIdx].size() < fFinalStateProductsPDG.size()) continue;

        /* First loop, to count how many times a MC Idx repeats */

        // DEBUG
        /*
        for (Int_t finalStateIdx = 0; finalStateIdx < (Int_t)getEsdIdx_fromReactionIdx[reactionIdx].size(); finalStateIdx++) {
            esdIdxFinalStateParticle = getEsdIdx_fromReactionIdx[reactionIdx][finalStateIdx];
            mcIdxFinalStateParticle = getMcIdx_fromEsdIdx[esdIdxFinalStateParticle];
            counter_mc_idx[mcIdxFinalStateParticle]++;
        }
        */

        /* Second loop, to get the reconstructed final-state particles info */

        for (Int_t finalStateIdx = 0; finalStateIdx < (Int_t)getEsdIdx_fromReactionIdx[reactionIdx].size(); finalStateIdx++) {

            esdIdxFinalStateParticle = getEsdIdx_fromReactionIdx[reactionIdx][finalStateIdx];
            finalStateTrack = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxFinalStateParticle));

            n_sigma_proton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(finalStateTrack, AliPID::kProton));
            n_sigma_kaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(finalStateTrack, AliPID::kKaon));
            n_sigma_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(finalStateTrack, AliPID::kPion));

            mcIdxFinalStateParticle = getMcIdx_fromEsdIdx[esdIdxFinalStateParticle];
            mcFinalStateParticle = (AliMCParticle*)fMC->GetTrack(mcIdxFinalStateParticle);

            /*
            // DEBUG
            if (counter_mc_idx[mcIdxFinalStateParticle] > 1) {
                AliInfoF("reaction %i", reactionIdx);
                AliInfoF("=> \033[1;31mmcIdx %i\033[0m, pdg %i", mcIdxFinalStateParticle, mcFinalStateParticle->PdgCode());
                AliInfoF("=> mc_px %.2f, mc_py %.2f, mc_pz %.2f, \033[1;31mmc_pt %.2f\033[0m", mcFinalStateParticle->Px(), mcFinalStateParticle->Py(),
                        mcFinalStateParticle->Pz(), mcFinalStateParticle->Pt());
                AliInfoF("=> esdIdx %i, n_sigma_proton %.0f, n_sigma_kaon %.0f, n_sigma_pion %.0f", esdIdxFinalStateParticle, n_sigma_proton,
                        n_sigma_kaon, n_sigma_pion);
                AliInfoF("=> px %.2f, py %.2f, pz %.2f, \033[1;31mpt %.2f\033[0m", finalStateTrack->Px(), finalStateTrack->Py(),
                        finalStateTrack->Pz(), finalStateTrack->Pt());
                AliInfoF("=> global chi2 %.2f, tpc_mom %.2f", finalStateTrack->GetGlobalChi2(), finalStateTrack->GetTPCmomentum());
            } else {
                AliInfoF("reaction %i", reactionIdx);
                AliInfoF("=> mcIdx %i, pdg %i", mcIdxFinalStateParticle, mcFinalStateParticle->PdgCode());
                AliInfoF("=> mc_px %.2f, mc_py %.2f, mc_pz %.2f, mc_pt %.2f", mcFinalStateParticle->Px(), mcFinalStateParticle->Py(),
                        mcFinalStateParticle->Pz(), mcFinalStateParticle->Pt());
                AliInfoF("=> esdIdx %i, n_sigma_proton %.0f, n_sigma_kaon %.0f, n_sigma_pion %.0f", esdIdxFinalStateParticle, n_sigma_proton,
                        n_sigma_kaon, n_sigma_pion);
                AliInfoF("=> px %.2f, py %.2f, pz %.2f, pt %.2f", finalStateTrack->Px(), finalStateTrack->Py(), finalStateTrack->Pz(),
                        finalStateTrack->Pt());
                AliInfoF("=> global chi2 %.2f, tpc_mom %.2f", finalStateTrack->GetGlobalChi2(), finalStateTrack->GetTPCmomentum());
            }
            */

            lvFinalStateParticle.SetPxPyPzE(mcFinalStateParticle->Px(), mcFinalStateParticle->Py(), mcFinalStateParticle->Pz(),
                                            mcFinalStateParticle->E());

            lvAntiSexaquark = lvAntiSexaquark + lvFinalStateParticle;
        }

        lvAntiSexaquark = lvAntiSexaquark - lvStruckNucleon;

        /* Fill histograms  */

        fHist_Findable_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
        fHist_Findable_AntiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
        fHist_Findable_AntiSexaquark_DCAwrtPV->Fill(  //
            LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), mcFinalStateParticle->Xv(), mcFinalStateParticle->Yv(),
                         mcFinalStateParticle->Zv(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(),
                         fMC_PrimaryVertex->GetZ()));  // BUG: shouldn't be the origin of the final-state particles! it should be the origin of the
                                                       // products!
        fHist_Findable_AntiSexaquark_CPAwrtPV->Fill(   //
            CosinePointingAngle(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), mcFinalStateParticle->Xv(),
                                mcFinalStateParticle->Yv(), mcFinalStateParticle->Zv(), fMC_PrimaryVertex->GetX(), fMC_PrimaryVertex->GetY(),
                                fMC_PrimaryVertex->GetZ()));  // BUG: shouldn't be the origin of the final-state particles! it should be the origin of
                                                              // the products!
        fHist_Findable_AntiSexaquark_Radius->Fill(            //
            TMath::Sqrt(TMath::Power(mcFinalStateParticle->Xv(), 2) +
                        TMath::Power(mcFinalStateParticle->Yv(), 2)));  // BUG: shouldn't be the origin of the final-state particles! it should be the
                                                                        // origin of the products!
    }
}

/*                             */
/**  V0s -- ALICE V0 Finders  **/
/*** ======================= ***/

/*
 Loop over all the V0s that were stored in the ESD file by the Official V0 Finders.
 - Uses: `fESD`, `fMC`, `fPDG`, `fPrimaryVertex`
 - Containers: `esdAntiLambdas`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`, `isAntiLambdaIdxATrueV0`,
 `isEsdIdxDaughterOfTrueV0`, `getMcIdxOfTrueV0_fromEsdIdx`, `isMcIdxSignal`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `isKaonZeroShortIdxATrueV0`
*/
void AliAnalysisTaskSexaquark::OfficialV0Finder(Bool_t online) {

    AliESDv0* V0;

    Int_t esdIdxPos;
    Int_t esdIdxNeg;

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;
    Float_t nsigma_pos_pion;
    Float_t nsigma_neg_pion;
    Float_t nsigma_antiproton;

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    Double_t pos_energy;
    Double_t neg_energy;

    Float_t mass;
    Float_t radius;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;

    AliMCParticle* mcV0;

    AliESDv0 outputV0;

    /* Loop over V0s that are stored in the ESD */

    for (Int_t esdIdxV0 = 0; esdIdxV0 < fESD->GetNumberOfV0s(); esdIdxV0++) {

        V0 = fESD->GetV0(esdIdxV0);

        /* Choose between offline and online (on-the-fly) V0s */

        if (!online && V0->GetOnFlyStatus()) {
            continue;
        }

        if (online && !V0->GetOnFlyStatus()) {
            continue;
        }

        esdIdxNeg = V0->GetNindex();
        esdIdxPos = V0->GetPindex();

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

        if (!PassesAntiLambdaCuts_OfficialCustom(V0, neg_track, pos_track) && !PassesKaonZeroShortCuts_OfficialCustom(V0, neg_track, pos_track)) {
            continue;
        }

        /* Get V0 properties */

        V0->GetXYZ(V0_X, V0_Y, V0_Z);
        V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
        V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

        radius = TMath::Sqrt(V0_X * V0_X + V0_Y * V0_Y);
        cpa_wrt_pv = V0->GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
        dca_wrt_pv =
            LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

        /* Case: anti-lambda -> anti-proton + pi+ */

        if (PassesAntiLambdaCuts_OfficialCustom(V0, neg_track, pos_track)) {

            pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(211)->Mass() * fPDG.GetParticle(211)->Mass());
            neg_energy = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(-2212)->Mass() * fPDG.GetParticle(-2212)->Mass());
            mass =
                TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) - V0->Px() * V0->Px() - V0->Py() * V0->Py() - V0->Pz() * V0->Pz());

            /* Store V0 */

            outputV0 = *V0;
            outputV0.ChangeMassHypothesis(-3122);
            esdAntiLambdas.push_back(outputV0);

            getEsdIdxOfNegDau_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = esdIdxNeg;
            getEsdIdxOfPosDau_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = esdIdxPos;

            isAntiLambdaIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                             getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

            /* Fill histograms */

            fHist_Found_AntiLambda_Mass->Fill(mass);
            fHist_Found_AntiLambda_Radius->Fill(radius);
            fHist_Found_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Found_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

            if (isAntiLambdaIdxATrueV0.back()) {

                getMcIdx_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                if (mcV0->PdgCode() == -3122) {

                    fHist_Found_True_AntiLambda_Mass->Fill(mass);
                    fHist_Found_True_AntiLambda_Radius->Fill(radius);
                    fHist_Found_True_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                    fHist_Found_True_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

                    if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                        fHist_Found_Signal_AntiLambda_Mass->Fill(mass);
                        fHist_Found_Signal_AntiLambda_Radius->Fill(radius);
                        fHist_Found_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_Signal_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
                    }
                }
            }
        }

        /* Case: K0S -> pi- + pi+ */

        if (PassesKaonZeroShortCuts_OfficialCustom(V0, neg_track, pos_track)) {

            pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(211)->Mass() * fPDG.GetParticle(211)->Mass());
            neg_energy = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(-211)->Mass() * fPDG.GetParticle(-211)->Mass());
            mass = TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) -  //
                               V0->Px() * V0->Px() - V0->Py() * V0->Py() - V0->Pz() * V0->Pz());

            /* Store V0 */

            outputV0 = *V0;
            outputV0.ChangeMassHypothesis(310);
            esdKaonsZeroShort.push_back(outputV0);

            getEsdIdxOfNegDau_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = esdIdxNeg;
            getEsdIdxOfPosDau_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = esdIdxPos;

            isKaonZeroShortIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                                getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

            /* Fill histograms */

            fHist_Found_KaonZeroShort_Mass->Fill(mass);
            fHist_Found_KaonZeroShort_Radius->Fill(radius);
            fHist_Found_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Found_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

            if (isKaonZeroShortIdxATrueV0.back()) {

                getMcIdx_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                if (mcV0->PdgCode() == 310) {

                    fHist_Found_True_KaonZeroShort_Mass->Fill(mass);
                    fHist_Found_True_KaonZeroShort_Radius->Fill(radius);
                    fHist_Found_True_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                    fHist_Found_True_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

                    if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                        fHist_Found_Signal_KaonZeroShort_Mass->Fill(mass);
                        fHist_Found_Signal_KaonZeroShort_Radius->Fill(radius);
                        fHist_Found_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_Signal_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
                    }
                }
            }
        }

    }  // end of loop over V0s
}

/*
 Custom V0 Finder.
 (Copied and adapted from `AliRoot/STEER/ESD/AliV0vertexer.cxx`)
 - Input: PENDING
 - Output: PENDING
*/
void AliAnalysisTaskSexaquark::CustomV0Finder(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {

    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Int_t nentr = fESD->GetNumberOfTracks();
    TArrayI neg(nentr);
    TArrayI pos(nentr);

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    Double_t radius;
    Double_t transv_momentum;
    Double_t pos_energy;
    Double_t neg_energy;
    Double_t mass;
    Double_t cpa_wrt_pv;
    Double_t dca_wrt_pv;

    AliMCParticle* mcV0;

    // define cuts
    Double_t COV0F_DNmin = 0.1;            // (d=0.1) min imp parameter for the negative daughter (depends on PV!)
    Double_t COV0F_DPmin = 0.1;            // (d=0.1) min imp parameter for the positive daughter (depends on PV!)
    Double_t COV0F_Chi2max = 33.;          // (d=33.) max chi2
    Double_t COV0F_DCAmax = 0.6;           // (d=1.) max DCA between the daughter tracks
    Double_t COV0F_CPAmin = 0.4;           // (d=0.998) max cosine of V0's pointing angle (depends on PV!)
    Double_t COV0F_Rmin = 50.;             // (d=0.9) min radius of the V0 radius
    Double_t COV0F_Rmax = 200.;            // (d=100) max radius of the V0 radius
    Double_t COV0F_Etamax_V0 = 1.0;        // (d=1.) max eta of V0
    Double_t COV0F_DCAmin_V0_wrtPV = 5.;   // min DCA between V0 and PV
    Double_t COV0F_Etamax_NegDau = 1.5;    // max eta of negative daughter
    Double_t COV0F_Etamax_PosDau = 1.5;    // max eta of positive daughter
    Int_t COV0F_NClustersTPCmin_Dau = 50;  // minimun number of clusters in the TPC
    Double_t COV0F_Ptmin = 1.0;            // min Pt of V0

    Double_t b = fESD->GetMagneticField();

    if (nentr < 2) {
        return;
    }

    // (addition)
    Bool_t UseImprovedFinding = kTRUE;

    // (addition)
    // [based on `AliV0HypSel::AccountBField(b)` and `AliV0HypSel::SetBFieldCoef(v)`]
    // account effect of B-field on pT resolution, ignoring the fact that the V0 mass resolution
    // is only partially determined by the prongs pT resolution
    const Float_t kNomField = 5.00668e+00;
    Float_t babs = TMath::Abs(b);
    Float_t BFieldCoef;
    if (babs > 1E-3) {
        BFieldCoef = kNomField / babs > 0. ? kNomField / babs : 1.0;
    }

    // (addition)
    const Int_t nHypSel = 6;
    TString v0_hyp_name[nHypSel] = {"gamma", "K0", "Lambda", "antiLambda", "HyperTriton", "antiHyperTriton"};
    Float_t v0_hyp_m0[nHypSel] = {0.5486e-3, 139.570e-3, 938.272e-3, 139.570e-3, 2.8092, 139.570e-3};
    Float_t v0_hyp_m1[nHypSel] = {0.5486e-3, 139.570e-3, 139.570e-3, 938.272e-3, 139.570e-3, 2.8092};
    Float_t v0_hyp_mass[nHypSel] = {1.099e-3, 497.7e-3, 1115.683e-3, 1115.683e-3, 2.992, 2.992};
    Float_t v0_hyp_sigma[nHypSel] = {0.001, 0.003, 0.001, 0.001, 0.0025, 0.0025};
    Float_t v0_hyp_cf0[nHypSel] = {20, 20, 20, 20, 14, 14};
    Float_t v0_hyp_cf1[nHypSel] = {0.6, 0.07, 0.07, 0.07, 0.07, 0.07};
    Float_t v0_hyp_nsig[nHypSel] = {0.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    Float_t v0_hyp_margin[nHypSel] = {0.0, 0.5, 0.5, 0.5, 0.5, 0.5};

    /*  Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    }

    /* Loop over all possible pairs of tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            /* Get tracks */

            neg_track = fESD->GetTrack(esdIdxNeg);
            pos_track = fESD->GetTrack(esdIdxPos);

            /* Cut on eta of daughters */

            if (TMath::Abs(pos_track->Eta()) > COV0F_Etamax_PosDau) {
                continue;
            }

            if (TMath::Abs(neg_track->Eta()) > COV0F_Etamax_NegDau) {
                continue;
            }

            AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;

            Double_t xn, xp;
            if (UseImprovedFinding && Preoptimize(ntp, ptp, &xn, &xp, b)) {
                // Move tracks to a better position if that helps
                nt.PropagateTo(xn, b);
                pt.PropagateTo(xp, b);
            }

            /* Cut on DCA between daughters */

            Double_t dca = ntp->GetDCA(ptp, b, xn, xp);
            if (dca > COV0F_DCAmax) {
                continue;
            }

            // idk what's this...
            if ((xn + xp) > 2 * COV0F_Rmax) {
                continue;
            }
            if ((xn + xp) < 2 * COV0F_Rmin) {
                continue;
            }

            nt.PropagateTo(xn, b);
            pt.PropagateTo(xp, b);

            /* Prepare the AliESDv0 */

            AliESDv0 vertex(nt, esdIdxNeg, pt, esdIdxPos);

            if (UseImprovedFinding) {
                vertex.Refit();  // imp pos + cov mat
            }

            /* Cut on Chi2 */

            if (vertex.GetChi2V0() > COV0F_Chi2max) {
                continue;
            }

            vertex.SetDcaV0Daughters(dca);

            /* Cut on CPA w.r.t. PV */

            cpa_wrt_pv = vertex.GetV0CosineOfPointingAngle(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            vertex.SetV0CosineOfPointingAngle(cpa_wrt_pv);
            if (cpa_wrt_pv < COV0F_CPAmin) {
                continue;
            }

            /* Cut on radius */

            radius = TMath::Sqrt(vertex.Xv() * vertex.Xv() + vertex.Yv() * vertex.Yv());
            if (radius < COV0F_Rmin) {
                continue;
            }
            if (radius > COV0F_Rmax) {
                continue;
            }

            /* Load coordinates and momentum of each component of the V0 */

            vertex.GetXYZ(V0_X, V0_Y, V0_Z);

            /* Cut on DCA w.r.t. PV */

            dca_wrt_pv = LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                      fPrimaryVertex->GetZ());
            if (dca_wrt_pv < COV0F_DCAmin_V0_wrtPV) {
                continue;
            }

            // (addition)
            if (nHypSel) {  // do we select particular hypotheses?
                Bool_t reject = kTRUE;
                Float_t v0_pt = vertex.Pt();
                for (int hh = 0; hh < nHypSel; hh++) {
                    Double_t v0_eff_mass = vertex.GetEffMassExplicit(v0_hyp_m0[hh], v0_hyp_m1[hh]);
                    Float_t v0_mass_margin = v0_hyp_nsig[hh] * BFieldCoef * v0_hyp_sigma[hh] * (v0_hyp_cf0[hh] + v0_pt * v0_hyp_cf1[hh])  //
                                             + v0_hyp_margin[hh];
                    if (TMath::Abs(v0_eff_mass - v0_hyp_mass[hh]) < v0_mass_margin) {
                        reject = kFALSE;
                        break;
                    }
                }
                if (reject) {
                    continue;
                }
            }

            /* Cut on Eta */

            if (TMath::Abs(vertex.Eta()) > COV0F_Etamax_V0) {
                continue;
            }

            /* Cut on Pt */

            transv_momentum = vertex.Pt();
            if (transv_momentum < COV0F_Ptmin) {
                continue;
            }

            /*
            // get NSigmas for PID
            nsigma_pos_pion = fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion);
            nsigma_neg_pion = fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kPion);
            nsigma_antiproton = fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kProton);

            // (cut) on PID
            // positive daughter must be a positive pion
            if (TMath::Abs(nsigma_pos_pion) > kMaxNSigma_Pion) {
                continue;
            }
            // negative daughter must be a negative pion or an anti-lambda
            if (TMath::Abs(nsigma_neg_pion) > kMaxNSigma_Pion && TMath::Abs(nsigma_antiproton) > kMaxNSigma_Proton) {
                continue;
            }

            // (cut) on sign
            if (pos_track->GetSign() == neg_track->GetSign()) {
                continue;
            }
            */

            /* Calculate energy */

            // get momentum of negative and positive components at V0 vertex
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            pos_energy =
                TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(pdgTrackPos)->Mass() * fPDG.GetParticle(pdgTrackPos)->Mass());
            neg_energy =
                TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(pdgTrackNeg)->Mass() * fPDG.GetParticle(pdgTrackNeg)->Mass());
            mass = TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) - vertex.Px() * vertex.Px() - vertex.Py() * vertex.Py() -
                               vertex.Pz() * vertex.Pz());

            /* Case: anti-lambda -> anti-proton + pi+ */

            if (pdgV0 == -3122) {

                /* Store V0 */

                vertex.ChangeMassHypothesis(pdgV0);
                esdAntiLambdas.push_back(vertex);

                getEsdIdxOfNegDau_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = esdIdxNeg;
                getEsdIdxOfPosDau_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = esdIdxPos;

                isAntiLambdaIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                                 getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

                /* Fill histograms */

                fHist_Found_AntiLambda_Mass->Fill(mass);
                fHist_Found_AntiLambda_Radius->Fill(radius);
                fHist_Found_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_Found_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

                if (isAntiLambdaIdxATrueV0.back()) {

                    getMcIdx_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                    mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                    if (mcV0->PdgCode() == pdgV0) {

                        fHist_Found_True_AntiLambda_Mass->Fill(mass);
                        fHist_Found_True_AntiLambda_Radius->Fill(radius);
                        fHist_Found_True_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_True_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

                        if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                            fHist_Found_Signal_AntiLambda_Mass->Fill(mass);
                            fHist_Found_Signal_AntiLambda_Radius->Fill(radius);
                            fHist_Found_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                            fHist_Found_Signal_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
                        }
                    }
                }
            }

            /* Case: K0S -> pi- + pi+ */

            if (pdgV0 == 310) {

                /* Store V0 */

                vertex.ChangeMassHypothesis(pdgV0);
                esdKaonsZeroShort.push_back(vertex);

                getEsdIdxOfNegDau_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = esdIdxNeg;
                getEsdIdxOfPosDau_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = esdIdxPos;

                isKaonZeroShortIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                                    getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

                /* Fill histograms */

                fHist_Found_KaonZeroShort_Mass->Fill(mass);
                fHist_Found_KaonZeroShort_Radius->Fill(radius);
                fHist_Found_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_Found_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

                if (isKaonZeroShortIdxATrueV0.back()) {

                    getMcIdx_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                    mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                    if (mcV0->PdgCode() == pdgV0) {

                        fHist_Found_True_KaonZeroShort_Mass->Fill(mass);
                        fHist_Found_True_KaonZeroShort_Radius->Fill(radius);
                        fHist_Found_True_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_True_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

                        if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                            fHist_Found_Signal_KaonZeroShort_Mass->Fill(mass);
                            fHist_Found_Signal_KaonZeroShort_Radius->Fill(radius);
                            fHist_Found_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                            fHist_Found_Signal_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
                        }
                    }
                }
            }

        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*
 Apply cuts to an (anti-lambda) candidate, reconstructed from an (anti-proton, pi+) pair. Relevant for channel A and D.
 - Input: `v0`, `neg_track`, `pos_track`
*/
Bool_t AliAnalysisTaskSexaquark::PassesAntiLambdaCuts_OfficialCustom(AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track) {
    // PENDING
    return kTRUE;
}

/*
 Apply cuts to a (K0S) candidate, reconstructed from a (pi-, pi+) pair. Relevant for channel A.
 - Input: `v0`, `neg_track`, `pos_track`
*/
Bool_t AliAnalysisTaskSexaquark::PassesKaonZeroShortCuts_OfficialCustom(AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track) {
    // PENDING
    return kTRUE;
}

/*                                     */
/**  Sexaquark -- Geometrical Finder  **/
/*** =============================== ***/

/*
 Using two AliESDv0 objects, reconstruct the anti-sexaquark candidate.
 - Containers: `esdAntiLambdas`, `esdKaonsZeroShort`, `getEsdIdxOfNegDau_fromAntiLambdaIdx`, `getEsdIdxOfNegDau_fromKaonZeroShortIdx`,
 `getEsdIdxOfPosDau_fromKaonZeroShortIdx`, `getEsdIdxOfPosDau_fromAntiLambdaIdx`
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelA_Geo() {

    /* Declare AliESDv0 objects */

    AliESDv0 AntiLambda;
    AliESDv0 KaonZeroShort;

    /* Declare TVector3s  */

    TVector3 AntiLambda_momentum, KaonZeroShort_momentum;
    TVector3 AntiLambda_vertex, KaonZeroShort_vertex;

    TVector3 AntiSexaquark_vertex;

    /* Declare TLorentzVectors */

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvAntiSexaquark;

    TLorentzVector lvStruckNucleon(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

    /* Declare auxiliar variables */

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    Float_t dca_between_daughters;
    TVector3 pca1, pca2;
    Float_t dca_wrt_pv;
    Float_t cpa_wrt_pv;
    Bool_t is_signal;

    /* Loop over all possible pairs of V0s */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)esdAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)esdKaonsZeroShort.size(); idxKaonZeroShort++) {

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Get AliESDv0 objects */

            AntiLambda = esdAntiLambdas[idxAntiLambda];
            KaonZeroShort = esdKaonsZeroShort[idxKaonZeroShort];

            AntiLambda_momentum.SetXYZ(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz());
            KaonZeroShort_momentum.SetXYZ(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz());

            AntiLambda_vertex.SetXYZ(AntiLambda.Xv(), AntiLambda.Yv(), AntiLambda.Zv());
            KaonZeroShort_vertex.SetXYZ(KaonZeroShort.Xv(), KaonZeroShort.Yv(), KaonZeroShort.Zv());

            dca_between_daughters =
                Calculate_TwoLinesDCA_v2(AntiLambda_momentum, AntiLambda_vertex, KaonZeroShort_momentum, KaonZeroShort_vertex, pca1, pca2);

            // PassesSexaquarkCuts_ChannelA_Geo(AntiLambda, KaonZeroShort);  // PENDING

            /* Reconstruct anti-sexaquark candidate */

            lvAntiLambda.SetPxPyPzE(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz(), AntiLambda.E());
            lvKaonZeroShort.SetPxPyPzE(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz(), KaonZeroShort.E());

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;
            AntiSexaquark_vertex = 0.5 * (pca1 + pca2);

            dca_wrt_pv =
                LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), AntiSexaquark_vertex.X(), AntiSexaquark_vertex.Y(),
                             AntiSexaquark_vertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), AntiSexaquark_vertex.X(),
                                             AntiSexaquark_vertex.Y(), AntiSexaquark_vertex.Z(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                             fPrimaryVertex->GetZ());

            /* Fill histograms */

            fHist_Found_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
            fHist_Found_AntiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
            fHist_Found_AntiSexaquark_DCAbtwDau->Fill(dca_between_daughters);
            fHist_Found_AntiSexaquark_DCAwrtPV->Fill(dca_wrt_pv);
            fHist_Found_AntiSexaquark_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Found_AntiSexaquark_Radius->Fill(AntiSexaquark_vertex.Perp());

            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            if (is_signal) {
                fHist_Found_Signal_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
                fHist_Found_Signal_AntiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
                fHist_Found_Signal_AntiSexaquark_DCAbtwDau->Fill(dca_between_daughters);
                fHist_Found_Signal_AntiSexaquark_DCAwrtPV->Fill(dca_wrt_pv);
                fHist_Found_Signal_AntiSexaquark_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_Found_Signal_AntiSexaquark_Radius->Fill(AntiSexaquark_vertex.Perp());
            }
        }
    }

    return;
}

/*
 PENDING
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelD_Geo() {
    // PENDING
    return;
}

/*
 PENDING
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelE_Geo() {
    // PENDING
    return;
}

/*
 Apply anti-sexaquark candidate cuts.
 - Input: `AntiLambda`, `KaonZeroShort`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA_Geo(AliESDv0 AntiLambda, AliESDv0 KaonZeroShort) {

    /* this function is PENDING, only use for printing and filling hists for the moment being... */

    TVector3 AntiLambda_momentum(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz());
    TVector3 KaonZeroShort_momentum(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz());

    TVector3 AntiLambda_vertex(AntiLambda.Xv(), AntiLambda.Yv(), AntiLambda.Zv());
    TVector3 KaonZeroShort_vertex(KaonZeroShort.Xv(), KaonZeroShort.Yv(), KaonZeroShort.Zv());

    Double_t dca;
    TVector3 pca1, pca2;

    AliInfoF("anti_lambda mom vertex radius = (%f, %f, %f) (%f, %f, %f) %f", AntiLambda_momentum.X(), AntiLambda_momentum.Y(),
             AntiLambda_momentum.Z(), AntiLambda_vertex.X(), AntiLambda_vertex.Y(), AntiLambda_vertex.Z(), AntiLambda_vertex.Perp());
    AliInfoF("kaon_zero_short mom vertex radius = (%f, %f, %f) (%f, %f, %f) %f", KaonZeroShort_momentum.X(), KaonZeroShort_momentum.Y(),
             KaonZeroShort_momentum.Z(), KaonZeroShort_vertex.X(), KaonZeroShort_vertex.Y(), KaonZeroShort_vertex.Z(), KaonZeroShort_vertex.Perp());

    dca = Calculate_TwoLinesDCA_v1(AntiLambda_momentum, AntiLambda_vertex, KaonZeroShort_momentum, KaonZeroShort_vertex, pca1, pca2);
    dca = Calculate_TwoLinesDCA_v2(AntiLambda_momentum, AntiLambda_vertex, KaonZeroShort_momentum, KaonZeroShort_vertex, pca1, pca2);

    return kTRUE;
}

/*                          */
/**  V0s -- Kalman Filter  **/
/*** ==================== ***/

/*
 Find all true V0s which had both of their daughters reconstructed.
 - Uses: `fESD`, `fMC`, `fMagneticField`
 - Input: `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 - Output: `kfV0s`, `idxDaughters`
*/
void AliAnalysisTaskSexaquark::ReconstructV0s_KF(Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos) {

    AliESDtrack* esdTrackNeg;
    AliESDtrack* esdTrackPos;

    AliMCParticle* mcV0;

    Bool_t successful_daughters_check;  // PENDING: I still don't understand this...

    Double_t impar_neg[2], impar_pos[2];

    Double_t radius;
    Double_t cpa_wrt_pv;
    Double_t dca_wrt_pv;

    /* Declare TLorentzVectors */

    TLorentzVector lvTrackNeg;
    TLorentzVector lvTrackPos;
    TLorentzVector lvV0;

    /* Define cuts */

    Double_t COV0F_DNmin = 0.1;            // PENDING
    Double_t COV0F_DPmin = 0.1;            // PENDING
    Double_t COV0F_Chi2max = 33.;          // PENDING
    Double_t COV0F_DCAmax = 0.6;           // PENDING
    Double_t COV0F_CPAmin = 0.4;           //
    Double_t COV0F_Rmin = 50.;             //
    Double_t COV0F_Rmax = 200.;            //
    Double_t COV0F_Etamax_V0 = 1.0;        //
    Double_t COV0F_DCAmin_V0_wrtPV = 5.;   //
    Double_t COV0F_Etamax_NegDau = 1.5;    //
    Double_t COV0F_Etamax_PosDau = 1.5;    //
    Int_t COV0F_NClustersTPCmin_Dau = 50;  //
    Double_t COV0F_Ptmin = 1.0;            //

    /*  Choose between anti-lambda and K0S */

    std::vector<Int_t> esdIndicesNegTracks;
    std::vector<Int_t> esdIndicesPosTracks;

    if (pdgV0 == -3122) {
        esdIndicesNegTracks = esdIndicesOfAntiProtonTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    } else if (pdgV0 == 310) {
        esdIndicesNegTracks = esdIndicesOfPiMinusTracks;
        esdIndicesPosTracks = esdIndicesOfPiPlusTracks;
    }

    /* Loop over all possible pairs of tracks */

    for (Int_t& esdIdxNeg : esdIndicesNegTracks) {
        for (Int_t& esdIdxPos : esdIndicesPosTracks) {

            if (esdIdxNeg == esdIdxPos) continue;

            esdTrackNeg = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg));
            esdTrackPos = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos));

            /* Kalman Filter */

            KFParticle kfDaughterNeg = CreateKFParticle(*esdTrackNeg, fPDG.GetParticle(pdgTrackNeg)->Mass(), (Int_t)esdTrackNeg->Charge());
            KFParticle kfDaughterPos = CreateKFParticle(*esdTrackPos, fPDG.GetParticle(pdgTrackPos)->Mass(), (Int_t)esdTrackPos->Charge());

            // (protection) PENDING
            successful_daughters_check = kTRUE;
            /*  */ KFParticleMother kfCheckPos;
            kfCheckPos.AddDaughter(kfDaughterNeg);
            if (!kfCheckPos.CheckDaughter(kfDaughterPos)) successful_daughters_check = kFALSE;
            /*  */ KFParticleMother kfCheckNeg;
            if (!kfCheckNeg.CheckDaughter(kfDaughterNeg)) successful_daughters_check = kFALSE;
            kfCheckNeg.AddDaughter(kfDaughterPos);

            KFParticleMother kfV0;
            if (successful_daughters_check) {
                kfV0.SetMassHypo(fPDG.GetParticle(pdgV0)->Mass());
                kfV0.SetConstructMethod(2);
                kfV0.AddDaughter(kfDaughterNeg);
                kfV0.AddDaughter(kfDaughterPos);
            } else {
                kfV0.AddDaughter(kfDaughterNeg);
                kfV0.AddDaughter(kfDaughterPos);
            }
            kfV0.TransportToDecayVertex();

            // get V0 vertex
            Double_t V0Vertex[3] = {kfV0.GetX(), kfV0.GetY(), kfV0.GetZ()};
            Double_t V0VertexErr[3] = {kfV0.GetErrX(), kfV0.GetErrY(), kfV0.GetErrZ()};
            AliESDVertex* esdV0Vertex = new AliESDVertex(V0Vertex, V0VertexErr);

            /* Clone ESD tracks to be able to propagate (~ modify) them. */

            AliESDtrack cloneTrackNeg(*esdTrackNeg);
            AliESDtrack cloneTrackPos(*esdTrackPos);

            /* Propagate V0 daughters to their DCA w.r.t. V0 vertex. Use ESDs for kinematics. */

            cloneTrackNeg.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar_neg);
            cloneTrackPos.PropagateToDCA(esdV0Vertex, fMagneticField, 10., impar_pos);

            /* Reconstruct V0 */

            lvTrackNeg.SetXYZM(cloneTrackNeg.Px(), cloneTrackNeg.Py(), cloneTrackNeg.Pz(), fPDG.GetParticle(pdgTrackNeg)->Mass());
            lvTrackPos.SetXYZM(cloneTrackPos.Px(), cloneTrackPos.Py(), cloneTrackPos.Pz(), fPDG.GetParticle(pdgTrackPos)->Mass());
            lvV0 = lvTrackNeg + lvTrackPos;

            /* Calculate variables */

            radius = TMath::Sqrt(kfV0.GetX() * kfV0.GetX() + kfV0.GetY() * kfV0.GetY());
            cpa_wrt_pv = CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(),
                                             fPrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(), fPrimaryVertex->GetX(),
                                      fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

            /* Apply cuts */
            // PENDING

            /*
            if (pdgV0 == -3122 && !PassesAntiLambdaCuts_KF(kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;
            if (pdgV0 == 310 && !PassesKaonZeroShortCuts_KF(kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;
            if (pdgV0 == 0 && !PassesPionPairCuts_KF(kfV0, kfDaughterNeg, kfDaughterPos, lvV0, lvTrackNeg, lvTrackPos)) continue;
            */

            if (cpa_wrt_pv < COV0F_CPAmin) continue;
            if (radius < COV0F_Rmin || radius > COV0F_Rmax) continue;
            if (dca_wrt_pv < COV0F_DCAmin_V0_wrtPV) continue;
            if (lvTrackNeg.Eta() > COV0F_Etamax_NegDau) continue;
            if (lvTrackPos.Eta() > COV0F_Etamax_PosDau) continue;
            if (lvV0.Eta() > COV0F_Etamax_V0) continue;
            if (lvV0.Pt() < COV0F_Ptmin) continue;

            /* Case: anti-lambda -> anti-proton + pi+ */

            if (pdgV0 == -3122) {

                /* Store V0 */

                kfAntiLambdas.push_back(kfV0);

                getEsdIdxOfNegDau_fromAntiLambdaIdx[kfAntiLambdas.size() - 1] = esdIdxNeg;
                getEsdIdxOfPosDau_fromAntiLambdaIdx[kfAntiLambdas.size() - 1] = esdIdxPos;

                isAntiLambdaIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                                 getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

                /* Fill histograms */

                fHist_Found_AntiLambda_Mass->Fill(lvV0.M());
                fHist_Found_AntiLambda_Radius->Fill(radius);
                fHist_Found_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_Found_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

                if (isAntiLambdaIdxATrueV0.back()) {

                    getMcIdx_fromAntiLambdaIdx[esdAntiLambdas.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                    mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                    if (mcV0->PdgCode() == pdgV0) {

                        fHist_Found_True_AntiLambda_Mass->Fill(lvV0.M());
                        fHist_Found_True_AntiLambda_Radius->Fill(radius);
                        fHist_Found_True_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_True_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);

                        if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                            fHist_Found_Signal_AntiLambda_Mass->Fill(lvV0.M());
                            fHist_Found_Signal_AntiLambda_Radius->Fill(radius);
                            fHist_Found_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
                            fHist_Found_Signal_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
                        }
                    }
                }
            }

            /* Case: K0S -> pi- + pi+ */

            if (pdgV0 == 310) {

                /* Store V0 */

                kfKaonsZeroShort.push_back(kfV0);

                getEsdIdxOfNegDau_fromKaonZeroShortIdx[kfKaonsZeroShort.size() - 1] = esdIdxNeg;
                getEsdIdxOfPosDau_fromKaonZeroShortIdx[kfKaonsZeroShort.size() - 1] = esdIdxPos;

                isKaonZeroShortIdxATrueV0.push_back(isEsdIdxDaughterOfTrueV0[esdIdxNeg] && isEsdIdxDaughterOfTrueV0[esdIdxPos] &&
                                                    getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg] == getMcIdxOfTrueV0_fromEsdIdx[esdIdxPos]);

                /* Fill histograms */

                fHist_Found_KaonZeroShort_Mass->Fill(lvV0.M());
                fHist_Found_KaonZeroShort_Radius->Fill(radius);
                fHist_Found_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                fHist_Found_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

                if (isKaonZeroShortIdxATrueV0.back()) {

                    getMcIdx_fromKaonZeroShortIdx[esdKaonsZeroShort.size() - 1] = getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg];

                    mcV0 = (AliMCParticle*)fMC->GetTrack(getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]);

                    if (mcV0->PdgCode() == pdgV0) {

                        fHist_Found_True_KaonZeroShort_Mass->Fill(lvV0.M());
                        fHist_Found_True_KaonZeroShort_Radius->Fill(radius);
                        fHist_Found_True_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                        fHist_Found_True_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);

                        if (isMcIdxSignal[getMcIdxOfTrueV0_fromEsdIdx[esdIdxNeg]]) {

                            fHist_Found_Signal_KaonZeroShort_Mass->Fill(lvV0.M());
                            fHist_Found_Signal_KaonZeroShort_Radius->Fill(radius);
                            fHist_Found_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
                            fHist_Found_Signal_KaonZeroShort_DCAwrtPV->Fill(dca_wrt_pv);
                        }
                    }
                }
            }
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*
 Apply cuts to an (anti-lambda) candidate, reconstructed from an (anti-proton, pi+) pair. Relevant for channel A and D.
 */
Bool_t AliAnalysisTaskSexaquark::PassesAntiLambdaCuts_KF(KFParticleMother kfV0, KFParticle kfNegDau, KFParticle kfPosDau, TLorentzVector lvV0,
                                                         TLorentzVector lvNegDau, TLorentzVector lvPosDau) {

    // define geometrical cuts
    // Float_t kfDCA_tracks = TMath::Abs(kfNegDau.GetDistanceFromParticleXY(kfPosDau));
    // Float_t kfDCA_neg_V0 =
    // TMath::Abs(kfNegDau.GetDistanceFromVertexXY(kfV0));  // = AliESDtrack::GetD(vertex_x, vertex_y, mag_field) = impar_neg[0]
    // Float_t kfDCA_pos_V0 = TMath::Abs(kfPosDau.GetDistanceFromVertexXY(kfV0));

    // define kinematic cuts
    /*
    Double_t armAlpha = ArmenterosAlpha(lvV0.Px(), lvV0.Py(), lvV0.Pz(),              //
                                        lvNegDau.Px(), lvNegDau.Py(), lvNegDau.Pz(),  //
                                        lvPosDau.Px(), lvPosDau.Py(), lvPosDau.Pz());
    Double_t armPt = ArmenterosQt(lvV0.Px(), lvV0.Py(), lvV0.Pz(), kfNegDau.Px(), kfNegDau.Py(), kfNegDau.Pz());
    Float_t mass = lvV0.M();
    Double_t cpa = CosinePointingAngle(lvV0, kfV0.GetX(), kfV0.GetY(), kfV0.GetZ(),  //
                                       fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());

    // apply them
    if ((armPt / TMath::Abs(armAlpha)) > kArmLambdaSlope) return kFALSE;
    if (mass > 1.13 || mass < 1.10) return kFALSE;  // (pending) refactorize
    if (cpa < 0.9) return kFALSE;                   // (pending)
    if (kfDCA_tracks > kMaxDCATracks) return kFALSE;
    */

    Double_t COV0F_DNmin = 0.1;            // (d=0.1) min imp parameter for the negative daughter (depends on PV!)
    Double_t COV0F_DPmin = 0.1;            // (d=0.1) min imp parameter for the positive daughter (depends on PV!)
    Double_t COV0F_Chi2max = 33.;          // (d=33.) max chi2
    Double_t COV0F_DCAmax = 0.6;           // (d=1.) max DCA between the daughter tracks
    Double_t COV0F_CPAmin = 0.4;           // (d=0.998) max cosine of V0's pointing angle (depends on PV!)
    Double_t COV0F_Rmin = 50.;             // (d=0.9) min radius of the V0 radius
    Double_t COV0F_Rmax = 200.;            // (d=100) max radius of the V0 radius
    Double_t COV0F_Etamax_V0 = 1.0;        // (d=1.) max eta of V0
    Double_t COV0F_DCAmin_V0_wrtPV = 5.;   // min DCA between V0 and PV
    Double_t COV0F_Etamax_NegDau = 1.5;    // max eta of negative daughter
    Double_t COV0F_Etamax_PosDau = 1.5;    // max eta of positive daughter
    Int_t COV0F_NClustersTPCmin_Dau = 50;  // minimun number of clusters in the TPC
    Double_t COV0F_Ptmin = 1.0;            // min Pt of V0

    // PENDING
    return kTRUE;
}

/*
 Apply cuts to a (K0S) candidate, reconstructed from a (pi-, pi+) pair. Relevant for channel A.
 */
Bool_t AliAnalysisTaskSexaquark::PassesKaonZeroShortCuts_KF(KFParticleMother kfV0, KFParticle kfNegDau, KFParticle kfPosDau,  //
                                                            TLorentzVector lvV0, TLorentzVector lvNegDau, TLorentzVector lvPosDau) {
    /* PENDING */
    return kTRUE;
}

/*
 Apply cuts to a (pi-, pi+) pair. Relevant for channel E.
 */
Bool_t AliAnalysisTaskSexaquark::PassesPionPairCuts_KF(KFParticleMother kfV0, KFParticle kfNegDau, KFParticle kfPosDau,  //
                                                       TLorentzVector lvV0, TLorentzVector lvNegDau, TLorentzVector lvPosDau) {
    /* PENDING */
    return kTRUE;
}

/*
 Apply cuts to a (pi-, K+) pair. Relevant for channel E.
*/
Bool_t AliAnalysisTaskSexaquark::PassesKaonPionPairCuts_KF(KFParticleMother kfV0, KFParticle kfDaughter1, KFParticle kfDaughter2,  //
                                                           TLorentzVector lvV0, TLorentzVector lvTrack1, TLorentzVector lvTrack2) {
    /* PENDING */
    return kTRUE;
}

/*
 Apply cuts to a (K+, K+) pair. Relevant for channel H.
 - Input: `kfV0`, `kfDaughter1`, `kfDaughter2`, `lvV0`, `lvTrack1`, `lvTrack2`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesPosKaonPairCuts_KF(KFParticleMother kfV0, KFParticle kfDaughter1, KFParticle kfDaughter2,  //
                                                          TLorentzVector lvV0, TLorentzVector lvTrack1, TLorentzVector lvTrack2) {
    /* PENDING */
    return kTRUE;
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel A.
 - Containers: `kf[AntiLambdas|KaonsZeroShort]`, `getEsdIdxOf[Neg|Pos]Dau_from[AntiLambda|KaonZeroShort]Idx`
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelA_KF() {

    Double_t impar[2];
    Bool_t fStatusT0PropToDCA;
    Bool_t fStatusT1PropToDCA;
    Bool_t fStatusT2PropToDCA;
    Bool_t fStatusT3PropToDCA;

    Int_t esdIdxNeg_AntiLambda, esdIdxPos_AntiLambda;
    Int_t esdIdxNeg_KaonZeroShort, esdIdxPos_KaonZeroShort;

    TLorentzVector lvNucleon(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());  // assume struck nucleon at rest

    AliESDVertex* esdVertex_AntiLambda;
    AliESDVertex* esdVertex_KaonZeroShort;

    AliESDtrack *esdTrack0, *esdTrack1;
    AliESDtrack *esdTrack2, *esdTrack3;

    TLorentzVector lvTrack0, lvTrack1;
    TLorentzVector lvTrack2, lvTrack3;

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;
    TLorentzVector lvAntiSexaquark;

    Float_t radius;
    Float_t cpa_wrt_pv;
    Float_t dca_wrt_pv;
    Bool_t is_signal;

    Int_t N_Comb = 0;
    Int_t N_PossibleComb = (Int_t)kfAntiLambdas.size() * (Int_t)kfKaonsZeroShort.size();

    AliInfoF("N_PossibleComb = %i", N_PossibleComb);

    /* Loop over all possible pairs of anti-lambdas and K0S */

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)kfAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)kfKaonsZeroShort.size(); idxKaonZeroShort++) {

            N_Comb++;
            if (N_Comb % (Int_t)TMath::Floor((Float_t)N_PossibleComb / 20.) == 0) AliInfoF("N_Comb = %i", N_Comb);

            /* Check that no daughter is repeated between the V0s */

            esdIdxNeg_AntiLambda = getEsdIdxOfNegDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxPos_AntiLambda = getEsdIdxOfPosDau_fromAntiLambdaIdx[idxAntiLambda];
            esdIdxNeg_KaonZeroShort = getEsdIdxOfNegDau_fromKaonZeroShortIdx[idxKaonZeroShort];
            esdIdxPos_KaonZeroShort = getEsdIdxOfPosDau_fromKaonZeroShortIdx[idxKaonZeroShort];

            if (esdIdxNeg_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxNeg_AntiLambda == esdIdxPos_KaonZeroShort ||
                esdIdxPos_AntiLambda == esdIdxNeg_KaonZeroShort || esdIdxPos_AntiLambda == esdIdxPos_KaonZeroShort) {
                continue;
            }

            /* Kalman Filter */

            KFParticleMother AntiLambda = kfAntiLambdas[idxAntiLambda];
            KFParticleMother KaonZeroShort = kfKaonsZeroShort[idxKaonZeroShort];

            Double_t vertex_AntiLambda[3] = {AntiLambda.GetX(), AntiLambda.GetY(), AntiLambda.GetZ()};
            Double_t vertexErr_AntiLambda[3] = {AntiLambda.GetErrX(), AntiLambda.GetErrY(), AntiLambda.GetErrZ()};
            esdVertex_AntiLambda = new AliESDVertex(vertex_AntiLambda, vertexErr_AntiLambda);

            Double_t vertex_KaonZeroShort[3] = {KaonZeroShort.GetX(), KaonZeroShort.GetY(), KaonZeroShort.GetZ()};
            Double_t vertexErr_KaonZeroShort[3] = {KaonZeroShort.GetErrX(), KaonZeroShort.GetErrY(), KaonZeroShort.GetErrZ()};
            esdVertex_KaonZeroShort = new AliESDVertex(vertex_KaonZeroShort, vertexErr_KaonZeroShort);

            /*
                AliInfoF("anti_lambda: vertex, radius, charge = (%f, %f, %f), (%f, %f, %f), %f, %i", AntiLambda.Px(), AntiLambda.Py(),
                AntiLambda.Pz(), AntiLambda.GetX(), AntiLambda.GetY(), AntiLambda.GetZ(), al_radius,
                (Int_t)AntiLambda.GetQ()); AliInfoF("kaon_zero_short: vertex, radius, charge = (%f, %f, %f), (%f, %f, %f), %f, %i",
               KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz(), KaonZeroShort.GetX(),
               KaonZeroShort.GetY(), KaonZeroShort.GetZ(), k0s_radius, (Int_t)KaonZeroShort.GetQ());
             */

            /* Add protection against ??? */

            // KFParticleMother CheckALAfterEmpty;

            // KFParticleMother CheckK0SAfterAL;
            // CheckK0SAfterAL.AddDaughter(AntiLambda);

            KFParticleMother AntiSexaquark;
            // if (CheckALAfterEmpty.CheckDaughter(AntiLambda) && CheckK0SAfterAL.CheckDaughter(KaonZeroShort)) {
            // AntiSexaquark.SetConstructMethod(2);
            // AntiSexaquark.AddDaughter(AntiLambda);
            // AntiSexaquark.AddDaughter(KaonZeroShort);
            // } else {
            AntiSexaquark.AddDaughter(AntiLambda);
            AntiSexaquark.AddDaughter(KaonZeroShort);
            // }
            AntiSexaquark.TransportToDecayVertex();

            /*
            // get Secondary Vertex
            Double_t SecVertex[3] = {AntiSexaquark.GetX(), AntiSexaquark.GetY(), AntiSexaquark.GetZ()};
            Double_t SecVertexErr[3] = {AntiSexaquark.GetErrX(), AntiSexaquark.GetErrY(), AntiSexaquark.GetErrZ()};
            AliESDVertex* esdSecVertex = new AliESDVertex(SecVertex, SecVertexErr);
            */

            /* Transport V0s to the secondary vertex */

            AntiLambda.SetProductionVertex(AntiSexaquark);
            KaonZeroShort.SetProductionVertex(AntiSexaquark);

            AntiLambda.TransportToProductionVertex();
            KaonZeroShort.TransportToProductionVertex();

            /* Get, clone and propagate the final-state AliESDtracks */

            esdTrack0 = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_AntiLambda));
            esdTrack1 = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_AntiLambda));
            esdTrack2 = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxNeg_KaonZeroShort));
            esdTrack3 = static_cast<AliESDtrack*>(fESD->GetTrack(esdIdxPos_KaonZeroShort));

            AliESDtrack cloneTrack0(*esdTrack0);
            AliESDtrack cloneTrack1(*esdTrack1);
            AliESDtrack cloneTrack2(*esdTrack2);
            AliESDtrack cloneTrack3(*esdTrack3);

            cloneTrack0.PropagateToDCA(esdVertex_AntiLambda, fMagneticField, 10., impar);
            cloneTrack1.PropagateToDCA(esdVertex_AntiLambda, fMagneticField, 10., impar);
            cloneTrack2.PropagateToDCA(esdVertex_KaonZeroShort, fMagneticField, 10., impar);
            cloneTrack3.PropagateToDCA(esdVertex_KaonZeroShort, fMagneticField, 10., impar);

            lvTrack0.SetPxPyPzE(cloneTrack0.Px(), cloneTrack0.Py(), cloneTrack0.Pz(),
                                TMath::Sqrt(cloneTrack0.P() * cloneTrack0.P() + fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass() *
                                                                                    fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass()));
            lvTrack1.SetPxPyPzE(cloneTrack1.Px(), cloneTrack1.Py(), cloneTrack1.Pz(),
                                TMath::Sqrt(cloneTrack1.P() * cloneTrack1.P() + fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass() *
                                                                                    fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass()));
            lvTrack2.SetPxPyPzE(cloneTrack2.Px(), cloneTrack2.Py(), cloneTrack2.Pz(),
                                TMath::Sqrt(cloneTrack2.P() * cloneTrack2.P() + fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass() *
                                                                                    fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass()));
            lvTrack3.SetPxPyPzE(cloneTrack3.Px(), cloneTrack3.Py(), cloneTrack3.Pz(),
                                TMath::Sqrt(cloneTrack3.P() * cloneTrack3.P() + fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass() *
                                                                                    fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass()));

            /* Reconstruct AntiSexaquark */

            lvAntiLambda = lvTrack0 + lvTrack1;
            lvKaonZeroShort = lvTrack2 + lvTrack3;

            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvNucleon;

            /* Calculate variables */

            radius = TMath::Sqrt(AntiSexaquark.GetX() * AntiSexaquark.GetX() + AntiSexaquark.GetY() * AntiSexaquark.GetY());
            cpa_wrt_pv = CosinePointingAngle(lvAntiSexaquark, AntiSexaquark.GetX(), AntiSexaquark.GetY(), AntiSexaquark.GetZ(),
                                             fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            dca_wrt_pv = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), AntiSexaquark.GetX(), AntiSexaquark.GetY(),
                                      AntiSexaquark.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            is_signal = isEsdIdxSignal[esdIdxNeg_AntiLambda] && isEsdIdxSignal[esdIdxPos_AntiLambda]  //
                        && isEsdIdxSignal[esdIdxNeg_KaonZeroShort] && isEsdIdxSignal[esdIdxPos_KaonZeroShort];
            is_signal = is_signal && getReactionIdx_fromEsdIdx[esdIdxNeg_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] &&
                        getReactionIdx_fromEsdIdx[esdIdxPos_AntiLambda] == getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] &&
                        getReactionIdx_fromEsdIdx[esdIdxNeg_KaonZeroShort] == getReactionIdx_fromEsdIdx[esdIdxPos_KaonZeroShort];

            /* Apply cuts */

            if (!PassesSexaquarkCuts_ChannelA_KF(AntiSexaquark, AntiLambda, KaonZeroShort, lvAntiSexaquark, lvAntiLambda, lvKaonZeroShort)) {
                continue;
            }

            /* Fill histograms */

            fHist_Found_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
            fHist_Found_AntiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
            fHist_Found_AntiSexaquark_DCAbtwDau->Fill(AntiLambda.GetDistanceFromVertexXY(KaonZeroShort));
            fHist_Found_AntiSexaquark_DCAwrtPV->Fill(dca_wrt_pv);
            fHist_Found_AntiSexaquark_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Found_AntiSexaquark_Radius->Fill(radius);

            if (!is_signal) continue;

            fHist_Found_Signal_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
            fHist_Found_Signal_AntiSexaquark_Pt->Fill(lvAntiSexaquark.Pt());
            fHist_Found_Signal_AntiSexaquark_DCAbtwDau->Fill(AntiLambda.GetDistanceFromVertexXY(KaonZeroShort));
            fHist_Found_Signal_AntiSexaquark_DCAwrtPV->Fill(dca_wrt_pv);
            fHist_Found_Signal_AntiSexaquark_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_Found_Signal_AntiSexaquark_Radius->Fill(radius);
        }  // end of loop over K0S
    }      // end of loop over anti-lambdas
}

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel D.
 - Uses: `fStruckNucleonPDG`, `fFinalStateProductsPDG`
 - Input: `kfAntiLambdas`, `idxAntiLambdaDaughters`, `idxPositiveKaons`
 - Output: `kfAntiSexaquarks`
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelD_KF(std::vector<KFParticleMother> kfAntiLambdas,                  //
                                                           std::vector<std::pair<Int_t, Int_t>> idxAntiLambdaDaughters,  //
                                                           std::vector<Int_t> idxPositiveKaons,                          //
                                                           std::vector<KFParticleMother>& kfAntiSexaquarks) {
    /* PENDING */
}

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel E.
 This one is trickier because it needs three components:
 (1) the rec. anti-lambda (anti-proton + pos. pion),
 (2) a pair of charged particles (neg. pion + pos. pion, or neg. pion + pos. kaon)
 (3) the remaining single pos. particle (pos. kaon or pos. pion)
 - Uses: `fStruckNucleonPDG`
 - Input: `kfAntiLambdas`, `idxAntiLambdaDaughters`, `idxPositiveKaons`
 - Output: `kfAntiSexaquarks`
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelE_KF(std::vector<KFParticleMother> kfAntiLambdas,                  //
                                                           std::vector<std::pair<Int_t, Int_t>> idxAntiLambdaDaughters,  //
                                                           std::vector<KFParticleMother> kfChargedPairs,                 //
                                                           std::vector<std::pair<Int_t, Int_t>> idxChargedPairsTracks,   //
                                                           std::vector<Int_t> idxSinglePositiveTrack,                    //
                                                           std::vector<KFParticleMother>& kfAntiSexaquarks) {
    /* PENDING */
}

/*
 Apply anti-sexaquark candidate cuts.
 - Input: `kf[AntiSexaquark|AntiLambda|KaonZeroShort]`, `lv[AntiSexaquark|AntiLambda|KaonZeroShort]`
 - Uses: `fPrimaryVertex`
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA_KF(KFParticleMother kfAntiSexaquark, KFParticle kfAntiLambda,
                                                                 KFParticle kfKaonZeroShort, TLorentzVector lvAntiSexaquark,
                                                                 TLorentzVector lvAntiLambda, TLorentzVector lvKaonZeroShort) {

    // >> DCA w.r.t. PV
    Double_t dca = LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(), kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                                kfAntiSexaquark.GetZ(), fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    if (dca > 0.1) {
        AliInfo("Found Sexaquark that had DCA w.r.t. PV > 0.1");
        return kFALSE;
    }

    return kTRUE;
}

/*
 Apply anti-sexaquark candidate cuts.
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelD_KF() {
    /* PENDING */
    return kTRUE;
}

/*
 Apply anti-sexaquark candidate cuts.
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelE_KF() {
    /* PENDING */
    return kTRUE;
}

/*                            */
/**  Mathematical Functions  **/
/*** ====================== ***/

/*
 Calculate the cosine of the pointing angle of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 - Input: `lvParticle`, `X`, `Y`, `Z`, `refPointX`, `refPointY`, `refPointZ`
 - Return: `cos(theta)`
 (Based on `AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()`)
*/
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(TLorentzVector lvParticle, Double_t X, Double_t Y, Double_t Z,  //
                                                       Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    return TMath::Cos(lvParticle.Angle(posRelativeToRef));
}

/*
 Overloaded version of `CosinePointingAngle(...)`, using Px, Py, Pz instead of the particle's TLorentzVector.
*/
Double_t AliAnalysisTaskSexaquark::CosinePointingAngle(Double_t Px, Double_t Py, Double_t Pz, Double_t X, Double_t Y, Double_t Z,  //
                                                       Double_t refPointX, Double_t refPointY, Double_t refPointZ) {
    TVector3 posRelativeToRef(X - refPointX, Y - refPointY, Z - refPointZ);
    TVector3 momParticle(Px, Py, Pz);
    return TMath::Cos(momParticle.Angle(posRelativeToRef));
}

/*
 This gives the Armenteros-Podolanski alpha.
 (Based on `AliRoot/STEER/ESD/AliESDv0::AlphaV0()`)
*/
Double_t AliAnalysisTaskSexaquark::Calculate_ArmAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
                                                      Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz,  //
                                                      Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momPos(Pos_Px, Pos_Py, Pos_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

    Double_t lQlNeg = momNeg.Dot(momTot) / momTot.Mag();
    Double_t lQlPos = momPos.Dot(momTot) / momTot.Mag();

    // (protection)
    if (lQlPos + lQlNeg == 0.) {
        return 2;
    }  // closure
    return (lQlPos - lQlNeg) / (lQlPos + lQlNeg);
}

/*
 This gives the Armenteros-Podolanski qT
 (Based on `AliRoot/STEER/ESD/AliESDv0::PtArmV0()`)
*/
Double_t AliAnalysisTaskSexaquark::Calculate_ArmPt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                   Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz) {
    TVector3 momTot(V0_Px, V0_Py, V0_Pz);
    TVector3 momNeg(Neg_Px, Neg_Py, Neg_Pz);

    return momNeg.Perp(momTot);
}

/*
 Find the distance of closest approach to the Primary Vertex, after backtracking a V0.
 This function stores the point of closest approach, and returns its distance to the PV.
 (PENDING: are we sure of that?)
 - Input: `V0_Px`, `V0_Py`, `V0_Pz`, `V0_X`, `V0_Y`, `V0_Z`, `PV_X`, `PV_Y`, `PV_Z`
 - Return: the distance of closest approach to the PV
*/
Double_t AliAnalysisTaskSexaquark::LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                Double_t refPointX, Double_t refPointY, Double_t refPointZ) {

    TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
    TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
    TVector3 RefVertex(refPointX, refPointY, refPointZ);

    TVector3 CrossProduct = (RefVertex - V0Vertex).Cross(V0Momentum);

    return CrossProduct.Mag() / V0Momentum.Mag();

    /*
    // convert input vectors to arrays
    Double_t ref[3] = {refPointX, refPointY, refPointZ};
    Double_t pos0[3] = {V0_X, V0_Y, V0_Z};
    Double_t dir0[3] = {V0_Px, V0_Py, V0_Pz};

    // lambda function
    auto func = [this, &ref, &pos0, &dir0](const Double_t* t) { return SquaredDistancePointToLine(t, ref, pos0, dir0); };
    ROOT::Math::Functor f(func, 1);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", 0., 0.01);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();

    return TMath::Sqrt(SquaredDistancePointToLine(xs, ref, pos0, dir0));
    */
}

/*
 Calculate transverse impact parameter w.r.t. Primary Vertex.
*/
Float_t AliAnalysisTaskSexaquark::MCRec_GetImpactParameter(AliESDtrack* track) {
    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

    return (Float_t)track->GetD(PV[0], PV[1], fESD->GetMagneticField());
}

/*
 Find a common vertex for two neutral particles, given their position and momentum, assuming they propagate in straight lines.
 If a common vertex cannot be found, then return the DCA between the two lines.
 (Based on Marty Cohen's answer in https://math.stackexchange.com/questions/2213165/find-shortest-distance-between-lines-in-3d)
*/
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v1(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    if (v3_dir0.Mag() == 0. || v3_dir1.Mag() == 0.) {
        AliWarning("One of the directions is null!");
        return -1.;
    }

    /* Require perpendicularity, that means both lines must intersect eventually */

    if (v3_dir0.Cross(v3_dir1).Mag() <= 0.) {
        AliWarning("TwoLinesDCA :: Not possible to intersect!");
        return -1.;
    }

    TVector3 diff = v3_pos0 - v3_pos1;
    Double_t determinant = v3_dir0.Dot(v3_dir1) * v3_dir0.Dot(v3_dir1) - v3_dir0.Mag2() * v3_dir1.Mag2();
    if (determinant == 0.) {
        AliWarning("TwoLinesDCA :: Lines are parallel!");
        return -1.;
    }
    Double_t t0 = (v3_dir1.Mag2() * v3_dir0.Dot(diff) - v3_dir1.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;
    Double_t t1 = (-1 * v3_dir0.Mag2() * v3_dir1.Dot(diff) + v3_dir0.Dot(diff) * v3_dir1.Dot(v3_dir0)) / determinant;

    TString aux = "";
    if (t0 >= 0. || t1 >= 0.) {
        aux = "*";
    }

    Double_t dca = (diff + v3_dir0 * t0 - v3_dir1 * t1).Mag();
    PCA0 = v3_pos0 + t0 * v3_dir0;
    PCA1 = v3_pos1 + t1 * v3_dir1;
    TVector3 middle_point((PCA0.X() + PCA1.X()) * 0.5, (PCA0.Y() + PCA1.Y()) * 0.5, (PCA0.Z() + PCA1.Z()) * 0.5);

    AliInfoF("dca pca1 pca2 radius = %f (%f, %f, %f) (%f, %f, %f) %f %s", dca, PCA0.X(), PCA0.Y(), PCA0.Z(), PCA1.X(), PCA1.Y(), PCA1.Z(),
             middle_point.Perp(), aux.Data());

    return dca;
}

/*
 Given the parametric form of a line, find the distance between it and a point.
 - Variable: `t`, the parameter of the line.
 - Parameters: `point`, the point to which we want to calculate the distance.
 - Parameters: `linePoint`, `lineDir`, the position and direction of the line.
 - Return: the square of the distance between the point and the line.
*/
Double_t AliAnalysisTaskSexaquark::SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]) {
    return TMath::Power(point[0] - linePoint[0] - t[0] * lineDir[0], 2) +  //
           TMath::Power(point[1] - linePoint[1] - t[0] * lineDir[1], 2) +  //
           TMath::Power(point[2] - linePoint[2] - t[0] * lineDir[2], 2);
}

/*
 Given the parametric form of two lines, find their distance.
 - Variable: `t`, array of size 2, where `t[0]` and `t[1] are the parameters for the first and second line, respectively.
 - Parameters: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 - Return: the square of the distance between the two lines.
*/
Double_t AliAnalysisTaskSexaquark::SquaredDistanceBetweenLines(const Double_t* t, Double_t pos0[], Double_t dir0[], Double_t pos1[],
                                                               Double_t dir1[]) {
    return TMath::Power(dir1[0] * t[1] + pos1[0] - (dir0[0] * t[0] + pos0[0]), 2) +
           TMath::Power(dir1[1] * t[1] + pos1[1] - (dir0[1] * t[0] + pos0[1]), 2) +
           TMath::Power(dir1[2] * t[1] + pos1[2] - (dir0[2] * t[0] + pos0[2]), 2);
}

/*
 Calculate the distance of closest approach between two lines, given their position and momentum, assuming they propagate in straight lines.
 - Input: `pos0`, `dir0`, `pos1`, `dir1`, the position and direction of the two lines.
 - Output: `CPA0`, `CPA1`, the closest point of approach for the two lines.
 - Return: the distance of closest approach between the two lines.
 */
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& PCA0,
                                                            TVector3& PCA1) {

    // convert input vectors to arrays
    Double_t pos0[3] = {v3_pos0.X(), v3_pos0.Y(), v3_pos0.Z()};
    Double_t dir0[3] = {v3_dir0.X(), v3_dir0.Y(), v3_dir0.Z()};

    Double_t pos1[3] = {v3_pos1.X(), v3_pos1.Y(), v3_pos1.Z()};
    Double_t dir1[3] = {v3_dir1.X(), v3_dir1.Y(), v3_dir1.Z()};

    // lambda function
    auto func = [this, &pos0, &dir0, &pos1, &dir1](const Double_t* t) { return SquaredDistanceBetweenLines(t, pos0, dir0, pos1, dir1); };
    ROOT::Math::Functor f(func, 2);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetLimitedVariable(0, "t0", 0.0, 0.1, -999., 0.);
    min->SetLimitedVariable(1, "t1", 0.0, 0.1, -999., 0.);
    min->Minimize();

    // print results
    const Double_t* xs = min->X();
    Double_t t0 = xs[0];
    Double_t t1 = xs[1];

    Double_t dca = TMath::Sqrt(SquaredDistanceBetweenLines(xs, pos0, dir0, pos1, dir1));
    PCA0.SetXYZ(dir0[0] * t0 + pos0[0], dir0[1] * t0 + pos0[1], dir0[2] * t0 + pos0[2]);
    PCA1.SetXYZ(dir1[0] * t1 + pos1[0], dir1[1] * t1 + pos1[1], dir1[2] * t1 + pos1[2]);
    TVector3 middle_point((PCA0.X() + PCA1.X()) * 0.5, (PCA0.Y() + PCA1.Y()) * 0.5, (PCA0.Z() + PCA1.Z()) * 0.5);

    AliInfoF("dca pca1 pca2 radius = %f (%f, %f, %f) (%f, %f, %f) %f", dca, PCA0.X(), PCA0.Y(), PCA0.Z(), PCA1.X(), PCA1.Y(), PCA1.Z(),
             middle_point.Perp());

    return dca;
}

/*** FUNCTIONS FOR THE CUSTOM OFFLINE V0 FINDER ***/

/*
 This function pre-optimizes a two-track pair in the XY plane
 and provides two X values for the tracks if successful
 [copied from AliRoot/STEER/ESD/AliV0vertexer.cxx]
*/
Bool_t AliAnalysisTaskSexaquark::Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn,
                                             Double_t* lPreprocessxp, const Double_t b) {
    Double_t lMinimumX = -3.0;
    Double_t lMaximumX = 300.0;

    Double_t nhelix[6], phelix[6];
    nt->GetHelixParameters(nhelix, b);
    pt->GetHelixParameters(phelix, b);
    Double_t lNegCenterR[2], lPosCenterR[2];

    // Negative track parameters in XY
    GetHelixCenter(nt, lNegCenterR, b);
    Double_t xNegCenter = lNegCenterR[0];
    Double_t yNegCenter = lNegCenterR[1];
    Double_t NegRadius = TMath::Abs(1. / nhelix[4]);

    // Positive track parameters in XY
    GetHelixCenter(pt, lPosCenterR, b);
    Double_t xPosCenter = lPosCenterR[0];
    Double_t yPosCenter = lPosCenterR[1];
    Double_t PosRadius = TMath::Abs(1. / phelix[4]);

    // Define convenient coordinate system
    // Logical zero: position of negative center
    Double_t ux = xPosCenter - xNegCenter;
    Double_t uy = yPosCenter - yNegCenter;

    // Check center-to-center distance
    Double_t lDist = TMath::Sqrt(TMath::Power(xNegCenter - xPosCenter, 2) + TMath::Power(yNegCenter - yPosCenter, 2));
    // Normalize ux, uz to unit vector
    ux /= lDist;
    uy /= lDist;

    // Calculate perpendicular vector (normalized)
    Double_t vx = -uy;
    Double_t vy = +ux;

    Double_t lPreprocessDCAxy = 1e+3;  // define outside scope
    *lPreprocessxp = pt->GetX();       // start at current location
    *lPreprocessxn = nt->GetX();       // start at current location

    //============================================================
    // Pre-optimization in the XY plane: cases considered here
    //============================================================
    //
    //  Case 1: Circles do not touch, centers far away
    //          (D > R1 + R2)
    //
    //  Case 2: Circles touch, centers at reasonable distance wrt D
    //          (D < R1 + R2) && (D > |R1-R2|)
    //
    //  Case 3: Circles do not touch, one inside the other
    //          (D < |R1-R2|)
    //
    //  Cases 1 and 2 are treated. Case 3 is not treated (unlikely
    //  to be a problem with unlike-sign charged tracks): brute
    //  force minimization takes place in any case
    //
    //============================================================

    //______________________
    // CASE 1
    if ((lDist > NegRadius + PosRadius)) {
        //================================================================
        // Case 1: distance bigger than sum of radii ("gamma-like")
        //        re-position tracks along the center-to-center axis
        // Re-position negative track
        Double_t xNegOptPosition = xNegCenter + NegRadius * ux;
        Double_t yNegOptPosition = yNegCenter + NegRadius * uy;
        Double_t csNeg = TMath::Cos(nt->GetAlpha());
        Double_t snNeg = TMath::Sin(nt->GetAlpha());
        Double_t xThisNeg = xNegOptPosition * csNeg + yNegOptPosition * snNeg;

        // Re-position positive track
        Double_t xPosOptPosition = xPosCenter - PosRadius * ux;
        Double_t yPosOptPosition = yPosCenter - PosRadius * uy;
        Double_t csPos = TMath::Cos(pt->GetAlpha());
        Double_t snPos = TMath::Sin(pt->GetAlpha());
        Double_t xThisPos = xPosOptPosition * csPos + yPosOptPosition * snPos;

        if (xThisNeg < lMaximumX && xThisPos < lMaximumX && xThisNeg > lMinimumX && xThisPos > lMinimumX) {
            Double_t lCase1NegR[3] = {0.}, lCase1PosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg, b, lCase1NegR) && pt->GetXYZAt(xThisPos, b, lCase1PosR)) {
                lPreprocessDCAxy = TMath::Sqrt(TMath::Power(lCase1NegR[0] - lCase1PosR[0], 2) + TMath::Power(lCase1NegR[1] - lCase1PosR[1], 2) +
                                               TMath::Power(lCase1NegR[2] - lCase1PosR[2], 2));
                // Pass coordinates
                if (lPreprocessDCAxy < 999) {
                    *lPreprocessxp = xThisPos;
                    *lPreprocessxn = xThisNeg;
                }
            }
        }
        //================================================================
    }

    //______________________
    // CASE 2
    if ((lDist > TMath::Abs(NegRadius - PosRadius)) && (lDist < NegRadius + PosRadius)) {
        //================================================================
        // Case 2: distance smaller than sum of radii (cowboy/sailor configs)

        // Calculate coordinate for radical line
        Double_t lRadical = (lDist * lDist - PosRadius * PosRadius + NegRadius * NegRadius) / (2 * lDist);

        // Calculate absolute displacement from center-to-center axis
        Double_t lDisplace = (0.5 / lDist) * TMath::Sqrt((-lDist + PosRadius - NegRadius) * (-lDist - PosRadius + NegRadius) *
                                                         (-lDist + PosRadius + NegRadius) * (lDist + PosRadius + NegRadius));

        // 3D distances in the two cases studied (prefer smallest)
        Double_t lCase2aDCA = 1e+3;
        Double_t lCase2bDCA = 1e+3;

        // 2 cases: positive and negative displacement
        Double_t xNegOptPosition[2], yNegOptPosition[2], xPosOptPosition[2], yPosOptPosition[2];
        Double_t csNeg, snNeg, csPos, snPos;
        Double_t xThisNeg[2], xThisPos[2];

        csNeg = TMath::Cos(nt->GetAlpha());
        snNeg = TMath::Sin(nt->GetAlpha());
        csPos = TMath::Cos(pt->GetAlpha());
        snPos = TMath::Sin(pt->GetAlpha());

        // Case 2a: Positive displacement along v vector
        // Re-position negative track
        xNegOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yNegOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisNeg[0] = xNegOptPosition[0] * csNeg + yNegOptPosition[0] * snNeg;
        // Re-position positive track
        xPosOptPosition[0] = xNegCenter + lRadical * ux + lDisplace * vx;
        yPosOptPosition[0] = yNegCenter + lRadical * uy + lDisplace * vy;
        xThisPos[0] = xPosOptPosition[0] * csPos + yPosOptPosition[0] * snPos;

        // Case 2b: Negative displacement along v vector
        // Re-position negative track
        xNegOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yNegOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisNeg[1] = xNegOptPosition[1] * csNeg + yNegOptPosition[1] * snNeg;
        // Re-position positive track
        xPosOptPosition[1] = xNegCenter + lRadical * ux - lDisplace * vx;
        yPosOptPosition[1] = yNegCenter + lRadical * uy - lDisplace * vy;
        xThisPos[1] = xPosOptPosition[1] * csPos + yPosOptPosition[1] * snPos;

        // Test the two cases, please

        // Case 2a
        if (xThisNeg[0] < lMaximumX && xThisPos[0] < lMaximumX && xThisNeg[0] > lMinimumX && xThisPos[0] > lMinimumX) {
            Double_t lCase2aNegR[3] = {0.}, lCase2aPosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg[0], b, lCase2aNegR) && pt->GetXYZAt(xThisPos[0], b, lCase2aPosR)) {
                lCase2aDCA = TMath::Sqrt(TMath::Power(lCase2aNegR[0] - lCase2aPosR[0], 2) + TMath::Power(lCase2aNegR[1] - lCase2aPosR[1], 2) +
                                         TMath::Power(lCase2aNegR[2] - lCase2aPosR[2], 2));
            }
        }

        // Case 2b
        if (xThisNeg[1] < lMaximumX && xThisPos[1] < lMaximumX && xThisNeg[1] > lMinimumX && xThisPos[1] > lMinimumX) {
            Double_t lCase2bNegR[3] = {0.}, lCase2bPosR[3] = {0.};
            if (nt->GetXYZAt(xThisNeg[1], b, lCase2bNegR) && pt->GetXYZAt(xThisPos[1], b, lCase2bPosR)) {
                lCase2bDCA = TMath::Sqrt(TMath::Power(lCase2bNegR[0] - lCase2bPosR[0], 2) + TMath::Power(lCase2bNegR[1] - lCase2bPosR[1], 2) +
                                         TMath::Power(lCase2bNegR[2] - lCase2bPosR[2], 2));
            }
        }
        // Minor detail: all things being equal, prefer closest X
        Double_t lCase2aSumX = xThisPos[0] + xThisNeg[0];
        Double_t lCase2bSumX = xThisPos[1] + xThisNeg[1];

        Double_t lDCAxySmallestR = lCase2aDCA;
        Double_t lxpSmallestR = xThisPos[0];
        Double_t lxnSmallestR = xThisNeg[0];

        Double_t lDCAxyLargestR = lCase2bDCA;
        Double_t lxpLargestR = xThisPos[1];
        Double_t lxnLargestR = xThisNeg[1];

        if (lCase2bSumX + 1e-6 < lCase2aSumX) {
            lDCAxySmallestR = lCase2bDCA;
            lxpSmallestR = xThisPos[1];
            lxnSmallestR = xThisNeg[1];
            lDCAxyLargestR = lCase2aDCA;
            lxpLargestR = xThisPos[0];
            lxnLargestR = xThisNeg[0];
        }

        // Pass conclusion to lPreprocess variables, please
        lPreprocessDCAxy = lDCAxySmallestR;
        *lPreprocessxp = lxpSmallestR;
        *lPreprocessxn = lxnSmallestR;
        if (lDCAxyLargestR + 1e-6 < lDCAxySmallestR) {  // beware epsilon: numerical calculations are unstable here
            lPreprocessDCAxy = lDCAxyLargestR;
            *lPreprocessxp = lxpLargestR;
            *lPreprocessxn = lxnLargestR;
        }
        // Protection against something too crazy, please
        if (lPreprocessDCAxy > 999) {
            *lPreprocessxp = pt->GetX();  // start at current location
            *lPreprocessxn = nt->GetX();  // start at current location
        }
    }
    // End of preprocessing stage!
    // at this point lPreprocessxp, lPreprocessxn are already good starting points: update helixparams
    Bool_t lWorked = kFALSE;
    if (lPreprocessDCAxy < 999) {
        lWorked = kTRUE;
    }
    return lWorked;
}

/*
 Copied from AliV0ReaderV1::GetHelixCenter
 Get Center of the helix track parametrization
 (based on AliRoot/STEER/ESD/AliV0vertexer.cxx)
*/
void AliAnalysisTaskSexaquark::GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b) {

    Int_t charge = track->Charge();

    Double_t helix[6];
    track->GetHelixParameters(helix, b);

    Double_t xpos = helix[5];
    Double_t ypos = helix[0];
    Double_t radius = TMath::Abs(1. / helix[4]);
    Double_t phi = helix[2];
    if (phi < 0) {
        phi = phi + 2 * TMath::Pi();
    }

    phi -= TMath::Pi() / 2.;
    Double_t xpoint = radius * TMath::Cos(phi);
    Double_t ypoint = radius * TMath::Sin(phi);

    if (b < 0 && charge > 0) {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }

    if (b > 0 && charge < 0) {
        xpoint = -xpoint;
        ypoint = -ypoint;
    }

    center[0] = xpos + xpoint;
    center[1] = ypos + ypoint;
    return;
}

/*
 Clear all the containers.
 - Uses:
*/
void AliAnalysisTaskSexaquark::ClearContainers() {

    mcIndicesOfTrueV0s.clear();

    isMcIdxSignal.clear();
    isMcIdxDaughterOfSecondaryV0.clear();
    isMcIdxDaughterOfTrueV0.clear();
    isMcIdxDaughterOfSignalV0.clear();

    getMcIdxOfTrueV0_fromMcIdxOfDau.clear();
    getMcIdxOfNegDau_fromMcIdxOfTrueV0.clear();
    getMcIdxOfPosDau_fromMcIdxOfTrueV0.clear();

    getReactionIdx_fromMcIdx.clear();
    getMcIdx_fromReactionIdx.clear();

    getReactionIdx_fromEsdIdx.clear();
    getEsdIdx_fromReactionIdx.clear();

    esdIndicesOfAntiProtonTracks.clear();
    esdIndicesOfPosKaonTracks.clear();
    esdIndicesOfPiPlusTracks.clear();
    esdIndicesOfPiMinusTracks.clear();

    getMcIdx_fromEsdIdx.clear();

    isEsdIdxSignal.clear();
    isEsdIdxDaughterOfSecondaryV0.clear();
    isEsdIdxDaughterOfTrueV0.clear();

    getMcIdxOfTrueV0_fromEsdIdx.clear();

    getEsdIdxOfRecNegDau_fromMcIdxOfTrueV0.clear();
    getEsdIdxOfRecPosDau_fromMcIdxOfTrueV0.clear();

    esdAntiLambdas.clear();
    esdKaonsZeroShort.clear();

    getEsdIdxOfNegDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfPosDau_fromAntiLambdaIdx.clear();
    getEsdIdxOfNegDau_fromKaonZeroShortIdx.clear();
    getEsdIdxOfPosDau_fromKaonZeroShortIdx.clear();

    isAntiLambdaIdxATrueV0.clear();
    isKaonZeroShortIdxATrueV0.clear();

    getMcIdx_fromAntiLambdaIdx.clear();
    getMcIdx_fromKaonZeroShortIdx.clear();

    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
}

/*** TREE OPERATIONS ***/

/*
 Define fTree's branches, and load them into fEvent (which is an Event_tt object)
*/
void AliAnalysisTaskSexaquark::SetBranches() {
    /* MC particles */
    fTree->Branch("N_MCGen", &fEvent.N_MCGen);
    fTree->Branch("MC_Px", &fEvent.MC_Px);
    fTree->Branch("MC_Py", &fEvent.MC_Py);
    fTree->Branch("MC_Pz", &fEvent.MC_Pz);
    fTree->Branch("MC_X", &fEvent.MC_X);
    fTree->Branch("MC_Y", &fEvent.MC_Y);
    fTree->Branch("MC_Z", &fEvent.MC_Z);
    fTree->Branch("MC_Xf", &fEvent.MC_Xf);
    fTree->Branch("MC_Yf", &fEvent.MC_Yf);
    fTree->Branch("MC_Zf", &fEvent.MC_Zf);
    fTree->Branch("MC_PID", &fEvent.MC_PID);
    fTree->Branch("MC_Mother", &fEvent.MC_Mother);
    fTree->Branch("MC_PID_Mother", &fEvent.MC_PID_Mother);
    fTree->Branch("MC_PID_GrandMother", &fEvent.MC_PID_GrandMother);
    fTree->Branch("MC_NDaughters", &fEvent.MC_NDaughters);
    fTree->Branch("MC_FirstDau", &fEvent.MC_FirstDau);
    fTree->Branch("MC_LastDau", &fEvent.MC_LastDau);
    fTree->Branch("MC_Status", &fEvent.MC_Status);
    fTree->Branch("MC_isSignal", &fEvent.MC_isSignal);
    /* MC Rec. Tracks */
    fTree->Branch("N_MCRec", &fEvent.N_MCRec);
    fTree->Branch("Idx_True", &fEvent.Idx_True);
    fTree->Branch("Rec_Px", &fEvent.Rec_Px);
    fTree->Branch("Rec_Py", &fEvent.Rec_Py);
    fTree->Branch("Rec_Pz", &fEvent.Rec_Pz);
    fTree->Branch("Rec_Charge", &fEvent.Rec_Charge);
    fTree->Branch("Rec_IP_wrtPV", &fEvent.Rec_IP_wrtPV);
    fTree->Branch("Rec_NSigmaPion", &fEvent.Rec_NSigmaPion);
    fTree->Branch("Rec_NSigmaProton", &fEvent.Rec_NSigmaProton);
    fTree->Branch("Rec_NClustersTPC", &fEvent.Rec_NClustersTPC);
    fTree->Branch("Rec_isDuplicate", &fEvent.Rec_isDuplicate);
    fTree->Branch("Rec_isSimilar", &fEvent.Rec_isSimilar);
    fTree->Branch("Rec_isSignal", &fEvent.Rec_isSignal);
    fTree->Branch("Rec_HelixParam0", &fEvent.Rec_HelixParam0);
    fTree->Branch("Rec_HelixParam1", &fEvent.Rec_HelixParam1);
    fTree->Branch("Rec_HelixParam2", &fEvent.Rec_HelixParam2);
    fTree->Branch("Rec_HelixParam3", &fEvent.Rec_HelixParam3);
    fTree->Branch("Rec_HelixParam4", &fEvent.Rec_HelixParam4);
    fTree->Branch("Rec_HelixParam5", &fEvent.Rec_HelixParam5);
    /* V0s */
    fTree->Branch("N_V0s", &fEvent.N_V0s);
    fTree->Branch("Idx_Pos", &fEvent.Idx_Pos);
    fTree->Branch("Idx_Neg", &fEvent.Idx_Neg);
    fTree->Branch("V0_Px", &fEvent.V0_Px);
    fTree->Branch("V0_Py", &fEvent.V0_Py);
    fTree->Branch("V0_Pz", &fEvent.V0_Pz);
    fTree->Branch("V0_X", &fEvent.V0_X);
    fTree->Branch("V0_Y", &fEvent.V0_Y);
    fTree->Branch("V0_Z", &fEvent.V0_Z);
    fTree->Branch("Pos_Px", &fEvent.Pos_Px);
    fTree->Branch("Pos_Py", &fEvent.Pos_Py);
    fTree->Branch("Pos_Pz", &fEvent.Pos_Pz);
    fTree->Branch("Neg_Px", &fEvent.Neg_Px);
    fTree->Branch("Neg_Py", &fEvent.Neg_Py);
    fTree->Branch("Neg_Pz", &fEvent.Neg_Pz);
    fTree->Branch("V0_isSignal", &fEvent.V0_isSignal);
    fTree->Branch("V0_E_asK0", &fEvent.V0_E_asK0);
    fTree->Branch("V0_E_asAL", &fEvent.V0_E_asAL);
    fTree->Branch("V0_couldBeK0", &fEvent.V0_couldBeK0);
    fTree->Branch("V0_couldBeAL", &fEvent.V0_couldBeAL);
    fTree->Branch("V0_Chi2", &fEvent.V0_Chi2);
    fTree->Branch("V0_DCA_Daughters", &fEvent.V0_DCA_Daughters);
    fTree->Branch("V0_IP_wrtPV", &fEvent.V0_IP_wrtPV);
    fTree->Branch("V0_CPA_wrtPV", &fEvent.V0_CPA_wrtPV);
    fTree->Branch("V0_ArmAlpha", &fEvent.V0_ArmAlpha);
    fTree->Branch("V0_ArmPt", &fEvent.V0_ArmPt);
    fTree->Branch("V0_DecayLength", &fEvent.V0_DecayLength);
    fTree->Branch("V0_DCA_wrtPV", &fEvent.V0_DCA_wrtPV);
}

/*
 Push one of the MC generated particles into the Event_tt object
*/
void AliAnalysisTaskSexaquark::MCGen_PushBack(Int_t evt_mc) {

    Int_t idx_mc = 0;
    AliMCParticle* mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);

    Int_t evt_mother;
    Int_t mother_pid;
    Int_t grandmother_pid;
    // MCGen_GetMotherInfo(mcPart, evt_mother, mother_pid, grandmother_pid);

    Int_t n_daughters;
    Int_t evt_first_dau;
    Int_t evt_last_dau;
    // MCGen_GetDaughtersInfo(mcPart, n_daughters, evt_first_dau, evt_last_dau);

    Float_t x_f = -999.;
    Float_t y_f = -999.;
    Float_t z_f = -999.;
    if (n_daughters) {
        AliMCParticle* mcDaughter = (AliMCParticle*)fMC->GetTrack(mcPart->GetDaughterFirst());
        x_f = mcDaughter->Xv();
        y_f = mcDaughter->Yv();
        z_f = mcDaughter->Zv();
    }

    Bool_t is_signal = 0;  // MCGen_IsSignal(mcPart);

#if DEBUG_MODE == 1
    AliInfoF("%6i %6i %6i %6i %6i %6i %6i %6i %6i", evt_mc, idx_mc, mcPart->PdgCode(), is_signal,  //
             evt_mother, mother_pid, grandmother_pid,                                              //
             evt_first_dau, evt_last_dau);
#endif

    fEvent.MC_Px.push_back(mcPart->Px());
    fEvent.MC_Py.push_back(mcPart->Py());
    fEvent.MC_Pz.push_back(mcPart->Pz());
    fEvent.MC_X.push_back(mcPart->Xv());
    fEvent.MC_Y.push_back(mcPart->Yv());
    fEvent.MC_Z.push_back(mcPart->Zv());
    fEvent.MC_Xf.push_back(x_f);
    fEvent.MC_Yf.push_back(y_f);
    fEvent.MC_Zf.push_back(z_f);
    fEvent.MC_PID.push_back(mcPart->PdgCode());
    fEvent.MC_Mother.push_back(evt_mother);
    fEvent.MC_PID_Mother.push_back(mother_pid);
    fEvent.MC_PID_GrandMother.push_back(grandmother_pid);
    fEvent.MC_NDaughters.push_back(n_daughters);
    fEvent.MC_FirstDau.push_back(evt_first_dau);
    fEvent.MC_LastDau.push_back(evt_last_dau);
    fEvent.MC_isSignal.push_back(is_signal);
    fEvent.MC_Status.push_back(mcPart->MCStatusCode());
}

/*
 Push one of the MC reconstructed tracks into the Event_tt object
*/
void AliAnalysisTaskSexaquark::MCRec_PushBack(Int_t evt_track) {

    Int_t idx_track = 0;  // fVec_Idx_MCRec[evt_track];
    AliESDtrack* track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

    Int_t idx_true = TMath::Abs(track->GetLabel());
    AliMCParticle* mcTrue = (AliMCParticle*)fMC->GetTrack(idx_true);

    Int_t evt_true = 0;          // fMap_Evt_MCGen[idx_true];
    Bool_t is_signal = 0;        // MCGen_IsSignal(mcTrue);
    Bool_t is_duplicate = 0;     // fVec_MCRec_IsDuplicate[evt_track];
    Bool_t is_similar = kFALSE;  // PENDING! fVec_MCRec_IsSimilar[evt_track];
    Float_t impact_param = MCRec_GetImpactParameter(track);
    Double_t helix_param[6];
    track->GetHelixParameters(helix_param, fESD->GetMagneticField());

#if DEBUG_MODE > 0
    AliInfoF("%6i %6i %6i %6i %6i", evt_track, idx_track, evt_true, idx_true, is_duplicate);
#endif

    fEvent.Idx_True.push_back(evt_true);
    fEvent.Rec_Px.push_back(track->Px());
    fEvent.Rec_Py.push_back(track->Py());
    fEvent.Rec_Pz.push_back(track->Pz());
    fEvent.Rec_Charge.push_back(track->Charge());
    fEvent.Rec_IP_wrtPV.push_back(impact_param);
    fEvent.Rec_NSigmaPion.push_back(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion));
    fEvent.Rec_NSigmaProton.push_back(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    fEvent.Rec_NClustersTPC.push_back(track->GetTPCncls());
    fEvent.Rec_isDuplicate.push_back(is_duplicate);
    fEvent.Rec_isSimilar.push_back(is_similar);
    fEvent.Rec_isSignal.push_back(is_signal);
    fEvent.Rec_HelixParam0.push_back(helix_param[0]);
    fEvent.Rec_HelixParam1.push_back(helix_param[1]);
    fEvent.Rec_HelixParam2.push_back(helix_param[2]);
    fEvent.Rec_HelixParam3.push_back(helix_param[3]);
    fEvent.Rec_HelixParam4.push_back(helix_param[4]);
    fEvent.Rec_HelixParam5.push_back(helix_param[5]);
}

/*
 Push one of the found V0s into the Event_tt object
*/
void AliAnalysisTaskSexaquark::V0_PushBack(Int_t evt_v0) {

    V0_tt this_V0;  // = fVec_V0s[evt_v0];

#if DEBUG_MODE == 2
    Int_t evt_pos = this_V0.Idx_Pos;
    Int_t evt_neg = this_V0.Idx_Neg;
    Int_t idx_pos = 0;  // fVec_Idx_MCRec[evt_pos];
    Int_t idx_neg = 0;  // fVec_Idx_MCRec[evt_neg];
    AliInfoF("%6i %6i %6i %6i", evt_pos, evt_neg, idx_pos, idx_neg);
#endif

    fEvent.Idx_Pos.push_back(this_V0.Idx_Pos);
    fEvent.Idx_Neg.push_back(this_V0.Idx_Neg);
    fEvent.V0_Px.push_back(this_V0.Px);
    fEvent.V0_Py.push_back(this_V0.Py);
    fEvent.V0_Pz.push_back(this_V0.Pz);
    fEvent.V0_X.push_back(this_V0.X);
    fEvent.V0_Y.push_back(this_V0.Y);
    fEvent.V0_Z.push_back(this_V0.Z);
    fEvent.Pos_Px.push_back(this_V0.Pos_Px);
    fEvent.Pos_Py.push_back(this_V0.Pos_Py);
    fEvent.Pos_Pz.push_back(this_V0.Pos_Pz);
    fEvent.Neg_Px.push_back(this_V0.Neg_Px);
    fEvent.Neg_Py.push_back(this_V0.Neg_Py);
    fEvent.Neg_Pz.push_back(this_V0.Neg_Pz);
    fEvent.V0_isSignal.push_back(this_V0.isSignal);
    fEvent.V0_E_asK0.push_back(this_V0.E_asK0);
    fEvent.V0_E_asAL.push_back(this_V0.E_asAL);
    fEvent.V0_couldBeK0.push_back(this_V0.couldBeK0);
    fEvent.V0_couldBeAL.push_back(this_V0.couldBeAL);
    fEvent.V0_Chi2.push_back(this_V0.Chi2);
    fEvent.V0_DCA_Daughters.push_back(this_V0.DCA_Daughters);
    fEvent.V0_IP_wrtPV.push_back(this_V0.IP_wrtPV);
    fEvent.V0_CPA_wrtPV.push_back(this_V0.CPA_wrtPV);
    fEvent.V0_ArmAlpha.push_back(this_V0.ArmAlpha);
    fEvent.V0_ArmPt.push_back(this_V0.ArmPt);
    fEvent.V0_DecayLength.push_back(this_V0.DecayLength);
    fEvent.V0_DCA_wrtPV.push_back(this_V0.DCA_wrtPV);
}

/*

 */
void AliAnalysisTaskSexaquark::FillTree() {
    AliInfo("MONTE CARLO TREE");
    AliInfo("");
    AliInfo("         EVT    IDX    PID SIGNAL MOTHER  M_PID GM_PID FIRDAU LASDAU");
    /*
    fEvent.N_MCGen = (Int_t)fVec_Idx_MCGen.size();
    for (Int_t ii = 0; ii < fEvent.N_MCGen; ii++) {
        MCGen_PushBack(ii);
    }

    AliInfo("RECONSTRUCTED TRACKS");
    AliInfo("");
    AliInfo("         EVT    IDX EVTTRU IDXTRU DUPLI?");
    fEvent.N_MCRec = (Int_t)fVec_Idx_MCRec.size();
    for (Int_t ii = 0; ii < fEvent.N_MCRec; ii++) {
        MCRec_PushBack(ii);
    }

    AliInfo("FOUND V0s");
    AliInfo("");
    AliInfo("   EVTPOS EVTNEG IDXPOS IDXNEG");
    fEvent.N_V0s = (Int_t)fVec_V0s.size();
    for (Int_t ii = 0; ii < fEvent.N_V0s; ii++) {
        V0_PushBack(ii);
    }
 */
    // (tree operations) after all variables have been assigned, fill tree
    fTree->Fill();
}

/*
 Clear all indices and vectors loaded in the Event_tt object
*/
void AliAnalysisTaskSexaquark::ClearEvent() {
    /* MC particles */
    fEvent.N_MCGen = 0;
    fEvent.MC_Px.clear();
    fEvent.MC_Py.clear();
    fEvent.MC_Pz.clear();
    fEvent.MC_X.clear();
    fEvent.MC_Y.clear();
    fEvent.MC_Z.clear();
    fEvent.MC_Xf.clear();
    fEvent.MC_Yf.clear();
    fEvent.MC_Zf.clear();
    fEvent.MC_PID.clear();
    fEvent.MC_Mother.clear();
    fEvent.MC_PID_Mother.clear();
    fEvent.MC_PID_GrandMother.clear();
    fEvent.MC_NDaughters.clear();
    fEvent.MC_FirstDau.clear();
    fEvent.MC_LastDau.clear();
    fEvent.MC_Status.clear();
    fEvent.MC_isSignal.clear();
    /* MC Rec. Tracks */
    fEvent.N_MCRec = 0;
    fEvent.Idx_True.clear();
    fEvent.Rec_Px.clear();
    fEvent.Rec_Py.clear();
    fEvent.Rec_Pz.clear();
    fEvent.Rec_Charge.clear();
    fEvent.Rec_IP_wrtPV.clear();
    fEvent.Rec_NSigmaPion.clear();
    fEvent.Rec_NSigmaProton.clear();
    fEvent.Rec_NClustersTPC.clear();
    fEvent.Rec_isDuplicate.clear();
    fEvent.Rec_isSimilar.clear();
    fEvent.Rec_isSignal.clear();
    fEvent.Rec_HelixParam0.clear();
    fEvent.Rec_HelixParam1.clear();
    fEvent.Rec_HelixParam2.clear();
    fEvent.Rec_HelixParam3.clear();
    fEvent.Rec_HelixParam4.clear();
    fEvent.Rec_HelixParam5.clear();
    /* V0s */
    fEvent.N_V0s = 0;
    fEvent.Idx_Pos.clear();
    fEvent.Idx_Neg.clear();
    fEvent.V0_Px.clear();
    fEvent.V0_Py.clear();
    fEvent.V0_Pz.clear();
    fEvent.V0_X.clear();
    fEvent.V0_Y.clear();
    fEvent.V0_Z.clear();
    fEvent.Pos_Px.clear();
    fEvent.Pos_Py.clear();
    fEvent.Pos_Pz.clear();
    fEvent.Neg_Px.clear();
    fEvent.Neg_Py.clear();
    fEvent.Neg_Pz.clear();
    fEvent.V0_isSignal.clear();
    fEvent.V0_E_asK0.clear();
    fEvent.V0_E_asAL.clear();
    fEvent.V0_couldBeK0.clear();
    fEvent.V0_couldBeAL.clear();
    fEvent.V0_Chi2.clear();
    fEvent.V0_DCA_Daughters.clear();
    fEvent.V0_IP_wrtPV.clear();
    fEvent.V0_CPA_wrtPV.clear();
    fEvent.V0_ArmAlpha.clear();
    fEvent.V0_ArmPt.clear();
    fEvent.V0_DecayLength.clear();
    fEvent.V0_DCA_wrtPV.clear();
}

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

/*
 Correct initialization of a KFParticle.
 (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
*/
KFParticle AliAnalysisTaskSexaquark::CreateKFParticle(AliExternalTrackParam& track, Double_t mass, Int_t charge) {

    Double_t fP[6];
    track.GetXYZ(fP);
    track.PxPyPz(fP + 3);

    Int_t fQ = track.Charge() * TMath::Abs(charge);
    fP[3] *= TMath::Abs(charge);
    fP[4] *= TMath::Abs(charge);
    fP[5] *= TMath::Abs(charge);

    Double_t pt = 1. / TMath::Abs(track.GetParameter()[4]) * TMath::Abs(charge);
    Double_t cs = TMath::Cos(track.GetAlpha());
    Double_t sn = TMath::Sin(track.GetAlpha());
    Double_t r = TMath::Sqrt((1. - track.GetParameter()[2]) * (1. + track.GetParameter()[2]));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + track.GetParameter()[2] * cs / r);
    Double_t m43 = -pt * pt * (r * cs - track.GetParameter()[2] * sn);
    Double_t m24 = pt * (cs - track.GetParameter()[2] * sn / r);
    Double_t m44 = -pt * pt * (r * sn + track.GetParameter()[2] * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * track.GetParameter()[3];

    m43 *= track.GetSign();
    m44 *= track.GetSign();
    m45 *= track.GetSign();

    const Double_t* cTr = track.GetCovariance();
    Double_t fC[21];
    fC[0] = cTr[0] * m00 * m00;
    fC[1] = cTr[0] * m00 * m10;
    fC[2] = cTr[0] * m10 * m10;
    fC[3] = cTr[1] * m00;
    fC[4] = cTr[1] * m10;
    fC[5] = cTr[2];
    fC[6] = m00 * (cTr[3] * m23 + cTr[10] * m43);
    fC[7] = m10 * (cTr[3] * m23 + cTr[10] * m43);
    fC[8] = cTr[4] * m23 + cTr[11] * m43;
    fC[9] = m23 * (cTr[5] * m23 + cTr[12] * m43) + m43 * (cTr[12] * m23 + cTr[14] * m43);
    fC[10] = m00 * (cTr[3] * m24 + cTr[10] * m44);
    fC[11] = m10 * (cTr[3] * m24 + cTr[10] * m44);
    fC[12] = cTr[4] * m24 + cTr[11] * m44;
    fC[13] = m23 * (cTr[5] * m24 + cTr[12] * m44) + m43 * (cTr[12] * m24 + cTr[14] * m44);
    fC[14] = m24 * (cTr[5] * m24 + cTr[12] * m44) + m44 * (cTr[12] * m24 + cTr[14] * m44);
    fC[15] = m00 * (cTr[6] * m35 + cTr[10] * m45);
    fC[16] = m10 * (cTr[6] * m35 + cTr[10] * m45);
    fC[17] = cTr[7] * m35 + cTr[11] * m45;
    fC[18] = m23 * (cTr[8] * m35 + cTr[12] * m45) + m43 * (cTr[13] * m35 + cTr[14] * m45);
    fC[19] = m24 * (cTr[8] * m35 + cTr[12] * m45) + m44 * (cTr[13] * m35 + cTr[14] * m45);
    fC[20] = m35 * (cTr[9] * m35 + cTr[13] * m45) + m45 * (cTr[13] * m35 + cTr[14] * m45);

    KFParticle part;
    part.Create(fP, fC, fQ, mass);

    return part;
}

/*
 Correct initialization of a KFVertex.
 (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
*/
KFVertex AliAnalysisTaskSexaquark::CreateKFVertex(const AliVVertex& vertex) {

    Double_t param[6];
    vertex.GetXYZ(param);

    Double_t cov[6];
    vertex.GetCovarianceMatrix(cov);

    KFPVertex kfpVtx;
    Float_t paramF[3] = {(Float_t)param[0], (Float_t)param[1], (Float_t)param[2]};
    kfpVtx.SetXYZ(paramF);
    Float_t covF[6] = {(Float_t)cov[0], (Float_t)cov[1], (Float_t)cov[2], (Float_t)cov[3], (Float_t)cov[4], (Float_t)cov[5]};
    kfpVtx.SetCovarianceMatrix(covF);

    KFVertex KFVtx(kfpVtx);

    return KFVtx;
}
