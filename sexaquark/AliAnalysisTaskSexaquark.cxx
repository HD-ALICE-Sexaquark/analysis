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

#define HomogeneousField  // homogenous field in z direction, required by KFParticle
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
      fEvent() {}

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
      fEvent() {
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

    /* That will be filled with MC Gen. (or True) information */

    TString def_True_Signal_SecVertex_Radius = "2D Radius of Signal Products";
    TString def_True_InjBkg_SecVertex_Radius = "2D Radius of Injected AntiNeutrons";

    fHist_True_Signal_SecVertex_Radius = new TH1F("True_Signal_SecVertex_Radius", def_True_Signal_SecVertex_Radius, 100, 0., 200.);
    fHist_True_InjBkg_SecVertex_Radius = new TH1F("True_InjBkg_SecVertex_Radius", def_True_InjBkg_SecVertex_Radius, 100, 0., 200.);

    fOutputListOfHists->Add(fHist_True_Signal_SecVertex_Radius);
    fOutputListOfHists->Add(fHist_True_InjBkg_SecVertex_Radius);

    TString def_True_FermiSexaquark_Pt = "p_T of injected anti-sexaquark (affected by Fermi Motion)";
    TString def_True_FermiSexaquark_Mass = "Mass of injected anti-sexaquark (affected by Fermi Motion)";

    fHist_True_FermiSexaquark_Pt = new TH1F("True_FermiSexaquark_Pt", def_True_FermiSexaquark_Pt, 140, -7., 7.);
    fHist_True_FermiSexaquark_Mass = new TH1F("True_FermiSexaquark_Mass", def_True_FermiSexaquark_Mass, 100, -5., 5.);

    fOutputListOfHists->Add(fHist_True_FermiSexaquark_Pt);
    fOutputListOfHists->Add(fHist_True_FermiSexaquark_Mass);

    TString def_True_Signal_AntiLambda_CPAwrtPV = "CPA w.r.t. PV of Signal AntiLambda";
    TString def_True_Signal_KaonZeroShort_CPAwrtPV = "CPA w.r.t. PV of Signal KaonZeroShort";

    fHist_True_Signal_AntiLambda_CPAwrtPV = new TH1F("True_Signal_AntiLambda_CPAwrtPV", def_True_Signal_AntiLambda_CPAwrtPV, 100, -5., 5.);
    fHist_True_Signal_KaonZeroShort_CPAwrtPV = new TH1F("True_Signal_KaonZeroShort_CPAwrtPV", def_True_Signal_KaonZeroShort_CPAwrtPV, 100, -5., 5.);

    fOutputListOfHists->Add(fHist_True_Signal_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_True_Signal_KaonZeroShort_CPAwrtPV);

    TString def_True_AntiLambda_Bookkeep = "Bookkeeping of AntiLambdas";
    TString def_True_KaonZeroShort_Bookkeep = "Bookkeeping of KaonZeroShort";
    TString def_True_AntiProton_Bookkeep = "Bookkeeping of AntiProton";
    TString def_True_PosKaon_Bookkeep = "Bookkeeping of PosKaon";
    TString def_True_PiPlus_Bookkeep = "Bookkeeping of PiPlus";
    TString def_True_PiMinus_Bookkeep = "Bookkeeping of PiMinus";

    fHist_True_AntiLambda_Bookkeep = new TH1F("True_AntiLambda_Bookkeep", def_True_AntiLambda_Bookkeep, 10, 0., 10.);
    fHist_True_KaonZeroShort_Bookkeep = new TH1F("True_KaonZeroShort_Bookkeep", def_True_KaonZeroShort_Bookkeep, 10, 0., 10.);
    fHist_True_AntiProton_Bookkeep = new TH1F("True_AntiProton_Bookkeep", def_True_AntiProton_Bookkeep, 10, 0., 10.);
    fHist_True_PosKaon_Bookkeep = new TH1F("True_PosKaon_Bookkeep", def_True_PosKaon_Bookkeep, 10, 0., 10.);
    fHist_True_PiPlus_Bookkeep = new TH1F("True_PiPlus_Bookkeep", def_True_PiPlus_Bookkeep, 10, 0., 10.);
    fHist_True_PiMinus_Bookkeep = new TH1F("True_PiMinus_Bookkeep", def_True_PiMinus_Bookkeep, 10, 0., 10.);

    fOutputListOfHists->Add(fHist_True_AntiLambda_Bookkeep);
    fOutputListOfHists->Add(fHist_True_KaonZeroShort_Bookkeep);
    fOutputListOfHists->Add(fHist_True_AntiProton_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PosKaon_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PiPlus_Bookkeep);
    fOutputListOfHists->Add(fHist_True_PiMinus_Bookkeep);

    TString def_True_AntiNeutron_Pt = "p_T of generated anti-neutrons";
    TString def_True_AntiNeutronThatInteracted_Pt = "p_T of generated anti-neutrons that interacted with material";
    TString def_True_AntiNeutron_NDaughters = "";
    TString def_True_AntiNeutron_PDGDaughters = "";
    TString def_True_AntiNeutron_Bookkeep = "";

    fHist_True_AntiNeutron_Pt = new TH1F("True_AntiNeutron_Pt", def_True_AntiNeutron_Pt, 100, 0., 10.);
    fHist_True_AntiNeutronThatInteracted_Pt = new TH1F("True_AntiNeutronThatInteracted_Pt", def_True_AntiNeutronThatInteracted_Pt, 100, 0., 10.);
    fHist_True_AntiNeutron_NDaughters = new TH1F("True_AntiNeutron_NDaughters", def_True_AntiNeutron_NDaughters, 100, 0., 100.);
    fHist_True_AntiNeutron_PDGDaughters = new TH1F("True_AntiNeutron_PDGDaughters", def_True_AntiNeutron_PDGDaughters, 7000, -3500., 3500.);
    fHist_True_AntiNeutron_Bookkeep = new TH1F("True_AntiNeutron_Bookkeep", def_True_AntiNeutron_Bookkeep, 10, 0., 10.);

    fOutputListOfHists->Add(fHist_True_AntiNeutron_Pt);
    fOutputListOfHists->Add(fHist_True_AntiNeutronThatInteracted_Pt);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_NDaughters);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_PDGDaughters);
    fOutputListOfHists->Add(fHist_True_AntiNeutron_Bookkeep);

    TString def_True_InjBkgProducts_Bookkeep = "";
    TString def_True_AllSecondaries_MothersPDG = "PDG code of all particles that come from a sec. vertex";

    fHist_True_InjBkgProducts_Bookkeep = new TH1F("True_InjBkgProducts_Bookkeep", def_True_InjBkgProducts_Bookkeep, 50, 0., 50.);
    fHist_True_AllSecondaries_MothersPDG = new TH1F("True_AllSecondaries_MothersPDG", def_True_AllSecondaries_MothersPDG, 7000, -3500, 3500.);

    fOutputListOfHists->Add(fHist_True_InjBkgProducts_Bookkeep);
    fOutputListOfHists->Add(fHist_True_AllSecondaries_MothersPDG);

    /* That will be filled with rec. tracks information (that require MC Gen. information) */

    TString def_Signal_AntiProton_Pt = "p_T of identified signal anti-protons";
    TString def_Signal_PosKaon_Pt = "p_T of identified signal K+";
    TString def_Signal_PiPlus_Pt = "p_T of identified signal pi+";
    TString def_Signal_PiMinus_Pt = "p_T of identified signal pi-";

    fHist_Signal_AntiProton_Pt = new TH1F("Signal_AntiProton_Pt", def_Signal_AntiProton_Pt, 100, 0., 10.);
    fHist_Signal_PosKaon_Pt = new TH1F("Signal_PosKaon_Pt", def_Signal_PosKaon_Pt, 100, 0., 10.);
    fHist_Signal_PiPlus_Pt = new TH1F("Signal_PiPlus_Pt", def_Signal_PiPlus_Pt, 100, 0., 10.);
    fHist_Signal_PiMinus_Pt = new TH1F("Signal_PiMinus_Pt", def_Signal_PiMinus_Pt, 100, 0., 10.);

    fOutputListOfHists->Add(fHist_Signal_AntiProton_Pt);
    fOutputListOfHists->Add(fHist_Signal_PosKaon_Pt);
    fOutputListOfHists->Add(fHist_Signal_PiPlus_Pt);
    fOutputListOfHists->Add(fHist_Signal_PiMinus_Pt);

    /* That will be filled with rec. tracks information */

    TString def_AntiProton_Pt = "p_T of identified anti-protons";
    TString def_PosKaon_Pt = "p_T of identified K+";
    TString def_PiPlus_Pt = "p_T of identified pi+";
    TString def_PiMinus_Pt = "p_T of identified pi-";

    fHist_AntiProton_Pt = new TH1F("AntiProton_Pt", def_AntiProton_Pt, 100, 0., 10.);
    fHist_PosKaon_Pt = new TH1F("PosKaon_Pt", def_PosKaon_Pt, 100, 0., 10.);
    fHist_PiPlus_Pt = new TH1F("PiPlus_Pt", def_PiPlus_Pt, 100, 0., 10.);
    fHist_PiMinus_Pt = new TH1F("PiMinus_Pt", def_PiMinus_Pt, 100, 0., 10.);

    fOutputListOfHists->Add(fHist_AntiProton_Pt);
    fOutputListOfHists->Add(fHist_PosKaon_Pt);
    fOutputListOfHists->Add(fHist_PiPlus_Pt);
    fOutputListOfHists->Add(fHist_PiMinus_Pt);

    TString def_TPC_Signal = "TPC dE/dx vs p/z";

    f2DHist_TPC_Signal = new TH2F("TPC_Signal", def_TPC_Signal, 500, 0., 10., 500, 0., 1000.);

    fOutputListOfHists->Add(f2DHist_TPC_Signal);

    /* That will be filled with V0s information */

    TString def_AntiLambda_Mass = "m(#bar{p},#pi^{+})";
    TString def_AntiLambda_CPAwrtPV = "CPA(#bar{p},#pi^{+}) w.r.t. PV";
    TString def_AntiLambda_DCAwrtPV = "DCA(#bar{p},#pi^{+}) w.r.t. PV";
    TString def_KaonZeroShort_Mass = "m(#pi^{-},#pi^{+})";
    TString def_KaonZeroShort_CPAwrtPV = "CPA(#pi^{-},#pi^{+}) w.r.t. PV";
    TString def_KaonZeroShort_DCAwrtPV = "DCA(#pi^{-},#pi^{+}) w.r.t. PV";

    fHist_AntiLambda_Mass = new TH1F("AntiLambda_Mass", def_AntiLambda_Mass, 100, 0., 10.);
    fHist_AntiLambda_CPAwrtPV = new TH1F("AntiLambda_CPAwrtPV", def_AntiLambda_CPAwrtPV, 100, -5., 5.);
    fHist_AntiLambda_DCAwrtPV = new TH1F("AntiLambda_DCAwrtPV", def_AntiLambda_DCAwrtPV, 100, 0., 10.);
    fHist_KaonZeroShort_Mass = new TH1F("KaonZeroShort_Mass", def_KaonZeroShort_Mass, 100, 0., 10.);
    fHist_KaonZeroShort_CPAwrtPV = new TH1F("KaonZeroShort_CPAwrtPV", def_KaonZeroShort_CPAwrtPV, 100, -5., 5.);
    fHist_KaonZeroShort_DCAwrtPV = new TH1F("KaonZeroShort_DCAwrtPV", def_KaonZeroShort_DCAwrtPV, 100, 0., 10.);

    fOutputListOfHists->Add(fHist_AntiLambda_Mass);
    fOutputListOfHists->Add(fHist_AntiLambda_CPAwrtPV);
    fOutputListOfHists->Add(fHist_AntiLambda_DCAwrtPV);
    fOutputListOfHists->Add(fHist_KaonZeroShort_Mass);
    fOutputListOfHists->Add(fHist_KaonZeroShort_CPAwrtPV);
    fOutputListOfHists->Add(fHist_KaonZeroShort_DCAwrtPV);

    /* This will be filled with the reconstructed anti-sexaquark information */

    TString def_AntiSexaquark_Mass = "m(#bar{#Lambda},K^{0}_{S})";
    TString def_AntiSexaquark_DCA = "DCA(#bar{#Lambda},K^{0}_{S})";

    fHist_AntiSexaquark_Mass = new TH1F("AntiSexaquark_Mass", def_AntiSexaquark_Mass, 150, -5., 10.);
    fHist_AntiSexaquark_DCA = new TH1F("AntiSexaquark_DCA", def_AntiSexaquark_DCA, 100, 0., 10.);

    fOutputListOfHists->Add(fHist_AntiSexaquark_Mass);
    fOutputListOfHists->Add(fHist_AntiSexaquark_DCA);

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
 Main function, called per each event at RUNTIME ~ execution on Grid
*/
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {

    /* Declare containers */

    std::set<Int_t> Indices_MCGen_FS_Signal;

    std::vector<Int_t> idxAntiProtonTracks;
    std::vector<Int_t> idxPosKaonTracks;
    std::vector<Int_t> idxPiPlusTracks;
    std::vector<Int_t> idxPiMinusTracks;

    std::vector<std::pair<Int_t, Int_t>> idxAntiLambdaDaughters;
    std::vector<std::pair<Int_t, Int_t>> idxKaonZeroShortDaughters;
    std::vector<std::pair<Int_t, Int_t>> idxPionPairs;
    std::vector<std::pair<Int_t, Int_t>> idxPionKaonPairs;
    std::vector<std::pair<Int_t, Int_t>> idxPosKaonPairs;

    std::vector<KFParticleMother> kfAntiLambdas;
    std::vector<KFParticleMother> kfKaonsZeroShort;
    std::vector<KFParticleMother> kfPionPairs;
    std::vector<KFParticleMother> kfPionKaonPairs;
    std::vector<KFParticleMother> kfPosKaonPairs;

    std::vector<KFParticleMother> kfAntiSexaquarks;

    std::map<Int_t, Int_t> Map_TrueDaughter_TrueV0;
    std::map<Int_t, std::pair<Int_t, Int_t>> Map_TrueV0_RecDaughters;
    std::vector<std::pair<Int_t, Int_t>> RecDaughters_SignalV0;

    std::vector<AliESDv0> esdAntiLambdas;
    std::vector<AliESDv0> esdKaonsZeroShort;

    std::vector<Int_t> idxAntiSexaquarkDaughters;

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

    if (fIsMC) ProcessMCGen(Indices_MCGen_FS_Signal, Map_TrueDaughter_TrueV0);

    ProcessTracks(Indices_MCGen_FS_Signal, Map_TrueDaughter_TrueV0, Map_TrueV0_RecDaughters, idxAntiProtonTracks, idxPosKaonTracks, idxPiPlusTracks,
                  idxPiMinusTracks);

    if (fIsMC) ReconstructV0s_True(Map_TrueV0_RecDaughters, RecDaughters_SignalV0);

    if (fSourceOfV0s == "offline") {
        ReconstructV0s_Official(kFALSE, idxAntiLambdaDaughters, idxKaonZeroShortDaughters, RecDaughters_SignalV0);
    }

    if (fSourceOfV0s == "online") {
        ReconstructV0s_Official(kTRUE, idxAntiLambdaDaughters, idxKaonZeroShortDaughters, RecDaughters_SignalV0);
    }

    if (fSourceOfV0s == "custom") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            ReconstructV0s_Custom(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, esdAntiLambdas, idxAntiLambdaDaughters,
                                  RecDaughters_SignalV0);
            ReconstructV0s_Custom(idxPiMinusTracks, idxPiPlusTracks, 310, -211, 211, esdKaonsZeroShort, idxKaonZeroShortDaughters,
                                  RecDaughters_SignalV0);
        }
    }

    if (fSourceOfV0s == "kalman") {
        if (fReactionChannel == "AntiSexaquark,N->AntiLambda,K0S") {
            ReconstructV0s_KF(idxAntiProtonTracks, idxPiPlusTracks, -3122, -2212, 211, kfAntiLambdas, idxAntiLambdaDaughters);
            ReconstructV0s_KF(idxPiMinusTracks, idxPiPlusTracks, 310, -211, 211, kfKaonsZeroShort, idxKaonZeroShortDaughters);
        }
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
    }

    if (fSourceOfV0s == "custom") {
        SexaquarkFinder_ChannelA_Geo(esdAntiLambdas, esdKaonsZeroShort, idxAntiLambdaDaughters, idxKaonZeroShortDaughters, idxAntiSexaquarkDaughters);
    }

    // FillTree();

    // stream the results the analysis of this event to the output manager
    PostData(1, fOutputListOfTrees);
    PostData(2, fOutputListOfHists);
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
 - Uses: `fMC`, `fPDG`, `fHist_True_*`
 - Output: `Indices_MCGen_FS_Signal`, `Map_TrueDaughter_TrueV0`
*/
void AliAnalysisTaskSexaquark::ProcessMCGen(std::set<Int_t>& Indices_MCGen_FS_Signal, std::map<Int_t, Int_t>& Map_TrueDaughter_TrueV0) {

    AliMCParticle* mcPart;
    Int_t pdg_mc;

    AliMCParticle* mcMother;
    Int_t pdg_mother;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    Int_t n_antilambda = 0;
    Int_t n_k0s = 0;
    Int_t n_antineutrons = 0;

    Int_t n_piplus = 0;
    Int_t n_piminus = 0;
    Int_t n_kplus = 0;
    Int_t n_antiprotons = 0;

    Bool_t is_secondary;
    Bool_t is_signal;

    Float_t antineutron_radius;

    // `key` = interaction id, `value` = vector of indices of first gen. daughters
    std::map<Int_t, std::vector<Int_t>> vec_idx_signal_fg_daughters;

    Float_t cpa_wrt_pv;
    Double_t PV[3];
    fMC_PrimaryVertex->GetXYZ(PV);

    /* Loop over MC gen. particles in a single event */

    AliInfo("   IDX    STATUS    PID MOMIDX MOMPID           ORIGIN           ISPRIM");
    for (Int_t idx_mc = 0; idx_mc < fMC->GetNumberOfTracks(); idx_mc++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(idx_mc);
        pdg_mc = mcPart->PdgCode();

        /* Fill bookkeeping histograms */
        /* PENDING: it can be simplified */

        if (pdg_mc == -3122) {
            fHist_True_AntiLambda_Bookkeep->Fill(0);
            for (Int_t idx_dau = mcPart->GetDaughterFirst(); mcPart->GetNDaughters() && idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {
                mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);
                pdg_dau = mcDaughter->PdgCode();
                if (pdg_dau == -2212)
                    Map_TrueDaughter_TrueV0[idx_dau] = idx_mc;
                else if (pdg_dau == 211)
                    Map_TrueDaughter_TrueV0[idx_dau] = idx_mc;
            }
        }
        if (pdg_mc == 310) {
            fHist_True_KaonZeroShort_Bookkeep->Fill(0);
            for (Int_t idx_dau = mcPart->GetDaughterFirst(); mcPart->GetNDaughters() && idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {
                mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);
                pdg_dau = mcDaughter->PdgCode();
                if (TMath::Abs(pdg_dau) == 211) Map_TrueDaughter_TrueV0[idx_dau] = idx_mc;
            }
        }
        if (pdg_mc == -2212) fHist_True_AntiProton_Bookkeep->Fill(0);
        if (pdg_mc == 321) fHist_True_PosKaon_Bookkeep->Fill(0);
        if (pdg_mc == 211) fHist_True_PiPlus_Bookkeep->Fill(0);
        if (pdg_mc == -211) fHist_True_PiMinus_Bookkeep->Fill(0);
        if (pdg_mc == -2112) {
            fHist_True_AntiNeutron_Bookkeep->Fill(0);
            fHist_True_AntiNeutron_Pt->Fill(mcPart->Pt());
        }

        /* If they're signal, their MC status number should be among 600 and 619. */
        /* This was done to recognize each interaction, assigning a different number to each of */
        /* the 20 antisexaquark-nucleon reactions that were injected in a single event */

        is_signal = mcPart->MCStatusCode() >= 600 && mcPart->MCStatusCode() < 620;

        /* Count secondary particles. */
        /* -- Note 1: Include signal particles, as well, because they were injected as primaries. */
        /* -- Note 2: The K0S (PDG Code = 310) that transformed from K0 (PDG Code = 311) are considered primaries, */
        /*            that's why it's important to exclude them from this count. */

        is_secondary = mcPart->GetMother() != -1;

        if (pdg_mc == -3122 && (is_secondary || is_signal)) fHist_True_AntiLambda_Bookkeep->Fill(1);
        if (pdg_mc == 310) {
            if (is_secondary) {
                mcMother = (AliMCParticle*)fMC->GetTrack(mcPart->GetMother());
                pdg_mother = mcMother->PdgCode();
                if (TMath::Abs(pdg_mother) != 311) fHist_True_KaonZeroShort_Bookkeep->Fill(1);
            } else if (is_signal) {
                fHist_True_KaonZeroShort_Bookkeep->Fill(1);
            }
        }
        if (pdg_mc == -2212 && (is_secondary || is_signal)) fHist_True_AntiProton_Bookkeep->Fill(1);
        if (pdg_mc == 321 && (is_secondary || is_signal)) fHist_True_PosKaon_Bookkeep->Fill(1);
        if (pdg_mc == 211 && (is_secondary || is_signal)) fHist_True_PiPlus_Bookkeep->Fill(1);
        if (pdg_mc == -211 && (is_secondary || is_signal)) fHist_True_PiMinus_Bookkeep->Fill(1);

        /* Count signal particles */
        /* -- Note: Only the first generation daughters right after the reaction */
        /* (PENDING) I should also check if the branching ratio of the decay of the AL and K0S correspond to the real one */

        if (pdg_mc == -3122 && is_signal) {
            cpa_wrt_pv = Calculate_CPA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->Xv(), mcPart->Yv(), mcPart->Zv(), PV[0], PV[1], PV[2]);
            fHist_True_Signal_AntiLambda_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_True_AntiLambda_Bookkeep->Fill(2);
            for (Int_t idx_dau = mcPart->GetDaughterFirst(); idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {
                mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);
                pdg_dau = mcDaughter->PdgCode();
                if (pdg_dau == -2212) {
                    fHist_True_AntiProton_Bookkeep->Fill(2);
                    Indices_MCGen_FS_Signal.insert(idx_dau);
                } else if (pdg_dau == 211) {
                    fHist_True_PiPlus_Bookkeep->Fill(2);
                    Indices_MCGen_FS_Signal.insert(idx_dau);
                }
            }
            if (fReactionID == 'A' || fReactionID == 'D' || fReactionID == 'E') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == 310 && is_signal) {
            cpa_wrt_pv = Calculate_CPA(mcPart->Px(), mcPart->Py(), mcPart->Pz(), mcPart->Xv(), mcPart->Yv(), mcPart->Zv(), PV[0], PV[1], PV[2]);
            fHist_True_Signal_KaonZeroShort_CPAwrtPV->Fill(cpa_wrt_pv);
            fHist_True_KaonZeroShort_Bookkeep->Fill(2);
            for (Int_t idx_dau = mcPart->GetDaughterFirst(); idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {
                mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);
                pdg_dau = mcDaughter->PdgCode();
                if (pdg_dau == -211) {
                    fHist_True_PiMinus_Bookkeep->Fill(2);
                    Indices_MCGen_FS_Signal.insert(idx_dau);
                } else if (pdg_dau == 211) {
                    fHist_True_PiPlus_Bookkeep->Fill(2);
                    Indices_MCGen_FS_Signal.insert(idx_dau);
                }
            }
            if (fReactionID == 'A') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == -2212 && is_signal) {
            fHist_True_AntiProton_Bookkeep->Fill(2);
            Indices_MCGen_FS_Signal.insert(idx_mc);
            if (fReactionID == 'H') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == 321 && is_signal) {
            fHist_True_PosKaon_Bookkeep->Fill(2);
            Indices_MCGen_FS_Signal.insert(idx_mc);
            if (fReactionID == 'D' || fReactionID == 'E' || fReactionID == 'H') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == 211 && is_signal) {
            fHist_True_PiPlus_Bookkeep->Fill(2);
            Indices_MCGen_FS_Signal.insert(idx_mc);
            if (fReactionID == 'E') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == -211 && is_signal) {
            fHist_True_PiMinus_Bookkeep->Fill(2);
            Indices_MCGen_FS_Signal.insert(idx_mc);
            if (fReactionID == 'E') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        if (pdg_mc == 111 && is_signal) {
            if (fReactionID == 'H') {
                vec_idx_signal_fg_daughters[(Int_t)mcPart->MCStatusCode()].push_back(idx_mc);
            }
        }

        /* Count anti-neutrons that had a daughter or more */

        if (pdg_mc == -2112 && mcPart->GetNDaughters() > 0) {

            fHist_True_AntiNeutron_Bookkeep->Fill(1);
            fHist_True_AntiNeutron_NDaughters->Fill(mcPart->GetNDaughters());
            fHist_True_AntiNeutronThatInteracted_Pt->Fill(mcPart->Pt());

            /* Fill, then, the PDG code of its daughters */

            for (Int_t idx_dau = mcPart->GetDaughterFirst(); idx_dau <= mcPart->GetDaughterLast(); idx_dau++) {

                mcDaughter = (AliMCParticle*)fMC->GetTrack(idx_dau);
                pdg_dau = mcDaughter->PdgCode();

                if (TMath::Abs(pdg_dau) < 3500) {
                    fHist_True_AntiNeutron_PDGDaughters->Fill(pdg_dau);
                } else {
                    fHist_True_AntiNeutron_PDGDaughters->Fill(0);
                }

                if (idx_dau == mcPart->GetDaughterFirst()) {
                    antineutron_radius = TMath::Sqrt(TMath::Power(mcDaughter->Xv(), 2) + TMath::Power(mcDaughter->Yv(), 2));
                    fHist_True_InjBkg_SecVertex_Radius->Fill(antineutron_radius);
                }
            }
        }
    }  // end of loop over MC particles

    /* Loop over antisexaquark-nucleon interactions */

    Float_t antisexaquark_pt;
    Float_t antisexaquark_px;
    Float_t antisexaquark_py;
    Float_t antisexaquark_pz;
    Float_t antisexaquark_p;
    Float_t antisexaquark_e;
    Float_t antisexaquark_mass;
    Float_t signal_radius;

    for (Int_t interaction_id = 600; interaction_id < 620; interaction_id++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(vec_idx_signal_fg_daughters[interaction_id][0]);
        signal_radius = TMath::Sqrt(TMath::Power(mcPart->Xv(), 2) + TMath::Power(mcPart->Yv(), 2));

        fHist_True_Signal_SecVertex_Radius->Fill(signal_radius);

        antisexaquark_px = 0.;
        antisexaquark_py = 0.;
        antisexaquark_pz = 0.;
        antisexaquark_e = 0.;

        for (Int_t& signal_prod_idx : vec_idx_signal_fg_daughters[interaction_id]) {

            mcPart = (AliMCParticle*)fMC->GetTrack(signal_prod_idx);
            pdg_mc = mcPart->PdgCode();

            antisexaquark_px += mcPart->Px();
            antisexaquark_py += mcPart->Py();
            antisexaquark_pz += mcPart->Pz();
            antisexaquark_e += TMath::Sqrt(mcPart->P() * mcPart->P() + fPDG.GetParticle(pdg_mc)->Mass() * fPDG.GetParticle(pdg_mc)->Mass());
        }
        antisexaquark_e -= fPDG.GetParticle(fStruckNucleonPDG)->Mass();

        antisexaquark_pt = TMath::Sqrt(TMath::Power(antisexaquark_px, 2) + TMath::Power(antisexaquark_py, 2));
        fHist_True_FermiSexaquark_Pt->Fill(antisexaquark_pt);

        antisexaquark_p = TMath::Sqrt(TMath::Power(antisexaquark_px, 2) + TMath::Power(antisexaquark_py, 2) + TMath::Power(antisexaquark_pz, 2));

        antisexaquark_mass = TMath::Sqrt(TMath::Power(antisexaquark_e, 2) - TMath::Power(antisexaquark_p, 2));
        fHist_True_FermiSexaquark_Mass->Fill(antisexaquark_mass);
    }
}

/*                */
/**    Tracks    **/
/*** ========== ***/

/*
 Determine if current AliESDtrack passes track selection,
 and fill the bookkeeping histograms to measure the effect of the cuts.
 - Input: `track`, `is_signal`
 - Output: `n_sigma_pion`, `n_sigma_kaon`, `n_sigma_proton`
 */
Bool_t AliAnalysisTaskSexaquark::PassesTrackSelection(AliESDtrack* track, Bool_t is_signal, Float_t& n_sigma_proton, Float_t& n_sigma_kaon,
                                                      Float_t& n_sigma_pion) {

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

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(3);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(3);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(3);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(3);

    // >> eta
    if (TMath::Abs(track->Eta()) > 0.9) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(4);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(4);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(4);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(4);

    // >> TPC clusters
    if (track->GetTPCNcls() < 50) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(5);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(5);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(5);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(5);

    // >> chi2 per TPC cluster
    if (track->GetTPCchi2() / (Float_t)track->GetTPCNcls() > 7.0) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(6);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(6);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(6);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(6);

    // >> require TPC refit
    if (!(track->GetStatus() & AliESDtrack::kTPCrefit)) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(7);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(7);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(7);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(7);

    // >> reject kinks
    if (track->GetKinkIndex(0) != 0) {
        return kFALSE;
    }

    if (n_sigma_proton < kMaxNSigma_Proton && is_signal && track->Charge() < 0.) fHist_True_AntiProton_Bookkeep->Fill(8);
    if (n_sigma_kaon < kMaxNSigma_Kaon && is_signal && track->Charge() > 0.) fHist_True_PosKaon_Bookkeep->Fill(8);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() > 0.) fHist_True_PiPlus_Bookkeep->Fill(8);
    if (n_sigma_pion < kMaxNSigma_Pion && is_signal && track->Charge() < 0.) fHist_True_PiMinus_Bookkeep->Fill(8);

    return kTRUE;
}

/*
 Loop over the reconstructed tracks in a single event.
 - Input: `Indices_MCGen_FS_Signal`, `Map_TrueV0_TrueDaughters`
 - Output: `Map_TrueV0_RecDaughters`, `idxAntiProtonTracks`, `idxPosKaonTracks`, `idxPiPlusTracks`, `idxPiMinusTracks`
*/
void AliAnalysisTaskSexaquark::ProcessTracks(std::set<Int_t> Indices_MCGen_FS_Signal, std::map<Int_t, Int_t> Map_TrueDaughter_TrueV0,
                                             std::map<Int_t, std::pair<Int_t, Int_t>>& Map_TrueV0_RecDaughters,
                                             std::vector<Int_t>& idxAntiProtonTracks, std::vector<Int_t>& idxPosKaonTracks,
                                             std::vector<Int_t>& idxPiPlusTracks, std::vector<Int_t>& idxPiMinusTracks) {

    AliESDtrack* track;

    Int_t idx_true;
    Bool_t is_signal;

    Float_t n_sigma_pion;
    Float_t n_sigma_kaon;
    Float_t n_sigma_proton;

    for (Int_t idx_track = 0; idx_track < fESD->GetNumberOfTracks(); idx_track++) {

        track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_track));

        idx_true = TMath::Abs(track->GetLabel());
        is_signal = (Int_t)Indices_MCGen_FS_Signal.count(idx_true) > 0;  // PENDING: check if binary search is faster than linear search

        if (!PassesTrackSelection(track, is_signal, n_sigma_proton, n_sigma_kaon, n_sigma_pion)) continue;

        f2DHist_TPC_Signal->Fill(track->P(), track->GetTPCsignal());

        if (n_sigma_proton < kMaxNSigma_Proton && track->Charge() < 0.) {
            idxAntiProtonTracks.push_back(idx_track);
            fHist_AntiProton_Pt->Fill(track->Pt());
            if (is_signal) fHist_Signal_AntiProton_Pt->Fill(track->Pt());
        }

        if (n_sigma_kaon < kMaxNSigma_Kaon && track->Charge() > 0.) {
            idxPosKaonTracks.push_back(idx_track);
            fHist_PosKaon_Pt->Fill(track->Pt());
            if (is_signal) fHist_Signal_PosKaon_Pt->Fill(track->Pt());
        }

        if (n_sigma_pion < kMaxNSigma_Pion) {
            if (track->Charge() > 0.) {
                idxPiPlusTracks.push_back(idx_track);
                fHist_PiPlus_Pt->Fill(track->Pt());
                if (is_signal) fHist_Signal_PiPlus_Pt->Fill(track->Pt());
            } else {
                idxPiMinusTracks.push_back(idx_track);
                fHist_PiMinus_Pt->Fill(track->Pt());
                if (is_signal) fHist_Signal_PiMinus_Pt->Fill(track->Pt());
            }
        }

        /* Assign this track to a true V0. Convoluted, but it works. */

        if (track->Charge() < 0.)
            Map_TrueV0_RecDaughters[Map_TrueDaughter_TrueV0[idx_true]].first = idx_track;
        else
            Map_TrueV0_RecDaughters[Map_TrueDaughter_TrueV0[idx_true]].second = idx_track;
    }  // end of loop over tracks
}

/*                             */
/**  V0s -- ALICE V0 Finders  **/
/*** ======================= ***/

/*
 Apply cuts to an (anti-lambda) candidate, reconstructed from an (anti-proton, pi+) pair. Relevant for channel A and D.
 - Input: `v0`, `neg_track`, `pos_track`
*/
Bool_t AliAnalysisTaskSexaquark::PassesAntiLambdaCuts_OfficialCustom(AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track) {

    Float_t neg_nsigma_proton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kProton));
    Float_t pos_nsigma_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion));

    // >> particle identification
    if (neg_nsigma_proton > kMaxNSigma_Proton || pos_nsigma_pion > kMaxNSigma_Pion) {
        return kFALSE;
    }

    return kTRUE;
}

/*
 Apply cuts to a (K0S) candidate, reconstructed from a (pi-, pi+) pair. Relevant for channel A.
 - Input: `v0`, `neg_track`, `pos_track`
*/
Bool_t AliAnalysisTaskSexaquark::PassesKaonZeroShortCuts_OfficialCustom(AliESDv0* v0, AliESDtrack* neg_track, AliESDtrack* pos_track) {

    Float_t neg_nsigma_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(neg_track, AliPID::kPion));
    Float_t pos_nsigma_pion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pos_track, AliPID::kPion));

    // >> particle identification
    if (neg_nsigma_pion > kMaxNSigma_Pion || pos_nsigma_pion > kMaxNSigma_Pion) {
        return kFALSE;
    }

    return kTRUE;
}

/*
 Find all true V0s that had both of their daughters reconstructed. Useful for bookkeeping.
 - Input: `Map_TrueV0_RecDaughters`
 - Output: `RecDaughters_SignalV0`
 - Notes: to get any rec. value from the true V0s, I need to do vertexing first. That means:
   (a) true vertexing & true values => trivial and already done by `ProcessMCGen()`
   (b) true vertexing & rec. values => feels artificial...
   (c) standard vertexing (a-la official V0 finders)
   (d) KF vertexing
   Considering (c) and (d), then this function should be used only as bookkeeping, and for calling the real "reconstructers"
*/
void AliAnalysisTaskSexaquark::ReconstructV0s_True(std::map<Int_t, std::pair<Int_t, Int_t>> Map_TrueV0_RecDaughters,
                                                   std::vector<std::pair<Int_t, Int_t>>& RecDaughters_SignalV0) {

    AliMCParticle* mcV0;

    Int_t idx_mc_V0;

    Bool_t is_signal;

    Int_t idx_rec_neg_dau;
    Int_t idx_rec_pos_dau;

    /* Loop over map */

    for (auto it = Map_TrueV0_RecDaughters.begin(); it != Map_TrueV0_RecDaughters.end(); it++) {

        idx_mc_V0 = it->first;

        mcV0 = (AliMCParticle*)fMC->GetTrack(idx_mc_V0);

        idx_rec_neg_dau = Map_TrueV0_RecDaughters[idx_mc_V0].first;
        idx_rec_pos_dau = Map_TrueV0_RecDaughters[idx_mc_V0].second;

        if (idx_rec_neg_dau && idx_rec_pos_dau) {

            is_signal = mcV0->MCStatusCode() >= 600 && mcV0->MCStatusCode() < 620;

            if (mcV0->PdgCode() == -3122) {
                fHist_True_AntiLambda_Bookkeep->Fill(3);
                if (is_signal) {
                    fHist_True_AntiLambda_Bookkeep->Fill(4);
                    RecDaughters_SignalV0.push_back(Map_TrueV0_RecDaughters[idx_mc_V0]);
                    AliInfoF("-3122 %i %i", Map_TrueV0_RecDaughters[idx_mc_V0].first, Map_TrueV0_RecDaughters[idx_mc_V0].second);  // debug
                }
            }
            if (mcV0->PdgCode() == 310) {
                fHist_True_KaonZeroShort_Bookkeep->Fill(3);
                if (is_signal) {
                    fHist_True_KaonZeroShort_Bookkeep->Fill(4);
                    RecDaughters_SignalV0.push_back(Map_TrueV0_RecDaughters[idx_mc_V0]);
                    AliInfoF(" 310 %i %i", Map_TrueV0_RecDaughters[idx_mc_V0].first, Map_TrueV0_RecDaughters[idx_mc_V0].second);  // debug
                }
            }
        }
    }

    return;
}

/*
 Loop over all the V0s found by the Official V0 Finders that are stored in the ESD file.
 - Input: `online`
 - Output: `idxAntiLambdaDaughters`, `idxKaonZeroShortDaughters`
*/
void AliAnalysisTaskSexaquark::ReconstructV0s_Official(Bool_t online, std::vector<std::pair<Int_t, Int_t>>& idxAntiLambdaDaughters,
                                                       std::vector<std::pair<Int_t, Int_t>>& idxKaonZeroShortDaughters,
                                                       std::vector<std::pair<Int_t, Int_t>> RecDaughters_SignalV0) {

    // define variables & objects
    AliESDv0* V0;

    Int_t idx_pos;
    Int_t idx_neg;

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

    Int_t idx_pos_true;
    Int_t idx_neg_true;
    AliMCParticle* mcPosTrue;
    AliMCParticle* mcNegTrue;

    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

    std::pair<Int_t, Int_t> aux_pair;
    Float_t aux_mass;

    /* Loop over V0s that are stored in the ESD */

    for (Int_t idx_v0 = 0; idx_v0 < fESD->GetNumberOfV0s(); idx_v0++) {

        V0 = fESD->GetV0(idx_v0);

        /* Choose between offline and online (on-the-fly) V0s */

        if (!online && V0->GetOnFlyStatus()) {
            continue;
        }

        if (online && !V0->GetOnFlyStatus()) {
            continue;
        }

        idx_neg = V0->GetNindex();
        idx_pos = V0->GetPindex();

        neg_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_neg));
        pos_track = static_cast<AliESDtrack*>(fESD->GetTrack(idx_pos));

        /* Store: anti-lambda -> anti-proton + pi+ */

        if (PassesAntiLambdaCuts_OfficialCustom(V0, neg_track, pos_track)) {

            aux_pair.first = idx_neg;
            aux_pair.second = idx_pos;
            idxAntiLambdaDaughters.push_back(aux_pair);

            V0->GetXYZ(V0_X, V0_Y, V0_Z);
            V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
            V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

            pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(211)->Mass() * fPDG.GetParticle(211)->Mass());
            neg_energy = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(-2212)->Mass() * fPDG.GetParticle(-2212)->Mass());
            aux_mass = TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) -  //
                                   V0->Px() * V0->Px() - V0->Py() * V0->Py() - V0->Pz() * V0->Pz());

            /* Fill histograms */

            fHist_AntiLambda_Mass->Fill(aux_mass);
            fHist_AntiLambda_CPAwrtPV->Fill(V0->GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]));
            fHist_AntiLambda_DCAwrtPV->Fill(Calculate_LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]));

            /* DEBUG */

            if (std::any_of(RecDaughters_SignalV0.begin(), RecDaughters_SignalV0.end(),
                            [&](std::pair<Int_t, Int_t> signal_pair) { return signal_pair == aux_pair; })) {
                AliInfoF("-3122 %i %i", aux_pair.first, aux_pair.second);
            }
        }

        /* Store: K0S -> pi- + pi+ */

        if (PassesKaonZeroShortCuts_OfficialCustom(V0, neg_track, pos_track)) {

            aux_pair.first = idx_neg;
            aux_pair.second = idx_pos;
            idxKaonZeroShortDaughters.push_back(aux_pair);

            V0->GetXYZ(V0_X, V0_Y, V0_Z);
            V0->GetPPxPyPz(P_Px, P_Py, P_Pz);
            V0->GetNPxPyPz(N_Px, N_Py, N_Pz);

            pos_energy = TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(211)->Mass() * fPDG.GetParticle(211)->Mass());
            neg_energy = TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(-211)->Mass() * fPDG.GetParticle(-211)->Mass());
            aux_mass = TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) -  //
                                   V0->Px() * V0->Px() - V0->Py() * V0->Py() - V0->Pz() * V0->Pz());

            /* Fill histograms */

            fHist_KaonZeroShort_Mass->Fill(aux_mass);
            fHist_KaonZeroShort_CPAwrtPV->Fill(V0->GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]));
            fHist_KaonZeroShort_DCAwrtPV->Fill(Calculate_LinePointDCA(V0->Px(), V0->Py(), V0->Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]));

            /* DEBUG */

            if (std::any_of(RecDaughters_SignalV0.begin(), RecDaughters_SignalV0.end(),
                            [&](std::pair<Int_t, Int_t> signal_pair) { return signal_pair == aux_pair; })) {
                AliInfoF("310 %i %i", aux_pair.first, aux_pair.second);
            }
        }

    }  // end of loop over V0s
}

/*
 Custom V0 Finder.
 (Copied and adapted from `AliRoot/STEER/ESD/AliV0vertexer.cxx`)
 - Input: `idxNegativeTracks`, `idxPositiveTracks`, `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 - Output: `esdV0s`, `idxDaughters`
*/
void AliAnalysisTaskSexaquark::ReconstructV0s_Custom(std::vector<Int_t> idxNegativeTracks, std::vector<Int_t> idxPositiveTracks, Int_t pdgV0,
                                                     Int_t pdgTrackNeg, Int_t pdgTrackPos, std::vector<AliESDv0>& esdV0s,
                                                     std::vector<std::pair<Int_t, Int_t>>& idxDaughters,
                                                     std::vector<std::pair<Int_t, Int_t>> RecDaughters_SignalV0) {

    // AliESDtrack* track;
    AliESDtrack* neg_track;
    AliESDtrack* pos_track;

    Int_t nentr = fESD->GetNumberOfTracks();
    TArrayI neg(nentr);
    TArrayI pos(nentr);

    Double_t V0_X, V0_Y, V0_Z;
    Double_t N_Px, N_Py, N_Pz;
    Double_t P_Px, P_Py, P_Pz;

    Double_t pos_energy;
    Double_t neg_energy;
    Double_t aux_mass;

    Double_t radius;
    Double_t dca_wrt_pv;

    std::pair<Int_t, Int_t> aux_pair;

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

    Double_t PV[3];
    fPrimaryVertex->GetXYZ(PV);

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
    if (babs > 1e-3) {
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

    /*
    // number of charged tracks found
    Int_t nneg = 0;
    Int_t npos = 0;
        // first loop: find neg and pos tracks
        for (Int_t idx_track = 0; idx_track < nentr; idx_track++) {

            track = fESD->GetTrack(idx_track);

            // (cut) on status
            ULong64_t status = track->GetStatus();
            if ((status & AliESDtrack::kTPCrefit) == 0) {
                continue;
            }

            // (cut) min. number of TPC clusters
            if ((Int_t)track->GetTPCncls() < COV0F_NClustersTPCmin_Dau) {
                continue;
            }

            // cut on transverse impact parameter
            Double_t d = 0.0;
            if (track->GetSign() < 0.) {
                d = track->GetD(PV[0], PV[1], b);
                if (TMath::Abs(d) < COV0F_DNmin) {
                    continue;
                }
                if (TMath::Abs(d) > COV0F_Rmax) {
                    continue;
                }
                neg[nneg++] = idx_track;
            } else {
                d = track->GetD(PV[0], PV[1], b);
                if (TMath::Abs(d) < COV0F_DPmin) {
                    continue;
                }
                if (TMath::Abs(d) > COV0F_Rmax) {
                    continue;
                }
                pos[npos++] = idx_track;
            }
        }  // end of first loop
     */

    // second loop: test all neg-pos pairs
    // for (Int_t i = 0; i < nneg; i++) {
    for (Int_t& nidx : idxNegativeTracks) {

        // Int_t nidx = neg[i];
        neg_track = fESD->GetTrack(nidx);

        // for (Int_t k = 0; k < npos; k++) {
        for (Int_t& pidx : idxPositiveTracks) {

            // Int_t pidx = pos[k];
            pos_track = fESD->GetTrack(pidx);

            // (cut) -1. < eta(pi+) < 1.
            if (TMath::Abs(pos_track->Eta()) > COV0F_Etamax_PosDau) {
                continue;
            }

            // (cut) -1.5 < eta(pi-,anti-p) < 1.5
            if (TMath::Abs(neg_track->Eta()) > COV0F_Etamax_NegDau) {
                continue;
            }

            /* PENDING: could be uplifted as well */

            AliExternalTrackParam nt(*neg_track), pt(*pos_track), *ntp = &nt, *ptp = &pt;

            Double_t xn, xp;
            if (UseImprovedFinding && Preoptimize(ntp, ptp, &xn, &xp, b)) {
                // Move tracks to a better position if that helps
                nt.PropagateTo(xn, b);
                pt.PropagateTo(xp, b);
            }

            // (cut) DCA between daughters
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

            /* PENDING: could be uplifted outside of the loop... could be, though? */

            AliESDv0 vertex(nt, nidx, pt, pidx);

            if (UseImprovedFinding) {
                vertex.Refit();  // imp pos + cov mat
            }

            // (cut) chi2 ~ quality of fit
            if (vertex.GetChi2V0() > COV0F_Chi2max) {
                continue;
            }

            vertex.SetDcaV0Daughters(dca);

            // (cut) CPA
            Float_t cpa = vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]);
            vertex.SetV0CosineOfPointingAngle(cpa);
            if (cpa < COV0F_CPAmin) {
                continue;
            }

            // (cut) radius(V0)
            Double_t x = vertex.Xv();
            Double_t y = vertex.Yv();
            Double_t v0_radius = TMath::Sqrt(x * x + y * y);
            if (v0_radius < COV0F_Rmin) {
                continue;
            }
            if (v0_radius > COV0F_Rmax) {
                continue;
            }

            /* Load coordinates and momentum of each component of the V0 */

            vertex.GetXYZ(V0_X, V0_Y, V0_Z);

            // (cut) dca(V0,PV)
            dca_wrt_pv = Calculate_LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]);
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

            // (addition) (cut) remove v0s "not contributing to physics"
            if (TMath::Abs(vertex.Eta()) > COV0F_Etamax_V0) {
                continue;
            }

            // get momentum of negative and positive components at V0 vertex
            vertex.GetNPxPyPz(N_Px, N_Py, N_Pz);
            vertex.GetPPxPyPz(P_Px, P_Py, P_Pz);

            // (cut) on Pt
            if (TMath::Sqrt(vertex.Px() * vertex.Px() + vertex.Py() * vertex.Py()) < COV0F_Ptmin) {
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

            // calculate energy
            pos_energy =
                TMath::Sqrt(P_Px * P_Px + P_Py * P_Py + P_Pz * P_Pz + fPDG.GetParticle(pdgTrackPos)->Mass() * fPDG.GetParticle(pdgTrackPos)->Mass());
            neg_energy =
                TMath::Sqrt(N_Px * N_Px + N_Py * N_Py + N_Pz * N_Pz + fPDG.GetParticle(pdgTrackNeg)->Mass() * fPDG.GetParticle(pdgTrackNeg)->Mass());
            aux_mass = TMath::Sqrt((pos_energy + neg_energy) * (pos_energy + neg_energy) - vertex.Px() * vertex.Px() - vertex.Py() * vertex.Py() -
                                   vertex.Pz() * vertex.Pz());

            /* Fill histograms */

            if (pdgV0 == -3122) {
                fHist_AntiLambda_Mass->Fill(aux_mass);
                fHist_AntiLambda_CPAwrtPV->Fill(vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]));
                fHist_AntiLambda_DCAwrtPV->Fill(dca_wrt_pv);
            }

            if (pdgV0 == 310) {
                fHist_KaonZeroShort_Mass->Fill(aux_mass);
                fHist_KaonZeroShort_CPAwrtPV->Fill(vertex.GetV0CosineOfPointingAngle(PV[0], PV[1], PV[2]));
                fHist_KaonZeroShort_DCAwrtPV->Fill(
                    Calculate_LinePointDCA(vertex.Px(), vertex.Py(), vertex.Pz(), V0_X, V0_Y, V0_Z, PV[0], PV[1], PV[2]));
            }

            esdV0s.push_back(vertex);
            aux_pair.first = nidx;
            aux_pair.second = pidx;
            idxDaughters.push_back(aux_pair);

            if (std::any_of(RecDaughters_SignalV0.begin(), RecDaughters_SignalV0.end(),
                            [&](std::pair<Int_t, Int_t> signal_pair) { return signal_pair == aux_pair; })) {
                AliInfoF("%i %i %i", pdgV0, aux_pair.first, aux_pair.second);  // debug
            }
        }  // end of loop over pos. tracks
    }      // end of loop over neg. tracks
}

/*                          */
/**  V0s -- Kalman Filter  **/
/*** ==================== ***/

/*
 Apply cuts to an (anti-lambda) candidate, reconstructed from an (anti-proton, pi+) pair. Relevant for channel A and D.
 */
Bool_t AliAnalysisTaskSexaquark::PassesAntiLambdaCuts_KF(KFParticleMother kfV0, KFParticle kfNegDau, KFParticle kfPosDau,  //
                                                         TLorentzVector lvV0, TLorentzVector lvNegDau, TLorentzVector lvPosDau) {

    // define geometrical cuts
    Float_t kfDCA_tracks = TMath::Abs(kfNegDau.GetDistanceFromParticleXY(kfPosDau));
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

/*
 Find all true V0s which had both of their daughters reconstructed.
 - Uses: `fESD`, `fMC`, `fMagneticField`
 - Input: `idxNegativeTracks`, `idxPositiveTracks`, `pdgV0`, `pdgTrackNeg`, `pdgTrackPos`
 - Output: `kfV0s`, `idxDaughters`
*/
void AliAnalysisTaskSexaquark::ReconstructV0s_KF(std::vector<Int_t> idxNegativeTracks,               //
                                                 std::vector<Int_t> idxPositiveTracks,               //
                                                 Int_t pdgV0, Int_t pdgTrackNeg, Int_t pdgTrackPos,  //
                                                 std::vector<KFParticleMother>& kfV0s,               //
                                                 std::vector<std::pair<Int_t, Int_t>>& idxDaughters) {
    /* PENDING */
    return;
}

/*                                     */
/**  Sexaquark -- Geometrical Finder  **/
/*** =============================== ***/

/*
 Apply anti-sexaquark candidate cuts.
 - Input: `AntiLambda`, `KaonZeroShort`
 - Return: `kTRUE` if the candidate passes the cuts, `kFALSE` otherwise
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA_Geo(AliESDv0 AntiLambda, AliESDv0 KaonZeroShort) {

    TVector3 AntiLambda_momentum(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz());
    TVector3 KaonZeroShort_momentum(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz());

    TVector3 AntiLambda_vertex(AntiLambda.Xv(), AntiLambda.Yv(), AntiLambda.Zv());
    TVector3 KaonZeroShort_vertex(KaonZeroShort.Xv(), KaonZeroShort.Yv(), KaonZeroShort.Zv());

    Double_t dca;
    TVector3 cpa1, cpa2;

    dca = Calculate_TwoLinesDCA_v2(AntiLambda_momentum, AntiLambda_vertex, KaonZeroShort_momentum, KaonZeroShort_vertex, cpa1, cpa2);
    AliInfoF("dca cpa1 cpa2 = %f (%f, %f, %f) (%f, %f, %f)", dca, cpa1.X(), cpa1.Y(), cpa1.Z(), cpa2.X(), cpa2.Y(), cpa2.Z());
    fHist_AntiSexaquark_DCA->Fill(dca);

    return kTRUE;
}

/*
 Using two AliESDv0 objects, reconstruct the anti-sexaquark candidate.
 - Input: `esdAntiLambdas`, `esdKaonsZeroShort`
 - Output: ????
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelA_Geo(std::vector<AliESDv0> esdAntiLambdas, std::vector<AliESDv0> esdKaonsZeroShort,
                                                            std::vector<std::pair<Int_t, Int_t>> idxAntiLambdaDaughters,
                                                            std::vector<std::pair<Int_t, Int_t>> idxKaonZeroShortDaughters,
                                                            std::vector<Int_t>& idxAntiSexaquarkDaughters) {

    AliESDv0 AntiLambda;
    AliESDv0 KaonZeroShort;

    Float_t antilambda_energy;
    Float_t kaonzero_energy;

    TLorentzVector lvAntiLambda;
    TLorentzVector lvKaonZeroShort;

    TLorentzVector lvStruckNucleon(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());
    TLorentzVector lvAntiSexaquark;

    for (Int_t idxAntiLambda = 0; idxAntiLambda < (Int_t)esdAntiLambdas.size(); idxAntiLambda++) {
        for (Int_t idxKaonZeroShort = 0; idxKaonZeroShort < (Int_t)esdKaonsZeroShort.size(); idxKaonZeroShort++) {

            AntiLambda = esdAntiLambdas[idxAntiLambda];
            KaonZeroShort = esdKaonsZeroShort[idxKaonZeroShort];

            AliInfoF("idxAL idxKZS = %i %i", idxAntiLambda, idxKaonZeroShort);
            PassesSexaquarkCuts_ChannelA_Geo(AntiLambda, KaonZeroShort);

            // fill tlorentzvectors with the information of the v0s
            antilambda_energy = TMath::Sqrt(AntiLambda.Px() * AntiLambda.Px() + AntiLambda.Py() * AntiLambda.Py() +
                                            AntiLambda.Pz() * AntiLambda.Pz() + fPDG.GetParticle(-3122)->Mass() * fPDG.GetParticle(-3122)->Mass());
            kaonzero_energy = TMath::Sqrt(KaonZeroShort.Px() * KaonZeroShort.Px() + KaonZeroShort.Py() * KaonZeroShort.Py() +
                                          KaonZeroShort.Pz() * KaonZeroShort.Pz() + fPDG.GetParticle(310)->Mass() * fPDG.GetParticle(310)->Mass());
            lvAntiLambda.SetPxPyPzE(AntiLambda.Px(), AntiLambda.Py(), AntiLambda.Pz(), antilambda_energy);
            lvKaonZeroShort.SetPxPyPzE(KaonZeroShort.Px(), KaonZeroShort.Py(), KaonZeroShort.Pz(), kaonzero_energy);

            // calculate the invariant mass of the anti-sexaquark candidate
            lvAntiSexaquark = lvAntiLambda + lvKaonZeroShort - lvStruckNucleon;

            fHist_AntiSexaquark_Mass->Fill(lvAntiSexaquark.M());
        }
    }

    return;
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/*
 Apply anti-sexaquark candidate cuts.
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelA_KF(KFParticleMother kfAntiSexaquark, KFParticle kfAntiLambda,
                                                                 KFParticle kfKaonZeroShort,  //
                                                                 TLorentzVector lvAntiSexaquark, TLorentzVector lvAntiLambda,
                                                                 TLorentzVector lvKaonZeroShort) {
    /*
        Float_t kfDCA_Lambda_Pion2 = TMath::Abs(kfLambda.GetDistanceFromParticleXY(kfPion2));
        Float_t kfDCA_Lambda_Pion3 = TMath::Abs(kfLambda.GetDistanceFromParticleXY(kfPion3));

        Float_t kfDCA_Pion2_SV = TMath::Abs(kfPion2.GetDistanceFromVertexXY(kfLambda1520));

        Float_t PrimVertex[3] = {(Float_t)fPrimaryVertex->GetX(), (Float_t)fPrimaryVertex->GetY(), (Float_t)fPrimaryVertex->GetZ()};
        Float_t kfDCA_Lambda1520_PV = TMath::Abs(kfLambda1520.GetDistanceFromVertexXY(PrimVertex));

        Float_t kfCPA_Lambda_SV = CosinePointingAngle(lvLambda, kfLambda.GetX(), kfLambda.GetY(), kfLambda.GetZ(),  //
                                                      kfLambda1520.GetX(), kfLambda1520.GetY(), kfLambda1520.GetZ());

        Float_t mass = lvLambda1520.M();

        if (kfDCA_Lambda_Pion2 > kMaxDCATracks) return kFALSE;
        if (kfDCA_Lambda_Pion3 > kMaxDCATracks) return kFALSE;
        if (kfDCA_Pion2_SV > kMaxDCATracks) return kFALSE;  // (pending) should assign to another variable, later, probably
        if (kfCPA_Lambda_SV < kMinCPA_Lambda_SV) return kFALSE;
        if (mass > 1.62 || mass < 1.42) return kFALSE;  // (pending) refactorize
        if (kfDCA_Lambda1520_PV > 0.5) return kFALSE;   // (pending) refactorize
     */
    /* PENDING */
    return kTRUE;
}

/*
 Using Kalman Filter, reconstruct anti-sexaquark candidates from the reaction channel A.
 - Uses: `fStruckNucleonPDG`, `fFinalStateProductsPDG`
 - Input: `kfFirstV0s`, `idxFirstV0Daughters`, `kfSecondV0s`, `idxSecondV0Daughters`
 - Output: `kfAntiSexaquarks`
*/
void AliAnalysisTaskSexaquark::SexaquarkFinder_ChannelA_KF(std::vector<KFParticleMother> kfAntiLambdas,                     //
                                                           std::vector<std::pair<Int_t, Int_t>> idxAntiLambdaDaughters,     //
                                                           std::vector<KFParticleMother> kfKaonsZeroShort,                  //
                                                           std::vector<std::pair<Int_t, Int_t>> idxKaonZeroShortDaughters,  //
                                                           std::vector<KFParticleMother>& kfAntiSexaquarks) {

    Double_t impar[2];
    Bool_t fStatusT0PropToDCA;
    Bool_t fStatusT1PropToDCA;
    Bool_t fStatusT2PropToDCA;
    Bool_t fStatusT3PropToDCA;

    Int_t idxV0A_NegDau, idxV0A_PosDau;
    Int_t idxV0B_NegDau, idxV0B_PosDau;

    TLorentzVector lvNucleon;
    TLorentzVector lvTrack0, lvTrack1, lvTrack2, lvTrack3;
    TLorentzVector lvAntiSexaquark;

    AliInfo("               Mass       Vx       Vy       Vz     Chi2");
    for (Int_t idxV0A = 0; idxV0A < (Int_t)kfAntiLambdas.size(); idxV0A++) {
        for (Int_t idxV0B = 0; idxV0B < (Int_t)kfKaonsZeroShort.size(); idxV0B++) {

            idxV0A_NegDau = idxAntiLambdaDaughters[idxV0A].first;
            idxV0A_PosDau = idxAntiLambdaDaughters[idxV0A].first;
            idxV0B_NegDau = idxKaonZeroShortDaughters[idxV0B].second;
            idxV0B_PosDau = idxKaonZeroShortDaughters[idxV0B].second;

            AliInfoF("%i %i %i %i", idxV0A_NegDau, idxV0A_PosDau, idxV0B_NegDau, idxV0B_PosDau);

            // (protection) avoid repetition of daughters
            if (idxV0A_NegDau == idxV0B_NegDau || idxV0A_PosDau == idxV0B_PosDau) continue;

            AliESDtrack* esdTrack0 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_NegDau));
            AliESDtrack* esdTrack1 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0A_PosDau));
            AliESDtrack* esdTrack2 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_NegDau));
            AliESDtrack* esdTrack3 = static_cast<AliESDtrack*>(fESD->GetTrack(idxV0B_PosDau));

            /* Kalman Filter */

            KFParticleMother kfV0A = kfAntiLambdas[idxV0A];
            KFParticleMother kfV0B = kfKaonsZeroShort[idxV0B];

            // KFParticle kfTrack0 = CreateKFParticle(*esdTrack0, fFinalStateProductsPDG[0], -1);  // not used
            // KFParticle kfTrack1 = CreateKFParticle(*esdTrack1, fFinalStateProductsPDG[1], 1);   // not used
            KFParticle kfTrack2 = CreateKFParticle(*esdTrack2, fFinalStateProductsPDG[2], -1);
            KFParticle kfTrack3 = CreateKFParticle(*esdTrack3, fFinalStateProductsPDG[3], 1);

            KFParticleMother kfAntiSexaquark;
            kfAntiSexaquark.SetConstructMethod(2);
            kfAntiSexaquark.AddDaughter(kfV0A);
            kfAntiSexaquark.AddDaughter(kfV0B);
            kfAntiSexaquark.TransportToDecayVertex();

            // clone ESD tracks
            AliESDtrack cloneTrack0(*esdTrack0);
            AliESDtrack cloneTrack1(*esdTrack1);
            AliESDtrack cloneTrack2(*esdTrack2);
            AliESDtrack cloneTrack3(*esdTrack3);

            // get V0A vertex
            Double_t V0A_Vertex[3] = {kfV0A.GetX(), kfV0A.GetY(), kfV0A.GetZ()};
            Double_t V0A_VertexErr[3] = {kfV0A.GetErrX(), kfV0A.GetErrY(), kfV0A.GetErrZ()};
            AliESDVertex* esdV0A_Vertex = new AliESDVertex(V0A_Vertex, V0A_VertexErr);
            // propagate V0 daughters to their DCA w.r.t. V0 vertex
            fStatusT0PropToDCA = cloneTrack0.PropagateToDCA(esdV0A_Vertex, fMagneticField, 10., impar) == 1;
            fStatusT1PropToDCA = cloneTrack1.PropagateToDCA(esdV0A_Vertex, fMagneticField, 10., impar) == 1;

            // get V0B vertex
            Double_t V0B_Vertex[3] = {kfV0B.GetX(), kfV0B.GetY(), kfV0B.GetZ()};
            Double_t V0B_VertexErr[3] = {kfV0B.GetErrX(), kfV0B.GetErrY(), kfV0B.GetErrZ()};
            AliESDVertex* esdV0B_Vertex = new AliESDVertex(V0B_Vertex, V0B_VertexErr);
            // propagate V0 daughters to their DCA w.r.t. V0 vertex
            fStatusT2PropToDCA = cloneTrack0.PropagateToDCA(esdV0B_Vertex, fMagneticField, 10., impar) == 1;
            fStatusT3PropToDCA = cloneTrack1.PropagateToDCA(esdV0B_Vertex, fMagneticField, 10., impar) == 1;

            // get Secondary Vertex
            Double_t SecVertex[3] = {kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ()};
            Double_t SecVertexErr[3] = {kfAntiSexaquark.GetErrX(), kfAntiSexaquark.GetErrY(), kfAntiSexaquark.GetErrZ()};
            AliESDVertex* esdSecVertex = new AliESDVertex(SecVertex, SecVertexErr);

            // (question) should I find a way to propagate the V0s to the SV?

            /* Reconstruct AntiSexaquark */

            lvNucleon.SetXYZM(0., 0., 0., fPDG.GetParticle(fStruckNucleonPDG)->Mass());
            lvTrack0.SetXYZM(cloneTrack0.Px(), cloneTrack0.Py(), cloneTrack0.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[0])->Mass());
            lvTrack1.SetXYZM(cloneTrack1.Px(), cloneTrack1.Py(), cloneTrack1.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[1])->Mass());
            lvTrack2.SetXYZM(cloneTrack2.Px(), cloneTrack2.Py(), cloneTrack2.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[2])->Mass());
            lvTrack3.SetXYZM(cloneTrack3.Px(), cloneTrack3.Py(), cloneTrack3.Pz(), fPDG.GetParticle(fFinalStateProductsPDG[3])->Mass());

            lvAntiSexaquark = lvTrack0 + lvTrack1 + lvTrack2 + lvTrack3 - lvNucleon;

            // (cut) (pending) DCA w.r.t. PV
            // dca = Calculate_LinePointDCA(lvAntiSexaquark.Px(), lvAntiSexaquark.Py(), lvAntiSexaquark.Pz(),        //
            //  kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(), kfAntiSexaquark.GetZ(),  //
            //  fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
            // // if (dca > 0.1) continue;

            AliInfoF("   -     - %8.3f %8.3f %8.3f %8.3f %8.3f", lvAntiSexaquark.M(), kfAntiSexaquark.GetX(), kfAntiSexaquark.GetY(),
                     kfAntiSexaquark.GetZ(), kfAntiSexaquark.Chi2() / (Double_t)kfAntiSexaquark.GetNDF());

            kfAntiSexaquarks.push_back(kfAntiSexaquark);
        }
    }
}

/*
 Apply anti-sexaquark candidate cuts.
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelD_KF() {
    /* PENDING */
    return kTRUE;
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
 Apply anti-sexaquark candidate cuts.
*/
Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_ChannelE_KF() {
    /* PENDING */
    return kTRUE;
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

/*** Mathematical Functions ***/

/*
 Calculate the cosine of the pointing angle
 of a particle with momentum Px,Py,Pz and vertex X,Y,Z w.r.t. to a reference point
 [based on AliRoot/STEER/ESD/AliESDv0::GetV0CosineOfPointingAngle()]
*/
Double_t AliAnalysisTaskSexaquark::Calculate_CPA(Double_t Px, Double_t Py, Double_t Pz,  //
                                                 Double_t X, Double_t Y, Double_t Z,     //
                                                 Double_t refPointX, Double_t refPointY, Double_t refPointZ) {

    Double_t deltaPos[3];
    deltaPos[0] = X - refPointX;
    deltaPos[1] = Y - refPointY;
    deltaPos[2] = Z - refPointZ;

    Double_t mom2 = Px * Px + Py * Py + Pz * Pz;
    Double_t deltaPos2 = deltaPos[0] * deltaPos[0] + deltaPos[1] * deltaPos[1] + deltaPos[2] * deltaPos[2];

    // (protection)
    if (deltaPos2 == 0. || mom2 == 0.) {
        return -2.;
    }

    return (deltaPos[0] * Px + deltaPos[1] * Py + deltaPos[2] * Pz) / TMath::Sqrt(mom2 * deltaPos2);
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
*/
Double_t AliAnalysisTaskSexaquark::Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                                          Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                                          Double_t PV_X, Double_t PV_Y, Double_t PV_Z) {

    TVector3 V0Momentum(V0_Px, V0_Py, V0_Pz);
    TVector3 V0Vertex(V0_X, V0_Y, V0_Z);
    TVector3 PrimVertex(PV_X, PV_Y, PV_Z);

    TVector3 CrossProduct = (PrimVertex - V0Vertex).Cross(V0Momentum);

    return CrossProduct.Mag() / V0Momentum.Mag();
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
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v1(TVector3 pos1, TVector3 dir1, TVector3 pos2, TVector3 dir2, TVector3& P1, TVector3& P2) {

    if (dir1.Mag() == 0. || dir2.Mag() == 0.) {
        AliWarning("One of the directions is null!");
        return -1.;
    }

    /* Require perpendicularity, that means both lines must intersect eventually */

    if (dir1.Cross(dir2).Mag() <= 0.) {
        AliWarning("TwoLinesDCA :: Not possible to intersect!");
        return -1.;
    }

    TVector3 diff = pos1 - pos2;
    Double_t determinant = dir1.Dot(dir2) * dir1.Dot(dir2) - dir1.Mag2() * dir2.Mag2();
    if (determinant == 0.) {
        AliWarning("TwoLinesDCA :: Lines are parallel!");
        return -1.;
    }
    Double_t sol1 = (dir2.Mag2() * dir1.Dot(diff) - dir2.Dot(diff) * dir2.Dot(dir1)) / determinant;
    Double_t sol2 = (-1 * dir1.Mag2() * dir2.Dot(diff) + dir1.Dot(diff) * dir2.Dot(dir1)) / determinant;

    // endpoints of the closest line
    P1 = pos1 + sol1 * dir1;
    P2 = pos2 + sol2 * dir2;

    // the distance of the closest line
    return (diff + dir1 * sol1 - dir2 * sol2).Mag();
}

/*
 Given the parametric form of two lines, find their closest distance.
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
Double_t AliAnalysisTaskSexaquark::Calculate_TwoLinesDCA_v2(TVector3 v3_pos0, TVector3 v3_dir0, TVector3 v3_pos1, TVector3 v3_dir1, TVector3& CPA0,
                                                            TVector3& CPA1) {

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
    CPA0.SetXYZ(dir0[0] * t0 + pos0[0], dir0[1] * t0 + pos0[1], dir0[2] * t0 + pos0[2]);
    CPA1.SetXYZ(dir1[0] * t1 + pos1[0], dir1[1] * t1 + pos1[1], dir1[2] * t1 + pos1[2]);

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