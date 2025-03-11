#include "AliAnalysisTaskSexaquark.h"

ClassImp(AliAnalysisTaskSexaquark);

/*
 * Empty I/O constructor. Non-persistent members are initialized to their default values from here.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark()
    : AliAnalysisTaskSE(),
      /*  */
      fIsMC(false),
      fIsSignalMC(false),
      /*  */
      fPIDResponse(nullptr),
      /*  */
      fMC(nullptr),
      fMC_PrimaryVertex(nullptr),
      /*  */
      fESD(nullptr),
      fPrimaryVertex(nullptr),
      v3_pv(),
      fEventCuts(),
      /*  */
      fAliEnPath(""),
      fSignalLog_NewBasename(""),
      /*  */
      fReactionID_(),
      fSexaquark_Px_(),
      fSexaquark_Py_(),
      fSexaquark_Pz_(),
      fNucleon_Px_(),
      fNucleon_Py_(),
      fNucleon_Pz_(),
      /*  */
      fOutputTree(nullptr),
      /*  */
      fRunNumber(0),
      fDirNumber(0),
      fDirNumberB(0),
      fEventNumber(0),
      fCentrality(0.),
      fMagneticField(0.),
      tEvents_MC_PV_Xv(0),
      tEvents_MC_PV_Yv(0),
      tEvents_MC_PV_Zv(0),
      tEvents_PV_Xv(0),
      tEvents_PV_Yv(0),
      tEvents_PV_Zv(0),
      tEvent_PV_CovMatrix(),
      tEvents_NTracks(0),
      tEvents_NTPCClusters(0),
      /*  */
      kf_pv(),
      /*  */
      tInjected_Nucleon_PdgCode(0),
      tInjected_Mass(0.),
      tInjected_ReactionID(),
      tInjected_Px(),
      tInjected_Py(),
      tInjected_Pz(),
      tInjected_Nucleon_Px(),
      tInjected_Nucleon_Py(),
      tInjected_Nucleon_Pz(),
      tInjected_Xv(),
      tInjected_Yv(),
      tInjected_Zv(),
      tInjected_Post_Px(),
      tInjected_Post_Py(),
      tInjected_Post_Pz(),
      /*  */
      tV0_Idx(),
      tV0_Px(),
      tV0_Py(),
      tV0_Pz(),
      tV0_E(),
      tV0_Xv(),
      tV0_Yv(),
      tV0_Zv(),
      tV0_CPAwrtPV(),
      tV0_DCAwrtPV(),
      tV0_ArmQt(),
      tV0_ArmAlpha(),
      tV0_DCA_Daughters(),
      tV0_Neg_EsdIdx(),
      tV0_Neg_Px(),
      tV0_Neg_Py(),
      tV0_Neg_Pz(),
      tV0_Neg_IPxy_V0(),
      tV0_Neg_IPz_V0(),
      tV0_Neg_TrackParamD_V0(),
      tV0_Neg_DCA_V0(),
      tV0_Neg_DCAxy_V0(),
      tV0_Pos_EsdIdx(),
      tV0_Pos_Px(),
      tV0_Pos_Py(),
      tV0_Pos_Pz(),
      tV0_Pos_IPxy_V0(),
      tV0_Pos_IPz_V0(),
      tV0_Pos_TrackParamD_V0(),
      tV0_Pos_DCA_V0(),
      tV0_Pos_DCAxy_V0(),
      /*  */
      tV0_McIdx(),
      tV0_PdgCode(),
      tV0_IsSignal(),
      tV0_ReactionID(),
      tV0_IsHybrid(),
      /*  */
      tTypeA_Px(),
      tTypeA_Py(),
      tTypeA_Pz(),
      tTypeA_E(),
      tTypeA_E_asDecay(),
      tTypeA_Xv(),
      tTypeA_Yv(),
      tTypeA_Zv(),
      tTypeA_V0a_Idx(),
      tTypeA_V0a_Px(),
      tTypeA_V0a_Py(),
      tTypeA_V0a_Pz(),
      tTypeA_V0a_E(),
      tTypeA_V0a_DecayLength(),
      tTypeA_DCAV0aSV(),
      tTypeA_DCAV0aNegSV(),
      tTypeA_DCAV0aPosSV(),
      tTypeA_V0b_Idx(),
      tTypeA_V0b_Px(),
      tTypeA_V0b_Py(),
      tTypeA_V0b_Pz(),
      tTypeA_V0b_E(),
      tTypeA_V0b_DecayLength(),
      tTypeA_DCAV0bSV(),
      tTypeA_DCAV0bNegSV(),
      tTypeA_DCAV0bPosSV(),
      tTypeA_DCAbtwV0s(),
      /*  */
      tTypeA_IsSignal(),
      tTypeA_ReactionID(),
      tTypeA_IsHybrid(),
      /*  */
      tTypeD_Px(),
      tTypeD_Py(),
      tTypeD_Pz(),
      tTypeD_E(),
      tTypeD_E_asDecay(),
      tTypeD_Xv(),
      tTypeD_Yv(),
      tTypeD_Zv(),
      tTypeD_V0_Idx(),
      tTypeD_V0_Px(),
      tTypeD_V0_Py(),
      tTypeD_V0_Pz(),
      tTypeD_V0_E(),
      tTypeD_V0_DecayLength(),
      tTypeD_DCAV0SV(),
      tTypeD_DCAV0NegSV(),
      tTypeD_DCAV0PosSV(),
      tTypeD_DCAV0NegKa(),
      tTypeD_DCAV0PosKa(),
      tTypeD_Ka_EsdIdx(),
      tTypeD_Ka_Px(),
      tTypeD_Ka_Py(),
      tTypeD_Ka_Pz(),
      tTypeD_DCAKaSV(),
      tTypeD_DCAKaV0(),
      /*  */
      tTypeD_IsSignal(),
      tTypeD_ReactionID(),
      tTypeD_IsHybrid(),
      /*  */
      fMC_PdgCode_(),
      fMC_Mother_McIdx_(),
      fMC_IsSignal_(),
      fMC_ReactionID_(),
      fReactionProducts_McIdx_(),
      fLinked_McIdx_(),
      /*  */
      fAntiProton_Indices(),
      fProton_Indices(),
      fNegKaon_Indices(),
      fPosKaon_Indices(),
      fPiMinus_Indices(),
      fPiPlus_Indices(),
      /*  */
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      mcAntiLambdas(),
      mcKaonsZeroShort() {}

/*
 * Constructor, called locally.
 */
AliAnalysisTaskSexaquark::AliAnalysisTaskSexaquark(const char* name)
    : AliAnalysisTaskSE(name),
      /*  */
      fIsMC(false),
      fIsSignalMC(false),
      /*  */
      fPIDResponse(nullptr),
      /*  */
      fMC(nullptr),
      fMC_PrimaryVertex(nullptr),
      /*  */
      fESD(nullptr),
      fPrimaryVertex(nullptr),
      v3_pv(),
      fEventCuts(),
      /*  */
      fAliEnPath(""),
      fSignalLog_NewBasename(""),
      /*  */
      fReactionID_(),
      fSexaquark_Px_(),
      fSexaquark_Py_(),
      fSexaquark_Pz_(),
      fNucleon_Px_(),
      fNucleon_Py_(),
      fNucleon_Pz_(),
      /*  */
      fOutputTree(nullptr),
      /*  */
      fRunNumber(0),
      fDirNumber(0),
      fDirNumberB(0),
      fEventNumber(0),
      fCentrality(0.),
      fMagneticField(0.),
      tEvents_MC_PV_Xv(0),
      tEvents_MC_PV_Yv(0),
      tEvents_MC_PV_Zv(0),
      tEvents_PV_Xv(0),
      tEvents_PV_Yv(0),
      tEvents_PV_Zv(0),
      tEvent_PV_CovMatrix(),
      tEvents_NTracks(0),
      tEvents_NTPCClusters(0),
      /*  */
      kf_pv(),
      /*  */
      tInjected_Nucleon_PdgCode(0),
      tInjected_Mass(0.),
      tInjected_ReactionID(),
      tInjected_Px(),
      tInjected_Py(),
      tInjected_Pz(),
      tInjected_Nucleon_Px(),
      tInjected_Nucleon_Py(),
      tInjected_Nucleon_Pz(),
      tInjected_Xv(),
      tInjected_Yv(),
      tInjected_Zv(),
      tInjected_Post_Px(),
      tInjected_Post_Py(),
      tInjected_Post_Pz(),
      /*  */
      tV0_Idx(),
      tV0_Px(),
      tV0_Py(),
      tV0_Pz(),
      tV0_E(),
      tV0_Xv(),
      tV0_Yv(),
      tV0_Zv(),
      tV0_CPAwrtPV(),
      tV0_DCAwrtPV(),
      tV0_ArmQt(),
      tV0_ArmAlpha(),
      tV0_DCA_Daughters(),
      tV0_Neg_EsdIdx(),
      tV0_Neg_Px(),
      tV0_Neg_Py(),
      tV0_Neg_Pz(),
      tV0_Neg_IPxy_V0(),
      tV0_Neg_IPz_V0(),
      tV0_Neg_TrackParamD_V0(),
      tV0_Neg_DCA_V0(),
      tV0_Neg_DCAxy_V0(),
      tV0_Pos_EsdIdx(),
      tV0_Pos_Px(),
      tV0_Pos_Py(),
      tV0_Pos_Pz(),
      tV0_Pos_IPxy_V0(),
      tV0_Pos_IPz_V0(),
      tV0_Pos_TrackParamD_V0(),
      tV0_Pos_DCA_V0(),
      tV0_Pos_DCAxy_V0(),
      /*  */
      tV0_McIdx(),
      tV0_PdgCode(),
      tV0_IsSignal(),
      tV0_ReactionID(),
      tV0_IsHybrid(),
      /*  */
      tTypeA_Px(),
      tTypeA_Py(),
      tTypeA_Pz(),
      tTypeA_E(),
      tTypeA_E_asDecay(),
      tTypeA_Xv(),
      tTypeA_Yv(),
      tTypeA_Zv(),
      tTypeA_V0a_Idx(),
      tTypeA_V0a_Px(),
      tTypeA_V0a_Py(),
      tTypeA_V0a_Pz(),
      tTypeA_V0a_E(),
      tTypeA_V0a_DecayLength(),
      tTypeA_DCAV0aSV(),
      tTypeA_DCAV0aNegSV(),
      tTypeA_DCAV0aPosSV(),
      tTypeA_V0b_Idx(),
      tTypeA_V0b_Px(),
      tTypeA_V0b_Py(),
      tTypeA_V0b_Pz(),
      tTypeA_V0b_E(),
      tTypeA_V0b_DecayLength(),
      tTypeA_DCAV0bSV(),
      tTypeA_DCAV0bNegSV(),
      tTypeA_DCAV0bPosSV(),
      tTypeA_DCAbtwV0s(),
      /*  */
      tTypeA_IsSignal(),
      tTypeA_ReactionID(),
      tTypeA_IsHybrid(),
      /*  */
      tTypeD_Px(),
      tTypeD_Py(),
      tTypeD_Pz(),
      tTypeD_E(),
      tTypeD_E_asDecay(),
      tTypeD_Xv(),
      tTypeD_Yv(),
      tTypeD_Zv(),
      tTypeD_V0_Idx(),
      tTypeD_V0_Px(),
      tTypeD_V0_Py(),
      tTypeD_V0_Pz(),
      tTypeD_V0_E(),
      tTypeD_V0_DecayLength(),
      tTypeD_DCAV0SV(),
      tTypeD_DCAV0NegSV(),
      tTypeD_DCAV0PosSV(),
      tTypeD_DCAV0NegKa(),
      tTypeD_DCAV0PosKa(),
      tTypeD_Ka_EsdIdx(),
      tTypeD_Ka_Px(),
      tTypeD_Ka_Py(),
      tTypeD_Ka_Pz(),
      tTypeD_DCAKaSV(),
      tTypeD_DCAKaV0(),
      /*  */
      tTypeD_IsSignal(),
      tTypeD_ReactionID(),
      tTypeD_IsHybrid(),
      /*  */
      fMC_PdgCode_(),
      fMC_Mother_McIdx_(),
      fMC_IsSignal_(),
      fMC_ReactionID_(),
      fReactionProducts_McIdx_(),
      fLinked_McIdx_(),
      /*  */
      fAntiProton_Indices(),
      fProton_Indices(),
      fNegKaon_Indices(),
      fPosKaon_Indices(),
      fPiMinus_Indices(),
      fPiPlus_Indices(),
      /*  */
      kfAntiLambdas(),
      kfKaonsZeroShort(),
      mcAntiLambdas(),
      mcKaonsZeroShort() {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TTree::Class());  // fOutputTree
}

AliAnalysisTaskSexaquark::~AliAnalysisTaskSexaquark() { delete fOutputTree; }

/*
 * Initialize analysis task. Needs to be called within the `AddTaskSexaquark.C` macro.
 */
void AliAnalysisTaskSexaquark::Initialize(Bool_t is_mc, Bool_t is_signal_mc) {
    //
    fIsMC = is_mc;
    fIsSignalMC = is_signal_mc;
    //
    AliInfo("Settings:");
    AliInfoF(">> IsMC       = %i", (Int_t)fIsMC);
    AliInfoF(">> IsSignalMC = %i", (Int_t)fIsSignalMC);
}

void AliAnalysisTaskSexaquark::PrintCuts() {
    /*  */
    AliInfoF("Track::Min_Pt                  = %f", SexaCuts::Track::Min_Pt);
    AliInfoF("Track::Max_Pt                  = %f", SexaCuts::Track::Max_Pt);
    AliInfoF("Track::AbsMax_PID_NSigma       = %f", SexaCuts::Track::AbsMax_PID_NSigma);
    AliInfoF("Track::AbsMax_Eta              = %f", SexaCuts::Track::AbsMax_Eta);
    AliInfoF("Track::Min_NTPCClusters        = %hu", SexaCuts::Track::Min_NTPCClusters);
    AliInfoF("Track::Max_Chi2PerNTPCClusters = %f", SexaCuts::Track::Max_Chi2PerNTPCClusters);
    AliInfoF("Track::TurnedOn_StatusCuts     = %i", SexaCuts::Track::TurnedOn_StatusCuts);
    AliInfoF("Track::TurnedOn_RejectKinks    = %i", SexaCuts::Track::TurnedOn_RejectKinks);
    AliInfoF("Track::AbsMin_DCAxy_wrtPV      = %f", SexaCuts::Track::AbsMin_DCAxy_wrtPV);
    /*  */
    AliInfoF("Lambda::Min_Pt                = %f", SexaCuts::Lambda::Min_Pt);
    AliInfoF("Lambda::Min_Mass              = %f", SexaCuts::Lambda::Min_Mass);
    AliInfoF("Lambda::Max_Mass              = %f", SexaCuts::Lambda::Max_Mass);
    AliInfoF("Lambda::AbsMax_Eta            = %f", SexaCuts::Lambda::AbsMax_Eta);
    AliInfoF("Lambda::Min_CPAwrtPV          = %f", SexaCuts::Lambda::Min_CPAwrtPV);
    AliInfoF("Lambda::Max_CPAwrtPV          = %f", SexaCuts::Lambda::Max_CPAwrtPV);
    AliInfoF("Lambda::Min_DCAwrtPV          = %f", SexaCuts::Lambda::Min_DCAwrtPV);
    AliInfoF("Lambda::AbsMax_ArmQtOverAlpha = %f", SexaCuts::Lambda::AbsMax_ArmQtOverAlpha);
    AliInfoF("Lambda::AbsMax_Zv             = %f", SexaCuts::Lambda::AbsMax_Zv);
    AliInfoF("Lambda::Min_Radius            = %f", SexaCuts::Lambda::Min_Radius);
    AliInfoF("Lambda::Max_Radius            = %f", SexaCuts::Lambda::Max_Radius);
    AliInfoF("Lambda::Max_DCAbtwDau         = %f", SexaCuts::Lambda::Max_DCAbtwDau);
    AliInfoF("Lambda::Max_DCAnegV0          = %f", SexaCuts::Lambda::Max_DCAnegV0);
    AliInfoF("Lambda::Max_DCAposV0          = %f", SexaCuts::Lambda::Max_DCAposV0);
    /*  */
    AliInfoF("KaonZeroShort::Min_Pt        = %f", SexaCuts::KaonZeroShort::Min_Pt);
    AliInfoF("KaonZeroShort::Min_Mass      = %f", SexaCuts::KaonZeroShort::Min_Mass);
    AliInfoF("KaonZeroShort::Max_Mass      = %f", SexaCuts::KaonZeroShort::Max_Mass);
    AliInfoF("KaonZeroShort::AbsMax_Eta    = %f", SexaCuts::KaonZeroShort::AbsMax_Eta);
    AliInfoF("KaonZeroShort::Min_CPAwrtPV  = %f", SexaCuts::KaonZeroShort::Min_CPAwrtPV);
    AliInfoF("KaonZeroShort::Max_CPAwrtPV  = %f", SexaCuts::KaonZeroShort::Max_CPAwrtPV);
    AliInfoF("KaonZeroShort::Min_DCAwrtPV  = %f", SexaCuts::KaonZeroShort::Min_DCAwrtPV);
    AliInfoF("KaonZeroShort::AbsMax_Zv     = %f", SexaCuts::KaonZeroShort::AbsMax_Zv);
    AliInfoF("KaonZeroShort::Min_Radius    = %f", SexaCuts::KaonZeroShort::Min_Radius);
    AliInfoF("KaonZeroShort::Max_Radius    = %f", SexaCuts::KaonZeroShort::Max_Radius);
    AliInfoF("KaonZeroShort::Max_DCAbtwDau = %f", SexaCuts::KaonZeroShort::Max_DCAbtwDau);
    AliInfoF("KaonZeroShort::Max_DCAnegV0  = %f", SexaCuts::KaonZeroShort::Max_DCAnegV0);
    AliInfoF("KaonZeroShort::Max_DCAposV0  = %f", SexaCuts::KaonZeroShort::Max_DCAposV0);
}

/*
 * Create output objects, called once at RUNTIME ~ execution on Grid.
 */
void AliAnalysisTaskSexaquark::UserCreateOutputObjects() {
    //
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (man == nullptr) AliFatal("AliAnalysisManager couldn't be found.");
    auto* inputHandler = dynamic_cast<AliESDInputHandler*>(man->GetInputEventHandler());
    if (inputHandler == nullptr) AliFatal("AliESDInputHandler couldn't be found.");
    /* Add mandatory routines */
    fPIDResponse = inputHandler->GetPIDResponse();
    /* Debug */
    PrintCuts();
    /* Prepare output tree */
    fOutputTree = new TTree("Events", "Events");
    AssociateBranches_Events();
    if (fIsMC && fIsSignalMC) AssociateBranches_Injected();
    AssociateBranches_V0s();
    // AssociateBranches_TypeA();
    // AssociateBranches_TypeD();
    /* Post data */
    PostData(1, fOutputTree);
}

/*
 * User implementation of `Notify()`. Needed for reading the AliEn path.
 * This function is loaded during `AliAnalysisManager::Notify()`.
 * It's called after `UserCreateOutputObjects()`, for each new file, and before the first `UserExec()`.
 */
Bool_t AliAnalysisTaskSexaquark::UserNotify() {
    //
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) AliFatal("Analysis Manager not found");
    TTree* man_tree = man->GetTree();
    if (!man_tree) AliFatal("Analysis Manager Tree not found");
    TFile* man_file = man_tree->GetCurrentFile();
    if (!man_file) AliFatal("Analysis Manager File not found");
    /* get AliEn path and tokenize it */
    fAliEnPath = man_file->GetName();
    if (fAliEnPath == "") AliWarning("fAliEnPath couldn't be found.");
    AliInfoF("AliEn Path : %s", fAliEnPath.Data());
    /*  */
    TObjArray* tokens = fAliEnPath.Tokenize("/");
    if (fIsMC) {
        /* path of general purpose MC ends as `.../LHC20e3a/297595/001/AliESDs.root` */
        fDirNumber = (dynamic_cast<TObjString*>(tokens->At(tokens->GetEntries() - 2)))->GetString().Atof();
        AliInfoF("Dir Number : %04i", (Int_t)fDirNumber);
        if (fIsSignalMC) {
            /* path of signal MC ends as `.../LHC23l1a3/A1.73/297595/001/AliESDs.root` */
            TString SimSet = (dynamic_cast<TObjString*>(tokens->At(tokens->GetEntries() - 4)))->GetString();
            AliInfoF("Simulation Set : %s", SimSet.Data());
            if (SimSet[0] == 'A')
                tInjected_Nucleon_PdgCode = PdgCode::Neutron;
            else
                tInjected_Nucleon_PdgCode = PdgCode::Proton;
            tInjected_Mass = ((TString)SimSet(1, 4)).Atof();
            ClearSignalLogs();
            BringSignalLogs();
            LoadSignalLogs();
        }
    } else {  // data
        /* path of data ends as `.../LHC15o/000245232/pass2/15000245232039.914/AliESDs.root` */
        TString aux_dir_nr = (dynamic_cast<TObjString*>(tokens->At(tokens->GetEntries() - 2)))->GetString();
        aux_dir_nr = TString(aux_dir_nr(2 + 3 + 6, 10));  // = "039.914"
        fDirNumber = TString(aux_dir_nr(0, 3)).Atoi();    // = 39
        fDirNumberB = TString(aux_dir_nr(4, 5)).Atoi();   // = 914
        AliInfoF("Dir Number (String): %s", aux_dir_nr.Data());
        AliInfoF("Dir Number A : %i", fDirNumber);
        AliInfoF("Dir Number B : %i", fDirNumberB);
    }
    return kTRUE;
}

/*
 * Main function, called per each event at RUNTIME ~ execution on Grid.
 */
void AliAnalysisTaskSexaquark::UserExec(Option_t* option) {
    /* Events */
    if (!ProcessEvent()) return;
    /* Set global properties */
    KFParticle::SetField(fMagneticField);
    kf_pv = CreateKFVertex(fPrimaryVertex);
    v3_pv.SetCoordinates(fPrimaryVertex->GetX(), fPrimaryVertex->GetY(), fPrimaryVertex->GetZ());
    /* MC Particles */
    if (fIsMC) {
        ProcessMCParticles();
        if (fIsSignalMC) ProcessInjected();
    }
    /* Tracks */
    ProcessTracks();
    /* V0s */
    KF_FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
    KF_FindV0s(PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);
    /* Sexaquarks */
    /* -- `AntiSexaquark,Neutron -> AntiLambda,K0S` */
    // KF_FindSexaquarks_TypeA(PdgCode::Neutron, {PdgCode::AntiLambda, PdgCode::KaonZeroShort});
    /* -- `AntiSexaquark,Proton -> AntiLambda,K+,(pi-,pi+)` */
    // KF_FindSexaquarks_TypeD(PdgCode::AntiNeutron, {PdgCode::Lambda, PdgCode::KaonZeroShort});
    /* Fill tree */
    fOutputTree->Fill();
    /* End of event */
    if (fIsMC && fIsSignalMC) ClearBranches_Injected();
    ClearBranches_V0s();
    // ClearBranches_TypeA();
    // ClearBranches_TypeD();
    ClearContainers();
    PostData(1, fOutputTree);
}

/*          */
/**  Tree  **/
/*** ==== ***/

void AliAnalysisTaskSexaquark::AssociateBranches_Events() {
    //
    fOutputTree->Branch("RunNumber", &fRunNumber);
    fOutputTree->Branch("DirNumber", &fDirNumber);
    if (!fIsMC) fOutputTree->Branch("DirNumberB", &fDirNumberB);
    fOutputTree->Branch("EventNumber", &fEventNumber);
    fOutputTree->Branch("Centrality", &fCentrality);
    fOutputTree->Branch("MagneticField", &fMagneticField);
    if (fIsMC) {
        fOutputTree->Branch("MC_PV_Xv", &tEvents_MC_PV_Xv);
        fOutputTree->Branch("MC_PV_Yv", &tEvents_MC_PV_Yv);
        fOutputTree->Branch("MC_PV_Zv", &tEvents_MC_PV_Zv);
    }
    fOutputTree->Branch("PV_Xv", &tEvents_PV_Xv);
    fOutputTree->Branch("PV_Yv", &tEvents_PV_Yv);
    fOutputTree->Branch("PV_Zv", &tEvents_PV_Zv);
    fOutputTree->Branch("PV_CovMatrix", &tEvent_PV_CovMatrix);
    fOutputTree->Branch("NTracks", &tEvents_NTracks);
    fOutputTree->Branch("NTPCClusters", &tEvents_NTPCClusters);
}

void AliAnalysisTaskSexaquark::AssociateBranches_Injected() {
    fOutputTree->Branch("Nucleon_PdgCode", &tInjected_Nucleon_PdgCode);
    fOutputTree->Branch("Injected_Mass", &tInjected_Mass);
    /*  */
    fOutputTree->Branch("ReactionID", &tInjected_ReactionID);
    fOutputTree->Branch("Injected_Px", &tInjected_Px);
    fOutputTree->Branch("Injected_Py", &tInjected_Py);
    fOutputTree->Branch("Injected_Pz", &tInjected_Pz);
    fOutputTree->Branch("Fermi_Px", &tInjected_Nucleon_Px);
    fOutputTree->Branch("Fermi_Py", &tInjected_Nucleon_Py);
    fOutputTree->Branch("Fermi_Pz", &tInjected_Nucleon_Pz);
    fOutputTree->Branch("Injected_Xv", &tInjected_Xv);
    fOutputTree->Branch("Injected_Yv", &tInjected_Yv);
    fOutputTree->Branch("Injected_Zv", &tInjected_Zv);
    fOutputTree->Branch("Injected_Post_Px", &tInjected_Post_Px);
    fOutputTree->Branch("Injected_Post_Py", &tInjected_Post_Py);
    fOutputTree->Branch("Injected_Post_Pz", &tInjected_Post_Pz);
}

void AliAnalysisTaskSexaquark::AssociateBranches_V0s() {
    //
    std::vector<Short_t> v0_pdg_codes = {PdgCode::AntiLambda, PdgCode::KaonZeroShort};
    std::vector<TString> v0_names = {"AL", "K0S"};
    for (size_t i = 0; i < v0_pdg_codes.size(); i++) {
        fOutputTree->Branch(v0_names[i] + "_Idx", &tV0_Idx[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Px", &tV0_Px[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Py", &tV0_Py[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pz", &tV0_Pz[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_E", &tV0_E[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Xv", &tV0_Xv[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Yv", &tV0_Yv[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Zv", &tV0_Zv[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_CPAwrtPV", &tV0_CPAwrtPV[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_DCAwrtPV", &tV0_DCAwrtPV[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_ArmQt", &tV0_ArmQt[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_ArmAlpha", &tV0_ArmAlpha[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_DCA_Daughters", &tV0_DCA_Daughters[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_EsdIdx", &tV0_Neg_EsdIdx[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Px", &tV0_Neg_Px[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Py", &tV0_Neg_Py[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Pz", &tV0_Neg_Pz[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_IPxy_V0", &tV0_Neg_IPxy_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_IPz_V0", &tV0_Neg_IPz_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_TrackParamD_V0", &tV0_Neg_TrackParamD_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_DCA_V0", &tV0_Neg_DCA_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_DCAxy_V0", &tV0_Neg_DCAxy_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_EsdIdx", &tV0_Pos_EsdIdx[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Px", &tV0_Pos_Px[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Py", &tV0_Pos_Py[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Pz", &tV0_Pos_Pz[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_IPxy_V0", &tV0_Pos_IPxy_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_IPz_V0", &tV0_Pos_IPz_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_TrackParamD_V0", &tV0_Pos_TrackParamD_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_DCA_V0", &tV0_Pos_DCA_V0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_DCAxy_V0", &tV0_Pos_DCAxy_V0[v0_pdg_codes[i]]);
        if (fIsMC) {
            fOutputTree->Branch(v0_names[i] + "_McIdx", &tV0_McIdx[v0_pdg_codes[i]]);
            fOutputTree->Branch(v0_names[i] + "_PdgCode", &tV0_PdgCode[v0_pdg_codes[i]]);
            fOutputTree->Branch(v0_names[i] + "_IsSignal", &tV0_IsSignal[v0_pdg_codes[i]]);
            fOutputTree->Branch(v0_names[i] + "_ReactionID", &tV0_ReactionID[v0_pdg_codes[i]]);
            fOutputTree->Branch(v0_names[i] + "_IsHybrid", &tV0_IsHybrid[v0_pdg_codes[i]]);
        }
    }
}

void AliAnalysisTaskSexaquark::AssociateBranches_TypeA() {
    //
    fOutputTree->Branch("ASA_Px", &tTypeA_Px);
    fOutputTree->Branch("ASA_Py", &tTypeA_Py);
    fOutputTree->Branch("ASA_Pz", &tTypeA_Pz);
    fOutputTree->Branch("ASA_E", &tTypeA_E);
    fOutputTree->Branch("ASA_E_asDecay", &tTypeA_E_asDecay);
    fOutputTree->Branch("ASA_Xv", &tTypeA_Xv);
    fOutputTree->Branch("ASA_Yv", &tTypeA_Yv);
    fOutputTree->Branch("ASA_Zv", &tTypeA_Zv);
    fOutputTree->Branch("ASA_V0a_Idx", &tTypeA_V0a_Idx);
    fOutputTree->Branch("ASA_V0a_Px", &tTypeA_V0a_Px);
    fOutputTree->Branch("ASA_V0a_Py", &tTypeA_V0a_Py);
    fOutputTree->Branch("ASA_V0a_Pz", &tTypeA_V0a_Pz);
    fOutputTree->Branch("ASA_V0a_E", &tTypeA_V0a_E);
    fOutputTree->Branch("ASA_V0a_DecayLength", &tTypeA_V0a_DecayLength);
    fOutputTree->Branch("ASA_DCAV0aSV", &tTypeA_DCAV0aSV);
    fOutputTree->Branch("ASA_DCAV0aNegSV", &tTypeA_DCAV0aNegSV);
    fOutputTree->Branch("ASA_DCAV0aPosSV", &tTypeA_DCAV0aPosSV);
    fOutputTree->Branch("ASA_V0b_Idx", &tTypeA_V0b_Idx);
    fOutputTree->Branch("ASA_V0b_Px", &tTypeA_V0b_Px);
    fOutputTree->Branch("ASA_V0b_Py", &tTypeA_V0b_Py);
    fOutputTree->Branch("ASA_V0b_Pz", &tTypeA_V0b_Pz);
    fOutputTree->Branch("ASA_V0b_E", &tTypeA_V0b_E);
    fOutputTree->Branch("ASA_V0b_DecayLength", &tTypeA_V0b_DecayLength);
    fOutputTree->Branch("ASA_DCAV0bSV", &tTypeA_DCAV0bSV);
    fOutputTree->Branch("ASA_DCAV0bNegSV", &tTypeA_DCAV0bNegSV);
    fOutputTree->Branch("ASA_DCAV0bPosSV", &tTypeA_DCAV0bPosSV);
    fOutputTree->Branch("ASA_DCAbtwV0s", &tTypeA_DCAbtwV0s);
    if (fIsMC) {
        fOutputTree->Branch("ASA_IsSignal", &tTypeA_IsSignal);
        fOutputTree->Branch("ASA_ReactionID", &tTypeA_ReactionID);
        fOutputTree->Branch("ASA_IsHybrid", &tTypeA_IsHybrid);
    }
}

void AliAnalysisTaskSexaquark::AssociateBranches_TypeD() {
    fOutputTree->Branch("ASD_Px", &tTypeD_Px);
    fOutputTree->Branch("ASD_Py", &tTypeD_Py);
    fOutputTree->Branch("ASD_Pz", &tTypeD_Pz);
    fOutputTree->Branch("ASD_E", &tTypeD_E);
    fOutputTree->Branch("ASD_E_asDecay", &tTypeD_E_asDecay);
    fOutputTree->Branch("ASD_Xv", &tTypeD_Xv);
    fOutputTree->Branch("ASD_Yv", &tTypeD_Yv);
    fOutputTree->Branch("ASD_Zv", &tTypeD_Zv);
    fOutputTree->Branch("ASD_V0_Idx", &tTypeD_V0_Idx);
    fOutputTree->Branch("ASD_V0_Px", &tTypeD_V0_Px);
    fOutputTree->Branch("ASD_V0_Py", &tTypeD_V0_Py);
    fOutputTree->Branch("ASD_V0_Pz", &tTypeD_V0_Pz);
    fOutputTree->Branch("ASD_V0_E", &tTypeD_V0_E);
    fOutputTree->Branch("ASD_V0_DecayLength", &tTypeD_V0_DecayLength);
    fOutputTree->Branch("ASD_DCAV0SV", &tTypeD_DCAV0SV);
    fOutputTree->Branch("ASD_DCAV0NegSV", &tTypeD_DCAV0NegSV);
    fOutputTree->Branch("ASD_DCAV0PosSV", &tTypeD_DCAV0PosSV);
    fOutputTree->Branch("ASD_DCAV0NegKa", &tTypeD_DCAV0NegKa);
    fOutputTree->Branch("ASD_DCAV0PosKa", &tTypeD_DCAV0PosKa);
    fOutputTree->Branch("ASD_Ka_EsdIdx", &tTypeD_Ka_EsdIdx);
    fOutputTree->Branch("ASD_Ka_Px", &tTypeD_Ka_Px);
    fOutputTree->Branch("ASD_Ka_Py", &tTypeD_Ka_Py);
    fOutputTree->Branch("ASD_Ka_Pz", &tTypeD_Ka_Pz);
    fOutputTree->Branch("ASD_DCAKaSV", &tTypeD_DCAKaSV);
    fOutputTree->Branch("ASD_DCAKaV0", &tTypeD_DCAKaV0);
    if (fIsMC) {
        fOutputTree->Branch("ASD_IsSignal", &tTypeD_IsSignal);
        fOutputTree->Branch("ASD_ReactionID", &tTypeD_ReactionID);
        fOutputTree->Branch("ASD_IsHybrid", &tTypeD_IsHybrid);
    }
}

/*            */
/**  Events  **/
/*** ====== ***/

Bool_t AliAnalysisTaskSexaquark::ProcessEvent() {
    /* Load MC event */
    if (fIsMC) {
        fMC = MCEvent();
        if (fMC == nullptr) AliFatal("AliMCEvent couldn't be found.");
        fPIDResponse->SetCurrentMCEvent(fMC);
    }
    /* Load reconstructed event */
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESD == nullptr) AliFatal("AliESDEvent couldn't be found.");
    /* Assign some branches (1) */
    fRunNumber = fESD->GetRunNumber();
    // `fDirNumber` should be set in `UserNotify()`
    // `fDirNumberB` should be set in `UserNotify()`
    fEventNumber = fESD->GetEventNumberInFile();
    fPrimaryVertex = fESD->GetPrimaryVertex();
    /* Apply selection */
    if (!PassesEventSelection()) return kFALSE;
    /* Assign rest of branches (2) */
    /*  */ auto* MultSelection = dynamic_cast<AliMultSelection*>(fESD->FindListObject("MultSelection"));
    /*  */ if (MultSelection == nullptr)
        /*  */ AliFatal("ERROR: AliMultSelection couldn't be found.");
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    fMagneticField = (Float_t)fESD->GetMagneticField();
    if (fIsMC) {
        fMC_PrimaryVertex = fMC->GetPrimaryVertex();
        tEvents_MC_PV_Xv = (Float_t)fMC_PrimaryVertex->GetX();
        tEvents_MC_PV_Yv = (Float_t)fMC_PrimaryVertex->GetY();
        tEvents_MC_PV_Zv = (Float_t)fMC_PrimaryVertex->GetZ();
    }
    tEvents_PV_Xv = (Float_t)fPrimaryVertex->GetX();
    tEvents_PV_Yv = (Float_t)fPrimaryVertex->GetY();
    tEvents_PV_Zv = (Float_t)fPrimaryVertex->GetZ();
    /*  */ Double_t PV_CovMatrix[6];
    /*  */ fPrimaryVertex->GetCovarianceMatrix(PV_CovMatrix);
    for (size_t i = 0; i < SexaConst::PV_CovMatrix_Size; i++) tEvent_PV_CovMatrix.at(i) = (Float_t)PV_CovMatrix[i];
    tEvents_NTracks = fESD->GetNumberOfTracks();
    tEvents_NTPCClusters = fESD->GetNumberOfTPCClusters();

    return kTRUE;
}

Bool_t AliAnalysisTaskSexaquark::PassesEventSelection() {
    //
    if (!fEventCuts.AcceptEvent(fESD)) return kFALSE;
    /* Reference: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/AliDPGRunList18r1 */
    if (!fIsMC && (fRunNumber == 296749 || fRunNumber == 296750 || fRunNumber == 296849 || fRunNumber == 296890 || fRunNumber == 297029 ||
                   fRunNumber == 297194 || fRunNumber == 297219 || fRunNumber == 297481)) {
        fEventCuts.UseTimeRangeCut();
        fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7);
        if (!fEventCuts.AcceptEvent(fESD)) return kFALSE;
    }
    /* Pileup Events */
    if (!fEventCuts.PassedCut(AliEventCuts::kPileUp)) return kFALSE;
    /* TPC Pileup Events */
    if (!fEventCuts.PassedCut(AliEventCuts::kTPCPileUp)) return kFALSE;
    /* Important for data? */
    if (!fIsMC && fESD->GetHeader()->GetEventType() != 7) return kFALSE;
    /* Trigger Selection */
    Bool_t IsMB = (fInputHandler->IsEventSelected() & AliVEvent::kINT7) != 0U;
    Bool_t IsHighMultV0 = (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0) != 0U;
    Bool_t IsHighMultSPD = (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) != 0U;
    Bool_t IsCentral = (fInputHandler->IsEventSelected() & AliVEvent::kCentral) != 0U;
    Bool_t IsSemiCentral = (fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral) != 0U;
    if (!IsMB && !IsHighMultV0 && !IsHighMultSPD && !IsCentral && !IsSemiCentral) return kFALSE;
    /* PV z-component */
    if (TMath::Abs(fPrimaryVertex->GetZ()) > SexaCuts::Event::AbsMax_PV_Zv) return kFALSE;

    return kTRUE;
}

/*                  */
/**  MC Generated  **/
/*** ============ ***/

void AliAnalysisTaskSexaquark::ProcessMCParticles() {
    //
    AliMCParticle* mcPart = nullptr;
    /* Auxiliary variables */
    Int_t pdg_code;
    Int_t mother_mc_idx;
    Bool_t is_signal;
    UInt_t reaction_id;
    for (Int_t mc_idx = 0; mc_idx < fMC->GetNumberOfTracks(); mc_idx++) {
        mcPart = dynamic_cast<AliMCParticle*>(fMC->GetTrack(mc_idx));
        if (mcPart == nullptr) continue;
        if (mcPart->P() < 0.01) continue;
        /* Get properties */
        pdg_code = mcPart->PdgCode();
        mother_mc_idx = mcPart->GetMother();
        is_signal = mcPart->GetGeneratorIndex() == 2;
        reaction_id = 0;
        if (mother_mc_idx == -1) {
            if (is_signal) {
                reaction_id = mcPart->MCStatusCode();                     // if no mother and signal, extract reaction id from status code
                fReactionProducts_McIdx_[reaction_id].push_back(mc_idx);  // store the reaction products for `ProcessInjected()`
            }
        } else {
            reaction_id = fMC_ReactionID_[mother_mc_idx];
        }
        /* Store */
        fMC_Mother_McIdx_[mc_idx] = mother_mc_idx;
        fMC_PdgCode_[mc_idx] = pdg_code;
        fMC_IsSignal_[mc_idx] = is_signal;
        fMC_ReactionID_[mc_idx] = reaction_id;
    }  // end of loop over MC particles
}

/*              */
/**  Injected  **/
/*** ======== ***/

/*
 * Assign the in-memory values to the tree branches.
 */
void AliAnalysisTaskSexaquark::ProcessInjected() {
    //
    AliMCParticle* mcPart = nullptr;
    /* Auxiliary variables */
    UInt_t reaction_id;
    Int_t mc_idx;
    Double_t sum_px, sum_py, sum_pz;
    for (size_t i = 0; i < fReactionID_[fEventNumber].size(); i++) {
        reaction_id = fReactionID_[fEventNumber][i];
        //
        tInjected_ReactionID.push_back(reaction_id);
        tInjected_Px.push_back(fSexaquark_Px_[fEventNumber][i]);
        tInjected_Py.push_back(fSexaquark_Py_[fEventNumber][i]);
        tInjected_Pz.push_back(fSexaquark_Pz_[fEventNumber][i]);
        tInjected_Nucleon_Px.push_back(fNucleon_Px_[fEventNumber][i]);
        tInjected_Nucleon_Py.push_back(fNucleon_Py_[fEventNumber][i]);
        tInjected_Nucleon_Pz.push_back(fNucleon_Pz_[fEventNumber][i]);
        /* Loop over reaction products */
        sum_px = 0.;
        sum_py = 0.;
        sum_pz = 0.;
        for (size_t j = 0; j < fReactionProducts_McIdx_[reaction_id].size(); j++) {
            mc_idx = fReactionProducts_McIdx_[reaction_id][j];
            mcPart = dynamic_cast<AliMCParticle*>(fMC->GetTrack(mc_idx));
            if (j == 0) {
                tInjected_Xv.push_back((Float_t)mcPart->Xv());
                tInjected_Yv.push_back((Float_t)mcPart->Yv());
                tInjected_Zv.push_back((Float_t)mcPart->Zv());
            }
            sum_px += mcPart->Px();
            sum_py += mcPart->Py();
            sum_pz += mcPart->Pz();
        }
        tInjected_Post_Px.push_back((Float_t)sum_px);
        tInjected_Post_Py.push_back((Float_t)sum_py);
        tInjected_Post_Pz.push_back((Float_t)sum_pz);
    }
}

/*
 * Open the respective `sim.log` that corresponds to the `RunNumber+DirNumber` that's being analyzed.
 * From it, read the injected anti-sexaquark and struck nucleon kinematics and store them into a tree.
 */
void AliAnalysisTaskSexaquark::BringSignalLogs() {
    //
    TGrid* alien = nullptr;
    if (gGrid == nullptr) {
        alien = TGrid::Connect("alien://");
        if (alien == nullptr) return;
    }

    TString AliEn_Dir = fAliEnPath(0, fAliEnPath.Last('/'));

    TString orig_path = Form("%s/sim.log", AliEn_Dir.Data());
    AliInfoF("Copying file %s ...", orig_path.Data());

    /* assuming path ends with format `.../LHC23l1a3/A1.73/297595/001/sim.log` */
    auto AliEn_DirNumber = static_cast<Int_t>(fDirNumber);
    TObjArray* tokens = fAliEnPath.Tokenize("/");
    Int_t AliEn_RunNumber = (dynamic_cast<TObjString*>(tokens->At(tokens->GetEntries() - 3)))->GetString().Atoi();
    TString AliEn_SimSubSet = (dynamic_cast<TObjString*>(tokens->At(tokens->GetEntries() - 4)))->GetString();

    fSignalLog_NewBasename = Form("sim_%s_%i_%03i.log", AliEn_SimSubSet.Data(), AliEn_RunNumber, AliEn_DirNumber);

    if (AliEn_Dir.BeginsWith("alien://")) {
        gSystem->Exec(Form("alien.py cp %s file://./%s", orig_path.Data(), fSignalLog_NewBasename.Data()));
    } else {
        gSystem->Exec(Form("cp %s ./%s", orig_path.Data(), fSignalLog_NewBasename.Data()));
    }

    TString new_path = Form("%s/%s", gSystem->pwd(), fSignalLog_NewBasename.Data());
    AliInfoF("Signal log file ready at %s ...", new_path.Data());
}

/*
 * Load injected anti-sexaquark and struck nucleon info.
 * From the `sim.log` file that corresponds to an entire dir number into memory.
 */
Bool_t AliAnalysisTaskSexaquark::LoadSignalLogs() {
    TString new_path = Form("%s/%s", gSystem->pwd(), fSignalLog_NewBasename.Data());
    AliInfoF("Opening file %s ...", new_path.Data());
    std::ifstream SignalLog(new_path);
    if (!SignalLog.is_open()) {
        AliWarningF("Unable to open file %s", new_path.Data());
        return kFALSE;
    }
    /* Read file */
    Int_t CurrentEventNumber = -1;
    std::string cstr_line;
    TString tstr_line, csv;
    TObjArray* csv_arr = nullptr;
    while (std::getline(SignalLog, cstr_line)) {
        tstr_line = cstr_line;
        /* A new event has appeared */
        if (tstr_line.Contains("I-AliGenCocktail::Generate: Generator 1: AliGenHijing")) CurrentEventNumber++;
        if (!tstr_line.Contains("I-AliGenSexaquarkReaction::GenerateN: 6")) continue;
        csv = static_cast<TString>(tstr_line(38, tstr_line.Length() - 1));
        csv_arr = csv.Tokenize(",");
        /* Load content to memory */
        fReactionID_[CurrentEventNumber].push_back(dynamic_cast<TObjString*>(csv_arr->At(0))->String().Atoi());
        fSexaquark_Px_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(1))->String().Atof());
        fSexaquark_Py_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(2))->String().Atof());
        fSexaquark_Pz_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(3))->String().Atof());
        fNucleon_Px_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(4))->String().Atof());
        fNucleon_Py_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(5))->String().Atof());
        fNucleon_Pz_[CurrentEventNumber].push_back((Float_t) dynamic_cast<TObjString*>(csv_arr->At(6))->String().Atof());
    }  // end of loop over lines
    AliInfoF("Closing file %s ...", new_path.Data());
    SignalLog.close();
    return kTRUE;
}

/*
 * Clear the maps from memory.
 */
void AliAnalysisTaskSexaquark::ClearSignalLogs() {
    fReactionID_.clear();
    fSexaquark_Px_.clear();
    fSexaquark_Py_.clear();
    fSexaquark_Pz_.clear();
    fNucleon_Px_.clear();
    fNucleon_Py_.clear();
    fNucleon_Pz_.clear();
}

/*            */
/**  Tracks  **/
/*** ====== ***/

void AliAnalysisTaskSexaquark::ProcessTracks() {
    //
    AliESDtrack* track = nullptr;
    /* Loop over tracks */
    for (Int_t esd_idx = 0; esd_idx < fESD->GetNumberOfTracks(); esd_idx++) {
        track = fESD->GetTrack(esd_idx);
        /* Track selection */
        if (!PassesTrackSelection(track)) continue;
        /* PID */
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fAntiProton_Indices.push_back(esd_idx);
            if (track->Charge() > 0) fProton_Indices.push_back(esd_idx);
        }
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fNegKaon_Indices.push_back(esd_idx);
            if (track->Charge() > 0) fPosKaon_Indices.push_back(esd_idx);
        }
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fPiMinus_Indices.push_back(esd_idx);
            if (track->Charge() > 0) fPiPlus_Indices.push_back(esd_idx);
        }
        /* Fill container */
        fLinked_McIdx_[esd_idx] = TMath::Abs(track->GetLabel());
    }  // end of loop over tracks
}

Bool_t AliAnalysisTaskSexaquark::PassesTrackSelection(const AliESDtrack* track) {
    //
    const AliExternalTrackParam* track_param = track->GetInnerParam();
    if (track_param == nullptr) return kFALSE;

    if (track_param->Pt() < SexaCuts::Track::Min_Pt || track_param->Pt() > SexaCuts::Track::Max_Pt) return kFALSE;
    /*  */ Float_t n_sigma_proton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    /*  */ Float_t n_sigma_kaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
    /*  */ Float_t n_sigma_pion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    if (TMath::Abs(n_sigma_proton) > SexaCuts::Track::AbsMax_PID_NSigma && TMath::Abs(n_sigma_kaon) > SexaCuts::Track::AbsMax_PID_NSigma &&
        TMath::Abs(n_sigma_pion) > SexaCuts::Track::AbsMax_PID_NSigma) {
        return kFALSE;
    }
    if (TMath::Abs(track_param->Eta()) > SexaCuts::Track::AbsMax_Eta) return kFALSE;
    /*  */ UShort_t NTPCClusters = track->GetTPCNcls();
    if (NTPCClusters < SexaCuts::Track::Min_NTPCClusters) return kFALSE;
    /*  */ Double_t Chi2PerNTPCClusters = NTPCClusters > 0 ? track->GetTPCchi2() / (Double_t)NTPCClusters : 999.;
    if (Chi2PerNTPCClusters > SexaCuts::Track::Max_Chi2PerNTPCClusters) return kFALSE;
    /*  */ Bool_t tpc_status = ((track->GetStatus() & AliESDtrack::kTPCout) != 0U) && ((track->GetStatus() & AliESDtrack::kTPCrefit) != 0U);
    /*  */ Bool_t its_status = ((track->GetStatus() & AliESDtrack::kITSin) == 0U) && ((track->GetStatus() & AliESDtrack::kITSout) == 0U) &&
                               ((track->GetStatus() & AliESDtrack::kITSrefit) == 0U);
    if (SexaCuts::Track::TurnedOn_StatusCuts && !tpc_status) return kFALSE;
    if (SexaCuts::Track::TurnedOn_StatusCuts && !its_status) return kFALSE;
    if (SexaCuts::Track::TurnedOn_RejectKinks && track->GetKinkIndex(0) > 0) return kFALSE;
    /*  */ Float_t DCAxy_wrtPV, DCAz_wrtPV;
    /*  */ track->GetImpactParameters(DCAxy_wrtPV, DCAz_wrtPV);
    if (TMath::Abs(DCAxy_wrtPV) < SexaCuts::Track::AbsMin_DCAxy_wrtPV) return kFALSE;

    return kTRUE;
}

/*         */
/**  V0s  **/
/*** === ***/

void AliAnalysisTaskSexaquark::KF_FindV0s(Short_t pdg_code_v0, Short_t pdg_code_neg, Short_t pdg_code_pos) {
    //
    const auto mass_neg = TDatabasePDG::Instance()->GetParticle(pdg_code_neg)->Mass();
    const auto mass_pos = TDatabasePDG::Instance()->GetParticle(pdg_code_pos)->Mass();
    //
    std::unique_ptr<AliExternalTrackParam> neg_param = std::make_unique<AliExternalTrackParam>();
    std::unique_ptr<AliExternalTrackParam> pos_param = std::make_unique<AliExternalTrackParam>();
    std::unique_ptr<AliESDVertex> esd_v0;
    Double_t arr_v0[3], arr_v0_err[3];
    /* Choose tracks species to loop over */
    const std::vector<Int_t>& neg_indices = (pdg_code_neg == PdgCode::AntiProton) ? fAntiProton_Indices : fPiMinus_Indices;
    const std::vector<Int_t>& pos_indices = (pdg_code_pos == PdgCode::Proton) ? fProton_Indices : fPiPlus_Indices;
    /* Loop over all possible pairs of tracks */
    for (auto esd_idx_neg : neg_indices) {
        neg_param->Reset();
        neg_param->CopyFromVTrack(fESD->GetTrack(esd_idx_neg)->GetInnerParam());
        for (auto esd_idx_pos : pos_indices) {
            pos_param->Reset();
            pos_param->CopyFromVTrack(fESD->GetTrack(esd_idx_pos)->GetInnerParam());
            /* sanity check */
            if (esd_idx_neg == esd_idx_pos) continue;
            /* fit */
            KF_V0 this_v0;
            this_v0.kf_neg = CreateKFParticle(neg_param.get(), mass_neg, (Int_t)neg_param->Charge());
            this_v0.kf_pos = CreateKFParticle(pos_param.get(), mass_pos, (Int_t)pos_param->Charge());
            /*  */
            this_v0.kf.AddDaughter(this_v0.kf_neg);
            this_v0.kf.AddDaughter(this_v0.kf_pos);
            this_v0.kf.TransportToDecayVertex();
            /*  */
            this_v0.v3.SetCoordinates(this_v0.kf.GetX(), this_v0.kf.GetY(), this_v0.kf.GetZ());
            arr_v0[0] = this_v0.kf.GetX();
            arr_v0[1] = this_v0.kf.GetY();
            arr_v0[2] = this_v0.kf.GetZ();
            arr_v0_err[0] = this_v0.kf.GetErrX();
            arr_v0_err[1] = this_v0.kf.GetErrY();
            arr_v0_err[2] = this_v0.kf.GetErrZ();
            /*  */
            esd_v0 = std::make_unique<AliESDVertex>(arr_v0, arr_v0_err);
            neg_param->PropagateToDCA(esd_v0.get(), fMagneticField, 10., this_v0.impar_neg);
            this_v0.param_d_neg_v0 = TMath::Abs(neg_param->GetD(arr_v0[0], arr_v0[1], fMagneticField));
            pos_param->PropagateToDCA(esd_v0.get(), fMagneticField, 10., this_v0.impar_pos);
            this_v0.param_d_pos_v0 = TMath::Abs(pos_param->GetD(arr_v0[0], arr_v0[1], fMagneticField));
            /*  */
            this_v0.dca_neg_v0 = TMath::Abs(this_v0.kf_neg.GetDistanceFromVertex(this_v0.kf));
            this_v0.dcaxy_neg_v0 = TMath::Abs(this_v0.kf_neg.GetDistanceFromVertexXY(this_v0.kf));
            this_v0.dca_pos_v0 = TMath::Abs(this_v0.kf_pos.GetDistanceFromVertex(this_v0.kf));
            this_v0.dcaxy_pos_v0 = TMath::Abs(this_v0.kf_pos.GetDistanceFromVertexXY(this_v0.kf));
            /* kinematics */
            this_v0.lv_neg.SetCoordinates(neg_param->Px(), neg_param->Py(), neg_param->Pz(), mass_neg);
            this_v0.lv_pos.SetCoordinates(pos_param->Px(), pos_param->Py(), pos_param->Pz(), mass_pos);
            this_v0.lv = this_v0.lv_neg + this_v0.lv_pos;
            /*  */
            this_v0.cpa_wrt_pv = CosinePointingAngle(this_v0.lv.Vect(), this_v0.v3, v3_pv);
            this_v0.dca_wrt_pv = LinePointDCA(this_v0.lv.Vect(), this_v0.v3, v3_pv);
            this_v0.arm_qt = ArmenterosQt(this_v0.lv.Vect(), this_v0.lv_neg.Vect());
            this_v0.arm_alpha = ArmenterosAlpha(this_v0.lv.Vect(), this_v0.lv_neg.Vect(), this_v0.lv_pos.Vect());
            this_v0.dca_btw_dau = TMath::Abs(this_v0.kf_neg.GetDistanceFromParticle(this_v0.kf_pos));
            /* apply cuts and store V0 */
            if (!PassesV0CutsAs(this_v0, pdg_code_v0)) continue;
            this_v0.idx = (pdg_code_v0 == PdgCode::AntiLambda) ? (Int_t)kfAntiLambdas.size() : (Int_t)kfKaonsZeroShort.size();
            this_v0.idx_neg = esd_idx_neg;
            this_v0.idx_pos = esd_idx_pos;
            StoreV0As(this_v0, pdg_code_v0, pdg_code_neg, pdg_code_pos);
        }  // end of loop over pos
    }      // end of loop over neg
}

Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(const KF_V0& this_v0, Short_t pdg_code_v0) {
    //
    if (pdg_code_v0 == PdgCode::AntiLambda) {
        /* anti-lambdas */
        /* -- kinematics cuts */
        if (this_v0.lv.Pt() < SexaCuts::Lambda::Min_Pt) return kFALSE;
        if (this_v0.lv.M() < SexaCuts::Lambda::Min_Mass || this_v0.lv.M() > SexaCuts::Lambda::Max_Mass) return kFALSE;
        if (TMath::Abs(this_v0.lv.Eta()) > SexaCuts::Lambda::AbsMax_Eta) return kFALSE;
        if (this_v0.cpa_wrt_pv < SexaCuts::Lambda::Min_CPAwrtPV || this_v0.cpa_wrt_pv > SexaCuts::Lambda::Max_CPAwrtPV) return kFALSE;
        if (this_v0.dca_wrt_pv < SexaCuts::Lambda::Min_DCAwrtPV) return kFALSE;
        if (this_v0.arm_qt / TMath::Abs(this_v0.arm_alpha) > SexaCuts::Lambda::AbsMax_ArmQtOverAlpha) return kFALSE;
        /* -- geometric cuts */
        if (TMath::Abs(this_v0.v3.Z()) > SexaCuts::Lambda::AbsMax_Zv) return kFALSE;
        if (this_v0.v3.Rho() < SexaCuts::Lambda::Min_Radius || this_v0.v3.Rho() > SexaCuts::Lambda::Max_Radius) return kFALSE;
        if (this_v0.dca_neg_v0 > SexaCuts::Lambda::Max_DCAnegV0) return kFALSE;
        if (this_v0.dca_pos_v0 > SexaCuts::Lambda::Max_DCAposV0) return kFALSE;
        if (this_v0.dca_btw_dau > SexaCuts::Lambda::Max_DCAbtwDau) return kFALSE;
    } else if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        /* kaons zero short */
        /* -- kinematics cuts */
        if (this_v0.lv.Pt() < SexaCuts::KaonZeroShort::Min_Pt) return kFALSE;
        if (this_v0.lv.M() < SexaCuts::KaonZeroShort::Min_Mass || this_v0.lv.M() > SexaCuts::KaonZeroShort::Max_Mass) return kFALSE;
        if (TMath::Abs(this_v0.lv.Eta()) > SexaCuts::KaonZeroShort::AbsMax_Eta) return kFALSE;
        if (this_v0.cpa_wrt_pv < SexaCuts::KaonZeroShort::Min_CPAwrtPV || this_v0.cpa_wrt_pv > SexaCuts::KaonZeroShort::Max_CPAwrtPV) return kFALSE;
        if (this_v0.dca_wrt_pv < SexaCuts::KaonZeroShort::Min_DCAwrtPV) return kFALSE;
        /* -- geometric cuts */
        if (TMath::Abs(this_v0.v3.Z()) > SexaCuts::Lambda::AbsMax_Zv) return kFALSE;
        if (this_v0.v3.Rho() < SexaCuts::KaonZeroShort::Min_Radius || this_v0.v3.Rho() > SexaCuts::KaonZeroShort::Max_Radius) return kFALSE;
        if (this_v0.dca_neg_v0 > SexaCuts::KaonZeroShort::Max_DCAnegV0) return kFALSE;
        if (this_v0.dca_pos_v0 > SexaCuts::KaonZeroShort::Max_DCAposV0) return kFALSE;
        if (this_v0.dca_btw_dau > SexaCuts::KaonZeroShort::Max_DCAbtwDau) return kFALSE;
    } else {
        AliWarning("Unknown V0 type.");
        return kFALSE;
    }

    return kTRUE;
}

void AliAnalysisTaskSexaquark::StoreV0As(const KF_V0& this_v0, Short_t pdg_code_v0, Short_t pdg_code_neg, Short_t pdg_code_pos) {
    //
    tV0_Idx[pdg_code_v0].push_back(this_v0.idx);
    tV0_Px[pdg_code_v0].push_back((Float_t)this_v0.lv.Px());
    tV0_Py[pdg_code_v0].push_back((Float_t)this_v0.lv.Py());
    tV0_Pz[pdg_code_v0].push_back((Float_t)this_v0.lv.Pz());
    tV0_E[pdg_code_v0].push_back((Float_t)this_v0.lv.E());
    tV0_Xv[pdg_code_v0].push_back(this_v0.kf.GetX());
    tV0_Yv[pdg_code_v0].push_back(this_v0.kf.GetY());
    tV0_Zv[pdg_code_v0].push_back(this_v0.kf.GetZ());
    tV0_CPAwrtPV[pdg_code_v0].push_back((Float_t)this_v0.cpa_wrt_pv);
    tV0_DCAwrtPV[pdg_code_v0].push_back((Float_t)this_v0.dca_wrt_pv);
    tV0_ArmQt[pdg_code_v0].push_back((Float_t)this_v0.arm_qt);
    tV0_ArmAlpha[pdg_code_v0].push_back((Float_t)this_v0.arm_alpha);
    tV0_DCA_Daughters[pdg_code_v0].push_back(this_v0.dca_btw_dau);
    /*  */
    tV0_Neg_EsdIdx[pdg_code_v0].push_back(this_v0.idx_neg);
    tV0_Neg_Px[pdg_code_v0].push_back((Float_t)this_v0.lv_neg.Px());
    tV0_Neg_Py[pdg_code_v0].push_back((Float_t)this_v0.lv_neg.Py());
    tV0_Neg_Pz[pdg_code_v0].push_back((Float_t)this_v0.lv_neg.Pz());
    tV0_Neg_IPxy_V0[pdg_code_v0].push_back((Float_t)this_v0.impar_neg[0]);
    tV0_Neg_IPz_V0[pdg_code_v0].push_back((Float_t)this_v0.impar_neg[1]);
    tV0_Neg_TrackParamD_V0[pdg_code_v0].push_back((Float_t)this_v0.param_d_neg_v0);
    tV0_Neg_DCA_V0[pdg_code_v0].push_back(this_v0.dca_neg_v0);
    tV0_Neg_DCAxy_V0[pdg_code_v0].push_back(this_v0.dcaxy_neg_v0);
    /*  */
    tV0_Pos_EsdIdx[pdg_code_v0].push_back(this_v0.idx_pos);
    tV0_Pos_Px[pdg_code_v0].push_back((Float_t)this_v0.lv_pos.Px());
    tV0_Pos_Py[pdg_code_v0].push_back((Float_t)this_v0.lv_pos.Py());
    tV0_Pos_Pz[pdg_code_v0].push_back((Float_t)this_v0.lv_pos.Pz());
    tV0_Pos_IPxy_V0[pdg_code_v0].push_back((Float_t)this_v0.impar_pos[0]);
    tV0_Pos_IPz_V0[pdg_code_v0].push_back((Float_t)this_v0.impar_pos[1]);
    tV0_Pos_TrackParamD_V0[pdg_code_v0].push_back((Float_t)this_v0.param_d_pos_v0);
    tV0_Pos_DCA_V0[pdg_code_v0].push_back(this_v0.dca_pos_v0);
    tV0_Pos_DCAxy_V0[pdg_code_v0].push_back(this_v0.dcaxy_pos_v0);
    /* Fill container */
    if (pdg_code_v0 == PdgCode::AntiLambda) {
        kfAntiLambdas.push_back(this_v0);
    } else if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        kfKaonsZeroShort.push_back(this_v0);
    }
    /* True information */
    if (fIsMC) {
        Int_t mc_idx_v0 = -1;
        Int_t mc_pdg_code_v0 = 0;
        Bool_t is_signal = false;
        UInt_t reaction_id = 0;
        /*  */
        Int_t mc_idx_neg = fLinked_McIdx_[this_v0.idx_neg];
        Int_t mc_idx_pos = fLinked_McIdx_[this_v0.idx_pos];
        Bool_t has_mc = fMC_Mother_McIdx_[mc_idx_neg] != -1 && fMC_Mother_McIdx_[mc_idx_neg] == fMC_Mother_McIdx_[mc_idx_pos];
        if (has_mc) {
            mc_idx_v0 = fMC_Mother_McIdx_[mc_idx_neg];
            mc_pdg_code_v0 = fMC_PdgCode_[mc_idx_v0];
            Bool_t is_true = fMC_PdgCode_[mc_idx_neg] == pdg_code_neg && fMC_PdgCode_[mc_idx_pos] == pdg_code_pos && mc_pdg_code_v0 == pdg_code_v0;
            if (is_true) {
                is_signal = fMC_IsSignal_[mc_idx_v0];
                reaction_id = fMC_ReactionID_[mc_idx_v0];
            }
        }
        Bool_t is_hybrid = !is_signal &&  //
                           ((fMC_IsSignal_[mc_idx_neg] && !fMC_IsSignal_[mc_idx_pos]) || (!fMC_IsSignal_[mc_idx_neg] && fMC_IsSignal_[mc_idx_pos]));
        /* Update vectors */
        tV0_McIdx[pdg_code_v0].push_back(mc_idx_v0);
        tV0_PdgCode[pdg_code_v0].push_back(mc_pdg_code_v0);
        tV0_IsSignal[pdg_code_v0].push_back(is_signal);
        tV0_ReactionID[pdg_code_v0].push_back(reaction_id);
        tV0_IsHybrid[pdg_code_v0].push_back(is_hybrid);
        /* -- Fill true container */
        if (pdg_code_v0 == PdgCode::AntiLambda) {
            mcAntiLambdas.push_back(MC_V0{is_signal, reaction_id, is_hybrid});
        } else if (pdg_code_v0 == PdgCode::KaonZeroShort) {
            mcKaonsZeroShort.push_back(MC_V0{is_signal, reaction_id, is_hybrid});
        }
    }
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/** Channel A **/

void AliAnalysisTaskSexaquark::KF_FindSexaquarks_TypeA(Short_t pdg_struck_nucleon, const std::vector<Short_t>& pdg_reaction_products) {
    //
    const Double_t neutron_mass = TDatabasePDG::Instance()->GetParticle(pdg_struck_nucleon)->Mass();
    if (!kfAntiLambdas.size() || !kfKaonsZeroShort.size()) return;
    /*  */
    KF_V0 cp_v0a, cp_v0b;
    PxPyPzEVector lv_v0a, lv_v0b;
    PxPyPzEVector lv_sexa, lv_sexa_asdecay;
    /* loop over all pairs */
    for (const auto& v0a : kfAntiLambdas) {
        for (const auto& v0b : kfKaonsZeroShort) {
            /* sanity check */
            std::set<Int_t> unique_track_entries = {v0a.idx_neg, v0a.idx_pos, v0b.idx_neg, v0b.idx_pos};
            if (unique_track_entries.size() < 4) continue;
            cp_v0a = v0a;
            cp_v0b = v0b;
            /* fit */
            KFParticle kf_sexa(cp_v0a.kf, cp_v0b.kf);
            kf_sexa.SetProductionVertex(kf_pv);
            //
            kf_sexa.TransportToDecayVertex();
            cp_v0a.kf.SetProductionVertex(kf_sexa);
            cp_v0b.kf.SetProductionVertex(kf_sexa);
            //
            cp_v0a.kf.TransportToDecayVertex();
            cp_v0a.kf_neg.SetProductionVertex(cp_v0a.kf);
            cp_v0a.kf_pos.SetProductionVertex(cp_v0a.kf);
            //
            cp_v0b.kf.TransportToDecayVertex();
            cp_v0b.kf_neg.SetProductionVertex(cp_v0b.kf);
            cp_v0b.kf_pos.SetProductionVertex(cp_v0b.kf);
            /* transport tracks to V0s vertices */
            cp_v0a.kf_neg.TransportToProductionVertex();
            cp_v0a.kf_pos.TransportToProductionVertex();
            cp_v0b.kf_neg.TransportToProductionVertex();
            cp_v0b.kf_pos.TransportToProductionVertex();
            /*
            lv_v0a_neg.SetCoordinates(cp_v0a.kf_neg.Px(), cp_v0a.kf_neg.Py(), cp_v0a.kf_neg.Pz(), lv_v0a_neg.M());
            lv_v0a_pos.SetCoordinates(cp_v0a.kf_pos.Px(), cp_v0a.kf_pos.Py(), cp_v0a.kf_pos.Pz(), lv_v0a_pos.M());
            lv_v0b_neg.SetCoordinates(cp_v0b.kf_neg.Px(), cp_v0b.kf_neg.Py(), cp_v0b.kf_neg.Pz(), lv_v0b_neg.M());
            lv_v0b_pos.SetCoordinates(cp_v0b.kf_pos.Px(), cp_v0b.kf_pos.Py(), cp_v0b.kf_pos.Pz(), lv_v0b_pos.M());
            */
            /* transport V0s to secondary vertex */
            cp_v0a.kf.TransportToProductionVertex();
            cp_v0b.kf.TransportToProductionVertex();
            lv_v0a.SetCoordinates(cp_v0a.kf.Px(), cp_v0a.kf.Py(), cp_v0a.kf.Pz(), cp_v0a.kf.E());
            lv_v0b.SetCoordinates(cp_v0b.kf.Px(), cp_v0b.kf.Py(), cp_v0b.kf.Pz(), cp_v0b.kf.E());
            /* fill struct */
            lv_sexa.SetCoordinates(lv_v0a.Px() + lv_v0b.Px(), lv_v0a.Py() + lv_v0b.Py(), lv_v0a.Pz() + lv_v0b.Pz(),
                                   lv_v0a.E() + lv_v0b.E() - neutron_mass);
            lv_sexa_asdecay = lv_v0a + lv_v0b;
            /* apply cuts and store */
            if (PassesSexaquarkCuts_TypeA(kf_sexa, lv_sexa, lv_sexa_asdecay, cp_v0a.kf, cp_v0a.kf_neg, cp_v0a.kf_pos, cp_v0b.kf, cp_v0b.kf_neg,
                                          cp_v0b.kf_pos)) {
                StoreSexaquark_TypeA(cp_v0a.idx, cp_v0b.idx, kf_sexa, lv_sexa, lv_sexa_asdecay, cp_v0a.kf, lv_v0a, cp_v0a.kf_neg, cp_v0a.kf_pos,
                                     cp_v0b.kf, lv_v0b, cp_v0b.kf_neg, cp_v0b.kf_pos);
            }
        }  // end of loop over k0s
    }      // end of loop over (anti)lambdas
}

Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_TypeA(const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa,
                                                           const PxPyPzEVector& lv_sexa_asdecay,                                                  //
                                                           const KFParticle& kf_v0a, const KFParticle& kf_v0a_neg, const KFParticle& kf_v0a_pos,  //
                                                           const KFParticle& kf_v0b, const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos) {
    XYZPoint v3_sexa(kf_sexa.GetX(), kf_sexa.GetY(), kf_sexa.GetZ());
    /* -- kinematics-dependent cuts */
    if (TMath::Abs(lv_sexa.Rapidity()) > SexaCuts::ChannelA::AbsMax_Rapidity) return kFALSE;
    /*  */ Double_t cpa_wrt_pv = CosinePointingAngle(lv_sexa.Vect(), v3_sexa, v3_pv);
    if (cpa_wrt_pv < SexaCuts::ChannelA::Min_CPAwrtPV || cpa_wrt_pv > SexaCuts::ChannelA::Max_CPAwrtPV) return kFALSE;
    if (lv_sexa_asdecay.M() < SexaCuts::ChannelA::Min_MassAsDecay || lv_sexa_asdecay.M() > SexaCuts::ChannelA::Max_MassAsDecay) return kFALSE;
    /* -- geometry-exclusive cuts */
    if (v3_sexa.Rho() < SexaCuts::ChannelA::Min_Radius || v3_sexa.Rho() > SexaCuts::ChannelA::Max_Radius) return kFALSE;
    /*  */ Double_t dca_la_sv = TMath::Abs((Double_t)kf_v0a.GetDistanceFromVertex(kf_sexa));
    if (dca_la_sv > SexaCuts::ChannelA::Max_DCALaSV) return kFALSE;
    /*  */ Double_t dca_la_neg_sv = TMath::Abs((Double_t)kf_v0a_neg.GetDistanceFromVertex(kf_sexa));
    if (dca_la_neg_sv > SexaCuts::ChannelA::Max_DCALaNegSV) return kFALSE;
    /*  */ Double_t dca_la_pos_sv = TMath::Abs((Double_t)kf_v0a_pos.GetDistanceFromVertex(kf_sexa));
    if (dca_la_pos_sv > SexaCuts::ChannelA::Max_DCALaPosSV) return kFALSE;
    /*  */ Double_t dca_k0s_sv = TMath::Abs((Double_t)kf_v0b.GetDistanceFromVertex(kf_sexa));
    if (dca_k0s_sv > SexaCuts::ChannelA::Max_DCAK0SV) return kFALSE;
    /*  */ Double_t dca_k0s_neg_sv = TMath::Abs((Double_t)kf_v0b_neg.GetDistanceFromVertex(kf_sexa));
    if (dca_k0s_neg_sv > SexaCuts::ChannelA::Max_DCAK0NegSV) return kFALSE;
    /*  */ Double_t dca_k0s_pos_sv = TMath::Abs((Double_t)kf_v0b_pos.GetDistanceFromVertex(kf_sexa));
    if (dca_k0s_pos_sv > SexaCuts::ChannelA::Max_DCAK0PosSV) return kFALSE;
    /*  */ Double_t dca_btw_v0s = TMath::Abs((Double_t)kf_v0a.GetDistanceFromParticle(kf_v0b));
    if (dca_btw_v0s > SexaCuts::ChannelA::Max_DCAbtwV0s) return kFALSE;
    /*  */ Double_t decay_length_la = TMath::Abs((Double_t)kf_v0a.GetDecayLength());
    if (decay_length_la > SexaCuts::ChannelA::Max_DecayLengthLa) return kFALSE;
    /*  */ Double_t decay_length_k0s = TMath::Abs((Double_t)kf_v0b.GetDecayLength());
    if (decay_length_k0s > SexaCuts::ChannelA::Max_DecayLengthK0) return kFALSE;
    return kTRUE;
}

void AliAnalysisTaskSexaquark::StoreSexaquark_TypeA(Int_t idx_v0a, Int_t idx_v0b, const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa,
                                                    const PxPyPzEVector& lv_sexa_asdecay, const KFParticle& kf_v0a, const PxPyPzEVector& lv_v0a,
                                                    const KFParticle& kf_v0a_neg, const KFParticle& kf_v0a_pos, const KFParticle& kf_v0b,
                                                    const PxPyPzEVector& lv_v0b, const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos) {
    tTypeA_Px.push_back((Float_t)lv_sexa.Px());
    tTypeA_Py.push_back((Float_t)lv_sexa.Py());
    tTypeA_Pz.push_back((Float_t)lv_sexa.Pz());
    tTypeA_E.push_back((Float_t)lv_sexa.E());
    tTypeA_E_asDecay.push_back((Float_t)lv_sexa_asdecay.E());
    tTypeA_Xv.push_back(kf_sexa.GetX());
    tTypeA_Yv.push_back(kf_sexa.GetY());
    tTypeA_Zv.push_back(kf_sexa.GetZ());
    tTypeA_V0a_Idx.push_back(idx_v0a);
    tTypeA_V0a_Px.push_back((Float_t)lv_v0a.Px());
    tTypeA_V0a_Py.push_back((Float_t)lv_v0a.Py());
    tTypeA_V0a_Pz.push_back((Float_t)lv_v0a.Pz());
    tTypeA_V0a_E.push_back((Float_t)lv_v0a.E());
    tTypeA_V0a_DecayLength.push_back(TMath::Abs(kf_v0a.GetDecayLength()));
    tTypeA_DCAV0aSV.push_back(TMath::Abs(kf_v0a.GetDistanceFromVertex(kf_sexa)));
    tTypeA_DCAV0aNegSV.push_back(TMath::Abs(kf_v0a_neg.GetDistanceFromVertex(kf_sexa)));
    tTypeA_DCAV0aPosSV.push_back(TMath::Abs(kf_v0a_pos.GetDistanceFromVertex(kf_sexa)));
    tTypeA_V0b_Idx.push_back(idx_v0b);
    tTypeA_V0b_Px.push_back((Float_t)lv_v0b.Px());
    tTypeA_V0b_Py.push_back((Float_t)lv_v0b.Py());
    tTypeA_V0b_Pz.push_back((Float_t)lv_v0b.Pz());
    tTypeA_V0b_E.push_back((Float_t)lv_v0b.E());
    tTypeA_V0b_DecayLength.push_back(TMath::Abs(kf_v0b.GetDecayLength()));
    tTypeA_DCAV0bSV.push_back(TMath::Abs(kf_v0b.GetDistanceFromVertex(kf_sexa)));
    tTypeA_DCAV0bNegSV.push_back(TMath::Abs(kf_v0b_neg.GetDistanceFromVertex(kf_sexa)));
    tTypeA_DCAV0bPosSV.push_back(TMath::Abs(kf_v0b_pos.GetDistanceFromVertex(kf_sexa)));
    tTypeA_DCAbtwV0s.push_back(TMath::Abs(kf_v0a.GetDistanceFromParticle(kf_v0b)));
    if (fIsMC) {
        /*  */
        Bool_t is_signal = kFALSE;
        UInt_t reaction_id = 0;
        /* fill values */
        const MC_V0& mc_v0a = mcAntiLambdas[idx_v0a];
        const MC_V0& mc_v0b = mcKaonsZeroShort[idx_v0b];
        if (mc_v0a.reaction_id == mc_v0b.reaction_id) {
            reaction_id = mc_v0a.reaction_id;
            is_signal = mc_v0a.is_signal && mc_v0b.is_signal;
        }
        Bool_t is_hybrid = !is_signal && ((mc_v0a.is_signal && !mc_v0b.is_signal) || (!mc_v0a.is_signal && mc_v0b.is_signal) ||  //
                                          mc_v0a.is_hybrid || mc_v0b.is_hybrid);
        /*  */
        tTypeA_IsSignal.push_back(is_signal);
        tTypeA_ReactionID.push_back(reaction_id);
        tTypeA_IsHybrid.push_back(is_hybrid);
    }
};

/** Channel D **/

void AliAnalysisTaskSexaquark::KF_FindSexaquarks_TypeD(Short_t pdg_struck_nucleon, const std::vector<Short_t>& pdg_reaction_products) {
    //
    // PENDING
    //
}

Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_TypeD(const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa, const KFParticle& kf_v0,
                                                           const KFParticle& kf_v0_neg, const KFParticle& kf_v0_pos, const KFParticle& kf_ka) {
    XYZPoint v3_sexa(kf_sexa.GetX(), kf_sexa.GetY(), kf_sexa.GetZ());

    /* -- kinematics-dependent cuts */
    if (TMath::Abs(lv_sexa.Rapidity()) > SexaCuts::ChannelD::AbsMax_Rapidity) return kFALSE;
    /*  */ Double_t cpa_wrt_pv = CosinePointingAngle(lv_sexa.Vect(), v3_sexa, v3_pv);
    if (cpa_wrt_pv < SexaCuts::ChannelD::Min_CPAwrtPV || cpa_wrt_pv > SexaCuts::ChannelD::Max_CPAwrtPV) return kFALSE;
    /* -- geometry-exclusive cuts */
    if (v3_sexa.Rho() < SexaCuts::ChannelD::Min_Radius || v3_sexa.Rho() > SexaCuts::ChannelD::Max_Radius) return kFALSE;
    /*  */ Double_t dca_la_sv = TMath::Abs((Double_t)kf_v0.GetDistanceFromVertex(kf_sexa));
    if (dca_la_sv > SexaCuts::ChannelD::Max_DCALaSV) return kFALSE;
    /*  */ Double_t dca_la_neg_sv = TMath::Abs((Double_t)kf_v0_neg.GetDistanceFromVertex(kf_sexa));
    if (dca_la_neg_sv > SexaCuts::ChannelD::Max_DCALaNegSV) return kFALSE;
    /*  */ Double_t dca_la_pos_sv = TMath::Abs((Double_t)kf_v0_pos.GetDistanceFromVertex(kf_sexa));
    if (dca_la_pos_sv > SexaCuts::ChannelD::Max_DCALaPosSV) return kFALSE;
    /*  */ Double_t dca_ka_sv = TMath::Abs((Double_t)kf_ka.GetDistanceFromVertex(kf_sexa));
    if (dca_ka_sv > SexaCuts::ChannelD::Max_DCAKaSV) return kFALSE;
    /*  */ Double_t dca_ka_la = TMath::Abs((Double_t)kf_ka.GetDistanceFromVertex(kf_v0));
    if (dca_ka_la > SexaCuts::ChannelD::Max_DCAKaLa) return kFALSE;
    /*  */ Double_t dca_la_neg_ka = TMath::Abs((Double_t)kf_v0_neg.GetDistanceFromParticle(kf_ka));
    if (dca_la_neg_ka > SexaCuts::ChannelD::Max_DCALaNegKa) return kFALSE;
    /*  */ Double_t dca_la_pos_ka = TMath::Abs((Double_t)kf_v0_pos.GetDistanceFromParticle(kf_ka));
    if (dca_la_pos_ka > SexaCuts::ChannelD::Max_DCALaPosKa) return kFALSE;

    return kTRUE;
}

void AliAnalysisTaskSexaquark::StoreSexaquark_TypeD(Int_t idx_v0, Int_t esd_idx_ka, const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa,
                                                    const KFParticle& kf_v0, const PxPyPzEVector& lv_v0, const KFParticle& kf_v0_neg,
                                                    const KFParticle& kf_v0_pos, const KFParticle& kf_ka, const PxPyPzEVector& lv_ka) {
    PxPyPzEVector lv_sexa_asdecay = lv_v0 + lv_ka;

    tTypeD_Px.push_back((Float_t)lv_sexa.Px());
    tTypeD_Py.push_back((Float_t)lv_sexa.Py());
    tTypeD_Pz.push_back((Float_t)lv_sexa.Pz());
    tTypeD_E.push_back((Float_t)lv_sexa.E());
    tTypeD_E_asDecay.push_back((Float_t)lv_sexa_asdecay.E());
    tTypeD_Xv.push_back(kf_sexa.GetX());
    tTypeD_Yv.push_back(kf_sexa.GetY());
    tTypeD_Zv.push_back(kf_sexa.GetZ());
    tTypeD_V0_Idx.push_back(idx_v0);
    tTypeD_V0_Px.push_back((Float_t)lv_v0.Px());
    tTypeD_V0_Py.push_back((Float_t)lv_v0.Py());
    tTypeD_V0_Pz.push_back((Float_t)lv_v0.Pz());
    tTypeD_V0_E.push_back((Float_t)lv_v0.E());
    tTypeD_V0_DecayLength.push_back(TMath::Abs(kf_v0.GetDecayLength()));
    tTypeD_DCAV0SV.push_back(TMath::Abs(kf_v0.GetDistanceFromVertex(kf_sexa)));
    tTypeD_DCAV0NegSV.push_back(TMath::Abs(kf_v0_neg.GetDistanceFromVertex(kf_sexa)));
    tTypeD_DCAV0PosSV.push_back(TMath::Abs(kf_v0_pos.GetDistanceFromVertex(kf_sexa)));
    tTypeD_DCAV0NegKa.push_back(TMath::Abs(kf_v0_neg.GetDistanceFromParticle(kf_ka)));
    tTypeD_DCAV0PosKa.push_back(TMath::Abs(kf_v0_pos.GetDistanceFromParticle(kf_ka)));
    tTypeD_Ka_EsdIdx.push_back(esd_idx_ka);
    tTypeD_Ka_Px.push_back((Float_t)lv_ka.Px());
    tTypeD_Ka_Py.push_back((Float_t)lv_ka.Py());
    tTypeD_Ka_Pz.push_back((Float_t)lv_ka.Pz());
    tTypeD_DCAKaSV.push_back(TMath::Abs(kf_ka.GetDistanceFromVertex(kf_sexa)));
    tTypeD_DCAKaV0.push_back(TMath::Abs(kf_ka.GetDistanceFromVertex(kf_v0)));
    if (fIsMC) {
        tTypeD_IsSignal.push_back(false);  // PENDING
        tTypeD_ReactionID.push_back(0);    // PENDING
        tTypeD_IsHybrid.push_back(false);  // PENDING
    }
}

/*                */
/**  Containers  **/
/*** ========== ***/

void AliAnalysisTaskSexaquark::ClearBranches_Injected() {
    /* `tInjected_Nucleon_PdgCode` not cleared on purpose, as an event branch */
    /* `tInjected_Mass` not cleared on purpose, as an event branch */
    tInjected_ReactionID.clear();
    tInjected_Px.clear();
    tInjected_Py.clear();
    tInjected_Pz.clear();
    tInjected_Nucleon_Px.clear();
    tInjected_Nucleon_Py.clear();
    tInjected_Nucleon_Pz.clear();
    tInjected_Xv.clear();
    tInjected_Yv.clear();
    tInjected_Zv.clear();
    tInjected_Post_Px.clear();
    tInjected_Post_Py.clear();
    tInjected_Post_Pz.clear();
}

void AliAnalysisTaskSexaquark::ClearBranches_V0s() {
    for (auto v0 : {PdgCode::AntiLambda, PdgCode::KaonZeroShort}) {
        tV0_Idx[v0].clear();
        tV0_Px[v0].clear();
        tV0_Py[v0].clear();
        tV0_Pz[v0].clear();
        tV0_E[v0].clear();
        tV0_Xv[v0].clear();
        tV0_Yv[v0].clear();
        tV0_Zv[v0].clear();
        tV0_CPAwrtPV[v0].clear();
        tV0_DCAwrtPV[v0].clear();
        tV0_ArmQt[v0].clear();
        tV0_ArmAlpha[v0].clear();
        tV0_DCA_Daughters[v0].clear();
        tV0_Neg_EsdIdx[v0].clear();
        tV0_Neg_Px[v0].clear();
        tV0_Neg_Py[v0].clear();
        tV0_Neg_Pz[v0].clear();
        tV0_Neg_IPxy_V0[v0].clear();
        tV0_Neg_IPz_V0[v0].clear();
        tV0_Neg_TrackParamD_V0[v0].clear();
        tV0_Neg_DCA_V0[v0].clear();
        tV0_Neg_DCAxy_V0[v0].clear();
        tV0_Pos_EsdIdx[v0].clear();
        tV0_Pos_Px[v0].clear();
        tV0_Pos_Py[v0].clear();
        tV0_Pos_Pz[v0].clear();
        tV0_Pos_IPxy_V0[v0].clear();
        tV0_Pos_IPz_V0[v0].clear();
        tV0_Pos_TrackParamD_V0[v0].clear();
        tV0_Pos_DCA_V0[v0].clear();
        tV0_Pos_DCAxy_V0[v0].clear();
        if (fIsMC) {
            tV0_McIdx[v0].clear();
            tV0_PdgCode[v0].clear();
            tV0_IsSignal[v0].clear();
            tV0_ReactionID[v0].clear();
            tV0_IsHybrid[v0].clear();
        }
    }
}

void AliAnalysisTaskSexaquark::ClearBranches_TypeA() {
    tTypeA_Px.clear();
    tTypeA_Py.clear();
    tTypeA_Pz.clear();
    tTypeA_E.clear();
    tTypeA_E_asDecay.clear();
    tTypeA_Xv.clear();
    tTypeA_Yv.clear();
    tTypeA_Zv.clear();
    tTypeA_V0a_Idx.clear();
    tTypeA_V0a_Px.clear();
    tTypeA_V0a_Py.clear();
    tTypeA_V0a_Pz.clear();
    tTypeA_V0a_E.clear();
    tTypeA_V0a_DecayLength.clear();
    tTypeA_DCAV0aSV.clear();
    tTypeA_DCAV0aNegSV.clear();
    tTypeA_DCAV0aPosSV.clear();
    tTypeA_V0b_Idx.clear();
    tTypeA_V0b_Px.clear();
    tTypeA_V0b_Py.clear();
    tTypeA_V0b_Pz.clear();
    tTypeA_V0b_E.clear();
    tTypeA_V0b_DecayLength.clear();
    tTypeA_DCAV0bSV.clear();
    tTypeA_DCAV0bNegSV.clear();
    tTypeA_DCAV0bPosSV.clear();
    tTypeA_DCAbtwV0s.clear();
    if (fIsMC) {
        tTypeA_IsSignal.clear();
        tTypeA_ReactionID.clear();
        tTypeA_IsHybrid.clear();
    }
}

void AliAnalysisTaskSexaquark::ClearBranches_TypeD() {
    tTypeD_Px.clear();
    tTypeD_Py.clear();
    tTypeD_Pz.clear();
    tTypeD_E.clear();
    tTypeD_E_asDecay.clear();
    tTypeD_Xv.clear();
    tTypeD_Yv.clear();
    tTypeD_Zv.clear();
    tTypeD_V0_Idx.clear();
    tTypeD_V0_Px.clear();
    tTypeD_V0_Py.clear();
    tTypeD_V0_Pz.clear();
    tTypeD_V0_E.clear();
    tTypeD_V0_DecayLength.clear();
    tTypeD_DCAV0SV.clear();
    tTypeD_DCAV0NegSV.clear();
    tTypeD_DCAV0PosSV.clear();
    tTypeD_DCAV0NegKa.clear();
    tTypeD_DCAV0PosKa.clear();
    tTypeD_Ka_EsdIdx.clear();
    tTypeD_Ka_Px.clear();
    tTypeD_Ka_Py.clear();
    tTypeD_Ka_Pz.clear();
    tTypeD_DCAKaSV.clear();
    tTypeD_DCAKaV0.clear();
    if (fIsMC) {
        tTypeD_IsSignal.clear();
        tTypeD_ReactionID.clear();
        tTypeD_IsHybrid.clear();
    }
}

void AliAnalysisTaskSexaquark::ClearContainers() {
    //
    fMC_PdgCode_.clear();
    fMC_Mother_McIdx_.clear();
    fMC_IsSignal_.clear();
    fMC_ReactionID_.clear();
    /* */
    fReactionProducts_McIdx_.clear();
    /* */
    fLinked_McIdx_.clear();
    /* */
    fAntiProton_Indices.clear();
    fProton_Indices.clear();
    fNegKaon_Indices.clear();
    fPosKaon_Indices.clear();
    fPiMinus_Indices.clear();
    fPiPlus_Indices.clear();
    /*  */
    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
    mcAntiLambdas.clear();
    mcKaonsZeroShort.clear();
}
