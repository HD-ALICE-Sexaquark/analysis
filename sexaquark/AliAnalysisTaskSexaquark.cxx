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
      tV0_Neg_EsdIdx(),
      tV0_Neg_Px(),
      tV0_Neg_Py(),
      tV0_Neg_Pz(),
      tV0_Pos_EsdIdx(),
      tV0_Pos_Px(),
      tV0_Pos_Py(),
      tV0_Pos_Pz(),
      tV0_DCAnegV0(),
      tV0_DCAposV0(),
      tV0_DCAbtwDau(),
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
      fAntiProton_EsdIdx(),
      fProton_EsdIdx(),
      fNegKaon_EsdIdx(),
      fPosKaon_EsdIdx(),
      fPiMinus_EsdIdx(),
      fPiPlus_EsdIdx(),
      /*  */
      kfAntiLambdas(),
      kfKaonsZeroShort() {}

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
      tV0_Neg_EsdIdx(),
      tV0_Neg_Px(),
      tV0_Neg_Py(),
      tV0_Neg_Pz(),
      tV0_Pos_EsdIdx(),
      tV0_Pos_Px(),
      tV0_Pos_Py(),
      tV0_Pos_Pz(),
      tV0_DCAnegV0(),
      tV0_DCAposV0(),
      tV0_DCAbtwDau(),
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
      fAntiProton_EsdIdx(),
      fProton_EsdIdx(),
      fNegKaon_EsdIdx(),
      fPosKaon_EsdIdx(),
      fPiMinus_EsdIdx(),
      fPiPlus_EsdIdx(),
      /*  */
      kfAntiLambdas(),
      kfKaonsZeroShort() {
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
    /* Prepare output tree */
    fOutputTree = new TTree("Events", "Events");
    AssociateBranches_Events();
    if (fIsMC && fIsSignalMC) AssociateBranches_Injected();
    // AssociateBranches_V0s();
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
void AliAnalysisTaskSexaquark::UserExec(Option_t*) {
    /* Events */
    if (!ProcessEvent()) return;
    /* Set KFParticle */
    KFParticle::SetField(fMagneticField);
    kf_pv = CreateKFVertex(fPrimaryVertex);
    /* MC Particles */
    if (fIsMC) {
        ProcessMCParticles();
        if (fIsSignalMC) ProcessInjected();
    }
    /* Tracks */
    ProcessTracks();
    /* V0s */
    // KF_FindV0s(PdgCode::AntiLambda, PdgCode::AntiProton, PdgCode::PiPlus);
    // KF_FindV0s(PdgCode::KaonZeroShort, PdgCode::PiMinus, PdgCode::PiPlus);
    /* Sexaquarks */
    /* -- `AntiSexaquark,Neutron -> AntiLambda,K0S` */
    // KF_FindSexaquarks_TypeA(PdgCode::Neutron, {PdgCode::AntiLambda, PdgCode::KaonZeroShort});
    /* -- `AntiSexaquark,Proton -> AntiLambda,K+,(pi-,pi+)` */
    // KF_FindSexaquarks_TypeD(PdgCode::AntiNeutron, {PdgCode::Lambda, PdgCode::KaonZeroShort});
    /* Fill tree */
    fOutputTree->Fill();
    /* End of event */
    if (fIsMC && fIsSignalMC) ClearBranches_Injected();
    // ClearBranches_V0s();
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
        fOutputTree->Branch(v0_names[i] + "_Neg_EsdIdx", &tV0_Neg_EsdIdx[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Px", &tV0_Neg_Px[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Py", &tV0_Neg_Py[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Neg_Pz", &tV0_Neg_Pz[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_EsdIdx", &tV0_Pos_EsdIdx[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Px", &tV0_Pos_Px[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Py", &tV0_Pos_Py[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_Pos_Pz", &tV0_Pos_Pz[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_DCAnegV0", &tV0_DCAnegV0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_DCAposV0", &tV0_DCAposV0[v0_pdg_codes[i]]);
        fOutputTree->Branch(v0_names[i] + "_DCAbtwDau", &tV0_DCAbtwDau[v0_pdg_codes[i]]);
        if (fIsMC) {
            fOutputTree->Branch(v0_names[i] + "_Idx_True", &tV0_McIdx[v0_pdg_codes[i]]);
            fOutputTree->Branch(v0_names[i] + "_True_PdgCode", &tV0_PdgCode[v0_pdg_codes[i]]);
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
    const AliExternalTrackParam* trackInnerParam = nullptr;
    /* Loop over tracks */
    for (Int_t esd_idx = 0; esd_idx < fESD->GetNumberOfTracks(); esd_idx++) {
        track = fESD->GetTrack(esd_idx);
        trackInnerParam = track->GetInnerParam();
        if (trackInnerParam == nullptr) continue;
        /* Track selection */
        if (!PassesTrackSelection(track)) continue;
        /* PID */
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fAntiProton_EsdIdx.push_back(esd_idx);
            if (track->Charge() > 0) fProton_EsdIdx.push_back(esd_idx);
        }
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fNegKaon_EsdIdx.push_back(esd_idx);
            if (track->Charge() > 0) fPosKaon_EsdIdx.push_back(esd_idx);
        }
        if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) < SexaCuts::Track::AbsMax_PID_NSigma) {
            if (track->Charge() < 0) fPiMinus_EsdIdx.push_back(esd_idx);
            if (track->Charge() > 0) fPiPlus_EsdIdx.push_back(esd_idx);
        }
    }  // end of loop over tracks
}

Bool_t AliAnalysisTaskSexaquark::PassesTrackSelection(const AliESDtrack* track) {
    //
    const AliExternalTrackParam* track_param = track->GetInnerParam();
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
    AliESDtrack *neg_track, *pos_track;
    /* Declare KFParticle objects */
    KFParticle kf_v0, kf_neg, kf_pos;
    /* Declare 4-momentum vectors */
    PxPyPzMVector lv_v0, lv_neg, lv_pos;
    /* Choose tracks species to loop over */
    const std::vector<Int_t>& neg_indices = pdg_code_v0 == PdgCode::AntiLambda ? fAntiProton_EsdIdx : fPiMinus_EsdIdx;
    const std::vector<Int_t>& pos_indices = fPiPlus_EsdIdx;
    /* Loop over all possible pairs of tracks */
    for (auto esd_idx_neg : neg_indices) {
        for (auto esd_idx_pos : pos_indices) {
            /* Sanity check */
            if (esd_idx_neg == esd_idx_pos) continue;
            /* Get tracks */
            neg_track = fESD->GetTrack(esd_idx_neg);
            pos_track = fESD->GetTrack(esd_idx_pos);
            /* Kalman Filter */
            kf_neg = CreateKFParticle(neg_track, mass_neg, (Int_t)neg_track->Charge());
            kf_pos = CreateKFParticle(pos_track, mass_pos, (Int_t)pos_track->Charge());
            kf_v0.AddDaughter(kf_neg);
            kf_v0.AddDaughter(kf_pos);
            /* Transport V0 and daughters */
            kf_v0.TransportToDecayVertex();
            /* Reconstruct V0 */
            lv_neg.SetCoordinates(kf_neg.Px(), kf_neg.Py(), kf_neg.Pz(), mass_neg);
            lv_pos.SetCoordinates(kf_pos.Px(), kf_pos.Py(), kf_pos.Pz(), mass_pos);
            lv_v0 = lv_neg + lv_pos;
            /* Apply cuts and store V0 */
            if (!PassesV0CutsAs(pdg_code_v0, kf_v0, kf_neg, kf_pos, lv_v0, lv_neg, lv_pos)) continue;
            StoreV0As(pdg_code_v0, esd_idx_neg, esd_idx_pos, kf_v0, kf_neg, kf_pos, lv_v0, lv_neg, lv_pos);
        }  // end of loop over pos
    }      // end of loop over neg
}

Bool_t AliAnalysisTaskSexaquark::PassesV0CutsAs(Short_t pdg_code_v0, const KFParticle& kf_v0, const KFParticle& kf_neg, const KFParticle& kf_pos,
                                                const PxPyPzMVector& lv_v0, const PxPyPzMVector& lv_neg, const PxPyPzMVector& lv_pos) {
    //
    XYZPoint v3_v0(kf_v0.GetX(), kf_v0.GetY(), kf_v0.GetZ());

    if (pdg_code_v0 == PdgCode::AntiLambda) {
        /* anti-lambdas */
        /* -- kinematics cuts */
        if (lv_v0.Pt() < SexaCuts::Lambda::Min_Pt) return kFALSE;
        if (lv_v0.M() < SexaCuts::Lambda::Min_Mass || lv_v0.M() > SexaCuts::Lambda::Max_Mass) return kFALSE;
        if (TMath::Abs(lv_v0.Eta()) > SexaCuts::Lambda::AbsMax_Eta) return kFALSE;
        /*  */ Double_t CPAwrtPV = CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
        if (CPAwrtPV < SexaCuts::Lambda::Min_CPAwrtPV || CPAwrtPV > SexaCuts::Lambda::Max_CPAwrtPV) return kFALSE;
        /*  */ Double_t DCAwrtPV = LinePointDCA(lv_v0.Vect(), v3_v0, v3_pv);
        if (DCAwrtPV < SexaCuts::Lambda::Min_DCAwrtPV || DCAwrtPV > SexaCuts::Lambda::Max_DCAwrtPV) return kFALSE;
        /*  */ Double_t ArmQt = ArmenterosQt(lv_v0.Vect(), lv_neg.Vect());
        /*  */ Double_t ArmAlpha = ArmenterosAlpha(lv_v0.Vect(), lv_neg.Vect(), lv_pos.Vect());
        /*  */ Double_t ArmQtOverAlpha = ArmQt / TMath::Abs(ArmAlpha);
        if (ArmQtOverAlpha > SexaCuts::Lambda::AbsMax_ArmQtOverAlpha) return kFALSE;
        /* -- geometry cuts */
        if (v3_v0.Rho() < SexaCuts::Lambda::Min_Radius || v3_v0.Rho() > SexaCuts::Lambda::Max_Radius) return kFALSE;
        /*  */ Double_t DCAnegV0 = TMath::Abs((Double_t)kf_neg.GetDistanceFromVertex(kf_v0));
        if (DCAnegV0 > SexaCuts::Lambda::Max_DCAnegV0) return kFALSE;
        /*  */ Double_t DCAposV0 = TMath::Abs((Double_t)kf_pos.GetDistanceFromVertex(kf_v0));
        if (DCAposV0 > SexaCuts::Lambda::Max_DCAposV0) return kFALSE;
        /*  */ Double_t DCAbtwDau = TMath::Abs((Double_t)kf_neg.GetDistanceFromParticle(kf_pos));
        if (DCAbtwDau > SexaCuts::Lambda::Max_DCAbtwDau) return kFALSE;
    } else if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        /* kaons zero short */
        /* -- kinematics cuts */
        if (lv_v0.Pt() < SexaCuts::KaonZeroShort::Min_Pt) return kFALSE;
        if (lv_v0.M() < SexaCuts::KaonZeroShort::Min_Mass || lv_v0.M() > SexaCuts::KaonZeroShort::Max_Mass) return kFALSE;
        if (TMath::Abs(lv_v0.Eta()) > SexaCuts::KaonZeroShort::AbsMax_Eta) return kFALSE;
        /*  */ Double_t CPAwrtPV = CosinePointingAngle(lv_v0.Vect(), v3_v0, v3_pv);
        if (CPAwrtPV < SexaCuts::KaonZeroShort::Min_CPAwrtPV || CPAwrtPV > SexaCuts::KaonZeroShort::Max_CPAwrtPV) return kFALSE;
        /*  */ Double_t DCAwrtPV = LinePointDCA(lv_v0.Vect(), v3_v0, v3_pv);
        if (DCAwrtPV < SexaCuts::KaonZeroShort::Min_DCAwrtPV || DCAwrtPV > SexaCuts::KaonZeroShort::Max_DCAwrtPV) return kFALSE;
        /* -- geometry cuts */
        if (v3_v0.Rho() < SexaCuts::KaonZeroShort::Min_Radius || v3_v0.Rho() > SexaCuts::KaonZeroShort::Max_Radius) return kFALSE;
        /*  */ Double_t DCAnegV0 = TMath::Abs((Double_t)kf_neg.GetDistanceFromVertex(kf_v0));
        if (DCAnegV0 > SexaCuts::KaonZeroShort::Max_DCAnegV0) return kFALSE;
        /*  */ Double_t DCAposV0 = TMath::Abs((Double_t)kf_pos.GetDistanceFromVertex(kf_v0));
        if (DCAposV0 > SexaCuts::KaonZeroShort::Max_DCAposV0) return kFALSE;
        /*  */ Double_t DCAbtwDau = TMath::Abs((Double_t)kf_neg.GetDistanceFromParticle(kf_pos));
        if (DCAbtwDau > SexaCuts::KaonZeroShort::Max_DCAbtwDau) return kFALSE;
    } else {
        AliWarning("Unknown V0 type.");
        return kFALSE;
    }

    return kTRUE;
}

void AliAnalysisTaskSexaquark::StoreV0As(Short_t pdg_code_v0, Int_t esd_idx_neg, Int_t esd_idx_pos, const KFParticle& kf_v0, const KFParticle& kf_neg,
                                         const KFParticle& kf_pos, const PxPyPzMVector& lv_v0, const PxPyPzMVector& lv_neg,
                                         const PxPyPzMVector& lv_pos) {
    /* Determine container */
    std::vector<KFParticle>* kfFoundV0s;
    if (pdg_code_v0 == PdgCode::KaonZeroShort) {
        kfFoundV0s = &kfKaonsZeroShort;
    } else if (pdg_code_v0 == PdgCode::AntiLambda) {
        kfFoundV0s = &kfAntiLambdas;
    }
    /* Assign branches and fill tree */
    tV0_Idx[pdg_code_v0].push_back((Int_t)kfFoundV0s->size());
    tV0_Px[pdg_code_v0].push_back((Float_t)lv_v0.Px());
    tV0_Py[pdg_code_v0].push_back((Float_t)lv_v0.Py());
    tV0_Pz[pdg_code_v0].push_back((Float_t)lv_v0.Pz());
    tV0_E[pdg_code_v0].push_back((Float_t)lv_v0.E());
    tV0_Xv[pdg_code_v0].push_back(kf_v0.GetX());
    tV0_Yv[pdg_code_v0].push_back(kf_v0.GetY());
    tV0_Zv[pdg_code_v0].push_back(kf_v0.GetZ());
    tV0_Neg_EsdIdx[pdg_code_v0].push_back(esd_idx_neg);
    tV0_Neg_Px[pdg_code_v0].push_back((Float_t)lv_neg.Px());
    tV0_Neg_Py[pdg_code_v0].push_back((Float_t)lv_neg.Py());
    tV0_Neg_Pz[pdg_code_v0].push_back((Float_t)lv_neg.Pz());
    tV0_Pos_EsdIdx[pdg_code_v0].push_back(esd_idx_pos);
    tV0_Pos_Px[pdg_code_v0].push_back((Float_t)lv_pos.Px());
    tV0_Pos_Py[pdg_code_v0].push_back((Float_t)lv_pos.Py());
    tV0_Pos_Pz[pdg_code_v0].push_back((Float_t)lv_pos.Pz());
    tV0_DCAnegV0[pdg_code_v0].push_back(TMath::Abs(kf_neg.GetDistanceFromVertex(kf_v0)));      // OPTIMIZE
    tV0_DCAposV0[pdg_code_v0].push_back(TMath::Abs(kf_pos.GetDistanceFromVertex(kf_v0)));      // OPTIMIZE
    tV0_DCAbtwDau[pdg_code_v0].push_back(TMath::Abs(kf_neg.GetDistanceFromParticle(kf_pos)));  // OPTIMIZE
    if (fIsMC) {
        tV0_McIdx[pdg_code_v0].push_back(0);         // PENDING
        tV0_PdgCode[pdg_code_v0].push_back(0);       // PENDING
        tV0_IsSignal[pdg_code_v0].push_back(false);  // PENDING
        tV0_ReactionID[pdg_code_v0].push_back(0);    // PENDING
        tV0_IsHybrid[pdg_code_v0].push_back(false);  // PENDING
    }
    /* Fill container */
    kfFoundV0s->push_back(kf_v0);
}

/*                                */
/**  Sexaquark -- Kalman Filter  **/
/*** ========================== ***/

/** Channel A **/

void AliAnalysisTaskSexaquark::KF_FindSexaquarks_TypeA(Short_t pdg_struck_nucleon, const std::vector<Short_t>& pdg_reaction_products) {
    //
    // PENDING
    //
}

Bool_t AliAnalysisTaskSexaquark::PassesSexaquarkCuts_TypeA(const KFParticle& kf_sexa, const PxPyPzEVector& lv_sexa,
                                                           const PxPyPzEVector& lv_sexa_asdecay,                                                  //
                                                           const KFParticle& kf_v0a, const KFParticle& kf_v0a_neg, const KFParticle& kf_v0a_pos,  //
                                                           const KFParticle& kf_v0b, const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos) {
    //
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
                                                    const KFParticle& kf_v0a, const PxPyPzEVector& lv_v0a, const KFParticle& kf_v0a_neg,
                                                    const KFParticle& kf_v0a_pos, const KFParticle& kf_v0b, const PxPyPzEVector& lv_v0b,
                                                    const KFParticle& kf_v0b_neg, const KFParticle& kf_v0b_pos) {

    PxPyPzEVector lv_sexa_asdecay = lv_v0a + lv_v0b;

    /* Assign branches and fill tree */
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
        tTypeA_IsSignal.push_back(0);    // PENDING
        tTypeA_ReactionID.push_back(0);  // PENDING
        tTypeA_IsHybrid.push_back(0);    // PENDING
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

/*                             */
/**  Kalman Filter Utilities  **/
/*** ======================= ***/

/*
 * Correct initialization of a KFParticle.
 * (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFParticle AliAnalysisTaskSexaquark::CreateKFParticle(const AliExternalTrackParam* track_param, Double_t mass, Int_t charge) {

    Double_t fP[6];
    track_param->GetXYZ(fP);
    track_param->PxPyPz(fP + 3);

    Int_t fQ = track_param->Charge() * TMath::Abs(charge);
    fP[3] *= TMath::Abs(charge);
    fP[4] *= TMath::Abs(charge);
    fP[5] *= TMath::Abs(charge);

    Double_t pt = 1. / TMath::Abs(track_param->GetParameter()[4]) * TMath::Abs(charge);
    Double_t cs = TMath::Cos(track_param->GetAlpha());
    Double_t sn = TMath::Sin(track_param->GetAlpha());
    Double_t r = TMath::Sqrt((1. - track_param->GetParameter()[2]) * (1. + track_param->GetParameter()[2]));

    Double_t m00 = -sn;
    Double_t m10 = cs;
    Double_t m23 = -pt * (sn + track_param->GetParameter()[2] * cs / r);
    Double_t m43 = -pt * pt * (r * cs - track_param->GetParameter()[2] * sn);
    Double_t m24 = pt * (cs - track_param->GetParameter()[2] * sn / r);
    Double_t m44 = -pt * pt * (r * sn + track_param->GetParameter()[2] * cs);
    Double_t m35 = pt;
    Double_t m45 = -pt * pt * track_param->GetParameter()[3];

    m43 *= track_param->GetSign();
    m44 *= track_param->GetSign();
    m45 *= track_param->GetSign();

    const Double_t* cTr = track_param->GetCovariance();
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
 * Correct initialization of a KFVertex.
 * (Copied from `AliPhysics/PWGLF/.../AliAnalysisTaskDoubleHypNucTree.cxx`)
 */
KFVertex AliAnalysisTaskSexaquark::CreateKFVertex(const AliVVertex* vertex) {

    Double_t param[6];
    vertex->GetXYZ(param);

    Double_t cov[6];
    vertex->GetCovarianceMatrix(cov);

    KFPVertex kfpVtx;
    Float_t paramF[3] = {(Float_t)param[0], (Float_t)param[1], (Float_t)param[2]};
    kfpVtx.SetXYZ(paramF);
    Float_t covF[6] = {(Float_t)cov[0], (Float_t)cov[1], (Float_t)cov[2], (Float_t)cov[3], (Float_t)cov[4], (Float_t)cov[5]};
    kfpVtx.SetCovarianceMatrix(covF);

    KFVertex KFVtx(kfpVtx);

    return KFVtx;
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
    tV0_Idx.clear();
    tV0_Px.clear();
    tV0_Py.clear();
    tV0_Pz.clear();
    tV0_E.clear();
    tV0_Xv.clear();
    tV0_Yv.clear();
    tV0_Zv.clear();
    tV0_Neg_EsdIdx.clear();
    tV0_Neg_Px.clear();
    tV0_Neg_Py.clear();
    tV0_Neg_Pz.clear();
    tV0_Pos_EsdIdx.clear();
    tV0_Pos_Px.clear();
    tV0_Pos_Py.clear();
    tV0_Pos_Pz.clear();
    tV0_DCAnegV0.clear();
    tV0_DCAposV0.clear();
    tV0_DCAbtwDau.clear();
    if (fIsMC) {
        tV0_McIdx.clear();
        tV0_PdgCode.clear();
        tV0_IsSignal.clear();
        tV0_ReactionID.clear();
        tV0_IsHybrid.clear();
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
    fAntiProton_EsdIdx.clear();
    fProton_EsdIdx.clear();
    fNegKaon_EsdIdx.clear();
    fPosKaon_EsdIdx.clear();
    fPiMinus_EsdIdx.clear();
    fPiPlus_EsdIdx.clear();
    /*  */
    kfAntiLambdas.clear();
    kfKaonsZeroShort.clear();
}
