#ifndef AliAnalysisTaskSexaquark_H
#define AliAnalysisTaskSexaquark_H

#include "AliAnalysisTaskSE.h"

struct V0_tt {
    //
    // Default common structure for the V0s
    // (since I want to study 3 alternatives, let's keep the same format for each)
    //
    Int_t Idx_Pos;          // index of positive track
    Int_t Idx_Neg;          // index of negative track
    Float_t Px;             // x-component of V0 momentum
    Float_t Py;             // y-component of V0 momentum
    Float_t Pz;             // z-component of V0 momentum
    Float_t X;              // x-coordinate of V0
    Float_t Y;              // y-coordinate of V0
    Float_t Z;              // z-coordinate of V0
    Float_t Pos_Px;         // x-component of positive track momentum at V0 position
    Float_t Pos_Py;         // y-component of positive track momentum at V0 position
    Float_t Pos_Pz;         // z-component of positive track momentum at V0 position
    Float_t Neg_Px;         // x-component of negative track momentum at V0 position
    Float_t Neg_Py;         // y-component of negative track momentum at V0 position
    Float_t Neg_Pz;         // z-component of negative track momentum at V0 position
    Bool_t isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    Float_t E_asK0;         // energy of V0, assuming it's a K0
    Float_t E_asAL;         // energy of V0, assuming it's an anti-lambda
    Bool_t couldBeK0;       // kTRUE if K0 candidate, kFALSE if not
    Bool_t couldBeAL;       // kTRUE if anti-lambda candidate, kFALSE if not
    Float_t Chi2;           // chi2 value from the Kalman Filter (?)
    Float_t DCA_Daughters;  // distance of closest approach between daughters
    Float_t IP_wrtPV;       // impact parameter w.r.t. Primary Vertex
    Float_t CPA_wrtPV;      // cosine of pointing angle w.r.t. Primary Vertex
    Float_t ArmAlpha;       // Armenteros-Podolanski variable alpha
    Float_t ArmPt;          // Armenteros-Podolanski variable Pt
    Float_t DecayLength;    // distance between PV and V0
    Float_t DCA_wrtPV;      // distance of closest approach to Primary Vertex
};

struct Event_tt {
    //
    // This is the structure of the output TTree
    //

    /* MC particles */

    Int_t N_MCGen;                          // number of MC particles
    std::vector<Float_t> MC_Px;             // x-component of true momentum
    std::vector<Float_t> MC_Py;             // y-component of true momentum
    std::vector<Float_t> MC_Pz;             // z-component of true momentum
    std::vector<Float_t> MC_X;              // initial x-coordinate of generation vertex
    std::vector<Float_t> MC_Y;              // initial y-coordinate of generation vertex
    std::vector<Float_t> MC_Z;              // initial z-coordinate of generation vertex
    std::vector<Float_t> MC_Xf;             // final x-coordinate of generation vertex
    std::vector<Float_t> MC_Yf;             // final y-coordinate of generation vertex
    std::vector<Float_t> MC_Zf;             // final z-coordinate of generation vertex
    std::vector<Int_t> MC_PID;              // PDG code
    std::vector<Int_t> MC_Mother;           // index of mother
    std::vector<Int_t> MC_PID_Mother;       // PID of mother
    std::vector<Int_t> MC_PID_GrandMother;  // PID of grandmother
    std::vector<Int_t> MC_NDaughters;       // number of daughters
    std::vector<Int_t> MC_FirstDau;         // index of first daughter
    std::vector<Int_t> MC_LastDau;          // index of last daughter
    std::vector<Int_t> MC_Status;           // MC status code
    std::vector<Bool_t> MC_isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background

    /* MC Rec. Tracks */

    Int_t N_MCRec;                          // number of MC reconstructed tracks
    std::vector<Int_t> Idx_True;            // index of true MC particle
    std::vector<Float_t> Rec_Px;            // x-component of reconstructed momentum
    std::vector<Float_t> Rec_Py;            // y-component of reconstructed momentum
    std::vector<Float_t> Rec_Pz;            // z-component of reconstructed momentum
    std::vector<Short_t> Rec_Charge;        // measured charge
    std::vector<Float_t> Rec_IP_wrtPV;      // transverse impact parameter w.r.t. Primary Vertex
    std::vector<Float_t> Rec_NSigmaPion;    // absolute value of likeness to be a charged pion (closer to 0, the most likely)
    std::vector<Float_t> Rec_NSigmaProton;  // absolute value of likeness to be a proton (closer to 0, the most likely)
    std::vector<Short_t> Rec_NClustersTPC;  // number of clusters assigned in the TPC
    std::vector<Bool_t> Rec_isDuplicate;    // kTRUE if track is a duplicate, kFALSE if not (v1: compared to true particle)
    std::vector<Bool_t> Rec_isSimilar;      // kTRUE if track is a duplicate, kFALSE if not (v2: compared to similar particles)
    std::vector<Bool_t> Rec_isSignal;       // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    std::vector<Float_t> Rec_HelixParam0;   // helix param. 0
    std::vector<Float_t> Rec_HelixParam1;   // helix param. 1
    std::vector<Float_t> Rec_HelixParam2;   // helix param. 2
    std::vector<Float_t> Rec_HelixParam3;   // helix param. 3
    std::vector<Float_t> Rec_HelixParam4;   // helix param. 4
    std::vector<Float_t> Rec_HelixParam5;   // helix param. 5

    /* V0s */

    Int_t N_V0s;                            // number of formed V0s
    std::vector<Int_t> Idx_Pos;             // index of positive daughter
    std::vector<Int_t> Idx_Neg;             // index of negative daughter
    std::vector<Float_t> V0_Px;             // x-component of V0 momentum
    std::vector<Float_t> V0_Py;             // y-component of V0 momentum
    std::vector<Float_t> V0_Pz;             // z-component of V0 momentum
    std::vector<Float_t> V0_X;              // x-coordinate of V0
    std::vector<Float_t> V0_Y;              // y-coordinate of V0
    std::vector<Float_t> V0_Z;              // z-coordinate of V0
    std::vector<Float_t> Pos_Px;            // x-component of positive track momentum at V0 position
    std::vector<Float_t> Pos_Py;            // y-component of positive track momentum at V0 position
    std::vector<Float_t> Pos_Pz;            // z-component of positive track momentum at V0 position
    std::vector<Float_t> Neg_Px;            // x-component of negative track momentum at V0 position
    std::vector<Float_t> Neg_Py;            // y-component of negative track momentum at V0 position
    std::vector<Float_t> Neg_Pz;            // z-component of negative track momentum at V0 position
    std::vector<Bool_t> V0_isSignal;        // kTRUE if it belongs to anti-sexaquark signal, kFALSE if background
    std::vector<Float_t> V0_E_asK0;         // energy of V0, assuming it's a K0
    std::vector<Float_t> V0_E_asAL;         // energy of V0, assuming it's an anti-lambda
    std::vector<Bool_t> V0_couldBeK0;       // kTRUE if K0 candidate, kFALSE if not
    std::vector<Bool_t> V0_couldBeAL;       // kTRUE if anti-lambda candidate, kFALSE if not
    std::vector<Float_t> V0_Chi2;           // chi2 value from the Kalman Filter (?)
    std::vector<Float_t> V0_DCA_Daughters;  // distance of closest approach between daughters
    std::vector<Float_t> V0_IP_wrtPV;       // impact parameter w.r.t. Primary Vertex
    std::vector<Float_t> V0_CPA_wrtPV;      // cosine of pointing angle w.r.t. Primary Vertex
    std::vector<Float_t> V0_ArmAlpha;       // Armenteros-Podolanski variable alpha
    std::vector<Float_t> V0_ArmPt;          // Armenteros-Podolanski variable Pt
    std::vector<Float_t> V0_DecayLength;    // distance between PV and V0
    std::vector<Float_t> V0_DCA_wrtPV;      // distance of closest approach to Primary Vertex
};

class AliPIDResponse;

class AliAnalysisTaskSexaquark : public AliAnalysisTaskSE {
   public:
    AliAnalysisTaskSexaquark();
    AliAnalysisTaskSexaquark(const char* name, Bool_t IsMC, TString SourceOfV0s, Char_t SimulationSet);
    virtual ~AliAnalysisTaskSexaquark();

   public:
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t* option);
    virtual void Terminate(Option_t* option) { return; }

   public:
    virtual void CheckForInputErrors();

    // MC Gen.
    virtual void MCGen_FindRelevant();
    virtual void MCGen_GetMotherInfo(AliMCParticle* mcPart, Int_t& evt_mother, Int_t& mother_pid, Int_t& grandmother_pid);
    virtual Bool_t MCGen_IsSignal(AliMCParticle* mcPart);
    virtual void MCGen_AddNonRelevant();
    virtual void MCGen_GetDaughtersInfo(AliMCParticle* mcPart, Int_t& n_dau, Int_t& n_first_dau, Int_t& n_last_dau);

    // MC Rec.
    virtual void MCRec_FindRelevant();
    virtual void MCRec_FindDuplicates();
    virtual void MCRec_FindSimilarTracks();
    virtual Float_t MCRec_GetImpactParameter(AliESDtrack* track);

    // V0 Finders
    virtual void V0_ProcessESD();
    virtual void V0_ProcessTrue();
    virtual void V0_ProcessCustom();

    // custom offline V0 finder (based on AliV0vertexer.cxx)
    virtual Bool_t Preoptimize(const AliExternalTrackParam* nt, AliExternalTrackParam* pt, Double_t* lPreprocessxn, Double_t* lPreprocessxp,
                               const Double_t b);
    virtual void GetHelixCenter(const AliExternalTrackParam* track, Double_t center[2], const Double_t b);

    // (mathematical functions)
    Double_t Calculate_CPA(Double_t Px, Double_t Py, Double_t Pz,  //
                           Double_t X, Double_t Y, Double_t Z,     //
                           Double_t refPointX, Double_t refPointY, Double_t refPointZ);
    Double_t Calculate_ArmAlpha(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,     //
                                Double_t Pos_Px, Double_t Pos_Py, Double_t Pos_Pz,  //
                                Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz);
    Double_t Calculate_ArmPt(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                             Double_t Neg_Px, Double_t Neg_Py, Double_t Neg_Pz);
    Double_t Calculate_LinePointDCA(Double_t V0_Px, Double_t V0_Py, Double_t V0_Pz,  //
                                    Double_t V0_X, Double_t V0_Y, Double_t V0_Z,     //
                                    Double_t PV_X, Double_t PV_Y, Double_t PV_Z);

   public:
    /* Tree operations */
    virtual void SetBranches();
    virtual void MCGen_PushBack(Int_t evt_mc);
    virtual void MCRec_PushBack(Int_t evt_track);
    virtual void V0_PushBack(Int_t evt_v0);
    virtual void ClearEvent();

   private:
    /* Input options */
    Bool_t fIsMC;
    // source of V0 reconstruction:
    // - "online"
    // - "offline"
    // - "custom"
    // - "true"
    // - "kalman"
    TString fSourceOfV0s;
    // reaction channel, could be:
    // - 'A' = AntiS + N -> AntiL + K0
    // - 'B' = AntiS + N -> AntiL + K0 + pim + pip
    // - 'C' = AntiS + N -> AntiP + K0 + K0 + pip
    // - 'D' = AntiS + P -> AntiL + Kp
    // - 'E' = AntiS + P -> AntiL + Kp + pim + pip
    // - 'F' = AntiS + P -> AntiP + Kp + K0 + pip
    // - 'G' = AntiS + N -> Xip + pim
    // - 'H' = AntiS + P -> AntiP + Kp + Kp + pi0
    Char_t fSimulationSet;

    Event_tt fEvent;

    std::unordered_map<Int_t, Int_t> fMap_Evt_MCGen;  // [Idx, Evt]
    std::vector<Int_t> fVec_Idx_MCGen;

    std::unordered_map<Int_t, Int_t> fMap_Evt_MCRec;  // [Idx, Evt]
    std::vector<Int_t> fVec_Idx_MCRec;

    std::vector<Bool_t> fVec_MCRec_IsDuplicate;
    std::vector<Bool_t> fVec_MCRec_IsSimilar;

    std::vector<Int_t> fVec_Idx_MCGen_NonRel;             // (non.rel.->rel.)
    std::map<Int_t, std::vector<Int_t>> fMap_Duplicates;  // (duplicated tracks) [Idx_MCGen, {Idx_MCRec_0, ...}]

    std::vector<V0_tt> fVec_V0s;

   private:
    /* AliRoot objects */
    AliMCEvent* fMC;               // corresponding MC event
    AliESDEvent* fESD;             // input event
    AliPIDResponse* fPIDResponse;  // pid response object
    AliESDVertex* fPrimaryVertex;  // primary vertex

   private:
    /* ROOT objects */
    TDatabasePDG fPDG;
    TList* fOutputListOfTrees;
    TTree* fTree;

    AliAnalysisTaskSexaquark(const AliAnalysisTaskSexaquark&);             // not implemented
    AliAnalysisTaskSexaquark& operator=(const AliAnalysisTaskSexaquark&);  // not implemented

    ClassDef(AliAnalysisTaskSexaquark, 1);
};

#endif
