/**
 * Tree Name: Injected
 * Branch Name: Injected
 */
struct Injected_tt {
    Int_t RunNumber;
    Int_t DirNumber;
    Int_t EventNumber;
    Int_t ReactionID;
    Float_t Px;
    Float_t Py;
    Float_t Pz;
    Int_t PdgCode_Nucleon;
    Float_t Px_Nucleon;
    Float_t Py_Nucleon;
    Float_t Pz_Nucleon;
    Char_t ReactionChannel;
    Float_t M;
};

/**
 * Tree Name: Events
 * Branch Name: Event
 */
struct Event_tt {
    Int_t Number;        //
    Int_t DirNumber;     //
    Int_t RunNumber;     //
    Float_t Centrality;  //
    Float_t MC_Xv_PV;    //
    Float_t MC_Yv_PV;    //
    Float_t MC_Zv_PV;    //
    Float_t Rec_Xv_PV;   //
    Float_t Rec_Yv_PV;   //
    Float_t Rec_Zv_PV;   //
};

/**
 * Tree Name: MCParticles
 * Branch Name: MCParticle
 */
struct MC_tt {
    Int_t EventNumber;   //
    Int_t DirNumber;     //
    Int_t RunNumber;     //
    Int_t ReactionID;    //
    Int_t PdgCode;       // pdg code
    Int_t Idx;           //
    Int_t Idx_Mother;    // index of mother (within MC_tt)
    Int_t NDaughters;    // number of daughters
    Int_t Idx_FirstDau;  // index of first daughter (within MC_tt)
    Int_t Idx_LastDau;   // index of last daughter (within MC_tt)
    Float_t Px;          // momentum
    Float_t Py;          //
    Float_t Pz;          //
    Float_t Xv_i;        // origin vertex
    Float_t Yv_i;        //
    Float_t Zv_i;        //
    Float_t Xv_f;        // decay vertex
    Float_t Yv_f;        //
    Float_t Zv_f;        //
    Int_t Status;        // MC status code
    Bool_t IsSecondary;  //
    Bool_t IsSignal;     //
};

/**
 * Tree Name: Tracks
 * Branch Name: Track
 */
struct Track_tt {
    Int_t ReactionID;      //
    Int_t EventNumber;     //
    Int_t DirNumber;       //
    Int_t RunNumber;       //
    Int_t Idx;             //
    Int_t Idx_True;        // index of true particle (within MC_tt)
    Float_t Px;            // momentum
    Float_t Py;            //
    Float_t Pz;            //
    Short_t Charge;        //
    Float_t NSigmaPion;    //
    Float_t NSigmaKaon;    //
    Float_t NSigmaProton;  //
    Bool_t IsSecondary;    //
    Bool_t IsSignal;       //
};

/**
 * Tree Name: V0s
 * Branch Name: V0
 */
struct V0_tt {
    Int_t ReactionID;    //
    Int_t EventNumber;   //
    Int_t DirNumber;     //
    Int_t RunNumber;     //
    Int_t Idx;           //
    Int_t Idx_Pos;       // index of positive track (within Track_tt)
    Int_t Idx_Neg;       // index of negative track (within Track_tt)
    Int_t Idx_True;      // index of true particle (within MC_tt)
    Float_t Px;          // momentum
    Float_t Py;          //
    Float_t Pz;          //
    Float_t E;           // energy
    Float_t Xv_f;        // decay vertex
    Float_t Yv_f;        //
    Float_t Zv_f;        //
    Float_t Px_Pos;      // momentum of positive daughter at v0 position
    Float_t Py_Pos;      //
    Float_t Pz_Pos;      //
    Float_t Px_Neg;      //
    Float_t Py_Neg;      //
    Float_t Pz_Neg;      //
    Bool_t IsTrue;       //
    Bool_t IsSecondary;  //
    Bool_t IsSignal;     //
};

/**
 * Tree Name: Sexaquarks
 * Branch Name: Sexaquark
 */
struct RecSexaquark_aa {
    Int_t ReactionID;   //
    Int_t EventNumber;  //
    Int_t DirNumber;    //
    Int_t RunNumber;    //
    Int_t Idx_AL;       // index of anti-lambda (within AntiLambdas' V0_tt)
    Int_t Idx_K0S;      // index of kaon-zero-short (within KaonsZeroShort' V0_tt)
    Int_t Idx_AL_Neg;   // index of anti-lambda negative track (within Track_tt)
    Int_t Idx_AL_Pos;   // index of anti-lambda positive track (within Track_tt)
    Int_t Idx_K0S_Neg;  // index of kaon-zero-short negative track (within Track_tt)
    Int_t Idx_K0S_Pos;  // index of kaon-zero-short positive track (within Track_tt)
    Float_t Px;         // momentum
    Float_t Py;         //
    Float_t Pz;         //
    Float_t E;          // energy
    Float_t Xv_SV;      // secondary vertex
    Float_t Yv_SV;      //
    Float_t Zv_SV;      //
    Float_t Chi2ndf;    //
};
