struct MCEvent_tt {
    Int_t Idx;        //
    Int_t DirNumber;  //
    Int_t RunNumber;  //
    Float_t Xv_PV;    //
    Float_t Yv_PV;    //
    Float_t Zv_PV;    //
};

struct MC_tt {
    Int_t PdgCode;         // pdg code
    Int_t Idx;             //
    Int_t Idx_Mother;      // index of mother (within MC_tt)
    Int_t Idx_NDaughters;  // number of daughters
    Int_t Idx_FirstDau;    // index of first daughter (within MC_tt)
    Int_t Idx_LastDau;     // index of last daughter (within MC_tt)
    Float_t Px;            // momentum
    Float_t Py;            //
    Float_t Pz;            //
    Float_t Xv_i;          // origin vertex
    Float_t Yv_i;          //
    Float_t Zv_i;          //
    Int_t Status;          // MC status code
    Bool_t IsSignal;       //
    Bool_t IsSecondary;    //
};

struct RecEvent_tt {
    Int_t Idx;           //
    Int_t DirNumber;     //
    Int_t RunNumber;     //
    Float_t Centrality;  //
    Float_t Xv_PV;       //
    Float_t Yv_PV;       //
    Float_t Zv_PV;       //
};

struct Track_tt {
    Int_t Idx;             //
    Int_t Idx_True;        // index of true particle (within MC_tt)
    Float_t Px;            // momentum
    Float_t Py;            //
    Float_t Pz;            //
    Short_t Charge;        //
    Float_t NSigmaPion;    //
    Float_t NSigmaKaon;    //
    Float_t NSigmaProton;  //
};

struct V0_tt {
    Int_t Idx;        //
    Int_t Idx_Pos;    // index of positive track (within Track_tt)
    Int_t Idx_Neg;    // index of negative track (within Track_tt)
    Int_t Idx_True;   // index of true particle (within MC_tt)
    Float_t Px;       // momentum
    Float_t Py;       //
    Float_t Pz;       //
    Float_t E;        // energy
    Float_t Xv_f;     // decay vertex
    Float_t Yv_f;     //
    Float_t Zv_f;     //
    Float_t Px_Pos;   // momentum of positive daughter at v0 position
    Float_t Py_Pos;   //
    Float_t Pz_Pos;   //
    Float_t Px_Neg;   //
    Float_t Py_Neg;   //
    Float_t Pz_Neg;   //
    Bool_t IsSignal;  // ktrue if it belongs to anti-sexaquark signal, kfalse if background
};

struct Sxqk_A {
    Int_t Idx;              //
    Int_t DirNumber;        //
    Int_t RunNumber;        //
    Int_t ReactionIdx;      //
    Int_t Idx_AntiLambda;   // index of anti-lambda (within V0_tt)
    Int_t Idx_K0S;          // index of kaon-zero-short (within V0_tt)
    Float_t Px;             // momentum
    Float_t Py;             //
    Float_t Pz;             //
    Float_t E;              // energy
    Float_t Xv_SV;          // secondary vertex
    Float_t Yv_SV;          //
    Float_t Zv_SV;          //
    Float_t Chi2ndf;        //
    Float_t Weight_Pt;      //
    Float_t Weight_Radius;  //
};
