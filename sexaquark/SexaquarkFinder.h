#include "macros/include/Headers.hxx"
#include "macros/include/Math.hxx"
#include "macros/include/TreeFunctions.hxx"
#include "macros/include/Utilities.hxx"

struct Sexa_tt {

    /* Identifiers */

    Int_t RunNumber;
    Int_t DirNumber;
    Int_t Event;

    /* Anti-Sexaquark Candidate */

    Int_t Idx_V0A;         // index of first daughter
    Int_t Idx_V0B;         // index of second daughter
    Float_t Sexa_E;        // energy of anti-sexaquark candidate
    Float_t Sexa_Px;       // x-component of anti-sexaquark candidate
    Float_t Sexa_Py;       // y-component of anti-sexaquark candidate
    Float_t Sexa_Pz;       // z-component of anti-sexaquark candidate
    Float_t Sexa_M;        // inv. mass of anti-sexaquark candidate
    Float_t Sexa_X;        // x-coordinate of anti-sexaquark candidate (secondary vertex)
    Float_t Sexa_Y;        // y-coordinate of anti-sexaquark candidate (secondary vertex)
    Float_t Sexa_Z;        // z-coordinate of anti-sexaquark candidate (secondary vertex)
    Float_t Sexa_DCA;      // distance of closest approach after extrapolation of V0s
    Bool_t Sexa_isSignal;  // kTRUE if signal, kFALSE if background

    /* First V0: "V0A" */

    Int_t V0A_Idx_Pos;             //
    Int_t V0A_Idx_Neg;             //
    Float_t V0A_Px;                //
    Float_t V0A_Py;                //
    Float_t V0A_Pz;                //
    Float_t V0A_X;                 //
    Float_t V0A_Y;                 //
    Float_t V0A_Z;                 //
    Float_t V0A_Pos_Px;            // x component of pos. daughter's momentum, in V0 vertex
    Float_t V0A_Pos_Py;            //
    Float_t V0A_Pos_Pz;            //
    Float_t V0A_Neg_Px;            //
    Float_t V0A_Neg_Py;            //
    Float_t V0A_Neg_Pz;            //
    Bool_t V0A_isSignal;           //
    Float_t V0A_E_asK0;            //
    Float_t V0A_E_asAL;            //
    Bool_t V0A_couldBeK0;          //
    Bool_t V0A_couldBeAL;          //
    Float_t V0A_Chi2;              //
    Float_t V0A_DCA_Daughters;     //
    Float_t V0A_IP_wrtPV;          //
    Float_t V0A_CPA_wrtPV;         //
    Float_t V0A_ArmAlpha;          //
    Float_t V0A_ArmPt;             //
    Float_t V0A_DecayLength;       //
    Bool_t V0A_Pos_isDuplicate;    //
    Bool_t V0A_Neg_isDuplicate;    //
    Float_t V0A_Pos_Rec_Px;        // x component of pos. daughter's momentum, as in AliESD
    Float_t V0A_Pos_Rec_Py;        //
    Float_t V0A_Pos_Rec_Pz;        //
    Float_t V0A_Neg_Rec_Px;        //
    Float_t V0A_Neg_Rec_Py;        //
    Float_t V0A_Neg_Rec_Pz;        //
    Float_t V0A_Pos_NSigmaPion;    //
    Float_t V0A_Pos_NSigmaProton;  //
    Float_t V0A_Neg_NSigmaPion;    //
    Float_t V0A_Neg_NSigmaProton;  //
    Int_t V0A_Mother_PID;          // PID of mother of V0, if available
    Int_t V0A_PID;                 // PID of V0, if available
    Int_t V0A_Pos_PID;             // PID taken for true/MC information
    Int_t V0A_Neg_PID;             //

    /* Second V0: "V0B" */

    Int_t V0B_Idx_Pos;             //
    Int_t V0B_Idx_Neg;             //
    Float_t V0B_Px;                //
    Float_t V0B_Py;                //
    Float_t V0B_Pz;                //
    Float_t V0B_X;                 //
    Float_t V0B_Y;                 //
    Float_t V0B_Z;                 //
    Float_t V0B_Pos_Px;            //
    Float_t V0B_Pos_Py;            //
    Float_t V0B_Pos_Pz;            //
    Float_t V0B_Neg_Px;            //
    Float_t V0B_Neg_Py;            //
    Float_t V0B_Neg_Pz;            //
    Bool_t V0B_isSignal;           //
    Float_t V0B_E_asK0;            //
    Float_t V0B_E_asAL;            //
    Bool_t V0B_couldBeK0;          //
    Bool_t V0B_couldBeAL;          //
    Float_t V0B_Chi2;              //
    Float_t V0B_DCA_Daughters;     //
    Float_t V0B_IP_wrtPV;          //
    Float_t V0B_CPA_wrtPV;         //
    Float_t V0B_ArmAlpha;          //
    Float_t V0B_ArmPt;             //
    Float_t V0B_DecayLength;       //
    Bool_t V0B_Pos_isDuplicate;    //
    Bool_t V0B_Neg_isDuplicate;    //
    Float_t V0B_Pos_Rec_Px;        //
    Float_t V0B_Pos_Rec_Py;        //
    Float_t V0B_Pos_Rec_Pz;        //
    Float_t V0B_Neg_Rec_Px;        //
    Float_t V0B_Neg_Rec_Py;        //
    Float_t V0B_Neg_Rec_Pz;        //
    Float_t V0B_Pos_NSigmaPion;    //
    Float_t V0B_Pos_NSigmaProton;  //
    Float_t V0B_Neg_NSigmaPion;    //
    Float_t V0B_Neg_NSigmaProton;  //
    Int_t V0B_Mother_PID;          // PID of mother of V0, if available
    Int_t V0B_PID;                 // PID of V0, if available
    Int_t V0B_Pos_PID;             // PID taken for true/MC information
    Int_t V0B_Neg_PID;             //
};

void SetSexaBranches(TTree* output_tree, Sexa_tt& this_sexa) {

    /* Identifiers */

    output_tree->Branch("RunNumber", &this_sexa.RunNumber);
    output_tree->Branch("DirNumber", &this_sexa.DirNumber);
    output_tree->Branch("Event", &this_sexa.Event);

    /* Anti-Sexaquark Candidate */

    output_tree->Branch("Idx_V0A", &this_sexa.Idx_V0A);
    output_tree->Branch("Idx_V0B", &this_sexa.Idx_V0B);
    output_tree->Branch("E", &this_sexa.Sexa_E);
    output_tree->Branch("Px", &this_sexa.Sexa_Px);
    output_tree->Branch("Py", &this_sexa.Sexa_Py);
    output_tree->Branch("Pz", &this_sexa.Sexa_Pz);
    output_tree->Branch("M", &this_sexa.Sexa_M);
    output_tree->Branch("X", &this_sexa.Sexa_X);
    output_tree->Branch("Y", &this_sexa.Sexa_Y);
    output_tree->Branch("Z", &this_sexa.Sexa_Z);
    output_tree->Branch("DCA", &this_sexa.Sexa_DCA);
    output_tree->Branch("isSignal", &this_sexa.Sexa_isSignal);

    /* First V0: "V0A" */

    output_tree->Branch("V0A_Idx_Pos", &this_sexa.V0A_Idx_Pos);
    output_tree->Branch("V0A_Idx_Neg", &this_sexa.V0A_Idx_Neg);
    output_tree->Branch("V0A_Px", &this_sexa.V0A_Px);
    output_tree->Branch("V0A_Py", &this_sexa.V0A_Py);
    output_tree->Branch("V0A_Pz", &this_sexa.V0A_Pz);
    output_tree->Branch("V0A_X", &this_sexa.V0A_X);
    output_tree->Branch("V0A_Y", &this_sexa.V0A_Y);
    output_tree->Branch("V0A_Z", &this_sexa.V0A_Z);
    output_tree->Branch("V0A_Pos_Px", &this_sexa.V0A_Pos_Px);
    output_tree->Branch("V0A_Pos_Py", &this_sexa.V0A_Pos_Py);
    output_tree->Branch("V0A_Pos_Pz", &this_sexa.V0A_Pos_Pz);
    output_tree->Branch("V0A_Neg_Px", &this_sexa.V0A_Neg_Px);
    output_tree->Branch("V0A_Neg_Py", &this_sexa.V0A_Neg_Py);
    output_tree->Branch("V0A_Neg_Pz", &this_sexa.V0A_Neg_Pz);
    output_tree->Branch("V0A_isSignal", &this_sexa.V0A_isSignal);
    output_tree->Branch("V0A_E_asK0", &this_sexa.V0A_E_asK0);
    output_tree->Branch("V0A_E_asAL", &this_sexa.V0A_E_asAL);
    output_tree->Branch("V0A_couldBeK0", &this_sexa.V0A_couldBeK0);
    output_tree->Branch("V0A_couldBeAL", &this_sexa.V0A_couldBeAL);
    output_tree->Branch("V0A_Chi2", &this_sexa.V0A_Chi2);
    output_tree->Branch("V0A_DCA_Daughters", &this_sexa.V0A_DCA_Daughters);
    output_tree->Branch("V0A_IP_wrtPV", &this_sexa.V0A_IP_wrtPV);
    output_tree->Branch("V0A_CPA_wrtPV", &this_sexa.V0A_CPA_wrtPV);
    output_tree->Branch("V0A_ArmAlpha", &this_sexa.V0A_ArmAlpha);
    output_tree->Branch("V0A_ArmPt", &this_sexa.V0A_ArmPt);
    output_tree->Branch("V0A_DecayLength", &this_sexa.V0A_DecayLength);

    output_tree->Branch("V0A_Pos_isDuplicate", &this_sexa.V0A_Pos_isDuplicate);
    output_tree->Branch("V0A_Neg_isDuplicate", &this_sexa.V0A_Neg_isDuplicate);
    output_tree->Branch("V0A_Pos_Rec_Px", &this_sexa.V0A_Pos_Rec_Px);
    output_tree->Branch("V0A_Pos_Rec_Py", &this_sexa.V0A_Pos_Rec_Py);
    output_tree->Branch("V0A_Pos_Rec_Pz", &this_sexa.V0A_Pos_Rec_Pz);
    output_tree->Branch("V0A_Neg_Rec_Px", &this_sexa.V0A_Neg_Rec_Px);
    output_tree->Branch("V0A_Neg_Rec_Py", &this_sexa.V0A_Neg_Rec_Py);
    output_tree->Branch("V0A_Neg_Rec_Pz", &this_sexa.V0A_Neg_Rec_Pz);
    output_tree->Branch("V0A_Pos_NSigmaPion", &this_sexa.V0A_Pos_NSigmaPion);
    output_tree->Branch("V0A_Pos_NSigmaProton", &this_sexa.V0A_Pos_NSigmaProton);
    output_tree->Branch("V0A_Neg_NSigmaPion", &this_sexa.V0A_Neg_NSigmaPion);
    output_tree->Branch("V0A_Neg_NSigmaProton", &this_sexa.V0A_Neg_NSigmaProton);

    output_tree->Branch("V0A_Mother_PID", &this_sexa.V0A_Mother_PID);
    output_tree->Branch("V0A_PID", &this_sexa.V0A_PID);
    output_tree->Branch("V0A_Pos_PID", &this_sexa.V0A_Pos_PID);
    output_tree->Branch("V0A_Neg_PID", &this_sexa.V0A_Neg_PID);

    /* Second V0: "V0B" */

    output_tree->Branch("V0B_Idx_Pos", &this_sexa.V0B_Idx_Pos);
    output_tree->Branch("V0B_Idx_Neg", &this_sexa.V0B_Idx_Neg);
    output_tree->Branch("V0B_Px", &this_sexa.V0B_Px);
    output_tree->Branch("V0B_Py", &this_sexa.V0B_Py);
    output_tree->Branch("V0B_Pz", &this_sexa.V0B_Pz);
    output_tree->Branch("V0B_X", &this_sexa.V0B_X);
    output_tree->Branch("V0B_Y", &this_sexa.V0B_Y);
    output_tree->Branch("V0B_Z", &this_sexa.V0B_Z);
    output_tree->Branch("V0B_Pos_Px", &this_sexa.V0B_Pos_Px);
    output_tree->Branch("V0B_Pos_Py", &this_sexa.V0B_Pos_Py);
    output_tree->Branch("V0B_Pos_Pz", &this_sexa.V0B_Pos_Pz);
    output_tree->Branch("V0B_Neg_Px", &this_sexa.V0B_Neg_Px);
    output_tree->Branch("V0B_Neg_Py", &this_sexa.V0B_Neg_Py);
    output_tree->Branch("V0B_Neg_Pz", &this_sexa.V0B_Neg_Pz);
    output_tree->Branch("V0B_isSignal", &this_sexa.V0B_isSignal);
    output_tree->Branch("V0B_E_asK0", &this_sexa.V0B_E_asK0);
    output_tree->Branch("V0B_E_asAL", &this_sexa.V0B_E_asAL);
    output_tree->Branch("V0B_couldBeK0", &this_sexa.V0B_couldBeK0);
    output_tree->Branch("V0B_couldBeAL", &this_sexa.V0B_couldBeAL);
    output_tree->Branch("V0B_Chi2", &this_sexa.V0B_Chi2);
    output_tree->Branch("V0B_DCA_Daughters", &this_sexa.V0B_DCA_Daughters);
    output_tree->Branch("V0B_IP_wrtPV", &this_sexa.V0B_IP_wrtPV);
    output_tree->Branch("V0B_CPA_wrtPV", &this_sexa.V0B_CPA_wrtPV);
    output_tree->Branch("V0B_ArmAlpha", &this_sexa.V0B_ArmAlpha);
    output_tree->Branch("V0B_ArmPt", &this_sexa.V0B_ArmPt);
    output_tree->Branch("V0B_DecayLength", &this_sexa.V0B_DecayLength);

    output_tree->Branch("V0B_Pos_isDuplicate", &this_sexa.V0B_Pos_isDuplicate);
    output_tree->Branch("V0B_Neg_isDuplicate", &this_sexa.V0B_Neg_isDuplicate);
    output_tree->Branch("V0B_Pos_Rec_Px", &this_sexa.V0B_Pos_Rec_Px);
    output_tree->Branch("V0B_Pos_Rec_Py", &this_sexa.V0B_Pos_Rec_Py);
    output_tree->Branch("V0B_Pos_Rec_Pz", &this_sexa.V0B_Pos_Rec_Pz);
    output_tree->Branch("V0B_Neg_Rec_Px", &this_sexa.V0B_Neg_Rec_Px);
    output_tree->Branch("V0B_Neg_Rec_Py", &this_sexa.V0B_Neg_Rec_Py);
    output_tree->Branch("V0B_Neg_Rec_Pz", &this_sexa.V0B_Neg_Rec_Pz);
    output_tree->Branch("V0B_Pos_NSigmaPion", &this_sexa.V0B_Pos_NSigmaPion);
    output_tree->Branch("V0B_Pos_NSigmaProton", &this_sexa.V0B_Pos_NSigmaProton);
    output_tree->Branch("V0B_Neg_NSigmaPion", &this_sexa.V0B_Neg_NSigmaPion);
    output_tree->Branch("V0B_Neg_NSigmaProton", &this_sexa.V0B_Neg_NSigmaProton);

    output_tree->Branch("V0B_Mother_PID", &this_sexa.V0B_Mother_PID);
    output_tree->Branch("V0B_PID", &this_sexa.V0B_PID);
    output_tree->Branch("V0B_Pos_PID", &this_sexa.V0B_Pos_PID);
    output_tree->Branch("V0B_Neg_PID", &this_sexa.V0B_Neg_PID);
}
