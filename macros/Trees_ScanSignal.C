#include "include/Headers.hxx"

// Indices
//     MajorIndex
//         All Trees : "RunNumber * 1000 + DirNumber"
//     MinorIndex
//         Events                   : "EventNumber"
//         Injected, Sexaquarks     : "EventNumber * 1000 + ReactionID"
//         MCParticles, Tracks, V0s : "EventNumber * 10000000 + Idx"

/*
 * Process an indexed `AnalysisResults.root` file. Scan reconstructed anti-sexaquark candidates and their decay products.
 */
void Trees_ScanSignal(TString InputFileName = "AnalysisResults_indexed.root") {

    /** Input **/

    TFile *InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << InputFileName << " !!" << std::endl;
        return;
    }

    std::map<TString, TTree *> Tree;
    std::vector<TString> TreeNames = {"Injected",    "Events",      "MCParticles",    "PiPluses",  "PiMinuses",
                                      "AntiProtons", "AntiLambdas", "KaonsZeroShort", "Sexaquarks"};

    for (const auto &treeName : TreeNames) {
        Tree[treeName] = (TTree *)InputFile->Get(treeName);
        if (!Tree[treeName]) {
            std::cout << "!! ERROR !! Couldn't find TTree " << treeName << " in " << InputFileName << " !!" << std::endl;
            InputFile->Close();
            return;
        }
    }

    /** Set Branches **/

    /* Sexaquarks */

    Int_t Sexaquark_RunNumber;
    Int_t Sexaquark_DirNumber;
    Int_t Sexaquark_EventNumber;
    Int_t Sexaquark_ReactionID;
    Int_t Sexaquark_Idx_AL;
    Int_t Sexaquark_Idx_K0S;
    Int_t Sexaquark_Idx_AL_Neg;
    Int_t Sexaquark_Idx_AL_Pos;
    Int_t Sexaquark_Idx_K0S_Neg;
    Int_t Sexaquark_Idx_K0S_Pos;
    Float_t Sexaquark_Px;
    Float_t Sexaquark_Py;

    Tree["Sexaquarks"]->SetBranchAddress("RunNumber", &Sexaquark_RunNumber);
    Tree["Sexaquarks"]->SetBranchAddress("DirNumber", &Sexaquark_DirNumber);
    Tree["Sexaquarks"]->SetBranchAddress("EventNumber", &Sexaquark_EventNumber);
    Tree["Sexaquarks"]->SetBranchAddress("ReactionID", &Sexaquark_ReactionID);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_AL", &Sexaquark_Idx_AL);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_K0S", &Sexaquark_Idx_K0S);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_AL_Neg", &Sexaquark_Idx_AL_Neg);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_AL_Pos", &Sexaquark_Idx_AL_Pos);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_K0S_Neg", &Sexaquark_Idx_K0S_Neg);
    Tree["Sexaquarks"]->SetBranchAddress("Idx_K0S_Pos", &Sexaquark_Idx_K0S_Pos);
    Tree["Sexaquarks"]->SetBranchAddress("Px", &Sexaquark_Px);
    Tree["Sexaquarks"]->SetBranchAddress("Py", &Sexaquark_Py);

    /* Events */

    Int_t Event_RunNumber;
    Int_t Event_DirNumber;
    Int_t Event_EventNumber;
    Float_t Event_MC_Xv_PV;
    Float_t Event_MC_Yv_PV;
    Float_t Event_MC_Zv_PV;
    Float_t Event_Rec_Xv_PV;
    Float_t Event_Rec_Yv_PV;
    Float_t Event_Rec_Zv_PV;
    Float_t Event_Centrality;

    Tree["Events"]->SetBranchAddress("RunNumber", &Event_RunNumber);
    Tree["Events"]->SetBranchAddress("DirNumber", &Event_DirNumber);
    Tree["Events"]->SetBranchAddress("EventNumber", &Event_EventNumber);
    Tree["Events"]->SetBranchAddress("MC_Xv_PV", &Event_MC_Xv_PV);
    Tree["Events"]->SetBranchAddress("MC_Yv_PV", &Event_MC_Yv_PV);
    Tree["Events"]->SetBranchAddress("MC_Zv_PV", &Event_MC_Zv_PV);
    Tree["Events"]->SetBranchAddress("Rec_Xv_PV", &Event_Rec_Xv_PV);
    Tree["Events"]->SetBranchAddress("Rec_Yv_PV", &Event_Rec_Yv_PV);
    Tree["Events"]->SetBranchAddress("Rec_Zv_PV", &Event_Rec_Zv_PV);
    Tree["Events"]->SetBranchAddress("Centrality", &Event_Centrality);

    /* Injected */

    Int_t Injected_RunNumber;
    Int_t Injected_DirNumber;
    Int_t Injected_EventNumber;
    Int_t Injected_ReactionID;
    Float_t Injected_Px;
    Float_t Injected_Py;

    Tree["Injected"]->SetBranchAddress("RunNumber", &Injected_RunNumber);
    Tree["Injected"]->SetBranchAddress("DirNumber", &Injected_DirNumber);
    Tree["Injected"]->SetBranchAddress("EventNumber", &Injected_EventNumber);
    Tree["Injected"]->SetBranchAddress("ReactionID", &Injected_ReactionID);
    Tree["Injected"]->SetBranchAddress("Px", &Injected_Px);
    Tree["Injected"]->SetBranchAddress("Py", &Injected_Py);

    /* MC Particles */

    Int_t MC_Idx;
    Int_t MC_PdgCode;
    Int_t MC_Status;
    Int_t MC_Idx_Mother;
    Int_t MC_ReactionID;

    Tree["MCParticles"]->SetBranchAddress("Idx", &MC_Idx);
    Tree["MCParticles"]->SetBranchAddress("PdgCode", &MC_PdgCode);
    Tree["MCParticles"]->SetBranchAddress("Status", &MC_Status);
    Tree["MCParticles"]->SetBranchAddress("Idx_Mother", &MC_Idx_Mother);
    Tree["MCParticles"]->SetBranchAddress("ReactionID", &MC_ReactionID);

    /* AntiLambdas */

    Int_t AL_Idx;
    Int_t AL_RunNumber;
    Int_t AL_DirNumber;
    Int_t AL_EventNumber;
    Int_t AL_ReactionID;
    Int_t AL_Idx_True;
    Float_t AL_Px;
    Float_t AL_Py;
    Float_t AL_Pz;
    Float_t AL_E;

    Tree["AntiLambdas"]->SetBranchAddress("Idx", &AL_Idx);
    Tree["AntiLambdas"]->SetBranchAddress("RunNumber", &AL_RunNumber);
    Tree["AntiLambdas"]->SetBranchAddress("DirNumber", &AL_DirNumber);
    Tree["AntiLambdas"]->SetBranchAddress("EventNumber", &AL_EventNumber);
    Tree["AntiLambdas"]->SetBranchAddress("ReactionID", &AL_ReactionID);
    Tree["AntiLambdas"]->SetBranchAddress("Idx_True", &AL_Idx_True);
    Tree["AntiLambdas"]->SetBranchAddress("Px", &AL_Px);
    Tree["AntiLambdas"]->SetBranchAddress("Py", &AL_Py);
    Tree["AntiLambdas"]->SetBranchAddress("Pz", &AL_Pz);
    Tree["AntiLambdas"]->SetBranchAddress("E", &AL_E);

    /* KaonsZeroShort */

    Int_t K0S_Idx;
    Int_t K0S_RunNumber;
    Int_t K0S_DirNumber;
    Int_t K0S_EventNumber;
    Int_t K0S_ReactionID;
    Int_t K0S_Idx_True;
    Float_t K0S_Px;
    Float_t K0S_Py;
    Float_t K0S_Pz;
    Float_t K0S_E;

    Tree["KaonsZeroShort"]->SetBranchAddress("Idx", &K0S_Idx);
    Tree["KaonsZeroShort"]->SetBranchAddress("RunNumber", &K0S_RunNumber);
    Tree["KaonsZeroShort"]->SetBranchAddress("DirNumber", &K0S_DirNumber);
    Tree["KaonsZeroShort"]->SetBranchAddress("EventNumber", &K0S_EventNumber);
    Tree["KaonsZeroShort"]->SetBranchAddress("ReactionID", &K0S_ReactionID);
    Tree["KaonsZeroShort"]->SetBranchAddress("Idx_True", &K0S_Idx_True);
    Tree["KaonsZeroShort"]->SetBranchAddress("Px", &K0S_Px);
    Tree["KaonsZeroShort"]->SetBranchAddress("Py", &K0S_Py);
    Tree["KaonsZeroShort"]->SetBranchAddress("Pz", &K0S_Pz);
    Tree["KaonsZeroShort"]->SetBranchAddress("E", &K0S_E);

    /* PiPluses */

    Int_t PiPlus_Idx;
    Int_t PiPlus_RunNumber;
    Int_t PiPlus_DirNumber;
    Int_t PiPlus_EventNumber;
    Int_t PiPlus_ReactionID;
    Int_t PiPlus_Idx_True;
    Float_t PiPlus_Px;
    Float_t PiPlus_Py;
    Float_t PiPlus_Pz;
    Float_t PiPlus_NSigmaPion;
    Float_t PiPlus_NSigmaKaon;
    Float_t PiPlus_NSigmaProton;

    Tree["PiPluses"]->SetBranchAddress("Idx", &PiPlus_Idx);
    Tree["PiPluses"]->SetBranchAddress("RunNumber", &PiPlus_RunNumber);
    Tree["PiPluses"]->SetBranchAddress("DirNumber", &PiPlus_DirNumber);
    Tree["PiPluses"]->SetBranchAddress("EventNumber", &PiPlus_EventNumber);
    Tree["PiPluses"]->SetBranchAddress("ReactionID", &PiPlus_ReactionID);
    Tree["PiPluses"]->SetBranchAddress("Idx_True", &PiPlus_Idx_True);
    Tree["PiPluses"]->SetBranchAddress("Px", &PiPlus_Px);
    Tree["PiPluses"]->SetBranchAddress("Py", &PiPlus_Py);
    Tree["PiPluses"]->SetBranchAddress("Pz", &PiPlus_Pz);
    Tree["PiPluses"]->SetBranchAddress("NSigmaPion", &PiPlus_NSigmaPion);
    Tree["PiPluses"]->SetBranchAddress("NSigmaKaon", &PiPlus_NSigmaKaon);
    Tree["PiPluses"]->SetBranchAddress("NSigmaProton", &PiPlus_NSigmaProton);

    /* PiMinuses */

    Int_t PiMinus_Idx;
    Int_t PiMinus_RunNumber;
    Int_t PiMinus_DirNumber;
    Int_t PiMinus_EventNumber;
    Int_t PiMinus_ReactionID;
    Int_t PiMinus_Idx_True;
    Float_t PiMinus_Px;
    Float_t PiMinus_Py;
    Float_t PiMinus_Pz;

    Tree["PiMinuses"]->SetBranchAddress("Idx", &PiMinus_Idx);
    Tree["PiMinuses"]->SetBranchAddress("RunNumber", &PiMinus_RunNumber);
    Tree["PiMinuses"]->SetBranchAddress("DirNumber", &PiMinus_DirNumber);
    Tree["PiMinuses"]->SetBranchAddress("EventNumber", &PiMinus_EventNumber);
    Tree["PiMinuses"]->SetBranchAddress("ReactionID", &PiMinus_ReactionID);
    Tree["PiMinuses"]->SetBranchAddress("Idx_True", &PiMinus_Idx_True);
    Tree["PiMinuses"]->SetBranchAddress("Px", &PiMinus_Px);
    Tree["PiMinuses"]->SetBranchAddress("Py", &PiMinus_Py);
    Tree["PiMinuses"]->SetBranchAddress("Pz", &PiMinus_Pz);

    /* AntiProtons */

    Int_t AntiProton_Idx;
    Int_t AntiProton_RunNumber;
    Int_t AntiProton_DirNumber;
    Int_t AntiProton_EventNumber;
    Int_t AntiProton_ReactionID;
    Int_t AntiProton_Idx_True;
    Float_t AntiProton_Px;
    Float_t AntiProton_Py;
    Float_t AntiProton_Pz;
    Float_t AntiProton_E;

    Tree["AntiProtons"]->SetBranchAddress("Idx", &AntiProton_Idx);
    Tree["AntiProtons"]->SetBranchAddress("RunNumber", &AntiProton_RunNumber);
    Tree["AntiProtons"]->SetBranchAddress("DirNumber", &AntiProton_DirNumber);
    Tree["AntiProtons"]->SetBranchAddress("EventNumber", &AntiProton_EventNumber);
    Tree["AntiProtons"]->SetBranchAddress("ReactionID", &AntiProton_ReactionID);
    Tree["AntiProtons"]->SetBranchAddress("Idx_True", &AntiProton_Idx_True);
    Tree["AntiProtons"]->SetBranchAddress("Px", &AntiProton_Px);
    Tree["AntiProtons"]->SetBranchAddress("Py", &AntiProton_Py);
    Tree["AntiProtons"]->SetBranchAddress("Pz", &AntiProton_Pz);

    /*** Loop over Sexaquarks candidates***/

    Int_t nEntries = Tree["Sexaquarks"]->GetEntries();

    Int_t RunNumber;
    Int_t DirNumber;
    Int_t EventNumber;
    Int_t ReactionID;

    for (Int_t i = 0; i < nEntries; i++) {

        /* Reconstructed Sexaquark Info */

        Tree["Sexaquarks"]->GetEntry(i);

        RunNumber = Sexaquark_RunNumber;
        DirNumber = Sexaquark_DirNumber;
        EventNumber = Sexaquark_EventNumber;
        ReactionID = Sexaquark_ReactionID;

        Int_t Idx_AL = Sexaquark_Idx_AL;
        Int_t Idx_K0S = Sexaquark_Idx_K0S;
        Int_t Idx_AL_Neg = Sexaquark_Idx_AL_Neg;
        Int_t Idx_AL_Pos = Sexaquark_Idx_AL_Pos;
        Int_t Idx_K0S_Neg = Sexaquark_Idx_K0S_Neg;
        Int_t Idx_K0S_Pos = Sexaquark_Idx_K0S_Pos;

        std::cout << "Sexaquark : (" << i << ") " << RunNumber << " " << DirNumber << " " << EventNumber << " " << ReactionID << std::endl;
        std::cout << "            " << Idx_AL << " " << Idx_K0S << " " << Idx_AL_Neg << " " << Idx_AL_Pos << " " << Idx_K0S_Neg << " " << Idx_K0S_Pos
                  << std::endl;
        std::cout << "            " << Sexaquark_Px << " " << Sexaquark_Py << std::endl;

        /* Event Info */

        Int_t eventEntry = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);

        if (eventEntry >= 0) {
            Tree["Events"]->GetEntry(eventEntry);
            std::cout << ">> Event : " << Event_RunNumber << " " << Event_DirNumber << " " << Event_EventNumber << std::endl;
            std::cout << "           " << Event_MC_Xv_PV << " " << Event_MC_Yv_PV << " " << Event_MC_Zv_PV << std::endl;
            std::cout << "           " << Event_Rec_Xv_PV << " " << Event_Rec_Yv_PV << " " << Event_Rec_Zv_PV << std::endl;
            std::cout << "           " << Event_Centrality << std::endl;
        } else {
            std::cout << ">> No matching event found" << std::endl;
        }

        /* Injected Info */

        Int_t injectedEntry = Tree["Injected"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 1000 + ReactionID);

        if (injectedEntry >= 0) {
            Tree["Injected"]->GetEntry(injectedEntry);
            std::cout << ">> Injected : " << Injected_RunNumber << " " << Injected_DirNumber << " " << Injected_EventNumber << " "
                      << Injected_ReactionID << std::endl;
            std::cout << "              " << Injected_Px << " " << Injected_Py << std::endl;
        } else {
            std::cout << ">> No matching injected found" << std::endl;
        }

        /* Event Info */

        RunNumber = Injected_RunNumber;
        DirNumber = Injected_DirNumber;
        EventNumber = Injected_EventNumber;
        eventEntry = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);

        if (eventEntry >= 0) {
            Tree["Events"]->GetEntry(eventEntry);
            std::cout << ">> Event : " << Event_RunNumber << " " << Event_DirNumber << " " << Event_EventNumber << std::endl;
            std::cout << "           " << Event_MC_Xv_PV << " " << Event_MC_Yv_PV << " " << Event_MC_Zv_PV << std::endl;
            std::cout << "           " << Event_Rec_Xv_PV << " " << Event_Rec_Yv_PV << " " << Event_Rec_Zv_PV << std::endl;
            std::cout << "           " << Event_Centrality << std::endl;
        } else {
            std::cout << ">> No matching event found" << std::endl;
        }

        /* Anti-Lambda Info */

        Int_t antilambdaEntry = Tree["AntiLambdas"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_AL);
        if (antilambdaEntry >= 0) {
            Tree["AntiLambdas"]->GetEntry(antilambdaEntry);
            std::cout << ">> Rec AntiLambda : (" << AL_Idx << ") " << AL_RunNumber << " " << AL_DirNumber << " " << AL_EventNumber << " "
                      << AL_ReactionID << std::endl;
            Float_t AL_Mass = TMath::Sqrt(AL_E * AL_E - AL_Px * AL_Px - AL_Py * AL_Py - AL_Pz * AL_Pz);
            std::cout << "                     " << AL_Mass << std::endl;
            std::cout << "                     " << AL_Idx_True << std::endl;

        } else {
            std::cout << ">> No matching AntiLambda found" << std::endl;
        }

        /* Kaon Zero Short Info */

        Int_t k0sEntry = Tree["KaonsZeroShort"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_K0S);
        if (k0sEntry >= 0) {
            Tree["KaonsZeroShort"]->GetEntry(k0sEntry);
            std::cout << ">> Rec KaonZeroShort : (" << K0S_Idx << ") " << K0S_RunNumber << " " << K0S_DirNumber << " " << K0S_EventNumber << " "
                      << K0S_ReactionID << std::endl;
            Float_t K0_Mass = TMath::Sqrt(K0S_E * K0S_E - K0S_Px * K0S_Px - K0S_Py * K0S_Py - K0S_Pz * K0S_Pz);
            std::cout << "                     " << K0_Mass << std::endl;
            std::cout << "                     " << K0S_Idx_True << std::endl;
        } else {
            std::cout << ">> No matching KaonZeroShort found" << std::endl;
        }

        /* Anti-Lambda Pos. Dau. */

        Int_t alposdauEntry = Tree["PiPluses"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_AL_Pos);
        if (alposdauEntry >= 0) {
            Tree["PiPluses"]->GetEntry(alposdauEntry);
            std::cout << ">> Rec PiPlus (from AL) : (" << PiPlus_Idx << ") " << PiPlus_RunNumber << " " << PiPlus_DirNumber << " "
                      << PiPlus_EventNumber << " " << PiPlus_ReactionID << std::endl;
            std::cout << "                          " << PiPlus_Idx_True << std::endl;
            /*  */
            Int_t mc_alposdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + PiPlus_Idx_True);
            if (mc_alposdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_alposdauEntry);
                std::cout << "   >> True: (" << MC_Idx << ") " << MC_PdgCode << " " << MC_Status << " " << MC_Idx_Mother << " " << MC_ReactionID
                          << std::endl;
            } else {
                std::cout << "   >> No matching true MC particle found" << std::endl;
            }
        } else {
            std::cout << ">> No matching PiPlus (AL) found" << std::endl;
        }

        /* Anti-Lambda Neg. Dau. */

        Int_t alnegdauEntry = Tree["AntiProtons"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_AL_Neg);
        if (alnegdauEntry >= 0) {
            Tree["AntiProtons"]->GetEntry(alnegdauEntry);
            std::cout << ">> Rec AntiProton (from AL) : (" << AntiProton_Idx << ") " << AntiProton_RunNumber << " " << AntiProton_DirNumber << " "
                      << AntiProton_EventNumber << " " << AntiProton_ReactionID << std::endl;
            std::cout << "                              " << AntiProton_Idx_True << std::endl;
            /*  */
            Int_t mc_alnegdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + AntiProton_Idx_True);
            if (mc_alnegdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_alnegdauEntry);
                std::cout << "   >> True: (" << MC_Idx << ") " << MC_PdgCode << " " << MC_Status << " " << MC_Idx_Mother << " " << MC_ReactionID
                          << std::endl;
            } else {
                std::cout << "   >> No matching true MC particle found" << std::endl;
            }
        } else {
            std::cout << ">> No matching AntiProton (AL) found" << std::endl;
        }

        /* Kaon-Zero-Short Pos. Dau. */

        Int_t k0sposdauEntry = Tree["PiPluses"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_K0S_Pos);
        if (k0sposdauEntry >= 0) {
            Tree["PiPluses"]->GetEntry(k0sposdauEntry);
            std::cout << ">> Rec PiPlus (from K0S) : (" << PiPlus_Idx << ") " << PiPlus_RunNumber << " " << PiPlus_DirNumber << " "
                      << PiPlus_EventNumber << " " << PiPlus_ReactionID << std::endl;
            std::cout << "                           " << PiPlus_Idx_True << std::endl;
            std::cout << "                           " << PiPlus_Px << " " << PiPlus_Py << " " << PiPlus_Pz << std::endl;
            std::cout << "                           " << PiPlus_NSigmaPion << " " << PiPlus_NSigmaKaon << " " << PiPlus_NSigmaProton << std::endl;
            /*  */
            Int_t mc_k0sposdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + PiPlus_Idx_True);
            if (mc_k0sposdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_k0sposdauEntry);
                std::cout << "   >> True: (" << MC_Idx << ") " << MC_PdgCode << " " << MC_Status << " " << MC_Idx_Mother << " " << MC_ReactionID
                          << std::endl;
            } else {
                std::cout << "   >> No matching true MC particle found" << std::endl;
            }
        } else {
            std::cout << ">> No matching PiPlus (K0S) found" << std::endl;
        }

        /* Kaon-Zero-Short Neg. Dau. */

        Int_t k0snegdauEntry = Tree["PiMinuses"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_K0S_Neg);
        if (k0snegdauEntry >= 0) {
            Tree["PiMinuses"]->GetEntry(k0snegdauEntry);
            std::cout << ">> Rec PiMinus (from K0S) : (" << PiMinus_Idx << ") " << PiMinus_RunNumber << " " << PiMinus_DirNumber << " "
                      << PiMinus_EventNumber << " " << PiMinus_ReactionID << std::endl;
            std::cout << "                            " << PiMinus_Idx_True << std::endl;
            /*  */
            Int_t mc_k0snegdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + PiMinus_Idx_True);
            if (mc_k0snegdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_k0snegdauEntry);
                std::cout << "   >> True: (" << MC_Idx << ") " << MC_PdgCode << " " << MC_Status << " " << MC_Idx_Mother << " " << MC_ReactionID
                          << std::endl;
            } else {
                std::cout << "   >> No matching true MC particle found" << std::endl;
            }
        } else {
            std::cout << ">> No matching PiMinus (K0S) found" << std::endl;
        }

        std::cout << std::endl;
    }

    InputFile->Close();
}
