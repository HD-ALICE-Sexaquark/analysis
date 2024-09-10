#include "include/Headers.hxx"

#include "../sexaquark/AliAnalysisTaskSexaquark_Structs.h"

// clang-format off
/*
 * Indices
 *     MajorIndex
 *         Injected:    "Injected.RunNumber * 1000 + Injected.DirNumber"
 *         Events:      "Event.RunNumber * 1000 + Event.DirNumber"
 *         MCParticles: "MCParticle.RunNumber * 1000 + MCParticle.DirNumber"
 *         Tracks:      "Track.RunNumber * 1000 + Track.DirNumber"
 *         V0s:         "V0.RunNumber * 1000 + V0.DirNumber"
 *         Sexaquarks:  "Sexaquark.RunNumber * 1000 + Sexaquark.DirNumber"
 *     MinorIndex
 *         Injected:    "Injected.EventNumber * 1000 + Injected.ReactionID"
 *         Events:      "Event.Number"
 *         MCParticles: "MCParticle.EventNumber * 10000000 + MCParticle.Idx"
 *         Tracks:      "Track.EventNumber * 10000000 + Track.Idx"
 *         V0s:         "V0.EventNumber * 10000000 + V0.Idx"
 *         Sexaquarks:  "Sexaquark.EventNumber * 1000 + Sexaquark.ReactionID"
 */
// clang-format on

void Trees_ScanSignal(TString InputFileName = "AnalysisResults_indexed.root") {

    /** Input **/

    TFile* InputFile = TFile::Open(InputFileName, "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! Couldn't open file " << InputFileName << " !!" << std::endl;
        return;
    }

    std::map<TString, TTree*> Tree;
    std::vector<TString> TreeNames = {"Injected",    "Events",      "MCParticles",    "PiPluses",  "PiMinuses",
                                      "AntiProtons", "AntiLambdas", "KaonsZeroShort", "Sexaquarks"};

    for (const auto& treeName : TreeNames) {
        Tree[treeName] = (TTree*)InputFile->Get(treeName);
        if (!Tree[treeName]) {
            std::cout << "!! ERROR !! Couldn't find TTree " << treeName << " in " << InputFileName << " !!" << std::endl;
            InputFile->Close();
            return;
        }
    }

    /** Set Branches **/

    Injected_tt* thisInjected = nullptr;
    Tree["Injected"]->SetBranchAddress("Injected", &thisInjected);

    Event_tt* thisEvent = nullptr;
    Tree["Events"]->SetBranchAddress("Event", &thisEvent);

    MC_tt* thisMC = nullptr;
    Tree["MCParticles"]->SetBranchAddress("MCParticle", &thisMC);

    Track_tt* thisPiPlus = nullptr;
    Tree["PiPluses"]->SetBranchAddress("PiPlus", &thisPiPlus);

    Track_tt* thisPiMinus = nullptr;
    Tree["PiMinuses"]->SetBranchAddress("PiMinus", &thisPiMinus);

    Track_tt* thisAntiProton = nullptr;
    Tree["AntiProtons"]->SetBranchAddress("AntiProton", &thisAntiProton);

    V0_tt* thisAL = nullptr;
    Tree["AntiLambdas"]->SetBranchAddress("AntiLambda", &thisAL);

    V0_tt* thisK0S = nullptr;
    Tree["KaonsZeroShort"]->SetBranchAddress("KaonZeroShort", &thisK0S);

    RecSexaquark_aa* thisSexaquark = nullptr;
    Tree["Sexaquarks"]->SetBranchAddress("Sexaquark", &thisSexaquark);

    /*** Loop over Sexaquarks candidates***/

    Int_t nEntries = Tree["Sexaquarks"]->GetEntries();

    Int_t RunNumber;
    Int_t DirNumber;
    Int_t EventNumber;
    Int_t ReactionID;

    for (Int_t i = 0; i < nEntries; i++) {

        /* Reconstructed Sexaquark Info */

        Tree["Sexaquarks"]->GetEntry(i);

        RunNumber = thisSexaquark->RunNumber;
        DirNumber = thisSexaquark->DirNumber;
        EventNumber = thisSexaquark->EventNumber;
        ReactionID = thisSexaquark->ReactionID;

        Int_t Idx_AL = thisSexaquark->Idx_AL;
        Int_t Idx_K0S = thisSexaquark->Idx_K0S;
        Int_t Idx_AL_Neg = thisSexaquark->Idx_AL_Neg;
        Int_t Idx_AL_Pos = thisSexaquark->Idx_AL_Pos;
        Int_t Idx_K0S_Neg = thisSexaquark->Idx_K0S_Neg;
        Int_t Idx_K0S_Pos = thisSexaquark->Idx_K0S_Pos;

        std::cout << "Sexaquark : (" << i << ") " << RunNumber << " " << DirNumber << " " << EventNumber << " " << ReactionID << std::endl;
        std::cout << "            " << Idx_AL << " " << Idx_K0S << " " << Idx_AL_Neg << " " << Idx_AL_Pos << " " << Idx_K0S_Neg << " " << Idx_K0S_Pos
                  << std::endl;
        std::cout << "            " << thisSexaquark->Px << " " << thisSexaquark->Py << std::endl;

        /* Event Info */

        Int_t eventEntry = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);

        if (eventEntry >= 0) {
            Tree["Events"]->GetEntry(eventEntry);
            std::cout << ">> Event : " << thisEvent->RunNumber << " " << thisEvent->DirNumber << " " << thisEvent->Number << std::endl;
            std::cout << "           " << thisEvent->MC_Xv_PV << " " << thisEvent->MC_Yv_PV << " " << thisEvent->MC_Zv_PV << std::endl;
            std::cout << "           " << thisEvent->Rec_Xv_PV << " " << thisEvent->Rec_Yv_PV << " " << thisEvent->Rec_Zv_PV << std::endl;
            std::cout << "           " << thisEvent->Centrality << std::endl;
        } else {
            std::cout << ">> No matching event found" << std::endl;
        }

        /* Injected Info */

        Int_t injectedEntry = Tree["Injected"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 1000 + ReactionID);

        if (injectedEntry >= 0) {
            Tree["Injected"]->GetEntry(injectedEntry);
            std::cout << ">> Injected : " << thisInjected->RunNumber << " " << thisInjected->DirNumber << " " << thisInjected->EventNumber << " "
                      << thisInjected->ReactionID << std::endl;
            std::cout << "              " << thisInjected->Px << " " << thisInjected->Py << std::endl;
        } else {
            std::cout << ">> No matching injected found" << std::endl;
        }

        /* Event Info */

        RunNumber = thisInjected->RunNumber;
        DirNumber = thisInjected->DirNumber;
        EventNumber = thisInjected->EventNumber;
        eventEntry = Tree["Events"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber);

        if (eventEntry >= 0) {
            Tree["Events"]->GetEntry(eventEntry);
            std::cout << ">> Event : " << thisEvent->RunNumber << " " << thisEvent->DirNumber << " " << thisEvent->Number << std::endl;
            std::cout << "           " << thisEvent->MC_Xv_PV << " " << thisEvent->MC_Yv_PV << " " << thisEvent->MC_Zv_PV << std::endl;
            std::cout << "           " << thisEvent->Rec_Xv_PV << " " << thisEvent->Rec_Yv_PV << " " << thisEvent->Rec_Zv_PV << std::endl;
            std::cout << "           " << thisEvent->Centrality << std::endl;
        } else {
            std::cout << ">> No matching event found" << std::endl;
        }

        /* Anti-Lambda Info */

        Int_t antilambdaEntry = Tree["AntiLambdas"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_AL);
        if (antilambdaEntry >= 0) {
            Tree["AntiLambdas"]->GetEntry(antilambdaEntry);
            std::cout << ">> Rec AntiLambda : (" << thisAL->Idx << ") " << thisAL->RunNumber << " " << thisAL->DirNumber << " " << thisAL->EventNumber
                      << " " << thisAL->ReactionID << std::endl;
            Float_t AL_Mass = TMath::Sqrt(thisAL->E * thisAL->E - thisAL->Px * thisAL->Px - thisAL->Py * thisAL->Py - thisAL->Pz * thisAL->Pz);
            std::cout << "                     " << AL_Mass << std::endl;
            std::cout << "                     " << thisAL->Idx_True << std::endl;

        } else {
            std::cout << ">> No matching AntiLambda found" << std::endl;
        }

        /* Kaon Zero Short Info */

        Int_t k0sEntry = Tree["KaonsZeroShort"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_K0S);
        if (k0sEntry >= 0) {
            Tree["KaonsZeroShort"]->GetEntry(k0sEntry);
            std::cout << ">> Rec KaonZeroShort : (" << thisK0S->Idx << ") " << thisK0S->RunNumber << " " << thisK0S->DirNumber << " "
                      << thisK0S->EventNumber << " " << thisK0S->ReactionID << std::endl;
            Float_t K0_Mass =
                TMath::Sqrt(thisK0S->E * thisK0S->E - thisK0S->Px * thisK0S->Px - thisK0S->Py * thisK0S->Py - thisK0S->Pz * thisK0S->Pz);
            std::cout << "                     " << K0_Mass << std::endl;
            std::cout << "                     " << thisK0S->Idx_True << std::endl;
        } else {
            std::cout << ">> No matching KaonZeroShort found" << std::endl;
        }

        /* Anti-Lambda Pos. Dau. */

        Int_t alposdauEntry = Tree["PiPluses"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + Idx_AL_Pos);
        if (alposdauEntry >= 0) {
            Tree["PiPluses"]->GetEntry(alposdauEntry);
            std::cout << ">> Rec PiPlus (from AL) : (" << thisPiPlus->Idx << ") " << thisPiPlus->RunNumber << " " << thisPiPlus->DirNumber << " "
                      << thisPiPlus->EventNumber << " " << thisPiPlus->ReactionID << std::endl;
            std::cout << "                          " << thisPiPlus->Idx_True << std::endl;
            /*  */
            Int_t mc_alposdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + thisPiPlus->Idx_True);
            if (mc_alposdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_alposdauEntry);
                std::cout << "   >> True: (" << thisMC->Idx << ") " << thisMC->PdgCode << " " << thisMC->Status << " " << thisMC->Idx_Mother << " "
                          << thisMC->ReactionID << std::endl;
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
            std::cout << ">> Rec AntiProton (from AL) : (" << thisAntiProton->Idx << ") " << thisAntiProton->RunNumber << " "
                      << thisAntiProton->DirNumber << " " << thisAntiProton->EventNumber << " " << thisAntiProton->ReactionID << std::endl;
            std::cout << "                              " << thisAntiProton->Idx_True << std::endl;
            /*  */
            Int_t mc_alnegdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + thisAntiProton->Idx_True);
            if (mc_alnegdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_alnegdauEntry);
                std::cout << "   >> True: (" << thisMC->Idx << ") " << thisMC->PdgCode << " " << thisMC->Status << " " << thisMC->Idx_Mother << " "
                          << thisMC->ReactionID << std::endl;
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
            std::cout << ">> Rec PiPlus (from K0S) : (" << thisPiPlus->Idx << ") " << thisPiPlus->RunNumber << " " << thisPiPlus->DirNumber << " "
                      << thisPiPlus->EventNumber << " " << thisPiPlus->ReactionID << std::endl;
            std::cout << "                           " << thisPiPlus->Idx_True << std::endl;
            std::cout << "                           " << thisPiPlus->Px << " " << thisPiPlus->Py << " " << thisPiPlus->Pz << std::endl;
            std::cout << "                           " << thisPiPlus->NSigmaPion << " " << thisPiPlus->NSigmaKaon << " " << thisPiPlus->NSigmaProton
                      << std::endl;
            /*  */
            Int_t mc_k0sposdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + thisPiPlus->Idx_True);
            if (mc_k0sposdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_k0sposdauEntry);
                std::cout << "   >> True: (" << thisMC->Idx << ") " << thisMC->PdgCode << " " << thisMC->Status << " " << thisMC->Idx_Mother << " "
                          << thisMC->ReactionID << std::endl;
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
            std::cout << ">> Rec PiMinus (from K0S) : (" << thisPiMinus->Idx << ") " << thisPiMinus->RunNumber << " " << thisPiMinus->DirNumber << " "
                      << thisPiMinus->EventNumber << " " << thisPiMinus->ReactionID << std::endl;
            std::cout << "                            " << thisPiMinus->Idx_True << std::endl;
            /*  */
            Int_t mc_k0snegdauEntry =
                Tree["MCParticles"]->GetEntryNumberWithIndex(RunNumber * 1000 + DirNumber, EventNumber * 10000000 + thisPiMinus->Idx_True);
            if (mc_k0snegdauEntry >= 0) {
                Tree["MCParticles"]->GetEntry(mc_k0snegdauEntry);
                std::cout << "   >> True: (" << thisMC->Idx << ") " << thisMC->PdgCode << " " << thisMC->Status << " " << thisMC->Idx_Mother << " "
                          << thisMC->ReactionID << std::endl;
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
