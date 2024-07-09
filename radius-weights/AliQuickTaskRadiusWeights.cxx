#include "AliQuickTaskRadiusWeights.h"

ClassImp(AliQuickTaskRadiusWeights);

/*
 Empty I/O constructor. Non-persistent members are initialized to their default values from here.
*/
AliQuickTaskRadiusWeights::AliQuickTaskRadiusWeights() : AliAnalysisTaskSE(), fOutputListOfHists(0), fMC(0), fMC_PrimaryVertex(0) {}

/*
 Constructor, called locally.
*/
AliQuickTaskRadiusWeights::AliQuickTaskRadiusWeights(const char* name)
    : AliAnalysisTaskSE(name), fOutputListOfHists(0), fMC(0), fMC_PrimaryVertex(0) {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

/*
 Destructor.
*/
AliQuickTaskRadiusWeights::~AliQuickTaskRadiusWeights() {
    if (fOutputListOfHists) {
        delete fOutputListOfHists;
    }
}

/*
 Create output objects, called once at RUNTIME ~ execution on Grid
*/
void AliQuickTaskRadiusWeights::UserCreateOutputObjects() {

    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    if (!man) {
        AliFatal("ERROR: AliAnalysisManager couldn't be found.");
    }

    AliESDInputHandler* inputHandler = (AliESDInputHandler*)(man->GetInputEventHandler());
    if (!inputHandler) {
        AliFatal("ERROR: AliESDInputHandler couldn't be found.");
    }

    /** Prepare output histograms **/

    fOutputListOfHists = new TList();
    fOutputListOfHists->SetOwner(kTRUE); /* the TList destructor should delete all objects added to it */

    PrepareHistograms();

    PostData(1, fOutputListOfHists);
}

/*
 */
void AliQuickTaskRadiusWeights::PrepareHistograms() {

    TString stage = "MCGen_PhotonConversions";

    fHist_MCGen_SFM_Radius = new TH1F(stage + "_Radius", "", 360, 0., 180.);
    fHist_MCGen_SFM_Radius->Sumw2();
    fOutputListOfHists->Add(fHist_MCGen_SFM_Radius);

    fHist_MCGen_SFM_YvsX = new TH2F(stage + "_YvsX", "", 600, -600., 600., 600, -600., 600.);
    fOutputListOfHists->Add(fHist_MCGen_SFM_YvsX);
}

/*
 Main function, called per each event at RUNTIME ~ execution on Grid
*/
void AliQuickTaskRadiusWeights::UserExec(Option_t*) {

    fMC = MCEvent();
    if (!fMC) AliFatal("ERROR: AliMCEvent couldn't be found.");
    fMC_PrimaryVertex = const_cast<AliVVertex*>(fMC->GetPrimaryVertex());

    if (TMath::Abs(fMC_PrimaryVertex->GetZ()) > 10.) {
        PostData(1, fOutputListOfHists);
        return;
    }

    ProcessMCGen();

    PostData(1, fOutputListOfHists);
}

/*
 Loop over MC particles in a single event.
 - Uses: `fMC`, `fMC_PrimaryVertex`
*/
void AliQuickTaskRadiusWeights::ProcessMCGen() {

    AliMCParticle* mcPart;
    Int_t pdg_mc;

    AliMCParticle* mcDaughter;
    Int_t pdg_dau;

    Int_t n_daughters;
    Int_t mcIdxNegDaughter;
    Int_t mcIdxPosDaughter;

    /* Loop over MC gen. particles in a single event */

    for (Int_t mcIdx = 0; mcIdx < fMC->GetNumberOfTracks(); mcIdx++) {

        mcPart = (AliMCParticle*)fMC->GetTrack(mcIdx);

        /* Protection against particles with near-zero momentum */

        if (mcPart->P() < 1E-6) continue;

        pdg_mc = mcPart->PdgCode();

        /* Search for photon conversions to fill tomography histograms */

        Double_t conv_radius;
        if (pdg_mc == 22) {
            n_daughters = mcPart->GetNDaughters();
            mcIdxNegDaughter = 0;
            mcIdxPosDaughter = 0;
            if (n_daughters && mcPart->IsPhysicalPrimary() && mcPart->Pt() > 0.05 && TMath::Abs(mcPart->Eta()) < 0.9) {
                for (Int_t mcIdxDaughter = mcPart->GetDaughterFirst(); mcIdxDaughter <= mcPart->GetDaughterLast(); mcIdxDaughter++) {
                    mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxDaughter);
                    pdg_dau = mcDaughter->PdgCode();
                    if (pdg_dau == 11 && mcDaughter->IsSecondaryFromMaterial()) mcIdxNegDaughter = mcIdxDaughter;
                    if (pdg_dau == -11 && mcDaughter->IsSecondaryFromMaterial()) mcIdxPosDaughter = mcIdxDaughter;
                }
                if (n_daughters && mcIdxNegDaughter && mcIdxPosDaughter) {
                    mcDaughter = (AliMCParticle*)fMC->GetTrack(mcIdxNegDaughter);
                    conv_radius = TMath::Sqrt(mcDaughter->Xv() * mcDaughter->Xv() + mcDaughter->Yv() * mcDaughter->Yv());
                    if (5. < conv_radius && conv_radius < 180.) {
                        fHist_MCGen_SFM_Radius->Fill(conv_radius);
                        fHist_MCGen_SFM_YvsX->Fill(mcDaughter->Xv(), mcDaughter->Yv());
                    }
                }
            }
        }

    }  // end of loop over MC particles
}
