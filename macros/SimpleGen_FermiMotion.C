#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

void SimpleGen_FermiMotion(Int_t N = 25000) {

    // prepare random numbers
    TRandom3* Rndm = new TRandom3(0);

    // define central value and error, obtained from:
    // -- Povh, Rith, Scholz, Zetsche, Rodejohann.
    //    "Particles and Nuclei: An Introduction to the Physical Concepts".
    //    (Springer, 2015, 7th edition)
    Double_t fFermiMomentum = 0.250;       // GeV
    Double_t fFermiMomentumError = 0.005;  // GeV

    // define function
    TF1* fFermiMomentumModel = new TF1("Fermi Momentum Model", "[0]*exp(-0.5*((x-[1])/[2])**2) ",  //
                                       fFermiMomentum - 5 * fFermiMomentumError, fFermiMomentum + 5 * fFermiMomentumError);
    fFermiMomentumModel->SetParameter(0, 1.);
    fFermiMomentumModel->SetParameter(1, fFermiMomentum);
    fFermiMomentumModel->SetParameter(2, fFermiMomentumError);

    Int_t i;
    Float_t Px, Py, Pz;

    // create a ROOT file
    TFile* file = new TFile("SimpleGen_FermiMotion.root", "RECREATE");

    // create a TTree
    TTree* tree = new TTree("IC", "SimpleGen_FermiMotion");
    tree->Branch("ID", &i, "ID/I");
    tree->Branch("NPx", &Px, "Px/F");
    tree->Branch("NPy", &Py, "Py/F");
    tree->Branch("NPz", &Pz, "Pz/F");

    for (i = 0; i < N; i++) {

        Float_t P_N = fFermiMomentumModel->GetRandom();
        // given the total momentum, set phi,theta as uniform and uncorrelated variables
        Float_t random_N[2];
        Rndm->RndmArray(2, random_N);
        Float_t Phi_N = 2 * TMath::Pi() * random_N[0];       // azimuthal angle (uniform distribution) (in radians)
        Float_t Theta_N = TMath::ACos(1 - 2 * random_N[1]);  // polar angle (based on inverse CDF) (in radians)
        // then, assign px,py,pz
        Px = P_N * TMath::Cos(Phi_N) * TMath::Sin(Theta_N);
        Py = P_N * TMath::Sin(Phi_N) * TMath::Sin(Theta_N);
        Pz = P_N * TMath::Cos(Theta_N);

        tree->Fill();
    }

    // write the tree to the file
    tree->Write();

    // close the file
    file->Save();
    file->Close();
}
