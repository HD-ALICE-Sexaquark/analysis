#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"

void SimpleGen_Reweighting(Int_t N = 25000) {

    // prepare random numbers
    TRandom3* Rndm = new TRandom3(0);

    // limits
    Float_t fSexaquarkMass = 1.8;
    Float_t fPtMin = 0.;
    Float_t fPtMax = 5.;
    Float_t fYMin = -0.8;
    Float_t fYMax = 0.8;
    Float_t fPhiMin = 0.;
    Float_t fPhiMax = 2 * TMath::Pi();
    Float_t fRadiusMin = 5.;
    Float_t fRadiusMax = 180.;

    Float_t fPV_Xv_Min = -0.5;
    Float_t fPV_Xv_Max = 0.5;
    Float_t fPV_Yv_Min = -0.5;
    Float_t fPV_Yv_Max = 0.5;
    Float_t fPV_Zv_Min = -30.;
    Float_t fPV_Zv_Max = 30.;

    // sexaquark and vertices
    Float_t Px, Py, Pz;
    Float_t Pt, Mt, Y, Phi, E, Theta;
    TLorentzVector Sexaquark;

    Float_t PV_Xv, PV_Yv, PV_Zv;
    TVector3 PrimaryVertex;

    Float_t SV_Xv_BD, SV_Yv_BD, SV_Zv_BD;
    Float_t SV_Radius3D_BD;
    Float_t SV_Radius_BD;
    TVector3 SecondaryVertex_BD;

    Float_t SV_Xv_AD, SV_Yv_AD, SV_Zv_AD;
    Float_t SV_Radius_AD;
    TVector3 SecondaryVertex_AD;

    // create a ROOT file
    TFile* file = new TFile("SimpleGen_Reweighting.root", "RECREATE");

    /* Prepare Tree */

    TTree* tree = new TTree("Tree", "Tree");

    tree->Branch("Px", &Px, "Px/F");
    tree->Branch("Py", &Py, "Py/F");
    tree->Branch("Pz", &Pz, "Pz/F");
    tree->Branch("PV_X", &PV_Xv, "PV_X/F");
    tree->Branch("PV_Y", &PV_Yv, "PV_Y/F");
    tree->Branch("PV_Z", &PV_Zv, "PV_Z/F");
    tree->Branch("SV_X_BD", &SV_Xv_BD, "SV_X/F");
    tree->Branch("SV_Y_BD", &SV_Yv_BD, "SV_Y/F");
    tree->Branch("SV_Z_BD", &SV_Zv_BD, "SV_Z/F");
    tree->Branch("SV_Radius_BD", &SV_Radius_BD, "SV_Radius_BD/F");
    tree->Branch("SV_X_AD", &SV_Xv_AD, "SV_X/F");
    tree->Branch("SV_Y_AD", &SV_Yv_AD, "SV_Y/F");
    tree->Branch("SV_Z_AD", &SV_Zv_AD, "SV_Z/F");
    tree->Branch("SV_Radius_AD", &SV_Radius_AD, "SV_Radius_AD/F");

    for (Int_t i = 0; i < N; i++) {

        Float_t random_N[7];
        Rndm->RndmArray(7, random_N);

        // sexaquark
        Pt = fPtMin + random_N[0] * (fPtMax - fPtMin);
        Mt = TMath::Sqrt(fSexaquarkMass * fSexaquarkMass + Pt * Pt);  // transverse mass (derived from inv. mass and Pt)
        Y = fYMin + random_N[1] * (fYMax - fYMin);                    // rapidity (uniform distribution)
        Phi = fPhiMin + random_N[2] * (fPhiMax - fPhiMin);            // azimuthal angle (uniform distribution) (in radians)
        Px = Pt * TMath::Cos(Phi);                                    // Px
        Py = Pt * TMath::Sin(Phi);                                    // Py
        Pz = Mt * TMath::SinH(Y);                                     // longitudinal momentum = Pz (derived from trans. mass and Y)
        E = Mt * TMath::CosH(Y);                                      // energy (derived from trans. mass and Y)

        Sexaquark.SetPxPyPzE(Px, Py, Pz, E);

        // primary vertex
        PV_Xv = fPV_Xv_Min + random_N[3] * (fPV_Xv_Max - fPV_Xv_Min);
        PV_Yv = fPV_Yv_Min + random_N[4] * (fPV_Yv_Max - fPV_Yv_Min);
        PV_Zv = fPV_Zv_Min + random_N[5] * (fPV_Zv_Max - fPV_Zv_Min);

        PrimaryVertex.SetXYZ(PV_Xv, PV_Yv, PV_Zv);

        // secondary vertex (before displacement)
        SV_Radius_BD = fRadiusMin + random_N[6] * (fRadiusMax - fRadiusMin);                      // 2D radius (uniform distribution) (in cm)
        Theta = Sexaquark.Theta();                                                                // polar angle (in radians)
        SV_Zv_BD = SV_Radius_BD / TMath::Tan(Theta);                                              // z, in both cartesian and cylindrical coordinates
        SV_Radius3D_BD = TMath::Sqrt(TMath::Power(SV_Radius_BD, 2) + TMath::Power(SV_Zv_BD, 2));  // 3D Radius -- radius in spherical coordinates
        SV_Xv_BD = SV_Radius3D_BD * TMath::Sin(Theta) * TMath::Cos(Phi);                          // x, in cartesian coordinates
        SV_Yv_BD = SV_Radius3D_BD * TMath::Sin(Theta) * TMath::Sin(Phi);                          // y, in cartesian coordinates

        SecondaryVertex_BD.SetXYZ(SV_Xv_BD, SV_Yv_BD, SV_Zv_BD);

        // displacement
        SecondaryVertex_AD = SecondaryVertex_BD + PrimaryVertex;
        SV_Xv_AD = SecondaryVertex_AD.X();
        SV_Yv_AD = SecondaryVertex_AD.Y();
        SV_Zv_AD = SecondaryVertex_AD.Z();
        SV_Radius_AD = SecondaryVertex_AD.Perp();

        tree->Fill();
    }

    tree->Write();

    /* Part 2... suboptimal */
    // (good example: PWGHF/vertexingHF/AliHFPtSpectrum.cxx)
    // (and stad_methods_ws2020_04_MC.pdf, slide 13)

    // prepare gaussian
    TF1* fPtModel = new TF1("fPtModel", "TMath::Gaus(x, 2.5, 0.7)", 0., 5);

    // convert to a histogram
    TH1F* hPtModel = new TH1F("hPtModel", "hPtModel", 100, 0., 5.);
    for (Int_t i = 0; i < 4 * N; i++) {
        hPtModel->Fill(fPtModel->GetRandom());
    }

    TH1F* hPtInjected = new TH1F("hPtInjected", "hPtInjected", 100, 0., 5.);
    for (Int_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        Pt = TMath::Sqrt(Px * Px + Py * Py);
        hPtInjected->Fill(Pt);
    }

    TH1F* hPtWeights = new TH1F("hPtWeights", "hPtWeights", 100, 0., 5.);
    hPtWeights->Divide(hPtModel, hPtInjected, hPtInjected->Integral(), hPtModel->Integral());

    TH1F* hPtReweighted = new TH1F("hPtReweighted", "hPtReweighted", 100, 0., 5.);
    hPtReweighted->Sumw2();
    for (Int_t i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        Pt = TMath::Sqrt(Px * Px + Py * Py);
        hPtReweighted->Fill(Pt, hPtWeights->GetBinContent(hPtWeights->FindBin(Pt)));
    }

    fPtModel->Write();
    hPtModel->Write();
    hPtInjected->Write();
    hPtWeights->Write();
    hPtReweighted->Write();

    /* The end */

    file->Save();
    file->Close();
}
