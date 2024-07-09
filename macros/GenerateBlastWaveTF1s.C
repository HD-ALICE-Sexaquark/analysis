// ROOT macro that plots the blastwave function

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom3.h"

#include "Math/Functor.h"
#include "Math/Integrator.h"

/*
  Define blast-wave function.
*/
Double_t BlastWaveFcn(Double_t *x, Double_t *par) {

    Double_t pt = x[0];

    Double_t m0 = par[0];
    Double_t T_kin = par[1];
    Double_t beta_s = par[2];
    Double_t R = par[3];
    Double_t np = par[4];

    Double_t mt = TMath::Sqrt(pt * pt + m0 * m0);

    auto BlastWaveIntegrand = [&](Double_t r) {
        Double_t beta = beta_s * TMath::Power(r / R, np);  // transverse velocity
        Double_t rho = TMath::ATanH(beta);                 // transverse rapidity
        return r * TMath::BesselI0(pt * TMath::SinH(rho) / T_kin) * TMath::BesselK1(mt * TMath::CosH(rho) / T_kin);
    };

    ROOT::Math::Functor1D blastWaveIntFunctor(BlastWaveIntegrand);
    ROOT::Math::Integrator integrator(blastWaveIntFunctor);

    return pt * mt * integrator.Integral(0., 1.);
}

void GenerateBlastWaveTF1s() {

    std::vector<Float_t> m0 = {1.73, 1.8, 1.87, 1.94, 2.01};  // mass of the anti-sexaquark (in GeV/c^2)
    Float_t min_pt = 0.;
    Float_t max_pt = 5.;

    Float_t kin_fo_temp = 0.107;  // kinetic freeze-out temperature (in GeV)
    Float_t beta_s = 0.6;         // transverse velocity at the surface
    Float_t radius = 1.;          // radius of the source
    Float_t np = 1.;              // radial flow power

    TString filename = "blastwave.root";
    TFile *file = new TFile(filename, "RECREATE");

    TRandom3 *rnd = new TRandom3(0);

    TF1 *f;

    TH1D *h;
    Int_t n_bins = 100;
    Int_t n_entries = 1E6;

    for (Float_t &inv_mass : m0) {

        /* Function */

        f = new TF1(Form("BlastWaveFcn_%.2f", inv_mass), BlastWaveFcn, min_pt, max_pt, 5);
        f->SetParameter(0, inv_mass);
        f->SetParameter(1, kin_fo_temp);
        f->SetParameter(2, beta_s);
        f->SetParameter(3, radius);
        f->SetParameter(4, np);

        f->Write();

        /* Histogram */

        h = new TH1D(Form("BlastWaveHist_%.2f", inv_mass), "", 100, min_pt, max_pt);
        h->Sumw2();

        for (Int_t i = 0; i < n_entries; i++) {
            h->Fill(f->GetRandom(rnd));
        }

        h->Write();

        delete h;
        delete f;
    }

    file->Close();

    std::cout << "File " << filename << " has been created." << std::endl;

    delete file;
    delete rnd;
}
