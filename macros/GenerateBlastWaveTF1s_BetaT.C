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
    Double_t beta_t = par[2];
    Double_t np = par[3];
    Double_t R = par[4];

    Double_t mt = TMath::Sqrt(pt * pt + m0 * m0);

    auto BlastWaveIntegrand = [&](Double_t r) {
        Double_t rho = TMath::ATanH(beta_t);  // transverse rapidity
        return r * TMath::BesselI0(pt * TMath::SinH(rho) / T_kin) * TMath::BesselK1(mt * TMath::CosH(rho) / T_kin);
    };

    ROOT::Math::Functor1D blastWaveIntFunctor(BlastWaveIntegrand);
    ROOT::Math::Integrator integrator(blastWaveIntFunctor);

    return pt * mt * integrator.Integral(0., 1.);
}

void GenerateBlastWaveTF1s_BetaT() {

    std::vector<Float_t> m0 = {1.73, 1.8, 1.87, 1.94, 2.01};  // mass of the anti-sexaquark (in GeV/c^2)
    Float_t min_pt = 0.;
    Float_t max_pt = 5.;

    /* Parameters taken from Production of charged pions, kaons, and (anti-)protons in Pb-Pb and inelastic pp collisions at √sNN = 5.02 TeV */
    /* -- ALICE Collab. PRC 101, 044907 (2020) */
    std::vector<TString> centr_str = {"0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90"};
    std::vector<Float_t> t_kin = {0.090, 0.091, 0.094, 0.097, 0.101, 0.108, 0.115, 0.129, 0.147, 0.161};   // kinetic freeze-out temperature (in GeV)
    std::vector<Float_t> beta_t = {0.663, 0.660, 0.655, 0.643, 0.622, 0.595, 0.557, 0.506, 0.435, 0.355};  // transverse velocity profile
    std::vector<Float_t> np = {0.735, 0.736, 0.739, 0.771, 0.828, 0.908, 1.052, 1.262, 1.678, 2.423};      // radial flow power
    Float_t radius = 1.;                                                                                   // radius of the source

    TString filename = "blastwave_bt.root";
    TFile *file = new TFile(filename, "RECREATE");

    TRandom3 *rnd = new TRandom3(0);

    TF1 *f;

    TH1D *h;
    Int_t n_bins = 100;
    Int_t n_entries = 1E6;

    for (Float_t &inv_mass : m0) {
        for (Int_t centr_i = 0; centr_i < (Int_t)centr_str.size(); centr_i++) {

            /* Function */

            f = new TF1(Form("BlastWaveFcn_%.2f_%s", inv_mass, centr_str[centr_i].Data()), BlastWaveFcn, min_pt, max_pt, 5);
            f->SetParameter(0, inv_mass);
            f->SetParameter(1, t_kin[centr_i]);
            f->SetParameter(2, beta_t[centr_i]);
            f->SetParameter(3, np[centr_i]);
            f->SetParameter(4, radius);

            f->Write();

            /* Histogram */

            h = new TH1D(Form("BlastWaveHist_%.2f_%s", inv_mass, centr_str[centr_i].Data()), "", 100, min_pt, max_pt);
            h->Sumw2();

            for (Int_t i = 0; i < n_entries; i++) {
                h->Fill(f->GetRandom(rnd));
            }

            h->Write();

            delete h;
            delete f;
        }
    }

    file->Close();

    std::cout << "File " << filename << " has been created." << std::endl;

    delete file;
    delete rnd;
}
