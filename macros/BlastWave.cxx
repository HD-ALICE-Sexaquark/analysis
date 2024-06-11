// ROOT macro that plots the blastwave function

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"

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

    // return pt * mt * integrator.Integral(0., 1.);
    return mt * integrator.Integral(0., 1.);
}

void BlastWave() {

    Float_t min_pt = 0.;
    Float_t max_pt = 5.;
    Float_t m0 = 1.73;     // mass of the anti-sexaquark (in GeV/c^2)
    Float_t Tkin = 0.107;  // kinetic freeze-out temperature (in GeV)
    Float_t beta_s = 0.6;  // transverse velocity at the surface
    Float_t R = 1.;        // radius of the source
    Float_t np = 1.;       // radial flow power

    /*
    // debug
        Double_t par[5] = {0.938, 0.1, 0.8, 1., 1.};
    for (Double_t pt = 0.; pt < 5.; pt += 0.5) {
        std::cout << pt << " " << BlastWaveFcn(&pt, par) << std::endl;
    }
    */

    TF1 *f = new TF1("f", BlastWaveFcn, min_pt, max_pt, 5);
    f->SetParameter(0, m0);
    f->SetParameter(1, Tkin);
    f->SetParameter(2, beta_s);
    f->SetParameter(3, R);
    f->SetParameter(4, np);

    f->SetLineColor(kBlue);

    // TF1 *f2 = new TF1("f2", BlastWaveFcn, min_pt, max_pt, 5);
    // f2->SetParameter(0, m0);
    // f2->SetParameter(1, 0.1);
    // f2->SetParameter(2, 0.8);
    // f2->SetParameter(3, R);
    // f2->SetParameter(4, np);

    // f2->SetLineColor(kRed);

    TCanvas *c = new TCanvas("c", "c", 720, 720);

    f->Draw();
    // f2->Draw("SAME");
}
