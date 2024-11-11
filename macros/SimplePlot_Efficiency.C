#include "include/Headers.hxx"
#include "include/Style.cxx"

/*
 *
 */
void SimplePlot_Efficiency(TString OutputDir = "gfx") {

    TFile *InputFile = TFile::Open("output_rdf/EfficiencyResults.root", "READ");
    if (!InputFile || InputFile->IsZombie()) {
        std::cout << "!! ERROR !! SimplePlot_Efficiency !! Couldn't open TFile ../output_rdf/EfficiencyResults.root !!" << std::endl;
        return;
    }

    /* Loop over all histograms on the file, draw them onto a canvas and save them with their name + ".png" */

    TIter next(InputFile->GetListOfKeys());
    TKey *key;

    TString HistName;
    TCanvas *c;
    while ((key = (TKey *)next())) {
        TObject *obj = key->ReadObj();

        if (!obj->InheritsFrom("TH1")) continue;

        TH1 *ThisHist = (TH1 *)obj;
        HistName = ThisHist->GetName();

        /* Info */

        if (HistName.Contains("Efficiency"))
            std::cout << "Avg " + HistName << " = " << ThisHist->Integral() / (Double_t)ThisHist->GetNbinsX() << std::endl;

        /* Draw */

        SetMyStyle();
        gStyle->SetOptStat(0);

        c = new TCanvas("c_" + HistName, "c_" + HistName, 1080, 1080);

        ThisHist->SetTitle("");
        ThisHist->SetLineWidth(4);
        ThisHist->SetLineColor(myBlue.GetNumber());
        ThisHist->SetFillStyle(0);

        ThisHist->GetYaxis()->SetTitleSize(0.04);
        ThisHist->GetYaxis()->SetTitleOffset(1.4);
        ThisHist->GetYaxis()->SetMaxDigits(3);
        ThisHist->GetYaxis()->SetTitleFont(62);

        ThisHist->GetXaxis()->SetTitleSize(0.04);
        ThisHist->GetXaxis()->SetTitleOffset(1.2);
        ThisHist->GetXaxis()->SetTitleFont(62);

        if (HistName.Contains("Radius_"))
            ThisHist->Draw("HIST");
        else
            ThisHist->Draw();

        c->Print(OutputDir + "/" + HistName + ".png");

        delete c;

        std::cout << std::endl;  // empty line
    }
}
