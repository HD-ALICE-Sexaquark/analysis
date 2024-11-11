#include "include/Headers.hxx"
#include "include/Style.cxx"

/*
 *
 */
void QA_OneDimHists(TString input_dir = "output", TString output_dir = "gfx") {

    std::vector<TString> InputSignalFilenames;
    std::vector<TString> InputSignalFilenamesA = {"QA_ChannelA_18qr.root", "QA_ChannelD_18qr.root", "QA_ChannelE_18qr.root", "QA_ChannelH_18qr.root"};
    std::vector<TString> InputSignalFilenamesB = {"QA_ChannelA_15o.root", "QA_ChannelD_15o.root", "QA_ChannelE_15o.root", "QA_ChannelH_15o.root"};
    TString InputBkgFilename;
    TString SignalProdName;
    TString BkgProdName;

    TString InputListName = "QA_Hists";

    std::vector<TString> InputHistNames = {"QA_MCGen_VertexZ", "QA_VertexZ",      "QA_SPD_NTracklets", "QA_NTracks",
                                           "QA_Tracks_Pt",     "QA_Tracks_DCAxy", "QA_Tracks_DCAz"};

    TFile *InputFile;
    TList *ListOfHists;
    TH1F *ThisHist;

    std::vector<TColor> TranspColors = {myTransparentBlue, myTransparentOrange, myTransparentGreen, myTransparentRed};
    std::vector<TColor> Colors = {myBlue, myOrange, myGreen, myRed};

    TString HistPath;
    TCanvas *c;
    TLegend *leg;

    for (Int_t hh = 0; hh < (Int_t)InputHistNames.size(); hh++) {

        c = new TCanvas("c", "c", 2 * 800, 800);
        c->Divide(2, 1);

        SetMyStyle();
        gStyle->SetOptStat(0);

        for (Int_t pp = 0; pp < 2; pp++) {

            leg = new TLegend(0.65, 0.72, 0.8, 0.87);
            leg->SetTextFont(62);
            leg->SetTextSize(0.03);
            leg->SetBorderSize(0);

            if (!pp) {
                InputSignalFilenames = InputSignalFilenamesA;
                InputBkgFilename = "QA_bkg_18qr_latest.root";
                SignalProdName = "23l1a3";
                BkgProdName = "20e3a";
                leg->SetHeader("MC anchored to 18qr");
            } else {
                InputSignalFilenames = InputSignalFilenamesB;
                InputBkgFilename = "QA_bkg_15o_latest.root";
                SignalProdName = "23l1b3";
                BkgProdName = "20j6a";
                leg->SetHeader("MC anchored to 15o");
            }

            /* Bkg */

            InputFile = TFile::Open(input_dir + "/" + InputBkgFilename, "READ");
            ListOfHists = (TList *)InputFile->Get(InputListName);
            ThisHist = (TH1F *)ListOfHists->FindObject(InputHistNames[hh]);

            c->cd(pp + 1);

            SetMyHistStyle(ThisHist);
            ThisHist->SetLineWidth(3);
            ThisHist->SetLineColor(myPink.GetNumber());

            ThisHist->GetYaxis()->SetTitle("Normalized Counts");
            ThisHist->GetXaxis()->SetLabelSize(0.03);
            ThisHist->Scale(1.0 / ThisHist->Integral());

            ThisHist->Draw("SAME");

            leg->AddEntry(ThisHist, BkgProdName, "l");
            ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(myPink.GetNumber());

            /* Signal */

            for (Int_t ff = 0; ff < (Int_t)InputSignalFilenames.size(); ff++) {

                InputFile = TFile::Open(input_dir + "/" + InputSignalFilenames[ff], "READ");
                ListOfHists = (TList *)InputFile->Get(InputListName);
                ThisHist = (TH1F *)ListOfHists->FindObject(InputHistNames[hh]);

                c->cd(pp + 1);

                SetMyHistStyle(ThisHist);
                ThisHist->SetLineWidth(3);
                ThisHist->SetLineColor(TranspColors[ff].GetNumber());

                ThisHist->GetYaxis()->SetTitle("Normalized Counts");
                ThisHist->GetXaxis()->SetLabelSize(0.03);
                if (InputHistNames[hh].Contains("QA_NTracks"))
                    ThisHist->Scale(1.0 / ThisHist->GetMaximum());
                else
                    ThisHist->Scale(1.0 / ThisHist->Integral());

                ThisHist->Draw("SAME");

                TString channel = InputSignalFilenames[ff](10, 1);
                leg->AddEntry(ThisHist, SignalProdName + " (ch." + channel + ")", "l");
                ((TLegendEntry *)leg->GetListOfPrimitives()->Last())->SetTextColor(Colors[ff].GetNumber());
            }

            leg->Draw();
            if (InputHistNames[hh] != "QA_NTracks") gPad->SetLogy();
        }

        system("mkdir -p " + output_dir);
        c->Print(output_dir + "/" + InputHistNames[hh] + ".png");

        delete leg;
        delete c;
    }
}
