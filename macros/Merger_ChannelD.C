#include "include/Headers.hxx"
#include "include/Utilities.hxx"

/*
 *
 */
void Merger_ChannelD(TString input_dir = "/home/ceres/borquez/some/output/D1.73_15o_latest", TString output_dir = "output/") {

    TString FilenamesPattern = input_dir.Contains("15o") ? "AnalysisResults_246991*.root" : "AnalysisResults_29759*.root";
    TString DirBaseName = gSystem->BaseName(input_dir);

    /* Prepare output */

    system("mkdir -p " + output_dir);
    TFile *OutputFile = new TFile(output_dir + DirBaseName + ".root", "RECREATE");

    TList *OutputList_V0s_Hists = new TList();
    TList *OutputList_Sexaquarks_Hists = new TList();

    /***                  ***/
    /*** Case 1: V0 Hists ***/
    /***                  ***/

    TString InputListName = "V0s_Hists";

    std::vector<TString> InputHistNames = {"Found_All_AntiLambda_Mass", "Found_Signal_AntiLambda_Mass"};
    std::vector<TString> InputHistXAxisTitle = {"m(#bar{p} #pi^{+}) (GeV/c^{2})", "m(#bar{p} #pi^{+}) (GeV/c^{2})"};

    TString ThisHistName;
    TString HistPath;

    TList *ListOfHists;
    TH1F *ThisHist;

    for (Int_t i = 0; i < (Int_t)InputHistNames.size(); i++) {

        ThisHistName = InputHistNames[i];
        ListOfHists = new TList();

        HistPath = input_dir + "/" + FilenamesPattern + "/" + InputListName + "/" + ThisHistName;
        AddObjects(ListOfHists, HistPath.Data());

        /* Add up histograms */

        ThisHist = new TH1F("ThisHist", ";;Counts", 1, 0, 1);  // dummy binning, replaced by first hist binning
        ThisHist->SetName(ThisHistName.Data());
        ThisHist->GetXaxis()->SetTitle(InputHistXAxisTitle[i].Data());

        TIter next(ListOfHists);
        while (TObject *obj = next()) {
            TH1F *CurrentHist = (TH1F *)obj;
            if (!ThisHist->GetEntries()) {
                ThisHist->SetBins(CurrentHist->GetNbinsX(), CurrentHist->GetXaxis()->GetXmin(), CurrentHist->GetXaxis()->GetXmax());
            }
            ThisHist->Add(CurrentHist);
        }

        /* Save histogram */

        OutputList_V0s_Hists->Add(ThisHist);

        /* Clean memory */

        delete ListOfHists;
    }

    /***                          ***/
    /*** Case 2: Sexaquarks Hists ***/
    /***                          ***/

    InputListName = "Sexaquarks_Hists";

    InputHistNames = {"AntiSexaquarks_Bookkeep"};
    InputHistXAxisTitle = {""};

    for (Int_t i = 0; i < (Int_t)InputHistNames.size(); i++) {

        ThisHistName = InputHistNames[i];
        ListOfHists = new TList();

        HistPath = input_dir + "/" + FilenamesPattern + "/" + InputListName + "/" + ThisHistName;
        AddObjects(ListOfHists, HistPath.Data());

        /* Add up histograms */

        ThisHist = new TH1F("ThisHist", ";;Counts", 1, 0, 1);  // dummy binning, replaced by first hist binning
        ThisHist->SetName(ThisHistName.Data());
        ThisHist->GetXaxis()->SetTitle(InputHistXAxisTitle[i].Data());

        TIter next(ListOfHists);
        while (TObject *obj = next()) {
            TH1F *CurrentHist = (TH1F *)obj;
            if (!ThisHist->GetEntries()) {
                ThisHist->SetBins(CurrentHist->GetNbinsX(), CurrentHist->GetXaxis()->GetXmin(), CurrentHist->GetXaxis()->GetXmax());
            }
            ThisHist->Add(CurrentHist);
        }

        /* Save histogram */

        OutputList_Sexaquarks_Hists->Add(ThisHist);

        /* Clean memory */

        delete ListOfHists;
    }

    /***         ***/
    /*** The end ***/
    /***         ***/

    /* Save lists */

    OutputFile->cd();

    OutputList_V0s_Hists->Write("V0s_Hists", 1);
    OutputList_Sexaquarks_Hists->Write("Sexaquarks_Hists", 1);

    /* Close file */

    OutputFile->Close();
    std::cout << ".. Output file " << output_dir + DirBaseName + ".root"
              << " has been created ." << std::endl;

    /* Clean memory */

    delete OutputFile;
    delete OutputList_V0s_Hists;
    delete OutputList_Sexaquarks_Hists;

    delete ThisHist;
}
