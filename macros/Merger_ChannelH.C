#include "include/Headers.hxx"
#include "include/Utilities.hxx"

/*
 *
 */
void Merger_ChannelH(TString input_dir = "/home/ceres/borquez/some/output/H1.73_15o_latest", TString output_dir = "output/") {

    TString FilenamesPattern = input_dir.Contains("15o") ? "AnalysisResults_246991*.root" : "AnalysisResults_29759*.root";
    TString DirBaseName = gSystem->BaseName(input_dir);

    /* Prepare output */

    system("mkdir -p " + output_dir);
    TFile *OutputFile = new TFile(output_dir + DirBaseName + ".root", "RECREATE");

    TList *OutputList_PosKaonPairs_Hists = new TList();

    /* Begin */

    TString InputListName = "PosKaonPairs_Hists";
    TString InputHistName = "PosKaonPairs_Bookkeep";

    TString HistPath = input_dir + "/" + FilenamesPattern + "/" + InputListName + "/" + InputHistName;

    TList *ListOfHists = new TList();
    AddObjects(ListOfHists, HistPath.Data());

    /* Add up histograms */

    TH1F *ThisHist = new TH1F("ThisHist", ";;Counts", 1, 0, 1);  // dummy binning, replaced by first hist binning
    ThisHist->SetName(InputHistName.Data());

    TIter next(ListOfHists);
    while (TObject *obj = next()) {
        TH1F *CurrentHist = (TH1F *)obj;
        if (!ThisHist->GetEntries()) {
            ThisHist->SetBins(CurrentHist->GetNbinsX(), CurrentHist->GetXaxis()->GetXmin(), CurrentHist->GetXaxis()->GetXmax());
        }
        ThisHist->Add(CurrentHist);
    }

    /* Save histogram */

    OutputList_PosKaonPairs_Hists->Add(ThisHist);

    /* Clean memory */

    delete ListOfHists;

    /* Save lists */

    OutputFile->cd();

    OutputList_PosKaonPairs_Hists->Write(InputListName, 1);

    /* Close file */

    OutputFile->Close();
    std::cout << ".. Output file " << output_dir + DirBaseName + ".root"
              << " has been created ." << std::endl;

    /* Clean memory */

    delete OutputFile;
    delete OutputList_PosKaonPairs_Hists;

    delete ThisHist;
}
