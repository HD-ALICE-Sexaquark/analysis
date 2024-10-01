#include "include/Headers.hxx"
#include "include/Utilities.hxx"

/*
 *
 */
void Merger_ChannelH(TString input_dir = "/home/ceres/borquez/some/output/H1.73_15o_latest", TString output_dir = "output") {

    TString FilenamesPattern = "AnalysisResults_2*.root";
    TString DirBaseName = gSystem->BaseName(input_dir);

    /* Prepare output */

    system("mkdir -p " + output_dir);
    TString OutputFilename = output_dir + "/" + DirBaseName + ".root";
    TFile *OutputFile = new TFile(OutputFilename, "RECREATE");

    TList *OutputList_PosKaonPairs_Hists = new TList();
    TList *OutputList_Sexaquarks_Hists = new TList();

    /***                           ***/
    /*** Case 1: PosKaonPair Hists ***/
    /***                           ***/

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

    /***                          ***/
    /*** Case 2: Sexaquarks Hists ***/
    /***                          ***/

    InputListName = "Sexaquarks_Hists";
    InputHistName = "MCGen_All_AntiSexaquark_Radius";
    TString InputHistXAxisTitle = "Radius (cm)";

    HistPath = input_dir + "/" + FilenamesPattern + "/" + InputListName + "/" + InputHistName;

    ListOfHists = new TList();
    AddObjects(ListOfHists, HistPath.Data());

    /* Add up histograms */

    ThisHist = new TH1F("ThisHist", ";;Counts", 1, 0, 1);  // dummy binning, replaced by first hist binning
    ThisHist->SetName(InputHistName.Data());
    ThisHist->GetXaxis()->SetTitle(InputHistXAxisTitle.Data());

    next = ListOfHists->MakeIterator();
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

    /***         ***/
    /*** The end ***/
    /***         ***/

    /* Save lists */

    OutputFile->cd();

    OutputList_PosKaonPairs_Hists->Write("PosKaonPairs_Hists", 1);
    OutputList_Sexaquarks_Hists->Write("Sexaquarks_Hists", 1);

    /* Close file */

    OutputFile->Close();
    std::cout << ".. Output file " << OutputFilename << " has been created ." << std::endl;

    /* Clean memory */

    delete OutputFile;
    delete OutputList_PosKaonPairs_Hists;
    delete OutputList_Sexaquarks_Hists;

    delete ThisHist;
}
