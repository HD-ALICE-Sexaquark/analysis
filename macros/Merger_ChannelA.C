#include "include/Headers.hxx"
#include "include/Utilities.hxx"

/*
 *
 */
void Merger_ChannelA(TString input_dir = "/home/ceres/borquez/some/output/A1.8_18qr",  //
                     TString pattern = "AnalysisResults_29759*.root",                  //
                     TString output_dir = "output_new") {

    TString DirBaseName = gSystem->BaseName(input_dir);

    /* Prepare output */

    system("mkdir -p " + output_dir);
    TString OutputFilename = output_dir + "/" + DirBaseName + ".root";
    TFile *OutputFile = new TFile(OutputFilename, "RECREATE");

    TList *OutputList_Trees = new TList();
    TList *OutputList_V0s_Hists = new TList();
    TList *OutputList_Sexaquarks_Hists = new TList();

    /***               ***/
    /*** Case 1: Trees ***/
    /***               ***/

    TString InputListName = "Trees";

    std::vector<TString> InputTreesNames = {"Injected"};

    TString TreePath;

    TList *ListOfTrees = nullptr;
    TTree *ThisTree = nullptr;

    for (TString &tree_name : InputTreesNames) {

        ListOfTrees = new TList();

        TreePath = input_dir + "/" + pattern + "/" + InputListName + "/" + tree_name;
        AddObjects(ListOfTrees, TreePath.Data());

        /* Add up trees */

        ThisTree = new TTree(tree_name, tree_name);

        TIter next(ListOfTrees);
        while (TObject *obj = next()) {
            TTree *CurrentTree = (TTree *)obj;
            if (ThisTree->GetEntries() == 0) {
                ThisTree = CurrentTree->CloneTree(0);
            } else {
                ThisTree->CopyEntries(CurrentTree);
            }
        }

        /* Save tree */

        OutputList_Trees->Add(ThisTree);

        /* Clean memory */

        delete ListOfTrees;
    }

    /***                   ***/
    /*** Case 2: V0s Hists ***/
    /***                   ***/

    InputListName = "V0s_Hists";

    std::vector<TString> InputV0s_V0sStageNames = {"Found"};
    std::vector<TString> InputV0s_V0sSetNames = {"All", "Signal"};
    std::vector<TString> InputV0s_V0sNames = {"AntiLambda", "KaonZeroShort"};
    std::vector<TString> InputV0s_VarNames = {"Mass", "Pt", "DecayRadius", "OpeningAngle", "CPAwrtPV", "DCAwrtPV", "DCAbtwDau"};

    TString ThisHistName;
    TString HistPath;

    TList *ListOfHists = nullptr;
    TH1F *ThisHist = nullptr;

    for (TString &stage_name : InputV0s_V0sStageNames) {
        for (TString &set_name : InputV0s_V0sSetNames) {
            for (TString &v0_name : InputV0s_V0sNames) {
                for (TString &var_name : InputV0s_VarNames) {

                    ThisHistName = stage_name + "_" + set_name + "_" + v0_name + "_" + var_name;
                    ListOfHists = new TList();

                    HistPath = input_dir + "/" + pattern + "/" + InputListName + "/" + ThisHistName;
                    AddObjects(ListOfHists, HistPath.Data());

                    /* Add up histograms */

                    ThisHist = new TH1F("ThisHist", ";;", 1, 0, 1);  // dummy name+binning, replaced by first hist binning
                    ThisHist->SetName(ThisHistName.Data());

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
            }
        }
    }

    /***                          ***/
    /*** Case 3: Sexaquarks Hists ***/
    /***                          ***/

    InputListName = "Sexaquarks_Hists";

    std::vector<TString> InputHistNames = {"MCGen_All_AntiSexaquark_Radius", "AntiSexaquarks_Bookkeep"};

    for (TString &hist_name : InputHistNames) {

        ListOfHists = new TList();

        HistPath = input_dir + "/" + pattern + "/" + InputListName + "/" + hist_name;
        AddObjects(ListOfHists, HistPath.Data());

        /* Add up histograms */

        ThisHist = new TH1F("ThisHist", ";;", 1, 0, 1);  // dummy name+binning, replaced by first hist binning
        ThisHist->SetName(hist_name.Data());

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

    OutputList_Trees->Write("Trees", 1);
    OutputList_V0s_Hists->Write("V0s_Hists", 1);
    OutputList_Sexaquarks_Hists->Write("Sexaquarks_Hists", 1);

    /* Close file */

    OutputFile->Close();
    std::cout << ".. Output file " << OutputFilename << " has been created ." << std::endl;

    /* Clean memory */

    delete OutputFile;
    delete OutputList_V0s_Hists;
    delete OutputList_Sexaquarks_Hists;

    delete ThisHist;
}
