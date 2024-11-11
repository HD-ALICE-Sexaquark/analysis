#include "include/Headers.hxx"

/*
 *
 */
void Test_HistogramContents(TString filename) {

    /* Open ROOT file */

    TFile* file = TFile::Open(filename);
    if (!file) {
        printf("!! ERROR !! Test_HistogramContents !! Couldn't open TFile %s !!\n", filename.Data());
        return;
    }

    /* Get lists of histograms */

    std::vector<TString> hist_containers = {"QA_Hists", "Tracks_Hists", "V0s_Hists", "Sexaquarks_Hists", "PosKaonPairs_Hists"};

    for (TString& hist_ctn : hist_containers) {

        TList* histList = dynamic_cast<TList*>(file->Get(hist_ctn));
        if (!histList) {
            printf("!! ERROR !! Test_HistogramContents !! Couldn't find TList %s !!\n", hist_ctn.Data());
            file->Close();
            return;
        }

        /* Loop over histograms */

        TIter next(histList);
        TObject* obj;

        Int_t n_entries;
        Int_t count_empty = 0;

        printf("!! Test_HistogramContents !! # %s !!\n\n", hist_ctn.Data());
        printf("!! Test_HistogramContents !! %-45s %10s !!\n", "Name", "N_Entries");
        printf("!! Test_HistogramContents !! %-45s %10s !!\n", "====", "=========");
        while ((obj = next())) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            Int_t n_entries = (Int_t)hist->GetEntries();
            /* Note: */
            /* - `\033[1;31m` is the ANSI escape code for red and bold */
            /* - `\033[0m` is the ANSI escape code to reset color and style to default */
            if (!n_entries) {
                printf("!! Test_HistogramContents !! %-45s \033[1;31m%10i\033[0m !!\n", hist->GetName(), n_entries);
                count_empty++;
            } else {
                printf("!! Test_HistogramContents !! %-45s %10i !!\n", hist->GetName(), n_entries);
            }
        }
        printf("!! Test_HistogramContents !! %56s !!\n", "==========");
        printf("!! Test_HistogramContents !! %45s %10i !!\n", "N_Histograms =", histList->GetEntries());
        printf("!! Test_HistogramContents !! %45s %10i !!\n\n", "N_Empty_Histograms =", count_empty);
    }

    file->Close();
}
