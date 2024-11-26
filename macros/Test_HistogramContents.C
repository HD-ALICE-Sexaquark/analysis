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

    TKey* single_key;
    TString list_name;
    TList* list_of_hists;

    UInt_t n_entries;
    UInt_t count_empty;

    /* Loop over file's objects */

    TList* file_content = file->GetListOfKeys();

    for (TObject* single_obj1 : *(file_content)) {

        /* Select only `TList` */

        single_key = dynamic_cast<TKey*>(single_obj1);
        if (static_cast<TString>(single_key->GetClassName()) != "TList") continue;

        /* Get `TList` */

        list_name = single_key->GetName();
        list_of_hists = file->Get<TList>(list_name);

        /* Loop over histograms */

        printf("!! Test_HistogramContents !! # %s !!\n\n", list_name.Data());
        printf("!! Test_HistogramContents !! %-45s %10s %10s %10s !!\n", "Name", "N_Entries", "Mean", "StdDev");
        printf("!! Test_HistogramContents !! %-45s %10s %10s %10s !!\n", "====", "=========", "====", "======");

        count_empty = 0;

        for (TObject* single_obj2 : *(list_of_hists)) {

            TH1* single_hist = dynamic_cast<TH1*>(single_obj2);
            n_entries = (UInt_t)single_hist->GetEntries();

            /* Note: */
            /* - `\033[1;31m` is the ANSI escape code for red and bold */
            /* - `\033[0m` is the ANSI escape code to reset color and style to default */

            if (!n_entries) {
                printf("!! Test_HistogramContents !! %-45s \033[1;31m%10i\033[0m %10s %10s !!\n", single_hist->GetName(), n_entries, "", "");
                count_empty++;
            } else {
                printf("!! Test_HistogramContents !! %-45s %10i %10.3f %10.3f !!\n", single_hist->GetName(), n_entries, single_hist->GetMean(),
                       single_hist->GetStdDev());
            }
        }
        printf("!! Test_HistogramContents !! %56s !!\n", "==========");
        printf("!! Test_HistogramContents !! %45s %10i !!\n", "N_Histograms =", list_of_hists->GetEntries());
        printf("!! Test_HistogramContents !! %45s %10i !!\n\n", "N_Empty_Histograms =", count_empty);
    }  // end of loop over file's objects
    file->Close();
}
