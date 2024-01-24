#include "TFile.h"
#include "TH1.h"
#include "TList.h"

Int_t Test_HistogramContents(TString filename) {

    /* Open ROOT file */

    TFile* file = TFile::Open(filename);
    if (!file) {
        printf("Error opening file: %s\n", filename.Data());
        return 1;
    }

    /* Get the list of histograms */

    TList* histList = dynamic_cast<TList*>(file->Get("Hists"));
    if (!histList) {
        printf("Error retrieving histogram list\n");
        file->Close();
        return 1;
    }

    /* Loop over the histograms */

    TIter next(histList);
    TObject* obj;

    Int_t n_entries;
    Int_t count_empty = 0;

    printf("%-40s %10s\n", "Name", "N_Entries");
    printf("%-40s %10s\n", "====", "=========");
    while ((obj = next())) {
        TH1* hist = dynamic_cast<TH1*>(obj);
        Int_t n_entries = (Int_t)hist->GetEntries();
        /* Note for myself: */
        /* - `\033[1;31m` is the ANSI escape code for red and bold */
        /* - `\033[0m` is the ANSI escape code to reset color and style to default */
        if (n_entries == 0) {
            printf("%-40s \033[1;31m%10i\033[0m\n", hist->GetName(), n_entries);
            count_empty++;
        } else {
            printf("%-40s %10i\n", hist->GetName(), n_entries);
        }
    }
    printf("%51s\n", "==========");
    printf("%40s %10i\n", "N_Histograms =", histList->GetEntries());
    printf("%40s %10i\n", "N_Empty_Histograms =", count_empty);

    file->Close();

    return 0;
}
