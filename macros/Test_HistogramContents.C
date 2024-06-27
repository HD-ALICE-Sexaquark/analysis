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

    std::vector<TString> hist_containers = {"QA_Hists", "Tracks_Hists", "V0s_Hists", "Sexaquarks_Hists", "PosKaonPairs_Hists"};

    for (TString& hist_ctn : hist_containers) {

        TList* histList = dynamic_cast<TList*>(file->Get(hist_ctn));
        if (!histList) {
            printf("Error retrieving histogram list %s\n", hist_ctn.Data());
            file->Close();
            return 1;
        }

        /* Loop over the histograms */

        TIter next(histList);
        TObject* obj;

        Int_t n_entries;
        Int_t count_empty = 0;

        printf("# %s\n\n", hist_ctn.Data());
        printf("%-45s %10s\n", "Name", "N_Entries");
        printf("%-45s %10s\n", "====", "=========");
        while ((obj = next())) {
            TH1* hist = dynamic_cast<TH1*>(obj);
            Int_t n_entries = (Int_t)hist->GetEntries();
            /* Note for myself: */
            /* - `\033[1;31m` is the ANSI escape code for red and bold */
            /* - `\033[0m` is the ANSI escape code to reset color and style to default */
            if (n_entries == 0) {
                printf("%-45s \033[1;31m%10i\033[0m\n", hist->GetName(), n_entries);
                count_empty++;
            } else {
                printf("%-45s %10i\n", hist->GetName(), n_entries);
            }
        }
        printf("%56s\n", "==========");
        printf("%45s %10i\n", "N_Histograms =", histList->GetEntries());
        printf("%45s %10i\n\n", "N_Empty_Histograms =", count_empty);
    }

    file->Close();

    return 0;
}
