#include "include/Headers.hxx"

/*
 *
 */
void V0s_Performance(TString path) {

    /* Find and open .root files ending with customV0s, kalmanV0s, offlineV0s and on-the-flyV0s on path */

    TFile *file_customV0s;
    TFile *file_kalmanV0s;
    TFile *file_offlineV0s;
    TFile *file_ontheflyV0s;

    TSystemDirectory dir(path, path);
    TList *files = dir.GetListOfFiles();
    files->Print();

    if (files) {
        TSystemFile *file;
        TString fname;

        TIter next(files);
        while ((file = (TSystemFile *)next())) {
            fname = file->GetName();

            if (!file->IsDirectory() && fname.EndsWith(".root")) {

                if (fname.Contains("customV0s")) {
                    file_customV0s = TFile::Open(path + fname);
                    if (!file_customV0s) printf("!! ERROR !! V0s_Performance !! Couldn't open TFile %s !!\n", fname.Data());
                }

                if (fname.Contains("kalmanV0s")) {
                    file_kalmanV0s = TFile::Open(path + fname);
                    if (!file_kalmanV0s) printf("!! ERROR !! V0s_Performance !! Couldn't open TFile %s !!\n", fname.Data());
                }

                if (fname.Contains("offlineV0s")) {
                    file_offlineV0s = TFile::Open(path + fname);
                    if (!file_offlineV0s) printf("!! ERROR !! V0s_Performance !! Couldn't open TFile %s !!\n", fname.Data());
                }

                if (fname.Contains("on-the-flyV0s")) {
                    file_ontheflyV0s = TFile::Open(path + fname);
                    if (!file_ontheflyV0s) printf("!! ERROR !! V0s_Performance !! Couldn't open TFile %s !!\n", fname.Data());
                }
            }
        }
    }

    /* Load the TH1 "Findable_All_AntiLambda_mass" within the TList "Hists" on every file */

    TList *lCustomV0s = (TList *)file_customV0s->Get("Hists");
    TList *lKalmanV0s = (TList *)file_kalmanV0s->Get("Hists");
    TList *lOfflineV0s = (TList *)file_offlineV0s->Get("Hists");
    TList *lOntheflyV0s = (TList *)file_ontheflyV0s->Get("Hists");

    TH1F *hFindable_AS_CustomV0s = (TH1F *)lCustomV0s->FindObject("Findable_Signal_AntiLambda_Radius");
    TH1F *hFindable_AS_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Findable_Signal_AntiLambda_Radius");
    TH1F *hFindable_AS_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Findable_Signal_AntiLambda_Radius");
    TH1F *hFindable_AS_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Findable_Signal_AntiLambda_Radius");

    TH1F *hFound_All_AS_CustomV0s = (TH1F *)lCustomV0s->FindObject("Found_All_AntiLambda_Mass");
    TH1F *hFound_All_AS_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Found_All_AntiLambda_Mass");
    TH1F *hFound_All_AS_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Found_All_AntiLambda_Mass");
    TH1F *hFound_All_AS_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Found_All_AntiLambda_Mass");

    TH1F *hFound_Signal_AS_CustomV0s = (TH1F *)lCustomV0s->FindObject("Found_Signal_AntiLambda_Mass");
    TH1F *hFound_Signal_AS_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Found_Signal_AntiLambda_Mass");
    TH1F *hFound_Signal_AS_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Found_Signal_AntiLambda_Mass");
    TH1F *hFound_Signal_AS_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Found_Signal_AntiLambda_Mass");

    /*  */

    TH1F *hFindable_K0S_CustomV0s = (TH1F *)lCustomV0s->FindObject("Findable_Signal_KaonZeroShort_Radius");
    TH1F *hFindable_K0S_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Findable_Signal_KaonZeroShort_Radius");
    TH1F *hFindable_K0S_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Findable_Signal_KaonZeroShort_Radius");
    TH1F *hFindable_K0S_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Findable_Signal_KaonZeroShort_Radius");

    TH1F *hFound_All_K0S_CustomV0s = (TH1F *)lCustomV0s->FindObject("Found_All_KaonZeroShort_Mass");
    TH1F *hFound_All_K0S_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Found_All_KaonZeroShort_Mass");
    TH1F *hFound_All_K0S_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Found_All_KaonZeroShort_Mass");
    TH1F *hFound_All_K0S_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Found_All_KaonZeroShort_Mass");

    TH1F *hFound_Signal_K0S_CustomV0s = (TH1F *)lCustomV0s->FindObject("Found_Signal_KaonZeroShort_Mass");
    TH1F *hFound_Signal_K0S_KalmanV0s = (TH1F *)lKalmanV0s->FindObject("Found_Signal_KaonZeroShort_Mass");
    TH1F *hFound_Signal_K0S_OfflineV0s = (TH1F *)lOfflineV0s->FindObject("Found_Signal_KaonZeroShort_Mass");
    TH1F *hFound_Signal_K0S_OntheflyV0s = (TH1F *)lOntheflyV0s->FindObject("Found_Signal_KaonZeroShort_Mass");

    /* Print their number of entries */

    TH1F *hNFound_Signal_AS = new TH1F("", "", 4, 0, 4);
    hNFound_Signal_AS->SetBinContent(1, hFound_Signal_AS_OfflineV0s->GetEntries());
    hNFound_Signal_AS->SetBinError(1, TMath::Sqrt(hFound_Signal_AS_OfflineV0s->GetEntries()));
    hNFound_Signal_AS->SetBinContent(2, hFound_Signal_AS_OntheflyV0s->GetEntries());
    hNFound_Signal_AS->SetBinError(2, TMath::Sqrt(hFound_Signal_AS_OntheflyV0s->GetEntries()));
    hNFound_Signal_AS->SetBinContent(3, hFound_Signal_AS_CustomV0s->GetEntries());
    hNFound_Signal_AS->SetBinError(3, TMath::Sqrt(hFound_Signal_AS_CustomV0s->GetEntries()));
    hNFound_Signal_AS->SetBinContent(4, hFound_Signal_AS_KalmanV0s->GetEntries());
    hNFound_Signal_AS->SetBinError(4, TMath::Sqrt(hFound_Signal_AS_KalmanV0s->GetEntries()));

    TH1F *hNFound_All_AS = new TH1F("", "", 4, 0, 4);
    hNFound_All_AS->SetBinContent(1, hFound_All_AS_OfflineV0s->GetEntries());
    hNFound_All_AS->SetBinError(1, TMath::Sqrt(hFound_All_AS_OfflineV0s->GetEntries()));
    hNFound_All_AS->SetBinContent(2, hFound_All_AS_OntheflyV0s->GetEntries());
    hNFound_All_AS->SetBinError(2, TMath::Sqrt(hFound_All_AS_OntheflyV0s->GetEntries()));
    hNFound_All_AS->SetBinContent(3, hFound_All_AS_CustomV0s->GetEntries());
    hNFound_All_AS->SetBinError(3, TMath::Sqrt(hFound_All_AS_CustomV0s->GetEntries()));
    hNFound_All_AS->SetBinContent(4, hFound_All_AS_KalmanV0s->GetEntries());
    hNFound_All_AS->SetBinError(4, TMath::Sqrt(hFound_All_AS_KalmanV0s->GetEntries()));

    TH1F *hNFindable_Signal_AS = new TH1F("", "", 4, 0, 4);
    hNFindable_Signal_AS->SetBinContent(1, hFindable_AS_OfflineV0s->GetEntries());
    hNFindable_Signal_AS->SetBinError(1, TMath::Sqrt(hFindable_AS_OfflineV0s->GetEntries()));
    hNFindable_Signal_AS->SetBinContent(2, hFindable_AS_OntheflyV0s->GetEntries());
    hNFindable_Signal_AS->SetBinError(2, TMath::Sqrt(hFindable_AS_OntheflyV0s->GetEntries()));
    hNFindable_Signal_AS->SetBinContent(3, hFindable_AS_CustomV0s->GetEntries());
    hNFindable_Signal_AS->SetBinError(3, TMath::Sqrt(hFindable_AS_CustomV0s->GetEntries()));
    hNFindable_Signal_AS->SetBinContent(4, hFindable_AS_KalmanV0s->GetEntries());
    hNFindable_Signal_AS->SetBinError(4, TMath::Sqrt(hFindable_AS_KalmanV0s->GetEntries()));

    /* print hNFindable_Signal_AS bin contents */

    for (Int_t i = 1; i <= 4; i++) {
        printf("!! V0s_Performance !! hNFound_Signal_AS bin %d: %f !!\n", i, hNFound_Signal_AS->GetBinContent(i));
        printf("!! V0s_Performance !! hNFindable_Signal_AS bin %d: %f !!\n", i, hNFindable_Signal_AS->GetBinContent(i));
    }

    TH1F *hPurity = new TH1F("", "", 4, 0, 4);  // N(signal) / N(all)
    hPurity->Divide(hNFound_Signal_AS, hNFound_All_AS, 1, 1, "B");

    TH1F *hEfficiency = new TH1F("", "", 4, 0, 4);  // N(signal) / N(findable signal)
    hEfficiency->Divide(hNFound_Signal_AS, hNFindable_Signal_AS, 1, 1, "B");

    /*  */

    TCanvas *cPurity = new TCanvas("c", "c", 2560, 2560);

    hPurity->Draw("E");
    hPurity->SetStats(0);

    hPurity->GetXaxis()->SetBinLabel(1, "OfflineV0s");
    hPurity->GetXaxis()->SetBinLabel(2, "On-the-flyV0s");
    hPurity->GetXaxis()->SetBinLabel(3, "CustomV0s");
    hPurity->GetXaxis()->SetBinLabel(4, "KalmanV0s");

    hPurity->GetYaxis()->SetTitle("Number of entries");
    hPurity->GetYaxis()->SetTitleOffset(1.5);
    hPurity->GetYaxis()->CenterTitle();

    cPurity->SaveAs("V0s_Purity.png");

    /* print hPurity bin contents and errors */

    for (Int_t i = 1; i <= 4; i++) {
        printf("!! V0s_Performance !! hPurity bin %d: %f +- %f !!\n", i, hPurity->GetBinContent(i), hPurity->GetBinError(i));
    }

    /*  */

    TCanvas *cEfficiency = new TCanvas("c2", "c2", 2560, 2560);

    hEfficiency->Draw("E");
    hEfficiency->SetStats(0);

    hEfficiency->GetXaxis()->SetBinLabel(1, "OfflineV0s");
    hEfficiency->GetXaxis()->SetBinLabel(2, "On-the-flyV0s");
    hEfficiency->GetXaxis()->SetBinLabel(3, "CustomV0s");
    hEfficiency->GetXaxis()->SetBinLabel(4, "KalmanV0s");

    hEfficiency->GetYaxis()->SetTitle("Number of entries");
    hEfficiency->GetYaxis()->SetTitleOffset(1.5);
    hEfficiency->GetYaxis()->CenterTitle();

    cEfficiency->SaveAs("V0s_Efficiency.png");

    /* print hEfficiency bin contents and errors */

    for (Int_t i = 1; i <= 4; i++) {
        printf("!! V0s_Performance !! hEfficiency bin %d: %f +- %f !!\n", i, hEfficiency->GetBinContent(i), hEfficiency->GetBinError(i));
    }

    /* Close files */

    if (file_customV0s) file_customV0s->Close();
    if (file_kalmanV0s) file_kalmanV0s->Close();
    if (file_offlineV0s) file_offlineV0s->Close();
    if (file_ontheflyV0s) file_ontheflyV0s->Close();
}
