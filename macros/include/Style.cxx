#include "Style.hxx"

// colors based on matplotlib main colors
TColor myBlue(TColor::GetFreeColorIndex(), 31. / 255., 119. / 255., 180. / 255., "myBlue");
TColor myOrange(TColor::GetFreeColorIndex(), 255. / 255., 127. / 255., 14. / 255., "myOrange");
TColor myGreen(TColor::GetFreeColorIndex(), 44. / 255., 160. / 255., 44. / 255., "myGreen");
TColor myRed(TColor::GetFreeColorIndex(), 214. / 255., 39. / 255., 40. / 255., "myRed");
TColor myPurple(TColor::GetFreeColorIndex(), 148. / 255., 103. / 255., 189. / 255., "myPurple");
TColor myBrown(TColor::GetFreeColorIndex(), 140. / 255., 86. / 255., 75. / 255., "myBrown");
TColor myPink(TColor::GetFreeColorIndex(), 227. / 255., 119. / 255., 194. / 255., "myPink");
TColor myGray(TColor::GetFreeColorIndex(), 127. / 255., 127. / 255., 127. / 255., "myGray");
TColor myOlive(TColor::GetFreeColorIndex(), 188. / 255., 189. / 255., 34. / 255., "myOlive");
TColor myCyan(TColor::GetFreeColorIndex(), 23. / 255., 190. / 255., 207. / 255., "myCyan");

TColor myTransparentBlue(TColor::GetFreeColorIndex(), 31. / 255., 119. / 255., 180. / 255., "myTransparentBlue", 0.5);
TColor myTransparentOrange(TColor::GetFreeColorIndex(), 255. / 255., 127. / 255., 14. / 255., "myTransparentOrange", 0.5);
TColor myTransparentGreen(TColor::GetFreeColorIndex(), 44. / 255., 160. / 255., 44. / 255., "myTransparentGreen", 0.5);
TColor myTransparentRed(TColor::GetFreeColorIndex(), 214. / 255., 39. / 255., 40. / 255., "myTransparentRed", 0.5);
TColor myTransparentPurple(TColor::GetFreeColorIndex(), 148. / 255., 103. / 255., 189. / 255., "myTransparentPurple", 0.5);
TColor myTransparentBrown(TColor::GetFreeColorIndex(), 140. / 255., 86. / 255., 75. / 255., "myTransparentBrown", 0.5);
TColor myTransparentPink(TColor::GetFreeColorIndex(), 227. / 255., 119. / 255., 194. / 255., "myTransparentPink", 0.5);
TColor myTransparentGray(TColor::GetFreeColorIndex(), 127. / 255., 127. / 255., 127. / 255., "myTransparentGray", 0.5);
TColor myTransparentOlive(TColor::GetFreeColorIndex(), 188. / 255., 189. / 255., 34. / 255., "myTransparentOlive", 0.5);
TColor myTransparentCyan(TColor::GetFreeColorIndex(), 23. / 255., 190. / 255., 207. / 255., "myTransparentCyan", 0.5);

void DrawHorizontalLine(Double_t x, Color_t color, Style_t style, Int_t width) {
    gPad->Update();
    Double_t v = (x - gPad->GetY1()) / (gPad->GetY2() - gPad->GetY1());
    // v = (y - y1)/(y2 - y1);
    TLine *liney = new TLine(gPad->GetLeftMargin(), v,  //
                             1. - gPad->GetRightMargin(), v);
    liney->SetLineColor(color);  // kBlack, kBlue, kRed, kGreen+2, kOrange+7, kGray+{1,2,3}
    liney->SetLineStyle(style);  // kSolid, kDashed, kDotted
    liney->SetLineWidth(width);
    liney->SetNDC(kTRUE);
    liney->Draw();
}

void DrawVerticalLine(Double_t x, Color_t color, Style_t style, Int_t width) {
    gPad->Update();
    Double_t u = (x - gPad->GetX1()) / (gPad->GetX2() - gPad->GetX1());
    // u = (x - x1)/(x2 - x1);
    TLine *linex = new TLine(u, gPad->GetBottomMargin(),  //
                             u, 1. - gPad->GetTopMargin());
    linex->SetLineColor(color);  // kBlack, kBlue, kRed, kGreen+2, kOrange+7, kGray+{1,2,3}
    linex->SetLineStyle(style);  // kSolid, kDashed, kDotted
    linex->SetLineWidth(width);
    linex->SetNDC(kTRUE);
    linex->Draw();
}

void SetMyHistStyle(TH1 *theHist) {
    theHist->SetTitle("");
    theHist->SetLineWidth(2);
    theHist->SetFillStyle(0);

    theHist->GetYaxis()->SetTitleSize(0.04);
    theHist->GetYaxis()->SetTitleOffset(1.4);
    theHist->GetYaxis()->SetMaxDigits(3);
    theHist->GetYaxis()->SetTitleFont(62);

    theHist->GetXaxis()->SetTitleSize(0.04);
    theHist->GetXaxis()->SetTitleOffset(1.2);
    theHist->GetXaxis()->SetTitleFont(62);
}

void SetMyStyle() {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(2);
    // set margin sizes
    gStyle->SetPadTopMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
}

void DrawText(TString content, Color_t color, Float_t min_x, Float_t min_y, Short_t text_h_align, Short_t text_v_align) {
    TPaveText *pav = new TPaveText(min_x, min_y, min_x + 0.3, min_y + 0.1, "NDC");  // x1,y1,x2,y2
    pav->AddText(content);
    pav->SetTextSize(0.03);
    pav->SetTextColor(color);
    ((TText *)pav->GetListOfLines()->Last())->SetTextAlign(text_h_align + text_v_align);
    pav->SetBorderSize(0);
    pav->SetFillStyle(0);
    // pav->SetFillColor(kWhite);
    pav->SetTextAlign(12);
    pav->Draw();
}

void SetMy2DHistStyle(TH2F *theHist) {
    theHist->SetTitle("");

    theHist->GetYaxis()->SetTitleSize(0.04);
    theHist->GetYaxis()->SetTitleOffset(1.2);
    theHist->GetYaxis()->SetTitleFont(62);

    theHist->GetXaxis()->SetTitleSize(0.04);
    theHist->GetXaxis()->SetTitleOffset(1.2);
    theHist->GetXaxis()->SetTitleFont(62);

    theHist->GetZaxis()->SetMaxDigits(3);
    theHist->GetZaxis()->SetTickLength(0.02);
}

void SetMy2DStyle() {
    gStyle->SetLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    // gStyle->SetPalette(kDarkBodyRadiator);
    gStyle->SetNumberContours(255);
}
