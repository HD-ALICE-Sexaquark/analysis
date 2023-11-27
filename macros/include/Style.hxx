#ifndef STYLE_HXX
#define STYLE_HXX

Color_t myRed = kRed + 1;
Color_t myGreen = kSpring - 8;
Color_t myBlue = kAzure + 2;
Color_t myOrange = kOrange + 2;
Color_t myViolet = kViolet + 2;
Color_t myYellow = kYellow + 1;
Color_t myCyan = kCyan + 1;
Color_t myMagenta = kPink - 8;
Color_t myBlack = kGray + 3;

//________________________________________________________________________
void DrawVerticalLine(Double_t x, Color_t color = kBlack, Style_t style = kDashed, Int_t width = 3) {
    gPad->Update();  // necessary
    Double_t u = (x - gPad->GetX1()) / (gPad->GetX2() - gPad->GetX1());
    // u = (x - x1)/(x2 - x1);
    TLine *linex = new TLine(u, 0.15, u, 0.75);
    linex->SetLineColor(color);  // kBlack, kBlue, kRed, kGreen+2, kOrange+7, kGray+{1,2,3}
    linex->SetLineStyle(style);  // kSolid, kDashed, kDotted
    linex->SetLineWidth(width);
    linex->SetNDC(kTRUE);
    linex->Draw();
}

//________________________________________________________________________
void SetMyHistStyle(TH1F *theHist) {
    theHist->SetTitle("");
    theHist->SetLineWidth(2);
    theHist->SetFillStyle(0);
    theHist->GetYaxis()->SetTitleSize(0.04);
    theHist->GetYaxis()->SetTitleOffset(1.2);
    theHist->GetYaxis()->SetMaxDigits(3);
    theHist->GetXaxis()->SetTitleSize(0.04);
    theHist->GetXaxis()->SetTitleOffset(1.2);
}

//________________________________________________________________________
void SetMyStyle() {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(2);
    // set margin sizes
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
}

//________________________________________________________________________
void DrawText(TString content, Color_t color = kBlack, Float_t min_x = 0.55, Float_t min_y = 0.8,  //
              Short_t text_h_align = kHAlignCenter, Short_t text_v_align = kVAlignCenter) {
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

//________________________________________________________________________
void SetMy2DHistStyle(TH2F *theHist) {
    theHist->SetTitle("");
    theHist->GetYaxis()->SetTitleSize(0.04);
    theHist->GetYaxis()->SetTitleOffset(1.2);
    theHist->GetXaxis()->SetTitleSize(0.04);
    theHist->GetXaxis()->SetTitleOffset(1.2);
    theHist->GetZaxis()->SetMaxDigits(3);
    theHist->GetZaxis()->SetTickLength(0.02);
}

//________________________________________________________________________
void SetMy2DStyle() {
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLineWidth(2);
    // only modify bottom and top margins
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetOptStat(0);
    // set palette color
    gStyle->SetPalette(kDarkBodyRadiator);
    gStyle->SetNumberContours(255);
}

#endif
