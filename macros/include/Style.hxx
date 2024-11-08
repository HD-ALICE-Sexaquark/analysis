#ifndef Macros_Style_hxx
#define Macros_Style_hxx

#ifndef Macros_Headers_hxx
#include "Headers.hxx"
#endif

extern TColor myBlue;
extern TColor myOrange;
extern TColor myGreen;
extern TColor myRed;
extern TColor myPurple;
extern TColor myBrown;
extern TColor myPink;
extern TColor myGray;
extern TColor myOlive;
extern TColor myCyan;

extern TColor myTransparentBlue;
extern TColor myTransparentOrange;
extern TColor myTransparentGreen;
extern TColor myTransparentRed;
extern TColor myTransparentPurple;
extern TColor myTransparentBrown;
extern TColor myTransparentPink;
extern TColor myTransparentGray;
extern TColor myTransparentOlive;
extern TColor myTransparentCyan;

void DrawHorizontalLine(Double_t x, Color_t color = kBlack, Style_t style = kDashed, Int_t width = 3);
void DrawVerticalLine(Double_t x, Color_t color = kBlack, Style_t style = kDashed, Int_t width = 3);

void SetMyHistStyle(TH1 *theHist);
void SetMyStyle();
void DrawText(TString content, Color_t color = kBlack, Float_t min_x = 0.55, Float_t min_y = 0.8,  //
              Short_t text_h_align = kHAlignCenter, Short_t text_v_align = kVAlignCenter);
void SetMy2DHistStyle(TH2F *theHist);
void SetMy2DStyle();

#endif
