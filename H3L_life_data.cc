#include <iostream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <fstream>
#include <math.h>
#include <string>
#include <unistd.h>
using namespace std;

#include "TApplication.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"
#include "TCut.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TText.h"
#include <TLegend.h>
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TString.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TMath.h"
#include "TGaxis.h"

#include "TRandom.h"

void SetGrAsymmErr(TGraphAsymmErrors *gr, int LColor, int MColor, int MStyle, double Msize){
//  gr->SetTitle(hname);
//  gr->SetName(hname);
//  gr->GetXaxis()->SetTitle(xname);
//  gr->GetXaxis()->CenterTitle();
//  gr->GetYaxis()->SetTitle(yname);
//  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(Msize);
//  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}

void SetTH2(TH2F *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);

  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->SetTitleOffset(0.85);
//  h2->GetXaxis()->CenterTitle();
  h2->GetXaxis()->SetNdivisions(000);
  h2->GetXaxis()->SetLabelSize(0.05);

  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->SetTitleOffset(0.95);
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)h2->GetYaxis())->SetMaxDigits(3);

  h2->SetMinimum(0.8);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.06,"x");
  h2->SetTitleSize(0.06,"y");
}

////////////////
void ReadCSV(string ifname,  vector<double> &xaxis, vector<double> &yaxis,  vector<double> &er_x, vector<double> &er_y, int &count, double xaxismin,double xaxismax, double xaxisoffset){
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double x,y,ey,ex;
double runnum;

count=0;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
  istringstream sline(line);
//csv data structure expected
/* runnum xvalue yvalue yerror  */
//for each line
  sline >> runnum;
  sline >> x;
  sline >> y;
  sline >> ey;
ex=0;
if( (x + xaxisoffset)> xaxismin    && (x + xaxisoffset) < xaxismax){

   xaxis.push_back(x + xaxisoffset);
   er_x.push_back(ex);
   yaxis.push_back(y);
   er_y.push_back(ey);
    count++;
//    if(count%1000==0)cout<<"count = "<<count<<", xaxis = "<<t*1E9 + xaxisoffset<<", voltage = "<<v<<endl;
    if(count>100)break;
  }
 }
cout<<"count = "<<count<<" "<<xaxis.size()<<endl;
cout<<"end of file"<<endl;

}
////////////////

int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadTopMargin(0.04);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.04);
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(1);
double xaxismin1 =  0.5;
double xaxismax1 =  6.5;
double xaxismin2 =  5.5;
double xaxismax2 = 10.5;
double xaxismax3 =  11.5;
TH2F *h_frame = new TH2F("h_frame","h_frame",10,xaxismin1,xaxismax1,10,0,450.); //emul & bubble data
SetTH2(h_frame , "", "", "{}^{3}_{#Lambda} H lifetime [ps]");

TH2F *h_frame1 = new TH2F("h_frame1","h_frame1",10,xaxismin2,xaxismax2,10,0,450.); //heavy ion data
SetTH2(h_frame1 , "", "", "{}^{3}_{#Lambda} H lifetime [ps]");

TH2F *h_frame2 = new TH2F("h_frame2","h_frame2",10,xaxismin1,xaxismax2,10,0,450.); //all data
SetTH2(h_frame2 , "", "", "{}^{3}_{#Lambda} H lifetime [ps]");

TH2F *h_frame3 = new TH2F("h_frame3","h_frame3",10,xaxismin2,xaxismax2,10,0,350.); //heavy ion data
SetTH2(h_frame3 , "", "", "{}^{3}_{#Lambda} H lifetime [ps]");

TH2F *h_frame4 = new TH2F("h_frame4","h_frame4",10,xaxismin1,xaxismax3,10,0,450.); //all data + expected data
SetTH2(h_frame4 , "", "", "{}^{3}_{#Lambda} H lifetime [ps]");

TGraphAsymmErrors *tg1, *tg2, *tg3;
TGraphAsymmErrors *tg4, *tg5, *tg6;
TGraphAsymmErrors *tg7, *tg8, *tg9;
TGraphAsymmErrors *tg10, *tg11;
TGraphAsymmErrors *tg6_sys, *tg7_sys, *tg8_sys, *tg9_sys, *tg10_sys, *tg11_sys;

double xaxisoffset1=0;

TLegend *leg_Lam = new TLegend( 0.82, 0.80, 0.98, 0.95);
leg_Lam -> SetBorderSize(0);
leg_Lam -> SetFillColor(0); 
leg_Lam -> SetFillStyle(0); 
leg_Lam -> SetTextFont(22); 

 double xaxis1[1] ={1} ;double ex1[1] ={0};
 double xaxis2[1] ={2} ;double ex2[1] ={0};
 double xaxis3[1] ={3} ;double ex3[1] ={0};
 double xaxis4[1] ={4} ;double ex4[1] ={0};
 double xaxis5[1] ={5} ;double ex5[1] ={0};
 double xaxis6[1] ={6} ;double ex6[1] ={0};double ex6_sys[1] ={0.1};
 double xaxis7[1] ={7} ;double ex7[1] ={0};double ex7_sys[1] ={0.1};
 double xaxis8[1] ={8} ;double ex8[1] ={0};double ex8_sys[1] ={0.1};
 double xaxis9[1] ={9} ;double ex9[1] ={0};double ex9_sys[1] ={0.1};
 double xaxis10[1]={10};double ex10[1]={0};double ex10_sys[1]={0.1};
 double xaxis11[1]={11};double ex11[1]={0};double ex11_sys[1]={0.1};
/* world data
Author,                   age, paper,              method,             ¦Ó mean,upper error,lower error,syst error
"R.J.Prem, P.H.Steinberg",1964,PR136(1964)B401,    bubble chamber,     90 ,   220,        40 ,        ???
G.Keyes,                  1968,PRL20(198)819,      bubble chamber,     232,   45 ,        34 ,        ???  //should not be included??
"P.E.Phillips, J.Schneps",1969,PR180(1969)1307,    emulsion,           285,   127,        105,        ???
G.Bohm,                   1970,NPB16(1970)46-52,   emulsion,           128,   35 ,        26 ,        ???
G.Keyes,                  1970,PRD1(1970)66,       bubble chamber,     228,   46 ,        33 ,        ???
G.Keyes,                  1973,NPB67(1973)269,     bubble chamber,     246,   62 ,        41 ,        ???
STAR            ,         2010,Science328(2010)58, heavy-ion collision,182,   89 ,        45 ,        27
C.Rappold(HypHI),         2013,NPA913(2013)170-184,heavy-ion collision,183,   42 ,        32 ,        37
ALICE,                    2016,PLB754(2016)360-372,heavy-ion collision,181,   54 ,        39 ,        33
ALICE preliminary,        2017,EPS-HEP2017        ,heavy-ion collision,223,   41 ,        33 ,        20 
expected result by ELPHexp 2020,                   TDL                 200,   
*/
 double life1[1]={90 }        ;double e_life1_high[1]={220};  double e_life1_low[1]={40 }; double e_life1_sys[1]={0};
// double life2[1]={232}        ;double e_life2_high[1]={45 };  double e_life2_low[1]={34 }; double e_life2_sys[1]={0};
 double life2[1]={285}        ;double e_life2_high[1]={127};  double e_life2_low[1]={105}; double e_life2_sys[1]={0};
 double life3[1]={128}        ;double e_life3_high[1]={35 };  double e_life3_low[1]={26 }; double e_life3_sys[1]={0};
 double life4[1]={228}        ;double e_life4_high[1]={46 };  double e_life4_low[1]={33 }; double e_life4_sys[1]={0};
 double life5[1]={246}        ;double e_life5_high[1]={62 };  double e_life5_low[1]={41 }; double e_life5_sys[1]={0};
 double life6[1]={182}        ;double e_life6_high[1]={89 };  double e_life6_low[1]={45 }; double e_life6_sys[1]={27};
 double life7[1]={183}        ;double e_life7_high[1]={42 };  double e_life7_low[1]={32 }; double e_life7_sys[1]={37};
 double life8[1]={181}        ;double e_life8_high[1]={54 };  double e_life8_low[1]={39 }; double e_life8_sys[1]={33};
 double life9[1]={142}        ;double e_life9_high[1]={24};   double e_life9_low[1]={21};  double e_life9_sys[1]={31};
 double life10[1]={223}       ;double e_life10_high[1]={41 }; double e_life10_low[1]={33 };double e_life10_sys[1]={20};
 double life11[1]={200}       ;double e_life11_high[1]={12 }; double e_life11_low[1]={12 };double e_life11_sys[1]={9.4};
 
tg1    = new TGraphAsymmErrors(1, xaxis1, life1 , ex1    , ex1    , e_life1_low , e_life1_high); 
tg2    = new TGraphAsymmErrors(1, xaxis2, life2 , ex2    , ex2    , e_life2_low , e_life2_high); 
tg3    = new TGraphAsymmErrors(1, xaxis3, life3 , ex3    , ex3    , e_life3_low , e_life3_high); 
tg4    = new TGraphAsymmErrors(1, xaxis4, life4 , ex4    , ex4    , e_life4_low , e_life4_high); 
tg5    = new TGraphAsymmErrors(1, xaxis5, life5 , ex5    , ex5    , e_life5_low , e_life5_high); 
tg6    = new TGraphAsymmErrors(1, xaxis6, life6 , ex6    , ex6    , e_life6_low , e_life6_high); 
tg7    = new TGraphAsymmErrors(1, xaxis7, life7 , ex7    , ex7    , e_life7_low , e_life7_high); 
tg8    = new TGraphAsymmErrors(1, xaxis8, life8 , ex8    , ex8    , e_life8_low , e_life8_high); 
tg9    = new TGraphAsymmErrors(1, xaxis9, life9 , ex9    , ex9    , e_life9_low , e_life9_high); 
tg10   = new TGraphAsymmErrors(1, xaxis10, life10 , ex10    , ex10    , e_life10_low , e_life10_high); 
tg11   = new TGraphAsymmErrors(1, xaxis11, life11 , ex11    , ex11    , e_life11_low , e_life11_high); 
tg6_sys= new TGraphAsymmErrors(1, xaxis6, life6 , ex6_sys, ex6_sys, e_life6_sys , e_life6_sys ); 
tg7_sys= new TGraphAsymmErrors(1, xaxis7, life7 , ex7_sys, ex7_sys, e_life7_sys , e_life7_sys ); 
tg8_sys= new TGraphAsymmErrors(1, xaxis8, life8 , ex8_sys, ex8_sys, e_life8_sys , e_life8_sys ); 
tg9_sys= new TGraphAsymmErrors(1, xaxis9, life9 , ex9_sys, ex9_sys, e_life9_sys , e_life9_sys ); 
tg10_sys= new TGraphAsymmErrors(1, xaxis10, life10 , ex10_sys, ex10_sys, e_life10_sys , e_life10_sys ); 
tg11_sys= new TGraphAsymmErrors(1, xaxis11, life11 , ex11_sys, ex11_sys, e_life11_sys , e_life11_sys ); 
//SetGrErr(TGraphErrors *gr, int LColor, int MColor, int MStyle, double Msize){
SetGrAsymmErr(tg1, 1, 1, 20, 6.);
SetGrAsymmErr(tg2, 1, 1, 21, 6.);
SetGrAsymmErr(tg3, 1, 1, 26, 6.);
SetGrAsymmErr(tg4, 1, 1, 27, 6.);
SetGrAsymmErr(tg5, 1, 1, 33, 6.);
SetGrAsymmErr(tg6, 1, 1, 30, 6.);
SetGrAsymmErr(tg7, 1, 1, 27, 6.);
SetGrAsymmErr(tg8, 1, 1, 28, 6.);
SetGrAsymmErr(tg9, 1, 1, 30, 6.);
SetGrAsymmErr(tg10,1, 1, 32, 6.);
SetGrAsymmErr(tg11,1, 2, 29, 7.);
SetGrAsymmErr(tg6_sys,  1, 0, 30, 3.);
SetGrAsymmErr(tg7_sys,  1, 0, 27, 3.);
SetGrAsymmErr(tg8_sys,  1, 0, 28, 3.);
SetGrAsymmErr(tg9_sys,  1, 0, 30, 3.);
SetGrAsymmErr(tg10_sys, 1, 0, 32, 3.);
SetGrAsymmErr(tg11_sys, 1, 0, 29, 3.);
tg6_sys ->SetFillColor(0); 
tg7_sys ->SetFillColor(0); 
tg8_sys ->SetFillColor(0); 
tg9_sys ->SetFillColor(0); 
tg10_sys->SetFillColor(0); 
tg11_sys->SetFillColor(0); 
tg6_sys ->SetFillStyle(3000); 
tg7_sys ->SetFillStyle(3000); 
tg8_sys ->SetFillStyle(3000); 
tg9_sys ->SetFillStyle(3000); 
tg10_sys->SetFillStyle(3000); 
tg11_sys->SetFillStyle(3000); 

TF1 *Lam_life = new TF1("Lam_life","[0]",0,11.5);
Lam_life -> SetParameter(0,263.2);
Lam_life -> SetLineColor(4);
Lam_life -> SetLineStyle(2);
//  TLine *l_mppc1Qmax  = new TLine(-2.1,mppc1Qmax ,2.1,mppc1Qmax );l_mppc1Qmax ->SetLineWidth(1); l_mppc1Qmax ->SetLineColor(6);
leg_Lam -> AddEntry(Lam_life,"Free #Lambda","l");

TF1 *Kama_line = new TF1("Kama_line","[0]",0,10.5);
Kama_line -> SetParameter(0,256.);
Kama_line -> SetLineColor(2);
Kama_line -> SetLineStyle(1);

TCanvas *c[5];
for(int i=0;i<5;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),3000,1600);}
c[0]->Clear();
h_frame -> Draw();
tg1     -> Draw("sameP");
tg2     -> Draw("sameP");
tg3     -> Draw("sameP");
tg4     -> Draw("sameP");
tg5     -> Draw("sameP");
tg6     -> Draw("sameP");
tg7     -> Draw("sameP");
tg8     -> Draw("sameP");
tg9     -> Draw("sameP");
tg10    -> Draw("sameP");
tg11    -> Draw("sameP");
tg7_sys -> Draw("same[]5");
tg8_sys -> Draw("same[]5");
tg9_sys -> Draw("same[]5");
tg10_sys-> Draw("same[]5");
tg11_sys-> Draw("same[]5");
Lam_life ->Draw("same");
leg_Lam -> Draw("same");

c[1]->Clear();
h_frame2 -> Draw();
tg1     -> Draw("sameP");
tg2     -> Draw("sameP");
tg3     -> Draw("sameP");
tg4     -> Draw("sameP");
tg5     -> Draw("sameP");
tg6     -> Draw("sameP");
tg7     -> Draw("sameP");
tg8     -> Draw("sameP");
tg9     -> Draw("sameP");
tg10    -> Draw("sameP");
tg11    -> Draw("sameP");
tg6_sys -> Draw("same[]5");
tg7_sys -> Draw("same[]5");
tg8_sys -> Draw("same[]5");
tg9_sys -> Draw("same[]5");
tg10_sys-> Draw("same[]5");
tg11_sys-> Draw("same[]5");
Lam_life ->Draw("same");
Kama_line->Draw("same");
leg_Lam -> Draw("same");

c[2]->Clear();
h_frame2-> Draw();
tg1     -> Draw("sameP");
tg2     -> Draw("sameP");
tg3     -> Draw("sameP");
tg4     -> Draw("sameP");
tg5     -> Draw("sameP");
tg6     -> Draw("sameP");
tg7     -> Draw("sameP");
tg8     -> Draw("sameP");
tg9     -> Draw("sameP");
tg10    -> Draw("sameP");
tg6_sys -> Draw("same[]5");
tg7_sys -> Draw("same[]5");
tg8_sys -> Draw("same[]5");
tg9_sys -> Draw("same[]5");
tg10_sys-> Draw("same[]5");
Lam_life ->Draw("same");
leg_Lam -> Draw("same");

c[3]->Clear();
h_frame3-> Draw();
tg6      -> Draw("sameP");
tg7      -> Draw("sameP");
tg8      -> Draw("sameP");
tg9      -> Draw("sameP");
tg10     -> Draw("sameP");
tg6_sys  -> Draw("same[]5");
tg7_sys  -> Draw("same[]5");
tg8_sys  -> Draw("same[]5");
tg9_sys  -> Draw("same[]5");
tg10_sys -> Draw("same[]5");
Lam_life ->Draw("same");
Kama_line->Draw("same");

c[4]->Clear();
h_frame4-> Draw();
tg1     -> Draw("sameP");
tg2     -> Draw("sameP");
tg3     -> Draw("sameP");
tg4     -> Draw("sameP");
tg5     -> Draw("sameP");
tg6     -> Draw("sameP");
tg7     -> Draw("sameP");
tg8     -> Draw("sameP");
tg9     -> Draw("sameP");
tg10    -> Draw("sameP");
tg6_sys -> Draw("same[]5");
tg7_sys -> Draw("same[]5");
tg8_sys -> Draw("same[]5");
tg9_sys -> Draw("same[]5");
tg10_sys-> Draw("same[]5");
Lam_life ->Draw("same");
leg_Lam -> Draw("same");

c[0] ->Print("pdf/life_H3L.pdf["  );
c[0] ->Print("pdf/life_H3L.pdf"   );
c[1] ->Print("pdf/life_H3L.pdf"   );
c[2] ->Print("pdf/life_H3L.pdf"   );
c[3] ->Print("pdf/life_H3L.pdf"   );
c[4] ->Print("pdf/life_H3L.pdf"   );
c[4] ->Print("pdf/life_H3L.pdf]"  );

    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
