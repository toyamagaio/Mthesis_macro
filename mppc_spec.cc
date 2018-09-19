#include <iostream>
#include <sstream>
#include <iomanip>
#include <csignal>
#include <stdlib.h>
#include <climits>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
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
#include "TProfile.h"
#include "TString.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TGaxis.h"
#include "TSpline.h"

#include "TRandom.h"

void SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset, double min, double max){
  gr->SetTitle(hname);
  gr->SetName(hname);
  gr->GetXaxis()->SetTitle(xname);
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->SetTitle(yname);
  gr->GetYaxis()->CenterTitle();
  gr->SetLineColor(LColor);
  gr->SetMarkerStyle(MStyle);
  gr->SetMarkerColor(MColor);
  gr->SetMarkerSize(0.5);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}

void SetTH2(TH2F *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);

  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->SetTitleOffset(0.85);
//  h2->GetXaxis()->CenterTitle();
  h2->GetXaxis()->SetNdivisions(505);
  h2->GetXaxis()->SetLabelSize(0.05);

  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->SetTitleOffset(0.78);
  h2->GetYaxis()->CenterTitle();
  h2->GetYaxis()->SetLabelSize(0.045);
  ((TGaxis*)h2->GetYaxis())->SetMaxDigits(4);

  h2->SetMinimum(0.8);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.06,"x");
  h2->SetTitleSize(0.06,"y");
}

////////////////////////////////////////////////////////////////
void SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
//SetGraph(tg, linecolor, count, time, eff);
tg -> SetMarkerStyle(21);
tg -> SetMarkerSize(2.1);
tg -> SetMarkerColor(6);
tg -> SetLineWidth(linewidth);
tg -> SetLineStyle(linestyle);
tg -> SetLineColor(linecolor);
}
////////////////////////////////////////////////////////////////
void SetSpl(TSpline3 *spl3, TString name, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize){
  spl3->SetTitle(name);
  spl3->SetName(name);

  spl3->SetLineColor(LColor);
  spl3->SetLineWidth(LWidth);
  spl3->SetLineStyle(LStyle);
  spl3->SetMarkerStyle(MStyle);
  spl3->SetMarkerColor(MColor);
  spl3->SetMarkerSize(MSize);

}

////////////////////////////////////////////////////////////////
double GetTime(vector<double> &time,vector<double> &eff, double Vth){
if( time.size() != eff.size() )return -1;
int Tth=-9999;

for(int i=0;i<time.size();i++){
 if( fabs(Vth-eff[i]) < 0.015){Tth=i;break;}
 }

return Tth;
}


//++++++++++++++ class ana_Vover++++++++++++++++//
class ana_Vover
{
public:
ana_Vover();
~ana_Vover();

void MakeTGraph(); 
void SetFrame();
void SetParam();
void DrawGraph();

private:

vector<double> Vover1, Vover2, Vover3, Vover4;
vector<double> eff1, eff2, eff3, eff4;
TLegend *leg_mppc;
TGaxis *axis1;
TGraph *tg_MPPC_eff;
TGraph *tg_MPPC_cross;
TGraph *tg_MPPC_gain;
TSpline3 *ts3[5];
TH2F *h_frame[5];
TH2F *h_frame_zoom;
string ifname[7];
string ofname_pdf;
TLine *tl_zoom_min;
TLine *tl_zoom_max;
TCanvas *c[6];
double aa;
double Vovermin;
double Vovermax;
double Ymin;
double Ymax;
double Vovermin_zoom;
double Vovermax_zoom;
unsigned int count1,count2,count3,count4;

void ReadCSV(string ifname,  vector<double> &Vover,vector<double> &eff,  unsigned int &count, double Vovermin,double Vovermax, int thinout);

 };

////////////////
ana_Vover::ana_Vover()
:Vover1(), Vover2(), Vover3(), Vover4(),
 eff1(), eff2(), eff3(), eff4()
 {
  }
////////////////
ana_Vover::~ana_Vover(){
}
////////////////
////////////////
void ana_Vover::ReadCSV(string ifname,  vector<double> &Vover,vector<double> &eff,  unsigned int &count, double Vovermin,double Vovermax, int thinout){ 
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double w,e, offset;
int aa=0;
string name;
cout<<"file name :"<<ifname<<endl;
count=0;
cout<<"count = "<<count<<" "<<Vover.size()<<endl;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
    if(line[0]=='t') {
  istringstream sline(line);
  sline>>name; //reading Voveroffset parameter from file
    sline >> offset;
    continue;
    }
  istringstream sline(line);
  sline >> w;
  sline >> e;
  aa++;
if(aa%thinout==0){
   Vover.push_back(w);
   eff.push_back(e);
    count++;
    if(count%10==0)cout<<"count = "<<count<<", Vover = "<<w<<", eff/gain = "<<e<<endl;
//    if(count>100)break;
   }
  }
cout<<"count = "<<count<<" "<<Vover.size()<<endl;
cout<<"end of file"<<endl;

 } 

////////////////
void  ana_Vover::SetParam(){
 ifname[0]="./csv/MPPC_eff50um_450nm.csv" ;
 ifname[1]="./csv/MPPC_cross50um.csv" ;
 ifname[2]="./csv/MPPC_gain50um.csv"; 
 ofname_pdf = "./pdf/mppc_spec.pdf";
 Vovermin      =     0;
 Vovermax      =    10;
 Ymin         =      0;
 Ymax         =     60;
 Vovermin_zoom =   300;
 Vovermax_zoom =   600;

h_frame[0] = new TH2F("h_frame1","h_frame1",10,Vovermin,Vovermax,10,Ymin,Ymax);
SetTH2(h_frame[0] , "", "V over [V]", "PDE, cross talk prob. [%]");
h_frame[1] = new TH2F("h_frame2","h_frame2",10,Vovermin,Vovermax,10,0,6E6);
SetTH2(h_frame[1] , "", "V over [V]", "gain");

h_frame_zoom = new TH2F("h_frame_zoom","h_frame_zoom",10,Vovermin_zoom,Vovermax_zoom,10,Ymin,Ymax);
SetTH2(h_frame_zoom , "", "wave length [nm]", "efficiency [%]");

h_frame[2] = new TH2F("h_frame3","h_frame3",10,Vovermin_zoom,Vovermax_zoom,10,0,1.05);
SetTH2(h_frame[2] , "", "wave length [nm]", "Amplitude");

h_frame[3] = new TH2F("h_frame4","h_frame4",10,Vovermin_zoom,Vovermax_zoom,10,0,100);
SetTH2(h_frame[3] , "", "wave length [nm]", "Num. of P.E. [/MeV/nm]");

h_frame[4] = new TH2F("h_frame5","h_frame5",10,Vovermin_zoom,Vovermax_zoom,10,0,100);
SetTH2(h_frame[4] , "PE", "wave length [nm]", "Num. of Ph.");

axis1 = new TGaxis(10,0,10,60,0,6E6,505,"+L");//(xmin,ymin,xmax,ymax,wmin,wmax,ndiv,chopt,gridrength)
//axis1 -> SetMaxDigits(1);
axis1->SetTitleOffset(0.78);
axis1->SetLabelSize(0.045);

leg_mppc= new TLegend( 0.10, 0.70, 0.35, 0.90);
leg_mppc -> SetBorderSize(1);
leg_mppc -> SetFillColor(10); 
leg_mppc -> SetFillStyle(1001); 
leg_mppc -> SetTextFont(42); 

  tl_zoom_min = new TLine(Vovermin_zoom,Ymin,Vovermin_zoom,Ymax);
  tl_zoom_min->SetLineStyle(2);tl_zoom_min->SetLineWidth(1); tl_zoom_min->SetLineColor(6);
  tl_zoom_max = new TLine(Vovermax_zoom,Ymin,Vovermax_zoom,Ymax);
  tl_zoom_max->SetLineStyle(2);tl_zoom_max->SetLineWidth(1); tl_zoom_max->SetLineColor(6);

for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2400,1600);}

 Vover1.clear(); Vover2.clear(); Vover3.clear(); Vover4.clear();
 eff1.clear(); eff2.clear(); eff3.clear(); eff4.clear();

 }

////////////////////////////////
void ana_Vover::MakeTGraph(){

ReadCSV(ifname[0], Vover1,eff1,count1, Vovermin, Vovermax ,1);
ReadCSV(ifname[1], Vover2,eff2,count2, Vovermin, Vovermax ,1);
ReadCSV(ifname[2], Vover3,eff3,count3, Vovermin, Vovermax ,1);
tg_MPPC_eff    = new TGraph(count1, &Vover1[0], &eff1[0]); 
tg_MPPC_cross  = new TGraph(count2, &Vover2[0], &eff2[0]); 
tg_MPPC_gain   = new TGraph(count3, &Vover3[0], &eff3[0]); 
SetGr(tg_MPPC_eff    ,4,2 ,1);//linewidth,linestyle,linecolor
SetGr(tg_MPPC_cross  ,4,1 ,4);//linewidth,linestyle,linecolor
SetGr(tg_MPPC_gain   ,4,2 ,2);//linewidth,linestyle,linecolor

leg_mppc -> AddEntry(tg_MPPC_eff   ,"PDE(#lambda=450 nm)" ,"l");
leg_mppc -> AddEntry(tg_MPPC_cross ,"cross talk" ,"l");
//leg_mppc -> AddEntry(tg_MPPC_gain  ,"gain" ,"l");

DrawGraph();
}
////////////////////////////////

void ana_Vover::DrawGraph(){
cout<<"start DrawGraph"<<endl;
c[0]->Clear();
h_frame[0]->Draw();
//axis1  ->Draw("same");
tg_MPPC_eff   ->Draw("samel");
tg_MPPC_cross ->Draw("samel");
leg_mppc      ->Draw("same");

c[1]->Clear();
h_frame[1]->Draw();
tg_MPPC_gain  ->Draw("samel");

//  ofname_pdf.append(".pdf");
  c[0]->Print(Form("%s[",ofname_pdf.c_str()));
  c[0]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1]->Print(Form("%s]",ofname_pdf.c_str()));

  }

//+++++++++++++++++++++++++++++++++//
int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadGridX(0);
gStyle->SetPadGridY(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadBottomMargin(0.17);
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(1);
ana_Vover *V_over = new ana_Vover(); 
V_over -> SetParam();
V_over -> MakeTGraph();
cout<<"Yes!!"<<endl;
    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
