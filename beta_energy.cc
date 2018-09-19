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
//SetGraph(tg, linecolor, count, time, intensity);
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
double GetTime(vector<double> &time,vector<double> &intensity, double Vth){
if( time.size() != intensity.size() )return -1;
int Tth=-9999;

for(int i=0;i<time.size();i++){
 if( fabs(Vth-intensity[i]) < 0.015){Tth=i;break;}
 }

return Tth;
}


//++++++++++++++ class ana_Ebeta++++++++++++++++//
class ana_Ebeta
{
public:
ana_Ebeta();
~ana_Ebeta();

void MakeTGraph(); 
void SetFrame();
void SetParam();
void DrawGraph();

private:

vector<double> Ebeta1, Ebeta2, Ebeta3, Ebeta4;
vector<double> intensity1, intensity2, intensity3, intensity4;
TLegend *leg_beta;
TGaxis *axis1;
TGraph *tg_Sr_intensity;
TGraph *tg_Y_intensity;
TGraph *tg_SrY_intensity;
TSpline3 *ts3[5];
TH2F *h_frame[5];
TH2F *h_frame_zoom;
string ifname[7];
string ofname_pdf;
TLine *tl_zoom_min;
TLine *tl_zoom_max;
TCanvas *c[6];
double aa;
double Ebetamin;
double Ebetamax;
double EbetaSrmax;
double Ymin;
double Ymax;
double YSrmax;
double YYmax;
double Ebetamin_zoom;
double Ebetamax_zoom;
unsigned int count1,count2,count3,count4;

void ReadCSV(string ifname,  vector<double> &Ebeta,vector<double> &intensity,  unsigned int &count, double Ebetamin,double Ebetamax, int thinout);

 };

////////////////
ana_Ebeta::ana_Ebeta()
:Ebeta1(), Ebeta2(), Ebeta3(), Ebeta4(),
 intensity1(), intensity2(), intensity3(), intensity4()
 {
  }
////////////////
ana_Ebeta::~ana_Ebeta(){
}
////////////////
////////////////
void ana_Ebeta::ReadCSV(string ifname,  vector<double> &Ebeta,vector<double> &intensity,  unsigned int &count, double Ebetamin,double Ebetamax, int thinout){ 
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
cout<<"count = "<<count<<" "<<Ebeta.size()<<endl;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
    if(line[0]=='t') {
  istringstream sline(line);
  sline>>name; //reading Ebetaoffset parameter from file
    sline >> offset;
    continue;
    }
  istringstream sline(line);
  sline >> w;
  sline >> e;
  aa++;
if(aa%thinout==0){
   Ebeta.push_back(w);
   intensity.push_back(e);
    count++;
    if(count%10==0)cout<<"count = "<<count<<", Ebeta = "<<w<<", intensity/gain = "<<e<<endl;
//    if(count>100)break;
   }
  }
cout<<"count = "<<count<<" "<<Ebeta.size()<<endl;
cout<<"end of file"<<endl;

 } 

////////////////
void  ana_Ebeta::SetParam(){
 ifname[0]="./csv/Sr90.csv" ;
 ifname[1]="./csv/Y90.csv" ;
 ifname[2]="./csv/Sr_Y90.csv"; 
 ofname_pdf = "./pdf/beta_energy.pdf";
 Ebetamin      =     0;
 Ebetamax      =     3;
 EbetaSrmax    =   0.6;
 Ymin         =      0;
 Ymax         =    3.5;
 YSrmax       =      3;
 YYmax        =    0.8;
 Ebetamin_zoom =     0;
 Ebetamax_zoom =   2.5;

h_frame[0] = new TH2F("h_frame1","h_frame1",10,Ebetamin,EbetaSrmax,10,Ymin,YSrmax);//Sr
h_frame[1] = new TH2F("h_frame2","h_frame2",10,Ebetamin,Ebetamax  ,10,Ymin,YYmax);//Y
h_frame[2] = new TH2F("h_frame2","h_frame2",10,Ebetamin,Ebetamax  ,10,Ymin,Ymax);//Sr+Y
SetTH2(h_frame[0] , "", "Energy [MeV]", "Num. of beta particle/1MeV/trans.");
SetTH2(h_frame[1] , "", "Energy [MeV]", "Num. of beta particle/1MeV/trans.");
SetTH2(h_frame[2] , "", "Energy [MeV]", "Num. of e^{-}/1MeV/trans.");
//SetTH2(h_frame[2] , "", "Energy [MeV]", "Num. of beta particle/1MeV/trans.");

//axis1 = new TGaxis(10,0,10,60,0,6E6,505,"+L");//(xmin,ymin,xmax,ymax,wmin,wmax,ndiv,chopt,gridrength)
//axis1 -> SetMaxDigits(1);
//axis1->SetTitleOffset(0.78);
//axis1->SetLabelSize(0.045);

leg_beta= new TLegend( 0.70, 0.60, 0.96, 0.90);
leg_beta -> SetBorderSize(1);
leg_beta -> SetFillColor(10); 
leg_beta -> SetFillStyle(1001); 
leg_beta -> SetTextFont(42); 

  tl_zoom_min = new TLine(Ebetamin_zoom,Ymin,Ebetamin_zoom,Ymax);
  tl_zoom_min->SetLineStyle(2);tl_zoom_min->SetLineWidth(1); tl_zoom_min->SetLineColor(6);
  tl_zoom_max = new TLine(Ebetamax_zoom,Ymin,Ebetamax_zoom,Ymax);
  tl_zoom_max->SetLineStyle(2);tl_zoom_max->SetLineWidth(1); tl_zoom_max->SetLineColor(6);

for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2400,1600);}

 Ebeta1.clear(); Ebeta2.clear(); Ebeta3.clear(); Ebeta4.clear();
 intensity1.clear(); intensity2.clear(); intensity3.clear(); intensity4.clear();

 }

////////////////////////////////
void ana_Ebeta::MakeTGraph(){

ReadCSV(ifname[0], Ebeta1,intensity1,count1, Ebetamin, Ebetamax ,2);
ReadCSV(ifname[1], Ebeta2,intensity2,count2, Ebetamin, Ebetamax ,2);
ReadCSV(ifname[2], Ebeta3,intensity3,count3, Ebetamin, Ebetamax ,1);
tg_Sr_intensity   = new TGraph(count1, &Ebeta1[0], &intensity1[0]); 
tg_Y_intensity    = new TGraph(count2, &Ebeta2[0], &intensity2[0]); 
tg_SrY_intensity  = new TGraph(count3, &Ebeta3[0], &intensity3[0]); 
SetGr(tg_Sr_intensity  , 2, 2, 2);//linewidth,linestyle,linecolor
SetGr(tg_Y_intensity   , 2, 2, 4);//linewidth,linestyle,linecolor
SetGr(tg_SrY_intensity , 2, 1, 1);//linewidth,linestyle,linecolor
ts3[0] = new TSpline3("ts3_1",&Ebeta1[0], &intensity1[0],count1);
ts3[1] = new TSpline3("ts3_2",&Ebeta2[0], &intensity2[0],count2);
//ts3[2] = new TSpline3("ts3_3",&Ebeta3[0], &intensity3[0],count3);

leg_beta -> AddEntry( tg_Sr_intensity   ,"{}^{90}Sr"          ,"l");
leg_beta -> AddEntry( tg_Y_intensity    ,"{}^{90}Y"           ,"l");
leg_beta -> AddEntry( tg_SrY_intensity  ,"{}^{90}Sr+{}^{90}Y" ,"l");

DrawGraph();
}
////////////////////////////////

void ana_Ebeta::DrawGraph(){
cout<<"start DrawGraph"<<endl;
c[0]->Clear();
h_frame[0]->Draw();
//axis1  ->Draw("same");
tg_Sr_intensity   ->Draw("sameC");
//ts3[0]->Draw("samel");

c[1]->Clear();
h_frame[1]->Draw();
tg_Y_intensity    ->Draw("sameC");
//ts3[1]->Draw("samel");

c[2]->Clear();
h_frame[2]->Draw();
tg_SrY_intensity  ->Draw("sameC");
//ts3[2]->Draw("samel");

c[3]->Clear();
h_frame[2]->Draw();
tg_Sr_intensity   ->Draw("sameC");
tg_Y_intensity    ->Draw("sameC");
tg_SrY_intensity  ->Draw("sameC");
leg_beta          ->Draw("same"); 

//  ofname_pdf.append(".pdf");
  c[0]->Print(Form("%s[",ofname_pdf.c_str()));
  c[0]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[2]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[3]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[3]->Print(Form("%s]",ofname_pdf.c_str()));

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
ana_Ebeta *E_beta= new ana_Ebeta(); 
E_beta -> SetParam();
E_beta -> MakeTGraph();
cout<<"Yes!!"<<endl;
    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
