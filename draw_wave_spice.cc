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
  ((TGaxis*)h2->GetYaxis())->SetMaxDigits(3);

  h2->SetMinimum(0.8);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.06,"x");
  h2->SetTitleSize(0.06,"y");
}

////////////////////////////////////////////////////////////////
void SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
//SetGraph(tg, linecolor, count, time, volt);
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
double GetTime(vector<double> &time,vector<double> &volt, double Vth){
if( time.size() != volt.size() )return -1;
int Tth=-9999;

for(int i=0;i<time.size();i++){
 if( fabs(Vth-volt[i]) < 0.015){Tth=i;break;}
 }

return Tth;
}


//++++++++++++++ class ana_waveform++++++++++++++++//
class ana_waveform
{
public:
ana_waveform();
~ana_waveform();

void MakeTGraph(); 
void SetFrame();
void SetParam();
void DrawGraph();

private:

vector<double> time1, time2, time3, time4;
vector<double> volt1, volt2, volt3, volt4;
vector<double> dvdt1, dvdt2, dvdt3, dvdt4;
vector<double> volt_a1, volt_a2, volt_a3, volt_a4;
vector<double> volt_b1, volt_b2, volt_b3, volt_b4;
TGraph *tg[5];
TGraph *tg_dif[5];
TSpline3 *ts3[5];
TSpline3 *ts3_dif[5];
TH2F *h_frame[5];
TH2F *h_frame_zoom;
TH2F *h_frame_dif;
string ifname[5];
string ofname_pdf;
TLine *tl_zoom_min;
TLine *tl_zoom_max;
TCanvas *c[6];
double timemin;
double timemax;
double Vmin,Vmin1;
double Vmax;
double timemin_zoom;
double timemax_zoom;
unsigned int count1,count2,count3,count4;
double timeoffset1;
double timeoffset2;
double timeoffset3;
double timeoffset4;

void ReadCSV(string ifname,  vector<double> &time,vector<double> &volt, vector<double> &dvdt, unsigned int &count, double timemin,double timemax, double timeoffset); 
void ReadTxt(string ifname,  vector<double> &time,vector<double> &volt, vector<double> &dvdt, unsigned int &count, double timemin,double timemax, double timeoffset); 
void ReadTxt2(string ifname,  vector<double> &time,vector<double> &volt1, vector<double> &volt2, unsigned int &count, double timemin,double timemax, double timeoffset); 
void ReadTxt3(string ifname,  vector<double> &time,vector<double> &volt1, vector<double> &volt2, vector<double> &volt3, unsigned int &count, double timemin,double timemax, double timeoffset); 

 };

////////////////
ana_waveform::ana_waveform()
:time1(), time2(), time3(), time4(),
 volt1(), volt2(), volt3(), volt4(),
 dvdt1(), dvdt2(), dvdt3(), dvdt4()
 {
  }
////////////////
ana_waveform::~ana_waveform(){
}
////////////////
void ana_waveform::ReadCSV(string ifname,  vector<double> &time,vector<double> &volt, vector<double> &dvdt, unsigned int &count, double timemin,double timemax, double timeoffset){ 
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double t,v,DVDT, offset;
string name;
cout<<"file name :"<<ifname<<endl;
count=0;
cout<<"count = "<<count<<" "<<time.size()<<endl;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
    if(line[0]=='t') {
  istringstream sline(line);
  sline>>name; //reading timeoffset parameter from file
    sline >> offset;
    timeoffset = offset;
    continue;
    }
  istringstream sline(line);
  sline >> t;
  sline >> v;
if( (t*1E9 + timeoffset)> timemin    && (t*1E9 + timeoffset) < timemax){
  if(count!=0) DVDT =(v - volt.at(count-1)) / (t*1E9+timeoffset - time.at(count-1) );
   else DVDT = 0;
   time.push_back(t*1E9 + timeoffset);
   volt.push_back(v);
//cout<<"dV/dt="<<(v - volt.at(count-1))<<" / "<<(t*1E9+timeoffset - time.at(count-1) )<<" = "<<DVDT<<endl;
   dvdt.push_back(DVDT);
//cout<<"count = "<<count<<" dvdt.size="<<dvdt.size()<<endl;
    count++;
    if(count%10==0)cout<<"count = "<<count<<", time = "<<t*1E9 + timeoffset<<", voltage = "<<v<<endl;
//    if(count>100)break;
    }
   }
cout<<"count = "<<count<<" "<<time.size()<<endl;
cout<<"end of file"<<endl;

 } 

////////////////
void ana_waveform::ReadTxt(string ifname,  vector<double> &time,vector<double> &volt, vector<double> &dvdt, unsigned int &count, double timemin,double timemax, double timeoffset){
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double t,v,DVDT;

count=0;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
  istringstream sline(line);
  sline >> t;
  sline >> v;
if( (t*1E9 + timeoffset)> timemin    && (t*1E9 + timeoffset) < timemax){
  if(count!=0) DVDT = (v-volt[count-1])/(t*1E9+timeoffset - time[count-1]);
   else DVDT = 0;
//   DVDT = 0;
   time.push_back(t*1E9 + timeoffset);
   volt.push_back(v);
   dvdt.push_back(DVDT);
    count++;
    if(count%1000==0)cout<<"count = "<<count<<", time = "<<t*1E9 + timeoffset<<", voltage = "<<v<<endl;
//    if(count>100)break;
  }
 }
cout<<"count = "<<count<<" "<<time.size()<<endl;
cout<<"end of file"<<endl;
 } 

////////////////
void ana_waveform::ReadTxt2(string ifname,  vector<double> &time,vector<double> &volt1, vector<double> &volt2, unsigned int &count, double timemin,double timemax, double timeoffset){
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double t,v1,v2;

count=0;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
  istringstream sline(line);
  sline >> t;
  sline >> v1;
  sline >> v2;
//if( (t*1E9 + timeoffset)> -10    && (t*1E9 + timeoffset) < timemax){
if( (t*1E9 + timeoffset)> timemin    && (t*1E9 + timeoffset) < timemax){
   time.push_back(t*1E9 + timeoffset);
   volt1.push_back(v1);
   volt2.push_back(v2);
    count++;
    if(count%1000==0)cout<<"count = "<<count<<", time = "<<t*1E9 + timeoffset<<", voltage1 = "<<v1<<", voltage2 = "<<v2<<endl;
//    if(count>100)break;
  }
 }
cout<<"count = "<<count<<" "<<time.size()<<endl;
cout<<"offset = "<<timeoffset<<endl;
cout<<"timemin= "<<timemin<<endl;

cout<<"end of file"<<endl;
 } 

////////////////
void ana_waveform::ReadTxt3(string ifname,  vector<double> &time,vector<double> &volt1, vector<double> &volt2, vector<double> &volt3, unsigned int &count, double timemin,double timemax, double timeoffset){
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  string line;
double t,v1,v2,v3;

count=0;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
  istringstream sline(line);
  sline >> t;
  sline >> v1;
  sline >> v2;
  sline >> v3;
if( (t*1E9 + timeoffset)> timemin    && (t*1E9 + timeoffset) < timemax){
   time.push_back(t*1E9 + timeoffset);
   volt1.push_back(v1);
   volt2.push_back(v2);
   volt3.push_back(v3);
    count++;
    if(count%1000==0)cout<<"count = "<<count<<", time = "<<t*1E9 + timeoffset<<", voltage1 = "<<v1<<", voltage2 = "<<v2<<", voltage3 = "<<v3<<endl;
//    if(count>100)break;
  }
 }
cout<<"count = "<<count<<" "<<time.size()<<endl;
cout<<"end of file"<<endl;
 } 

////////////////

////////////////////////////////
void  ana_waveform::SetParam(){
 ifname[0]="./waveform/LTspice/integ_ad8009_pzc_20170126.txt" ;
 ifname[1]="./waveform/LTspice/original.txt" ; 
 ifname[2]="./waveform/LTspice/add_fbC.txt" ; 
 ofname_pdf = "./pdf/spice_pzc.pdf";
 timemin      =   -10;
 timemax      =   100;
 Vmin         =  -2.2;
 Vmin1        =  -3.6;
 Vmax         =   0.6;
 timemin_zoom =     0;
 timemax_zoom =   100;
 timeoffset1  =   -20;
 timeoffset2  =   -90;
 timeoffset3  =     0;
 timeoffset4  =     0;

h_frame[0] = new TH2F("h_frame1","h_frame1",10,timemin,timemax,10,Vmin,Vmax);
SetTH2(h_frame[0] , "", "time [ns]", "Voltage [V]");
h_frame[1] = new TH2F("h_frame2","h_frame2",10,timemin,timemax,10,Vmin1,Vmax);
SetTH2(h_frame[1] , "", "time [ns]", "Voltage [V]");

h_frame_zoom = new TH2F("h_frame_zoom","h_frame_zoom",10,timemin_zoom,timemax_zoom,10,-0.4,Vmax);
SetTH2(h_frame_zoom , "", "time [ns]", "Voltage [V]");

h_frame_dif = new TH2F("h_frame_dif","h_frame_dif",10,timemin_zoom,timemax_zoom,10,-0.5,0.5);
SetTH2(h_frame_dif , "wave form differencial", "time [ns]", "dV/dt [V/ns]");

  tl_zoom_min = new TLine(timemin_zoom,Vmin,timemin_zoom,Vmax);
  tl_zoom_min->SetLineStyle(2);tl_zoom_min->SetLineWidth(1); tl_zoom_min->SetLineColor(6);
  tl_zoom_max = new TLine(timemax_zoom,Vmin,timemax_zoom,Vmax);
  tl_zoom_max->SetLineStyle(2);tl_zoom_max->SetLineWidth(1); tl_zoom_max->SetLineColor(6);

for(int i=0;i<6;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800);}

 time1.clear(); time2.clear(); time3.clear(); time4.clear();
 volt1.clear(); volt2.clear(); volt3.clear(); volt4.clear();
 dvdt1.clear(); dvdt2.clear(); dvdt3.clear(); dvdt4.clear();
 volt_a1.clear(); volt_a2.clear(); volt_a3.clear(); volt_a4.clear();

 }
////////////////////////////////
void ana_waveform::MakeTGraph(){

//cout<<"time1 size"<<time1.size()<<endl;
ReadTxt2(ifname[0], time1,volt1,volt2,count1, -20, timemax, timeoffset1);
tg[0]  = new TGraph(count1, &time1[0], &volt1[0]); 
tg[1]  = new TGraph(count1, &time1[0], &volt2[0]); 
SetGr(tg[0],2,1,2);
SetGr(tg[1],2,1,4);
ReadTxt(ifname[1], time2,volt_a1,dvdt1,count2, timemin, timemax, timeoffset2);
tg[2]  = new TGraph(count2, &time2[0], &volt_a1[0]); 
SetGr(tg[2],2,1,1);

ReadTxt(ifname[2], time3,volt_b1,dvdt2,count3, timemin, timemax, timeoffset2);
tg[3]  = new TGraph(count3, &time3[0], &volt_b1[0]); 
SetGr(tg[3],2,1,2);


DrawGraph();
}
////////////////////////////////

void ana_waveform::DrawGraph(){
cout<<"start DrawGraph"<<endl;
c[0]->Clear();
h_frame[0]->Draw();
tg[2] ->Draw("samel");
tl_zoom_min  ->Draw("same");
tl_zoom_max  ->Draw("same");

c[1]->Clear();
h_frame[0]->Draw();
tg[3] ->Draw("samel");
//tl_zoom_min  ->Draw("same");
//tl_zoom_max  ->Draw("same");

c[2]->Clear();
h_frame_zoom->Draw();
tg[2] ->Draw("samel");
tg[3] ->Draw("samel");

c[3]->Clear();
h_frame[0]->Draw();
tg[2] ->Draw("samel");
tg[3] ->Draw("samel");

c[4]->Clear();
h_frame[1]->Draw();
tg[0] ->Draw("samel");
tg[1] ->Draw("samel");

//  ofname_pdf.append(".pdf");
  c[0]->Print(Form("%s[",ofname_pdf.c_str()));
  c[0]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[2]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[3]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[4]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[5]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[5]->Print(Form("%s]",ofname_pdf.c_str()));

  }

//+++++++++++++++++++++++++++++++++//
int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadGridX(1);
gStyle->SetPadGridY(1);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadBottomMargin(0.17);
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(1);
ana_waveform *waveform = new ana_waveform(); 
waveform -> SetParam();
waveform -> MakeTGraph();

    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
