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


//++++++++++++++ class ana_wavelength++++++++++++++++//
class ana_wavelength
{
public:
ana_wavelength();
~ana_wavelength();

void MakeTGraph(); 
void SetFrame();
void SetParam();
void DrawGraph();

private:

vector<double> wavelen1, wavelen2, wavelen3, wavelen4;
vector<double> wavelen_BC418;
vector<double> wavelen_BC422;
vector<double> wavelen_BC408;
vector<double> wavelen_BC404;
vector<double> wavelen_BC400;
vector<double> emit_BC418;
vector<double> emit_BC422;
vector<double> emit_BC408;
vector<double> emit_BC404;
vector<double> emit_BC400;
vector<double> yield_BC418_CS;
vector<double> yield_BC422_CS;
vector<double> yield_BC408_CS;
vector<double> yield_BC404_CS;
vector<double> yield_BC400_CS;
vector<double> yield_BC418_PE;
vector<double> yield_BC422_PE;
vector<double> yield_BC408_PE;
vector<double> yield_BC404_PE;
vector<double> yield_BC400_PE;
vector<double> eff1, eff2, eff3, eff4;
vector<double> dvdt1, dvdt2, dvdt3, dvdt4;
vector<double> eff_a1, eff_a2, eff_a3, eff_a4;
TLegend *leg_Scint;
TGraph *tg_MPPC_eff;
TGraph *tg_MPPCPE_eff;
TGraph *tg_BC418;
TGraph *tg_BC422;
TGraph *tg_BC408;
TGraph *tg_BC404;
TGraph *tg_BC400;
TGraph *tg_BC418_eff_CS;
TGraph *tg_BC422_eff_CS;
TGraph *tg_BC408_eff_CS;
TGraph *tg_BC404_eff_CS;
TGraph *tg_BC400_eff_CS;
TGraph *tg_BC418_eff_PE;
TGraph *tg_BC422_eff_PE;
TGraph *tg_BC408_eff_PE;
TGraph *tg_BC404_eff_PE;
TGraph *tg_BC400_eff_PE;
TSpline3 *ts3[5];
TSpline3 *ts_BC418;
TSpline3 *ts_BC422;
TSpline3 *ts_BC408;
TSpline3 *ts_BC404;
TSpline3 *ts_BC400;
TH2F *h_frame[5];
TH2F *h_frame_zoom;
string ifname[7];
string ofname_pdf;
TLine *tl_zoom_min;
TLine *tl_zoom_max;
TCanvas *c[7];
double aa;
double wavelenmin;
double wavelenmax;
double Ymin;
double Ymax;
double wavelenmin_zoom;
double wavelenmax_zoom;
double integ1,integ2,integ3,integ4;
unsigned int count1,count2,count3,count4;
unsigned int count_BC418;
unsigned int count_BC422;
unsigned int count_BC408;
unsigned int count_BC404;
unsigned int count_BC400;
double scale_BC418;
double scale_BC422;
double scale_BC408;
double scale_BC404;
double scale_BC400;
double integ_BC418;
double integ_BC422;
double integ_BC408;
double integ_BC404;
double integ_BC400;
double NPE_BC418_CS;
double NPE_BC422_CS;
double NPE_BC408_CS;
double NPE_BC404_CS;
double NPE_BC400_CS;
double NPE_BC418_PE;
double NPE_BC422_PE;
double NPE_BC408_PE;
double NPE_BC404_PE;
double NPE_BC400_PE;

void ReadCSV(string ifname,  vector<double> &wavelen,vector<double> &eff,  unsigned int &count, double &integ, double wavelenmin,double wavelenmax, int thinout);
void CalcNP( TSpline3 *ts3, vector<double> &wavelen,vector<double> &emit, vector<double> &yield, double scale, double &integ);

 };

////////////////
ana_wavelength::ana_wavelength()
:wavelen1(), wavelen2(), wavelen3(), wavelen4(),
 eff1(), eff2(), eff3(), eff4(),
 dvdt1(), dvdt2(), dvdt3(), dvdt4()
 {
  }
////////////////
ana_wavelength::~ana_wavelength(){
}
////////////////
////////////////
void ana_wavelength::ReadCSV(string ifname,  vector<double> &wavelen,vector<double> &eff,  unsigned int &count, double &integ, double wavelenmin,double wavelenmax, int thinout){ 
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
cout<<"count = "<<count<<" "<<wavelen.size()<<endl;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
    if(line[0]=='t') {
  istringstream sline(line);
  sline>>name; //reading wavelenoffset parameter from file
    sline >> offset;
    continue;
    }
  istringstream sline(line);
  sline >> w;
  sline >> e;
  aa++;
if(aa%thinout==0){
   wavelen.push_back(w);
   eff.push_back(e);
    count++;
    integ += e;
    if(count%10==0)cout<<"count = "<<count<<", wavelen = "<<w<<", effage = "<<e<<endl;
//    if(count>100)break;
   }
  }
cout<<"count = "<<count<<" "<<wavelen.size()<<endl;
cout<<"end of file"<<endl;

 } 

////////////////
void ana_wavelength::CalcNP( TSpline3 *ts3, vector<double> &wavelen,vector<double> &emit, vector<double> &yield, double scale, double &integ){
int n = wavelen.size();
cout<<n<<endl;
double eff;
for(int i=0;i<n;i++){
eff = ts3->Eval(wavelen[i]) ;
yield.push_back(eff*emit[i]*scale*0.01);
integ += eff*emit[i]*scale*0.01;
if(i%10==0) cout <<"wavelen:"<<wavelen[i] <<", eff:"<<eff<<", emit:"<<emit[i]<<", yield"<<yield[i] << endl;
  }

}

////////////////////////////////
void  ana_wavelength::SetParam(){
// ifname[0]="./csv/MPPC_eff50um.csv" ;
 ifname[0]="./csv/MPPC_eff50um_CS.csv" ;
 ifname[1]="./csv/MPPC_eff50um_PE.csv" ;
 ifname[2]="./csv/EJ228_BC418.csv"; 
 ifname[3]="./csv/EJ232_BC422.csv";
 ifname[4]="./csv/EJ200_BC408.csv";
 ifname[5]="./csv/EJ204_BC404.csv";
 ifname[6]="./csv/EJ212_BC400.csv";
 ofname_pdf = "./pdf/scinti_emit.pdf";
 wavelenmin      =   300;
 wavelenmax      =   1000;
 Ymin         =     0;
 Ymax         =    50;
 wavelenmin_zoom =   300;
 wavelenmax_zoom =   600;

h_frame[0] = new TH2F("h_frame1","h_frame1",10,wavelenmin,wavelenmax,10,Ymin,Ymax);
SetTH2(h_frame[0] , "", "wave length [nm]", "PDE [%]");
h_frame[1] = new TH2F("h_frame2","h_frame2",10,wavelenmin,wavelenmax,10,0,1.05);
SetTH2(h_frame[1] , "", "wave length [nm]", "Amplitude");

h_frame_zoom = new TH2F("h_frame_zoom","h_frame_zoom",10,wavelenmin_zoom,wavelenmax_zoom,10,Ymin,Ymax);
SetTH2(h_frame_zoom , "", "wave length [nm]", "efficiency [%]");

h_frame[2] = new TH2F("h_frame3","h_frame3",10,wavelenmin_zoom,wavelenmax_zoom,10,0,1.05);
SetTH2(h_frame[2] , "", "wave length [nm]", "Amplitude");

h_frame[3] = new TH2F("h_frame4","h_frame4",10,wavelenmin_zoom,wavelenmax_zoom,10,0,100);
SetTH2(h_frame[3] , "", "wave length [nm]", "Num. of P.E. [/MeV/nm]");

h_frame[4] = new TH2F("h_frame5","h_frame5",10,wavelenmin_zoom,wavelenmax_zoom,10,0,100);
SetTH2(h_frame[4] , "PE", "wave length [nm]", "Num. of Ph.");

leg_Scint= new TLegend( 0.65, 0.50, 0.96, 0.90);
leg_Scint -> SetBorderSize(1);
leg_Scint -> SetFillColor(10); 
leg_Scint -> SetFillStyle(1001); 
leg_Scint -> SetTextFont(42); 

  tl_zoom_min = new TLine(wavelenmin_zoom,Ymin,wavelenmin_zoom,Ymax);
  tl_zoom_min->SetLineStyle(2);tl_zoom_min->SetLineWidth(1); tl_zoom_min->SetLineColor(6);
  tl_zoom_max = new TLine(wavelenmax_zoom,Ymin,wavelenmax_zoom,Ymax);
  tl_zoom_max->SetLineStyle(2);tl_zoom_max->SetLineWidth(1); tl_zoom_max->SetLineColor(6);

for(int i=0;i<7;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2400,1600);}

scale_BC418=10200.;  //[photon/MeV]
scale_BC422=8400.; //[photon/MeV]
scale_BC408=10000.; //[photon/MeV]
scale_BC404=10400.; //[photon/MeV]
scale_BC400=10000.; //[photon/MeV]

integ_BC418=0.;
integ_BC422=0.;
integ_BC408=0.;
integ_BC404=0.;
integ_BC400=0.;

NPE_BC418_CS=0.;
NPE_BC422_CS=0.;
NPE_BC408_CS=0.;
NPE_BC404_CS=0.;
NPE_BC400_CS=0.;

NPE_BC418_PE=0.;
NPE_BC422_PE=0.;
NPE_BC408_PE=0.;
NPE_BC404_PE=0.;
NPE_BC400_PE=0.;

 wavelen1.clear(); wavelen2.clear(); wavelen3.clear(); wavelen4.clear();
 eff1.clear(); eff2.clear(); eff3.clear(); eff4.clear();
 dvdt1.clear(); dvdt2.clear(); dvdt3.clear(); dvdt4.clear();
 eff_a1.clear(); eff_a2.clear(); eff_a3.clear(); eff_a4.clear();

 }

////////////////////////////////
void ana_wavelength::MakeTGraph(){

ReadCSV(ifname[0], wavelen1,eff1,count1, integ1,wavelenmin, wavelenmax ,13);
tg_MPPC_eff  = new TGraph(count1, &wavelen1[0], &eff1[0]); 
SetGr(tg_MPPC_eff  ,4,1,2);//linewidth,linestyle,linecolor
ts3[0] = new TSpline3("ts3_1",&wavelen1[0], &eff1[0],count1);
//ts3[0] = new TSpline3("ts3_1",tg[0]);
SetSpl(ts3[0],"spline",2,2,1,2,24,2.0);
//SetSpl(TSpline3 *spl3, TString name, int LColor, int LWidth, int LStyle, int MColor, int MStyle, double MSize)

ReadCSV(ifname[1], wavelen2,eff2,count2, integ2, wavelenmin, wavelenmax ,13);
tg_MPPCPE_eff  = new TGraph(count2, &wavelen2[0], &eff2[0]); 
SetGr(tg_MPPCPE_eff  ,4,7,1);//linewidth,linestyle,linecolor
ts3[1] = new TSpline3("ts3_2",&wavelen2[0], &eff2[0],count2);
SetSpl(ts3[1],"spline",1,2,7,1,25,2.0);

ReadCSV(ifname[2], wavelen_BC418,emit_BC418,count_BC418,integ_BC418, wavelenmin, wavelenmax ,1);
ReadCSV(ifname[3], wavelen_BC422,emit_BC422,count_BC422,integ_BC422, wavelenmin, wavelenmax ,1);
ReadCSV(ifname[4], wavelen_BC408,emit_BC408,count_BC408,integ_BC408, wavelenmin, wavelenmax ,1);
ReadCSV(ifname[5], wavelen_BC404,emit_BC404,count_BC404,integ_BC404, wavelenmin, wavelenmax ,1);
ReadCSV(ifname[6], wavelen_BC400,emit_BC400,count_BC400,integ_BC400, wavelenmin, wavelenmax ,1);

CalcNP( ts3[0], wavelen_BC418,emit_BC418,  yield_BC418_CS, scale_BC418/integ_BC418 , NPE_BC418_CS);
CalcNP( ts3[0], wavelen_BC422,emit_BC422,  yield_BC422_CS, scale_BC422/integ_BC422 , NPE_BC422_CS);
CalcNP( ts3[0], wavelen_BC408,emit_BC408,  yield_BC408_CS, scale_BC408/integ_BC408 , NPE_BC408_CS);
CalcNP( ts3[0], wavelen_BC404,emit_BC404,  yield_BC404_CS, scale_BC404/integ_BC404 , NPE_BC404_CS);
CalcNP( ts3[0], wavelen_BC400,emit_BC400,  yield_BC400_CS, scale_BC400/integ_BC400 , NPE_BC400_CS);

CalcNP( ts3[1], wavelen_BC418,emit_BC418,  yield_BC418_PE, scale_BC418/integ_BC418 , NPE_BC418_PE);
CalcNP( ts3[1], wavelen_BC422,emit_BC422,  yield_BC422_PE, scale_BC422/integ_BC422 , NPE_BC422_PE);
CalcNP( ts3[1], wavelen_BC408,emit_BC408,  yield_BC408_PE, scale_BC408/integ_BC408 , NPE_BC408_PE);
CalcNP( ts3[1], wavelen_BC404,emit_BC404,  yield_BC404_PE, scale_BC404/integ_BC404 , NPE_BC404_PE);
CalcNP( ts3[1], wavelen_BC400,emit_BC400,  yield_BC400_PE, scale_BC400/integ_BC400 , NPE_BC400_PE);

tg_BC418 = new TGraph(count_BC418, &wavelen_BC418[0], &emit_BC418[0]);
tg_BC422 = new TGraph(count_BC422, &wavelen_BC422[0], &emit_BC422[0]);
tg_BC408 = new TGraph(count_BC408, &wavelen_BC408[0], &emit_BC408[0]);
tg_BC404 = new TGraph(count_BC404, &wavelen_BC404[0], &emit_BC404[0]);
tg_BC400 = new TGraph(count_BC400, &wavelen_BC400[0], &emit_BC400[0]);
SetGr(tg_BC400  ,2,2,1);
SetGr(tg_BC404  ,2,2,3);
SetGr(tg_BC408  ,2,2,2);
SetGr(tg_BC418  ,2,2,4);
SetGr(tg_BC422  ,2,2,6);

tg_BC418_eff_CS = new TGraph(count_BC418, &wavelen_BC418[0], &yield_BC418_CS[0]);
tg_BC422_eff_CS = new TGraph(count_BC422, &wavelen_BC422[0], &yield_BC422_CS[0]);
tg_BC408_eff_CS = new TGraph(count_BC408, &wavelen_BC408[0], &yield_BC408_CS[0]);
tg_BC404_eff_CS = new TGraph(count_BC404, &wavelen_BC404[0], &yield_BC404_CS[0]);
tg_BC400_eff_CS = new TGraph(count_BC400, &wavelen_BC400[0], &yield_BC400_CS[0]);
SetGr(tg_BC400_eff_CS  ,2,2,1);
SetGr(tg_BC404_eff_CS  ,2,2,3);
SetGr(tg_BC408_eff_CS  ,2,2,2);
SetGr(tg_BC418_eff_CS  ,2,2,4);
SetGr(tg_BC422_eff_CS  ,2,2,6);

tg_BC418_eff_PE = new TGraph(count_BC418, &wavelen_BC418[0], &yield_BC418_PE[0]);
tg_BC422_eff_PE = new TGraph(count_BC422, &wavelen_BC422[0], &yield_BC422_PE[0]);
tg_BC408_eff_PE = new TGraph(count_BC408, &wavelen_BC408[0], &yield_BC408_PE[0]);
tg_BC404_eff_PE = new TGraph(count_BC404, &wavelen_BC404[0], &yield_BC404_PE[0]);
tg_BC400_eff_PE = new TGraph(count_BC400, &wavelen_BC400[0], &yield_BC400_PE[0]);
SetGr(tg_BC400_eff_PE  ,2,3,1);
SetGr(tg_BC404_eff_PE  ,2,3,3);
SetGr(tg_BC408_eff_PE  ,2,3,2);
SetGr(tg_BC418_eff_PE  ,2,3,4);
SetGr(tg_BC422_eff_PE  ,2,3,6);

leg_Scint -> AddEntry(tg_BC400,"BC400(EJ212)" ,"l");//"BC400"
leg_Scint -> AddEntry(tg_BC404,"BC404(EJ204)" ,"l");//"BC404"
leg_Scint -> AddEntry(tg_BC408,"BC408(EJ200)" ,"l");//"BC408"
leg_Scint -> AddEntry(tg_BC418,"BC418(EJ228)" ,"l");//"BC418"
leg_Scint -> AddEntry(tg_BC422,"BC422(EJ232)" ,"l");//"BC422"


DrawGraph();
cout<<"integ BC418:"<<integ_BC418 <<endl;
cout<<"integ BC422:"<<integ_BC422 <<endl;
cout<<"integ BC408:"<<integ_BC408 <<endl;
cout<<"integ BC404:"<<integ_BC404 <<endl;
cout<<"integ BC400:"<<integ_BC400 <<endl;
cout<<"NPE BC418(CS):"<<NPE_BC418_CS <<endl;
cout<<"NPE BC422(CS):"<<NPE_BC422_CS <<endl;
cout<<"NPE BC408(CS):"<<NPE_BC408_CS <<endl;
cout<<"NPE BC404(CS):"<<NPE_BC404_CS <<endl;
cout<<"NPE BC400(CS):"<<NPE_BC400_CS <<endl;
cout<<"NPE BC418(PE):"<<NPE_BC418_PE <<endl;
cout<<"NPE BC422(PE):"<<NPE_BC422_PE <<endl;
cout<<"NPE BC408(PE):"<<NPE_BC408_PE <<endl;
cout<<"NPE BC404(PE):"<<NPE_BC404_PE <<endl;
cout<<"NPE BC400(PE):"<<NPE_BC400_PE <<endl;
}
////////////////////////////////

void ana_wavelength::DrawGraph(){
cout<<"start DrawGraph"<<endl;
c[0]->Clear();
h_frame[0]->Draw();
//tg[0]->Draw("samep");
ts3[0] ->Draw("lcsame");//tg1->Draw("samep");
ts3[1] ->Draw("lcsame");//tg1->Draw("samep");
//tl_zoom_min  ->Draw("same");
//tl_zoom_max  ->Draw("same");

c[1]->Clear();
h_frame[0]->Draw();
//tg[0]->Draw("samep");
ts3[0] ->Draw("lcsame");//tg1->Draw("samep");

c[2]->Clear();
h_frame[1]->Draw();
tg_BC418 ->Draw("samel");
tg_BC422 ->Draw("samel");
tg_BC408 ->Draw("samel");
tg_BC404 ->Draw("samel");
tg_BC400 ->Draw("samel");

c[3]->Clear();
h_frame_zoom->Draw();
//tg[0]->Draw("samep");
ts3[0] ->Draw("samelc");//tg1->Draw("samep");
ts3[1] ->Draw("samelc");//tg1->Draw("samep");

c[4]->Clear();
h_frame[2]->Draw();
tg_BC418 ->Draw("samel");
tg_BC422 ->Draw("samel");
tg_BC408 ->Draw("samel");
tg_BC404 ->Draw("samel");
tg_BC400 ->Draw("samel");
leg_Scint ->Draw("same");

c[5]->Clear();
h_frame[3]->Draw();
tg_BC418_eff_CS ->Draw("samel");
tg_BC422_eff_CS ->Draw("samel");
tg_BC408_eff_CS ->Draw("samel");
tg_BC404_eff_CS ->Draw("samel");
tg_BC400_eff_CS ->Draw("samel");
leg_Scint ->Draw("same");

c[6]->Clear();
h_frame[4]->Draw();
tg_BC418_eff_CS ->Draw("samel");
tg_BC422_eff_CS ->Draw("samel");
tg_BC408_eff_CS ->Draw("samel");
tg_BC404_eff_CS ->Draw("samel");
tg_BC400_eff_CS ->Draw("samel");
tg_BC418_eff_PE ->Draw("samel");
tg_BC422_eff_PE ->Draw("samel");
tg_BC408_eff_PE ->Draw("samel");
tg_BC404_eff_PE ->Draw("samel");
tg_BC400_eff_PE ->Draw("samel");
leg_Scint ->Draw("same");
//  ofname_pdf.append(".pdf");
  c[0]->Print(Form("%s[",ofname_pdf.c_str()));
  c[0]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[1]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[2]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[3]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[4]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[5]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[6]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[6]->Print(Form("%s]",ofname_pdf.c_str()));

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
ana_wavelength *waveform = new ana_wavelength(); 
waveform -> SetParam();
waveform -> MakeTGraph();
cout<<"Yes!!"<<endl;
    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
