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
#include "TProfile.h"
#include "TString.h"
#include "TPaveText.h"
#include "TSystem.h"

#include "TRandom.h"
const double Mp = 938.272046;          // proton       mass (MeV/c2)

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
  gr->SetMarkerSize(4.0);
  gr->GetYaxis()->SetTitleOffset(Yoffset);
//  gr->GetYaxis()->SetRangeUser(min,max);
}

void SetTH2(TH2F *h2, TString hname, TString xname, TString yname){
  h2->SetTitle(hname);
  h2->GetXaxis()->SetTitle(xname);
  h2->GetXaxis()->SetTitleOffset(0.65);
//  h2->GetXaxis()->CenterTitle();
  h2->GetXaxis()->SetNdivisions(505);
  h2->GetYaxis()->SetTitle(yname);
  h2->GetYaxis()->SetTitleOffset(0.52);
  h2->GetYaxis()->CenterTitle();
  h2->SetMinimum(0.8);
  h2->SetLineWidth(0);
  h2->SetTitleSize(0.06,"x");
  h2->SetTitleSize(0.06,"y");
}

double GetMom(double M,double Ekin){
double E = M+Ekin;
double mom = sqrt(E*E-M*M);
return mom;
}
////////////////
void ReadData(string ifname,  vector<double> &xaxis,vector<double> &yaxis,  vector<double> &er_x,vector<double> &er_y, int &count, double Mass){
  ifstream ifs( ifname.c_str() );
  if ( ifs.fail() ) {
    std::cout << "file open fail : " << ifname << std::endl;
    return;
  }
  std::cout<<ifname.c_str()<<std::endl;
  string line;
double x,y,ey,ex;
double Ekin,dEdx_e,dEdx_n,Range,mom;

count=0;
  while(!ifs.eof()){
    getline(ifs,line);
    if(line[0]=='#') continue;
  istringstream sline(line);
//data structure expected
/* Ekin dEdx(Elec) dEdx(Nucl.) Range*/
//for each line
  sline >> Ekin;
  sline >> dEdx_e;
  sline >> dEdx_n;
  sline >> Range;
mom = GetMom(Mass,Ekin);
      ex=ey=0;

   xaxis.push_back(mom);
   er_x.push_back(ex);
   yaxis.push_back(Range);
   er_y.push_back(ey);
    count++;
//    if(count%1000==0)cout<<"count = "<<count<<", xaxis = "<<t*1E9 + xaxisoffset<<", voltage = "<<v<<endl;
//    if(count>100)break;
 }
cout<<"count = "<<count<<" "<<xaxis.size()<<endl;
cout<<"end of file"<<endl;

}
////////////////

int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadGridX(0);
gStyle->SetPadGridY(0);
gStyle->SetPadTopMargin(0.1);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadBottomMargin(0.17);
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(1);
string ifname1 = "./Hydrogen_in_PilotB.txt";
double xaxismin =   0;
double xaxismax =   600;
double Ymin =   0;
double Ymax = 10;
TH2F *h_frame = new TH2F("h_frame","h_frame",10,xaxismin,xaxismax,10,Ymin,Ymax);
SetTH2(h_frame , "", "momentum[MeV/#it{c}]", "range[g/cm^{2}]");
TH2F *h_frame2 = new TH2F("h_frame2","h_frame2",10,xaxismin,xaxismax,10,60,90.);
SetTH2(h_frame2 , "ToF resolution", "Vb [V]", "RMS_{ToF}[ps]");
TGraphErrors *tg1, *tg2, *tg3, *tg4, *tg5, *tg6;

vector<double> xaxis1, xaxis2, xaxis3, xaxis4, xaxis5, xaxis6;
vector<double> ex1, ex2, ex3, ex4, ex5, ex6;
vector<double> range1, range2, range3, range4, range5, range6;
vector<double> e_range1, e_range2, e_range3, e_range4, e_range5, e_range6;
int count1,count2,count3,count4,count5,count6;
double xaxisoffset=0;

TLegend *leg_range = new TLegend( 0.70, 0.20, 0.95, 0.50);
leg_range -> SetBorderSize(1);
leg_range -> SetFillColor(0);
leg_range -> SetFillStyle(1);
leg_range -> SetTextFont(22);
TCanvas *c[5];
for(int i=0;i<5;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),2400,1600);}
ReadData(ifname1,  xaxis1, range1, ex1,e_range1, count1, Mp);
tg1 = new TGraphErrors(count1, &xaxis1[0], &range1[0] , &ex1[0], &e_range1[0] ); 
//SetGrErr(TGraphErrors *gr, TString hname, TString xname, TString yname, int LColor, int MColor, int MStyle, double Yoffset, double min, double max);
SetGrErr(tg1, "", "", "", 1, 2, 24, 0, 0, 0);

c[0]->Clear();
c[0]->SetTicks(0,1);
h_frame->Draw();
tg1  ->Draw("sameP");

c[1]->Clear();
c[1]->SetTicks(0,1);
h_frame->Draw();
tg1 ->Draw("sameL");

c[0]->Print("./pdf/proton_range.pdf["  );
c[0]->Print("./pdf/proton_range.pdf"   );
c[1]->Print("./pdf/proton_range.pdf"   );
c[1]->Print("./pdf/proton_range.pdf]"  );

    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
