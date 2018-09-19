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

const double PI = 4.*atan(1.);
const double deg_to_rad = PI / 180.;
const double rad_to_deg = 180. / PI;
const double sigma_to_fwhm = 2.*sqrt(2.*log(2.));
const double fwhm_to_sigma = 1./sigma_to_fwhm;
const double cm_to_barn = 1e+24;
const double alpha = 1./137.035999074; // fine structure constant
const double hbar = 6.58211928*1e-22;  // Planc constant (reduced) (MeV x s)
const double hbarc = 197.3269718;      // conversion constant (MeV x fm)
const double kb = 8.6173324*1e-5;      // Boltzmann constant
const double e = 1.602176565*1e-19;    // electron charge magnitude (C)
const double c = 0.299792458;          // speed of light in vacuum (m/ns)
const double re = 2.817e-13;           // classical electron radius (cm)
const double Na = 6.02214129*1e+23;    // Avogadro constant
const double Me = 0.510998928;         // electron     mass (MeV/c2)
const double Mmu = 105.6583715;        // muon         mass (MeV/c2)
const double Mpi = 139.57018;          // charged pion mass (MeV/c2)
const double Mpi0 = 134.9766;          // charged pion mass (MeV/c2)
const double MK = 493.677;             // charged Kaon mass (MeV/c2)
const double Mp = 938.272046;          // proton       mass (MeV/c2)
const double Mn = 939.565379;          // proton       mass (MeV/c2)
const double Mu = 931.494061;          // proton       mass (MeV/c2)
const double ML = 1115.683;            // Lambda       mass (MeV/c2)
const double MS0 = 1192.642;           // Sigma Zero   mass (MeV/c2)
const double MSm = 1197.449;           // Sigma Minus  mass (MeV/c2)
const double MSp = 1189.37;            // Sigma Plus   mass (MeV/c2)
const int step = 100;



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

//++++++++++++++ class Recoilmom++++++++++++++++//
class Recoilmom
{
public:
Recoilmom();
~Recoilmom();

void MakeTGraph(); 
void SetFrame();
void SetParam();
void DrawGraph();

private:

vector<double> mom1_0deg , mom2_0deg , mom3_0deg ;
vector<double> mom1_10deg, mom2_10deg, mom3_10deg;
vector<double> mom1_20deg, mom2_20deg, mom3_20deg;
vector<double> mom_r1_0deg , mom_r2_0deg , mom_r3_0deg ;
vector<double> mom_r1_10deg, mom_r2_10deg, mom_r3_10deg;
vector<double> mom_r1_20deg, mom_r2_20deg, mom_r3_20deg;
vector<double> volt1, volt2, volt3, volt4;
vector<double> dvdt1, dvdt2, dvdt3, dvdt4;
vector<double> volt_a1, volt_a2, volt_a3, volt_a4;
TGraph *tg_0deg[3] ;
TGraph *tg_10deg[3];
TGraph *tg_20deg[3];
TH2F *h_frame;
TH2F *h_frame_zoom;
TH2F *h_frame_dif;
string ifname[5];
string ofname_pdf;
TLine *tl_zoom_min;
TLine *tl_zoom_max;
TCanvas *c[4];
double mommin;
double mommax;
double Vmin;
double Vmax;
double timemin_zoom;
double timemax_zoom;
unsigned int count1_0deg , count2_0deg , count3_0deg ;
unsigned int count1_10deg, count2_10deg, count3_10deg;
unsigned int count1_20deg, count2_20deg, count3_20deg;
double mom_rec;
double mom_in;

double GetMom_0deg(const double M_in,const double M_detect,const double M_tar, double mom_in);
double GetMom(const double M_in,const double M_detect, double theta_detect,const double M_tar, double mom_in);
void SetMom();
 };

////////////////
Recoilmom::Recoilmom()
:mom1_0deg() , mom2_0deg() , mom3_0deg(),
 mom1_10deg(), mom2_10deg(), mom3_10deg(),
 mom1_20deg(), mom2_20deg(), mom3_20deg(),
 mom_r1_0deg()  , mom_r2_0deg()  , mom_r3_0deg() ,
 mom_r1_10deg() , mom_r2_10deg() , mom_r3_10deg() ,
 mom_r1_20deg() , mom_r2_20deg() , mom_r3_20deg() ,
 volt1(), volt2(), volt3(), volt4(),
 dvdt1(), dvdt2(), dvdt3(), dvdt4()
 {
  }
////////////////
Recoilmom::~Recoilmom(){
}


////////////////

double Recoilmom::GetMom_0deg(const double M_in,const double M_detect, const double M_tar, double mom_in){
double E_in = sqrt(mom_in*mom_in + M_in*M_in); //Energy of beam particle
double s_sq = M_in*M_in + M_tar*M_tar + 2* M_tar *E_in ; //invariant mass square
double mom_cm_sq =  ( (s_sq - ML*ML -M_detect*M_detect  )*(s_sq - ML*ML -M_detect*M_detect  ) -4*ML*ML*M_detect*M_detect  )/(4*s_sq) ;//momentum^2 of lambda or K+ or pi- @CM
double mom_cm;
if(mom_cm_sq>0){
mom_cm = sqrt( mom_cm_sq  );//momentum of lambda or K+ or pi- @CM
 }
else return -99; 
double E_detcm = sqrt(mom_cm*mom_cm + M_detect*M_detect); //energy of detected particle @CM frame
double beta = mom_in/(E_in + M_tar); //beta (Lorentz factor)
double gamma = 1/sqrt(1-beta*beta); //gamma (Lorentz factor)
double mom_detect = beta*gamma*E_detcm + gamma*mom_cm;//momentum of detected particle @Lab frame
//double mom_detect = sqrt( (beta*gamma*E_detcm)*(beta*gamma*E_detcm) +gamma*gamma*mom_cm*mom_cm +2*beta*gamma*gamma*E_detcm*mom_cm );//momentum of detected particle @Lab frame

//cout<<"beta="<<beta<<", mom_detect="<<mom_detect<<endl;

double recoil_mom = mom_in - mom_detect;
if(recoil_mom<0)recoil_mom=-1*recoil_mom;
//cout<<recoil_mom<<endl;
return recoil_mom;
}
////////////////

double Recoilmom::GetMom(const double M_in,const double M_detect, double theta_detect, const double M_tar, double mom_in){
   double theta_det = theta_detect * deg_to_rad; //convert degree to radian
   double E_in = sqrt(mom_in*mom_in + M_in*M_in); //Energy of beam particle
   double s_sq = M_in*M_in + M_tar*M_tar + 2* M_tar *E_in ; //invariant mass square
   double mom_cm_sq =  ( (s_sq - ML*ML -M_detect*M_detect  )*(s_sq - ML*ML -M_detect*M_detect  ) -4*ML*ML*M_detect*M_detect  )/(4*s_sq) ;//momentum^2 of lambda or K+ or pi- @CM
   double mom_cm;
   if(mom_cm_sq>0){
   mom_cm = sqrt( mom_cm_sq  );//momentum of lambda or K+ or pi- @CM
    }
   else return -99; 
   double E_detcm = sqrt(mom_cm*mom_cm + M_detect*M_detect); //energy of detected particle @CM frame
   double beta = mom_in/(E_in + M_tar); //beta (Lorentz factor)
   double gamma = 1/sqrt(1-beta*beta); //gamma (Lorentz factor)
   
   double a0 = (beta*gamma*E_detcm)*(beta*gamma*E_detcm);
   double a1 = a0 + gamma*gamma*mom_cm*mom_cm;
   double b1 = 1+(gamma*gamma-1)*sin(theta_det)*sin(theta_det);
   double c0 = 2*beta*gamma*gamma*E_detcm;
   
   double a = b1*b1;
   double b = c0*c0*sin(theta_det)*sin(theta_det)-2*a1*b1;
   double c = a1*a1-c0*c0*mom_cm*mom_cm;
   
   double detect_mom=sqrt( (-b + sqrt(b*b - 4*a*c)  )/(2*a) );
   
   
   double recoil_mom = sqrt(detect_mom*detect_mom + mom_in*mom_in - 2*detect_mom*mom_in*cos(theta_det) );
   return recoil_mom;
   
   }
   ////////////////
void Recoilmom::SetMom(){
  
  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom_0deg(Mpi,MK, Mp,mom_in);//pi K reaction
  if(mom_rec>-1){mom1_0deg.push_back(mom_in*0.001);mom_r1_0deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_pi="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count1_0deg = mom1_0deg.size();
  
  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom_0deg(MK,Mpi, Mp,mom_in);//K pi reaction
  if(mom_rec>-1){mom2_0deg.push_back(mom_in*0.001);mom_r2_0deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_K="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count2_0deg = mom2_0deg.size();
  
  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom_0deg(0.00,MK, Mp,mom_in);//gamma K reaction
  if(mom_rec>-1){mom3_0deg.push_back(mom_in*0.001);mom_r3_0deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_gamma="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count3_0deg = mom3_0deg.size();
  
//GetMom(const double M_in,const double M_detect, double theta_detect,const double M_tar, double mom_in);
  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(Mpi,MK,10.0, Mp,mom_in);//pi K reaction
  if(mom_rec>-1){mom1_10deg.push_back(mom_in*0.001);mom_r1_10deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_pi="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count1_10deg = mom1_10deg.size();

  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(MK,Mpi, 10.0,Mp,mom_in);//K pi reaction
  if(mom_rec>-1){mom2_10deg.push_back(mom_in*0.001);mom_r2_10deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_K="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count2_10deg = mom2_10deg.size();

  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(0.00,MK, 10.0,Mp,mom_in);//gamma K reaction
  if(mom_rec>-1){mom3_10deg.push_back(mom_in*0.001);mom_r3_10deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_gamma="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count3_10deg = mom3_10deg.size();

  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(Mpi,MK,20.0, Mp,mom_in);//pi K reaction
  if(mom_rec>-1){mom1_20deg.push_back(mom_in*0.001);mom_r1_20deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_pi="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count1_20deg = mom1_20deg.size();
  
  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(MK,Mpi, 20.0,Mp,mom_in);//K pi reaction
  if(mom_rec>-1){mom2_20deg.push_back(mom_in*0.001);mom_r2_20deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_K="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count2_20deg = mom2_20deg.size();

  for(int i=0;i<step;i++){
  mom_in = mommin+i*(mommax-mommin)/step;
  mom_rec = GetMom(0.00,MK, 20.0,Mp,mom_in);//gamma K reaction
  if(mom_rec>-1){mom3_20deg.push_back(mom_in*0.001);mom_r3_20deg.push_back(mom_rec*0.001);}
  //cout<<i<<" mom_gamma="<<mom_in<<", mom_Lambda="<<mom_rec<<endl;
   }
  count3_20deg = mom3_20deg.size();
  
}
////////////////
void Recoilmom::MakeTGraph(){
  SetMom();
  tg_0deg[0]  = new TGraph(count1_0deg, &mom1_0deg[0], &mom_r1_0deg[0]); 
  tg_0deg[1]  = new TGraph(count2_0deg, &mom2_0deg[0], &mom_r2_0deg[0]); 
  tg_0deg[2]  = new TGraph(count3_0deg, &mom3_0deg[0], &mom_r3_0deg[0]); 
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_0deg[0],2,1,4);
  SetGr(tg_0deg[1],2,1,2);
  SetGr(tg_0deg[2],2,1,6);
  
  tg_10deg[0]  = new TGraph(count1_10deg, &mom1_10deg[0], &mom_r1_10deg[0]); 
  tg_10deg[1]  = new TGraph(count2_10deg, &mom2_10deg[0], &mom_r2_10deg[0]); 
  tg_10deg[2]  = new TGraph(count3_10deg, &mom3_10deg[0], &mom_r3_10deg[0]); 
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_10deg[0],2,2,4);
  SetGr(tg_10deg[1],2,2,2);
  SetGr(tg_10deg[2],2,2,6);
  
  tg_20deg[0]  = new TGraph(count1_20deg, &mom1_20deg[0], &mom_r1_20deg[0]); 
  tg_20deg[1]  = new TGraph(count2_20deg, &mom2_20deg[0], &mom_r2_20deg[0]); 
  tg_20deg[2]  = new TGraph(count3_20deg, &mom3_20deg[0], &mom_r3_20deg[0]); 
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_20deg[0],2,10,4);
  SetGr(tg_20deg[1],2,10,2);
  SetGr(tg_20deg[2],2,10,6);
  
  DrawGraph();
}
////////////////

////////////////

////////////////////////////////
void  Recoilmom::SetParam(){
 ofname_pdf = "./pdf/recoil_mom.pdf";
 mommin       =     0;
 mommax       =  2500;
 Vmin         =  -2.0;
 Vmax         =   0.5;
 timemin_zoom =    50;
 timemax_zoom =    70;

h_frame = new TH2F("h_frame","h_frame",10,0,2.5,10,0,1.0);
SetTH2(h_frame , "", "beam momentum [GeV/#it{c}]", "recoil momentum [GeV/#it{c}]");


for(int i=0;i<4;i++){c[i] = new TCanvas(Form("c%d",i+1),Form("c%d",i+1),1200,800);}

 mom1_0deg.clear() , mom2_0deg.clear() , mom3_0deg.clear(),
 mom1_10deg.clear(), mom2_10deg.clear(), mom3_10deg.clear(),
 mom1_20deg.clear(), mom2_20deg.clear(), mom3_20deg.clear(),
 mom_r1_0deg.clear()  , mom_r2_0deg.clear()  , mom_r3_0deg.clear() ,
 mom_r1_10deg.clear() , mom_r2_10deg.clear() , mom_r3_10deg.clear() ,
 mom_r1_20deg.clear() , mom_r2_20deg.clear() , mom_r3_20deg.clear() ;

 }
////////////////////////////////

void Recoilmom::DrawGraph(){
cout<<"start DrawGraph"<<endl;
c[0]->Clear();
h_frame->Draw();
tg_0deg[0]->Draw("samel");
tg_0deg[1]->Draw("samel");
tg_0deg[2]->Draw("samel");
tg_10deg[0]->Draw("samel");
tg_10deg[1]->Draw("samel");
tg_10deg[2]->Draw("samel");
tg_20deg[0]->Draw("samel");
tg_20deg[1]->Draw("samel");
tg_20deg[2]->Draw("samel");

//  ofname_pdf.append(".pdf");
  c[0]->Print(Form("%s[",ofname_pdf.c_str()));
  c[0]->Print(Form("%s" ,ofname_pdf.c_str()));
  c[0]->Print(Form("%s]",ofname_pdf.c_str()));

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
Recoilmom *recoil= new Recoilmom(); 
recoil -> SetParam();
recoil -> MakeTGraph();

    gSystem->Exit(0);
  theApp.Run();
  return 0;

}
