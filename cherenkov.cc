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
double GetBeta(const double M_in, double mom_in){
  double beta;
  beta = mom_in/sqrt(mom_in*mom_in +M_in*M_in);
  return beta;
}
////////////////////////////////////////////////////////////////
int main(int argc, char** argv){
gStyle->SetOptStat(0);
gStyle->SetPadGridX(1);
gStyle->SetPadGridY(1);
gStyle->SetPadTopMargin(0.04);
gStyle->SetPadRightMargin(0.04);
gStyle->SetPadBottomMargin(0.17);
  TApplication theApp("App", &argc, argv);
  gROOT->SetBatch(1);

string ofname_pdf;
 ofname_pdf = "./pdf/cherenkov_thre.pdf";
const int step = 1000;
double mommin=0.01;//GeV/c
double mommax=2.;  //GeV/c
TH2F *h_frame = new TH2F("h_frame","h_frame",10,mommin,mommax,10,0,1.1);
SetTH2(h_frame , "", "momentum [GeV/c]", "#beta");
TGraph *tg_e, *tg_pi, *tg_k, *tg_p;

double n[2]={1.05, 1.5};//kussetsu-ritsu

double mom_e[step],mom_pi[step],mom_k[step],mom_p[step];
double beta_e[step],beta_pi[step],beta_k[step],beta_p[step];
  for(int i=0;i<step;i++){
  mom_e[i]  = mommin+i*(mommax-mommin)/step;
  beta_e[i] = GetBeta(Me,1000.*mom_e[i]);
   }
tg_e = new TGraph(step,mom_e,beta_e);
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_e,2,1,4);

  for(int i=0;i<step;i++){
  mom_pi[i]  = mommin+i*(mommax-mommin)/step;
  beta_pi[i] = GetBeta(Mpi,1000.*mom_pi[i]);
   }
tg_pi = new TGraph(step,mom_pi,beta_pi);
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_pi,2,1,6);

  for(int i=0;i<step;i++){
  mom_k[i]  = mommin+i*(mommax-mommin)/step;
  beta_k[i] = GetBeta(MK,1000.*mom_k[i]);
   }
tg_k = new TGraph(step,mom_k,beta_k);
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_k,2,1,3);

  for(int i=0;i<step;i++){
  mom_p[i]  = mommin+i*(mommax-mommin)/step;
  beta_p[i] = GetBeta(Mp,1000.*mom_p[i]);
   }
tg_p = new TGraph(step,mom_p,beta_p);
  //SetGr(TGraph *tg, double linewidth, int linestyle, int linecolor){
  SetGr(tg_p,2,1,2);

TF1 *f[2];
for(int i=0;i<2;i++){
   f[i] =new TF1(Form("f%d",i+1),"[0]",mommin,mommax);
   f[i] ->SetParameter(0, 1./n[i]);
   f[i] ->SetLineColor(1);
   f[i]->SetLineStyle(2);
   f[i]->SetLineWidth(2);
  }


TLegend *leg_beta = new TLegend( 0.7, 0.2, 0.95, 0.50);
leg_beta -> SetBorderSize(1);
leg_beta -> SetFillColor(10);
leg_beta -> SetFillStyle(1001);
leg_beta -> SetTextFont(22);
leg_beta -> AddEntry(tg_e  , "electron", "l");
leg_beta -> AddEntry(tg_pi , "pion"    , "l");
leg_beta -> AddEntry(tg_k  , "Kaon"    , "l");
leg_beta -> AddEntry(tg_p  , "proton"  , "l");

TCanvas *c1 =new TCanvas("c1","canvas",800,600);
//gPad->SetLogx(1);
h_frame ->Draw();
tg_e  ->Draw("same");
tg_pi ->Draw("same");
tg_k  ->Draw("same");
tg_p  ->Draw("same");
f[0]->Draw("same");
f[1]->Draw("same");
leg_beta ->Draw("same"); 

  c1->Print(Form("%s[",ofname_pdf.c_str()));
  c1->Print(Form("%s" ,ofname_pdf.c_str()));
  c1->Print(Form("%s]",ofname_pdf.c_str()));
    gSystem->Exit(0);
  theApp.Run();
  return 0;
}
