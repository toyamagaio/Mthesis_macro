double t2NPE(double thick){
double PDE =0.5;
double N0 = 1600;// /mm
double pixel = 7200.;
double NPE = pixel*( 1-TMath::Exp(-thick*PDE*N0/pixel) );
return NPE;
}
///*************?////
double N2NPE(double N0){
double PDE =0.5;
double pixel = 7200.;
double NPE = pixel*( 1-TMath::Exp(-PDE*N0/pixel) );
return NPE;
}
///*************?////
void NPE(){
gStyle->SetOptStat("");
const int step = 1000;
double Nd[step],thick[step],lin[step],N0[step],NPE[step],paw[step];
double thickmin=0;
double thickmax=20;
double Nmin=0;
double Nmax=20000;
TGraph *tg_lin;
TGraph *tg_nd;
TGraph *tg_npe;

TCanvas *c[2];
for(int i=0;i<2;i++){c[i] =new TCanvas(Form("c%d",i+1),Form("c%d",i+1),800,600);}
TH2F *h_frame = new TH2F("h_frame","h_frame",10,thickmin,thickmax,10,0,7200);
SetTH2(h_frame , "", "thickness [mm]", "N_{det}");

TH2F *h_frame1 = new TH2F("h_frame1","h_frame1",10,0,20000,10,0,7200);
SetTH2(h_frame1 , "", "N_{0}", "N_{det}");

  for(int i=0;i<step;i++){
  thick[i]  = 0+i*(thickmax-thickmin)/step;
  Nd[i]=t2NPE(thick[i]);
  lin[i]=thick[i]*0.5*1600;

  N0[i]  = 0+i*(Nmax-Nmin)/step;
  NPE[i] = N2NPE(N0[i]);
  paw[i]=N0[i]*0.5;
   }
tg_nd  = new TGraph(step,thick,Nd);
tg_lin = new TGraph(step,thick,lin);
tg_npe = new TGraph(step,N0,NPE);
tg_inf = new TGraph(step,N0,paw);
SetGr(tg_nd  ,"","","",2,3,1,1,1,1);
SetGr(tg_lin,"","","",1,3,1,1,1,1);
SetGr(tg_npe,"","","",2,3,1,1,1,1);
SetGr(tg_inf,"","","",1,3,1,1,1,1);
c[0]->Clear();
h_frame ->Draw();
tg_nd   ->Draw("samel");
tg_lin  ->Draw("samel");


c[1]->Clear();
h_frame1 ->Draw();
tg_npe    ->Draw("samel");
tg_inf    ->Draw("samel");


c[0]->Print("pdf/Nd_vs_thick.pdf[");
c[0]->Print("pdf/Nd_vs_thick.pdf");
c[1]->Print("pdf/Nd_vs_thick.pdf");
c[1]->Print("pdf/Nd_vs_thick.pdf]");
}
