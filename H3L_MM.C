void H3L_MM(){
gROOT->Reset();
gRandom -> SetSeed( time(NULL) ); //seed set by time
TH1F *hist[6];
//char histname[10],histtitle[12];
for(int k;k<6;k++){
//sprintf(histname,"hist%d",k+1);
//sprintf(histtitle,"SNsimulaton");
hist[k] = new TH1F(Form("hist%d",k+1),"",500,2.95,3.05);
SetTH1(hist[k],"simulation", "mass [GeV/c^{2}]", "count/(2MeV/c^{2})",1,3001,k+2);
 }
//hist[2]->SetFillColor(2);
//hist[2]->SetFillStyle(1001);
//hist[1]->SetFillColor(3);
//hist[1]->SetFillStyle(1001);
double M_He3 = 2.8;//GeV
double Eg    = 1.;//GeV
double cos   = 1.;
double momk    = 0.6;
double MK = 0.493677;// Kaon mass (GeV/c2)
double Ek;

double rEg;   //gamma energy include resolution
double rmomk; //K+ momentum include resolution
double rEk;
double sigma_g = 0.003;//gamma energy resolution[GeV]
double sigma_k = 0.001;//K+ momentum resolution[GeV/c]
double sigma_cos = 0.01;//cos resolution
double rcos;


double mass_x;
double mom_x;

for(int i=0;i<1000;i++){
  rEg   = Eg   + gRandom -> Gaus(0.,sigma_g);
  rmomk = momk + gRandom -> Gaus(0.,sigma_k);
  Ek = sqrt(momk*momk + MK*MK);
  rEk = sqrt(rmomk*rmomk + MK*MK);
  rcos = cos + gRandom -> Gaus(0.,sigma_cos);

  mass_x = sqrt( pow(Eg + M_He3 - Ek,2) - (Eg*Eg + momk*momk -2*Eg*momk*cos) ) ;//
  mom_x  = (Eg*Eg + momk*momk -2*Eg*momk*cos);
hist[0]->Fill(mass_x);

//include tagger resolution
  mass_x = sqrt( pow(rEg + M_He3 - Ek,2) - (rEg*rEg + momk*momk -2*rEg*momk*cos) ) ;//
  mom_x  = (rEg*rEg + momk*momk -2*rEg*momk*cos);
hist[1]->Fill(mass_x);

//include spektrometer resolution
  mass_x = sqrt( pow(rEg + M_He3 - rEk,2) - (rEg*rEg + rmomk*rmomk -2*rEg*rmomk*cos) ) ;//
  mom_x  = (rEg*rEg + rmomk*rmomk -2*rEg*rmomk*cos);
hist[2]->Fill(mass_x);

//include spektrometer angular resolution
  mass_x = sqrt( pow(rEg + M_He3 - rEk,2) - (rEg*rEg + rmomk*rmomk -2*rEg*rmomk*rcos) ) ;//
  mom_x  = (rEg*rEg + rmomk*rmomk -2*rEg*rmomk*rcos);
hist[3]->Fill(mass_x);
// cout<<"momk="<<momk<<"GeV/c, Mx="<<mass_x<<"GeV/c2"<<endl; 
}
TCanvas *c1 =new TCanvas("c1","canvas",800,600);

hist[1]->Draw("");
hist[2]->Draw("same");
//hist[3]->Draw("same");
hist[0]->Draw("same");
/*
int counter = 0;
double S[50],N[200];
double sigma =0.5;
for(int i;i<50;i++){
S[i] = gRandom -> Gaus(0.,sigma);
cout<<S[i]<<endl;
if(fabs(S[i])<log(2)*sigma)counter++;
if(counter>20)break;
hist[0]->Fill(S[i]);
}
cout<<"counter = "<<counter<<endl;
int counter1 = 0;
for(int j;j<200;j++){
N[j] = gRandom -> Uniform(-5*sigma, 5*sigma);
hist[1]->Fill(N[j]);
if(fabs(N[j])<log(2)*sigma)counter1++;
if(counter1>29)break;
}
cout << "counter1 = " <<counter1<<endl;
hist[2]->Add(hist[0],1);
hist[2]->Add(hist[1],1);


*/
}
