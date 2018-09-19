void Esimu(){
gROOT->Reset();
TH1F *hist[6];
//char histname[10],histtitle[12];
for(int k;k<6;k++){
//sprintf(histname,"hist%d",k+1);
//sprintf(histtitle,"SNsimulaton");
hist[k] = new TH1F(Form("hist%d",k+1),"",100,0,1);
//SetTH1(hist[k],"simulation", "beta", "count");
 }
hist[2]->SetFillColor(2);
hist[2]->SetFillStyle(1001);
hist[1]->SetFillColor(3);
hist[1]->SetFillStyle(1001);

int counter = 0;
double pi_mass        = 139.57;
double mom_pi=120;
double Epi  = sqrt(mom_pi*mom_pi+pi_mass*pi_mass);
double beta = mom_pi/Epi;
double Ekin = Epi - pi_mass;
double rEpi, rEkin, rmom_pi;
double sigma =5.; //Ekin resolution [%]
double rbeta;
cout<<beta<<endl;
cout<<Ekin<<endl;
gRandom -> SetSeed( time(NULL) ); //seed set by time
for(int i;i<10000;i++){
rEkin = gRandom -> Gaus(Ekin,Ekin*sigma*0.01);
rEpi = rEkin + pi_mass;
rmom_pi = sqrt(rEpi*rEpi-pi_mass*pi_mass);
rbeta = rmom_pi/rEpi;
hist[0]->Fill(rbeta);
counter ++;
}
cout<<"counter = "<<counter<<endl;


TCanvas *c1 =new TCanvas("c1","canvas",800,600);
hist[0]->Draw();
hist[0]->Fit("gaus");
}
