
{

int nbd = 6; // 6 bins per decade
int ndec = 6;

float arg1 = 1./float(nbd);
int noe = nbd*ndec;
cout << "Number of energy bin " << noe << endl;
float e0 = 10.;

float Ebin[noe+1];
Ebin[0] = e0;
for (int i=1; i<noe+1; i++) {
    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
}

const int nsetHe = 7; //number of different energy intervals used to produce the MC data
const int nsetP = 6;
// HELIUM
Double_t nwHe[nsetHe] = {
	(1./20380000.)*log(10.),          //10-100 GeV//
    (1./19556000.)*log(10.),          //100GeV-1TeV//
    (1./11346000.)*log(10.),          //1TeV -10TeV//
    (1./17694000.)*log(10.),         //10TeV-100TeV//
	(1./10305100.)*log(5.),       //100TeV-500TeV//
	(1./10281395.)*log(2.),			//500TeV - 1 PeV 
	(1./5064900.)*log(5.)      //1PeV-5PeV
};
// PROTON
Double_t nwP[nsetP] = {
	(1./1740700000.)*log(10.),      //10-100 GeV//      
    (1./519070000.)*log(10.),       //100GeV-1TeV//
    (1./201032000.)*log(10.),       //1TeV-10TeV//
	(1./138744900.)*log(10.),       //10TeV-100TeV
	(1./19441600.)*log(10.),      //100TeV-1PeV
	(1./5094200.)*log(5.)      //1PeV-5PeV
};

cout << "... HELIUM ..." << endl;

TH1D *h1he1 = new TH1D("h1he1", "h1he1", noe, Ebin); //10-100 GeV//
TH1D *h1he2 = new TH1D("h1he2", "h1he2", noe, Ebin); //100GeV-1TeV//
TH1D *h1he3 = new TH1D("h1he3", "h1he3", noe, Ebin); //1TeV-10TeV//
TH1D *h1he4 = new TH1D("h1he4", "h1he4", noe, Ebin); //10TeV-100TeV
TH1D *h1he5 = new TH1D("h1he5", "h1he5", noe, Ebin); //100TeV-500TeV
TH1D *h1he6 = new TH1D("h1he6", "h1he6", noe, Ebin); //500TeV-1PeV
TH1D *h1he7 = new TH1D("h1he7", "h1he7", noe, Ebin); //1PeV-5PeV

h1he1->Sumw2();
h1he2->Sumw2();
h1he3->Sumw2();
h1he4->Sumw2();
h1he5->Sumw2();
h1he6->Sumw2();
h1he7->Sumw2();

temp=20380000./6.;
for(int i1=1;i1<7;i1++){h1he1->SetBinContent(i1,temp/*1.25*/);h1he1->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<1<<endl;

h1he1->Scale(nwHe[0]);

temp=19556000./6.;
for(int i1=7;i1<13;i1++){h1he2->SetBinContent(i1,temp/*1.25*/);h1he2->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<2<<endl;

h1he2->Scale(nwHe[1]);

temp=11346000./6.;
for(int i1=13;i1<19;i1++){h1he3->SetBinContent(i1,temp/*1.25*/);h1he3->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<3<<endl;

h1he3->Scale(nwHe[2]);

temp=17694000./6.;
for(int i1=19;i1<25;i1++){h1he4->SetBinContent(i1,temp/*1.25*/);h1he4->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<4<<endl;

h1he4->Scale(nwHe[3]);

temp=10305100./4.; //33.55; //34
for(int i1=25;i1<29;i1++){h1he5->SetBinContent(i1,temp/*1.25*/);h1he5->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<5<<endl;

h1he5->Scale(nwHe[4]);

temp=10281395./2.; 
for(int i1=29;i1<31;i1++){h1he6->SetBinContent(i1,temp/*1.25*/);h1he6->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<6<<endl;

h1he6->Scale(nwHe[5]);

temp=5064900./4.;
for(int i1=31;i1<35;i1++){h1he7->SetBinContent(i1,temp/*1.25*/);h1he7->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<7<<endl;

h1he7->Scale(nwHe[6]);

//*********************************
cout << "... PROTON ..." << endl;

TH1D *h1p1 = new TH1D("h1p1", "h1p1", noe, Ebin); //10-100 GeV//
TH1D *h1p2 = new TH1D("h1p2", "h1p2", noe, Ebin); //100GeV-1TeV//
TH1D *h1p3 = new TH1D("h1p3", "h1p3", noe, Ebin); //1TeV-10TeV//
TH1D *h1p4 = new TH1D("h1p4", "h1p4", noe, Ebin); //10TeV-100TeV
TH1D *h1p5 = new TH1D("h1p5", "h1p5", noe, Ebin); //100TeV-1PeV
TH1D *h1p6 = new TH1D("h1p6", "h1p6", noe, Ebin); //1PeV-5PeV

h1p1->Sumw2();
h1p2->Sumw2();
h1p3->Sumw2();
h1p4->Sumw2();
h1p5->Sumw2();
h1p6->Sumw2();

temp=1740700000./6.;
for(int i1=1;i1<7;i1++){h1p1->SetBinContent(i1,temp/*1.25*/);h1p1->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<1<<endl;

h1p1->Scale(nwP[0]);

temp=519070000./6.;
for(int i1=7;i1<13;i1++){h1p2->SetBinContent(i1,temp/*1.25*/);h1p2->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<2<<endl;

h1p2->Scale(nwP[1]);

temp=201032000./6.;
for(int i1=13;i1<19;i1++){h1p3->SetBinContent(i1,temp/*1.25*/);h1p3->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<3<<endl;

h1p3->Scale(nwP[2]);

temp=138744900./6.;
for(int i1=19;i1<25;i1++){h1p4->SetBinContent(i1,temp/*1.25*/);h1p4->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<4<<endl;

h1p4->Scale(nwP[3]);

temp=19441600./6.; 
for(int i1=25;i1<31;i1++){h1p5->SetBinContent(i1,temp/*1.25*/);h1p5->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<5<<endl;

h1p5->Scale(nwP[4]);

temp=5094200./4.;
for(int i1=31;i1<35;i1++){h1p6->SetBinContent(i1,temp/*1.25*/);h1p6->SetBinError(i1,sqrt(temp/*1.25*/));}
cout<<6<<endl;

h1p6->Scale(nwP[5]);


//********************************** 
TH1 *h1Ngen=h1he1;
h1Ngen->Add(h1he2);
h1Ngen->Add(h1he3);
h1Ngen->Add(h1he4);
h1Ngen->Add(h1he5);
h1Ngen->Add(h1he6);
h1Ngen->Add(h1he7);

TH1 *h1NgenP=h1p1;
h1NgenP->Add(h1p2);
h1NgenP->Add(h1p3);
h1NgenP->Add(h1p4);
h1NgenP->Add(h1p5);
h1NgenP->Add(h1p6);

h1NgenP->Add(h1Ngen);

TCanvas can;
can.cd();
h1NgenP->Draw();
can.SetLogx();
can.SetGridx();

TFile *fout = new TFile("PHe_MC_FTFP_EPOSLHC_h1Ngen_6binsPerDecade.root","RECREATE");


fout->cd();
h1NgenP->Write();
fout->Close();


}