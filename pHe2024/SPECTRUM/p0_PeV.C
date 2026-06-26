{

int nbd = 48;

 int ndec = 6;

 float arg1 = 1./float(nbd);
 int noe = nbd*ndec;
 cout << "Number of energy bin " << noe << endl;
 float e0 = 1.;

 float Ebin[noe+1];
 Ebin[0] = e0;
 for (int i=1; i<noe+1; i++) {
    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
 }


const int nsetHe =5;
const int nsetP =4;


//
Double_t nwP[nsetP] = {
		  (1./2029540000.)*log(10.),           //10-100 GeV//      
		  (1./1186785000.)*log(10.),          //100GeV-1TeV//
                  (1./233180000.)*log(10.),             //1TeV-10TeV//
		  (1./22873000.)*log(10.)                 //10TeV-100TeV
			 };
//
Double_t nwHe[nsetHe] = {(1./168736517.)*log(10.),          //10-100 GeV//
                    (1./94071084.)*log(10.),          //100GeV-1TeV//
                    (1./86580000.)*log(10.),          //1TeV -10TeV//
                    (1./148863174.)*log(10.),          //10TeV-100TeV//
		    (1./19747200.)*log(5.)		//100 TeV - 500TeV
		   };

 TRandom3*R3=new TRandom3();


TH1F *h1he1 = new TH1F("h1he1", "h1he1", noe, Ebin);
TH1F *h1he2 = new TH1F("h1he2", "h1he2", noe, Ebin);
TH1F *h1he3 = new TH1F("h1he3", "h1he3", noe, Ebin);
TH1F *h1he4 = new TH1F("h1he4", "h1he4", noe, Ebin);
TH1F *h1he5 = new TH1F("h1he5", "h1he5", noe, Ebin);

TH1F *h1p1 = new TH1F("h1p1", "h1p1", noe, Ebin);
TH1F *h1p2 = new TH1F("h1p2", "h1p2", noe, Ebin);
TH1F *h1p3 = new TH1F("h1p3", "h1p3", noe, Ebin);
TH1F *h1p4 = new TH1F("h1p4", "h1p4", noe, Ebin);


 TF1*func1=new TF1("func1","1/x",1.,1e6);

 
 h1he1->Sumw2();
 h1he2->Sumw2();
 h1he3->Sumw2();
 h1he4->Sumw2();
 h1he5->Sumw2();


 temp=168736517./48.;
 for(int i1=49;i1<97;i1++){h1he1->SetBinContent(i1,temp/*1.25*/);h1he1->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<1<<endl;

 h1he1->Scale(nwHe[0]);

 temp=94071084./48.;
 for(int i1=97;i1<145;i1++){h1he2->SetBinContent(i1,temp/*1.25*/);h1he2->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<2<<endl;

 h1he2->Scale(nwHe[1]);

 temp=86580000./48.;
 for(int i1=145;i1<193;i1++){h1he3->SetBinContent(i1,temp/*1.25*/);h1he3->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<3<<endl;

 h1he3->Scale(nwHe[2]);

 temp=148863174./48.;
 for(int i1=193;i1<241;i1++){h1he4->SetBinContent(i1,temp/*1.25*/);h1he4->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<4<<endl;

 h1he4->Scale(nwHe[3]);

 temp=19747200./33.55; //34
 for(int i1=241;i1<275;i1++){h1he5->SetBinContent(i1,temp/*1.25*/);h1he5->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<5<<endl;

 h1he5->Scale(nwHe[4]);

//*********************************

 h1p1->Sumw2();
 h1p2->Sumw2();
 h1p3->Sumw2();
 h1p4->Sumw2();

 
temp=2029540000./48.;
 for(int i1=49;i1<97;i1++){h1p1->SetBinContent(i1,temp/*1.25*/);h1p1->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<1<<endl;

 h1p1->Scale(nwP[0]);

 temp=1186785000./48.;
 for(int i1=97;i1<145;i1++){h1p2->SetBinContent(i1,temp/*1.25*/);h1p2->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<2<<endl;

 h1p2->Scale(nwP[1]);

 temp=233180000./48.;
 for(int i1=145;i1<193;i1++){h1p3->SetBinContent(i1,temp/*1.25*/);h1p3->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<3<<endl;

 h1p3->Scale(nwP[2]);

 temp=22873000./48.;
 for(int i1=193;i1<241;i1++){h1p4->SetBinContent(i1,temp/*1.25*/);h1p4->SetBinError(i1,sqrt(temp/*1.25*/));}
 cout<<4<<endl;

 h1p4->Scale(nwP[3]);

//********************************** 
  TH1 *h1Ngen=h1he1;

  h1Ngen->Add(h1he2);
  h1Ngen->Add(h1he3);
  h1Ngen->Add(h1he4);
  h1Ngen->Add(h1he5);

  TH1 *h1NgenP=h1p1;

  h1NgenP->Add(h1p2);
  h1NgenP->Add(h1p3);
  h1NgenP->Add(h1p4);

h1NgenP->Add(h1Ngen);

  TCanvas can;
  can.cd();
  h1NgenP->Draw();
  can.SetLogx();



 TFile *fout = new TFile("PHe_MC_QGSP_FLUKA_h1Ngen_48binsPerDecade.root","RECREATE");


 fout->cd();
 h1NgenP->Write();
 fout->Close();

}
