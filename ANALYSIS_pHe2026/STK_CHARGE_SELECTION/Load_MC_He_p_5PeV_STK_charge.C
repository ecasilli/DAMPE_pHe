{
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
//
//

// =============================
// ----- Rilevamento automatico del sito ----
 
TString basePath;
TString hostname = gSystem->HostName();
 
if (hostname.Contains("cnaf")) {
    basePath = "/storage/gpfs_data/dampe/users/SKIM_2026_pHe";
} else if (hostname.Contains("le.infn.it")) {
    basePath = "/nfs/argo/dampe/SKIM_2026_pHe";
} else {
    cout << "ERROR: hostname non riconosciuto: " << hostname << endl;
    cout << "  Aggiungere il path corretto per questo sito." << endl;
    return;
}
 
cout << "Hostname: " << hostname << endl;
cout << "Base path: " << basePath << endl;

// =============================

const int nsetHe = 7; //number of different energy intervals used to produce the MC data
const int nsetP = 6;

std::vector<TString> filesP;
filesP.push_back(basePath + "/PROTON/Proton_10GeV_100GeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_100GeV_1TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_1TeV_10TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_10TeV_100TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_100TeV_1PeV_EPOSLHC_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_1PeV_5PeV_EPOSLHC_FTFP_merged.root");

std::vector<TString> filesHe;
filesHe.push_back(basePath + "/HELIUM/He4_10GeV_100GeV_FTFP_BGO_Quenching_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100GeV_1TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1TeV_10TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_10TeV_100TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100TeV_500TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_500TeV_1PeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1PeV_5PeV_EPOSLHC_FTFP_merged.root");


//
TChain *skim= new TChain("newtree");
TChain *sk_he[nsetHe];
TChain *sk_p[nsetP];
for (int i=0; i< nsetP; i++) { sk_p[i] = new TChain("newtree"); };
for (int i=0; i< nsetHe; i++) { sk_he[i] = new TChain("newtree"); };

//`````````````````````````````````````````````````````````````````````````` ******* PROTON
// -p-
// 10 GeV - 100 GeV
sk_p[0]->Add(filesP[0]);
// 100 GeV - 1 TeV 
sk_p[1]->Add(filesP[1]);
// 1 TeV - 10 TeV 
sk_p[2]->Add(filesP[2]);
// 10 TeV - 100 TeV 
sk_p[3]->Add(filesP[3]);
// 100 TeV - 1 PeV 
sk_p[4]->Add(filesP[4]);
// 100 TeV - 1 PeV 
sk_p[5]->Add(filesP[5]);


 cout<<" "<<endl;
 cout<<"  10 GeV - 100 GeV ...sk_p[0]-> "<<sk_p[0]->GetEntries()<<endl;
 cout<<" 100 GeV - 1 TeV   ...sk_p[1]-> "<<sk_p[1]->GetEntries()<<endl;
 cout<<"   1 TeV - 10 TeV  ...sk_p[2]-> "<<sk_p[2]->GetEntries()<<endl;
 cout<<"  10 TeV - 100 TeV ...sk_p[3]-> "<<sk_p[3]->GetEntries()<<endl;
 cout<<" 100 TeV - 1 PeV   ...sk_p[4]-> "<<sk_p[4]->GetEntries()<<endl;
 cout<<"   1 PeV - 5 PeV   ...sk_p[5]-> "<<sk_p[5]->GetEntries()<<endl;
 cout<<" "<<endl;

// -------------------------------------------------------- Evaluate Ngen:
/*
TH1F *h1Ngen_p = new TH1F("h1Ngen_p", "Ngen(Et)", noe, Ebin);
for(int i=0; i<nsetP; i++){

	TFile *p = TFile::Open(filesP[i], "READ");
	TH1D *h = (TH1D*) p->Get("h_energy_truth_6bins_weight_E2e7_cut_0");

	for (int j=1; j<=nbd; j++){
		int globalBin = j + i*nbd;
		h1Ngen_p->SetBinContent(globalBin, h->GetBinContent(globalBin));
	}
	p->Close();
}*/

//`````````````````````````````````````````````````````````````````````````` ******* HELIUM
// HE 10 GeV - 100 GeV
sk_he[0]->Add(filesHe[0]); 
// HE 100 GeV - 1 TeV
sk_he[1]->Add(filesHe[1]); 
// HE 1 TeV - 10 TeV 
sk_he[2]->Add(filesHe[2]);
// HE 10 TeV - 100 TeV
sk_he[3]->Add(filesHe[3]);
// HE 100 TeV - 500 TeV
sk_he[4]->Add(filesHe[4]);
// HE 500 TeV - 1 PeV 
sk_he[5]->Add(filesHe[5]);
// HE 1 PeV - 5 PeV 
sk_he[6]->Add(filesHe[6]);

cout<<" "<<endl;
cout<<"  10 GeV - 100 GeV ...sk_he[0]-> "<<sk_he[0]->GetEntries()<<endl;
cout<<" 100 GeV - 1 TeV   ...sk_he[1]-> "<<sk_he[1]->GetEntries()<<endl;
cout<<"   1 TeV - 10 TeV  ...sk_he[2]-> "<<sk_he[2]->GetEntries()<<endl;
cout<<"  10 TeV - 100 TeV ...sk_he[3]-> "<<sk_he[3]->GetEntries()<<endl;
cout<<" 100 TeV - 500 TeV ...sk_he[4]-> "<<sk_he[4]->GetEntries()<<endl;
cout<<" 500 TeV - 1 PeV   ...sk_he[5]-> "<<sk_he[5]->GetEntries()<<endl;
cout<<"   1 PeV - 5 PeV   ...sk_he[6]-> "<<sk_he[6]->GetEntries()<<endl;
cout<<" "<<endl;

// -------------------------------------------------------- Evaluate Ngen:
/*
TH1F *h1Ngen_he = new TH1F("h1Ngen_he", "Ngen(Et)", noe, Ebin);
double splitEnergy = 5e5; // 500 TeV
int overlapBin = h1Ngen_he->FindBin(splitEnergy);
double norm45 = 10305100./10281395.;
for(int i=0; i<nsetHe; i++){

	TFile *p = TFile::Open(filesHe[i], "READ");
	TH1D *h = (TH1D*) p->Get("h_energy_truth_6bins_weight_E2e7_cut_0");

	for (int j=1; j<=nbd; j++){
		int globalBin = j + i*nbd;
		double val = h->GetBinContent(globalBin);
		// add scaled contribution from second sample in overlap bin
        if (i == 5 && globalBin == overlapBin)
            val *= norm45;

        h1Ngen_he->AddBinContent(globalBin, val);
	}
	p->Close();
}*/
//
// -----------------------------------------------------------------------------------------


//== here we have a normalization factor (SAME NUMBERS, as before)
TCut wP[nsetP] = {
	"(1./1740700000.)*log(10.)",      //10-100 GeV//      
    "(1./519070000.)*log(10.)",       //100GeV-1TeV//
    "(1./201032000.)*log(10.)" ,       //1TeV-10TeV//
	"(1./138744900.)*log(10.)",       //10TeV-100TeV
	"(1./19441600.)*log(10.)",      //100TeV-1PeV
	"(1./5094200.)*log(5.)"      //1PeV-5PeV
};
//
TCut wHe[nsetHe] = {
	"(1./20380000.)*log(10.)",          //10-100 GeV//
    "(1./19556000.)*log(10.)",          //100GeV-1TeV//
    "(1./11346000.)*log(10.)",          //1TeV -10TeV//
    "(1./17694000.)*log(10.)",         //10TeV-100TeV//
	"(1./10305100.)*log(5.)",       //100TeV-500TeV//
	"(1./10281395.)*log(2.)",			//500TeV - 1 PeV 
	"(1./5064900.)*log(5.)"      //1PeV-5PeV
};

TCut wHeN1[nsetHe];
TCut wPN1[nsetP];

TCut wcomP = "2*3.1415*3.1415*1.38*1.38";                                                 // !!!!!!!!!!!!!!!!!!!!!!
TCut wcomHe = "2*3.14159*3.14159";  
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)"; 

TCut wspectAmsPs =  "(0.4544*(1./45.)**(-2.849)*(MC_EnergyT)**(-1.849)*(1.+(MC_EnergyT/336.)**(5.542))**(0.024))";
TCut wspectAmsHes =  "((0.0948/2.)*(1./(2*45.))**(-2.780)*(MC_EnergyT)**(-1.780)*(1.+(MC_EnergyT/(2*245.))**(4.407))**(0.027))";

TCut wEnergy = "(MC_EnergyT)**(-1.7)";

// MODIFICARE ???
for (int i=0; i< nsetP; i++) { wPN1[i] = wP[i]*wEnergy; }; //wspectAmsPs; };
for (int i=0; i< nsetHe; i++) { wHeN1[i] = wHe[i]*wEnergy*GeoCorr; }; //*wspectAmsHes; };

/*
TCut weightsInterval[noe];

int index = 2.7;
TString weights[noe];
TString PriE0[noe];
TString PriE1[noe];
TString TrueEnergy = "MC_EnergyT";
TString wNgen_p[noe];
TString wNgen_he[noe];

TCut Weights[noe];
TCut wPN[noe], wHeN[noe];

for (int i=0; i<noe; i++){
	weightsInterval[i]=Form("MC_EnergyT>=%f && MC_EnergyT<%f",Ebin[i],Ebin[i+1]);
	PriE0[i]=Form("%f",Ebin[i]);
	PriE1[i]=Form("%f",Ebin[i+1]);

	weights[i]="(1.7)*TMath::Log("+PriE1[i]+"/"+PriE0[i]+")/(TMath::Power("+PriE0[i]+",-(1.7))-TMath::Power("+PriE1[i]+",-(1.7)))";
	Weights[i] = (weights[i]*weightsInterval[i]);

	wNgen_p[i] = Form("%e", 1./h1Ngen_p->GetBinContent(i+1));
	wPN[i] = wNgen_p[i];
	wNgen_he[i] = Form("%e", 1./h1Ngen_he->GetBinContent(i+1));
	wHeN[i] = wNgen_he[i];

}

TString weights2 = "TMath::Power(MC_EnergyT,-(1.7))";
TCut Weights2 = weights2;

*/

// ------------------------------------------------------------------------------
//    LISTA TAGLI


// taglio BGO acceptance analitico --> DA AGGIUNGERE?
TCut bgo_valid = "!( (BGO_slopeXZ_analy==0 && BGO_interceptXZ_analy==0) || (BGO_slopeYZ_analy==0 && BGO_interceptYZ_analy==0) )";
TCut bgo_acceptance =
    "fabs(BGO_interceptXZ_analy + 448.*BGO_slopeXZ_analy) < 280. && "
    "fabs(BGO_interceptYZ_analy + 448.*BGO_slopeYZ_analy) < 280. && "
    "fabs(BGO_interceptXZ_analy + 46.*BGO_slopeXZ_analy)  < 280. && "
    "fabs(BGO_interceptYZ_analy + 46.*BGO_slopeYZ_analy)  < 280.";
TCut bgo_acc = bgo_valid*bgo_acceptance;

TCut Trig_HEP="BGO_HET>0.";
TCut cc204s = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3>20.)";

TCut cut00  = cc204s*Trig_HEP;
TCut cut01 = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)"; // un segnale su entrambe i due piani di PSD
TCut cut02 = "BGO_EneLay[0]+BGO_EneLay[1]<BGO_EneLay[2]+BGO_EneLay[3]";

/*
//TCut cut03="Etrack01>0.0 && Etrack23>0.0";                         // un segnale su entrambe i due piani di PSD
TCut cut03ups="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03dws="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03=cut03ups*cut03dws; 

//TCut cut04="(fabs(sqrt(Etrack01)-sqrt(epsdmax01/facslope))<5. && fabs(sqrt(Etrack23)-sqrt(epsdmax23/facslope))<5.)";
TCut cut04ups="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))-sqrt(PSD_psdY_Emax_layer*STK_theta_correction))<5.";
TCut cut04dws="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))-sqrt(PSD_psdX_Emax_layer*STK_theta_correction))<5.";
TCut cut04=cut04ups*cut04dws;
*/

//TCut cut05="(fabs(cbgomax01-cbgostk01)<30.) && (fabs(cbgomax02-cbgostk02)<30.)";
TCut cut05 ="TMath::Abs(BGO_cbgomax[0]-BGO_cbgostk[0])<30. && TMath::Abs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.";
TCut cut06 = "fabs(STKtrack_to_PSD_topY)<400. && fabs(STKtrack_to_PSD_topX)<400."; 

TCut SpCut = "BGO_xtr > 12."; // nuovo taglio per gli elettroni
 
TCut cutEt = "MC_EnergyT>20.";

TString STK_charge = "((((TMath::Sign(1.,STK_chargeY[0])+1.)/2.*STK_chargeY[0]+(TMath::Sign(1.,STK_chargeX[0])+1.)/2.*STK_chargeX[0])/((TMath::Sign(1.,STK_chargeY[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX[0])+1.)/2.)))";

// Simple average
//TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";

TCut bgo01 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 10.0)     && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 100.0)";
TCut bgo02 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 100.0)    && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 1000.0)";
TCut bgo03 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 1000.0)   && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 10000.0)";
TCut bgo04 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 10000.0)  && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 100000.0)";
TCut bgo05 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 100000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 1000000.0)";

TString name01 = "10 GeV < E_{BGO} < 100 GeV";
TString name02 = "100 GeV < E_{BGO} < 1 TeV";
TString name03 = "1 TeV < E_{BGO} < 10 TeV";
TString name04 = "10 TeV < E_{BGO} < 100 TeV";
TString name05 = "100 TeV < E_{BGO} < 1 PeV";

TH1F *h01_p=new TH1F("h01_p", name01, 200, 0., 1000.); h01_p->GetXaxis()->SetTitle("STK signal (ADC counts)"); h01_p->GetYaxis()->SetTitle("Normalized counts"); h01_p->SetLineColor(9); h01_p->SetMarkerColor(9); h01_p->Sumw2();
TH1F *h02_p=new TH1F("h02_p", name02, 200, 0., 1000.); h02_p->GetXaxis()->SetTitle("STK signal (ADC counts)"); h02_p->GetYaxis()->SetTitle("Normalized counts"); h02_p->SetLineColor(9); h02_p->SetMarkerColor(9); h02_p->Sumw2();
TH1F *h03_p=new TH1F("h03_p", name03, 200, 0., 1000.); h03_p->GetXaxis()->SetTitle("STK signal (ADC counts)"); h03_p->GetYaxis()->SetTitle("Normalized counts"); h03_p->SetLineColor(9); h03_p->SetMarkerColor(9); h03_p->Sumw2();
TH1F *h04_p=new TH1F("h04_p", name04, 200, 0., 1000.); h04_p->GetXaxis()->SetTitle("STK signal (ADC counts)"); h04_p->GetYaxis()->SetTitle("Normalized counts"); h04_p->SetLineColor(9); h04_p->SetMarkerColor(9); h04_p->Sumw2();
TH1F *h05_p=new TH1F("h05_p", name05, 200, 0., 1000.); h05_p->GetXaxis()->SetTitle("STK signal (ADC counts)"); h05_p->GetYaxis()->SetTitle("Normalized counts"); h05_p->SetLineColor(9); h05_p->SetMarkerColor(9); h05_p->Sumw2();

TH1F *h01_he=new TH1F("h01_he", name01, 200, 0., 1000.); h01_he->GetXaxis()->SetTitle("STK signal (ADC counts)"); h01_he->GetYaxis()->SetTitle("Normalized counts"); h01_he->SetLineColor(6); h01_he->SetMarkerColor(6); h01_he->Sumw2();
TH1F *h02_he=new TH1F("h02_he", name02, 200, 0., 1000.); h02_he->GetXaxis()->SetTitle("STK signal (ADC counts)"); h02_he->GetYaxis()->SetTitle("Normalized counts"); h02_he->SetLineColor(6); h02_he->SetMarkerColor(6); h02_he->Sumw2();
TH1F *h03_he=new TH1F("h03_he", name03, 200, 0., 1000.); h03_he->GetXaxis()->SetTitle("STK signal (ADC counts)"); h03_he->GetYaxis()->SetTitle("Normalized counts"); h03_he->SetLineColor(6); h03_he->SetMarkerColor(6); h03_he->Sumw2();
TH1F *h04_he=new TH1F("h04_he", name04, 200, 0., 1000.); h04_he->GetXaxis()->SetTitle("STK signal (ADC counts)"); h04_he->GetYaxis()->SetTitle("Normalized counts"); h04_he->SetLineColor(6); h04_he->SetMarkerColor(6); h04_he->Sumw2();
TH1F *h05_he=new TH1F("h05_he", name05, 200, 0., 1000.); h05_he->GetXaxis()->SetTitle("STK signal (ADC counts)"); h05_he->GetYaxis()->SetTitle("Normalized counts"); h05_he->SetLineColor(6); h05_he->SetMarkerColor(6); h05_he->Sumw2();

TH1F *h01=new TH1F("h01", name01, 200, 0., 1000.); h01->GetXaxis()->SetTitle("STK signal (ADC counts)"); h01->GetYaxis()->SetTitle("Normalized counts"); h01->SetLineColor(9); h01->SetMarkerColor(9); h01->Sumw2();
TH1F *h02=new TH1F("h02", name02, 200, 0., 1000.); h02->GetXaxis()->SetTitle("STK signal (ADC counts)"); h02->GetYaxis()->SetTitle("Normalized counts"); h02->SetLineColor(9); h02->SetMarkerColor(9); h02->Sumw2();
TH1F *h03=new TH1F("h03", name03, 200, 0., 1000.); h03->GetXaxis()->SetTitle("STK signal (ADC counts)"); h03->GetYaxis()->SetTitle("Normalized counts"); h03->SetLineColor(9); h03->SetMarkerColor(9); h03->Sumw2();
TH1F *h04=new TH1F("h04", name04, 200, 0., 1000.); h04->GetXaxis()->SetTitle("STK signal (ADC counts)"); h04->GetYaxis()->SetTitle("Normalized counts"); h04->SetLineColor(9); h04->SetMarkerColor(9); h04->Sumw2();
TH1F *h05=new TH1F("h05", name05, 200, 0., 1000.); h05->GetXaxis()->SetTitle("STK signal (ADC counts)"); h05->GetYaxis()->SetTitle("Normalized counts"); h05->SetLineColor(9); h05->SetMarkerColor(9); h05->Sumw2();

// ---------------  TOTAL CUTS
//TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe_MCCor6;
TCut ctot=cutEt*bgo_acc*cut00*cut01*cut02*cut05*cut06*SpCut;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~

TCanvas* c[10];

for (int i = 0; i < 10; ++i) {

    TString name  = Form("c_%02d", i+1);
    TString title = Form(" %02d", i+1);

    c[i] = new TCanvas(name, title, 800, 600);
}

// -------------- h01
c[0]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h01_p", ctot*wPN1[i]*bgo01, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h01_he",ctot*wHeN1[i]*bgo01,"");} cout<<"... h01 ..."<<endl;

// -------------- h02
c[1]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h02_p", ctot*wPN1[i]*bgo02, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h02_he",ctot*wHeN1[i]*bgo02,"");} cout<<"... h02 ..."<<endl;

// -------------- h03
c[2]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h03_p", ctot*wPN1[i]*bgo03, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h03_he",ctot*wHeN1[i]*bgo03,"");} cout<<"... h03 ..."<<endl;

// -------------- h04
c[3]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h04_p", ctot*wPN1[i]*bgo04, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h04_he",ctot*wHeN1[i]*bgo04,"");} cout<<"... h04 ..."<<endl;

// -------------- h05
c[4]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h05_p", ctot*wPN1[i]*bgo05, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h05_he",ctot*wHeN1[i]*bgo05,"");} cout<<"... h05 ..."<<endl;

// -------------- h01
c[5]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h01", ctot*wPN1[i]*bgo01, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h01",ctot*wHeN1[i]*bgo01,"");} cout<<"... h01 ..."<<endl;

// -------------- h02
c[6]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h02", ctot*wPN1[i]*bgo02, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h02",ctot*wHeN1[i]*bgo02,"");} cout<<"... h02 ..."<<endl;

// -------------- h03
c[7]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h03", ctot*wPN1[i]*bgo03, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h03",ctot*wHeN1[i]*bgo03,"");} cout<<"... h03 ..."<<endl;

// -------------- h04
c[8]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h04", ctot*wPN1[i]*bgo04, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h04",ctot*wHeN1[i]*bgo04,"");} cout<<"... h04 ..."<<endl;

// -------------- h05
c[9]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( STK_charge+">>+h05", ctot*wPN1[i]*bgo05, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(STK_charge+">>+h05",ctot*wHeN1[i]*bgo05,"");} cout<<"... h05 ..."<<endl;
	


TFile *fout = new TFile("PHe_MC_p_He_5PeV_STK_charge_200bins.root","RECREATE");

fout->cd();

h01_p->Write(); // save the histogram
h02_p->Write(); // save the histogram
h03_p->Write(); // save the histogram
h04_p->Write(); // save the histogram
h05_p->Write(); // save the histogram

h01_he->Write(); // save the histogram
h02_he->Write(); // save the histogram
h03_he->Write(); // save the histogram
h04_he->Write(); // save the histogram
h05_he->Write(); // save the histogram

h01->Write(); // save the histogram
h02->Write(); // save the histogram
h03->Write(); // save the histogram
h04->Write(); // save the histogram
h05->Write(); // save the histogram

for (int i=0; i<10; i++){
	c[i]->Write();
}


fout->Close();


cout<<"End script."<<endl;
};



