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

/*
// taglio BGO acceptance analitico --> DA AGGIUNGERE?
TCut bgo_valid = "!( (BGO_slopeXZ_analy==0 && BGO_interceptXZ_analy==0) || (BGO_slopeYZ_analy==0 && BGO_interceptYZ_analy==0) )";
TCut bgo_acceptance =
    "fabs(BGO_interceptXZ_analy + 448.*BGO_slopeXZ_analy) < 280. && "
    "fabs(BGO_interceptYZ_analy + 448.*BGO_slopeYZ_analy) < 280. && "
    "fabs(BGO_interceptXZ_analy + 46.*BGO_slopeXZ_analy)  < 280. && "
    "fabs(BGO_interceptYZ_analy + 46.*BGO_slopeYZ_analy)  < 280.";
*/

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

// Simple average
TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";

TCut bgo01 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 25.1189) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 39.8107)";
TCut bgo02 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 39.8107) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 63.0957)";
TCut bgo03 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 63.0957) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 100.0)";
TCut bgo04 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 100.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 158.489)";
TCut bgo05 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 158.489) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 251.189)";
TCut bgo06 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 251.189) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 398.107)";
TCut bgo07 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 398.107) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 630.957)";
TCut bgo08 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 630.957) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 1000.0)";
TCut bgo09 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 1000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 1584.89)";
TCut bgo10 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 1584.89) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 2511.89)";
TCut bgo11 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 2511.89) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 3981.07)";
TCut bgo12 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 3981.07) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 6309.57)";
TCut bgo13 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 6309.57) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 10000.0)";
TCut bgo14 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 10000.0) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 31622.8)";
TCut bgo15 = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3 > 31622.8) && (BGO_EnergyG_QuenchSatCorr_ML_ions_v3 < 100000.0)";

TString name01 = "25 GeV < E_{BGO} < 39 GeV";
TString name02 = "39 GeV < E_{BGO} < 63 GeV";
TString name03 = "63 GeV < E_{BGO} < 100 GeV";
TString name04 = "100 GeV < E_{BGO} < 158 GeV";
TString name05 = "158 GeV < E_{BGO} < 251 GeV";
TString name06 = "251 GeV < E_{BGO} < 398 GeV";
TString name07 = "398 GeV < E_{BGO} < 630) GeV";
TString name08 = "630 GeV < E_{BGO} < 1000) GeV";
TString name09 = "1000 GeV < E_{BGO} < 1584 GeV";
TString name10 = "1584 GeV < E_{BGO} < 2511 GeV";
TString name11 = "2511 GeV < E_{BGO} < 3981 GeV";
TString name12 = "3981 GeV < E_{BGO} < 6309 GeV";
TString name13 = "6309 GeV < E_{BGO} < 10000 GeV";
TString name14 = "10000 GeV < E_{BGO} < 31622 GeV";
TString name15 = "31622 GeV < E_{BGO} < 100000 GeV";

TH1F *h01_p=new TH1F("h01_p", name01, 240, 0., 6.); h01_p->GetXaxis()->SetTitle("PSD charge"); h01_p->GetYaxis()->SetTitle("Normalized counts"); h01_p->SetLineColor(9); h01_p->SetMarkerColor(9); h01_p->Sumw2();
TH1F *h02_p=new TH1F("h02_p", name02, 240, 0., 6.); h02_p->GetXaxis()->SetTitle("PSD charge"); h02_p->GetYaxis()->SetTitle("Normalized counts"); h02_p->SetLineColor(9); h02_p->SetMarkerColor(9); h02_p->Sumw2();
TH1F *h03_p=new TH1F("h03_p", name03, 240, 0., 6.); h03_p->GetXaxis()->SetTitle("PSD charge"); h03_p->GetYaxis()->SetTitle("Normalized counts"); h03_p->SetLineColor(9); h03_p->SetMarkerColor(9); h03_p->Sumw2();
TH1F *h04_p=new TH1F("h04_p", name04, 240, 0., 6.); h04_p->GetXaxis()->SetTitle("PSD charge"); h04_p->GetYaxis()->SetTitle("Normalized counts"); h04_p->SetLineColor(9); h04_p->SetMarkerColor(9); h04_p->Sumw2();
TH1F *h05_p=new TH1F("h05_p", name05, 240, 0., 6.); h05_p->GetXaxis()->SetTitle("PSD charge"); h05_p->GetYaxis()->SetTitle("Normalized counts"); h05_p->SetLineColor(9); h05_p->SetMarkerColor(9); h05_p->Sumw2();
TH1F *h06_p=new TH1F("h06_p", name06, 240, 0., 6.); h06_p->GetXaxis()->SetTitle("PSD charge"); h06_p->GetYaxis()->SetTitle("Normalized counts"); h06_p->SetLineColor(9); h06_p->SetMarkerColor(9); h06_p->Sumw2();
TH1F *h07_p=new TH1F("h07_p", name07, 240, 0., 6.); h07_p->GetXaxis()->SetTitle("PSD charge"); h07_p->GetYaxis()->SetTitle("Normalized counts"); h07_p->SetLineColor(9); h07_p->SetMarkerColor(9); h07_p->Sumw2();
TH1F *h08_p=new TH1F("h08_p", name08, 240, 0., 6.); h08_p->GetXaxis()->SetTitle("PSD charge"); h08_p->GetYaxis()->SetTitle("Normalized counts"); h08_p->SetLineColor(9); h08_p->SetMarkerColor(9); h08_p->Sumw2();
TH1F *h09_p=new TH1F("h09_p", name09, 240, 0., 6.); h09_p->GetXaxis()->SetTitle("PSD charge"); h09_p->GetYaxis()->SetTitle("Normalized counts"); h09_p->SetLineColor(9); h09_p->SetMarkerColor(9); h09_p->Sumw2();
TH1F *h10_p=new TH1F("h10_p", name10, 240, 0., 6.); h10_p->GetXaxis()->SetTitle("PSD charge"); h10_p->GetYaxis()->SetTitle("Normalized counts"); h10_p->SetLineColor(9); h10_p->SetMarkerColor(9); h10_p->Sumw2();
TH1F *h11_p=new TH1F("h11_p", name11, 240, 0., 6.); h11_p->GetXaxis()->SetTitle("PSD charge"); h11_p->GetYaxis()->SetTitle("Normalized counts"); h11_p->SetLineColor(9); h11_p->SetMarkerColor(9); h11_p->Sumw2();
TH1F *h12_p=new TH1F("h12_p", name12, 240, 0., 6.); h12_p->GetXaxis()->SetTitle("PSD charge"); h12_p->GetYaxis()->SetTitle("Normalized counts"); h12_p->SetLineColor(9); h12_p->SetMarkerColor(9); h12_p->Sumw2();
TH1F *h13_p=new TH1F("h13_p", name13, 240, 0., 6.); h13_p->GetXaxis()->SetTitle("PSD charge"); h13_p->GetYaxis()->SetTitle("Normalized counts"); h13_p->SetLineColor(9); h13_p->SetMarkerColor(9); h13_p->Sumw2();
TH1F *h14_p=new TH1F("h14_p", name14, 240, 0., 6.); h14_p->GetXaxis()->SetTitle("PSD charge"); h14_p->GetYaxis()->SetTitle("Normalized counts"); h14_p->SetLineColor(9); h14_p->SetMarkerColor(9); h14_p->Sumw2();
TH1F *h15_p=new TH1F("h15_p", name15, 240, 0., 6.); h15_p->GetXaxis()->SetTitle("PSD charge"); h15_p->GetYaxis()->SetTitle("Normalized counts"); h15_p->SetLineColor(9); h15_p->SetMarkerColor(9); h15_p->Sumw2();


TH1F *h01_he=new TH1F("h01_he", name01, 240, 0., 6.); h01_he->GetXaxis()->SetTitle("PSD charge"); h01_he->GetYaxis()->SetTitle("Normalized counts"); h01_he->SetLineColor(6); h01_he->SetMarkerColor(6); h01_he->Sumw2();
TH1F *h02_he=new TH1F("h02_he", name02, 240, 0., 6.); h02_he->GetXaxis()->SetTitle("PSD charge"); h02_he->GetYaxis()->SetTitle("Normalized counts"); h02_he->SetLineColor(6); h02_he->SetMarkerColor(6); h02_he->Sumw2();
TH1F *h03_he=new TH1F("h03_he", name03, 240, 0., 6.); h03_he->GetXaxis()->SetTitle("PSD charge"); h03_he->GetYaxis()->SetTitle("Normalized counts"); h03_he->SetLineColor(6); h03_he->SetMarkerColor(6); h03_he->Sumw2();
TH1F *h04_he=new TH1F("h04_he", name04, 240, 0., 6.); h04_he->GetXaxis()->SetTitle("PSD charge"); h04_he->GetYaxis()->SetTitle("Normalized counts"); h04_he->SetLineColor(6); h04_he->SetMarkerColor(6); h04_he->Sumw2();
TH1F *h05_he=new TH1F("h05_he", name05, 240, 0., 6.); h05_he->GetXaxis()->SetTitle("PSD charge"); h05_he->GetYaxis()->SetTitle("Normalized counts"); h05_he->SetLineColor(6); h05_he->SetMarkerColor(6); h05_he->Sumw2();
TH1F *h06_he=new TH1F("h06_he", name06, 240, 0., 6.); h06_he->GetXaxis()->SetTitle("PSD charge"); h06_he->GetYaxis()->SetTitle("Normalized counts"); h06_he->SetLineColor(6); h06_he->SetMarkerColor(6); h06_he->Sumw2();
TH1F *h07_he=new TH1F("h07_he", name07, 240, 0., 6.); h07_he->GetXaxis()->SetTitle("PSD charge"); h07_he->GetYaxis()->SetTitle("Normalized counts"); h07_he->SetLineColor(6); h07_he->SetMarkerColor(6); h07_he->Sumw2();
TH1F *h08_he=new TH1F("h08_he", name08, 240, 0., 6.); h08_he->GetXaxis()->SetTitle("PSD charge"); h08_he->GetYaxis()->SetTitle("Normalized counts"); h08_he->SetLineColor(6); h08_he->SetMarkerColor(6); h08_he->Sumw2();
TH1F *h09_he=new TH1F("h09_he", name09, 240, 0., 6.); h09_he->GetXaxis()->SetTitle("PSD charge"); h09_he->GetYaxis()->SetTitle("Normalized counts"); h09_he->SetLineColor(6); h09_he->SetMarkerColor(6); h09_he->Sumw2();
TH1F *h10_he=new TH1F("h10_he", name10, 240, 0., 6.); h10_he->GetXaxis()->SetTitle("PSD charge"); h10_he->GetYaxis()->SetTitle("Normalized counts"); h10_he->SetLineColor(6); h10_he->SetMarkerColor(6); h10_he->Sumw2();
TH1F *h11_he=new TH1F("h11_he", name11, 240, 0., 6.); h11_he->GetXaxis()->SetTitle("PSD charge"); h11_he->GetYaxis()->SetTitle("Normalized counts"); h11_he->SetLineColor(6); h11_he->SetMarkerColor(6); h11_he->Sumw2();
TH1F *h12_he=new TH1F("h12_he", name12, 240, 0., 6.); h12_he->GetXaxis()->SetTitle("PSD charge"); h12_he->GetYaxis()->SetTitle("Normalized counts"); h12_he->SetLineColor(6); h12_he->SetMarkerColor(6); h12_he->Sumw2();
TH1F *h13_he=new TH1F("h13_he", name13, 240, 0., 6.); h13_he->GetXaxis()->SetTitle("PSD charge"); h13_he->GetYaxis()->SetTitle("Normalized counts"); h13_he->SetLineColor(6); h13_he->SetMarkerColor(6); h13_he->Sumw2();
TH1F *h14_he=new TH1F("h14_he", name14, 240, 0., 6.); h14_he->GetXaxis()->SetTitle("PSD charge"); h14_he->GetYaxis()->SetTitle("Normalized counts"); h14_he->SetLineColor(6); h14_he->SetMarkerColor(6); h14_he->Sumw2();
TH1F *h15_he=new TH1F("h15_he", name15, 240, 0., 6.); h15_he->GetXaxis()->SetTitle("PSD charge"); h15_he->GetYaxis()->SetTitle("Normalized counts"); h15_he->SetLineColor(6); h15_he->SetMarkerColor(6); h15_he->Sumw2();



//************** CHARGE SELECTION ***********************

// 25/05/2026
// HELIUM
TString HeMPVf= "(1.88586+(0.185264*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0867115*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0200701*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.000810826*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString HeWidthf = "(-0.0824082+(0.201315*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0989824*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0218456*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00136192*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeGSigmaf = "(0.0396516)";

TString HeMcMPVf= "(1.24779+(1.25798*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.698468*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.169728*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.0123542*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcWidthf = "(-0.207322+(0.351724*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.173554*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0388396*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00242265*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcGSigmaf = "(0.0626831)";

TString HeMcFSig= "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
//TString HeMcFSigCor= "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

// PROTON
TString PMPVf= "(0.993603+(-0.00402039*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0408375*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.0227868*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00494042*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString PWidthf = "(-0.0611545+(0.140328*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0683705*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0143596*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.000527989*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PGSigmaf = "(5.5475e-08)";

TString PMcMPVf= "(0.526726+(0.806415*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.488731*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.129689*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00965891*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcWidthf = "(-0.151343+(0.289688*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.154242*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0383406*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00270484*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcGSigmaf = "(0.000127013)";

TString PMcFSig= "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
//TString PMcFSigCor= "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

//Questo è il valore di PSD energy corretto da usare per i nuovi istogrammi di carica
//TString HeliumCor ="(("+PSDcharge+"-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+")";
//TString ProtonCor ="(("+PSDcharge+"-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+")";
TString HeliumCor = "(PSD_Global_Charge)";
TString ProtonCor = "(PSD_Global_Charge)";
/*
TString p_low_mcCor= "("+PMcMPVf_Cor+"-"+ProtonCor +")<(2*"+PMcFSigCor+")";
TString he_high_mcCor6= "("+HeliumCor+"-"+HeMcMPVf_Cor+")<(6*"+HeMcFSigCor+")";


TString p_low_mc= "("+PMcMPVf+"-"+""+PSDcharge+")<(2*"+PMcFSig+")";
TString he_high_mc6= "("+PSDcharge+"-"+HeMcMPVf+")<(6*"+HeMcFSig+")";

TString p_low= "("+PMPVf+"-"+""+PSDcharge+")<(2*"+PFSig+")*(("+PSDcharge+")>0.)";
TString he_high6= "("+PSDcharge+"-"+HeMPVf+")<(6*"+HeFSig+")";



TCut PLow_MC = p_low_mc;
TCut HeHigh_MC6= he_high_mc6;

TCut PLow_MCCor = p_low_mcCor;
TCut HeHigh_MCCor6= he_high_mcCor6;

TCut PLow = p_low;
TCut HeHigh6= he_high6;

TCut Cut_PHe_MC6 =PLow_MC*HeHigh_MC6;
TCut Cut_PHe_MCCor6 =PLow_MCCor*HeHigh_MCCor6;
TCut Cut_PHe6 =PLow*HeHigh6;
*/
// ---------------  TOTAL CUTS
//TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe_MCCor6;
TCut ctot=cutEt*cut00*cut01*cut02*cut05*cut06*SpCut;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
TCanvas *c0=new TCanvas("c0","BGO-PSDglob",1200,900); c0->Divide(3,5);
c0_1->cd();  c0_1->SetTicks();  c0_2->cd();  c0_2->SetTicks();  c0_3->cd();  c0_3->SetTicks();  
c0_4->cd();  c0_4->SetTicks();  c0_5->cd();  c0_5->SetTicks();  c0_6->cd();  c0_6->SetTicks();
c0_7->cd();  c0_7->SetTicks();  c0_8->cd();  c0_8->SetTicks();  c0_9->cd();  c0_9->SetTicks(); 
c0_10->cd();  c0_10->SetTicks();  c0_11->cd();  c0_11->SetTicks();  c0_12->cd();  c0_12->SetTicks(); 
c0_13->cd();  c0_13->SetTicks();  c0_14->cd();  c0_14->SetTicks();  c0_15->cd();  c0_15->SetTicks(); */

TCanvas* c[15];

for (int i = 0; i < 15; ++i) {

    TString name  = Form("c_%02d", i+1);
    TString title = Form(" %02d", i+1);

    c[i] = new TCanvas(name, title, 800, 600);
}

// -------------- h01
c[0]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h01_p", ctot*wPN1[i]*bgo01, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h01_he",ctot*wHeN1[i]*bgo01,"");} cout<<"... h01 ..."<<endl;

// -------------- h02
c[1]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h02_p", ctot*wPN1[i]*bgo02, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h02_he",ctot*wHeN1[i]*bgo02,"");} cout<<"... h02 ..."<<endl;

// -------------- h03
c[2]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h03_p", ctot*wPN1[i]*bgo03, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h03_he",ctot*wHeN1[i]*bgo03,"");} cout<<"... h03 ..."<<endl;

// -------------- h04
c[3]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h04_p", ctot*wPN1[i]*bgo04, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h04_he",ctot*wHeN1[i]*bgo04,"");} cout<<"... h04 ..."<<endl;

// -------------- h05
c[4]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h05_p", ctot*wPN1[i]*bgo05, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h05_he",ctot*wHeN1[i]*bgo05,"");} cout<<"... h05 ..."<<endl;

// -------------- h06
c[5]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h06_p", ctot*wPN1[i]*bgo06, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h06_he",ctot*wHeN1[i]*bgo06,"");} cout<<"... h06 ..."<<endl;

// -------------- h07
c[6]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h07_p", ctot*wPN1[i]*bgo07, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h07_he",ctot*wHeN1[i]*bgo07,"");} cout<<"... h07 ..."<<endl;

// -------------- h08
c[7]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h08_p", ctot*wPN1[i]*bgo08, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h08_he",ctot*wHeN1[i]*bgo08,"");} cout<<"... h08 ..."<<endl;

// -------------- h09
c[8]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h09_p", ctot*wPN1[i]*bgo09, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h09_he",ctot*wHeN1[i]*bgo09,"");} cout<<"... h09 ..."<<endl;

// -------------- h10
c[9]->cd();for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h10_p", ctot*wPN1[i]*bgo10, "");} 
	 		for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h10_he",ctot*wHeN1[i]*bgo10,"");} cout<<"... h10 ..."<<endl;

// -------------- h11
c[10]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h11_p", ctot*wPN1[i]*bgo11, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h11_he",ctot*wHeN1[i]*bgo11,"");} cout<<"... h11 ..."<<endl;

// -------------- h12
c[11]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h12_p", ctot*wPN1[i]*bgo12, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h12_he",ctot*wHeN1[i]*bgo12,"");} cout<<"... h12 ..."<<endl;

// -------------- h13
c[12]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h13_p", ctot*wPN1[i]*bgo13, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h13_he",ctot*wHeN1[i]*bgo13,"");} cout<<"... h13 ..."<<endl;

// -------------- h14
c[13]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h14_p", ctot*wPN1[i]*bgo14, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h14_he",ctot*wHeN1[i]*bgo14,"");} cout<<"... h14 ..."<<endl;

// -------------- h15
c[14]->cd(); for(int i=0; i<nsetP; i++){  sk_p[i]->Draw( ProtonCor+">>+h15_p", ctot*wPN1[i]*bgo15, "");} 
			for(int i=0; i<nsetHe; i++){ sk_he[i]->Draw(HeliumCor+">>+h15_he",ctot*wHeN1[i]*bgo15,"");} cout<<"... h15 ..."<<endl;
	


TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_PSD_Global_Charge_240bins.root","RECREATE");

fout->cd();

h01_p->Write(); // save the histogram
h02_p->Write(); // save the histogram
h03_p->Write(); // save the histogram
h04_p->Write(); // save the histogram
h05_p->Write(); // save the histogram
h06_p->Write(); // save the histogram
h07_p->Write(); // save the histogram
h08_p->Write(); // save the histogram
h09_p->Write(); // save the histogram
h10_p->Write(); // save the histogram
h11_p->Write(); // save the histogram
h12_p->Write(); // save the histogram
h13_p->Write(); // save the histogram
h14_p->Write(); // save the histogram
h15_p->Write(); // save the histogram

h01_he->Write(); // save the histogram
h02_he->Write(); // save the histogram
h03_he->Write(); // save the histogram
h04_he->Write(); // save the histogram
h05_he->Write(); // save the histogram
h06_he->Write(); // save the histogram
h07_he->Write(); // save the histogram
h08_he->Write(); // save the histogram
h09_he->Write(); // save the histogram
h10_he->Write(); // save the histogram
h11_he->Write(); // save the histogram
h12_he->Write(); // save the histogram
h13_he->Write(); // save the histogram
h14_he->Write(); // save the histogram
h15_he->Write(); // save the histogram

for (int i=0; i<15; i++){
	c[i]->Write();
}


fout->Close();


cout<<"End script."<<endl;
};



