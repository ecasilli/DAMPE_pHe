{
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
//
//

// =============================

const int nsetHe = 7; //number of different energy intervals used to produce the MC data
const int nsetP = 6;

std::vector<TString> filesP = {
        "/nfs/argo/dampe/PATH/Proton_10GeV_100GeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/Proton_100GeV_1TeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/Proton_1TeV_10TeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/Proton_10TeV_100TeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/Proton_100TeV_1PeV_EPOSLHC_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/Proton_1PeV_5PeV_EPOSLHC_FTFP_merged.root"
    };

std::vector<TString> filesHe = {
        "/nfs/argo/dampe/PATH/He4_10GeV_100GeV_FTFP_BGO_Quenching_merged.root",
        "/nfs/argo/dampe/PATH/He4_100GeV_1TeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/He4_1TeV_10TeV_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/He4_10TeV_100TeV_EPOSLHC_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/He4_100TeV_500TeV_EPOSLHC_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/He4_500TeV_1PeV_EPOSLHC_FTFP_merged.root",
        "/nfs/argo/dampe/PATH/He4_1PeV_5PeV_EPOSLHC_FTFP_merged.root"
    };


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

TH1F *h1Ngen_p = new TH1F("h1Ngen_p", "Ngen(Et)", noe, Ebin);
for(int i=0; i<nsetP; i++){

	TFile *p = TFile::Open(filesP[i], "READ");
	TH1D *h = (TH1D*) p->Get("h_energy_truth_6bins_weight_E2e7_cut_0");

	for (int j=1; j<=nbd; j++){
		int globalBin = j + i*nbd;
		h1Ngen_p->SetBinContent(globalBin, h->GetBinContent(globalBin));
	}
	p->Close();
}

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

TH1F *h1Ngen_he = new TH1F("h1Ngen_he", "Ngen(Et)", noe, Ebin);
double splitEnergy = 5e5; // 500 TeV
int overlapBin = h1Ngen_he->FindBin(splitEnergy);
double norm45 = 00/00;
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
}
//
// -----------------------------------------------------------------------------------------

/* // vecchi tagli?
//== here we have a normalization factor (SAME NUMBERS, as before)
TCut wP[nsetP] = {
	"(1./00.)*log(100.)",      //1-100 GeV//      
    "(1./00.)*log(10.)",       //100GeV-1TeV//
    "(1./00.)*log(10.)" ,       //1TeV-10TeV//
	"(1./00.)*log(10.)",       //10TeV-100TeV
	"(1./00.)*log(10.)",      //100TeV-1PeV
	"(1./00.)*log(5.)"      //1PeV-5PeV
};
//
TCut wHe[nsetHe] = {
	"(1./00.)*log(10.)",          //10-100 GeV//
    "(1./00.)*log(10.)",          //100GeV-1TeV//
    "(1./00.)*log(10.)",          //1TeV -10TeV//
    "(1./00.)*log(10.)",         //10TeV-100TeV//
	"(1./00.)*log(5.)",       //100TeV-500TeV//
	"(1./00.)*log(2.)",			//500TeV - 1 PeV 
	"(1./00.)*log(5.)"      //1PeV-5PeV
};

TCut wHeN1[nsetHe];
TCut wPN1[nsetP];

TCut wcomP = "2*3.1415*3.1415*1.38*1.38";                                                 // !!!!!!!!!!!!!!!!!!!!!!
TCut wcomHe = "2*3.14159*3.14159";  
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)"; 

TCut wspectAmsPs =  "(0.4544*(1./45.)**(-2.849)*(MC_EnergyT)**(-1.849)*(1.+(MC_EnergyT/336.)**(5.542))**(0.024))";
TCut wspectAmsHes =  "((0.0948/2.)*(1./(2*45.))**(-2.780)*(MC_EnergyT)**(-1.780)*(1.+(MC_EnergyT/(2*245.))**(4.407))**(0.027))";


// MODIFICARE ???
for (int i=0; i< nsetP; i++) { wPN1[i] = wP[i]*wspectAmsPs; };
for (int i=0; i< nsetHe; i++) { wHeN1[i] = wHe[i]*wspectAmsHes; };
*/

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
TCut cc204s = "(BGO_EnergyG>20.)";

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

TCut SpCut = "BGO_xtr > 16."; // nuovo taglio per gli elettroni
 
TCut cutEt = "MC_EnergyT>20.";



//************** CHARGE SELECTION ***********************//
// nota: sono quelli vecchi, da modificare
TString HeMcMPVf="((12.69)+(-6.543)*log10(BGO_EnergyG)+(3.776)*log10(BGO_EnergyG)**2+(-0.9158)*log10(BGO_EnergyG)**3+(0.1075)*log10(BGO_EnergyG)**4)";
TString HeMcMPVf_Cor="((10.95)+(-4.175)*log10(BGO_EnergyG)+(2.653)*log10(BGO_EnergyG)**2+(-0.6897)*log10(BGO_EnergyG)**3+(0.07261)*log10(BGO_EnergyG)**4)";
TString HeMPVf="((7.852)+(1.035)*log10(BGO_EnergyG)+(-0.4836)*log10(BGO_EnergyG)**2+(0.1124)*log10(BGO_EnergyG)**3+(-0.001648)*log10(BGO_EnergyG)**4)";
TString HeMcWidthf="((2.04)+(-3.504)*log10(BGO_EnergyG)+(2.481)*log10(BGO_EnergyG)**2+(-0.7159)*log10(BGO_EnergyG)**3+(0.08252)*log10(BGO_EnergyG)**4)";
TString HeMcWidthf_Cor="((0.6497)+(-0.9196)*log10(BGO_EnergyG)+(0.7154)*log10(BGO_EnergyG)**2+(-0.1817)*log10(BGO_EnergyG)**3+(0.01912)*log10(BGO_EnergyG)**4)";
TString HeWidthf="((0.2505)+(0.02188)*log10(BGO_EnergyG)+(0.1554)*log10(BGO_EnergyG)**2+(-0.06573)*log10(BGO_EnergyG)**3+(0.01144)*log10(BGO_EnergyG)**4)";


TString HeMcGSigmaf="(0.5807)+(-1.46e-42)*log10(BGO_EnergyG)";
TString HeMcGSigmaf_Cor="(0.5533)+(-2.406e-22)*log10(BGO_EnergyG)";
TString HeGSigmaf="(0.4789)+(-4.384e-33)*log10(BGO_EnergyG)";	


TString HeMcFSig= "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
TString HeMcFSigCor= "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";


TString PMcMPVf="((3.945)+(-3.557)*log10(BGO_EnergyG)+(2.546)*log10(BGO_EnergyG)**2+(-0.7748)*log10(BGO_EnergyG)**3+(0.09923)*log10(BGO_EnergyG)**4)";
TString PMcMPVf_Cor="((1.948)+(-0.1462)*log10(BGO_EnergyG)+(0.3469)*log10(BGO_EnergyG)**2+(-0.1806)*log10(BGO_EnergyG)**3+(0.03454)*log10(BGO_EnergyG)**4)";
TString PMPVf="((2.071)+(-0.3114)*log10(BGO_EnergyG)+(0.4153)*log10(BGO_EnergyG)**2+(-0.1876)*log10(BGO_EnergyG)**3+(0.03373)*log10(BGO_EnergyG)**4)";
TString PMcWidthf="((0.8941)+(-1.455)*log10(BGO_EnergyG)+(1.07)*log10(BGO_EnergyG)**2+(-0.3192)*log10(BGO_EnergyG)**3+(0.03884)*log10(BGO_EnergyG)**4)";
TString PMcWidthf_Cor="((-0.1775)+(0.3985)*log10(BGO_EnergyG)+(-0.09441)*log10(BGO_EnergyG)**2+(-0.02121)*log10(BGO_EnergyG)**3+(0.01006)*log10(BGO_EnergyG)**4)";
TString PWidthf="((0.03606)+(0.09163)*log10(BGO_EnergyG)+(0.05757)*log10(BGO_EnergyG)**2+(-0.04983)*log10(BGO_EnergyG)**3+(0.01135)*log10(BGO_EnergyG)**4)";

TString PMcGSigmaf="(1.056E-09)+(-1.461E-40)*log10(BGO_EnergyG)";
TString PMcGSigmaf_Cor="(8.949E-10)+(-1.461E-40)*log10(BGO_EnergyG)";
TString PGSigmaf="(6.945E-09)+(-1.461E-40)*log10(BGO_EnergyG)";		


//********************** END NEW **********************************//

TString PMcFSig= "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
TString PMcFSigCor= "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";


// Simple average
//TString PSDcharge = "( (((TMath.Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath.Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath.Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath.Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath.Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath.Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath.Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath.Sign(1.,PSD_ChargeX1)+1.)/2.) )";

// Pathweighted average
TString PSDchargePW = "PSD_PathWeighted_Charge";

// Simple average Y-X
//TString PSDchargeY = "( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1)/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.) )";
//TString PSDchargeX = "( ((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1)/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";

// Pathweighted average Y-X
//TString PSDchargeAvX="( ((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_ChargeX0*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_ChargeX1*PSD_psdX_pathlength[1])/((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.*PSD_psdX_pathlength[0]+(TMath::Sign(1.,PSD_ChargeX1)+1.)/2.*PSD_psdX_pathlength[1]))";
//TString PSDchargeAvY="( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_ChargeY0*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_ChargeY1*PSD_psdY_pathlength[1])/((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.*PSD_psdY_pathlength[0]+(TMath::Sign(1.,PSD_ChargeY1)+1.)/2.*PSD_psdY_pathlength[1]))";


//Questo è il valore di PSD energy corretto da usare per i nuovi istogrammi di carica
TString HeliumCor ="(("+PSDcharge+"-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+")";
TString ProtonCor ="(("+PSDcharge+"-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+")";

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

// ---------------  TOTAL CUTS
TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe_MCCor6;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH2F *h2Ntrig_wgt = new TH2F("h2Ntrig_wgt", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2F *h2Ntrig_wgt_v3 = new TH2F("h2Ntrig_wgt_v3", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);

// =============================


for (int i=0; i<nsetHe; i++) {
	for (int j=0; j<noe; j++) {

		sk_he[i]->Draw("(BGO_EnergyG_QuenchSatCorr_ML_ions):(MC_EnergyT)>>+h2Ntrig_wgt",ctot*Weights[j]*Weights2*wHeN[j]*GeoCorr,"goff");
		sk_he[i]->Draw("(BGO_EnergyG_QuenchSatCorr_ML_ions_v3):(MC_EnergyT)>>+h2Ntrig_wgt_v3",ctot*Weights[j]*Weights2*wHeN[j]*GeoCorr,"goff");
	}
    cout<<"He...Cor in:"<<i<<endl;
}

for (int i=0; i<nsetP; i++) {
	for (int j=0; j<noe; j++) {

		sk_p[i]->Draw("(BGO_EnergyG_QuenchSatCorr_ML_ions):(MC_EnergyT)>>+h2Ntrig_wgt",ctot*Weights[j]*Weights2*wPN[j]*GeoCorr,"goff");
		sk_p[i]->Draw("(BGO_EnergyG_QuenchSatCorr_ML_ions_v3):(MC_EnergyT)>>+h2Ntrig_wgt_v3",ctot*Weights[j]*Weights2*wPN[j]*GeoCorr,"goff");
	}
    cout<<"P...Cor in:"<<i<<endl;
}


TFile *fout = new TFile("PHe_MC_p_He_5PeV_6binperdecade.root","RECREATE");

fout->cd();

h2Ntrig_wgt->Write();
h2Ntrig_wgt_v3->Write();

fout->Close();


cout<<"End script."<<endl;
};



