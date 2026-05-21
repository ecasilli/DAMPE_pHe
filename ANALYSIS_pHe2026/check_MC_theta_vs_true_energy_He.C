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

const int nsetHe = 7; 

std::vector<TString> filesHe;
filesHe.push_back(basePath + "/HELIUM/He4_10GeV_100GeV_FTFP_BGO_Quenching_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100GeV_1TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1TeV_10TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_10TeV_100TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100TeV_500TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_500TeV_1PeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1PeV_5PeV_EPOSLHC_FTFP_merged.root");


int nbd = 48; // 6 bins per decade
int ndec = 6;

Double_t  arg1 = 1./ Double_t(nbd);
int noe = nbd*ndec;
cout << "Number of energy bin " << noe << endl;
Double_t  e0 = 10.;

Double_t Ebin[noe+1];
Ebin[0] = e0;
for (int i=1; i<noe+1; i++) {
    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
}

//
TChain *sk_he[nsetHe];
for (int i=0; i< nsetHe; i++) { sk_he[i] = new TChain("newtree"); };

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



// -----------------------------------------------------------------------------------------

 // vecchi tagli?
//== here we have a normalization factor (SAME NUMBERS, as before)
/*TCut wP[nsetP] = {
	"(1./1740700000.)*log(10.)",      //10-100 GeV//      
    "(1./519070000.)*log(10.)",       //100GeV-1TeV//
    "(1./201032000.)*log(10.)" ,       //1TeV-10TeV//
	"(1./138744900.)*log(10.)",       //10TeV-100TeV
	"(1./19441600.)*log(10.)",      //100TeV-1PeV
	"(1./5094200.)*log(5.)"      //1PeV-5PeV
};
*/
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
//TCut wPN1[nsetP];

TCut wcomP = "2*3.1415*3.1415*1.38*1.38";                                                 // !!!!!!!!!!!!!!!!!!!!!!
TCut wcomHe = "2*3.14159*3.14159";  
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)"; 

TCut wspectAmsPs =  "(0.4544*(1./45.)**(-2.849)*(MC_EnergyT)**(-1.849)*(1.+(MC_EnergyT/336.)**(5.542))**(0.024))";
TCut wspectAmsHes =  "((0.0948/2.)*(1./(2*45.))**(-2.780)*(MC_EnergyT)**(-1.780)*(1.+(MC_EnergyT/(2*245.))**(4.407))**(0.027))";


TCut wEnergy = "pow(MC_EnergyT, -1.7)";

// MODIFICARE ???
//for (int i=0; i< nsetP; i++) { wPN1[i] = wP[i]*wEnergy; }; //wspectAmsPs; };
for (int i=0; i< nsetHe; i++) { wHeN1[i] = wHe[i]; };//*wEnergy*GeoCorr;  }; //*wspectAmsHes; };


// ------------------------------------------------------------------------------
//    LISTA TAGLI


// taglio BGO acceptance analitico --> DA AGGIUNGERE?
TCut bgo_valid = "!( (BGO_slopeXZ_analy==0 && BGO_interceptXZ_analy==0) || (BGO_slopeYZ_analy==0 && BGO_interceptYZ_analy==0) )";
TCut bgo_acceptance =
    "fabs(BGO_interceptXZ_analy + 448.*BGO_slopeXZ_analy) < 280. && "
    "fabs(BGO_interceptYZ_analy + 448.*BGO_slopeYZ_analy) < 280. && "
    "fabs(BGO_interceptXZ_analy + 46.*BGO_slopeXZ_analy)  < 280. && "
    "fabs(BGO_interceptYZ_analy + 46.*BGO_slopeYZ_analy)  < 280.";


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

TCut cutNtrack = "STK_ntrack>0";

//************** CHARGE SELECTION ***********************//


// Simple average
TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";


// ---------------  TOTAL CUTS
TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut;
//TCut ctot=bgo_valid*bgo_acceptance*cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe_MCCor6;
//TCut ctot=cutNtrack*cut00*cut01*cut02*cut05*cut06*SpCut;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~

TH2F *h2_theta_vs_trueE = new TH2F("h2_theta_vs_trueE", "MC theta vs true energy", noe, Ebin, 100, 0., 100.);

// =============================


for (int i=0; i<nsetHe; i++) {

		sk_he[i]->Draw("MC_theta:MC_EnergyT>>+h2_theta_vs_trueE", cutEt*ctot*wHeN1[i], "goff");

    cout<<"He...Cor in:"<<i<<endl;
}

TFile *fout = new TFile("He_MC_theta_vs_true_energy_2d.root","RECREATE");

fout->cd();

h2_theta_vs_trueE->Write();

fout->Close();


cout<<"End script."<<endl;
};





