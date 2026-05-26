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
for (int i=0; i< nsetP; i++) { wPN1[i] = wP[i];}; //*wEnergy; }; //wspectAmsPs; };
for (int i=0; i< nsetHe; i++) { wHeN1[i] = wHe[i]*GeoCorr; }; //*wEnergy; }; //*wspectAmsHes; };

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

//************** CHARGE SELECTION ***********************//


// Simple average
TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";

// HELIUM 05/2026
TString HeMPVf= "(1.81869+(0.290607*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.146116*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0338788*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00192039*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString HeMcMPVf= "(1.24779+(1.25798*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.698468*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.169728*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.0123542*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeWidthf = "(-0.404366+(0.742843*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.422432*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.10251*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00846052*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcWidthf = "(-0.207321+(0.351723*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.173553*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0388396*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00242266*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeGSigmaf = "(0.060102+(-1.43351e-39*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)))";
TString HeMcGSigmaf = "(0.0626831+(-2.19752e-40*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)))";

TString HeMcMPVf_Cor= "(2.00152+(-0.0127222*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0376776*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.0140209*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00265541*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcWidthf_Cor= "(-0.201235+(0.341636*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.160983*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0339152*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00215655*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcGSigmaf_Cor= "(0.0603583)";

TString HeMcFSig= "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
TString HeMcFSigCor= "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";


// PROTON 05/2026
TString PMPVf= "(1.01823+(-0.0471954*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0688881*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.0306054*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.0057345*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString PMcMPVf= "(0.526726+(0.806415*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.488731*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.12969*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00965891*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PWidthf = "(0.0115231+(0.0107847*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0168705*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.00982277*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00196218*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcWidthf = "(-0.151341+(0.289686*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.154241*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0383403*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00270482*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PGSigmaf = "(4.15872e-08+(-1.46051e-40*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)))";
TString PMcGSigmaf = "(0.000112937+(-1.43351e-40*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)))";

TString PMcMPVf_Cor= "(0.980513+(0.0274524*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0158863*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.0142169*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00389876*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcWidthf_Cor= "(-0.0221089+(0.0690865*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.019641*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.000380102*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.000891969*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcGSigmaf_Cor= "(0.000170636+(-2.71193e-08*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)))";

TString PMcFSig= "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
TString PMcFSigCor= "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";



//Questo è il valore di PSD energy corretto da usare per i nuovi istogrammi di carica
TString HeliumCor ="(("+PSDcharge+"-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+")";
TString ProtonCor ="(("+PSDcharge+"-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+")";

TString p_low_mcCor= "("+PMcMPVf_Cor+"-"+ProtonCor +")<(2*"+PMcFSigCor+")";
TString he_high_mcCor6= "("+HeliumCor+"-"+HeMcMPVf_Cor+")<(6*"+HeMcFSigCor+")";

/*
TString p_low_mc= "("+PMcMPVf+"-"+""+PSDcharge+")<(2*"+PMcFSig+")";
TString he_high_mc6= "("+PSDcharge+"-"+HeMcMPVf+")<(6*"+HeMcFSig+")";

TString p_low= "("+PMPVf+"-"+""+PSDcharge+")<(2*"+PFSig+")*(("+PSDcharge+")>0.)";
TString he_high6= "("+PSDcharge+"-"+HeMPVf+")<(6*"+HeFSig+")";



TCut PLow_MC = p_low_mc;
TCut HeHigh_MC6= he_high_mc6;
*/
TCut PLow_MCCor = p_low_mcCor;
TCut HeHigh_MCCor6= he_high_mcCor6;

//TCut PLow = p_low;
//TCut HeHigh6= he_high6;

//TCut Cut_PHe_MC6 =PLow_MC*HeHigh_MC6;
TCut Cut_PHe_MCCor6 =PLow_MCCor*HeHigh_MCCor6;
//TCut Cut_PHe6 =PLow*HeHigh6;

// ---------------  TOTAL CUTS
TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe_MCCor6;
//TCut ctot=cutEt*cut00*cut01*cut02*cut05*cut06*SpCut;


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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~
TH2F *h2_theta_vs_trueE = new TH2F("h2_theta_vs_trueE", "MC theta vs true energy", noe, Ebin, 100, 0., 100.);

// =============================

for (int i=0; i<nsetP; i++) {
	//for (int j=0; j<noe; j++) {

		sk_p[i]->Draw("(MC_theta):(MC_EnergyT)>>+h2_theta_vs_trueE",cutEt*ctot*wPN1[i],"goff");
	//}
    cout<<"P...Cor in:"<<i<<endl;
}

for (int i=0; i<nsetHe; i++) {
	//for (int j=0; j<noe; j++) {

		sk_he[i]->Draw("(MC_theta):(MC_EnergyT)>>+h2_theta_vs_trueE",cutEt*ctot*wHeN1[i],"goff");
	//}
    cout<<"He...Cor in:"<<i<<endl;
}

TFile *fout = new TFile("pHe_MC_theta_vs_true_energy_2d_wChargeSel.root","RECREATE");

fout->cd();

h2_theta_vs_trueE->Write();

fout->Close();


cout<<"End script."<<endl;
};



