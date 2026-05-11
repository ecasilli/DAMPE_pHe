{
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
//
//
TChain *skim= new TChain("newtree");

// =============================
//---------Orbital data------
// 2016 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2016_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2016_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2016_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2016_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2016_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2016_merged.root");

// 2017 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2017_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2017_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2017_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2017_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2017_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2017_merged.root");

// 2018 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2018_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2018_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2018_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2018_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2018_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2018_merged.root");

// 2019 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2019_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2019_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2019_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2019_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2019_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2019_merged.root");

// 2020 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2020_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2020_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2020_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2020_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2020_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2020_merged.root");

// 2021 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2021_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2021_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2021_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2021_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2021_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2021_merged.root");

// 2022 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2022_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2022_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2022_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2022_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2022_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2022_merged.root");

// 2023 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2023_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2023_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2023_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2023_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2023_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2023_merged.root");

// 2024 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2024_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2024_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2024_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2024_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2024_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2024_merged.root");

// 2025 (Jan-Dec)
skim->Add("/nfs/argo/dampe/PATH/skim_flight_002_010_2025_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_010_025_2025_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_025_050_2025_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_050_100_2025_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_100_500_2025_merged.root");
skim->Add("/nfs/argo/dampe/PATH/skim_flight_500_000_2025_merged.root");


cout<<"Orbital Data Enties:"<<skim -> GetEntries()<<endl;

/*
TCut wcomP = "(2.*3.14159*3.14159*1.38*1.38)";
TCut wcomHe = "(2.*3.14159*3.14159)";
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)";
// 
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


//**************  CHARGE SELECTION 2026 ***********************//
// nota: sono quelli vecchi, da modificare
TString HeMPVf="((7.852)+(1.035)*log10(BGO_EnergyG)+(-0.4836)*log10(BGO_EnergyG)**2+(0.1124)*log10(BGO_EnergyG)**3+(-0.001648)*log10(BGO_EnergyG)**4)";
TString HeWidthf="((0.2505)+(0.02188)*log10(BGO_EnergyG)+(0.1554)*log10(BGO_EnergyG)**2+(-0.06573)*log10(BGO_EnergyG)**3+(0.01144)*log10(BGO_EnergyG)**4)";
TString HeGSigmaf="(0.4789)+(-4.384e-33)*log10(BGO_EnergyG)";	
String HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

TString PMPVf="((2.071)+(-0.3114)*log10(BGO_EnergyG)+(0.4153)*log10(BGO_EnergyG)**2+(-0.1876)*log10(BGO_EnergyG)**3+(0.03373)*log10(BGO_EnergyG)**4)";
TString PWidthf="((0.03606)+(0.09163)*log10(BGO_EnergyG)+(0.05757)*log10(BGO_EnergyG)**2+(-0.04983)*log10(BGO_EnergyG)**3+(0.01135)*log10(BGO_EnergyG)**4)";
TString PGSigmaf="(6.867E-09)+(9.713E-15)*log10(BGO_EnergyG)";		
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



TString p_low= "("+PMPVf+"-"+PSDcharge+")<(2*"+PFSig+")";
TString he_high= "("+PSDcharge+"-"+HeMPVf+")<(6*"+HeFSig+")";

TCut PLow = p_low;
TCut HeHigh= he_high;


TCut Cut_PHe =PLow*HeHigh;

// ---------------  TOTAL CUTS
TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe;

// -----------------------------------------------------------------------------------------------------
// for spectrum

int nbd = 6; // meglio non rebinnare? 
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

cout<<"Now on Orbital Data..."<<endl;


// for spectrum
 
TH1F *h1SelBGO_orb = new TH1F("h1SelBGO_orb", "Selected(E_bgo) orbital", noe, Ebin);//Orbital Data
h1SelBGO_orb->Sumw2();
skim->Draw("(BGO_EnergyG_SatCorr_ML_ions)>>+h1SelBGO_orb",ctot,""); 

TH1F *h1SelBGO_orb_v3 = new TH1F("h1SelBGO_orb_v3", "Selected(E_bgo) orbital", noe, Ebin);//Orbital Data
h1SelBGO_orb_v3->Sumw2();
skim->Draw("(BGO_EnergyG_SatCorr_ML_ions_v3)>>+h1SelBGO_orb_v3",ctot,""); 

cout<<"finished..."<<endl;

TFile *fout1 = new TFile("PHe_skim_Orb120Month_6binperdecade.root","RECREATE");

// for spectrum

fout1->cd();
h1SelBGO_orb->Write();
fout1->Close();

cout<<"End script."<<endl;
};

