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
skim->Add("../../DATA/skim_flight_002_010_2016_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2016_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2016_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2016_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2016_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2016_merged.root");

// 2017 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2017_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2017_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2017_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2017_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2017_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2017_merged.root");

// 2018 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2018_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2018_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2018_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2018_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2018_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2018_merged.root");

// 2019 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2019_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2019_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2019_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2019_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2019_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2019_merged.root");

// 2020 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2020_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2020_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2020_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2020_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2020_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2020_merged.root");

// 2021 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2021_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2021_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2021_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2021_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2021_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2021_merged.root");

// 2022 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2022_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2022_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2022_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2022_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2022_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2022_merged.root");

// 2023 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2023_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2023_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2023_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2023_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2023_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2023_merged.root");

// 2024 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2024_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2024_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2024_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2024_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2024_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2024_merged.root");

// 2025 (Jan-Dec)
skim->Add("../../DATA/skim_flight_002_010_2025_partially_merged.root");
skim->Add("../../DATA/skim_flight_010_025_2025_merged.root");
skim->Add("../../DATA/skim_flight_025_050_2025_merged.root");
skim->Add("../../DATA/skim_flight_050_100_2025_merged.root");
skim->Add("../../DATA/skim_flight_100_500_2025_merged.root");
skim->Add("../../DATA/skim_flight_500_000_2025_merged.root");


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

TCut SpCut = "BGO_xtr > 12."; // nuovo taglio per gli elettroni


// Simple average
TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";

TCut bgo01 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 25.1189) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 39.8107)";
TCut bgo02 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 39.8107) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 63.0957)";
TCut bgo03 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 63.0957) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 100.0)";
TCut bgo04 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 100.0) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 158.489)";
TCut bgo05 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 158.489) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 251.189)";
TCut bgo06 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 251.189) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 398.107)";
TCut bgo07 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 398.107) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 630.957)";
TCut bgo08 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 630.957) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 1000.0)";
TCut bgo09 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 1000.0) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 1584.89)";
TCut bgo10 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 1584.89) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 2511.89)";
TCut bgo11 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 2511.89) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 3981.07)";
TCut bgo12 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 3981.07) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 6309.57)";
TCut bgo13 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 6309.57) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 10000.0)";
TCut bgo14 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 10000.0) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 31622.8)";
TCut bgo15 = "(BGO_EnergyG_SatCorr_ML_ions_v3 > 31622.8) && (BGO_EnergyG_SatCorr_ML_ions_v3 < 100000.0)";

TString name01 = "25 < E_{BGO} < 39";
TString name02 = "39 < E_{BGO} < 63";
TString name03 = "63 < E_{BGO} < 100";
TString name04 = "100 < E_{BGO} < 158";
TString name05 = "158 < E_{BGO} < 251";
TString name06 = "251 < E_{BGO} < 398";
TString name07 = "398 < E_{BGO} < 630)";
TString name08 = "630 < E_{BGO} < 1000)";
TString name09 = "1000 < E_{BGO} < 1584";
TString name10 = "1584 < E_{BGO} < 2511";
TString name11 = "2511 < E_{BGO} < 3981";
TString name12 = "3981 < E_{BGO} < 6309";
TString name13 = "6309 < E_{BGO} < 10000";
TString name14 = "10000 < E_{BGO} < 31622";
TString name15 = "31622 < E_{BGO} < 100000";

TH1F *h01=new TH1F("h01", name01, 240, 0. , 5.); h01->GetXaxis()->SetTitle("PSD charge"); h01->GetYaxis()->SetTitle("events"); h01->SetLineColor(kBlack); h01->SetMarkerColor(kBlack); h01->Sumw2();
TH1F *h02=new TH1F("h02", name02, 240, 0. , 5.); h02->GetXaxis()->SetTitle("PSD charge"); h02->GetYaxis()->SetTitle("events"); h02->SetLineColor(kBlack); h02->SetMarkerColor(kBlack); h02->Sumw2();
TH1F *h03=new TH1F("h03", name03, 240, 0. , 5.); h03->GetXaxis()->SetTitle("PSD charge"); h03->GetYaxis()->SetTitle("events"); h03->SetLineColor(kBlack); h03->SetMarkerColor(kBlack); h03->Sumw2();
TH1F *h04=new TH1F("h04", name04, 240, 0. , 5.); h04->GetXaxis()->SetTitle("PSD charge"); h04->GetYaxis()->SetTitle("events"); h04->SetLineColor(kBlack); h04->SetMarkerColor(kBlack); h04->Sumw2();
TH1F *h05=new TH1F("h05", name05, 240, 0. , 5.); h05->GetXaxis()->SetTitle("PSD charge"); h05->GetYaxis()->SetTitle("events"); h05->SetLineColor(kBlack); h05->SetMarkerColor(kBlack); h05->Sumw2();
TH1F *h06=new TH1F("h06", name06, 240, 0. , 5.); h06->GetXaxis()->SetTitle("PSD charge"); h06->GetYaxis()->SetTitle("events"); h06->SetLineColor(kBlack); h06->SetMarkerColor(kBlack); h06->Sumw2();
TH1F *h07=new TH1F("h07", name07, 240, 0. , 5.); h07->GetXaxis()->SetTitle("PSD charge"); h07->GetYaxis()->SetTitle("events"); h07->SetLineColor(kBlack); h07->SetMarkerColor(kBlack); h07->Sumw2();
TH1F *h08=new TH1F("h08", name08, 240, 0. , 5.); h08->GetXaxis()->SetTitle("PSD charge"); h08->GetYaxis()->SetTitle("events"); h08->SetLineColor(kBlack); h08->SetMarkerColor(kBlack); h08->Sumw2();
TH1F *h09=new TH1F("h09", name09, 240, 0. , 5.); h09->GetXaxis()->SetTitle("PSD charge"); h09->GetYaxis()->SetTitle("events"); h09->SetLineColor(kBlack); h09->SetMarkerColor(kBlack); h09->Sumw2();
TH1F *h10=new TH1F("h10", name10, 240, 0. , 5.); h10->GetXaxis()->SetTitle("PSD charge"); h10->GetYaxis()->SetTitle("events"); h10->SetLineColor(kBlack); h10->SetMarkerColor(kBlack); h10->Sumw2();
TH1F *h11=new TH1F("h11", name11, 240, 0. , 5.); h11->GetXaxis()->SetTitle("PSD charge"); h11->GetYaxis()->SetTitle("events"); h11->SetLineColor(kBlack); h11->SetMarkerColor(kBlack); h11->Sumw2();
TH1F *h12=new TH1F("h12", name12, 240, 0. , 5.); h12->GetXaxis()->SetTitle("PSD charge"); h12->GetYaxis()->SetTitle("events"); h12->SetLineColor(kBlack); h12->SetMarkerColor(kBlack); h12->Sumw2();
TH1F *h13=new TH1F("h13", name13, 240, 0. , 5.); h13->GetXaxis()->SetTitle("PSD charge"); h13->GetYaxis()->SetTitle("events"); h13->SetLineColor(kBlack); h13->SetMarkerColor(kBlack); h13->Sumw2();
TH1F *h14=new TH1F("h14", name14, 240, 0. , 5.); h14->GetXaxis()->SetTitle("PSD charge"); h14->GetYaxis()->SetTitle("events"); h14->SetLineColor(kBlack); h14->SetMarkerColor(kBlack); h14->Sumw2();
TH1F *h15=new TH1F("h15", name15, 240, 0. , 5.); h15->GetXaxis()->SetTitle("PSD charge"); h15->GetYaxis()->SetTitle("events"); h15->SetLineColor(kBlack); h15->SetMarkerColor(kBlack); h15->Sumw2();

/*
//**************  CHARGE SELECTION 2026 ***********************
// nota: sono quelli vecchi, da modificare
TString HeMPVf="((7.852)+(1.035)*log10(BGO_EnergyG)+(-0.4836)*log10(BGO_EnergyG)**2+(0.1124)*log10(BGO_EnergyG)**3+(-0.001648)*log10(BGO_EnergyG)**4)";
TString HeWidthf="((0.2505)+(0.02188)*log10(BGO_EnergyG)+(0.1554)*log10(BGO_EnergyG)**2+(-0.06573)*log10(BGO_EnergyG)**3+(0.01144)*log10(BGO_EnergyG)**4)";
TString HeGSigmaf="(0.4789)+(-4.384e-33)*log10(BGO_EnergyG)";	
String HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

TString PMPVf="((2.071)+(-0.3114)*log10(BGO_EnergyG)+(0.4153)*log10(BGO_EnergyG)**2+(-0.1876)*log10(BGO_EnergyG)**3+(0.03373)*log10(BGO_EnergyG)**4)";
TString PWidthf="((0.03606)+(0.09163)*log10(BGO_EnergyG)+(0.05757)*log10(BGO_EnergyG)**2+(-0.04983)*log10(BGO_EnergyG)**3+(0.01135)*log10(BGO_EnergyG)**4)";
TString PGSigmaf="(6.867E-09)+(9.713E-15)*log10(BGO_EnergyG)";		
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

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
*/
// ---------------  TOTAL CUTS
//TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut*Cut_PHe;
TCut ctot=cut00*cut01*cut02*cut05*cut06*SpCut;

// -----------------------------------------------------------------------------------------------------

cout<<"Now on Orbital Data..."<<endl;
TCanvas *c0=new TCanvas("c0","BGO-PSDglob",1200,900); c0->Divide(3,5);
c0_1->cd();  c0_1->SetTicks();  c0_2->cd();  c0_2->SetTicks();  c0_3->cd();  c0_3->SetTicks();  
c0_4->cd();  c0_4->SetTicks();  c0_5->cd();  c0_5->SetTicks();  c0_6->cd();  c0_6->SetTicks();
c0_7->cd();  c0_7->SetTicks();  c0_8->cd();  c0_8->SetTicks();  c0_9->cd();  c0_9->SetTicks(); 
c0_10->cd();  c0_10->SetTicks();  c0_11->cd();  c0_11->SetTicks();  c0_12->cd();  c0_12->SetTicks(); 
c0_13->cd();  c0_13->SetTicks();  c0_14->cd();  c0_14->SetTicks();  c0_15->cd();  c0_15->SetTicks(); 

c0_1->cd();  skim->Draw(""+PSDcharge+">>h01",ctot*bgo01); cout<<"... h01 ..."<<endl;
c0_2->cd();  skim->Draw(""+PSDcharge+">>h02",ctot*bgo02); cout<<"... h02 ..."<<endl;
c0_3->cd();  skim->Draw(""+PSDcharge+">>h03",ctot*bgo03); cout<<"... h03 ..."<<endl;
c0_4->cd();  skim->Draw(""+PSDcharge+">>h04",ctot*bgo04); cout<<"... h04 ..."<<endl;
c0_5->cd();  skim->Draw(""+PSDcharge+">>h05",ctot*bgo05); cout<<"... h05 ..."<<endl;
c0_6->cd();  skim->Draw(""+PSDcharge+">>h06",ctot*bgo06); cout<<"... h06 ..."<<endl;
c0_7->cd();  skim->Draw(""+PSDcharge+">>h07",ctot*bgo07); cout<<"... h07 ..."<<endl;
c0_8->cd();  skim->Draw(""+PSDcharge+">>h08",ctot*bgo08); cout<<"... h08 ..."<<endl;
c0_9->cd();  skim->Draw(""+PSDcharge+">>h09",ctot*bgo09); cout<<"... h09 ..."<<endl;
c0_10->cd(); skim->Draw(""+PSDcharge+">>h10",ctot*bgo10); cout<<"... h10 ..."<<endl;
c0_11->cd(); skim->Draw(""+PSDcharge+">>h11",ctot*bgo11); cout<<"... h11 ..."<<endl;
c0_12->cd(); skim->Draw(""+PSDcharge+">>h12",ctot*bgo12); cout<<"... h12 ..."<<endl;
c0_13->cd(); skim->Draw(""+PSDcharge+">>h13",ctot*bgo13); cout<<"... h13 ..."<<endl;
c0_14->cd(); skim->Draw(""+PSDcharge+">>h14",ctot*bgo14); cout<<"... h14 ..."<<endl;
c0_15->cd(); skim->Draw(""+PSDcharge+">>h15",ctot*bgo15); cout<<"... h15 ..."<<endl;


cout<<"finished..."<<endl;

TFile *fout1 = new TFile("PHe_charge_Orb120Month_240bins.root","RECREATE");

// for spectrum

fout1->cd();
h01->Write(); // save the histogram
h02->Write(); // save the histogram
h03->Write(); // save the histogram
h04->Write(); // save the histogram
h05->Write(); // save the histogram
h06->Write(); // save the histogram
h07->Write(); // save the histogram
h08->Write(); // save the histogram
h09->Write(); // save the histogram
h10->Write(); // save the histogram
h11->Write(); // save the histogram
h12->Write(); // save the histogram
h13->Write(); // save the histogram
h14->Write(); // save the histogram
h15->Write(); // save the histogram

c0->Write();
fout1->Close();

cout<<"End script."<<endl;
};

