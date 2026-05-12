void Spectra_Indirect(){

    gROOT->Reset();

    double  TotalEnergyDAMPE[29]={0}, FluxDAMPE[29]={0}, FluxEDAMPE[29]={0}, Stat_errDAMPE[29]={0}, Stat_errEDAMPE[29]={0}, 
            EnergyErrDAMPE[29]={0}, Sys_errDAMPE[29]={0}, Sys_errDAMPEFlux[27]={0}, Sys_errEDAMPE_minus[29]={0}, 
            Sys_errEDAMPE_plus[29]={0}, SysHad[27]={0}, SysHadFlux[27]={0};

    double  Energy_ARGOTeV[7]={0.}, Energy_ARGO[7]={0.}, Flux_ARGO_E265[7]={0.}, Stat_err_up_ARGO[7]={0.},
            Stat_err_down_ARGO[7]={0.}, Sys_err_up_ARGO[7]={0.}, Sys_err_down_ARGO[7]={0.},
            Flux_ARGO_E26[7]={0.}, Stat_err_up_ARGO_E26[7]={0.}, Stat_err_down_ARGO_E26[7]={0.},
            Sys_err_up_ARGO_E26[7]={0.}, Sys_err_down_ARGO_E26[7]={0.}, Cont[7]={0.};

    double  Energy_HAWC[14]={0.}, Flux_HAWC[14]={0.}, Stat_err_HAWC[14]={0.}, Sys_err_up_HAWC[14]={0.}, Sys_err_down_HAWC[14]={0.},
            Flux_HAWCE[14]={0.}, Stat_err_HAWCE[14]={0.}, Sys_err_up_HAWCE[14]={0.}, Sys_err_down_HAWCE[14]={0.},
            Energy_HAWC_Err[14]={0.};

    double  Energy_KASCADE_QGSJet_GeV[20]={0.}, Flux_KASCADE_Protons[20]={0.}, Err_up_KASCADE_Protons[20]={0.}, Err_down_KASCADE_Protons[20]={0.}, Flux_KASCADE_Helium[20]={0.},
            Err_up_KASCADE_Helium[20]={0.}, Err_down_KASCADE_Helium[20]={0.}, Err_up_KASCADE_Protons_26[20]={0.}, Err_down_KASCADE_Protons_26[20]={0.}, 
            Err_up_KASCADE_Helium_26[20]={0.}, Err_down_KASCADE_Helium_26[20]={0.}, Flux_KASCADE_Proton_Helium[20]={0.};

    double  Energy_KASCADE_SIBYLL_GeV[20]={0.}, Flux_KASCADE_SIBYLL_Protons[20]={0.}, Err_up_KASCADE_SIBYLL_Protons[20]={0.}, Err_down_KASCADE_SIBYLL_Protons[20]={0.}, 
            Flux_KASCADE_SIBYLL_Helium[20]={0.}, Err_up_KASCADE_SIBYLL_Helium[20]={0.}, Err_down_KASCADE_SIBYLL_Helium[20]={0.},
            Err_up_KASCADE_SIBYLL_Helium_26[20]={0.}, Err_down_KASCADE_SIBYLL_Helium_26[20]={0.},
            Err_up_KASCADE_SIBYLL_Protons_26[20]={0.}, Err_down_KASCADE_SIBYLL_Protons_26[20]={0.},
            Flux_KASCADE_SIBYLL_Proton_Helium[20]={0.}, Err_up_KASCADE_SIBYLL[20]={0.}, Err_down_KASCADE_SIBYLL[20]={0.};

    ifstream in1, in2, in3, in4, in5;
    in1.open("DAMPE_p+He_72M26.dat");
    //in2.open("ARGO2.txt");
    in2.open("ARGO2New.txt");
    in3.open("HAWC2022.txt");
    in4.open("KASCADE_QGSjet01.txt");
    in5.open("KASCADE_SIBYLL.txt");

 
    for(int j=0; j<26; j++){
       in1>>TotalEnergyDAMPE[j]>>FluxDAMPE[j]>>Stat_errDAMPE[j]>>EnergyErrDAMPE[j]>>Sys_errDAMPE[j]>>SysHad[j];
       FluxEDAMPE[j]=FluxDAMPE[j]*(TMath::Power(TotalEnergyDAMPE[j],2.6));
       Stat_errEDAMPE[j]=(Stat_errDAMPE[j]*FluxEDAMPE[j])/100;
       Sys_errDAMPEFlux[j]=(Sys_errDAMPE[j]*FluxEDAMPE[j]);
       SysHadFlux[j]=(SysHad[j]*FluxEDAMPE[j]);
       Sys_errEDAMPE_plus[j]=FluxEDAMPE[j]+(Sys_errDAMPE[j]*FluxEDAMPE[j]);
       Sys_errEDAMPE_minus[j]=FluxEDAMPE[j]-(Sys_errDAMPE[j]*FluxEDAMPE[j]);
    }


    for(int i=0; i<7; i++){
       in2>>Energy_ARGOTeV[i]>>Flux_ARGO_E265[i]>>Stat_err_up_ARGO[i]>>Stat_err_down_ARGO[i]>>Sys_err_up_ARGO[i]>>Sys_err_down_ARGO[i]>>Cont[i];
    
      
        Energy_ARGO[i] = (TMath::Power(10, Energy_ARGOTeV[i]))*1000;
        //Flux_ARGO_E26[i] = ((Flux_ARGO_E265[i])/(10.0E12))*(TMath::Power(Energy_ARGO[i],2.6));  // così è come dovrebbe essere
        Flux_ARGO_E26[i] = ((Flux_ARGO_E265[i])/(10.0E11))*(TMath::Power(Energy_ARGO[i],2.6));     // così è come torna con lo spettro sul paper
        Stat_err_up_ARGO_E26[i] = ((Stat_err_up_ARGO[i])/(10.0E11))*(TMath::Power(Energy_ARGO[i],2.6));
        Stat_err_down_ARGO_E26[i] = Stat_err_up_ARGO_E26[i];

        Sys_err_up_ARGO_E26[i] = Flux_ARGO_E26[i] + (Sys_err_up_ARGO[i]*Flux_ARGO_E26[i])/100.0;

        Cont[i]=(Cont[i]*Flux_ARGO_E26[i])/100.0;
        Sys_err_down_ARGO_E26[i] = (Flux_ARGO_E26[i] - (Sys_err_down_ARGO[i]*Flux_ARGO_E26[i])/100.0) - Cont[i];
                

    /*
       Flux_ARGO_E26[i] = Flux_ARGO_E265[i]/(TMath::Power(Energy_ARGO[i],0.05));
       Stat_err_up_ARGO_E26[i] = (Stat_err_up_ARGO[i]-Flux_ARGO_E265[i])/(TMath::Power(Energy_ARGO[i],0.05));
       Stat_err_down_ARGO_E26[i] = (-Stat_err_down_ARGO[i]+Flux_ARGO_E265[i])/(TMath::Power(Energy_ARGO[i],0.05));
       Sys_err_up_ARGO_E26[i] = (Sys_err_up_ARGO[i])/(TMath::Power(Energy_ARGO[i],0.05));
       Sys_err_down_ARGO_E26[i] = (Sys_err_down_ARGO[i])/(TMath::Power(Energy_ARGO[i],0.05));
    */
    }


    for(int i=0; i<7; i++){
       in3>>Energy_HAWC[i]>>Flux_HAWC[i]>>Stat_err_HAWC[i]>>Sys_err_up_HAWC[i]>>Sys_err_down_HAWC[i];
       Flux_HAWCE[i]=Flux_HAWC[i]*(TMath::Power(Energy_HAWC[i],2.6));
       Stat_err_HAWCE[i]=Stat_err_HAWC[i]*(TMath::Power(Energy_HAWC[i],2.6));
       Sys_err_down_HAWCE[i]=((-Sys_err_down_HAWC[i]*(TMath::Power(Energy_HAWC[i],2.6)))+Flux_HAWCE[i]);
       Sys_err_up_HAWCE[i]=((Sys_err_up_HAWC[i]*(TMath::Power(Energy_HAWC[i],2.6)))+Flux_HAWCE[i]);
       Energy_HAWC_Err[i]=Energy_HAWC[i]*1E-05;
      // Energy_HAWC_log[i] = TMath::Power(10,Energy_HAWC[i]);
    }


    for(int i=0; i<20; i++){
       in4>>Energy_KASCADE_QGSJet_GeV[i]>>Flux_KASCADE_Protons[i]>>Err_down_KASCADE_Protons[i]>>Err_up_KASCADE_Protons[i]>>Flux_KASCADE_Helium[i]>>Err_down_KASCADE_Helium[i]>>Err_up_KASCADE_Helium[i];
       //Energy_KASCADE_QGSJet_GeV[i] = Energy_KASCADE_QGSJet_eV[i]/TMath::Power(10,9);
       Flux_KASCADE_Proton_Helium[i] = (Flux_KASCADE_Protons[i] + Flux_KASCADE_Helium[i])*(TMath::Power(Energy_KASCADE_QGSJet_GeV[i],2.6));
       Err_up_KASCADE_Protons_26[i] = (Err_up_KASCADE_Protons[i])*(TMath::Power(Energy_KASCADE_QGSJet_GeV[i],2.6));
       Err_down_KASCADE_Protons_26[i] = (Err_down_KASCADE_Protons[i])*(TMath::Power(Energy_KASCADE_QGSJet_GeV[i],2.6));
       Err_up_KASCADE_Helium_26[i] = (Err_up_KASCADE_Helium[i])*(TMath::Power(Energy_KASCADE_QGSJet_GeV[i],2.6));
       Err_down_KASCADE_Helium_26[i] = (Err_down_KASCADE_Helium[i])*(TMath::Power(Energy_KASCADE_QGSJet_GeV[i],2.6));
    }


    for(int i=0; i<20; i++){
       in5>>Energy_KASCADE_SIBYLL_GeV[i]>>Flux_KASCADE_SIBYLL_Protons[i]>>Err_up_KASCADE_SIBYLL_Protons[i]>>Err_down_KASCADE_SIBYLL_Protons[i]>>
       Flux_KASCADE_SIBYLL_Helium[i]>>Err_up_KASCADE_SIBYLL_Helium[i]>>Err_down_KASCADE_SIBYLL_Helium[i];
       
       Flux_KASCADE_SIBYLL_Proton_Helium[i] = (Flux_KASCADE_SIBYLL_Protons[i] + Flux_KASCADE_SIBYLL_Helium[i])*(TMath::Power(Energy_KASCADE_SIBYLL_GeV[i],2.6));

       Err_up_KASCADE_SIBYLL_Helium_26[i] = (Err_up_KASCADE_SIBYLL_Helium[i])*(TMath::Power(Energy_KASCADE_SIBYLL_GeV[i],2.6));
       Err_down_KASCADE_SIBYLL_Helium_26[i] = (Err_down_KASCADE_SIBYLL_Helium[i])*(TMath::Power(Energy_KASCADE_SIBYLL_GeV[i],2.6));
       Err_up_KASCADE_SIBYLL_Protons_26[i] = (Err_up_KASCADE_SIBYLL_Protons[i])*(TMath::Power(Energy_KASCADE_SIBYLL_GeV[i],2.6));
       Err_down_KASCADE_SIBYLL_Protons_26[i] = (Err_down_KASCADE_SIBYLL_Protons[i])*(TMath::Power(Energy_KASCADE_SIBYLL_GeV[i],2.6));

      Err_up_KASCADE_SIBYLL[i]=sqrt((Err_up_KASCADE_SIBYLL_Helium_26[i])*(Err_up_KASCADE_SIBYLL_Helium_26[i])+(Err_up_KASCADE_SIBYLL_Protons_26[i])*(Err_up_KASCADE_SIBYLL_Protons_26[i]));
      Err_down_KASCADE_SIBYLL[i]=sqrt((Err_down_KASCADE_SIBYLL_Helium_26[i])*(Err_down_KASCADE_SIBYLL_Helium_26[i])+(Err_down_KASCADE_SIBYLL_Protons_26[i])*(Err_down_KASCADE_SIBYLL_Protons_26[i]));

    }

double F_EAS_TOP = 0, Err_EAS_TOP = 0;
const Double_t E_EAS_TOP=80E03, Err_E_EAS_TOP=0.01; //80 TeV

F_EAS_TOP = (1.80E-09)*(TMath::Power(E_EAS_TOP,2.6));
Err_EAS_TOP = (0.44E-09)*(TMath::Power(E_EAS_TOP,2.6)); //stat+sys

auto mg = new TMultiGraph();

//********************************************************************************************************//
TGraphErrors *EAS_TOP = new TGraphErrors(1,&E_EAS_TOP,&F_EAS_TOP,&Err_E_EAS_TOP, &Err_EAS_TOP);
//TGraph *EAS_TOP = new TGraphErrors(1,&E_EAS_TOP,&F_EAS_TOP);

EAS_TOP->SetMarkerStyle(30);
EAS_TOP->SetMarkerSize(2.6);
EAS_TOP->SetMarkerColor(28);
EAS_TOP->SetLineColor(28);
//EAS_TOP->Draw();

mg->Add(EAS_TOP);

//********************************************************************************************************//

double exl[7]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
double exh[7]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

TGraphAsymmErrors *ARGO = new TGraphAsymmErrors(7, Energy_ARGO, Flux_ARGO_E26, exl, exh, Stat_err_down_ARGO_E26, Stat_err_up_ARGO_E26);
ARGO->SetMarkerStyle(26);
ARGO->SetMarkerSize(1.4);
ARGO->SetMarkerColor(kGray+2);
ARGO->SetLineColor(kGray+2);
mg->Add(ARGO);

TGraph *ARGO_Sys_up = new TGraph(7, Energy_ARGO, Sys_err_up_ARGO_E26);
ARGO_Sys_up->SetMarkerStyle(26);
ARGO_Sys_up->SetMarkerSize(1.4);
ARGO_Sys_up->SetMarkerColor(kGray+2);
ARGO_Sys_up->SetLineColor(kGray+2);
mg->Add(ARGO_Sys_up,"l");

TGraph *ARGO_Sys_down = new TGraph(7, Energy_ARGO, Sys_err_down_ARGO_E26);
ARGO_Sys_down->SetMarkerStyle(26);
ARGO_Sys_down->SetMarkerSize(1.4);
ARGO_Sys_down->SetMarkerColor(kGray+2);
ARGO_Sys_down->SetLineColor(kGray+2);
mg->Add(ARGO_Sys_down,"l");

//******************************************************************************************************//

TGraphErrors *HAWC = new TGraphErrors(7, Energy_HAWC, Flux_HAWCE, Energy_HAWC_Err, Stat_err_HAWCE);
HAWC->SetMarkerStyle(32);
HAWC->SetMarkerSize(1.4);
HAWC->SetMarkerColor(kGray+3);
HAWC->SetLineColor(kGray+3);
mg->Add(HAWC);

TGraph *HAWC_Sys_up = new TGraph(7, Energy_HAWC, Sys_err_up_HAWCE);
HAWC_Sys_up->SetMarkerStyle(32);
HAWC_Sys_up->SetMarkerSize(1.4);
HAWC_Sys_up->SetMarkerColor(kGray+3);
HAWC_Sys_up->SetLineColor(kGray+3);
HAWC_Sys_up->SetLineStyle(2);
mg->Add(HAWC_Sys_up,"l");

TGraph *HAWC_Sys_down = new TGraph(7, Energy_HAWC, Sys_err_down_HAWCE);
HAWC_Sys_down->SetMarkerStyle(32);
HAWC_Sys_down->SetMarkerSize(1.4);
HAWC_Sys_down->SetMarkerColor(kGray+3);
HAWC_Sys_down->SetLineColor(kGray+3);
HAWC_Sys_down->SetLineStyle(2);
mg->Add(HAWC_Sys_down,"l");

//******************************************************************************************************//

double exl_K_QGSJet[10]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
double exh_K_QGSJet[10]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

TGraphAsymmErrors *KASCADE_QGSJet = new TGraphAsymmErrors(20, Energy_KASCADE_QGSJet_GeV, Flux_KASCADE_Proton_Helium, exl_K_QGSJet, exh_K_QGSJet, Err_down_KASCADE_Protons_26, Err_up_KASCADE_Protons_26);
KASCADE_QGSJet->SetMarkerStyle(27);
KASCADE_QGSJet->SetMarkerSize(2.0);
KASCADE_QGSJet->SetMarkerColor(27);
KASCADE_QGSJet->SetLineColor(27);
mg->Add(KASCADE_QGSJet);

//******************************************************************************************************//

double exl_K_SIBYLL[8]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
double exh_K_SIBYLL[8]={0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

TGraphAsymmErrors *KASCADE_SIBYLL = new TGraphAsymmErrors(20, Energy_KASCADE_SIBYLL_GeV, Flux_KASCADE_SIBYLL_Proton_Helium, exl_K_SIBYLL, exh_K_SIBYLL, Err_down_KASCADE_SIBYLL, Err_up_KASCADE_SIBYLL);
//TGraph *KASCADE_SIBYLL = new TGraph(8, Energy_KASCADE_SIBYLL_GeV, Flux_KASCADE_SIBYLL_Proton_Helium);
KASCADE_SIBYLL->SetMarkerStyle(46);
KASCADE_SIBYLL->SetMarkerSize(2.0);
KASCADE_SIBYLL->SetMarkerColor(27);
KASCADE_SIBYLL->SetLineColor(27);
mg->Add(KASCADE_SIBYLL);

//******************************************************************************************************//

TCanvas *C2 = new TCanvas("C2","C2",1000,800);
C2->SetLogx();
C2->SetLogy();
C2->SetTickx();
C2->SetTicky();


TGraphErrors *pHe_flux_DAMPE = new TGraphErrors(24, TotalEnergyDAMPE, FluxEDAMPE, EnergyErrDAMPE, Stat_errEDAMPE);
pHe_flux_DAMPE->SetMarkerStyle(20);
pHe_flux_DAMPE->SetMarkerSize(1.4);
pHe_flux_DAMPE->SetMarkerColor(kRed);
pHe_flux_DAMPE->SetLineColor(kRed);

TGraphErrors *pHe_flux_Sys_Had = new TGraphErrors(24, TotalEnergyDAMPE, FluxEDAMPE, 0, SysHadFlux);
pHe_flux_Sys_Had->SetMarkerStyle(20);
pHe_flux_Sys_Had->SetMarkerSize(0.8);
pHe_flux_Sys_Had->SetMarkerColor(kRed);
pHe_flux_Sys_Had->SetLineColor(kRed);
pHe_flux_Sys_Had->SetFillColor(kAzure-9);  
pHe_flux_Sys_Had->SetFillStyle(1001);  
//pHe_flux_Sys_Had->SetFillStyle(3001);  

TGraphErrors *pHe_flux_Sys_Up = new TGraphErrors(24, TotalEnergyDAMPE, FluxEDAMPE, 0, Sys_errDAMPEFlux);
pHe_flux_Sys_Up->SetMarkerStyle(20);
pHe_flux_Sys_Up->SetMarkerSize(0.8);
pHe_flux_Sys_Up->SetMarkerColor(kRed);
pHe_flux_Sys_Up->SetLineColor(kRed);
pHe_flux_Sys_Up->SetFillColor(kAzure-2);  
pHe_flux_Sys_Up->SetFillStyle(1001);  


pHe_flux_Sys_Had->SetTitle(" ;Primary energy (GeV); E^{2.6} #Phi(E)(GeV^{1.6}s^{-1}m^{-2}sr^{-1})");
//pHe_flux_Sys_Had->Draw("a3");

pHe_flux_Sys_Had->GetYaxis()->SetRangeUser(400,30000);
pHe_flux_Sys_Had->GetXaxis()->SetLimits(30,35E6);

pHe_flux_Sys_Had->GetXaxis()->SetTitleSize(0.042);
pHe_flux_Sys_Had->GetXaxis()->SetLabelSize(0.04);
pHe_flux_Sys_Had->GetYaxis()->SetTitleSize(0.042);
pHe_flux_Sys_Had->GetYaxis()->SetLabelSize(0.04);
pHe_flux_Sys_Had->GetXaxis()->SetTitleOffset(1.07);
pHe_flux_Sys_Had->GetYaxis()->SetTitleOffset(1.07);
pHe_flux_Sys_Had->GetYaxis()->SetDecimals();

pHe_flux_Sys_Had->Draw("a3");
pHe_flux_Sys_Up->Draw("3");    
pHe_flux_DAMPE->Draw("P");
mg->Draw("P");

TLegend* leg = new TLegend(0.1352705,0.6745306,0.3857715,0.8748195); //0.1352705,0.6745306,0.3857715,0.8748195,NULL,"brNDC"
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetHeader("p+He indirect measurements");
leg->AddEntry(EAS_TOP,"EAS-TOP+MACRO (2004)","ep");
leg->AddEntry(ARGO,"ARGO YBJ+WFCT (2015)","ep");
leg->AddEntry(HAWC,"HAWC (2022)","ep");
leg->AddEntry(KASCADE_QGSJet,"KASCADE QGSjet01 (2005)","ep");
leg->AddEntry(KASCADE_SIBYLL,"KASCADE SIBYLL-2.1 (2005)","ep");
leg->AddEntry(pHe_flux_DAMPE,"DAMPE (this work)","ep");
leg->Draw();

C2->Update();

}