void Spectra_Direct(){

    gROOT->Reset();

    double  TotalEnergyDAMPE[27]={0}, FluxDAMPE[27]={0}, FluxEDAMPE[27]={0}, Stat_errDAMPE[27]={0}, Stat_errEDAMPE[27]={0}, 
            EnergyErrDAMPE[27]={0}, Sys_errDAMPE[27]={0}, Sys_errDAMPEFlux[27]={0}, Sys_errEDAMPE_minus[27]={0}, Sys_errEDAMPE_plus[27]={0}, SysHad[27]={0}, SysHadFlux[27]={0};

    double  TotalEnergyDAMPE_2026[36]={0}, FluxDAMPE_2026[36]={0}, FluxEDAMPE_2026[36]={0}, Stat_errDAMPE_2026[36]={0}, Stat_errEDAMPE_2026[36]={0}, 
            EnergyErrDAMPE_2026[36]={0};

    double  TotalEnergyDAMPE_2026_1[36]={0}, FluxDAMPE_2026_1[36]={0}, FluxEDAMPE_2026_1[36]={0}, Stat_errDAMPE_2026_1[36]={0}, Stat_errEDAMPE_2026_1[36]={0}, 
            EnergyErrDAMPE_2026_1[36]={0};

    double  TotalEnergyDAMPE_2026_2[36]={0}, FluxDAMPE_2026_2[36]={0}, FluxEDAMPE_2026_2[36]={0}, Stat_errDAMPE_2026_2[36]={0}, Stat_errEDAMPE_2026_2[36]={0}, 
            EnergyErrDAMPE_2026_2[36]={0};

    double  TotalEnergyDAMPE_2026_3[36]={0}, FluxDAMPE_2026_3[36]={0}, FluxEDAMPE_2026_3[36]={0}, Stat_errDAMPE_2026_3[36]={0}, Stat_errEDAMPE_2026_3[36]={0}, 
            EnergyErrDAMPE_2026_3[36]={0};

    double EnergyATIC[15]={0.}, FluxProtonATIC[15]={0.}, ErrProtonATIC[15]={0.}, FluxHeliumATIC[15]={0.}, ErrHeliumATIC[15]={0.}, ErrEnergyATIC[15]={0.},
           Flux_Proton_Helium_ATIC_27[15]={0.}, Err_Proton_Helium_ATIC_27[15]={0.};

    double EnergyCREAM[10]={0.}, FluxProtonCREAM[10]={0.}, ErrProtonCREAM[10]={0.}, FluxHeliumCREAM[10]={0.}, ErrHeliumCREAM[10]={0.}, ErrEnergyCREAM[10]={0.},
           Flux_Proton_Helium_CREAM_27[10]={0.}, Err_Proton_Helium_CREAM_27[10]={0.};

    double EnergyNucleon[10]={0.}, FluxProtonNucleon[10]={0.}, ErrProtonNucleon[10]={0.}, FluxHeliumNucleon[10]={0.}, ErrHeliumNucleon[10]={0.}, ErrEnergyNucleon[10]={0.},
           Flux_Proton_Helium_Nucleon_27[10]={0.}, Err_Proton_Helium_Nucleon_27[10]={0.};


    ifstream in1, in2, in3, in4, in5, in6, in7, in8;
    in1.open("../pHe2024/SPECTRUM/Direct/DAMPE_p+He_72M26.dat");
    in2.open("../pHe2024/SPECTRUM/Direct/ATIC2_ProtonStat_HeliumStat.dat");
    in3.open("../pHe2024/SPECTRUM/Direct/CREAM_Proton_Helium.dat");
    in4.open("../pHe2024/SPECTRUM/Direct/Nucleon_p_He.dat");

    in5.open("flux_spectrum_pHe_2026_MLionsv3.dat");
    in6.open("flux_spectrum_pHe_2026_MLionsv3_ntrackCut.dat");
    in7.open("flux_spectrum_pHe_2026_MLionsv3_BGOAccCut_v0.dat");

    in8.open("flux_spectrum_pHe_2026_MLionsv3_smooth.dat");


    for(int j=0; j<26; j++){
       in1>>TotalEnergyDAMPE[j]>>FluxDAMPE[j]>>Stat_errDAMPE[j]>>EnergyErrDAMPE[j]>>Sys_errDAMPE[j]>>SysHad[j];
       FluxEDAMPE[j]=FluxDAMPE[j]*(TMath::Power(TotalEnergyDAMPE[j],2.7));
       Stat_errEDAMPE[j]=(Stat_errDAMPE[j]*FluxEDAMPE[j])/100;
       Sys_errDAMPEFlux[j]=(Sys_errDAMPE[j]*FluxEDAMPE[j]);
       SysHadFlux[j]=(SysHad[j]*FluxEDAMPE[j]);
       Sys_errEDAMPE_plus[j]=FluxEDAMPE[j]+(Sys_errDAMPE[j]*FluxEDAMPE[j]);
       Sys_errEDAMPE_minus[j]=FluxEDAMPE[j]-(Sys_errDAMPE[j]*FluxEDAMPE[j]);
 //      cout<<TotalEnergyDAMPE[j]<<"\t"<<Stat_errEDAMPE[j]/(TMath::Power(TotalEnergyDAMPE[j],2.7))<<endl;
 //       cout<<TotalEnergyDAMPE[j]<<"\t"<<Sys_errDAMPEFlux[j]/(TMath::Power(TotalEnergyDAMPE[j],2.7))<<endl;
       cout<<TotalEnergyDAMPE[j]<<"\t"<<SysHadFlux[j]/(TMath::Power(TotalEnergyDAMPE[j],2.7))<<endl;
    }

    for(int i=0; i<15; i++){
            in2>>EnergyATIC[i]>>FluxProtonATIC[i]>>ErrProtonATIC[i]>>FluxHeliumATIC[i]>>ErrHeliumATIC[i]>>ErrEnergyATIC[i];
            Flux_Proton_Helium_ATIC_27[i]=(FluxProtonATIC[i] + FluxHeliumATIC[i])*(TMath::Power(EnergyATIC[i],2.7));
            Err_Proton_Helium_ATIC_27[i]=(sqrt((ErrProtonATIC[i])*(ErrProtonATIC[i])+(ErrHeliumATIC[i])*(ErrHeliumATIC[i])))*(TMath::Power(EnergyATIC[i],2.7));
    }

    for(int i=0; i<10; i++){
            in3>>EnergyCREAM[i]>>FluxProtonCREAM[i]>>ErrProtonCREAM[i]>>FluxHeliumCREAM[i]>>ErrHeliumCREAM[i]>>ErrEnergyCREAM[i];
            Flux_Proton_Helium_CREAM_27[i]=(FluxProtonCREAM[i] + FluxHeliumCREAM[i])*(TMath::Power(EnergyCREAM[i],2.7));
            Err_Proton_Helium_CREAM_27[i]=(sqrt((ErrProtonCREAM[i])*(ErrProtonCREAM[i])+(ErrHeliumCREAM[i])*(ErrHeliumCREAM[i])))*(TMath::Power(EnergyCREAM[i],2.7));
    }

    for(int i=0; i<10; i++){
            in4>>EnergyNucleon[i]>>FluxProtonNucleon[i]>>ErrProtonNucleon[i]>>FluxHeliumNucleon[i]>>ErrHeliumNucleon[i]>>ErrEnergyNucleon[i];
            Flux_Proton_Helium_Nucleon_27[i]=(FluxProtonNucleon[i] + FluxHeliumNucleon[i])*(TMath::Power(EnergyNucleon[i],2.7));
            Err_Proton_Helium_Nucleon_27[i]=(sqrt((ErrProtonNucleon[i])*(ErrProtonNucleon[i])+(ErrHeliumNucleon[i])*(ErrHeliumNucleon[i])))*(TMath::Power(EnergyNucleon[i],2.7));
    }

    // outfile << "# E  lowE  upE  Flux  StatErr\n";
    for(int j=0; j<36; j++){
        in5>>TotalEnergyDAMPE_2026[j]>>FluxDAMPE_2026[j]>>Stat_errDAMPE_2026[j]>>EnergyErrDAMPE_2026[j];
        FluxEDAMPE_2026[j]=FluxDAMPE_2026[j]*(TMath::Power(TotalEnergyDAMPE_2026[j],2.7));
        Stat_errEDAMPE_2026[j]=Stat_errDAMPE_2026[j]*(TMath::Power(TotalEnergyDAMPE_2026[j],2.7));
    }

    // outfile << "# E  lowE  upE  Flux  StatErr\n";
    for(int j=0; j<36; j++){
        in6>>TotalEnergyDAMPE_2026_1[j]>>FluxDAMPE_2026_1[j]>>Stat_errDAMPE_2026_1[j]>>EnergyErrDAMPE_2026_1[j];
        FluxEDAMPE_2026_1[j]=FluxDAMPE_2026_1[j]*(TMath::Power(TotalEnergyDAMPE_2026_1[j],2.7));
        Stat_errEDAMPE_2026_1[j]=Stat_errDAMPE_2026_1[j]*(TMath::Power(TotalEnergyDAMPE_2026_1[j],2.7));
    }

    // outfile << "# E  lowE  upE  Flux  StatErr\n";
    for(int j=0; j<36; j++){
        in7>>TotalEnergyDAMPE_2026_2[j]>>FluxDAMPE_2026_2[j]>>Stat_errDAMPE_2026_2[j]>>EnergyErrDAMPE_2026_2[j];
        FluxEDAMPE_2026_2[j]=FluxDAMPE_2026_2[j]*(TMath::Power(TotalEnergyDAMPE_2026_2[j],2.7));
        Stat_errEDAMPE_2026_2[j]=Stat_errDAMPE_2026_2[j]*(TMath::Power(TotalEnergyDAMPE_2026_2[j],2.7));
    }

    // outfile << "# E  lowE  upE  Flux  StatErr\n";
    for(int j=0; j<36; j++){
        in8>>TotalEnergyDAMPE_2026_3[j]>>FluxDAMPE_2026_3[j]>>Stat_errDAMPE_2026_3[j]>>EnergyErrDAMPE_2026_3[j];
        FluxEDAMPE_2026_3[j]=FluxDAMPE_2026_3[j]*(TMath::Power(TotalEnergyDAMPE_2026_3[j],2.7));
        Stat_errEDAMPE_2026_3[j]=Stat_errDAMPE_2026_3[j]*(TMath::Power(TotalEnergyDAMPE_2026_3[j],2.7));
    }

TGaxis::SetMaxDigits(3);

auto mg = new TMultiGraph();


//********************************************************************************************************//

TGraphErrors *ATIC = new TGraphErrors(15, EnergyATIC, Flux_Proton_Helium_ATIC_27, ErrEnergyATIC, Err_Proton_Helium_ATIC_27);
ATIC->SetMarkerStyle(21);
ATIC->SetMarkerSize(1.2);
ATIC->SetMarkerColor(kGray+2);
ATIC->SetLineColor(kGray+2);
mg->Add(ATIC);

//*******************************************************************************************************//

TGraphErrors *CREAM = new TGraphErrors(10, EnergyCREAM, Flux_Proton_Helium_CREAM_27, ErrEnergyCREAM, Err_Proton_Helium_CREAM_27);
CREAM->SetMarkerStyle(23);
CREAM->SetMarkerSize(1.2);
CREAM->SetMarkerColor(kGray+3);
CREAM->SetLineColor(kGray+3);
mg->Add(CREAM);

//*******************************************************************************************************//

TGraphErrors *Nucleon = new TGraphErrors(10, EnergyNucleon, Flux_Proton_Helium_Nucleon_27, ErrEnergyNucleon, Err_Proton_Helium_Nucleon_27);
Nucleon->SetMarkerStyle(33);
Nucleon->SetMarkerSize(1.7);
Nucleon->SetMarkerColor(15);
Nucleon->SetLineColor(15);
mg->Add(Nucleon);

//*******************************************************************************************************//
TCanvas *C2 = new TCanvas("C2","C2",1000,800);
C2->SetLogx();
//C2->SetLogy();
C2->SetTickx();
C2->SetTicky();

TPad *pad = new TPad("pad","pad",0,0,1,1);
pad->SetLogx();
pad->Draw();
pad->cd();

TH1F *frame = pad->DrawFrame(2e1, 0, 1e6, 50E03);

frame->GetXaxis()->SetTitleSize(0.042);
frame->GetXaxis()->SetLabelSize(0.04);

frame->GetYaxis()->SetTitleSize(0.042);
frame->GetYaxis()->SetLabelSize(0.045);

frame->GetXaxis()->SetTitleOffset(1.25);
frame->GetYaxis()->SetTitleOffset(0.7);

frame->GetYaxis()->SetDecimals();

frame->SetTitle(";Primary energy (GeV);E^{2.7} #Phi(E) (GeV^{1.7} s^{-1} m^{-2} sr^{-1})");

TGraphErrors *pHe_flux_DAMPE = new TGraphErrors(26, TotalEnergyDAMPE, FluxEDAMPE, EnergyErrDAMPE, Stat_errEDAMPE);
pHe_flux_DAMPE->SetMarkerStyle(20);
pHe_flux_DAMPE->SetMarkerSize(1.4);
pHe_flux_DAMPE->SetMarkerColor(kRed);
pHe_flux_DAMPE->SetLineColor(kRed); 
//mg->Add(pHe_flux_DAMPE);     

TGraphErrors *pHe_flux_Sys_Had = new TGraphErrors(24, TotalEnergyDAMPE, FluxEDAMPE, 0, SysHadFlux);
pHe_flux_Sys_Had->SetMarkerStyle(20);
pHe_flux_Sys_Had->SetMarkerSize(0.8);
pHe_flux_Sys_Had->SetMarkerColor(kRed);
pHe_flux_Sys_Had->SetLineColor(kRed);
pHe_flux_Sys_Had->SetFillColor(kAzure-9);  
pHe_flux_Sys_Had->SetFillStyle(1001);  

TGraphErrors *pHe_flux_Sys_Up = new TGraphErrors(24, TotalEnergyDAMPE, FluxEDAMPE, 0, Sys_errDAMPEFlux);
pHe_flux_Sys_Up->SetMarkerStyle(20);
pHe_flux_Sys_Up->SetMarkerSize(0.8);
pHe_flux_Sys_Up->SetMarkerColor(kRed);
pHe_flux_Sys_Up->SetLineColor(kRed);
pHe_flux_Sys_Up->SetFillColor(kAzure-2);  
pHe_flux_Sys_Up->SetFillStyle(1001);  

pHe_flux_Sys_Had->SetTitle(" ;Primary energy (GeV); E^{2.7} #Phi(E)(GeV^{1.7}s^{-1}m^{-2}sr^{-1})");

   pHe_flux_Sys_Had->GetXaxis()->SetRangeUser(20,1.0E6);
   pHe_flux_Sys_Had->GetYaxis()->SetRangeUser(0,60E03);
   pHe_flux_Sys_Had->GetXaxis()->SetTitleSize(0.042);
   pHe_flux_Sys_Had->GetXaxis()->SetLabelSize(0.04);
   pHe_flux_Sys_Had->GetYaxis()->SetTitleSize(0.042);
   pHe_flux_Sys_Had->GetYaxis()->SetLabelSize(0.04);
   pHe_flux_Sys_Had->GetXaxis()->SetTitleOffset(1.07);
   pHe_flux_Sys_Had->GetYaxis()->SetTitleOffset(1.07);
   pHe_flux_Sys_Had->GetYaxis()->SetDecimals();


TGraphErrors *pHe_flux_DAMPE_2026 = new TGraphErrors(35, TotalEnergyDAMPE_2026, FluxEDAMPE_2026, EnergyErrDAMPE_2026, Stat_errEDAMPE_2026);
pHe_flux_DAMPE_2026->SetMarkerStyle(20);
pHe_flux_DAMPE_2026->SetMarkerSize(1.4);
pHe_flux_DAMPE_2026->SetMarkerColor(kBlack);
pHe_flux_DAMPE_2026->SetLineColor(kBlack); 
mg->Add(pHe_flux_DAMPE_2026);  

TGraphErrors *pHe_flux_DAMPE_2026_1 = new TGraphErrors(35, TotalEnergyDAMPE_2026_1, FluxEDAMPE_2026_1, EnergyErrDAMPE_2026_1, Stat_errEDAMPE_2026_1);
pHe_flux_DAMPE_2026_1->SetMarkerStyle(20);
pHe_flux_DAMPE_2026_1->SetMarkerSize(1.4);
pHe_flux_DAMPE_2026_1->SetMarkerColor(kMagenta);
pHe_flux_DAMPE_2026_1->SetLineColor(kMagenta); 
//mg->Add(pHe_flux_DAMPE_2026_1);  

TGraphErrors *pHe_flux_DAMPE_2026_2 = new TGraphErrors(35, TotalEnergyDAMPE_2026_2, FluxEDAMPE_2026_2, EnergyErrDAMPE_2026_2, Stat_errEDAMPE_2026_2);
pHe_flux_DAMPE_2026_2->SetMarkerStyle(20);
pHe_flux_DAMPE_2026_2->SetMarkerSize(1.4);
pHe_flux_DAMPE_2026_2->SetMarkerColor(kBlue);
pHe_flux_DAMPE_2026_2->SetLineColor(kBlue); 
//mg->Add(pHe_flux_DAMPE_2026_2);   

TGraphErrors *pHe_flux_DAMPE_2026_3 = new TGraphErrors(35, TotalEnergyDAMPE_2026_3, FluxEDAMPE_2026_3, EnergyErrDAMPE_2026_3, Stat_errEDAMPE_2026_3);
pHe_flux_DAMPE_2026_3->SetMarkerStyle(20);
pHe_flux_DAMPE_2026_3->SetMarkerSize(1.4);
pHe_flux_DAMPE_2026_3->SetMarkerColor(kBlue);
pHe_flux_DAMPE_2026_3->SetLineColor(kBlue); 
//mg->Add(pHe_flux_DAMPE_2026_3);  

pHe_flux_Sys_Had->Draw("3");
pHe_flux_Sys_Up->Draw("3");
mg->Draw("P");
pHe_flux_DAMPE->Draw("P");
//He_flux_DAMPE->Draw("P");

/*
TGraph *pHe_flux_Sys_Up = new TGraph(23, TotalEnergyDAMPE, Sys_errEDAMPE_plus);
pHe_flux_Sys_Up->SetMarkerStyle(20);
pHe_flux_Sys_Up->SetMarkerSize(0.8);
pHe_flux_Sys_Up->SetMarkerColor(kRed);
pHe_flux_Sys_Up->SetLineColor(kRed);
mg->Add(pHe_flux_Sys_Up, "l");

TGraph *pHe_flux_Sys_Down = new TGraph(23, TotalEnergyDAMPE, Sys_errEDAMPE_minus);
pHe_flux_Sys_Down->SetMarkerStyle(20);
pHe_flux_Sys_Down->SetMarkerSize(0.8);
pHe_flux_Sys_Down->SetMarkerColor(kRed);
pHe_flux_Sys_Down->SetLineColor(kRed);
mg->Add(pHe_flux_Sys_Down, "l");
*/

//TCanvas *C2 = new TCanvas("C2","C2",1000,800);
//C2->SetLogx();
//C2->SetTickx();
//C2->SetTicky();
//mg->SetTitle(" ;Primary energy (GeV); E^{2.7} #Phi(E)(GeV^{1.7}s^{-1}m^{-2}sr^{-1})");
//mg->Draw("AP");


TLegend* leg = new TLegend(0.1352705,0.6745306,0.3857715,0.8748195); //0.1352705,0.6745306,0.3857715,0.8748195,NULL,"brNDC"
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetHeader("p+He direct measurements");
leg->AddEntry(ATIC,"ATIC (2006)","ep");
leg->AddEntry(CREAM,"CREAM-III (2017) ","ep");
leg->AddEntry(Nucleon,"NUCLEON (KLEM; 2017)","ep");
leg->AddEntry(pHe_flux_DAMPE,"DAMPE (2024)","ep");
leg->AddEntry(pHe_flux_DAMPE_2026,"DAMPE (2026)","ep");
leg->AddEntry(pHe_flux_DAMPE_2026,"DAMPE (2026) w smoothing","ep");
//leg->AddEntry(pHe_flux_DAMPE_2026_1,"DAMPE (2026) w N track (analy)","ep");
//leg->AddEntry(pHe_flux_DAMPE_2026_2,"DAMPE (2026) w BGO acc (analy)","ep");

leg->Draw();
/*
	TText *label = new TText();
 label-> SetNDC();
 label -> SetTextFont(2);
 label -> SetTextColor(kRed+1);
 label -> SetTextSize(0.05);
 label -> SetTextAlign(23);
 label -> SetTextAngle(0);
 label -> DrawText(0.6, 0.82, "Preliminary");
 */
C2->Update();

}

