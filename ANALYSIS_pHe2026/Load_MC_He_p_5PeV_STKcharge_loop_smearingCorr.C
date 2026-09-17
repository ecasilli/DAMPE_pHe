#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLegend.h"

using namespace std;

Double_t GetSTKLayerSignal(Double_t qX, Double_t qY, Double_t minSignal=0.){
    const Bool_t hasX = (qX > minSignal);
    const Bool_t hasY = (qY > minSignal);

    if (hasX && hasY)
        return 0.5 * (qX + qY);
    if (hasX)
        return qX;
    if (hasY)
        return qY;

    return -999.;
}

Double_t CorrectSTKProton(Double_t q)
{
    return (q - 57.3422) * 1.0719792 + 58.8187;
}


Double_t CorrectSTKHelium(Double_t q)
{
    return (q - 245.048) * 1.0976715 + 250.715;
}

void ProcessMCSpecies(TChain **chains, Int_t nsets, TH1F *hTot, TH1F **hEnergy, const Double_t *norm, const TString &speciesName, const Double_t *eMin, const Double_t *eMax) {

    const Int_t nEnergyBins = 4;

    const Bool_t isProton = (speciesName == "PROTON");
    const Bool_t isHelium = (speciesName == "HELIUM");

     if (!isProton && !isHelium) {
        cout << "ERROR: unknown species " << speciesName << endl;
        return;
    }

    // Loop sui sample MC
    for (Int_t iset = 0; iset < nsets; ++iset) {
        
        TChain *chain = chains[iset];

        Double_t MC_EnergyT, BGO_E_corr, BGO_xtr;
        Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
        Double_t STK_Y, STK_X, STK_chargeY_etaCorr[6], STK_chargeX_etaCorr[6];
        Int_t BGO_HET;

        // =============================
        // ----- Branch selection ------
        
        chain->SetBranchStatus("*", 0);
        chain->SetBranchStatus("MC_EnergyT",                   1);
        chain->SetBranchStatus("BGO_HET",                      1);
        chain->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions2", 1);
        chain->SetBranchStatus("PSD_ChargeY0",                 1);
        chain->SetBranchStatus("PSD_ChargeY1",                 1);
        chain->SetBranchStatus("PSD_ChargeX0",                 1);
        chain->SetBranchStatus("PSD_ChargeX1",                 1);
        chain->SetBranchStatus("STKtrack_to_PSD_topY",         1);
        chain->SetBranchStatus("STKtrack_to_PSD_topX",         1);
        chain->SetBranchStatus("BGO_xtr",                      1);
        chain->SetBranchStatus("STK_chargeY_etaCorr",          1);
        chain->SetBranchStatus("STK_chargeX_etaCorr",          1);

        chain->SetCacheSize(200 * 1024 * 1024);
        chain->AddBranchToCache("*", kTRUE);

        // =============================
        // ----- Branch addresses ------

        chain->SetBranchAddress("MC_EnergyT",                   &MC_EnergyT);
        chain->SetBranchAddress("BGO_HET",                      &BGO_HET);
        chain->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions2", &BGO_E_corr);
        chain->SetBranchAddress("PSD_ChargeY0",                 &PSD_CY0);
        chain->SetBranchAddress("PSD_ChargeY1",                 &PSD_CY1);
        chain->SetBranchAddress("PSD_ChargeX0",                 &PSD_CX0);
        chain->SetBranchAddress("PSD_ChargeX1",                 &PSD_CX1);
        chain->SetBranchAddress("STKtrack_to_PSD_topY",         &STK_Y);
        chain->SetBranchAddress("STKtrack_to_PSD_topX",         &STK_X);
        chain->SetBranchAddress("BGO_xtr",                      &BGO_xtr);
        chain->SetBranchAddress("STK_chargeY_etaCorr",          STK_chargeY_etaCorr);
        chain->SetBranchAddress("STK_chargeX_etaCorr",          STK_chargeX_etaCorr);

        // =============================
        // ----- Loop unico su tutti gli eventi ------

        Long64_t nEntries = chain->GetEntries();

        cout << endl;
        cout << "==============================================" << endl;
        cout << speciesName << " sample " << iset << " : " << nEntries << " events" << endl;
        cout << "==============================================" << endl;

        for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {

            if (ievt % 1000000 == 0) cout << "  " << speciesName << "[" << iset << "] " << ievt << " / " << nEntries << endl;

            chain->GetEntry(ievt);

            if (MC_EnergyT <= 20.) continue;

            // cut00 = cc204s * Trig_HEP
            if (BGO_HET <= 0)   continue;
            if (BGO_E_corr <= 20.) continue;

            // cut01: almeno un segnale su entrambi i piani di PSD
            if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
                (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;


            // cut06: accettanza STK-PSD
            if (TMath::Abs(STK_Y) >= 400.) continue;
            if (TMath::Abs(STK_X) >= 400.) continue;

            // SpCut: taglio per elettroni
            if (BGO_xtr <= 12.) continue;

            Double_t stkLayer[6];
            for (Int_t ilayer = 0; ilayer < 1; ++ilayer) {
                stkLayer[ilayer] = GetSTKLayerSignal(STK_chargeX_etaCorr[ilayer], STK_chargeY_etaCorr[ilayer]);
            }

            Double_t stkCharge = stkLayer[0];

            // ================================================================
            // SMEARING CORRECTION
            // ================================================================

            Double_t stkChargeCorr = -999.;

            if (isProton) 
                stkChargeCorr = CorrectSTKProton(stkCharge);
            
            else if (isHelium) 
                stkChargeCorr = CorrectSTKHelium(stkCharge);
            

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = norm[iset] * TMath::Power(MC_EnergyT, -1.7);

            hTot->Fill(stkChargeCorr, weight);

            // ---- Fill nell'istogramma del bin energetico corretto ----
            for (Int_t j = 0; j < nEnergyBins; ++j) {
                if (BGO_E_corr > eMin[j] && BGO_E_corr < eMax[j]) {
                    hEnergy[j]->Fill(stkChargeCorr, weight);
                    break;
                }
            }            


        } // event loop 

        cout << "Finished " << speciesName << " sample " << iset << endl;

    } // MC loop
}

void Load_MC_He_p_5PeV_STKcharge_loop_smearingCorr(){

    TString basePath;
    TString hostname = gSystem->HostName();
     
    if (hostname.Contains("cnaf")) {
        basePath = "/storage/gpfs_data/dampe/users/SKIM_2026";
    } else if (hostname.Contains("le.infn.it")) {
        basePath = "/nfs/argo/dampe/SKIM_2026";
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
    filesHe.push_back(basePath + "/HELIUM/He4_10GeV_100GeV_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_100GeV_1TeV_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_1TeV_10TeV_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_10TeV_100TeV_EPOSLHC_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_100TeV_500TeV_EPOSLHC_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_500TeV_1PeV_EPOSLHC_FTFP_merged.root");
    filesHe.push_back(basePath + "/HELIUM/He4_1PeV_5PeV_EPOSLHC_FTFP_merged.root");


    // ================================
    // Build TChains
    
    TChain *sk_p[nsetP];
    TChain *sk_he[nsetHe];

    for (Int_t i = 0; i < nsetP; ++i) {
        sk_p[i] = new TChain("newtree");
        sk_p[i]->Add(filesP[i]);
    }


    for (Int_t i = 0; i < nsetHe; ++i) {
        sk_he[i] = new TChain("newtree");
        sk_he[i]->Add(filesHe[i]);
    }

    cout << endl;
    cout << "PROTON MC" << endl;

    for (Int_t i = 0; i < nsetP; ++i) cout << "  sample " << i << " : " << sk_p[i]->GetEntries() << endl;

    cout << endl;
    cout << "HELIUM MC" << endl;

    for (Int_t i = 0; i < nsetHe; ++i) cout << "  sample " << i << " : "<< sk_he[i]->GetEntries()<< endl;

    cout << endl;

    // Array di supporto per fill rapido
/*    Double_t eMin[15] = {25.1189,  39.8107,  63.0957,  100.0,   158.489,
                         251.189,  398.107,  630.957,  1000.0,  1584.89,
                         2511.89,  3981.07,  6309.57,  10000.0, 31622.8};
    Double_t eMax[15] = {39.8107,  63.0957,  100.0,   158.489,  251.189,
                         398.107,  630.957,  1000.0,  1584.89,  2511.89,
                         3981.07,  6309.57,  10000.0, 31622.8,  100000.0};
*/
    const Int_t nEnergyBins = 4;

    Double_t eMin[nEnergyBins] = {20.0,  100.0,  1000.0,  10000.0};
    Double_t eMax[nEnergyBins] = {100.0, 1000.0, 10000.0, 100000.0};

    // ==========================================
    // MC normalization factors
    // proton: wP  * Etrue^-1.7
    // helium: wHe * Etrue^-1.7 * GeoCorr
    
    Double_t normP[nsetP] = {
        TMath::Log(10.) / 1740700000.,
        TMath::Log(10.) / 519070000.,
        TMath::Log(10.) / 201032000.,
        TMath::Log(10.) / 138744900.,
        TMath::Log(10.) / 19441600.,
        TMath::Log(5.)  / 5094200.
    };


    Double_t GeoCorr = (2. * TMath::Pi() * TMath::Pi()) / (2. * TMath::Pi() * TMath::Pi() * 1.38 * 1.38);

    Double_t normHe[nsetHe] = {

        (TMath::Log(10.) / 20380000.) * GeoCorr,
        (TMath::Log(10.) / 19556000.) * GeoCorr,
        (TMath::Log(10.) / 11346000.) * GeoCorr,
        (TMath::Log(10.) / 17694000.) * GeoCorr,
        (TMath::Log(5.) / 10305100.) * GeoCorr,
        (TMath::Log(2.) / 10281395.) * GeoCorr,
        (TMath::Log(5.) / 5064900.) * GeoCorr
    };

    // ===============================================
    // Histograms

    TH1F *hTot_p = new TH1F("hTot_p", " ", 480, 0., 1000.);
    hTot_p->GetXaxis()->SetTitle("STK charge");
    hTot_p->GetYaxis()->SetTitle("Weighted counts");
    hTot_p->SetLineColor(9);
    hTot_p->SetMarkerColor(9);
    hTot_p->Sumw2();

    TH1F *hTot_He = new TH1F("hTot_He", " ", 480, 0., 1000.);
    hTot_He->GetXaxis()->SetTitle("STK charge");
    hTot_He->GetYaxis()->SetTitle("Weighted counts");
    hTot_He->SetLineColor(6);
    hTot_He->SetMarkerColor(6);
    hTot_He->Sumw2();
    
    TH1F *hP[nEnergyBins];
    TH1F *hHe[nEnergyBins];

    for (Int_t i = 0; i < nEnergyBins; ++i) {

        TString nameP  = Form("h%02d_p" , i+1);
        TString nameHe = Form("h%02d_he", i+1);

        TString title = Form("%.0f < E_{BGO} < %.0f GeV", eMin[i], eMax[i]);

        hP[i]  = new TH1F(nameP,  title, 480, 0., 1000.);
        hHe[i] = new TH1F(nameHe, title, 480, 0., 1000.);

        hP[i]->GetXaxis()->SetTitle("STK charge");
        hP[i]->GetYaxis()->SetTitle("Weighted counts");

        hHe[i]->GetXaxis()->SetTitle("STK charge");
        hHe[i]->GetYaxis()->SetTitle("Weighted counts");

        hP[i]->SetLineColor(9);
        hP[i]->SetMarkerColor(9);

        hHe[i]->SetLineColor(6);
        hHe[i]->SetMarkerColor(6);

        hP[i]->Sumw2();
        hHe[i]->Sumw2();
    }

    // ======================================
    // Process PROTON

    cout << endl;
    cout << "Starting PROTON processing..." << endl;

     ProcessMCSpecies(sk_p, nsetP, hTot_p, hP, normP, "PROTON", eMin, eMax);

    cout << "PROTON processing finished." << endl;

    // ======================================
    // Process HELIUM

    cout << endl;
    cout << "Starting HELIUM processing..." << endl;

    ProcessMCSpecies(sk_he, nsetHe, hTot_He, hHe, normHe, "HELIUM", eMin, eMax);

    cout << "HELIUM processing finished." << endl;

    // ======================================
    // Canvas
    TCanvas *cTot = new TCanvas("cTot", "STK corrected MC - total", 900, 650);
    cTot->cd();
    gPad->SetTicks();

    hTot_p->SetLineColor(9);
    hTot_He->SetLineColor(6);

    hTot_p->Draw("HIST");
    hTot_He->Draw("HIST SAME");

    TLegend *legTot = new TLegend( 0.70, 0.75, 0.88, 0.88 );
    legTot->AddEntry( hTot_p, "Proton", "l");
    legTot->AddEntry(hTot_He, "Helium", "l");
    legTot->Draw();

    TCanvas *c[nEnergyBins];

    for (Int_t i = 0; i < nEnergyBins; ++i) {

        TString cname  = Form("c_%02d", i+1);
        TString ctitle = Form("STK charge corrected %.0f-%.0f GeV", eMin[i], eMax[i]);
        c[i] = new TCanvas(cname, ctitle, 800, 600);
        c[i]->cd();
        gPad->SetTicks();

        hP[i]->Draw("HIST");
        hHe[i]->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.70, 0.75, 0.88, 0.88);
        leg->AddEntry(hP[i] , "Proton", "l");
        leg->AddEntry(hHe[i], "Helium", "l");
        leg->Draw();
    }

    // =======================================
    // Output

    TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_STKcharge_adc_480bins_16sett26_nocut05_CORRECTED.root", "RECREATE");
    fout->cd();
    hTot_p->Write();
    hTot_He->Write();

    for (Int_t i = 0; i < nEnergyBins; ++i) {
    
        hP[i]->Write();
        hHe[i]->Write();
        c[i]->Write();
    
    }
    
    fout->Close();

    cout << endl;
    cout << "End script." << endl;


}
