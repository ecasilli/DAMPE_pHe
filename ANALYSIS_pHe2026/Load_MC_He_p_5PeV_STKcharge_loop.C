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

Double_t GetPSDProgressiveCharge(const std::vector<Double_t>& q, Double_t dQmin, Double_t dQmax, Int_t &nLayersUsed) {
    nLayersUsed = 0;

    if (q.empty()) return -999.;

    // Il primo layer valido viene sempre preso
    Double_t sum = q[0];
    nLayersUsed = 1;

    for (size_t i=1; i < q.size(); ++i) {
        Double_t dQ = q[i] - q[i-1];

        if (dQ > dQmin && dQ < dQmax){
            sum += q[i];
            nLayersUsed++;
        } 
        else {
            // ci fermiamo alla prima incopatibilità
            break;
        }
    }

    return sum / ((Double_t)nLayersUsed);
}

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

void ProcessMCSpecies(TChain **chains, Int_t nsets, TH1F **hArr, const Double_t *norm, const TString &speciesName, const Double_t *eMin, const Double_t *eMax) {
    
    const Double_t dQmin = -0.3;
    const Double_t dQmax =  0.7;

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
            for (Int_t ilayer = 0; ilayer < 2; ++ilayer) {
                stkLayer[ilayer] = GetSTKLayerSignal(STK_chargeX_etaCorr[ilayer], STK_chargeY_etaCorr[ilayer]);
            }

            Double_t stkFirstLayer = stkLayer[0];
            Double_t stkCharge = stkFirstLayer;

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = norm[iset] * TMath::Power(MC_EnergyT, -1.7);

            // ---- Fill nell'istogramma del bin energetico corretto ----
            for (Int_t j = 0; j < 15; ++j) {
                if (BGO_E_corr > eMin[j] && BGO_E_corr < eMax[j]) {
                    hArr[j]->Fill(stkCharge, weight);
                    break;
                }
            }


        } // event loop 

        cout << "Finished " << speciesName << " sample " << iset << endl;

    } // MC loop
}

void Load_MC_He_p_5PeV_STKcharge_loop(){

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
    Double_t eMin[15] = {25.1189,  39.8107,  63.0957,  100.0,   158.489,
                         251.189,  398.107,  630.957,  1000.0,  1584.89,
                         2511.89,  3981.07,  6309.57,  10000.0, 31622.8};
    Double_t eMax[15] = {39.8107,  63.0957,  100.0,   158.489,  251.189,
                         398.107,  630.957,  1000.0,  1584.89,  2511.89,
                         3981.07,  6309.57,  10000.0, 31622.8,  100000.0};

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

    TH1F *hP[15];
    TH1F *hHe[15];

    for (Int_t i = 0; i < 15; ++i) {

        TString nameP  = Form("h%02d_p" , i+1);
        TString nameHe = Form("h%02d_he", i+1);

        TString title = Form("%.0f < E_{BGO} < %.0f GeV", eMin[i], eMax[i]);

        hP[i]  = new TH1F(nameP,  title, 480, 0., 1000.);
        hHe[i] = new TH1F(nameHe, title, 480, 0., 1000.);

        hP[i]->GetXaxis()->SetTitle("PSD charge");
        hP[i]->GetYaxis()->SetTitle("Weighted counts");

        hHe[i]->GetXaxis()->SetTitle("PSD charge");
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

    ProcessMCSpecies(sk_p, nsetP, hP, normP, "PROTON", eMin, eMax);

    cout << "PROTON processing finished." << endl;

    // ======================================
    // Process HELIUM

    cout << endl;
    cout << "Starting HELIUM processing..." << endl;

    ProcessMCSpecies(sk_he, nsetHe, hHe, normHe, "HELIUM", eMin, eMax);

    cout << "HELIUM processing finished." << endl;

    // ======================================
    // Canvas

    TCanvas *c[15];

    for (Int_t i = 0; i < 15; ++i) {

        TString cname  = Form("c_%02d", i+1);
        TString ctitle = Form("PSD charge %.0f-%.0f GeV", eMin[i], eMax[i]);
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

    TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_STKcharge_adc_480bins_14sett26_nocut05.root", "RECREATE");
    fout->cd();

    for (Int_t i = 0; i < 15; ++i) {
    
        hP[i]->Write();
        hHe[i]->Write();
        c[i]->Write();
    
    }
    
    fout->Close();

    cout << endl;
    cout << "End script." << endl;


}
