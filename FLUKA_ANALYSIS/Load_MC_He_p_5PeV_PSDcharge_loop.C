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

void ProcessMCSpecies(TChain **chains, Int_t nsets, TH1F **hArr, const Double_t *norm, const TString &speciesName, const Double_t *eMin, const Double_t *eMax) {
    
    const Double_t dQmin = -0.3;
    const Double_t dQmax =  0.7;

    // Loop sui sample MC
    for (Int_t iset = 0; iset < nsets; ++iset) {
        
        TChain *chain = chains[iset];

        Double_t MC_EnergyT, BGO_E_corr, BGO_xtr;
        Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
        Double_t STK_Y, STK_X;
        Int_t BGO_HET;
        Double_t MC_weight_p_FLUKA, MC_weight_he_FLUKA;

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
        chain->SetBranchStatus("BGO_xtr", 1);

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

            // NEW PROGRESSIVE CHARGE: asimmetrica, come nel caso adottato da Paul
            std::vector<Double_t> psdvec;

            Double_t PSDCharges[4] = { PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1};
            for (int ilay=0; ilay < 4; ilay++) {
                if (PSDCharges[ilay] > 0.2) // NOTA: non ho messo proprio 0 preciso
                    psdvec.push_back(PSDCharges[ilay]);
            }

            Int_t nPSDlayers = 0;

            Double_t charge = GetPSDProgressiveCharge(psdvec, dQmin, dQmax, nPSDlayers);
            if (charge < 0.2) continue;

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = norm[iset] * TMath::Power(MC_EnergyT, -1.7);

            // ---- Fill nell'istogramma del bin energetico corretto ----
            for (Int_t j = 0; j < 15; ++j) {
                if (BGO_E_corr > eMin[j] && BGO_E_corr < eMax[j]) {
                    hArr[j]->Fill(charge, weight);
                    break;
                }
            }


        } // event loop 

        cout << "Finished " << speciesName << " sample " << iset << endl;

    } // MC loop
}

void Load_MC_He_p_5PeV_PSDcharge_loop(){

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
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_10GeV_100GeV_FLUKA_merged.root");
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_100GeV_1TeV_FLUKA_merged.root");
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_1TeV_10TeV_FLUKA_merged.root");
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_10TeV_100TeV_FLUKA_merged.root");
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_100TeV_1PeV_FLUKA_merged.root");
    filesP.push_back(basePath + "/PROTON/FLUKA/Proton_1PeV_5PeV_FLUKA_merged.root");

    std::vector<TString> filesHe;
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p18_10GeV_100GeV_FLUKA_partial.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p18_100GeV_1TeV_FLUKA_merged.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p18_1TeV_10TeV_FLUKA_merged.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p18_10TeV_100TeV_FLUKA_merged.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p10_100TeV_500TeV_FLUKA_merged.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p18_500TeV_1PeV_FLUKA_merged.root");
    filesHe.push_back(basePath + "/HELIUM/FLUKA/He4-v6r0p20_1PeV_5PeV_FLUKA_merged.root");


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
        TMath::Log(10.) / 1055454555., // 10 GeV - 100 GeV 
        TMath::Log(10.) / 579900871., // 100 GeV - 1 TeV 
        TMath::Log(10.) / 274526589., // 1 TeV - 10 TeV
        TMath::Log(10.) / 208226000., // 10 TeV - 100 TeV
        TMath::Log(10.) / 9999400., // 100 TeV - 1 PeV
        TMath::Log(5.)  / 2564400. // 1 PeV - 5 PeV
    };


    Double_t GeoCorr = (2. * TMath::Pi() * TMath::Pi()) / (2. * TMath::Pi() * TMath::Pi() * 1.38 * 1.38);

    Double_t normHe[nsetHe] = {

        //(TMath::Log(10.) / 198586517.) * GeoCorr, // 10 GeV - 100 GeV, v6r0p10
        (TMath::Log(10.) / 474272213.) * GeoCorr, // 10 GeV - 100 GeV, v6r0p18
        (TMath::Log(10.) / 66009860.) * GeoCorr, // 100 GeV - 1 TeV, v6r0p18
        (TMath::Log(10.) / 77281442.) * GeoCorr, // 1 TeV - 10 TeV, v6r0p18
        (TMath::Log(10.) / 19993301.) * GeoCorr, // 10 TeV - 100 TeV, v6r0p18
        (TMath::Log(5.) / 20132200.) * GeoCorr, // 100 TeV - 500 TeV, v6r0p10
        (TMath::Log(2.) / 1246818.) * GeoCorr, // 500 TeV - 1 PeV, v6r0p18
        (TMath::Log(5.) / 2263134.) * GeoCorr // 1 PeV - 5 PeV
    };

    // ===============================================
    // Histograms

    TH1F *hP[15];
    TH1F *hHe[15];

    for (Int_t i = 0; i < 15; ++i) {

        TString nameP  = Form("h%02d_p" , i+1);
        TString nameHe = Form("h%02d_he", i+1);

        TString title = Form("%.0f < E_{BGO} < %.0f GeV", eMin[i], eMax[i]);

        hP[i]  = new TH1F(nameP,  title, 240, 0., 5.);
        hHe[i] = new TH1F(nameHe, title, 240, 0., 5.);

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

    TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_FLUKA_charge_240bins_NewProgrCharge_asym_11nov26.root", "RECREATE");
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
