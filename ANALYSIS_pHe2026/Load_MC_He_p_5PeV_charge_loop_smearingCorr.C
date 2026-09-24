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

// CORREZIONE DI SMEARING
struct SmearingPars {
    Double_t mpvData[5];
    Double_t mpvMC[5];
    Double_t widthData[5];
    Double_t widthMC[5];
    Double_t sigmaData;
    Double_t sigmaMC;
};

Double_t EvalPol4LogE(Double_t E, const Double_t p[5]) {
    const Double_t x = TMath::Log10(E);
    return p[0] + p[1]*x + p[2]*x*x + p[3]*x*x*x + p[4]*x*x*x*x;
}

Double_t ApplySmearingCorrection(Double_t chargeMC, Double_t BGOenergy, const SmearingPars &p) {
    if (BGOenergy <= 0.) return chargeMC;

    const Double_t mpvData = EvalPol4LogE(BGOenergy, p.mpvData);
    const Double_t mpvMC   = EvalPol4LogE(BGOenergy, p.mpvMC);

    const Double_t widthData = EvalPol4LogE(BGOenergy, p.widthData);
    const Double_t widthMC   = EvalPol4LogE(BGOenergy, p.widthMC);

    const Double_t sigmaData = TMath::Sqrt(widthData*widthData + p.sigmaData*p.sigmaData);
    const Double_t sigmaMC   = TMath::Sqrt( widthMC*widthMC + p.sigmaMC*p.sigmaMC );

    if (sigmaMC <= 0.) return chargeMC;

    return (chargeMC - mpvMC) * (sigmaData / sigmaMC) + mpvData;
}

const SmearingPars protonSmearing = {
    // MPV DATA
    { 0.988625, 0.0166979, 0.0136979, -0.0114983, 0.00313886 },
    // MPV MC
    { 0.885521, 0.192799, -0.115362,   0.0313466,-0.00102642 },

    // WIDTH DATA
    {-0.0152013, 0.0753179, -0.0344289,  0.00656681, 0.000105438 },
    // WIDTH MC
    { 0.0362995, 0.00614734, 0.00140585, 0.00113887, 0.000375804 },

    // sigma Gauss DATA
    2.02488e-08,
    // sigma Gauss MC
    0.031548
};

const SmearingPars heliumSmearing = {
    // MPV DATA
    { 1.97763, 0.0402723, 0.00636717, -0.00967592, 0.00268719 },
    // MPV MC
    //{ 2.4225, -0.653307, 0.418625, -0.108939,   0.0118722 },
    { 2.17123,-0.236375, 0.16997,  -0.0457528,  0.00609007 },

    // WIDTH DATA
    { -0.0105333, 0.129859, -0.0713011, 0.0169532, -0.00117928 },
    // WIDTH MC
    //{  0.236337, -0.345337,  0.231388, -0.064774,   0.00709489 },
    {  0.0471634,-0.0287278, 0.0407409,-0.0157212,  0.00254262 },

    // sigma Gauss DATA
    0.0227879,
    // sigma Gauss MC
    //0.0838493
    0.0835585
};

void ProcessMCSpecies(TChain **chains, Int_t nsets, TH1F **hArr, const Double_t *norm, const TString &speciesName, const Double_t *eMin, const Double_t *eMax, const SmearingPars &smearingPars) {
    
    const Double_t dQmin = -0.3;
    const Double_t dQmax =  0.7;

    // Loop sui sample MC
    for (Int_t iset = 0; iset < nsets; ++iset) {
        
        TChain *chain = chains[iset];

        Double_t MC_EnergyT, BGO_E_corr, BGO_xtr;
        Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
        Double_t STK_Y, STK_X;
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

            // SMEARING CORRECTION
            Double_t chargeCorr =
            ApplySmearingCorrection(charge, BGO_E_corr, smearingPars);

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = norm[iset] * TMath::Power(MC_EnergyT, -1.7);

            // ---- Fill nell'istogramma del bin energetico corretto ----
            for (Int_t j = 0; j < 15; ++j) {
                if (BGO_E_corr > eMin[j] && BGO_E_corr < eMax[j]) {
                    hArr[j]->Fill(chargeCorr, weight);
                    break;
                }
            }


        } // event loop 

        cout << "Finished " << speciesName << " sample " << iset << endl;

    } // MC loop
}

void Load_MC_He_p_5PeV_charge_loop_smearingCorr(){

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

    ProcessMCSpecies(sk_p, nsetP, hP, normP, "PROTON", eMin, eMax, protonSmearing);

    cout << "PROTON processing finished." << endl;

    // ======================================
    // Process HELIUM

    cout << endl;
    cout << "Starting HELIUM processing..." << endl;

    ProcessMCSpecies(sk_he, nsetHe, hHe, normHe, "HELIUM", eMin, eMax, heliumSmearing);

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

    TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_charge_240bins_NewProgrCharge_asym_24sett26_nocut05_CORRECTED.root", "RECREATE");
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
