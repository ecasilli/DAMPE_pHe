#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>

// =============================
// ----- Helper function -------

void addYear(TChain *skim, TString basePath, int year) {
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_%d_merged.root", year));
}

void Load_orb_120month_MLSATcor_STK_charge()
{

// =============================
// ----- Selezione del sito ----
TString basePath;
TString hostname = gSystem->HostName();
if (hostname.Contains("cnaf")) {
    basePath = "/storage/gpfs_data/dampe/users";
} else if (hostname.Contains("le.infn")) {
    basePath = "/nfs/argo/dampe";
} else {
    cout << "WARNING: hostname non riconosciuto: " << hostname << endl;
    return;
}

// =============================
// ----- Dichiarazioni variabili

// Variabili per SetBranchAddress
Double_t BGO_E, BGO_E_corr, BGO_xtr;
Double_t BGO_EneLay[14], BGO_cbgomax[14], BGO_cbgostk[14];
Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
Double_t STK_Y, STK_X, STK_Y_charge[6], STK_X_charge[6];
Int_t    BGO_HET;

// Array di supporto per fill rapido
Double_t eMin[5] = {10.0, 100.0, 1000.0, 10000.0, 100000.0};
Double_t eMax[5] = {100.0, 1000.0, 10000.0, 100000.0, 1000000.0};

// =============================
// --------- Orbital data ------

TChain *skim = new TChain("newtree");

// Anni 2016-2024 (tutti "merged")
for (int y = 2016; y <= 2025; y++) addYear(skim, basePath, y);
/*
// 2025: il primo bin ha suffix diverso
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2025_partially_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2025_merged.root");
*/
cout << "Orbital Data Entries: " << skim->GetEntries() << endl;

// =============================
// ----- Histograms ------------

    TH1F *h01 = new TH1F("h01", "10 GeV < E_{BGO} < 100 GeV", 240, 0., 1000.);
    h01->GetXaxis()->SetTitle("STK signal (ADC counts)");
    h01->GetYaxis()->SetTitle("Number of events");
    h01->SetLineColor(kBlack);
    h01->SetMarkerColor(kBlack);
    h01->Sumw2();
 
    TH1F *h02 = new TH1F("h02", "100 GeV < E_{BGO} < 1 TeV", 240, 0., 1000.);
    h02->GetXaxis()->SetTitle("STK signal (ADC counts)");
    h02->GetYaxis()->SetTitle("Number of events");
    h02->SetLineColor(kBlack);
    h02->SetMarkerColor(kBlack);
    h02->Sumw2();
 
    TH1F *h03 = new TH1F("h03", "1 TeV < E_{BGO} < 10 TeV", 240, 0., 1000.);
    h03->GetXaxis()->SetTitle("STK signal (ADC counts)");
    h03->GetYaxis()->SetTitle("Number of events");
    h03->SetLineColor(kBlack);
    h03->SetMarkerColor(kBlack);
    h03->Sumw2();
 
    TH1F *h04 = new TH1F("h04", "10 TeV < E_{BGO} < 100 TeV", 240, 0., 1000.);
    h04->GetXaxis()->SetTitle("STK signal (ADC counts)");
    h04->GetYaxis()->SetTitle("Number of events");
    h04->SetLineColor(kBlack);
    h04->SetMarkerColor(kBlack);
    h04->Sumw2();
 
    TH1F *h05 = new TH1F("h05", "100 TeV < E_{BGO} < 1 PeV", 240, 0., 1000.);
    h05->GetXaxis()->SetTitle("STK signal (ADC counts)");
    h05->GetYaxis()->SetTitle("Number of events");
    h05->SetLineColor(kBlack);
    h05->SetMarkerColor(kBlack);
    h05->Sumw2();

TH1F* hArr[5] = {h01, h02, h03, h04, h05};

// =============================
// ----- Branch selection ------

skim->SetBranchStatus("*", 0);
skim->SetBranchStatus("BGO_HET",                        1);
skim->SetBranchStatus("BGO_EnergyG",                    1);
skim->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions_v3", 1);
skim->SetBranchStatus("BGO_EneLay",                     1);
skim->SetBranchStatus("BGO_cbgomax",                    1);
skim->SetBranchStatus("BGO_cbgostk",                    1);
skim->SetBranchStatus("PSD_ChargeY0",                   1);
skim->SetBranchStatus("PSD_ChargeY1",                   1);
skim->SetBranchStatus("PSD_ChargeX0",                   1);
skim->SetBranchStatus("PSD_ChargeX1",                   1);
skim->SetBranchStatus("STKtrack_to_PSD_topY",           1);
skim->SetBranchStatus("STKtrack_to_PSD_topX",           1);
skim->SetBranchStatus("BGO_xtr",                        1);
skim->SetBranchStatus("STK_chargeX",                    1);
skim->SetBranchStatus("STK_chargeY",                    1);

// =============================
// ----- Branch addresses ------

skim->SetBranchAddress("BGO_HET",                         &BGO_HET);
skim->SetBranchAddress("BGO_EnergyG",                     &BGO_E);
skim->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions_v3",  &BGO_E_corr);
skim->SetBranchAddress("BGO_xtr",                         &BGO_xtr);
skim->SetBranchAddress("BGO_EneLay",                       BGO_EneLay);  // no & per array
skim->SetBranchAddress("BGO_cbgomax",                      BGO_cbgomax); // no & per array
skim->SetBranchAddress("BGO_cbgostk",                      BGO_cbgostk); // no & per array
skim->SetBranchAddress("PSD_ChargeY0",                    &PSD_CY0);
skim->SetBranchAddress("PSD_ChargeY1",                    &PSD_CY1);
skim->SetBranchAddress("PSD_ChargeX0",                    &PSD_CX0);
skim->SetBranchAddress("PSD_ChargeX1",                    &PSD_CX1);
skim->SetBranchAddress("STKtrack_to_PSD_topY",            &STK_Y);
skim->SetBranchAddress("STKtrack_to_PSD_topX",            &STK_X);
skim->SetBranchAddress("STK_chargeX",                      STK_X_charge);
skim->SetBranchAddress("STK_chargeY",                      STK_Y_charge);

// =============================
// ----- Loop unico su tutti gli eventi ------

Long64_t nEntries = skim->GetEntries();
cout << "Starting loop over " << nEntries << " entries..." << endl;

for (Long64_t i = 0; i < nEntries; i++) {

    if (i % 1000000 == 0) cout << "  Processing entry " << i << " / " << nEntries << endl;

    skim->GetEntry(i);

    // ---- Tagli ----

    // cut00 = cc204s * Trig_HEP
    if (BGO_HET <= 0)   continue;
    if (BGO_E   <= 20.) continue;

    // cut01: almeno un segnale su entrambi i piani di PSD
    if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
        (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;

    // cut02: profilo di shower
    if (BGO_EneLay[0] + BGO_EneLay[1] >= BGO_EneLay[2] + BGO_EneLay[3]) continue;

    // cut05: consistenza posizione BGO
    if (TMath::Abs(BGO_cbgomax[0] - BGO_cbgostk[0]) >= 30.) continue;
    if (TMath::Abs(BGO_cbgomax[1] - BGO_cbgostk[1]) >= 30.) continue;

    // cut06: accettanza STK-PSD
    if (TMath::Abs(STK_Y) >= 400.) continue;
    if (TMath::Abs(STK_X) >= 400.) continue;

    // SpCut: taglio per elettroni
    if (BGO_xtr <= 12.) continue;

    // ---- Calcolo STK charge (simple average sui segnali positivi) ----
    Double_t num = 0., den = 0.;
    if (STK_Y_charge[0] > 0.) { num += STK_Y_charge[0]; den += 1.; }
    if (STK_X_charge[0] > 0.) { num += STK_X_charge[0]; den += 1.; }
    if (den == 0.) continue;
    Double_t charge = num / den;

    // ---- Fill nell'istogramma del bin energetico corretto ----
    for (int j = 0; j < 5; j++) {
        if (BGO_E_corr > eMin[j] && BGO_E_corr < eMax[j]) {
            hArr[j]->Fill(charge);
            break;
        }
    }
}

cout << "Loop finished." << endl;

// =============================
// ----- Canvas ----------------

TCanvas *c0 = new TCanvas("c0", "BGO-PSDglob", 1200, 900);
c0->Divide(2, 3);
for (int j = 0; j < 5; j++) {
    c0->cd(j + 1);
    gPad->SetTicks();
    hArr[j]->Draw();
}

// =============================
// ----- Salvataggio -----------

TFile *fout1 = new TFile("PHe_STK_Orb120Month_240bins.root", "RECREATE");
fout1->cd();
h01->Write();
h02->Write();
h03->Write();
h04->Write();
h05->Write();
fout1->Close();

cout << "End script." << endl;
}
