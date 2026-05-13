#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>

void Load_orb_120month_MLSATcor_charge()
{
// =============================
// ----- Dichiarazioni variabili

// Variabili per SetBranchAddress
Double_t BGO_E, BGO_E_corr, BGO_xtr;
Double_t BGO_EneLay[14], BGO_cbgomax[14], BGO_cbgostk[14];
Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
Double_t STK_Y, STK_X;
Int_t    BGO_HET;

// Array di supporto per fill rapido
Double_t eMin[15] = {25.1189,  39.8107,  63.0957,  100.0,   158.489,
                     251.189,  398.107,  630.957,  1000.0,  1584.89,
                     2511.89,  3981.07,  6309.57,  10000.0, 31622.8};
Double_t eMax[15] = {39.8107,  63.0957,  100.0,   158.489,  251.189,
                     398.107,  630.957,  1000.0,  1584.89,  2511.89,
                     3981.07,  6309.57,  10000.0, 31622.8,  100000.0};

// =============================
// --------- Orbital data ------

TChain *skim = new TChain("newtree");

// 2016 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2016_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2016_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2016_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2016_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2016_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2016_merged.root");

// 2017 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2017_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2017_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2017_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2017_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2017_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2017_merged.root");

// 2018 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2018_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2018_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2018_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2018_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2018_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2018_merged.root");

// 2019 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2019_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2019_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2019_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2019_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2019_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2019_merged.root");

// 2020 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2020_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2020_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2020_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2020_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2020_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2020_merged.root");

// 2021 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2021_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2021_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2021_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2021_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2021_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2021_merged.root");

// 2022 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2022_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2022_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2022_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2022_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2022_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2022_merged.root");

// 2023 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2023_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2023_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2023_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2023_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2023_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2023_merged.root");

// 2024 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2024_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2024_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2024_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2024_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2024_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2024_merged.root");

// 2025 (Jan-Dec)
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2025_partially_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2025_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2025_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2025_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2025_merged.root");
skim->Add("/storage/gpfs_data/dampe/users/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2025_merged.root");

cout << "Orbital Data Entries: " << skim->GetEntries() << endl;

// =============================
// ----- Histograms ------------

    TH1F *h01 = new TH1F("h01", "25 < E_{BGO} < 39", 240, 0., 5.);
    h01->GetXaxis()->SetTitle("PSD charge");
    h01->GetYaxis()->SetTitle("events");
    h01->SetLineColor(kBlack);
    h01->SetMarkerColor(kBlack);
    h01->Sumw2();
 
    TH1F *h02 = new TH1F("h02", "39 < E_{BGO} < 63", 240, 0., 5.);
    h02->GetXaxis()->SetTitle("PSD charge");
    h02->GetYaxis()->SetTitle("events");
    h02->SetLineColor(kBlack);
    h02->SetMarkerColor(kBlack);
    h02->Sumw2();
 
    TH1F *h03 = new TH1F("h03", "63 < E_{BGO} < 100", 240, 0., 5.);
    h03->GetXaxis()->SetTitle("PSD charge");
    h03->GetYaxis()->SetTitle("events");
    h03->SetLineColor(kBlack);
    h03->SetMarkerColor(kBlack);
    h03->Sumw2();
 
    TH1F *h04 = new TH1F("h04", "100 < E_{BGO} < 158", 240, 0., 5.);
    h04->GetXaxis()->SetTitle("PSD charge");
    h04->GetYaxis()->SetTitle("events");
    h04->SetLineColor(kBlack);
    h04->SetMarkerColor(kBlack);
    h04->Sumw2();
 
    TH1F *h05 = new TH1F("h05", "158 < E_{BGO} < 251", 240, 0., 5.);
    h05->GetXaxis()->SetTitle("PSD charge");
    h05->GetYaxis()->SetTitle("events");
    h05->SetLineColor(kBlack);
    h05->SetMarkerColor(kBlack);
    h05->Sumw2();
 
    TH1F *h06 = new TH1F("h06", "251 < E_{BGO} < 398", 240, 0., 5.);
    h06->GetXaxis()->SetTitle("PSD charge");
    h06->GetYaxis()->SetTitle("events");
    h06->SetLineColor(kBlack);
    h06->SetMarkerColor(kBlack);
    h06->Sumw2();
 
    TH1F *h07 = new TH1F("h07", "398 < E_{BGO} < 630", 240, 0., 5.);
    h07->GetXaxis()->SetTitle("PSD charge");
    h07->GetYaxis()->SetTitle("events");
    h07->SetLineColor(kBlack);
    h07->SetMarkerColor(kBlack);
    h07->Sumw2();
 
    TH1F *h08 = new TH1F("h08", "630 < E_{BGO} < 1000", 240, 0., 5.);
    h08->GetXaxis()->SetTitle("PSD charge");
    h08->GetYaxis()->SetTitle("events");
    h08->SetLineColor(kBlack);
    h08->SetMarkerColor(kBlack);
    h08->Sumw2();
 
    TH1F *h09 = new TH1F("h09", "1000 < E_{BGO} < 1584", 240, 0., 5.);
    h09->GetXaxis()->SetTitle("PSD charge");
    h09->GetYaxis()->SetTitle("events");
    h09->SetLineColor(kBlack);
    h09->SetMarkerColor(kBlack);
    h09->Sumw2();
 
    TH1F *h10 = new TH1F("h10", "1584 < E_{BGO} < 2511", 240, 0., 5.);
    h10->GetXaxis()->SetTitle("PSD charge");
    h10->GetYaxis()->SetTitle("events");
    h10->SetLineColor(kBlack);
    h10->SetMarkerColor(kBlack);
    h10->Sumw2();
 
    TH1F *h11 = new TH1F("h11", "2511 < E_{BGO} < 3981", 240, 0., 5.);
    h11->GetXaxis()->SetTitle("PSD charge");
    h11->GetYaxis()->SetTitle("events");
    h11->SetLineColor(kBlack);
    h11->SetMarkerColor(kBlack);
    h11->Sumw2();
 
    TH1F *h12 = new TH1F("h12", "3981 < E_{BGO} < 6309", 240, 0., 5.);
    h12->GetXaxis()->SetTitle("PSD charge");
    h12->GetYaxis()->SetTitle("events");
    h12->SetLineColor(kBlack);
    h12->SetMarkerColor(kBlack);
    h12->Sumw2();
 
    TH1F *h13 = new TH1F("h13", "6309 < E_{BGO} < 10000", 240, 0., 5.);
    h13->GetXaxis()->SetTitle("PSD charge");
    h13->GetYaxis()->SetTitle("events");
    h13->SetLineColor(kBlack);
    h13->SetMarkerColor(kBlack);
    h13->Sumw2();
 
    TH1F *h14 = new TH1F("h14", "10000 < E_{BGO} < 31622", 240, 0., 5.);
    h14->GetXaxis()->SetTitle("PSD charge");
    h14->GetYaxis()->SetTitle("events");
    h14->SetLineColor(kBlack);
    h14->SetMarkerColor(kBlack);
    h14->Sumw2();
 
    TH1F *h15 = new TH1F("h15", "31622 < E_{BGO} < 100000", 240, 0., 5.);
    h15->GetXaxis()->SetTitle("PSD charge");
    h15->GetYaxis()->SetTitle("events");
    h15->SetLineColor(kBlack);
    h15->SetMarkerColor(kBlack);
    h15->Sumw2();


TH1F* hArr[15] = {h01, h02, h03, h04, h05, h06, h07, h08,
                  h09, h10, h11, h12, h13, h14, h15};

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

    // ---- Calcolo PSD charge (simple average sui segnali positivi) ----
    Double_t num = 0., den = 0.;
    if (PSD_CY0 > 0.) { num += PSD_CY0; den += 1.; }
    if (PSD_CY1 > 0.) { num += PSD_CY1; den += 1.; }
    if (PSD_CX0 > 0.) { num += PSD_CX0; den += 1.; }
    if (PSD_CX1 > 0.) { num += PSD_CX1; den += 1.; }
    if (den == 0.) continue;
    Double_t charge = num / den;

    // ---- Fill nell'istogramma del bin energetico corretto ----
    for (int j = 0; j < 15; j++) {
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
c0->Divide(3, 5);
for (int j = 0; j < 15; j++) {
    c0->cd(j + 1);
    gPad->SetTicks();
    hArr[j]->Draw();
}

// =============================
// ----- Salvataggio -----------

TFile *fout1 = new TFile("PHe_charge_Orb120Month_240bins.root", "RECREATE");
fout1->cd();
fout1->Close();

cout << "End script." << endl;
}
