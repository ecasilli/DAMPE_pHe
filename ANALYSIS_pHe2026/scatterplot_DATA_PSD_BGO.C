
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

void scatterplot_DATA_PSD_BGO()
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
Double_t STK_Y, STK_X;
Int_t    BGO_HET;

// =============================
// --------- Orbital data ------

TChain *skim = new TChain("newtree");


// Anni 2016-2024 (tutti "merged")
for (int y = 2016; y <= 2024; y++) addYear(skim, basePath, y);

// 2025: il primo bin ha suffix diverso
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2025_partially_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_010_025_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_025_050_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_050_100_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_100_500_2025_merged.root");
skim->Add(basePath + "/SKIM_2026_pHe/FLIGHT/skim_flight_500_000_2025_merged.root");

cout << "Orbital Data Entries: " << skim->GetEntries() << endl;

const int nbd = 48;
const int ndec = 6;
float arg1 = 1./float(nbd);
const int noe = nbd*ndec;
cout << "Number of energy bins: " << noe << endl;
Double_t e0 = 10.;
Double_t Ebin[noe+1];
Ebin[0] = e0;
for (int j = 1; j < noe+1; j++) {
    Ebin[j] = Ebin[j-1]*TMath::Power(10., arg1);
}

// =============================
// ----- Histograms ------------

TH2F *h2psdbgo = new TH2F("h2psdbgo","",noe,Ebin,1200,0.,30.); 
h2psdbgo->GetXaxis()->SetTitle("Deposited energy (GeV)"); 
h2psdbgo->GetYaxis()->SetTitle("PSD charge");
h2psdbgo->GetYaxis()->SetTitleOffset(1.3); 
h2psdbgo->GetXaxis()->SetTitleOffset(1.1); 
h2psdbgo->GetZaxis()->SetLabelSize(0.03); 
	

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

TFile *fout1 = new TFile("scatterplot_PHe_PSD_BGO_Orb120Month.root", "RECREATE");

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
    h2psdbgo->Fill(BGO_E_corr,charge);
}

cout << "Loop finished." << endl;

// =============================
// ----- Salvataggio -----------

fout1->cd();
h2psdbgo->Write();
fout1->Close();

cout << "End script." << endl;
}
