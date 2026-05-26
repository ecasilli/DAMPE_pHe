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

void Load_orb_120month_MLSATcor()
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
Int_t    BGO_HET, STK_ntrack;
Double_t BGO_E_corr_v1, BGO_EnergyG_SatCorr_ML_ions_v3;
Double_t BGO_slopeXZ_analy, BGO_interceptXZ_analy;
Double_t BGO_slopeYZ_analy, BGO_interceptYZ_analy;

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

// =============================
// =============================
// ----- Branch selection ------

skim->SetBranchStatus("*", 0);
skim->SetBranchStatus("BGO_HET",                        1);
skim->SetBranchStatus("BGO_EnergyG",                    1);
skim->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions",    1);
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
skim->SetBranchStatus("BGO_slopeXZ_analy",              1);
skim->SetBranchStatus("BGO_interceptXZ_analy",          1);
skim->SetBranchStatus("BGO_slopeYZ_analy",              1);
skim->SetBranchStatus("BGO_interceptYZ_analy",          1);
skim->SetBranchStatus("STK_ntrack",                     1);

// =============================
// ----- Branch addresses ------

skim->SetBranchAddress("BGO_HET",                         &BGO_HET);
skim->SetBranchAddress("BGO_EnergyG",                     &BGO_E);
skim->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions",     &BGO_E_corr_v1);
skim->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions_v3",  &BGO_EnergyG_SatCorr_ML_ions_v3);
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
skim->SetBranchAddress("BGO_slopeXZ_analy",               &BGO_slopeXZ_analy);
skim->SetBranchAddress("BGO_interceptXZ_analy",           &BGO_interceptXZ_analy);
skim->SetBranchAddress("BGO_slopeYZ_analy",               &BGO_slopeYZ_analy);
skim->SetBranchAddress("BGO_interceptYZ_analy",           &BGO_interceptYZ_analy);
skim->SetBranchAddress("STK_ntrack",                      &STK_ntrack);

// =============================

const int nbd = 6;
const int ndec = 6;
float arg1 = 1./float(nbd);
const int noe = nbd*ndec;
cout << "Number of energy bins: " << noe << endl;
float e0 = 10.;
float Ebin[noe+1];
Ebin[0] = e0;
for (int j = 1; j < noe+1; j++) {
    Ebin[j] = Ebin[j-1]*TMath::Power(10., arg1);
}

TFile *fout1 = new TFile("PHe_skim_Orb120Month_6binperdecade_2e5sigmaLow_6sigmaUp_new_noCut02.root", "RECREATE");

TH1F *h1SelBGO_orb    = new TH1F("h1SelBGO_orb",    "Selected(E_bgo) orbital", noe, Ebin);
h1SelBGO_orb->Sumw2();
TH1F *h1SelBGO_orb_v3 = new TH1F("h1SelBGO_orb_v3", "Selected(E_bgo) orbital", noe, Ebin);
h1SelBGO_orb_v3->Sumw2();


// ----- Loop su tutti gli eventi ------

Long64_t nEntries = skim->GetEntries();
cout << "Starting loop over " << nEntries << " entries..." << endl;

for (Long64_t i = 0; i < nEntries; i++) {

    if (i % 1000000 == 0) cout << "  Processing entry " << i << " / " << nEntries << endl;

    skim->GetEntry(i);

    // ---- Tagli ----
    /*
    // BGO analytical acceptance cut
    bool bgo_valid =
        !((BGO_slopeXZ_analy == 0. && BGO_interceptXZ_analy == 0.) ||
          (BGO_slopeYZ_analy == 0. && BGO_interceptYZ_analy == 0.));

    if (!bgo_valid) continue;

    bool bgo_acceptance =
        fabs(BGO_interceptXZ_analy + 448.*BGO_slopeXZ_analy) < 280. &&
        fabs(BGO_interceptYZ_analy + 448.*BGO_slopeYZ_analy) < 280. &&
        fabs(BGO_interceptXZ_analy + 46.*BGO_slopeXZ_analy)  < 280. &&
        fabs(BGO_interceptYZ_analy + 46.*BGO_slopeYZ_analy)  < 280.;

    if (!bgo_acceptance) continue;
    */
    // Ntrack 
    //if (STK_ntrack == 0) continue;

    // cut00 = cc204s * Trig_HEP
    if (BGO_HET <= 0)   continue;
    if (BGO_E   <= 20.) continue;

    // cut01: almeno un segnale su entrambi i piani di PSD
    if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
        (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;

    // cut02: profilo di shower
    //if (BGO_EneLay[0] + BGO_EneLay[1] >= BGO_EneLay[2] + BGO_EneLay[3]) continue;

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

    // 25/05/2026
    Double_t HeMPVf= (1.88586+(0.185264*log10(BGO_EnergyG_SatCorr_ML_ions_v3))+(-0.0867115*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**2)+(0.0200701*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**3)+(-0.000810826*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**4) );
    Double_t HeWidthf = (-0.0824082+(0.201315*log10(BGO_EnergyG_SatCorr_ML_ions_v3))+(-0.0989824*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**2)+(0.0218456*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**3)+(-0.00136192*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**4));
    Double_t HeGSigmaf = (0.0396516);
    Double_t HeFSig= TMath::Sqrt(HeWidthf*HeWidthf + HeGSigmaf*HeGSigmaf);

    Double_t PMPVf= (0.993603+(-0.00402039*log10(BGO_EnergyG_SatCorr_ML_ions_v3))+(0.0408375*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**2)+(-0.0227868*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**3)+(0.00494042*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**4) );
    Double_t PWidthf = (-0.0611545+(0.140328*log10(BGO_EnergyG_SatCorr_ML_ions_v3))+(-0.0683705*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**2)+(0.0143596*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**3)+(-0.000527989*log10(BGO_EnergyG_SatCorr_ML_ions_v3)**4));
    Double_t PGSigmaf = (5.5475e-08);
    Double_t PFSig= TMath::Sqrt(PWidthf*PWidthf+PGSigmaf*PGSigmaf);

    // PLow: charge < proton MPV + 2*sigma
    if (!((PMPVf - charge) < (2.5 * PFSig)))  continue;
    // HeHigh: charge - He MPV < 6*sigma
    if (!((charge - HeMPVf) < (6. * HeFSig))) continue;

    // Fill spettro energetico
    h1SelBGO_orb->Fill(BGO_E_corr_v1);
    h1SelBGO_orb_v3->Fill(BGO_EnergyG_SatCorr_ML_ions_v3);
}

cout << "Loop finished." << endl;


// =============================
// ----- Salvataggio -----------

fout1->cd();
h1SelBGO_orb->Write();
h1SelBGO_orb_v3->Write();
fout1->Close();
cout << "End script." << endl;
}
