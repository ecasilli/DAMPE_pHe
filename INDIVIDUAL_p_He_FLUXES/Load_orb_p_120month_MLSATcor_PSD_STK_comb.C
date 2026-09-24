#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

// =============================
// ----- Helper function -------

void addYear(TChain *skim, TString basePath, int year) {
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_002_010_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_010_025_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_025_050_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_050_100_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_100_500_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_500_000_%d_merged.root", year));
}

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

Double_t GetSTKLayerSignal(Double_t qX, Double_t qY, Double_t minSignal=0.) {
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

Double_t EvalPol4(Double_t x, const Double_t p[5]) {
    return p[0]
         + p[1] * x
         + p[2] * x*x
         + p[3] * x*x*x
         + p[4] * x*x*x*x;
}


// PSD DATA CHARGE BAND
// Qlow  = MPV_p  - 2.8 sigma_p
// Qhigh = MPV_He + 6.0 sigma_He
void GetPSDChargeLimits(Double_t BGOenergy, Double_t &qLow, Double_t &qHigh) {
    const Double_t x = TMath::Log10(BGOenergy);

    // PROTON DATA
    const Double_t pMPVpars[5] = { 0.988625, 0.0166979, 0.0136979, -0.0114983, 0.00313886 };
    const Double_t pWidthPars[5] = { -0.0152013, 0.0753179, -0.0344289, 0.00656681, 0.000105438 };
    const Double_t pGSigma = 2.02488e-08;

    // HELIUM DATA
    const Double_t heMPVpars[5] = { 1.92403, 0.1541, -0.082833, 0.0208528, -0.0011579 };
    const Double_t heWidthPars[5] = { -0.0105333, 0.129859, -0.0713011, 0.0169532, -0.00117928};
    const Double_t heGSigma = 0.0227879;

    // Evaluate fits
    const Double_t pMPV = EvalPol4(x, pMPVpars);
    const Double_t pWidth = EvalPol4(x, pWidthPars);
    const Double_t pSigma = TMath::Sqrt(pWidth * pWidth + pGSigma * pGSigma);

    const Double_t heMPV = EvalPol4(x, heMPVpars);
    const Double_t heWidth = EvalPol4(x, heWidthPars);
    const Double_t heSigma = TMath::Sqrt( heWidth * heWidth + heGSigma * heGSigma );

    // Final PSD band
    qLow = pMPV - 2.8 * pSigma;
    qHigh = heMPV + 6.0 * heSigma;
}

void Load_orb_120month_MLSATcor_PSD_STK_comb() 
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
Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1, PSD_Global_Charge;
Double_t STK_Y, STK_X, STK_chargeY_etaCorr[6], STK_chargeX_etaCorr[6], STK_vertexPrediction;
Int_t    BGO_HET, STK_ntrack;
Double_t BGO_E_corr_v1, BGO_EnergyG_SatCorr_ML_ions_v3;
Double_t BGO_slopeXZ_analy, BGO_interceptXZ_analy;
Double_t BGO_slopeYZ_analy, BGO_interceptYZ_analy;

// =============================
// --------- Orbital data ------

TChain *skim = new TChain("newtree");

// Anni 2016-2024 (tutti "merged")
for (int y = 2016; y <= 2025; y++) addYear(skim, basePath, y);
cout << "Orbital Data Entries: " << skim->GetEntries() << endl;

// =============================
// =============================
// ----- Branch selection ------

skim->SetBranchStatus("*", 0);
skim->SetBranchStatus("BGO_HET",                      1);
skim->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions2", 1);
skim->SetBranchStatus("PSD_ChargeY0",                 1);
skim->SetBranchStatus("PSD_ChargeY1",                 1);
skim->SetBranchStatus("PSD_ChargeX0",                 1);
skim->SetBranchStatus("PSD_ChargeX1",                 1);
skim->SetBranchStatus("STKtrack_to_PSD_topY",         1);
skim->SetBranchStatus("STKtrack_to_PSD_topX",         1);
skim->SetBranchStatus("BGO_xtr",                      1);
skim->SetBranchStatus("STK_chargeX_etaCorr",          1);
skim->SetBranchStatus("STK_chargeY_etaCorr",          1);
skim->SetBranchStatus("STK_vertexPrediction",         1);

// ROOT cache
skim->SetCacheSize(200 * 1024 * 1024);
skim->AddBranchToCache("*", kTRUE);

// =============================
// ----- Branch addresses ------

skim->SetBranchAddress("BGO_HET",                      &BGO_HET);
skim->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions2", &BGO_E_corr);
skim->SetBranchAddress("BGO_xtr",                      &BGO_xtr);
skim->SetBranchAddress("PSD_ChargeY0",                 &PSD_CY0);
skim->SetBranchAddress("PSD_ChargeY1",                 &PSD_CY1);
skim->SetBranchAddress("PSD_ChargeX0",                 &PSD_CX0);
skim->SetBranchAddress("PSD_ChargeX1",                 &PSD_CX1);
skim->SetBranchAddress("STKtrack_to_PSD_topY",         &STK_Y);
skim->SetBranchAddress("STKtrack_to_PSD_topX",         &STK_X);
skim->SetBranchAddress("STK_chargeX_etaCorr",           STK_chargeX_etaCorr);
skim->SetBranchAddress("STK_chargeY_etaCorr",           STK_chargeY_etaCorr);
skim->SetBranchAddress("STK_vertexPrediction",         &STK_vertexPrediction);


// =============================
// ENERGY BINNING

const int nbd = 5;
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

TFile *fout1 = new TFile("ROOT_FILES/PHe_skim_Orb120Month_5binperdecade_2e8sigmaLow_6sigmaUp_PSDprogr_STKcharge_comb_STKvert0e7_nocut06_17sett26.root", "RECREATE");

TH1D *h1SelBGO_orb    = new TH1D("h1SelBGO_orb",    "Selected(E_bgo) orbital", noe, Ebin);
h1SelBGO_orb->Sumw2();

// PSD progressive
const Double_t minPSDSignal = 0.2;

const Double_t dQmin = -0.3;
const Double_t dQmax =  0.7;

// Vertex decision
const Double_t vertexCut = 0.7;

// STK charge selection
const Double_t stkMin = 25.;
const Double_t stkMax = 420.;

// ----- Loop su tutti gli eventi ------

Long64_t nEntries = skim->GetEntries();
cout << "Starting loop over " << nEntries << " entries..." << endl;

for (Long64_t i = 0; i < nEntries; i++) {

    if (i % 1000000 == 0) cout << "  Processing entry " << i << " / " << nEntries << endl;

    skim->GetEntry(i);

    // ---- Tagli ----

    // cut00 = cc204s * Trig_HEP
    if (BGO_HET <= 0)   continue;
    if (BGO_E_corr <= 20.) continue;

    // cut01: almeno un segnale su entrambi i piani di PSD
    if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
        (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;

    // cut02: profilo di shower
    //if (BGO_EneLay[0] + BGO_EneLay[1] >= BGO_EneLay[2] + BGO_EneLay[3]) continue;

    // cut05: consistenza posizione BGO
    //if (TMath::Abs(BGO_cbgomax[0] - BGO_cbgostk[0]) >= 30.) continue;
    //if (TMath::Abs(BGO_cbgomax[1] - BGO_cbgostk[1]) >= 30.) continue;

    // cut06: accettanza STK-PSD
    //if (TMath::Abs(STK_Y) >= 400.) continue;
    //if (TMath::Abs(STK_X) >= 400.) continue;

    // SpCut: taglio per elettroni
    if (BGO_xtr <= 12.) continue;

    // NEW PROGRESSIVE CHARGE: asimmetrica, come nel caso adottato da Paul
    std::vector<Double_t> psdvec;

    Double_t PSDCharges[4] = { PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1};
    for (int ilay=0; ilay < 4; ilay++) {
        if (PSDCharges[ilay] > minPSDSignal) // NOTA: non ho messo proprio 0 preciso
            psdvec.push_back(PSDCharges[ilay]);
    }

    Int_t nPSDlayers = 0;

    Double_t psdCharge = GetPSDProgressiveCharge(psdvec, dQmin, dQmax, nPSDlayers);

    // STK charge (first layer)
    Double_t stkCharge = GetSTKLayerSignal(STK_chargeX_etaCorr[0], STK_chargeY_etaCorr[0]);

    // VERTEX-DEPENDENT CHARGE SELECTION
    if (STK_vertexPrediction < vertexCut) { // PSD selection
        if (psdCharge < 0.) continue;
        
        Double_t psdLow;
        Double_t psdHigh;
        GetPSDChargeLimits(BGO_E_corr, psdLow, psdHigh);
        
        if (psdCharge < psdLow) continue;
        if (psdCharge > psdHigh) continue;

    } else {                                // STK selection
            
        if (stkCharge < stkMin) continue;
        if (stkCharge > stkMax) continue;

    }
    
    h1SelBGO_orb->Fill(BGO_E_corr);

}

cout << endl;

// =============================
// ----- Salvataggio -----------

fout1->cd();
h1SelBGO_orb->Write();
fout1->ls();
fout1->Close();
cout << "End script." << endl;
}
