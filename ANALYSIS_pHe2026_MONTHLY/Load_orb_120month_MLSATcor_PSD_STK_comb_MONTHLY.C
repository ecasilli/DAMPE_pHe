#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "TChain.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TString.h"
#include "TSystem.h"

using namespace std;

// Unisce nella TChain i 6 file dei range energetici per un mese.
// Se un file manca, lo segnala; il main evita di processare quel mese.
void addYearMonth(TChain *skim, TString inputDir, int year, int month) {

    const char* energies[] = {"002_010", "010_025", "025_050", "050_100", "100_500", "500_000"};

    for (const auto& energy : energies) {
        TString filename = TString::Format(
            "%s/skim_flight_%s_%d_%02d.root",
            inputDir.Data(), energy, year, month
        );

        if (gSystem->AccessPathName(filename)) {
            cout << "WARNING: File not found: " << filename << endl;
            continue;
        }
        skim->Add(filename);
    }
}

// Funzioni per la selezione della carica: invariate rispetto alla macro originale.
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
// Qlow  = MPV_p  - 3.0 sigma_p
// Qhigh = MPV_He + 6.0 sigma_He
void GetPSDChargeLimits(Double_t BGOenergy, Double_t &qLow, Double_t &qHigh) {
    const Double_t x = TMath::Log10(BGOenergy);

    // PROTON DATA
    const Double_t pMPVpars[5] = { 0.988625, 0.0166979, 0.0136979, -0.0114983, 0.00313886 };
    const Double_t pWidthPars[5] = { -0.0152013, 0.0753179, -0.0344289, 0.00656681, 0.000105438 };
    const Double_t pGSigma = 2.02488e-08;

    // HELIUM DATA
    const Double_t heMPVpars[5] = { 1.97763, 0.0402723, 0.00636717, -0.00967592, 0.00268719 };
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
    qLow = pMPV - 3.0 * pSigma;
    qHigh = heMPV + 6.0 * heSigma;
}

// Parametri opzionali: per default elabora tutti i 24 mesi.
// Esempio per gennaio 2024:
// Load_orb_120month_MLSATcor_PSD_STK_comb_MONTHLY(2024, 2024, 1, 1)
void Load_orb_120month_MLSATcor_PSD_STK_comb_MONTHLY(
    int firstYear = 2016, int lastYear = 2025,
    int firstMonth = 1, int lastMonth = 12)
{
    if (firstYear > lastYear || firstMonth < 1 || lastMonth > 12 || firstMonth > lastMonth) {
        cerr << "ERROR: invalid year/month range" << endl;
        return;
    }

    // =============================
    // ----- Selezione del sito ----
    // Mantiene la logica originale per CNAF e Lecce.
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

    const TString inputDir = basePath + "/SKIM_2026/FLIGHT/MONTHLY";
    const TString outputDir = "ROOT_FILES";

    // Crea la directory dei risultati, se non esiste.
    if (gSystem->AccessPathName(outputDir.Data())) {
        if (gSystem->mkdir(outputDir.Data(), kTRUE) != 0) {
            cerr << "ERROR: cannot create output directory: " << outputDir << endl;
            return;
        }
    }

    // ENERGY BINNING: identico alla macro originale (5 bin/decade).
    const int nbd = 5;
    const int ndec = 6;
    const int noe = nbd * ndec;
    const float arg1 = 1. / float(nbd);
    Float_t Ebin[noe + 1];
    Ebin[0] = 10.;
    for (int j = 1; j <= noe; ++j) {
        Ebin[j] = Ebin[j-1] * TMath::Power(10., arg1);
    }

    // Parametri dei tagli: invariati.
    const Double_t minPSDSignal = 0.2;
    const Double_t dQmin = -0.3;
    const Double_t dQmax =  0.7;
    const Double_t vertexCut = 0.7;
    const Double_t stkMin = 25.;
    const Double_t stkMax = 450.;

    int success = 0;
    int skipped = 0;
    int failed = 0;

    for (int year = firstYear; year <= lastYear; ++year) {
        for (int month = firstMonth; month <= lastMonth; ++month) {

            cout << "\n======================================" << endl;
            cout << "Analisi: " << year << "-" << TString::Format("%02d", month) << endl;

            // Una TChain nuova per ogni mese; non mescola eventi di mesi diversi.
            TChain skim("newtree");
            addYearMonth(&skim, inputDir, year, month);

            // Il mese deve contenere tutti i 6 file energetici.
            if (skim.GetListOfFiles()->GetEntries() != 6) {
                cerr << "SKIPPED: " << year << "-" << TString::Format("%02d", month) << " (presenti " << skim.GetListOfFiles()->GetEntries() << " file energetici su 6)" << endl;
                ++skipped;
                continue;
            }

            const Long64_t nEntries = skim.GetEntries();
            // Controllo anche che in tutti i 6 ROOT ci sia l'albero newtree.
            if (skim.GetNtrees() != 6 || nEntries < 0) {
                cerr << "ERROR: missing/unreadable newtree for " << year << "-" << TString::Format("%02d", month) << " (trees=" << skim.GetNtrees() << ")" << endl;
                ++failed;
                continue;
            }
            cout << "Input files: 6 | Entries: " << nEntries << endl;

            // Variabili originali utilizzate dai tagli e dai branch.
            Double_t BGO_E_corr = 0., BGO_xtr = 0.;
            Double_t PSD_CY0 = 0., PSD_CY1 = 0., PSD_CX0 = 0., PSD_CX1 = 0.;
            Double_t STK_Y = 0., STK_X = 0.;
            Double_t STK_chargeY_etaCorr[6] = {0.};
            Double_t STK_chargeX_etaCorr[6] = {0.};
            Double_t STK_vertexPrediction = 0.;
            Int_t BGO_HET = 0;

            skim.SetBranchStatus("*", 0);
            skim.SetBranchStatus("BGO_HET",                      1);
            skim.SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions2", 1);
            skim.SetBranchStatus("PSD_ChargeY0",                 1);
            skim.SetBranchStatus("PSD_ChargeY1",                 1);
            skim.SetBranchStatus("PSD_ChargeX0",                 1);
            skim.SetBranchStatus("PSD_ChargeX1",                 1);
            skim.SetBranchStatus("STKtrack_to_PSD_topY",         1);
            skim.SetBranchStatus("STKtrack_to_PSD_topX",         1);
            skim.SetBranchStatus("BGO_xtr",                      1);
            skim.SetBranchStatus("STK_chargeX_etaCorr",          1);
            skim.SetBranchStatus("STK_chargeY_etaCorr",          1);
            skim.SetBranchStatus("STK_vertexPrediction",         1);

            skim.SetCacheSize(200 * 1024 * 1024);
            skim.AddBranchToCache("*", kTRUE);

            skim.SetBranchAddress("BGO_HET",                      &BGO_HET);
            skim.SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions2", &BGO_E_corr);
            skim.SetBranchAddress("BGO_xtr",                      &BGO_xtr);
            skim.SetBranchAddress("PSD_ChargeY0",                 &PSD_CY0);
            skim.SetBranchAddress("PSD_ChargeY1",                 &PSD_CY1);
            skim.SetBranchAddress("PSD_ChargeX0",                 &PSD_CX0);
            skim.SetBranchAddress("PSD_ChargeX1",                 &PSD_CX1);
            skim.SetBranchAddress("STKtrack_to_PSD_topY",         &STK_Y);
            skim.SetBranchAddress("STKtrack_to_PSD_topX",         &STK_X);
            skim.SetBranchAddress("STK_chargeX_etaCorr",           STK_chargeX_etaCorr);
            skim.SetBranchAddress("STK_chargeY_etaCorr",           STK_chargeY_etaCorr);
            skim.SetBranchAddress("STK_vertexPrediction",         &STK_vertexPrediction);

            // Istogramma nuovo per ogni mese.
            TH1D h1SelBGO_orb("h1SelBGO_orb", "Selected(E_bgo) orbital", noe, Ebin);
            h1SelBGO_orb.SetDirectory(nullptr);
            h1SelBGO_orb.Sumw2();

            Long64_t selected = 0;
            bool readError = false;
            std::vector<Double_t> psdvec;
            psdvec.reserve(4); // Evita riallocazioni per ogni evento.
            for (Long64_t i = 0; i < nEntries; ++i) {
                if (i % 1000000 == 0) {
                    cout << "  Processing entry " << i << " / " << nEntries << endl;
                }
                if (skim.GetEntry(i) <= 0) {
                    cerr << "ERROR: failed to read entry " << i << endl;
                    readError = true;
                    break;
                }

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
                psdvec.clear();

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

                h1SelBGO_orb.Fill(BGO_E_corr);

                ++selected;
            }

            // Evita di salvare un istogramma incompleto se una entry non e' leggibile.
            if (readError) {
                ++failed;
                continue;
            }

            const TString outputFile = TString::Format(
                "%s/PHe_skim_%d_%02d_5binperdecade_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_STKvert0e7_27ott26.root",
                outputDir.Data(), year, month
            );

            TFile fout(outputFile.Data(), "RECREATE");
            if (fout.IsZombie()) {
                cerr << "ERROR creating output: " << outputFile << endl;
                ++failed;
                continue;
            }
            fout.cd();
            const Int_t bytesWritten = h1SelBGO_orb.Write();
            fout.Close();

            if (bytesWritten <= 0) {
                cerr << "ERROR writing histogram: " << outputFile << endl;
                ++failed;
                continue;
            }
            cout << "Saved: " << outputFile << endl;
            cout << "Selected events: " << selected << endl;
            ++success;
        }
    }

    cout << "\n======================================" << endl;
    cout << "Monthly analysis summary" << endl;
    cout << "Successful: " << success << endl;
    cout << "Skipped:    " << skipped << endl;
    cout << "Failed:     " << failed << endl;
    cout << "======================================" << endl;
}
