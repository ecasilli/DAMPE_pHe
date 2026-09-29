
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
#include "TNamed.h"

#include "ChargeSystematicsConfig.h"

using namespace std;

void addYear(TChain *skim, TString basePath, int year) {
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_002_010_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_010_025_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_025_050_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_050_100_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_100_500_%d_merged.root", year));
    skim->Add(basePath + TString::Format("/SKIM_2026/FLIGHT/skim_flight_500_000_%d_merged.root", year));
}

Double_t GetPSDProgressiveCharge(const std::vector<Double_t>& q, Double_t dQmin,
                                 Double_t dQmax, Int_t &nLayersUsed) {
    nLayersUsed = 0;
    if (q.empty()) return -999.;
    Double_t sum = q[0];
    nLayersUsed = 1;
    for (size_t i=1; i < q.size(); ++i) {
        Double_t dQ = q[i] - q[i-1];
        if (dQ > dQmin && dQ < dQmax) {
            sum += q[i];
            nLayersUsed++;
        } else break;
    }
    return sum / ((Double_t)nLayersUsed);
}

Double_t GetSTKLayerSignal(Double_t qX, Double_t qY, Double_t minSignal=0.) {
    const Bool_t hasX = (qX > minSignal);
    const Bool_t hasY = (qY > minSignal);
    if (hasX && hasY) return 0.5 * (qX + qY);
    if (hasX) return qX;
    if (hasY) return qY;
    return -999.;
}

Double_t EvalPol4(Double_t x, const Double_t p[5]) {
    return p[0] + p[1]*x + p[2]*x*x + p[3]*x*x*x + p[4]*x*x*x*x;
}

void GetPSDChargeLimits(Double_t BGOenergy, Double_t &qLow, Double_t &qHigh) {
    const Double_t x = TMath::Log10(BGOenergy);
    const Double_t pMPVpars[5] = { 0.988625, 0.0166979, 0.0136979, -0.0114983, 0.00313886 };
    const Double_t pWidthPars[5] = { -0.0152013, 0.0753179, -0.0344289, 0.00656681, 0.000105438 };
    const Double_t pGSigma = 2.02488e-08;
    const Double_t heMPVpars[5] = { 1.97763, 0.0402723, 0.00636717, -0.00967592, 0.00268719 };
    const Double_t heWidthPars[5] = { -0.0105333, 0.129859, -0.0713011, 0.0169532, -0.00117928};
    const Double_t heGSigma = 0.0227879;
    const Double_t pMPV = EvalPol4(x, pMPVpars);
    const Double_t pWidth = EvalPol4(x, pWidthPars);
    const Double_t pSigma = TMath::Sqrt(pWidth*pWidth + pGSigma*pGSigma);
    const Double_t heMPV = EvalPol4(x, heMPVpars);
    const Double_t heWidth = EvalPol4(x, heWidthPars);
    const Double_t heSigma = TMath::Sqrt(heWidth*heWidth + heGSigma*heGSigma);
    qLow = pMPV - 3.0*pSigma;
    qHigh = heMPV + 6.0*heSigma;
}

void Load_data(const char* scenario="nominal") {
    const ChargeSysConfig cfg = MakeChargeSysConfig(scenario);
    // Correction variations modify MC only, so re-use nominal flight file.
    if (cfg.psdCorrectionStrength != 1.0 || cfg.stkCorrectionStrength != 1.0) {
        cerr << "For smearing scenarios, use the nominal flight histogram and vary MC only.\n";
        return;
    }

    TString basePath;
    TString hostname = gSystem->HostName();
    if (hostname.Contains("cnaf")) basePath = "/storage/gpfs_data/dampe/users";
    else if (hostname.Contains("le.infn")) basePath = "/nfs/argo/dampe";
    else {
        cerr << "WARNING: hostname non riconosciuto: " << hostname << endl;
        return;
    }

    Double_t BGO_E_corr, BGO_xtr;
    Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
    Double_t STK_chargeY_etaCorr[6], STK_chargeX_etaCorr[6], STK_vertexPrediction;
    Int_t BGO_HET;

    TChain *skim = new TChain("newtree");
    for (int y=2016; y<=2025; ++y) addYear(skim, basePath, y);
    cout << "Orbital Data Entries: " << skim->GetEntries() << endl;

    skim->SetBranchStatus("*", 0);
    skim->SetBranchStatus("BGO_HET", 1);
    skim->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions2", 1);
    skim->SetBranchStatus("PSD_ChargeY0", 1);
    skim->SetBranchStatus("PSD_ChargeY1", 1);
    skim->SetBranchStatus("PSD_ChargeX0", 1);
    skim->SetBranchStatus("PSD_ChargeX1", 1);
    skim->SetBranchStatus("BGO_xtr", 1);
    skim->SetBranchStatus("STK_chargeX_etaCorr", 1);
    skim->SetBranchStatus("STK_chargeY_etaCorr", 1);
    skim->SetBranchStatus("STK_vertexPrediction", 1);
    skim->SetCacheSize(200*1024*1024);
    skim->AddBranchToCache("*", kTRUE);

    skim->SetBranchAddress("BGO_HET", &BGO_HET);
    skim->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions2", &BGO_E_corr);
    skim->SetBranchAddress("BGO_xtr", &BGO_xtr);
    skim->SetBranchAddress("PSD_ChargeY0", &PSD_CY0);
    skim->SetBranchAddress("PSD_ChargeY1", &PSD_CY1);
    skim->SetBranchAddress("PSD_ChargeX0", &PSD_CX0);
    skim->SetBranchAddress("PSD_ChargeX1", &PSD_CX1);
    skim->SetBranchAddress("STK_chargeX_etaCorr", STK_chargeX_etaCorr);
    skim->SetBranchAddress("STK_chargeY_etaCorr", STK_chargeY_etaCorr);
    skim->SetBranchAddress("STK_vertexPrediction", &STK_vertexPrediction);

    const int nbd=5, ndec=6, noe=nbd*ndec;
    const Float_t arg1 = 1./Float_t(nbd);
    Float_t Ebin[noe+1];
    Ebin[0] = 10.;
    for (int j=1; j<noe+1; ++j) Ebin[j] = Ebin[j-1]*TMath::Power(10.,arg1);

    gSystem->mkdir("ROOT_FILES", kTRUE);
    const TString outputName = TString::Format("ROOT_FILES/PHe_skim_Orb120Month_5binperdecade_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_%s.root", cfg.tag.c_str());
    TFile fout(outputName, "RECREATE");
    if (fout.IsZombie()) { cerr << "Cannot create output: " << outputName << endl; return; }

    TH1D *h1SelBGO_orb = new TH1D("h1SelBGO_orb", "Selected(E_bgo) orbital", noe, Ebin);
    h1SelBGO_orb->Sumw2();
    TH1D *h1SelPSD_orb = new TH1D("h1SelPSD_orb", "Selected via PSD;E_{BGO} [GeV];Events", noe, Ebin);
    TH1D *h1SelSTK_orb = new TH1D("h1SelSTK_orb", "Selected via STK;E_{BGO} [GeV];Events", noe, Ebin);
    h1SelPSD_orb->Sumw2();
    h1SelSTK_orb->Sumw2();

    const Double_t minPSDSignal = 0.2;
    const Double_t dQmin=-0.3, dQmax=0.7;
    const Double_t vertexCut=0.7;
    const Double_t stkMin=25., stkMax=450.;
    Double_t stkLow=stkMin, stkHigh=stkMax;
    VaryChargeWindow(stkLow,stkHigh,cfg.stkWidthFraction);
    cout << "Scenario: " << cfg.tag << " | STK [" << stkLow << ", " << stkHigh << "]\n";

    const Long64_t nEntries = skim->GetEntries();
    for (Long64_t i=0; i<nEntries; ++i) {
        if (i%1000000==0) cout << "  Processing entry " << i << " / " << nEntries << endl;
        skim->GetEntry(i);

        if (BGO_HET<=0 || BGO_E_corr<=20.) continue;
        if ((PSD_CY0<=0. && PSD_CY1<=0.) || (PSD_CX0<=0. && PSD_CX1<=0.)) continue;
        if (BGO_xtr<=12.) continue;

        std::vector<Double_t> psdvec;
        const Double_t PSDCharges[4] = {PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1};
        for (int ilay=0; ilay<4; ++ilay)
            if (PSDCharges[ilay] > minPSDSignal) psdvec.push_back(PSDCharges[ilay]);
        Int_t nPSDlayers=0;
        const Double_t psdCharge=GetPSDProgressiveCharge(psdvec,dQmin,dQmax,nPSDlayers);
        const Double_t stkCharge=GetSTKLayerSignal(STK_chargeX_etaCorr[0],STK_chargeY_etaCorr[0]);

        if (STK_vertexPrediction < vertexCut) {
            if (psdCharge < 0.) continue;
            Double_t psdLow, psdHigh;
            GetPSDChargeLimits(BGO_E_corr, psdLow, psdHigh);
            VaryChargeWindow(psdLow, psdHigh, cfg.psdWidthFraction);
            if (psdCharge < psdLow || psdCharge > psdHigh) continue;
            h1SelPSD_orb->Fill(BGO_E_corr);
        } else {
            if (stkCharge < stkLow || stkCharge > stkHigh) continue;
            h1SelSTK_orb->Fill(BGO_E_corr);
        }
        h1SelBGO_orb->Fill(BGO_E_corr);
    }

    fout.cd();
    h1SelBGO_orb->Write();
    h1SelPSD_orb->Write();
    h1SelSTK_orb->Write();
    TNamed("ChargeScenario", cfg.tag.c_str()).Write();
    fout.Close();
    cout << "Wrote " << outputName << endl;
}
