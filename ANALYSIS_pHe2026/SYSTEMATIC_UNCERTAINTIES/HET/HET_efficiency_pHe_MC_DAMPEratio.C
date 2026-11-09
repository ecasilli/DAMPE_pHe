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
#include "TH2D.h"
#include "TH1D.h"
#include "TNamed.h"
#include <iostream>
#include "DAMPE_pHe_ratio.h"

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

Double_t EvalPol4LogE(Double_t E, const Double_t p[5]) {

    const Double_t x = TMath::Log10(E);

    return p[0]
         + p[1] * x
         + p[2] * x*x
         + p[3] * x*x*x
         + p[4] * x*x*x*x;
}

// CORREZIONE DI SMEARING PSD
struct PSDSmearingPars {
    Double_t mpvData[5];
    Double_t widthData[5];
    Double_t sigmaData;
    
    Double_t mpvMC[5];
    Double_t widthMC[5];
    Double_t sigmaMC;
};

PSDSmearingPars pPSDpars;

PSDSmearingPars HePSDpars;

Double_t ApplyPSDSmearingCorrection(Double_t chargeMC, Double_t BGOenergy, const PSDSmearingPars &p) {
    if (BGOenergy <= 0.) return chargeMC;

    const Double_t mpvData = EvalPol4LogE(BGOenergy, p.mpvData);
    const Double_t widthData = EvalPol4LogE(BGOenergy, p.widthData);
    const Double_t sigmaData = TMath::Sqrt(widthData*widthData + p.sigmaData*p.sigmaData);

    const Double_t mpvMC   = EvalPol4LogE(BGOenergy, p.mpvMC);
    const Double_t widthMC   = EvalPol4LogE(BGOenergy, p.widthMC);
    const Double_t sigmaMC   = TMath::Sqrt( widthMC*widthMC + p.sigmaMC*p.sigmaMC );

    if (sigmaMC <= 0.) return chargeMC;

    return (chargeMC - mpvMC) * (sigmaData / sigmaMC) + mpvData;
}

// PROTON MC corrected MPV
Double_t pMPVcorr[5];
// PROTON MC corrected WIDTH
Double_t pWidthCorr[5];
const Double_t pGSigmaCorr = 0.000224763;

// HELIUM MC corrected MPV
Double_t heMPVcorr[5];
// HELIUM MC corrected WIDTH
Double_t heWidthCorr[5];
const Double_t heGSigmaCorr = 0.0368612;

// Explicit assignments for ROOT 5 / CINT. Called before event processing.
void InitializePSDParameters_ROOT5() {
    pPSDpars.mpvData[0] = 0.988625;
    pPSDpars.mpvData[1] = 0.0166979;
    pPSDpars.mpvData[2] = 0.0136979;
    pPSDpars.mpvData[3] = -0.0114983;
    pPSDpars.mpvData[4] = 0.00313886;
    pPSDpars.widthData[0] = -0.0152013;
    pPSDpars.widthData[1] = 0.0753179;
    pPSDpars.widthData[2] = -0.0344289;
    pPSDpars.widthData[3] = 0.00656681;
    pPSDpars.widthData[4] = 0.000105438;
    pPSDpars.sigmaData = 2.02488e-08;
    pPSDpars.mpvMC[0] = 0.885521;
    pPSDpars.mpvMC[1] = 0.192799;
    pPSDpars.mpvMC[2] = -0.115362;
    pPSDpars.mpvMC[3] = 0.0313466;
    pPSDpars.mpvMC[4] = -0.00102642;
    pPSDpars.widthMC[0] = 0.0362995;
    pPSDpars.widthMC[1] = 0.00614734;
    pPSDpars.widthMC[2] = 0.00140585;
    pPSDpars.widthMC[3] = 0.00113887;
    pPSDpars.widthMC[4] = 0.000375804;
    pPSDpars.sigmaMC = 0.031548;
    HePSDpars.mpvData[0] = 1.92403;
    HePSDpars.mpvData[1] = 0.1541;
    HePSDpars.mpvData[2] = -0.082833;
    HePSDpars.mpvData[3] = 0.0208528;
    HePSDpars.mpvData[4] = -0.0011579;
    HePSDpars.widthData[0] = -0.0105333;
    HePSDpars.widthData[1] = 0.129859;
    HePSDpars.widthData[2] = -0.0713011;
    HePSDpars.widthData[3] = 0.0169532;
    HePSDpars.widthData[4] = -0.00117928;
    HePSDpars.sigmaData = 0.0227879;
    HePSDpars.mpvMC[0] = 2.17123;
    HePSDpars.mpvMC[1] = -0.236375;
    HePSDpars.mpvMC[2] = 0.16997;
    HePSDpars.mpvMC[3] = -0.0457528;
    HePSDpars.mpvMC[4] = 0.00609007;
    HePSDpars.widthMC[0] = 0.0471634;
    HePSDpars.widthMC[1] = -0.0287278;
    HePSDpars.widthMC[2] = 0.0407409;
    HePSDpars.widthMC[3] = -0.0157212;
    HePSDpars.widthMC[4] = 0.00254262;
    HePSDpars.sigmaMC = 0.0835585;
    pMPVcorr[0] = 1.07996;
    pMPVcorr[1] = -0.146468;
    pMPVcorr[2] = 0.119006;
    pMPVcorr[3] = -0.0403108;
    pMPVcorr[4] = 0.00593183;
    pWidthCorr[0] = 0.00440895;
    pWidthCorr[1] = 0.0276405;
    pWidthCorr[2] = -0.00368363;
    pWidthCorr[3] = -0.00124546;
    pWidthCorr[4] = 0.000830814;
    heMPVcorr[0] = 2.08278;
    heMPVcorr[1] = -0.10484;
    heMPVcorr[2] = 0.0723122;
    heMPVcorr[3] = -0.0183575;
    heMPVcorr[4] = 0.00239775;
    heWidthCorr[0] = -0.0530575;
    heWidthCorr[1] = 0.161719;
    heWidthCorr[2] = -0.0782506;
    heWidthCorr[3] = 0.0159716;
    heWidthCorr[4] = -0.000868304;
}

void GetPSDMCSelectionLimits(Double_t BGOenergy, Double_t &qLow, Double_t &qHigh) {
    // proton
    const Double_t pMPV = EvalPol4LogE(BGOenergy, pMPVcorr);
    const Double_t pWidth = EvalPol4LogE(BGOenergy, pWidthCorr);
    const Double_t pSigma = TMath::Sqrt(pWidth * pWidth + pGSigmaCorr * pGSigmaCorr);

    // helium
    const Double_t heMPV = EvalPol4LogE(BGOenergy, heMPVcorr);
    const Double_t heWidth = EvalPol4LogE(BGOenergy, heWidthCorr);
    const Double_t heSigma = TMath::Sqrt(heWidth * heWidth + heGSigmaCorr * heGSigmaCorr);

    qLow  = pMPV - 3.0 * pSigma;
    qHigh = heMPV + 6.0 * heSigma;
}

Double_t CorrectSTKProton(Double_t q) {
    const Double_t mpvMC = 57.3422;
    const Double_t mpvData = 58.8187;

    const Double_t sigmaMC  = TMath::Sqrt(5.47914 * 5.47914 + 5.65029 * 5.65029);
    const Double_t sigmaData =TMath::Sqrt(5.77431 * 5.77431 + 6.15165 * 6.15165);

    return
        (q - mpvMC) * (sigmaData / sigmaMC) + mpvData;
}


Double_t CorrectSTKHelium(Double_t q) {
    const Double_t mpvMC = 245.048;
    const Double_t mpvData = 250.715;

    const Double_t sigmaMC  = TMath::Sqrt(14.4618 * 14.4618 + 14.3083 * 14.3083);
    const Double_t sigmaData =TMath::Sqrt(18.4519 * 18.4519 + 12.5775 * 12.5775);

    return
        (q - mpvMC) * (sigmaData / sigmaMC) + mpvData;
}

void ProcessMCSpecies(TChain **chains,
                      Int_t nsets,
                      const Double_t *sampleNorm,
                      Bool_t isHelium,
                      TH1D *h1_Unb,
                      TH1D *h1Nobs,
                      TH1D *h1Nobs_And,
                      TH1D *h1_Unb_RatioOutside,
                      TH1D *h1Nobs_And_RatioOutside) {
    
    const Double_t minPSDSignal = 0.2;
    const Double_t dQmin =-0.3;
    const Double_t dQmax = 0.7;

    const Double_t vertexCut = 0.7;
    const Double_t stkMin = 25.;
    const Double_t stkMax = 450.;

    // helium geometrical normalization
    const Double_t GeoCorr = (2. * TMath::Pi() * TMath::Pi()) / (2. * TMath::Pi() * TMath::Pi() * 1.38 * 1.38);
    const TString species = isHelium ? "HELIUM" : "PROTON";

    // Loop sui sample MC
    for (Int_t iset = 0; iset < nsets; ++iset) {
        
        TChain *chain = chains[iset];

        Double_t MC_EnergyT, BGO_E_corr, BGO_xtr;
        Double_t PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1;
        Double_t STK_Y, STK_X;
        Double_t STK_chargeX_etaCorr[6];
        Double_t STK_chargeY_etaCorr[6];
        Int_t BGO_HET, BGO_UBT;
        Double_t STK_vertexPrediction;

        // =============================
        // ----- Branch selection ------
        
        chain->SetBranchStatus("*", 0);
        chain->SetBranchStatus("MC_EnergyT",                   1);
        chain->SetBranchStatus("BGO_HET",                      1);
        chain->SetBranchStatus("BGO_UBT",                      1);
        chain->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions2", 1);
        chain->SetBranchStatus("PSD_ChargeY0",                 1);
        chain->SetBranchStatus("PSD_ChargeY1",                 1);
        chain->SetBranchStatus("PSD_ChargeX0",                 1);
        chain->SetBranchStatus("PSD_ChargeX1",                 1);
        chain->SetBranchStatus("STKtrack_to_PSD_topY",         1);
        chain->SetBranchStatus("STKtrack_to_PSD_topX",         1);
        chain->SetBranchStatus("BGO_xtr",                      1);
        chain->SetBranchStatus("STK_chargeX_etaCorr",          1);
        chain->SetBranchStatus("STK_chargeY_etaCorr",          1);
        chain->SetBranchStatus("STK_vertexPrediction",         1);

        chain->SetCacheSize(200 * 1024 * 1024);
        chain->AddBranchToCache("*", kTRUE);

        // =============================
        // ----- Branch addresses ------

        chain->SetBranchAddress("MC_EnergyT",                   &MC_EnergyT);
        chain->SetBranchAddress("BGO_HET",                      &BGO_HET);
        chain->SetBranchAddress("BGO_UBT",                      &BGO_UBT);
        chain->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions2", &BGO_E_corr);
        chain->SetBranchAddress("PSD_ChargeY0",                 &PSD_CY0);
        chain->SetBranchAddress("PSD_ChargeY1",                 &PSD_CY1);
        chain->SetBranchAddress("PSD_ChargeX0",                 &PSD_CX0);
        chain->SetBranchAddress("PSD_ChargeX1",                 &PSD_CX1);
        chain->SetBranchAddress("STKtrack_to_PSD_topY",         &STK_Y);
        chain->SetBranchAddress("STKtrack_to_PSD_topX",         &STK_X);
        chain->SetBranchAddress("BGO_xtr",                      &BGO_xtr);
        chain->SetBranchAddress("STK_chargeX_etaCorr",           STK_chargeX_etaCorr);
        chain->SetBranchAddress("STK_chargeY_etaCorr",           STK_chargeY_etaCorr);
        chain->SetBranchAddress("STK_vertexPrediction",         &STK_vertexPrediction);

        // =============================
        // ----- Loop unico su tutti gli eventi ------

        Long64_t nEntries = chain->GetEntries();

        cout << endl;
        cout << "==============================================" << endl;
        cout << species << " sample " << iset << " : " << nEntries << " events" << endl;
        cout << "==============================================" << endl;

        for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {

            if (ievt % 1000000 == 0) cout << "  " << species << "[" << iset << "] " << ievt << " / " << nEntries << endl;

            chain->GetEntry(ievt);

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = sampleNorm[iset] * TMath::Power(MC_EnergyT, -1.7);

            if (isHelium)
                weight *= GeoCorr;

            // REQUIRED CONVENTION: MC_EnergyT is kinetic energy per nucleus
            // in GeV. If its definition differs, convert it here for the ratio.
            const Double_t TkinGeV = MC_EnergyT;
            weight *= DAMPEpHe::compositionWeight(TkinGeV, isHelium);
            const Bool_t ratioOutside = DAMPEpHe::outsideSupport(TkinGeV);

            //if (MC_EnergyT <= 20.) continue;

            // cut00 = cc204s * Trig_HEP
            //if (BGO_HET <= 0)   continue;
            if (BGO_E_corr <= 20.) continue;

            // cut01: almeno un segnale su entrambi i piani di PSD
            if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
                (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;


            // cut06: accettanza STK-PSD
            //if (TMath::Abs(STK_Y) >= 400.) continue;
            //if (TMath::Abs(STK_X) >= 400.) continue;

            // SpCut: taglio per elettroni
            if (BGO_xtr <= 12.) continue;

            // NEW PROGRESSIVE CHARGE: asimmetrica, come nel caso adottato da Paul
            std::vector<Double_t> psdvec;

            Double_t PSDCharges[4];
            PSDCharges[0] = PSD_CY0;
            PSDCharges[1] = PSD_CY1;
            PSDCharges[2] = PSD_CX0;
            PSDCharges[3] = PSD_CX1;
            for (int ilay=0; ilay < 4; ilay++) {
                if (PSDCharges[ilay] > minPSDSignal) // NOTA: non ho messo proprio 0 preciso
                    psdvec.push_back(PSDCharges[ilay]);
            }

            Int_t nPSDlayers = 0;

            const Double_t psdCharge = GetPSDProgressiveCharge(psdvec, dQmin, dQmax, nPSDlayers);

            // STK charge
            const Double_t stkCharge = GetSTKLayerSignal(STK_chargeX_etaCorr[0], STK_chargeY_etaCorr[0]);

            // VERTEX-DEPENDENT CHARGE SELECTION
            if (STK_vertexPrediction < vertexCut) {         // PSD branch

                if (psdCharge < 0.)
                    continue;

                const Double_t psdChargeCorrP = ApplyPSDSmearingCorrection(psdCharge, BGO_E_corr, pPSDpars);
                const Double_t psdChargeCorrHe = ApplyPSDSmearingCorrection(psdCharge, BGO_E_corr, HePSDpars);

                Double_t qLow;
                Double_t qHigh;

                GetPSDMCSelectionLimits(BGO_E_corr, qLow, qHigh);
                // MPV_p - 2.8 sigma_p
                if (psdChargeCorrP < qLow)
                    continue;

                // MPV_He + 6 sigma_He
                if (psdChargeCorrHe > qHigh)
                    continue;
            }
            else {                                          // STK branch

                if (stkCharge <= 0.)
                    continue;

                Double_t stkChargeCorr;

                if (isHelium) {
                    stkChargeCorr = CorrectSTKHelium(stkCharge);
                }
                else {
                    stkChargeCorr = CorrectSTKProton(stkCharge);
                }

                if (stkChargeCorr < stkMin)
                    continue;

                if (stkChargeCorr > stkMax)
                    continue;
            }

            // Tutti i tagli offline sono stati superati.
            const Bool_t passUBT = (BGO_UBT > 0);
            const Bool_t passHET = (BGO_HET > 0);

            if (passUBT)
                h1_Unb->Fill(BGO_E_corr, weight);

            if (passHET)
                h1Nobs->Fill(BGO_E_corr, weight);

            if (passUBT && passHET)
                h1Nobs_And->Fill(BGO_E_corr, weight);

            // Combined-species diagnostics: no p/He histogram separation.
            if (ratioOutside && passUBT)
                h1_Unb_RatioOutside->Fill(BGO_E_corr, weight);
            if (ratioOutside && passUBT && passHET)
                h1Nobs_And_RatioOutside->Fill(BGO_E_corr, weight);

        } // event loop 

        cout << "Finished " << species << " sample " << iset << endl;

    } // MC loop
}

void HET_efficiency_pHe_MC_DAMPEratio(Int_t binsPerDecade = 3){
    InitializePSDParameters_ROOT5();
    if (binsPerDecade < 1 || binsPerDecade > 100) {
        std::cerr << "Invalid binsPerDecade (allowed: 1..100)." << std::endl;
        return;
    }
    std::cout << "Composition: supplied DAMPE SEP2026 p/He central values."
              << "\nAssumption: MC_EnergyT = kinetic energy per nucleus [GeV]."
              << "\nRatio support: " << DAMPEpHe::minEnergyGeV() << " -- "
              << DAMPEpHe::maxEnergyGeV() << " GeV. Endpoint ratios frozen outside."
              << std::endl;

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
    //filesP.push_back(basePath + "/PROTON/Proton_10TeV_100TeV_FTFP_merged.root");
    filesP.push_back(basePath + "/PROTON/Proton_10TeV_100TeV_EPOSLHC_FTFP_merged.root");
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

    // =============================
    // ENERGY BINNING

    const int nbd = binsPerDecade;
    const int ndec = 6;
    Double_t arg1 = 1./Double_t(nbd);
    const int noe = nbd*ndec;
    cout << "Number of energy bins: " << noe << endl;
    Double_t e0 = 10.;
    std::vector<Double_t> Ebin(noe+1);
    Ebin[0] = e0;
    for (int j = 1; j < noe+1; j++) {
        Ebin[j] = Ebin[j-1]*TMath::Power(10., arg1);
    }

    // ==========================================
    // MC normalization factors
    // proton: wP  * Etrue^-1.7
    // helium: wHe * Etrue^-1.7 * GeoCorr
    
    Double_t normP[nsetP];
    normP[0] = TMath::Log(10.) / 1740700000.;
    normP[1] = TMath::Log(10.) / 519070000.;
    normP[2] = TMath::Log(10.) / 201432000.;
    normP[3] = TMath::Log(10.) / 107980000.;
    normP[4] = TMath::Log(10.) / 19431600.;
    normP[5] = TMath::Log(5.)  / 5094200.;

    Double_t normHe[nsetHe];
    normHe[0] = (TMath::Log(10.) / 20380000.);
    normHe[1] = (TMath::Log(10.) / 19556000.);
    normHe[2] = (TMath::Log(10.) / 11346000.);
    normHe[3] = (TMath::Log(10.) / 17694000.);
    normHe[4] = (TMath::Log(5.) / 10305100.);
    normHe[5] = (TMath::Log(2.) / 10281395.);
    normHe[6] = (TMath::Log(5.) / 5064900.);

    // =======================================
    // Output

    TString outputName = TString::Format("HET_eff_MC_pHe_DAMPEratio_%dbins.root", binsPerDecade);
    TFile *fout = new TFile(outputName, "RECREATE");
    if (fout->IsZombie()) {
        std::cerr << "Cannot create " << outputName << std::endl;
        return;
    }
    fout->cd();

    // ===============================================
    // Histograms


    TH1D *h1_Unb = new TH1D("h1_Unb", "Unbias-Trig (p+He); BGO deposited energy [GeV]; Events",noe, &Ebin[0]);
    h1_Unb->Sumw2();

    TH1D *h1Nobs = new TH1D("h1Nobs", "HE-Trig (p+He); BGO deposited energy [GeV]; Events",noe, &Ebin[0]);
    h1Nobs->Sumw2();

    TH1D *h1Nobs_And = new TH1D("h1Nobs_And", "HE*Unbias-Trig (p+He); BGO deposited energy [GeV]; Events",noe, &Ebin[0]);
    h1Nobs_And->Sumw2();

    TH1D *h1_HET_Eff = new TH1D("h1_HET_Eff", "HET-Eff vs Ebgo (p+He)", noe, &Ebin[0]);
    h1_HET_Eff->Sumw2();

    TH1D *h1_Unb_RatioOutside = new TH1D("h1_Unb_RatioOutside",
        "UNB with true energy outside ratio support; BGO deposited energy [GeV]; Weighted events", noe, &Ebin[0]);
    h1_Unb_RatioOutside->Sumw2();
    TH1D *h1Nobs_And_RatioOutside = new TH1D("h1Nobs_And_RatioOutside",
        "UNB and HET outside ratio support; BGO deposited energy [GeV]; Weighted events", noe, &Ebin[0]);
    h1Nobs_And_RatioOutside->Sumw2();

    // ======================================
    // Process PROTON

    cout << endl;
    cout << "Starting PROTON processing..." << endl;
    ProcessMCSpecies(
        sk_p,
        nsetP,
        normP,
        kFALSE,

        h1_Unb,
        h1Nobs,
        h1Nobs_And,
        h1_Unb_RatioOutside,
        h1Nobs_And_RatioOutside
    );

    cout << "PROTON processing finished." << endl;

    // ======================================
    // Process HELIUM

    cout << endl;
    cout << "Starting HELIUM processing..." << endl;
    ProcessMCSpecies(
        sk_he,
        nsetHe,
        normHe,
        kTRUE,

        h1_Unb,
        h1Nobs,
        h1Nobs_And,
        h1_Unb_RatioOutside,
        h1Nobs_And_RatioOutside
    );

    cout << "HELIUM processing finished." << endl;

    // Retained from the original macro. Use counts and TEfficiency in the
    // plotting macro for uncertainties; plain Divide ignores subset covariance.
    h1_HET_Eff->Divide(h1Nobs_And,h1_Unb);

    // ======================================
    // SAVE
    
    fout->cd();

    h1_Unb->Write();
    h1Nobs->Write();
    h1Nobs_And->Write();
    h1_HET_Eff->Write();
    h1_Unb_RatioOutside->Write();
    h1Nobs_And_RatioOutside->Write();
    TNamed compositionInfo("composition_model",
        "SEP2026 supplied flux tables; log-log interpolation; constant endpoint ratio; "
        "weights 2r/(1+r) for p, 2/(1+r) for He; total incident E^-2.7 retained; "
        "MC_EnergyT assumed kinetic energy per nucleus [GeV]");
    compositionInfo.Write();

    std::cout << "\nWeighted UNB fraction outside representative-energy support:" << std::endl;
    for (Int_t ib = 1; ib <= noe; ++ib) {
        const Double_t total = h1_Unb->GetBinContent(ib);
        if (total > 0.) {
            std::cout << h1_Unb->GetXaxis()->GetBinLowEdge(ib) << " -- "
                      << h1_Unb->GetXaxis()->GetBinUpEdge(ib) << " GeV: "
                      << h1_Unb_RatioOutside->GetBinContent(ib)/total << std::endl;
        }
    }

    fout->Close();


    cout << endl;
    cout << "========================================" << endl;
    cout << "End script." << endl;
    cout << "========================================" << endl;


}

