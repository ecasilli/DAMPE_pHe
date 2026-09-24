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

const PSDSmearingPars pPSDpars = {
    // MPV DATA
    { 0.988625, 0.0166979, 0.0136979, -0.0114983, 0.00313886 },
    // WIDTH DATA
    {-0.0152013, 0.0753179, -0.0344289,  0.00656681, 0.000105438 },
    // sigma Gauss DATA
    2.02488e-08,

    // MPV MC
    { 0.885521, 0.192799, -0.115362,   0.0313466,-0.00102642 },
    // WIDTH MC
    { 0.0362995, 0.00614734, 0.00140585, 0.00113887, 0.000375804 },
    // sigma Gauss MC
    0.031548
};

const PSDSmearingPars HePSDpars = {
    // MPV DATA
    { 1.97763, 0.0402723, 0.00636717, -0.00967592, 0.00268719 },
    // WIDTH DATA
    { -0.0105333, 0.129859, -0.0713011, 0.0169532, -0.00117928 },
    // sigma Gauss DATA
    0.0227879,

    // MPV MC
    //{ 2.4225, -0.653307, 0.418625, -0.108939,   0.0118722 },
    { 2.17123,-0.236375, 0.16997,  -0.0457528,  0.00609007 },
    // WIDTH MC
    //{  0.236337, -0.345337,  0.231388, -0.064774,   0.00709489 },
    {  0.0471634,-0.0287278, 0.0407409,-0.0157212,  0.00254262 },
    // sigma Gauss MC
    //0.0838493
    0.0835585
};

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
const Double_t pMPVcorr[5] = { 1.07996,-0.146468, 0.119006,-0.0403108, 0.00593183};
// PROTON MC corrected WIDTH
const Double_t pWidthCorr[5] = { 0.00440895, 0.0276405,-0.00368363, -0.00124546, 0.000830814};
const Double_t pGSigmaCorr = 0.000224763;

// HELIUM MC corrected MPV
//const Double_t heMPVcorr[5] = {2.08278, -0.10484, 0.0723122,-0.0183575, 0.00239775};
const Double_t heMPVcorr[5] = {2.09943, -0.177259, 0.141829, -0.0437065, 0.00563705};
// HELIUM MC corrected WIDTH
//const Double_t heWidthCorr[5] = { -0.0530575, 0.161719, -0.0782506, 0.0159716, -0.000868304};
//const Double_t heGSigmaCorr = 0.0368612;
const Double_t heWidthCorr[5] = { 0.0410602, 0.018532, -0.0020576, -0.000574209, 0.000384744};
const Double_t heGSigmaCorr = 0.0307387;

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
    qHigh = heMPV + 6.5 * heSigma;
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

// FIND TRUE-ENERGY BIN
// Equivalent to:
// MC_EnergyT >= Ebin[j] && MC_EnergyT < Ebin[j+1]

Int_t FindEnergyBin(Double_t E, const Double_t *Ebin, Int_t nBins) {
    if (E < Ebin[0] || E >= Ebin[nBins])
        return -1;

    for (Int_t j = 0; j < nBins; ++j) {
        if (E >= Ebin[j] && E <  Ebin[j+1])
            return j;
    }

    return -1;
}

void ProcessMCSpecies(TChain **chains,
                      Int_t nsets,
                      const Double_t *sampleNorm,
                      Bool_t isHelium,
                      TH2D *hAll,
                      TH2D *hCut00,
                      TH2D *hCut01,
                      //TH2D *hCut06,
                      TH2D *hSpCut,
                      TH2D *hFinal,
                      const Double_t *Ebin,
                      const Double_t *truthBinNorm,
                      Int_t nBins) {
    
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
        Int_t BGO_HET;
        Double_t STK_vertexPrediction;

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

            // TRUE ENERGY BIN
            const Int_t jTrue = FindEnergyBin( MC_EnergyT, Ebin, nBins);
            if (jTrue < 0)
                continue;

            // ================================================================
            // EVENT WEIGHT
            // Riproduce: weight = normalization * MC_EnergyT^(-1.7)

            Double_t weight = sampleNorm[iset] * truthBinNorm[jTrue] * TMath::Power( MC_EnergyT, -1.7);

            if (isHelium)
                weight *= GeoCorr;

            if (MC_EnergyT <= 20.) continue;

            hAll->Fill(MC_EnergyT, BGO_E_corr, weight);

            // cut00 = cc204s * Trig_HEP
            if (BGO_HET <= 0)   continue;
            if (BGO_E_corr <= 20.) continue;

            hCut00->Fill(MC_EnergyT, BGO_E_corr, weight);

            // cut01: almeno un segnale su entrambi i piani di PSD
            if ((PSD_CY0 <= 0. && PSD_CY1 <= 0.) ||
                (PSD_CX0 <= 0. && PSD_CX1 <= 0.)) continue;

            hCut01->Fill(MC_EnergyT, BGO_E_corr, weight);


            // cut06: accettanza STK-PSD
            //if (TMath::Abs(STK_Y) >= 400.) continue;
            //if (TMath::Abs(STK_X) >= 400.) continue;
            //hCut06->Fill( MC_EnergyT, BGO_E_corr, weight);

            // SpCut: taglio per elettroni
            if (BGO_xtr <= 12.) continue;
            hSpCut->Fill(MC_EnergyT, BGO_E_corr, weight);

            // NEW PROGRESSIVE CHARGE: asimmetrica, come nel caso adottato da Paul
            std::vector<Double_t> psdvec;

            Double_t PSDCharges[4] = { PSD_CY0, PSD_CY1, PSD_CX0, PSD_CX1};
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

            // RESPONSE MATRIX
            //
            // x = TRUE energy
            // y = reconstructed BGO energy

            hFinal->Fill(MC_EnergyT, BGO_E_corr, weight);

        } // event loop 

        cout << "Finished " << species << " sample " << iset << endl;

    } // MC loop
}

void Load_MC_He_p_5PeV_spectrum_PSD_STK_comb(){

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

    // =============================
    // ENERGY BINNING

    const int nbd = 5;
    const int ndec = 6;
    Double_t arg1 = 1./Double_t(nbd);
    const int noe = nbd*ndec;
    cout << "Number of energy bins: " << noe << endl;
    Double_t e0 = 10.;
    Double_t Ebin[noe+1];
    Ebin[0] = e0;
    for (int j = 1; j < noe+1; j++) {
        Ebin[j] = Ebin[j-1]*TMath::Power(10., arg1);
    }

    // ==========================================
    // MC normalization factors
    // proton: wP  * Etrue^-1.7
    // helium: wHe * Etrue^-1.7 * GeoCorr
    
    Double_t normP[nsetP] = {
        TMath::Log(10.) / 1740700000.,
        TMath::Log(10.) / 519070000.,
        TMath::Log(10.) / 201432000.,
        TMath::Log(10.) / 138744900.,
        TMath::Log(10.) / 19431600.,
        TMath::Log(5.)  / 5094200.
    };

    Double_t normHe[nsetHe] = {
        (TMath::Log(10.) / 20380000.),
        (TMath::Log(10.) / 19556000.),
        (TMath::Log(10.) / 11346000.),
        (TMath::Log(10.) / 17694000.),
        (TMath::Log(5.) / 10305100.),
        (TMath::Log(2.) / 10281395.),
        (TMath::Log(5.) / 5064900.)
    };

    // TRUE-ENERGY BIN NORMALIZATION

    Double_t truthBinNorm[noe];

    for (Int_t j = 0; j < noe; ++j) {
        const Double_t E0 = Ebin[j];
        const Double_t E1 = Ebin[j+1];

        truthBinNorm[j] = 1.7 * TMath::Log(E1 / E0) / (TMath::Power(E0, -1.7) - TMath::Power(E1,-1.7));
    }

    // =======================================
    // Output

    TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_5binperdecade_3sigmaLow_6e5sigmaUp_PSDprogr_STKcharge450_comb_STKvert0e7_24sett26.root", "RECREATE");
    fout->cd();

    // ===============================================
    // Histograms
    TH2D *h2Ntrig_wgt = new TH2D("h2Ntrig_wgt", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    h2Ntrig_wgt->Sumw2();

    TH2D *h2Ntrig_wgt_all = new TH2D("h2Ntrig_wgt_all", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    h2Ntrig_wgt_all->Sumw2();
    TH2D *h2Ntrig_wgt_cut00 = new TH2D("h2Ntrig_wgt_cut00", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    h2Ntrig_wgt_cut00->Sumw2();
    TH2D *h2Ntrig_wgt_cut01 = new TH2D("h2Ntrig_wgt_cut01", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    h2Ntrig_wgt_cut01->Sumw2();
    //TH2D *h2Ntrig_wgt_cut06 = new TH2D("h2Ntrig_wgt_cut06", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    //h2Ntrig_wgt_cut06->Sumw2();
    TH2D *h2Ntrig_wgt_SpCut = new TH2D("h2Ntrig_wgt_SpCut", "Selected MC; MC true energy [GeV]; BGO reconstructed energy [GeV]",noe, Ebin, noe, Ebin);
    h2Ntrig_wgt_SpCut->Sumw2();

    // ======================================
    // Process PROTON

    cout << endl;
    cout << "Starting PROTON processing..." << endl;
    ProcessMCSpecies(
        sk_p,
        nsetP,
        normP,
        kFALSE,

        h2Ntrig_wgt_all,
        h2Ntrig_wgt_cut00,
        h2Ntrig_wgt_cut01,
        //h2Ntrig_wgt_cut06,
        h2Ntrig_wgt_SpCut,
        h2Ntrig_wgt,

        Ebin,
        truthBinNorm,
        noe
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

        h2Ntrig_wgt_all,
        h2Ntrig_wgt_cut00,
        h2Ntrig_wgt_cut01,
        //h2Ntrig_wgt_cut06,
        h2Ntrig_wgt_SpCut,
        h2Ntrig_wgt,

        Ebin,
        truthBinNorm,
        noe
    );

    cout << "HELIUM processing finished." << endl;

    // ======================================
    // SAVE
    
    fout->cd();

    h2Ntrig_wgt->Write();
    h2Ntrig_wgt_all->Write();
    h2Ntrig_wgt_cut00->Write();
    h2Ntrig_wgt_cut01->Write();
    //h2Ntrig_wgt_cut06->Write();
    h2Ntrig_wgt_SpCut->Write();

    fout->Close();


    cout << endl;
    cout << "========================================" << endl;
    cout << "End script." << endl;
    cout << "========================================" << endl;


}
