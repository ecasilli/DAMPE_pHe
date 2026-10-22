// Variante a 3 bin/decade in Etrue e E_BGO; la somma 10-100 TeV si fa SOLO nel plot.
// TEST: composizione COSTANTE p/He, default 70/30 in ogni bin Etrue.
// Questa macro usa i singoli istogrammi DI GENERAZIONE (pre-trigger):
//   h_energy_truth_3bins_cut_0
// per ricostruire Gp, GHe con gli stessi sampleNorm, Etrue^-1.7 e GeoCorr
// utilizzati nell'efficienza. NON utilizza file SBPL ne' h1Ngen esterni.
// Istogrammi di generazione e bin Etrue dell'analisi devono avere gli STESSI bordi.
// Nessuna interpolazione/redistribuzione tra bin: per E^-1.7 si usa la media
// analitica nel singolo bin, assumendo generazione uniforme in log(E).
// Nessuna lettura dei TTree per il calcolo dei G: gli skim potrebbero essere
// gia' stati selezionati e NON rappresentare tutti gli eventi generati.
// L'efficienza viene sempre riempita in E_BGO (default 3 bin/decade).

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TH1.h"
#include "TH1D.h"

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
    { 1.92403, 0.1541,  -0.082833,  0.0208528, -0.0011579 },
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
const Double_t heMPVcorr[5] = {2.08278, -0.10484, 0.0723122,-0.0183575, 0.00239775};
// HELIUM MC corrected WIDTH
const Double_t heWidthCorr[5] = { -0.0530575, 0.161719, -0.0782506, 0.0159716, -0.000868304};
const Double_t heGSigmaCorr = 0.0368612;

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


// ============================================================================
// Normalizzazione specifica di OGNI sample e spettro di generazione.
// Ngen e intervallo energetico corrispondono alla lista dei sample MC.
// Il file ROOT di ogni sample contiene l'istogramma di GENERAZIONE
// 'h_energy_truth_3bins_cut_0' il cui integrale deve essere circa Ngen.
// ============================================================================
struct MCGenerationSample {
    Double_t Emin;
    Double_t Emax;
    Double_t Ngen;
};

Double_t SampleNorm(const MCGenerationSample& sample) {
    return TMath::Log(sample.Emax / sample.Emin) / sample.Ngen;
}

// h_energy_truth_3bins_cut_0 riporta conteggi NON pesati.
// Se la generazione e' piatta in log(E), in un bin sorgente [a,b]
// l'aspettazione di E^-1.7 integrata su [x0,x1] e':
//   nBin * (x0^-1.7 - x1^-1.7) / (1.7 * ln(b/a))
// dove x0,x1 e' l'intersezione tra bin sorgente, range del sample e bin Etrue.
// Se un bin attraversa Emin/Emax di un sample, si integra solo
// nella porzione effettivamente generata [a,b].
Bool_t AccumulateGeneratedSpectrum(
    const std::vector<TString>& sampleFiles,
    const std::vector<MCGenerationSample>& samples,
    Bool_t isHelium,
    const std::vector<Double_t>& trueEdges,
    std::vector<Double_t>& G) {

    if (sampleFiles.size() != samples.size()) {
        cerr << "ERROR: numero di sample e range incompatibili" << endl;
        return kFALSE;
    }
    const Double_t geoCorr = 1. / (1.38 * 1.38);
    const Double_t geo = isHelium ? geoCorr : 1.;
    const Int_t nTrue = static_cast<Int_t>(trueEdges.size()) - 1;
    G.assign(nTrue, 0.);

    for (size_t iset = 0; iset < samples.size(); ++iset) {
        const MCGenerationSample& s = samples[iset];
        if (!(s.Emin > 0. && s.Emax > s.Emin && s.Ngen > 0.)) {
            cerr << "ERROR: intervallo o Ngen non valido al sample " << iset << endl;
            return kFALSE;
        }
        TFile *fin = TFile::Open(sampleFiles[iset], "READ");
        if (!fin || fin->IsZombie()) {
            cerr << "ERROR: apertura file generazione " << sampleFiles[iset] << endl;
            if (fin) { fin->Close(); delete fin; }
            return kFALSE;
        }
        TH1 *h = dynamic_cast<TH1*>(fin->Get("h_energy_truth_8bins_cut_0"));
        if (!h) {
            cerr << "ERROR: istogramma h_energy_truth_3bins_cut_0 assente in "
                 << sampleFiles[iset] << endl;
            fin->Close(); delete fin;
            return kFALSE;
        }

        // Istogramma di generazione gia' a 3 bin/decade:
        // il suo binning DEVE coincidere con quello usato nei pesi.
        if (h->GetNbinsX() != nTrue) {
            cerr << "ERROR: istogramma di generazione con " << h->GetNbinsX()
                 << " bin anziche' " << nTrue << ": " << sampleFiles[iset] << endl;
            fin->Close(); delete fin;
            return kFALSE;
        }
        for (Int_t j = 0; j < nTrue; ++j) {
            const Double_t lo = h->GetXaxis()->GetBinLowEdge(j+1);
            const Double_t hi = h->GetXaxis()->GetBinUpEdge(j+1);
            if (std::fabs(lo / trueEdges[j] - 1.) > 1.e-4 ||
                std::fabs(hi / trueEdges[j+1] - 1.) > 1.e-4) {
                cerr << "ERROR: bordi Etrue incompatibili nel bin " << j+1
                     << " del file " << sampleFiles[iset] << endl;
                fin->Close(); delete fin;
                return kFALSE;
            }
        }

        const Double_t nHist = h->Integral(0, h->GetNbinsX()+1);
        const Double_t relDiff = TMath::Abs(nHist / s.Ngen - 1.);
        cout << (isHelium ? "He" : "p") << " gen sample " << iset
             << ": Ngen=" << s.Ngen << ", Nhist=" << nHist
             << ", relative difference=" << relDiff << endl;
        // Se h non contiene TUTTI gli eventi generati, la composizione
        // ricavata non e' attendibile. Non procediamo silenziosamente.
        if (relDiff > 0.01) {
            cerr << "ERROR: l'istogramma di generazione differisce da Ngen"
                 << " di oltre l'1%: " << sampleFiles[iset] << endl;
            fin->Close(); delete fin;
            return kFALSE;
        }
        if (relDiff > 1.e-5)
            cout << "WARNING: piccolo scostamento Nhist/Ngen nel sample " << iset << endl;
        if (h->GetBinContent(0) != 0. ||
            h->GetBinContent(h->GetNbinsX()+1) != 0.) {
            cerr << "ERROR: eventi in underflow/overflow dell'istogramma di generazione"
                 << " per " << sampleFiles[iset] << endl;
            fin->Close(); delete fin;
            return kFALSE;
        }

        const Double_t norm = SampleNorm(s);
        for (Int_t ib = 1; ib <= h->GetNbinsX(); ++ib) {
            const Double_t nGenerated = h->GetBinContent(ib);
            if (nGenerated == 0.) continue;
            if (nGenerated < 0. || !std::isfinite(nGenerated)) {
                cerr << "ERROR: contenuto di generazione non valido al bin " << ib << endl;
                fin->Close(); delete fin;
                return kFALSE;
            }

            const Double_t rawLow = h->GetXaxis()->GetBinLowEdge(ib);
            const Double_t rawHigh = h->GetXaxis()->GetBinUpEdge(ib);
            const Double_t low = std::max(rawLow, s.Emin);
            const Double_t high = std::min(rawHigh, s.Emax);
            if (!(high > low && low > 0.)) {
                cerr << "ERROR: conteggi non nulli fuori dal range del sample "
                     << iset << ", bin " << ib << endl;
                fin->Close(); delete fin;
                return kFALSE;
            }
            // Lo stesso bin ib e' il bin Etrue della composizione.
            // Il contributo del sample e' limitato al suo range generato.
            const Double_t meanEminus17 =
                (TMath::Power(low, -1.7) - TMath::Power(high, -1.7))
                / (1.7 * TMath::Log(high / low));
            G[ib-1] += nGenerated * norm * geo * meanEminus17;
        }
        fin->Close();
        delete fin;
    }
    return kTRUE;
}

Bool_t BuildFixedMixWeights(
    const std::vector<TString>& filesP,
    const std::vector<MCGenerationSample>& samplesP,
    const std::vector<TString>& filesHe,
    const std::vector<MCGenerationSample>& samplesHe,
    const std::vector<Double_t>& trueEdges,
    Double_t targetP,
    std::vector<Double_t>& GP,
    std::vector<Double_t>& GHe,
    std::vector<Double_t>& genFracP,
    std::vector<Double_t>& genFracHe,
    std::vector<Double_t>& weightP,
    std::vector<Double_t>& weightHe) {

    if (!AccumulateGeneratedSpectrum(filesP,samplesP,kFALSE,trueEdges,GP) ||
        !AccumulateGeneratedSpectrum(filesHe,samplesHe,kTRUE,trueEdges,GHe))
        return kFALSE;

    const Int_t nTrue = static_cast<Int_t>(trueEdges.size()) - 1;
    genFracP.assign(nTrue,0.);
    genFracHe.assign(nTrue,0.);
    weightP.assign(nTrue,0.);
    weightHe.assign(nTrue,0.);
    const Double_t targetHe = 1. - targetP;
    cout << "\n# Etrue_GeV Gp GHe gP gHe weightP weightHe checkTotal" << endl;
    for (Int_t j = 0; j < nTrue; ++j) {
        const Double_t gp0=GP[j], ghe0=GHe[j];
        const Double_t total=gp0+ghe0;
        if (total <= 0.) {
            cout << "Etrue bin " << j+1 << " senza eventi generati (p+He); pesi=0\n";
            continue;
        }
        if (gp0 <= 0. || ghe0 <= 0.) {
            cerr << "ERROR: una delle due specie non ha generazione nel bin Etrue "
                 << j+1 << "; impossibile imporre 70/30" << endl;
            return kFALSE;
        }
        genFracP[j] = gp0 / total;
        genFracHe[j] = ghe0 / total;
        weightP[j] = targetP / genFracP[j];
        weightHe[j] = targetHe / genFracHe[j];
        const Double_t closure=(weightP[j]*gp0+weightHe[j]*ghe0)/total;
        const Double_t E = std::sqrt(trueEdges[j]*trueEdges[j+1]);
        cout << E << " " << gp0 << " " << ghe0 << " " << genFracP[j]
             << " " << genFracHe[j] << " " << weightP[j] << " "
             << weightHe[j] << " " << closure << endl;
        if (TMath::Abs(closure-1.) > 1.e-10) {
            cerr << "ERROR: verifica di chiusura fallita nel bin " << j+1 << endl;
            return kFALSE;
        }
    }
    return kTRUE;
}


void ProcessMCSpecies(TChain **chains,
                      Int_t nsets,
                      const Double_t *sampleNorm,
                      Bool_t isHelium,
                      TH1D *h1_Unb,
                      TH1D *h1Nobs,
                      TH1D *h1Nobs_And,
                      const Double_t *EbinTrue,
                      Int_t nBinsTrue,
                      const Double_t *relativeMixWeight) {
    
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

            // La composizione e' definita in energia VERA, non in energia BGO.
            const Int_t jTrue = FindEnergyBin(MC_EnergyT, EbinTrue, nBinsTrue);
            if (jTrue < 0)
                continue;

            // Normalizzazione originale MC e spettro Etrue^-1.7.
            Double_t weight = sampleNorm[iset] * TMath::Power(MC_EnergyT, -1.7);

            if (isHelium)
                weight *= GeoCorr;

            // Ripesatura RELATIVA: w_s(Etrue) = f_s(costante) / g_s(MC).
            // Questa forma preserva la normalizzazione p+He nel singolo bin vero.
            weight *= relativeMixWeight[jTrue];

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

            // Tutti i tagli offline sono stati superati.
            const Bool_t passUBT = (BGO_UBT > 0);
            const Bool_t passHET = (BGO_HET > 0);

            if (passUBT)
                h1_Unb->Fill(BGO_E_corr, weight);

            if (passHET)
                h1Nobs->Fill(BGO_E_corr, weight);

            if (passUBT && passHET)
                h1Nobs_And->Fill(BGO_E_corr, weight);

        } // event loop 

        cout << "Finished " << species << " sample " << iset << endl;

    } // MC loop
}

void HET_efficiency_pHe_MC_fixedMix(Int_t nbdBGO = 8,
                                    Double_t fractionP = 0.70,
                                    Int_t nbdTrue = 8){

    if (!(fractionP > 0. && fractionP < 1.)) {
        cerr << "ERROR: fractionP deve essere compresa strettamente tra 0 e 1" << endl;
        return;
    }
    if (nbdTrue != 3) {
        cerr << "ERROR: questa macro legge h_energy_truth_3bins_cut_0:"
             << " controllare di aver modificato" << endl;
        //return;
    }
    const Double_t fractionHe = 1. - fractionP;

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

    // Binning vero: usato SOLO per i pesi p/He; default 3 bin/decade.
    // Binning BGO: usato per la vera efficienza HET (default 3 bin/decade).
    if (nbdBGO <= 0) {
        cerr << "ERROR: nbdBGO deve essere positivo" << endl;
        return;
    }
    const int ndec = 6;  // 10 GeV -> 10 PeV
    const int nTrue = nbdTrue * ndec;
    const int nBGO  = nbdBGO  * ndec;
    std::vector<Double_t> EbinTrue(nTrue + 1);
    std::vector<Double_t> EbinBGO(nBGO + 1);
    for (int j = 0; j <= nTrue; ++j)
        EbinTrue[j] = 10. * TMath::Power(10., Double_t(j) / nbdTrue);
    for (int j = 0; j <= nBGO; ++j)
        EbinBGO[j] = 10. * TMath::Power(10., Double_t(j) / nbdBGO);

    cout << "True-energy bins: " << nTrue << "  -  BGO-energy bins: " << nBGO << endl;
    cout << "FIXED MIX: p=" << fractionP << " He=" << fractionHe << endl;

    // I sample sono pesati con i loro REALI range ed Ngen,
    // identici a quelli usati nella macro originale di generazione.
    const std::vector<MCGenerationSample> samplesP = {
        {1.e1, 1.e2, 1740700000.},
        {1.e2, 1.e3, 519070000.},
        {1.e3, 1.e4, 201432000.},
        {1.e4, 1.e5, 107980000.},
        {1.e5, 1.e6, 19431600.},
        {1.e6, 5.e6, 5094200.}
    };
    const std::vector<MCGenerationSample> samplesHe = {
        {1.e1, 1.e2, 20380000.},
        {1.e2, 1.e3, 19556000.},
        {1.e3, 1.e4, 11346000.},
        {1.e4, 1.e5, 17694000.},
        {1.e5, 5.e5, 10305100.},
        {5.e5, 1.e6, 10281395.},
        {1.e6, 5.e6, 5064900.}
    };
    if (samplesP.size() != filesP.size() || samplesHe.size() != filesHe.size()) {
        cerr << "ERROR: sample e file non hanno lo stesso numero di elementi" << endl;
        return;
    }

    std::vector<Double_t> normP(nsetP), normHe(nsetHe);
    for (Int_t i=0; i<nsetP; ++i) normP[i]=SampleNorm(samplesP[i]);
    for (Int_t i=0; i<nsetHe; ++i) normHe[i]=SampleNorm(samplesHe[i]);

    // Calcolo della miscela MC sui CAMPIONI GENERATI, NON su eventi
    // gia' selezionati dalla skim. Stessa convenzione dei pesi dell'HET:
    // sampleNorm * Etrue^-1.7 * GeoCorr(He).
    std::vector<Double_t> Gp, GHe;
    std::vector<Double_t> genFracP, genFracHe;
    std::vector<Double_t> relativeWeightP, relativeWeightHe;
    if (!BuildFixedMixWeights(filesP, samplesP, filesHe, samplesHe,
                              EbinTrue, fractionP, Gp, GHe, genFracP, genFracHe,
                              relativeWeightP, relativeWeightHe))
        return;

    // Nessuna truthBinNorm: per l'efficienza rispetto a E_BGO
    // non vogliamo appiattire il peso relativo delle energie vere.

    // =======================================
    // Output

    const TString outputFile = TString::Format(
        "HET_eff_MC_pHe_fixedMix_p%g_He%g_%dbinBGO_%dbinTrue.root",
        100.*fractionP, 100.*fractionHe, nbdBGO, nbdTrue);
    TFile *fout = TFile::Open(outputFile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        cerr << "ERROR: impossibile creare output " << outputFile << endl;
        if (fout) { fout->Close(); delete fout; }
        return;
    }
    fout->cd();

    // ===============================================
    // Histograms


    TH1D *h1_Unb = new TH1D("h1_Unb", "Unbias-Trig (p+He); BGO deposited energy [GeV]; Events", nBGO, EbinBGO.data());
    h1_Unb->Sumw2();

    TH1D *h1Nobs = new TH1D("h1Nobs", "HE-Trig (p+He); BGO deposited energy [GeV]; Events", nBGO, EbinBGO.data());
    h1Nobs->Sumw2();

    TH1D *h1Nobs_And = new TH1D("h1Nobs_And", "HE*Unbias-Trig (p+He); BGO deposited energy [GeV]; Events", nBGO, EbinBGO.data());
    h1Nobs_And->Sumw2();

    TH1D *h1_HET_Eff = new TH1D("h1_HET_Eff", "HET-Eff vs Ebgo (p+He);BGO deposited energy [GeV];Efficiency", nBGO, EbinBGO.data());
    h1_HET_Eff->Sumw2();

    // Diagnostica dei pesi usati (tutti in energia VERA).
    TH1D *hGenWeightedP = new TH1D("h_gen_weighted_p", "Generated weighted p;MC true energy [GeV];G_{p}", nTrue, EbinTrue.data());
    TH1D *hGenWeightedHe = new TH1D("h_gen_weighted_He", "Generated weighted He;MC true energy [GeV];G_{He}", nTrue, EbinTrue.data());
    TH1D *hMixP = new TH1D("h_frac_target_p", "Fixed p fraction;MC true energy [GeV];f_{p}", nTrue, EbinTrue.data());
    TH1D *hMixHe = new TH1D("h_frac_target_He", "Fixed He fraction;MC true energy [GeV];f_{He}", nTrue, EbinTrue.data());
    TH1D *hGP = new TH1D("h_gen_frac_p", "MC p fraction;MC true energy [GeV];g_{p}", nTrue, EbinTrue.data());
    TH1D *hGHe = new TH1D("h_gen_frac_He", "MC He fraction;MC true energy [GeV];g_{He}", nTrue, EbinTrue.data());
    TH1D *hWP = new TH1D("h_mix_relative_weight_p", "Relative p weight;MC true energy [GeV];f_{p}/g_{p}", nTrue, EbinTrue.data());
    TH1D *hWHe = new TH1D("h_mix_relative_weight_He", "Relative He weight;MC true energy [GeV];f_{He}/g_{He}", nTrue, EbinTrue.data());
    for (Int_t j = 0; j < nTrue; ++j) {
        hGenWeightedP->SetBinContent(j+1, Gp[j]);
        hGenWeightedHe->SetBinContent(j+1, GHe[j]);
        hMixP->SetBinContent(j+1, fractionP);
        hMixHe->SetBinContent(j+1, fractionHe);
        hGP->SetBinContent(j+1, genFracP[j]);
        hGHe->SetBinContent(j+1, genFracHe[j]);
        hWP->SetBinContent(j+1, relativeWeightP[j]);
        hWHe->SetBinContent(j+1, relativeWeightHe[j]);
    }

    // ======================================
    // Process PROTON

    cout << endl;
    cout << "Starting PROTON processing..." << endl;
    ProcessMCSpecies(
        sk_p,
        nsetP,
        normP.data(),
        kFALSE,

        h1_Unb,
        h1Nobs,
        h1Nobs_And,
        EbinTrue.data(), nTrue,
        relativeWeightP.data()
    );

    cout << "PROTON processing finished." << endl;

    // ======================================
    // Process HELIUM

    cout << endl;
    cout << "Starting HELIUM processing..." << endl;
    ProcessMCSpecies(
        sk_he,
        nsetHe,
        normHe.data(),
        kTRUE,

        h1_Unb,
        h1Nobs,
        h1Nobs_And,
        EbinTrue.data(), nTrue,
        relativeWeightHe.data()
    );

    cout << "HELIUM processing finished." << endl;

    // Efficienza E_HET|UBT = N(HET && UBT) / N(UBT).
    // Con eventi pesati, i due conteggi sono correlati:
    // var(e) = [(1-e)^2 * sum_{pass}w^2 + e^2 * sum_{fail}w^2]
    //          / [sum_{UBT}w]^2.
    for (Int_t ib = 1; ib <= nBGO; ++ib) {
        const Double_t den = h1_Unb->GetBinContent(ib);
        const Double_t num = h1Nobs_And->GetBinContent(ib);
        if (den <= 0.) {
            h1_HET_Eff->SetBinContent(ib, 0.);
            h1_HET_Eff->SetBinError(ib, 0.);
            continue;
        }
        const Double_t eff = num / den;
        const Double_t sw2All = TMath::Power(h1_Unb->GetBinError(ib), 2);
        const Double_t sw2Pass = TMath::Power(h1Nobs_And->GetBinError(ib), 2);
        const Double_t sw2Fail = std::max(0., sw2All - sw2Pass);
        const Double_t variance =
            (TMath::Power(1.-eff, 2)*sw2Pass + eff*eff*sw2Fail) / (den*den);
        h1_HET_Eff->SetBinContent(ib, eff);
        h1_HET_Eff->SetBinError(ib, TMath::Sqrt(std::max(0., variance)));
    }

    // ======================================
    // SAVE
    
    fout->cd();

    h1_Unb->Write();
    h1Nobs->Write();
    h1Nobs_And->Write();
    h1_HET_Eff->Write();
    hGenWeightedP->Write();
    hGenWeightedHe->Write();
    hMixP->Write();
    hMixHe->Write();
    hGP->Write();
    hGHe->Write();
    hWP->Write();
    hWHe->Write();

    fout->Close();


    cout << endl;
    cout << "========================================" << endl;
    cout << "End script." << endl;
    cout << "========================================" << endl;


}
