#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TPad.h"

using namespace std;


void efficienza_HET_DATA()
{

    // ============================================================
    // INPUT FILE
    // ============================================================

    TFile *fin = TFile::Open("HETeff_skim_Orb120Month_2e8sigmaLow_6sigmaUp_PSDprogr_STKcharge_comb_STKvert0e7.root","READ");

    if (!fin || fin->IsZombie()) {
        cout << "ERROR: cannot open input file" << endl;
        return;
    }

    // ============================================================
    // READ HISTOGRAMS
    // ============================================================

    TH1D *h1_Unb = nullptr;
    TH1D *h1Nobs_And = nullptr;

    fin->GetObject("h1_Unb", h1_Unb);
    fin->GetObject("h1Nobs_And", h1Nobs_And);

    if (!h1_Unb || !h1Nobs_And) {
        cout << "ERROR: input histograms not found" << endl;
        fin->Close();
        return;
    }


    // ============================================================
    // CHECK HISTOGRAM CONSISTENCY
    // ============================================================

    if (!TEfficiency::CheckConsistency( *h1Nobs_And, *h1_Unb)) {
        cout << "ERROR: inconsistent histograms" << endl;
        fin->Close();
        return;
    }


    // ============================================================
    // HET EFFICIENCY
    // ============================================================

    TEfficiency *pEff = new TEfficiency( *h1Nobs_And, *h1_Unb);
    pEff->SetName("HET_Eff_DATA");
    pEff->SetTitle("HET Efficiency DATA (p+He); BGO deposited energy [GeV]; HET Efficiency");

    // ============================================================
    // STATISTICAL TREATMENT
    // ============================================================

    // Clopper-Pearson confidence intervals
    pEff->SetStatisticOption(TEfficiency::kFCP);

    // Confidence level corresponding approximately to 1 sigma
    pEff->SetConfidenceLevel(0.682689492);

    // ============================================================
    // PRINT EFFICIENCY
    // ============================================================

    const Int_t nBins = h1_Unb->GetNbinsX();
    cout << endl;

    cout << "========================================" << endl;
    cout << "HET EFFICIENCY DATA" << endl;
    cout << "========================================" << endl;


    for (Int_t ibin = 1; ibin <= nBins; ++ibin) {

        const Double_t E_low = h1_Unb->GetXaxis()->GetBinLowEdge(ibin);
        const Double_t E_up = h1_Unb->GetXaxis()->GetBinUpEdge(ibin);
        const Double_t N_Unb = h1_Unb->GetBinContent(ibin);
        const Double_t N_And = h1Nobs_And->GetBinContent(ibin);

        // Skip empty bins
        if (N_Unb <= 0.) {
            cout << "Bin " << ibin << " : empty denominator" << endl;
            continue;
        }

        const Double_t efficiency = pEff->GetEfficiency(ibin);
        const Double_t errLow = pEff->GetEfficiencyErrorLow(ibin);
        const Double_t errUp = pEff->GetEfficiencyErrorUp(ibin);

        cout << "Bin " << ibin << " [" << E_low << ", " << E_up << "] GeV" << endl;
        cout << "  N_Unb = " << N_Unb << "  N_And = " << N_And << endl;
        cout << "  Efficiency = " << efficiency << " -" << errLow << " +" << errUp << endl;

    }


    // ============================================================
    // DRAW EFFICIENCY
    // ============================================================

    TCanvas *c1 = new TCanvas("c_HET_DATA", "HET Efficiency DATA", 800, 600);
    c1->cd();
    gPad->SetLogx();
    pEff->Draw("AP");
    c1->Update();

    // ============================================================
    // OUTPUT FILE
    // ============================================================

    TFile *fout = TFile::Open("HET_Efficiency_DATA_pHe.root", "RECREATE");

    if (!fout || fout->IsZombie()) {
        cout << "ERROR: cannot create output file" << endl;
        fin->Close();
        return;
    }

    fout->cd();
    pEff->Write();
    c1->Write();
    fout->Close();
    fin->Close();

    cout << endl;
    cout << "HET DATA efficiency calculation completed." << endl;

}
