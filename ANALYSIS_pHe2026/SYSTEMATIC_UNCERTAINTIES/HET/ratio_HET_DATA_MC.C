#include <iostream>
#include <cmath>
#include <algorithm>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"

using namespace std;


// ================================================================
// HET EFFICIENCY DATA-MC COMPARISON
// ================================================================

void ratio_HET_DATA_MC()
{

    // ============================================================
    // INPUT FILES
    // ============================================================

    TFile *fDATA = TFile::Open("HET_Efficiency_DATA_pHe.root", "READ");
    TFile *fMC = TFile::Open("HET_Efficiency_MC_pHe.root", "READ");

    if (!fDATA || fDATA->IsZombie() || !fMC   || fMC->IsZombie()) {
        cout << "ERROR: cannot open input files" << endl;
        return;
    }


    // ============================================================
    // READ TEFFICIENCY OBJECTS
    // ============================================================

    TEfficiency *effDATA = nullptr;
    TEfficiency *effMC   = nullptr;

    fDATA->GetObject("HET_Eff_DATA", effDATA);
    fMC->GetObject("HET_Eff_MC", effMC);

    if (!effDATA || !effMC) {
        cout << "ERROR: TEfficiency objects not found" << endl;
        return;
    }


    // ============================================================
    // CHECK BINNING
    // ============================================================

    const TH1 *hDATA = effDATA->GetTotalHistogram();
    const TH1 *hMC = effMC->GetTotalHistogram();

    if (!hDATA || !hMC) {
        cout << "ERROR: total histograms not found" << endl;
        return;
    }

    const Int_t nBins = hDATA->GetNbinsX();

    if (nBins != hMC->GetNbinsX()) {
        cout << "ERROR: different number of energy bins" << endl;
        return;
    }

    // Controllo anche che i bordi dei bin coincidano.
    for (Int_t i = 1; i <= nBins + 1; ++i) {
        const Double_t edgeDATA = hDATA->GetXaxis()->GetBinLowEdge(i);
        const Double_t edgeMC = hMC->GetXaxis()->GetBinLowEdge(i);
        const Double_t tolerance = 1.e-10 * std::max(1., std::max(std::abs(edgeDATA),std::abs(edgeMC)));

        if (std::abs(edgeDATA - edgeMC) > tolerance) {
            cout << "ERROR: different energy binning" << endl;
            return;
        }

    }


    // ============================================================
    // ENERGY RANGE
    // ============================================================

    const Double_t Emin = hDATA->GetXaxis()->GetBinLowEdge(1);
    const Double_t Emax =hDATA->GetXaxis()->GetBinUpEdge(nBins);

    // ============================================================
    // CREATE GRAPHS
    // ============================================================

    TGraphAsymmErrors *gDATA = new TGraphAsymmErrors();
    TGraphAsymmErrors *gMC = new TGraphAsymmErrors();
    TGraphAsymmErrors *gRatio = new TGraphAsymmErrors();

    gDATA->SetName("g_HET_Eff_DATA");
    gMC->SetName("g_HET_Eff_MC");
    gRatio->SetName("g_HET_Ratio_DATA_MC");

    // ============================================================
    // FILL GRAPHS
    // ============================================================

    Double_t minRatio = 1.;
    Double_t maxRatio = 1.;

    for (Int_t ibin = 1; ibin <= nBins; ++ibin) {

        // --------------------------------------------------------
        // ENERGY BIN
        // --------------------------------------------------------

        const Double_t E_low = hDATA->GetXaxis()->GetBinLowEdge(ibin);
        const Double_t E_up = hDATA->GetXaxis()->GetBinUpEdge(ibin);

        // Centro geometrico, appropriato per il binning logaritmico
        const Double_t E_center = std::sqrt(E_low * E_up);
        const Double_t E_errLow = E_center - E_low;
        const Double_t E_errUp = E_up - E_center;

        // --------------------------------------------------------
        // TOTAL EVENTS
        // --------------------------------------------------------

        const Double_t N_DATA = hDATA->GetBinContent(ibin);
        const Double_t N_MC = hMC->GetBinContent(ibin);

        // --------------------------------------------------------
        // DATA EFFICIENCY
        // --------------------------------------------------------

        Double_t epsDATA = 0.;
        Double_t errDATA_low = 0.;
        Double_t errDATA_up = 0.;

        if (N_DATA > 0.) {
            epsDATA = effDATA->GetEfficiency(ibin);
            errDATA_low = effDATA->GetEfficiencyErrorLow(ibin);
            errDATA_up = effDATA->GetEfficiencyErrorUp(ibin);

            const Int_t pointDATA = gDATA->GetN();

            gDATA->SetPoint(pointDATA, E_center, epsDATA);
            gDATA->SetPointError(pointDATA, E_errLow, E_errUp, errDATA_low, errDATA_up);
        }


        // --------------------------------------------------------
        // MC EFFICIENCY
        // --------------------------------------------------------

        Double_t epsMC = 0.;
        Double_t errMC_low = 0.;
        Double_t errMC_up = 0.;

        if (N_MC > 0.) {
            epsMC = effMC->GetEfficiency(ibin);
            errMC_low = effMC->GetEfficiencyErrorLow(ibin);
            errMC_up = effMC->GetEfficiencyErrorUp(ibin);

            const Int_t pointMC = gMC->GetN();
            gMC->SetPoint(pointMC, E_center, epsMC);
            gMC->SetPointError(pointMC, E_errLow, E_errUp, errMC_low, errMC_up);

        }

        // --------------------------------------------------------
        // DATA / MC RATIO
        // --------------------------------------------------------

        // Il rapporto viene calcolato soltanto nei bin
        // con denominatore DATA non nullo e efficienza MC > 0.

        if (N_DATA <= 0. || N_MC <= 0.)
            continue;

        if (epsMC <= 0.)
            continue;

        const Double_t ratio = epsDATA / epsMC;


        // --------------------------------------------------------
        // ASYMMETRIC ERROR PROPAGATION
        // --------------------------------------------------------

        const Double_t ratioErrLow = std::sqrt(std::pow(errDATA_low / epsMC, 2.) + std::pow(epsDATA * errMC_up / (epsMC * epsMC), 2.));
        const Double_t ratioErrUp = std::sqrt(std::pow(errDATA_up / epsMC, 2.) + std::pow(epsDATA * errMC_low / (epsMC * epsMC), 2.));

        // --------------------------------------------------------
        // FILL RATIO GRAPH
        // --------------------------------------------------------

        const Int_t pointRatio = gRatio->GetN();
        gRatio->SetPoint(pointRatio, E_center, ratio);
        gRatio->SetPointError(pointRatio, E_errLow, E_errUp, ratioErrLow, ratioErrUp);

        // --------------------------------------------------------
        // AUTOMATIC RATIO RANGE
        // --------------------------------------------------------

        minRatio = std::min(minRatio, ratio - ratioErrLow);
        maxRatio = std::max(maxRatio, ratio + ratioErrUp);

        // --------------------------------------------------------
        // PRINT RESULTS
        // --------------------------------------------------------

        cout << "Bin " << ibin << "  E = [" << E_low << ", " << E_up << "] GeV" << endl;
        cout << "  DATA = " << epsDATA << " -" << errDATA_low << " +" << errDATA_up << endl;
        cout << "  MC   = " << epsMC << " -" << errMC_low << " +" << errMC_up << endl;
        cout << "  DATA/MC = " << ratio << " -" << ratioErrLow << " +" << ratioErrUp << endl;
    }


    // ============================================================
    // GRAPH STYLE
    // ============================================================

    // DATA
    gDATA->SetLineColor(kRed + 1);
    gDATA->SetMarkerColor(kRed + 1);

    gDATA->SetMarkerStyle(20);
    gDATA->SetMarkerSize(0.9);

    gDATA->SetLineWidth(2);


    // MC
    gMC->SetLineColor(kBlue + 1);
    gMC->SetMarkerColor(kBlue + 1);

    gMC->SetMarkerStyle(24);
    gMC->SetMarkerSize(0.9);

    gMC->SetLineWidth(2);


    // RATIO
    gRatio->SetLineColor(kBlack);
    gRatio->SetMarkerColor(kBlack);

    gRatio->SetMarkerStyle(20);
    gRatio->SetMarkerSize(0.8);

    gRatio->SetLineWidth(2);


    // ============================================================
    // CANVAS
    // ============================================================

    TCanvas *c = new TCanvas("c_HET_DATA_MC","HET efficiency DATA-MC comparison",900,850);

    // ============================================================
    // UPPER PAD
    // ============================================================

    TPad *pad1 = new TPad("pad1", "HET efficiency", 0., 0.30, 1., 1.);
    pad1->SetBottomMargin(0.03);
    pad1->SetTopMargin(0.07);

    pad1->SetLeftMargin(0.13);
    pad1->SetRightMargin(0.04);

    pad1->SetLogx();

    pad1->SetGridx();
    pad1->SetGridy();

    pad1->SetTickx();
    pad1->SetTicky();

    pad1->Draw();
    pad1->cd();

    // ============================================================
    // UPPER FRAME
    // ============================================================

    TH1F *frame1 = pad1->DrawFrame(Emin, 0., Emax, 1.12);
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle("HET efficiency");

    frame1->GetYaxis()->SetTitleSize(0.055);
    frame1->GetYaxis()->SetLabelSize(0.045);

    frame1->GetYaxis()->SetTitleOffset(1.10);

    frame1->GetYaxis()->SetNdivisions(505);

    // Nascondo i valori dell'asse X nel pannello superiore

    frame1->GetXaxis()->SetLabelSize(0.);
    frame1->GetXaxis()->SetTitleSize(0.);

    // ============================================================
    // DRAW DATA AND MC
    // ============================================================

    gMC->Draw("P SAME");
    gDATA->Draw("P SAME");

    // ============================================================
    // LEGEND
    // ============================================================

    TLegend *legend = new TLegend(0.54,0.18,0.88,0.33);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.040);
    legend->AddEntry(gDATA,"Flight data (p+He)","lep");
    legend->AddEntry(gMC,"MC (p+He)","lep");
    legend->Draw();

    // ============================================================
    // LOWER PAD
    // ============================================================

    c->cd();
    TPad *pad2 = new TPad("pad2","DATA / MC", 0., 0., 1., 0.30);

    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.32);

    pad2->SetLeftMargin(0.13);
    pad2->SetRightMargin(0.04);

    pad2->SetLogx();

    pad2->SetGridx();
    pad2->SetGridy();

    pad2->SetTickx();
    pad2->SetTicky();

    pad2->Draw();
    pad2->cd();


    // ============================================================
    // AUTOMATIC Y RANGE FOR RATIO
    // ============================================================

    const Double_t deltaRatio = std::max(0.1, maxRatio - minRatio);
    const Double_t ratioMin = std::max(0., minRatio - 0.1 * deltaRatio);
    const Double_t ratioMax = maxRatio + 0.1 * deltaRatio;

    // ============================================================
    // LOWER FRAME
    // ============================================================

    TH1F *frame2 = pad2->DrawFrame(Emin, ratioMin, Emax, ratioMax);
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle("DATA / MC");
    frame2->GetXaxis()->SetTitle("E_{BGO} [GeV]");

    // Y axis
    frame2->GetYaxis()->SetTitleSize(0.115);
    frame2->GetYaxis()->SetLabelSize(0.095);
    frame2->GetYaxis()->SetTitleOffset(0.47);
    frame2->GetYaxis()->SetNdivisions(505);

    // X axis
    frame2->GetXaxis()->SetTitleSize(0.12);
    frame2->GetXaxis()->SetLabelSize(0.095);
    frame2->GetXaxis()->SetTitleOffset(1.05);
    frame2->GetXaxis()->SetMoreLogLabels();
    frame2->GetXaxis()->SetNoExponent(kFALSE);

    // ============================================================
    // REFERENCE LINE: DATA / MC = 1
    // ============================================================

    TLine *line1 = new TLine(Emin, 1., Emax, 1.);
    line1->SetLineColor(kGray + 2);
    line1->SetLineStyle(2);
    line1->SetLineWidth(2);
    line1->Draw("SAME");

    // ============================================================
    // DRAW RATIO
    // ============================================================

    gRatio->Draw("P SAME");

    // ============================================================
    // UPDATE CANVAS
    // ============================================================

    c->cd();
    c->Update();

    // ============================================================
    // SAVE PLOT
    // ============================================================

    c->SaveAs("HET_Efficiency_DATA_MC_pHe.png");
    c->SaveAs("HET_Efficiency_DATA_MC_pHe.pdf");

    // ============================================================
    // SAVE ROOT OBJECTS
    // ============================================================

    TFile *fout = TFile::Open("HET_Efficiency_DATA_MC_Comparison.root", "RECREATE");

    if (fout && !fout->IsZombie()) {

        fout->cd();

        gDATA->Write();
        gMC->Write();
        gRatio->Write();

        c->Write();

        fout->Close();

    }


    // ============================================================
    // END
    // ============================================================

    fDATA->Close();
    fMC->Close();


    cout << endl;

    cout << "========================================" << endl;
    cout << "HET DATA-MC comparison completed." << endl;
    cout << "========================================" << endl;

}
