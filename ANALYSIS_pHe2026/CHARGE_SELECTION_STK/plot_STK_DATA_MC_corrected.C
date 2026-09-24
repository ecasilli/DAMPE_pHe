
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"

#include <iostream>

using namespace std;

void DrawDataMCShape(TH1F *hData, TH1F *hMCP, TH1F *hMCHe, const TString &title, const TString &tag){
    if (!hData || !hMCP || !hMCHe) {
        cout << "ERROR: missing histogram for " << tag << endl;
        return;
    }

    // Clone histograms
    TH1F *hD = (TH1F*)hData->Clone(Form("%s_data", tag.Data()));
    TH1F *hP = (TH1F*)hMCP->Clone(Form("%s_pmc", tag.Data()));
    TH1F *hHe = (TH1F*)hMCHe->Clone(Form("%s_hemc", tag.Data()));

    hD->SetDirectory(0);
    hP->SetDirectory(0);
    hHe->SetDirectory(0);

    // Normalization windows
    const Double_t pMin  = 30.;
    const Double_t pMax  = 90.;

    const Double_t heMin = 180.;
    const Double_t heMax = 320.;

    // PROTON normalization

    Int_t binPMin = hD->GetXaxis()->FindBin(pMin);
    Int_t binPMax = hD->GetXaxis()->FindBin(pMax);

    Double_t intDataP = hD->Integral(binPMin, binPMax);
    Double_t intMCP = hP->Integral(binPMin, binPMax);

    if (intMCP > 0.) {

        Double_t scaleP = intDataP / intMCP;
        hP->Scale(scaleP);

        cout << tag << " proton scale = " << scaleP << endl;
    }

    // HELIUM normalization

    Int_t binHeMin = hD->GetXaxis()->FindBin(heMin);
    Int_t binHeMax = hD->GetXaxis()->FindBin(heMax);

    Double_t intDataHe = hD->Integral(binHeMin, binHeMax);
    Double_t intMCHe = hHe->Integral(binHeMin, binHeMax);

    if (intMCHe > 0.) {

        Double_t scaleHe = intDataHe / intMCHe;
        hHe->Scale(scaleHe);

        cout << tag << " helium scale = " << scaleHe << endl;
    }

    // STYLE
    hD->SetLineColor(kBlack);
    hD->SetMarkerColor(kBlack);
    hD->SetMarkerStyle(20);
    hD->SetMarkerSize(0.45);

    hP->SetLineColor(kBlue + 1);
    hP->SetLineWidth(2);

    hHe->SetLineColor(kRed + 1);
    hHe->SetLineWidth(2);

    // X range
    const Double_t xmin = 0.;
    const Double_t xmax = 660.;

    hD->GetXaxis()->SetRangeUser(xmin, xmax);
    hP->GetXaxis()->SetRangeUser(xmin, xmax);
    hHe->GetXaxis()->SetRangeUser(xmin,xmax);

    // Y range
    Double_t ymax = TMath::Max(hD->GetMaximum(), TMath::Max(hP->GetMaximum(), hHe->GetMaximum() ));
    hD->SetMinimum(0.);
    hD->SetMaximum(1.20 * ymax);

    // AXES / TITLE
    hD->SetTitle(title);
    hD->GetXaxis()->SetTitle("STK charge signal (ADC)");
    hD->GetYaxis()->SetTitle("Counts");
    hD->GetXaxis()->SetTitleSize(0.045);
    hD->GetYaxis()->SetTitleSize(0.045);
    hD->GetXaxis()->SetLabelSize(0.040);
    hD->GetYaxis()->SetLabelSize(0.040);

    // DRAW
    hD->Draw("E");
    hP->Draw("HIST SAME");
    hHe->Draw("HIST SAME");

    // DATA again on top
    hD->Draw("E SAME");

    // Vertical STK cuts
    Double_t yMin = 0.;
    Double_t yMax = hD->GetMaximum() ;

    TLine *lineLow = new TLine(25., yMin, 25., yMax);
    lineLow->SetLineColor(kGreen + 2);
    lineLow->SetLineWidth(2);
    lineLow->SetLineStyle(2);
    lineLow->Draw("SAME");

    TLine *lineHigh = new TLine(450., yMin, 450., yMax);
    lineHigh->SetLineColor(kGreen + 2);
    lineHigh->SetLineWidth(2);
    lineHigh->SetLineStyle(2);
    lineHigh->Draw("SAME");

    // LEGEND
    TLegend *leg = new TLegend(0.35, 0.65, 0.72, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry(hD, "DATA", "lep");
    leg->AddEntry(hP,"Proton MC corrected", "l");
    leg->AddEntry(hHe,"Helium MC corrected","l");
    leg->Draw();

    gPad->SetTicks();
    gPad->Modified();
    gPad->Update();
}



// ============================================================================
// MAIN
// ============================================================================
void plot_STK_DATA_MC_corrected()
{
    gStyle->SetOptStat(0);

    // INPUT FILES
    TString dataFileName = "../ROOT_FILES/PHe_STKcharge_adc_Orb120Month_480bins_15sett26_nocut05.root";
    TString mcFileName = "../ROOT_FILES/PHe_MC_p_He_5PeV_STKcharge_adc_480bins_16sett26_nocut05_CORRECTED.root";


    TFile *fData = TFile::Open(dataFileName, "READ");
    TFile *fMC = TFile::Open(mcFileName, "READ");

    if (!fData || fData->IsZombie()) {
        cout << "ERROR opening DATA file" << endl;
        return;
    }

    if (!fMC || fMC->IsZombie()) {
        cout << "ERROR opening MC file" << endl;
        return;
    }

    // ENERGY INTERVALS
    const Int_t nE = 4;

    Double_t eMin[nE] = {20., 100., 1000., 10000. };
    Double_t eMax[nE] = {100., 1000., 10000., 100000.};

    // HISTOGRAMS
    TH1F *hData[nE];
    TH1F *hMCP[nE];
    TH1F *hMCHe[nE];

    for (Int_t i = 0; i < nE; ++i) {

        // DATA:
        // h01, h02, h03, h04
        hData[i] = (TH1F*)fData->Get(Form("h%02d", i+1));

        // corrected proton MC
        hMCP[i] = (TH1F*)fMC->Get(Form("h%02d_p", i+1));

        // corrected helium MC
        hMCHe[i] = (TH1F*)fMC->Get(Form("h%02d_he", i+1));

        if (!hData[i])
            cout << "Missing DATA h" << Form("%02d", i+1) << endl;

        if (!hMCP[i])
            cout << "Missing proton MC h" << Form("%02d_p", i+1) << endl;

        if (!hMCHe[i])
            cout << "Missing helium MC h" << Form("%02d_he", i+1) << endl;
    }

    // CANVAS
    TCanvas *c = new TCanvas( "c_STK_shape_comparison", "STK DATA vs corrected MC", 1400, 1000);
    c->Divide(2, 2, 0.01, 0.01);

    for (Int_t i = 0; i<nE; ++i) {

        c->cd(i+1);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.05);
        gPad->SetBottomMargin(0.12);
        gPad->SetTopMargin(0.10);

        TString title = Form("%.0f < E_{BGO} < %.0f GeV", eMin[i], eMax[i]);
        DrawDataMCShape(hData[i], hMCP[i], hMCHe[i], title, Form("energy_%d", i));
    }

    // SAVE
    c->SaveAs("PLOTS/STK_DATA_MC_corrected_shape_comparison_3.pdf");

    TFile *fout = new TFile( "STK_DATA_MC_corrected_shape_comparison_3.root", "RECREATE");

    if (fout && !fout->IsZombie()) {
        fout->cd();
        c->Write();
        fout->Close();
    }

    fData->Close();
    fMC->Close();

    cout << endl;
    cout << "==========================================" << endl;
    cout << " Comparison completed." << endl;
    cout << "==========================================" << endl;
}