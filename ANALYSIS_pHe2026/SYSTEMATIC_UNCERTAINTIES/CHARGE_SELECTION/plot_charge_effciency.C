// plot_charge_efficiency.C
//
// Diagnostic plots for the vertex-dependent p+He charge selection.
//
// DATA histograms expected:
//   h1PrePSD_orb, h1SelPSD_orb, h1PreSTK_orb, h1SelSTK_orb
// MC histograms expected:
//   h2PreCharge_PSD, h2Ntrig_wgt_PSD, h2PreCharge_STK, h2Ntrig_wgt_STK
//
// Usage:
//   root -l 'plot_charge_efficiency.C()'
// or
//   root -l 'plot_charge_efficiency.C("myData.root","myMC.root")'

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

namespace {

TH1D* MakeSubsetEfficiency(const TH1D *hPass, const TH1D *hPre, const char *name) {
    if (!hPass || !hPre) return nullptr;

    if (hPass->GetNbinsX() != hPre->GetNbinsX()) {
        cerr << "ERROR: incompatible binning for " << name << endl;
        return nullptr;
    }

    TH1D *hEff = (TH1D*)hPre->Clone(name);
    hEff->Reset("ICES");
    hEff->SetDirectory(nullptr);
    hEff->SetTitle("");

    for (int ib = 1; ib <= hPre->GetNbinsX(); ++ib) {
        const double pre  = hPre ->GetBinContent(ib);
        const double pass = hPass->GetBinContent(ib);

        if (pre <= 0.) {
            hEff->SetBinContent(ib, 0.);
            hEff->SetBinError(ib, 0.);
            continue;
        }

        double eps = pass / pre;

        if (eps < -1e-10 || eps > 1. + 1e-10) {
            cerr << "WARNING: efficiency outside [0,1] in bin " << ib
                 << " for " << name
                 << ": pass=" << pass << " pre=" << pre
                 << " eps=" << eps << endl;
        }

        eps = std::max(0.0, std::min(1.0, eps));

        // Because pass is a subset of pre, use the subset-aware variance.
        // This reduces to eps*(1-eps)/N for unweighted data.
        const double sumw2Pre  = std::pow(hPre ->GetBinError(ib), 2);
        const double sumw2Pass = std::pow(hPass->GetBinError(ib), 2);
        const double sumw2Fail = std::max(0.0, sumw2Pre - sumw2Pass);

        const double variance = ( sumw2Pass * std::pow(1. - eps, 2) + sumw2Fail * std::pow(eps, 2) ) / (pre * pre);

        hEff->SetBinContent(ib, eps);
        hEff->SetBinError(ib, std::sqrt(std::max(0.0, variance)));
    }

    return hEff;
}

TH1D* MakeRatio(const TH1D *hData, const TH1D *hMC, const char *name) {
    if (!hData || !hMC) return nullptr;

    TH1D *hRatio = (TH1D*)hData->Clone(name);
    hRatio->Reset("ICES");
    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle("");

    for (int ib = 1; ib <= hData->GetNbinsX(); ++ib) {
        const double d  = hData->GetBinContent(ib);
        const double ed = hData->GetBinError(ib);
        const double m  = hMC->GetBinContent(ib);
        const double em = hMC->GetBinError(ib);

        if (d <= 0. || m <= 0.) continue;

        const double r = d / m;
        const double er = r * std::sqrt(std::pow(ed/d, 2) + std::pow(em/m, 2));

        hRatio->SetBinContent(ib, r);
        hRatio->SetBinError(ib, er);
    }

    return hRatio;
}

void GetPlotRange(const TH1D *h, bool useFluxRange, int removeFirst, int removeLast, double &xmin, double &xmax) {
    const int nb = h->GetNbinsX();
    int first = 1;
    int last  = nb;

    if (useFluxRange) {
        first = std::min(nb, std::max(1, removeFirst + 1));
        last  = std::max(first, std::min(nb, nb - removeLast));
    }

    xmin = h->GetXaxis()->GetBinLowEdge(first);
    xmax = h->GetXaxis()->GetBinUpEdge(last);
}

double RatioHalfRange(const TH1D *hRatio, double xmin, double xmax) {
    if (!hRatio) return 0.05;

    double maxDev = 0.;
    for (int ib=1; ib<=hRatio->GetNbinsX(); ++ib) {
        const double x = hRatio->GetXaxis()->GetBinCenter(ib);
        if (x < xmin || x > xmax) continue;

        const double y = hRatio->GetBinContent(ib);
        const double e = hRatio->GetBinError(ib);
        if (y <= 0.) continue;

        maxDev = std::max(maxDev, std::fabs(y - 1.) + e);
    }

    return std::max(0.02, 1.25*maxDev);
}

void DrawDetectorCanvas(const char *detector, TH1D *hDataEff, TH1D *hMCEff, TH1D *hRatio, double xmin, double xmax, const char *outBase) {
    const bool haveData = (hDataEff != nullptr && hRatio != nullptr);

    TCanvas *c = new TCanvas(Form("cEff_%s",detector),Form("%s charge efficiency",detector),900,820);

    TPad *p1 = nullptr;
    TPad *p2 = nullptr;

    if (haveData) {
        p1 = new TPad(Form("p1_%s",detector),"",0.,0.31,1.,1.);
        p2 = new TPad(Form("p2_%s",detector),"",0.,0.,1.,0.31);
        p1->SetBottomMargin(0.02);
        p2->SetTopMargin(0.03);
        p2->SetBottomMargin(0.32);
        p1->Draw();
        p2->Draw();
    } else {
        p1 = new TPad(Form("p1_%s",detector),"",0.,0.,1.,1.);
        p1->SetBottomMargin(0.13);
        p1->Draw();
    }

    p1->cd();
    p1->SetLogx();
    p1->SetGridx();
    p1->SetGridy();

    TH1D *frame = (TH1D*)hMCEff->Clone(Form("frame_%s",detector));
    frame->Reset("ICES");
    frame->SetDirectory(nullptr);
    frame->GetXaxis()->SetRangeUser(xmin,xmax);
    frame->GetYaxis()->SetRangeUser(0.,1.05);
    frame->GetYaxis()->SetTitle(Form("%s charge efficiency",detector));
    frame->GetXaxis()->SetTitle("BGO reconstructed energy [GeV]");

    if (haveData) {
        frame->GetXaxis()->SetLabelSize(0.);
        frame->GetXaxis()->SetTitleSize(0.);
        frame->GetYaxis()->SetTitleSize(0.060);
        frame->GetYaxis()->SetLabelSize(0.050);
        frame->GetYaxis()->SetTitleOffset(0.90);
    }

    frame->Draw("AXIS");

    hMCEff->SetMarkerStyle(24);
    hMCEff->SetMarkerSize(0.85);
    hMCEff->SetMarkerColor(kRed+1);
    hMCEff->SetLineColor(kRed+1);
    hMCEff->SetLineWidth(2);
    hMCEff->Draw("E1 SAME");

    if (hDataEff) {
        hDataEff->SetMarkerStyle(20);
        hDataEff->SetMarkerSize(0.85);
        hDataEff->SetMarkerColor(kBlack);
        hDataEff->SetLineColor(kBlack);
        hDataEff->SetLineWidth(2);
        hDataEff->Draw("E1 SAME");
    }

    TLegend *leg = new TLegend(0.64,0.18,0.88,0.33);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    if (hDataEff) leg->AddEntry(hDataEff,"Flight data","lep");
    leg->AddEntry(hMCEff,"MC p+He","lep");
    leg->Draw();

    if (haveData) {
        p2->cd();
        p2->SetLogx();
        p2->SetGridx();
        p2->SetGridy();

        const double halfRange = RatioHalfRange(hRatio,xmin,xmax);
        TH1D *rframe = (TH1D*)hRatio->Clone(Form("rframe_%s",detector));
        rframe->Reset("ICES");
        rframe->SetDirectory(nullptr);
        rframe->GetXaxis()->SetRangeUser(xmin,xmax);
        rframe->GetYaxis()->SetRangeUser(1.-halfRange,1.+halfRange);
        rframe->GetXaxis()->SetTitle("BGO reconstructed energy [GeV]");
        rframe->GetYaxis()->SetTitle("Data / MC");
        rframe->GetXaxis()->SetTitleSize(0.115);
        rframe->GetXaxis()->SetLabelSize(0.095);
        rframe->GetYaxis()->SetTitleSize(0.100);
        rframe->GetYaxis()->SetLabelSize(0.082);
        rframe->GetYaxis()->SetTitleOffset(0.52);
        rframe->GetYaxis()->SetNdivisions(505);
        rframe->Draw("AXIS");

        hRatio->SetMarkerStyle(20);
        hRatio->SetMarkerSize(0.80);
        hRatio->SetMarkerColor(kBlack);
        hRatio->SetLineColor(kBlack);
        hRatio->Draw("E1 SAME");

        TLine *l1 = new TLine(xmin,1.,xmax,1.);
        l1->SetLineStyle(2);
        l1->SetLineWidth(2);
        l1->Draw();
    }

    c->SaveAs(Form("%s.pdf",outBase));
    c->SaveAs(Form("%s.png",outBase));
}

void WriteTable(const char *filename, TH1D *dPSD, TH1D *mPSD, TH1D *rPSD, TH1D *dSTK, TH1D *mSTK, TH1D *rSTK) {
    ofstream out(filename);
    out << "# E_GeV"
        << " dataPSD errDataPSD mcPSD errMcPSD ratioPSD errRatioPSD"
        << " dataSTK errDataSTK mcSTK errMcSTK ratioSTK errRatioSTK\n";

    const int nb = mPSD->GetNbinsX();
    out << setprecision(12);

    for (int ib=1; ib<=nb; ++ib) {
        const double E = mPSD->GetXaxis()->GetBinCenter(ib);
        auto c = [ib](TH1D *h) { return h ? h->GetBinContent(ib) : -1.; };
        auto e = [ib](TH1D *h) { return h ? h->GetBinError(ib)   : -1.; };

        out << E << " "
            << c(dPSD) << " " << e(dPSD) << " "
            << c(mPSD) << " " << e(mPSD) << " "
            << c(rPSD) << " " << e(rPSD) << " "
            << c(dSTK) << " " << e(dSTK) << " "
            << c(mSTK) << " " << e(mSTK) << " "
            << c(rSTK) << " " << e(rSTK) << "\n";
    }
}

} // namespace

void plot_charge_efficiency(
    const char *dataFile =
      "ROOT_FILES/PHe_skim_Orb120Month_5binperdecade_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_nominal.root",
    const char *mcFile =
      "ROOT_FILES/PHe_MC_p_He_5PeV_5binperdecade_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_nominal.root")
{
    gStyle->SetOptStat(0);

    // Same range used in your flux comparisons.
    // Set false to inspect all 30 bins.
    const bool useFluxRange = true;
    const int removeFirst = 3;
    const int removeLast  = 5;

    gSystem->mkdir("PLOTS",kTRUE);
    gSystem->mkdir("TXT_FILES",kTRUE);

    // =========================
    // MC -- mandatory
    // =========================
    TFile *fMC = TFile::Open(mcFile,"READ");
    if (!fMC || fMC->IsZombie()) {
        cerr << "ERROR: cannot open MC file: " << mcFile << endl;
        return;
    }

    TH2D *h2PrePSD = dynamic_cast<TH2D*>(fMC->Get("h2PreCharge_PSD"));
    TH2D *h2SelPSD = dynamic_cast<TH2D*>(fMC->Get("h2Ntrig_wgt_PSD"));
    TH2D *h2PreSTK = dynamic_cast<TH2D*>(fMC->Get("h2PreCharge_STK"));
    TH2D *h2SelSTK = dynamic_cast<TH2D*>(fMC->Get("h2Ntrig_wgt_STK"));

    if (!h2PrePSD || !h2SelPSD || !h2PreSTK || !h2SelSTK) {
        cerr << "ERROR: missing MC histograms. Expected:" << endl;
        cerr << "  h2PreCharge_PSD" << endl;
        cerr << "  h2Ntrig_wgt_PSD" << endl;
        cerr << "  h2PreCharge_STK" << endl;
        cerr << "  h2Ntrig_wgt_STK" << endl;
        fMC->Close();
        return;
    }

    // Projection Y = reconstructed BGO energy, for direct Data/MC comparison.
    TH1D *mcPrePSD_reco = h2PrePSD->ProjectionY("mcPrePSD_reco",1,-1,"e");
    TH1D *mcSelPSD_reco = h2SelPSD->ProjectionY("mcSelPSD_reco",1,-1,"e");
    TH1D *mcPreSTK_reco = h2PreSTK->ProjectionY("mcPreSTK_reco",1,-1,"e");
    TH1D *mcSelSTK_reco = h2SelSTK->ProjectionY("mcSelSTK_reco",1,-1,"e");

    mcPrePSD_reco->SetDirectory(nullptr);
    mcSelPSD_reco->SetDirectory(nullptr);
    mcPreSTK_reco->SetDirectory(nullptr);
    mcSelSTK_reco->SetDirectory(nullptr);

    TH1D *mcEffPSD_reco = MakeSubsetEfficiency(mcSelPSD_reco,mcPrePSD_reco,"mcEffPSD_reco");
    TH1D *mcEffSTK_reco = MakeSubsetEfficiency(mcSelSTK_reco,mcPreSTK_reco,"mcEffSTK_reco");

    // Projection X = true energy, useful MC-only diagnostic.
    TH1D *mcPrePSD_true = h2PrePSD->ProjectionX("mcPrePSD_true",1,-1,"e");
    TH1D *mcSelPSD_true = h2SelPSD->ProjectionX("mcSelPSD_true",1,-1,"e");
    TH1D *mcPreSTK_true = h2PreSTK->ProjectionX("mcPreSTK_true",1,-1,"e");
    TH1D *mcSelSTK_true = h2SelSTK->ProjectionX("mcSelSTK_true",1,-1,"e");

    mcPrePSD_true->SetDirectory(nullptr);
    mcSelPSD_true->SetDirectory(nullptr);
    mcPreSTK_true->SetDirectory(nullptr);
    mcSelSTK_true->SetDirectory(nullptr);

    TH1D *mcEffPSD_true = MakeSubsetEfficiency(mcSelPSD_true,mcPrePSD_true,"mcEffPSD_true");
    TH1D *mcEffSTK_true = MakeSubsetEfficiency(mcSelSTK_true,mcPreSTK_true,"mcEffSTK_true");

    fMC->Close();

    // =========================
    // DATA -- optional
    // =========================
    TH1D *dataEffPSD = nullptr;
    TH1D *dataEffSTK = nullptr;

    TFile *fData = TFile::Open(dataFile,"READ");

    if (!fData || fData->IsZombie()) {
        cout << "INFO: Data file not available yet: " << dataFile << endl;
        cout << "Proceeding with MC-only plots." << endl;
        if (fData) { fData->Close(); delete fData; fData=nullptr; }
    } else {
        TH1D *hPrePSD = dynamic_cast<TH1D*>(fData->Get("h1PrePSD_orb"));
        TH1D *hSelPSD = dynamic_cast<TH1D*>(fData->Get("h1SelPSD_orb"));
        TH1D *hPreSTK = dynamic_cast<TH1D*>(fData->Get("h1PreSTK_orb"));
        TH1D *hSelSTK = dynamic_cast<TH1D*>(fData->Get("h1SelSTK_orb"));

        if (!hPrePSD || !hSelPSD || !hPreSTK || !hSelSTK) {
            cout << "INFO: Data file exists but the new pre-charge histograms are missing." << endl;
            cout << "Proceeding with MC-only plots." << endl;
        } else {
            TH1D *dPrePSD = (TH1D*)hPrePSD->Clone("dPrePSD");
            TH1D *dSelPSD = (TH1D*)hSelPSD->Clone("dSelPSD");
            TH1D *dPreSTK = (TH1D*)hPreSTK->Clone("dPreSTK");
            TH1D *dSelSTK = (TH1D*)hSelSTK->Clone("dSelSTK");

            dPrePSD->SetDirectory(nullptr);
            dSelPSD->SetDirectory(nullptr);
            dPreSTK->SetDirectory(nullptr);
            dSelSTK->SetDirectory(nullptr);

            dataEffPSD = MakeSubsetEfficiency(dSelPSD,dPrePSD,"dataEffPSD");
            dataEffSTK = MakeSubsetEfficiency(dSelSTK,dPreSTK,"dataEffSTK");
        }

        fData->Close();
    }

    // =========================
    // Data / MC ratios
    // =========================
    TH1D *ratioPSD = nullptr;
    TH1D *ratioSTK = nullptr;

    if (dataEffPSD && dataEffSTK) {
        ratioPSD = MakeRatio(dataEffPSD,mcEffPSD_reco,"ratioPSD");
        ratioSTK = MakeRatio(dataEffSTK,mcEffSTK_reco,"ratioSTK");
    }

    // =========================
    // Plot range
    // =========================
    double xmin=0., xmax=0.;
    GetPlotRange(mcEffPSD_reco,useFluxRange,removeFirst,removeLast,xmin,xmax);

    // =========================
    // PSD / STK reconstructed-energy plots
    // =========================
    DrawDetectorCanvas("PSD",dataEffPSD,mcEffPSD_reco,ratioPSD, xmin,xmax,"PLOTS/pHe_charge_efficiency_PSD");
    DrawDetectorCanvas("STK",dataEffSTK,mcEffSTK_reco,ratioSTK, xmin,xmax,"PLOTS/pHe_charge_efficiency_STK");

    // =========================
    // MC-only true-energy plot
    // =========================
    TCanvas *cTrue = new TCanvas("cChargeEffMCTrue","MC charge efficiency vs true energy",900,700);
    cTrue->SetLogx();
    cTrue->SetGridx();
    cTrue->SetGridy();

    mcEffPSD_true->SetMarkerStyle(20);
    mcEffPSD_true->SetMarkerColor(kBlue+1);
    mcEffPSD_true->SetLineColor(kBlue+1);
    mcEffPSD_true->SetLineWidth(2);

    mcEffSTK_true->SetMarkerStyle(21);
    mcEffSTK_true->SetMarkerColor(kRed+1);
    mcEffSTK_true->SetLineColor(kRed+1);
    mcEffSTK_true->SetLineWidth(2);

    TH1D *trueFrame = (TH1D*)mcEffPSD_true->Clone("trueFrame");
    trueFrame->Reset("ICES");
    trueFrame->SetDirectory(nullptr);
    trueFrame->GetXaxis()->SetRangeUser(xmin,xmax);
    trueFrame->GetYaxis()->SetRangeUser(0.,1.05);
    trueFrame->GetXaxis()->SetTitle("MC true energy [GeV]");
    trueFrame->GetYaxis()->SetTitle("Charge efficiency");
    trueFrame->Draw("AXIS");

    mcEffPSD_true->Draw("E1 SAME");
    mcEffSTK_true->Draw("E1 SAME");

    TLegend *legTrue = new TLegend(0.65,0.18,0.88,0.31);
    legTrue->SetBorderSize(0);
    legTrue->SetFillStyle(0);
    legTrue->AddEntry(mcEffPSD_true,"PSD branch","lep");
    legTrue->AddEntry(mcEffSTK_true,"STK branch","lep");
    legTrue->Draw();

    cTrue->SaveAs("PLOTS/pHe_charge_efficiency_MC_true.pdf");
    cTrue->SaveAs("PLOTS/pHe_charge_efficiency_MC_true.png");

    // =========================
    // Text output
    // =========================
    WriteTable("TXT_FILES/pHe_charge_efficiency_reco.dat",dataEffPSD,mcEffPSD_reco,ratioPSD,dataEffSTK,mcEffSTK_reco,ratioSTK);

    cout << endl;
    cout << "==============================================" << endl;
    cout << "Charge-efficiency diagnostic completed" << endl;
    cout << "==============================================" << endl;
    cout << "MC file:   " << mcFile << endl;
    cout << "Data file: " << dataFile << endl;

    if (dataEffPSD && dataEffSTK)
        cout << "Data/MC comparison: AVAILABLE" << endl;
    else
        cout << "Data/MC comparison: not available yet (MC-only)." << endl;

    cout << "Outputs:" << endl;
    cout << "  PLOTS/pHe_charge_efficiency_PSD.pdf/.png" << endl;
    cout << "  PLOTS/pHe_charge_efficiency_STK.pdf/.png" << endl;
    cout << "  PLOTS/pHe_charge_efficiency_MC_true.pdf/.png" << endl;
    cout << "  TXT_FILES/pHe_charge_efficiency_reco.dat" << endl;
}

