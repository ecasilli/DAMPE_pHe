#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"

// Parte dai TH1D ORIGINALI (48 bin, 8/decade, 10 GeV--10 PeV).
// NON rielabora eventi e NON ribinna le efficienze gia' divise.
// Binning: 8 bin/decade fino a 1 TeV; 4 tra 1 e 10 TeV;
//          2 tra 10 e 100 TeV; 1 tra 100 TeV e 10 PeV.
void HET_pHe_rebin_and_plot_DATA_MC()
{
    const char *dataName = "HETeff_skim_Orb108Month_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7.root";
    const char *mcName   = "HET_eff_MC_pHe_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC.root";

    TFile *fData = TFile::Open(dataName, "READ");
    TFile *fMC   = TFile::Open(mcName, "READ");
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        std::cerr << "ERRORE: file DATA/MC non apribili." << std::endl;
        return;
    }

    TH1D *hDU=nullptr, *hDA=nullptr, *hMU=nullptr, *hMA=nullptr;
    fData->GetObject("h1_Unb", hDU);
    fData->GetObject("h1Nobs_And", hDA);
    fMC->GetObject("h1_Unb", hMU);
    fMC->GetObject("h1Nobs_And", hMA);
    if (!hDU || !hDA || !hMU || !hMA) {
        std::cerr << "ERRORE: manca h1_Unb o h1Nobs_And." << std::endl;
        return;
    }

    // Indici dei BORDI del vecchio istogramma: da 0 a 48.
    // 0..16 = 10 GeV..1 TeV, un bin originale per volta;
    // 18,20,22,24 = 4 bin tra 1 e 10 TeV;
    // 28,32 = 2 bin tra 10 e 100 TeV;
    // 48 = 1 bin tra 100 TeV e 10 PeV.
    std::vector<int> edgeIndex;
    for (int i=0; i<=16; ++i) edgeIndex.push_back(i);
    for (int i=18; i<=24; i+=2) edgeIndex.push_back(i);
    edgeIndex.push_back(28);
    edgeIndex.push_back(32);
    edgeIndex.push_back(48);

    const TH1D *orig[4] = {hDU,hDA,hMU,hMA};
    for (auto *h : orig) {
        if (h->GetNbinsX() != 48) {
            std::cerr << "ERRORE: richiesto l'istogramma ORIGINALE da 48 bin: "
                      << h->GetName() << std::endl;
            return;
        }
    }
    // Costruisce i bordi esattamente da quelli originali, senza approssimazioni.
    std::vector<std::vector<double>> allEdges(4);
    for (int k=0; k<4; ++k) {
        for (int idx : edgeIndex) {
            allEdges[k].push_back(orig[k]->GetXaxis()->GetBinLowEdge(idx+1));
        }
        for (std::size_t j=0; j<edgeIndex.size(); ++j) {
            const double x = allEdges[0][j], y = allEdges[k][j];
            if (std::abs(x-y)>1e-9*std::max(1.,std::abs(x))) {
                std::cerr << "ERRORE: bordi originari DATA/MC incompatibili." << std::endl;
                return;
            }
        }
    }
    const int nNew = static_cast<int>(edgeIndex.size())-1;

    // ROOT Rebin somma i contenuti e propaga correttamente Sumw2.
    TH1D *dTot = static_cast<TH1D*>(hDU->Rebin(nNew,"hDATA_Unb_rebin",allEdges[0].data()));
    TH1D *dPass= static_cast<TH1D*>(hDA->Rebin(nNew,"hDATA_And_rebin",allEdges[1].data()));
    TH1D *mTot = static_cast<TH1D*>(hMU->Rebin(nNew,"hMC_Unb_rebin",allEdges[2].data()));
    TH1D *mPass= static_cast<TH1D*>(hMA->Rebin(nNew,"hMC_And_rebin",allEdges[3].data()));
    if (!dTot || !dPass || !mTot || !mPass) {
        std::cerr << "ERRORE durante Rebin." << std::endl;
        return;
    }
    for (TH1D *h : {dTot,dPass,mTot,mPass}) h->SetDirectory(nullptr);

    if (!TEfficiency::CheckConsistency(*dPass,*dTot)) {
        std::cerr << "ERRORE: istogrammi DATA incompatibili per TEfficiency." << std::endl;
        return;
    }
    if (!TEfficiency::CheckConsistency(*mPass,*mTot,"w")) {
        std::cerr << "ERRORE: istogrammi MC pesati incompatibili per TEfficiency." << std::endl;
        return;
    }

    TEfficiency *effD = new TEfficiency(*dPass,*dTot);
    effD->SetName("HET_Eff_DATA_rebin");
    effD->SetStatisticOption(TEfficiency::kFCP);
    effD->SetConfidenceLevel(0.682689492);

    TEfficiency *effM = new TEfficiency(*mPass,*mTot);
    effM->SetName("HET_Eff_MC_rebin");
    effM->SetUseWeightedEvents();
    effM->SetStatisticOption(TEfficiency::kFNormal);
    effM->SetConfidenceLevel(0.682689492);

    TGraphAsymmErrors *gD = new TGraphAsymmErrors();
    TGraphAsymmErrors *gM = new TGraphAsymmErrors();
    TGraphAsymmErrors *gR = new TGraphAsymmErrors();
    gD->SetName("g_HET_DATA_rebin");
    gM->SetName("g_HET_MC_rebin");
    gR->SetName("g_DATA_over_MC_rebin");

    double ratioMin=1., ratioMax=1.;
    for (int ibin=1; ibin<=nNew; ++ibin) {
        const double elo=dTot->GetXaxis()->GetBinLowEdge(ibin);
        const double ehi=dTot->GetXaxis()->GetBinUpEdge(ibin);
        const double x=std::sqrt(elo*ehi); // centro geometrico in GeV
        const double exl=x-elo, exh=ehi-x;
        const bool okD = dTot->GetBinContent(ibin)>0.;
        const bool okM = mTot->GetBinContent(ibin)>0.;

        double d=0., dl=0., du=0., m=0., ml=0., mu=0.;
        if (okD) {
            d=effD->GetEfficiency(ibin);
            dl=effD->GetEfficiencyErrorLow(ibin);
            du=effD->GetEfficiencyErrorUp(ibin);
            if (std::isfinite(d) && std::isfinite(dl) && std::isfinite(du)) {
                int n=gD->GetN();
                gD->SetPoint(n,x,d);
                gD->SetPointError(n,exl,exh,dl,du);
            }
        }
        if (okM) {
            m=effM->GetEfficiency(ibin);
            ml=effM->GetEfficiencyErrorLow(ibin);
            mu=effM->GetEfficiencyErrorUp(ibin);
            if (std::isfinite(m) && std::isfinite(ml) && std::isfinite(mu)) {
                int n=gM->GetN();
                gM->SetPoint(n,x,m);
                gM->SetPointError(n,exl,exh,ml,mu);
            }
        }
        if (!okD || !okM || m<=0. || !std::isfinite(d) || !std::isfinite(m)
            || !std::isfinite(dl) || !std::isfinite(du)
            || !std::isfinite(ml) || !std::isfinite(mu)) continue;

        const double r=d/m;
        const double rl=std::sqrt(std::pow(dl/m,2.)+std::pow(d*mu/(m*m),2.));
        const double ru=std::sqrt(std::pow(du/m,2.)+std::pow(d*ml/(m*m),2.));
        if (!std::isfinite(r) || !std::isfinite(rl) || !std::isfinite(ru)) continue;
        const int n=gR->GetN();
        gR->SetPoint(n,x,r);
        gR->SetPointError(n,exl,exh,rl,ru);
        ratioMin=std::min(ratioMin,r-rl);
        ratioMax=std::max(ratioMax,r+ru);
    }

    // 
    // ============================================================
    // MANUAL Y-AXIS RANGES
    // ============================================================

    // Upper plot: HET efficiencies
    const double yEffMin = 0.05;
    const double yEffMax = 1.06;

    // Lower plot: DATA / MC ratio
    const double yRatioMin = 0.78;
    const double yRatioMax = 1.22;

    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gD->SetLineColor(kRed+1);  gD->SetMarkerColor(kRed+1);
    gD->SetMarkerStyle(20);    gD->SetMarkerSize(1.1); gD->SetLineWidth(2);
    gM->SetLineColor(kBlue+1); gM->SetMarkerColor(kBlue+1);
    gM->SetMarkerStyle(24);    gM->SetMarkerSize(1.1); gM->SetLineWidth(2);
    gR->SetLineColor(kBlack);  gR->SetMarkerColor(kBlack);
    gR->SetMarkerStyle(21);    gR->SetMarkerSize(1.0); gR->SetLineWidth(2);

    const double xmin=17., xmax=1.e5;
    TCanvas *c=new TCanvas("c_HET_pHe_rebin","HET efficiency p+He",1000,950);
    TPad *pad1=new TPad("pad1_pHe","efficiencies",0.,0.30,1.,1.);
    pad1->SetLeftMargin(0.10); pad1->SetRightMargin(0.03);
    pad1->SetBottomMargin(0.0); pad1->SetTopMargin(0.04);
    pad1->SetLogx(); pad1->SetGridx(); pad1->SetGridy();
    pad1->SetTickx(); pad1->SetTicky(); pad1->Draw(); pad1->cd();
    TH1F *frame1=pad1->DrawFrame(xmin,yEffMin,xmax,yEffMax);
    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0.);
    frame1->GetYaxis()->SetLabelSize(0.034);
    frame1->GetYaxis()->SetNdivisions(505);
    gM->Draw("P SAME"); gD->Draw("P SAME");
    TLatex y1; y1.SetTextAngle(90); y1.SetTextFont(42);
    y1.SetTextSize(0.037); y1.DrawLatexNDC(0.035,0.55,"HET efficiency");
    
    TLegend *leg=new TLegend(0.19,0.56,0.4,0.76);
    leg->SetLineColor(kBlack); leg->SetFillColor(0);
    leg->SetTextSize(0.037); leg->SetHeader("#scale[1.1]{p + He}","C");
    leg->AddEntry(gD,"Flight data","lep");
    leg->AddEntry(gM,"MC","lep"); leg->Draw();

    c->cd();
    TPad *pad2=new TPad("pad2_pHe","DATA / MC",0.,0.,1.,0.30);
    pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.03);
    pad2->SetTopMargin(0.0); pad2->SetBottomMargin(0.30);
    pad2->SetLogx(); pad2->SetGridx(); pad2->SetGridy();
    pad2->SetTickx(); pad2->SetTicky(); pad2->Draw(); pad2->cd();
    const double span=std::max(0.1,ratioMax-ratioMin);
    TH1F *frame2=pad2->DrawFrame(xmin, yRatioMin, xmax, yRatioMax);
    frame2->SetTitle("");
    frame2->GetYaxis()->SetLabelSize(0.073);
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetXaxis()->SetTitle("Deposited energy (GeV)");
    frame2->GetXaxis()->CenterTitle();
    frame2->GetXaxis()->SetTitleFont(42);
    frame2->GetXaxis()->SetTitleSize(0.10);
    frame2->GetXaxis()->SetTitleOffset(1.00);
    frame2->GetXaxis()->SetLabelSize(0.076);
    frame2->GetXaxis()->SetMoreLogLabels(kFALSE);
    TLine unity(xmin,1.,xmax,1.);
    unity.SetLineColor(kGray+2); unity.SetLineStyle(2);
    unity.SetLineWidth(3); unity.DrawClone("SAME");
    TLine diff(xmin,0.95,xmax,0.95);
    diff.SetLineColor(kRed); diff.SetLineStyle(2);
    diff.SetLineWidth(3); diff.DrawClone("SAME");
    gR->Draw("P SAME");
    TLatex y2; y2.SetTextAngle(90); y2.SetTextFont(42);
    y2.SetTextSize(0.076); y2.DrawLatexNDC(0.035,0.45,"DATA / MC");
    c->cd(); c->Update();

    c->SaveAs("HET_pHe_rebin_DATA_MC_10TeV_EPOSLHC.png");
    c->SaveAs("HET_pHe_rebin_DATA_MC_10TeV_EPOSLHC.pdf");
    c->SaveAs("HET_pHe_rebin_DATA_MC_10TeV_EPOSLHC.eps");
    TFile *out=TFile::Open("HET_pHe_rebin_DATA_MC_10TeV_EPOSLHC.root","RECREATE");
    if (out && !out->IsZombie()) {
        out->cd();
        for (TH1D *h : {dTot,dPass,mTot,mPass}) h->Write();
        effD->Write(); effM->Write();
        gD->Write(); gM->Write(); gR->Write(); c->Write();
        out->Close();
    }
    fData->Close(); fMC->Close();
    std::cout << "Binning nuovo: " << nNew
              << " bin. Output: HET_pHe_rebin_DATA_MC_10TeV_EPOSLHC.*" << std::endl;
}

