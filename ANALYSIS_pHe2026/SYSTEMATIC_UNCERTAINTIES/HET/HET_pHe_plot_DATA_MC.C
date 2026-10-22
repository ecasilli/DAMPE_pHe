#include <algorithm>
#include <cmath>
#include <utility>
#include <iomanip>
#include <iostream>
#include <limits>

#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"

// Confronto efficienza HET DATA/MC in E_BGO, 3 bin/decade.
// INPUT: 18 bin identici tra 10 GeV e 10 PeV.
// NON applica ulteriori rebinning; usa i conteggi originali.
// MC originale senza reweighting SBPL della composizione p/He.

namespace HET_PLOT {

bool SameBinning(const TH1D* a, const TH1D* b, double tol=1.e-8) {
    if (!a || !b || a->GetNbinsX()!=b->GetNbinsX()) return false;
    const int n=a->GetNbinsX();
    for (int k=1; k<=n+1; ++k) {
        const double ea=a->GetXaxis()->GetBinLowEdge(k);
        const double eb=b->GetXaxis()->GetBinLowEdge(k);
        if (std::fabs(ea-eb)>tol*std::max(1.,std::fabs(ea))) return false;
    }
    return true;
}

TH1D* MergeCountsByIndex(const TH1D* orig,
                         const std::vector<int>& borderIndices,
                         const std::vector<double>& newEdges,
                         const char* name) {
    const int nNew = static_cast<int>(borderIndices.size()) - 1;
    TH1D* h = new TH1D(name, orig->GetTitle(), nNew, newEdges.data());
    h->SetDirectory(nullptr);
    h->Sumw2();
    for (int j = 0; j < nNew; ++j) {
        double sum = 0., sumErr2 = 0.;
        for (int ib = borderIndices[j] + 1; ib <= borderIndices[j+1]; ++ib) {
            sum += orig->GetBinContent(ib);
            const double err = orig->GetBinError(ib);
            sumErr2 += err * err;
        }
        h->SetBinContent(j+1, sum);
        h->SetBinError(j+1, std::sqrt(sumErr2));
    }
    return h;
}


} // namespace HET_PLOT

void HET_pHe_plot_DATA_MC(
    const char *dataName="HETeff_skim_Orb118Month_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_3bins.root",
    //const char *mcName="HET_eff_MC_pHe_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_3bins.root",
    const char *mcName="HET_eff_MC_pHe_fixedMix_p70_He30_3binBGO_3binTrue.root",
    Double_t xmin=45., Double_t xmax=1.e5,
    Double_t yRatioMin=0.88, Double_t yRatioMax=1.18)
{
    if (!(xmin>0.) || !(xmax>xmin) || !(yRatioMax>yRatioMin)) {
        std::cerr << "ERRORE: intervalli del grafico non validi.\n";
        return;
    }

    TFile *fData=TFile::Open(dataName,"READ");
    TFile *fMC=TFile::Open(mcName,"READ");
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        std::cerr << "ERRORE: impossibile aprire i file:\n  DATA: " << dataName
                  << "\n  MC:   " << mcName << "\n";
        if (fData) { fData->Close(); delete fData; }
        if (fMC) { fMC->Close(); delete fMC; }
        return;
    }

    TH1D *dTot=nullptr,*dPass=nullptr,*mTot=nullptr,*mPass=nullptr;
    fData->GetObject("h1_Unb",dTot);
    fData->GetObject("h1Nobs_And",dPass);
    fMC->GetObject("h1_Unb",mTot);
    fMC->GetObject("h1Nobs_And",mPass);
    if (!dTot || !dPass || !mTot || !mPass) {
        std::cerr << "ERRORE: manca h1_Unb o h1Nobs_And.\n";
        fData->Close(); fMC->Close();
        return;
    }

    const int nOriginalBins =18; // 3 bin/decade * 6 decadi
    for (TH1D* h: {dTot,dPass,mTot,mPass}) {
        if (h->GetNbinsX()!=nOriginalBins || !HET_PLOT::SameBinning(dTot,h)) {
            std::cerr << "ERRORE: richiesti 18 bin BGO coincidenti DATA/MC; problema in "
                      << h->GetName() << " (" << h->GetNbinsX() << " bin).\n";
            fData->Close(); fMC->Close();
            return;
        }
    }
    const double eMin=dTot->GetXaxis()->GetBinLowEdge(1);
    const double eMax=dTot->GetXaxis()->GetBinUpEdge(nOriginalBins);
    if (std::fabs(eMin/10.-1.)>1.e-4 || std::fabs(eMax/1.e7-1.)>1.e-4) {
        std::cerr << "ERRORE: intervallo istogrammi diverso da 10 GeV - 10 PeV.\n";
        fData->Close(); fMC->Close();
        return;
    }

    // 3 bin/decade: 10 TeV = bordo 9, 100 TeV = bordo 12.
    // Eliminiamo SOLO i bordi interni 10 e 11: 3 bin --> 1 bin.
    // Restano 16 bin (18 - 2), di cui uno unico tra 10 e 100 TeV.
    const std::vector<int> edgeIndex = {
        0,1,2,3,4,5,6,7,8,9,12,13,14,15,16,17,18
    };
    const auto* ax = dTot->GetXaxis();
    const double edge10TeV = ax->GetBinLowEdge(10); // bordo indice 9
    const double edge100TeV = ax->GetBinLowEdge(13); // bordo indice 12
    if (std::fabs(edge10TeV / 1.e4 - 1.) > 1.e-4 ||
        std::fabs(edge100TeV / 1.e5 - 1.) > 1.e-4) {
        std::cerr << "ERRORE: bordi 10 e 100 TeV incompatibili.\n";
        fData->Close(); fMC->Close();
        return;
    }
    std::vector<double> edges;
    edges.reserve(edgeIndex.size());
    for (int ix : edgeIndex) edges.push_back(ax->GetBinLowEdge(ix+1));
    TH1D *dTotMerged  = HET_PLOT::MergeCountsByIndex(dTot, edgeIndex, edges, "hDATA_Unb_merge");
    TH1D *dPassMerged = HET_PLOT::MergeCountsByIndex(dPass,edgeIndex, edges, "hDATA_And_merge");
    TH1D *mTotMerged  = HET_PLOT::MergeCountsByIndex(mTot, edgeIndex, edges, "hMC_Unb_merge");
    TH1D *mPassMerged = HET_PLOT::MergeCountsByIndex(mPass,edgeIndex, edges, "hMC_And_merge");
    dTot = dTotMerged;
    dPass = dPassMerged;
    mTot = mTotMerged;
    mPass = mPassMerged;
    const int nBins = static_cast<int>(edgeIndex.size()) - 1;

    if (!TEfficiency::CheckConsistency(*dPass,*dTot)) {
        std::cerr << "ERRORE: conteggi DATA incompatibili con TEfficiency.\n";
        fData->Close(); fMC->Close();
        return;
    }
    if (!TEfficiency::CheckConsistency(*mPass,*mTot,"w")) {
        std::cerr << "ERRORE: conteggi MC pesati incompatibili con TEfficiency.\n";
        fData->Close(); fMC->Close();
        return;
    }

    TEfficiency* effD=new TEfficiency(*dPass,*dTot);
    effD->SetName("HET_Eff_DATA_3bins");
    effD->SetStatisticOption(TEfficiency::kFCP);
    effD->SetConfidenceLevel(0.682689492);

    TEfficiency* effM=new TEfficiency(*mPass,*mTot);
    effM->SetName("HET_Eff_MC_original_3bins");
    effM->SetUseWeightedEvents();
    effM->SetStatisticOption(TEfficiency::kFNormal);
    effM->SetConfidenceLevel(0.682689492);

    TGraphAsymmErrors* gD=new TGraphAsymmErrors();
    TGraphAsymmErrors* gM=new TGraphAsymmErrors();
    TGraphAsymmErrors* gR=new TGraphAsymmErrors();
    gD->SetName("g_HET_DATA");
    gM->SetName("g_HET_MC");
    gR->SetName("g_DATA_over_MC");

    TH1D* hNeffMC=(TH1D*)mTot->Clone("h_MC_Neff_UBT");
    hNeffMC->SetDirectory(nullptr);
    hNeffMC->Reset("ICES");
    hNeffMC->SetTitle("Effective MC UBT statistics;Deposited BGO energy [GeV];N_{eff}");

    std::cout << "\nDATA: " << dataName << "\nMC:   " << mcName
              << "\nBinning BGO: " << nBins << " bin = 3 bin/decade, nessun rebin"
              << "\nIntervallo plot: " << xmin << " - " << xmax << " GeV\n\n";
    std::cout << std::setw(12) << "Elo(GeV)" << std::setw(12) << "Ehi(GeV)"
              << std::setw(12) << "DATA_UBT" << std::setw(13) << "MC_Neff"
              << std::setw(13) << "eff_DATA" << std::setw(13) << "eff_MC"
              << std::setw(13) << "D/MC" << "\n";

    int nD=0,nM=0,nR=0;
    for (int ib=1; ib<=nBins; ++ib) {
        const double lo=dTot->GetXaxis()->GetBinLowEdge(ib);
        const double hi=dTot->GetXaxis()->GetBinUpEdge(ib);
        const double x=std::sqrt(lo*hi);
        const double exl=x-lo, exh=hi-x;
        const double nd=dTot->GetBinContent(ib);
        const double nm=mTot->GetBinContent(ib);
        const double sumW2=std::pow(mTot->GetBinError(ib),2.);
        const double neff=(nm>0. && sumW2>0.) ? nm*nm/sumW2 : 0.;
        hNeffMC->SetBinContent(ib,neff);

        const bool okD=nd>0., okM=nm>0.;
        double d=0.,dl=0.,du=0.,m=0.,ml=0.,mu=0.;
        if (okD) {
            d=effD->GetEfficiency(ib);
            dl=effD->GetEfficiencyErrorLow(ib);
            du=effD->GetEfficiencyErrorUp(ib);
        }
        if (okM) {
            m=effM->GetEfficiency(ib);
            ml=effM->GetEfficiencyErrorLow(ib);
            mu=effM->GetEfficiencyErrorUp(ib);
        }

        double r=std::numeric_limits<double>::quiet_NaN();
        if (okD && okM && m>0.) r=d/m;
        // Primo bin escluso dal grafico: centro sotto 17 GeV e taglio BGO >20 GeV.
        // Eventuali bin senza statistica restano esclusi (no punti artificiali).
        const bool draw=(x>=xmin && x<=xmax);
        if (draw && okD && std::isfinite(d) && std::isfinite(dl) && std::isfinite(du)) {
            gD->SetPoint(nD,x,d);
            gD->SetPointError(nD,exl,exh,dl,du);
            ++nD;
        }
        if (draw && okM && std::isfinite(m) && std::isfinite(ml) && std::isfinite(mu)) {
            gM->SetPoint(nM,x,m);
            gM->SetPointError(nM,exl,exh,ml,mu);
            ++nM;
        }
        if (draw && std::isfinite(r) && std::isfinite(dl) && std::isfinite(du)
                 && std::isfinite(ml) && std::isfinite(mu)) {
            // Errori DATA/MC indipendenti, propagazione asimmetrica di primo ordine.
            const double rl=std::hypot(dl/m,d*mu/(m*m));
            const double ru=std::hypot(du/m,d*ml/(m*m));
            if (std::isfinite(rl) && std::isfinite(ru)) {
                gR->SetPoint(nR,x,r);
                gR->SetPointError(nR,exl,exh,rl,ru);
                ++nR;
                if (r-rl<yRatioMin || r+ru>yRatioMax)
                    std::cout << "[INFO] Rapporto o errore fuori range Y per "
                              << lo << " - " << hi << " GeV."
                              << " Regolare yRatioMin/yRatioMax se necessario.\n";
            }
        }

        if (draw && okD && nd<20.)
            std::cerr << "[WARNING] Statistica UBT DATA bassa tra "
                      << lo << " e " << hi << " GeV: " << nd << " eventi.\n";
        if (draw && neff>0. && neff<20.)
            std::cerr << "[WARNING] MC Neff basso tra "
                      << lo << " e " << hi << " GeV: " << neff << ".\n";

        if (hi<=xmax*(1.+1.e-8)) {
            std::cout << std::scientific << std::setprecision(3)
                      << std::setw(12) << lo << std::setw(12) << hi
                      << std::fixed << std::setprecision(1)
                      << std::setw(12) << nd << std::setw(13) << neff
                      << std::setprecision(4)
                      << std::setw(13) << (okD?d:-1.)
                      << std::setw(13) << (okM?m:-1.)
                      << std::setw(13) << r << "\n";
        }
    }
    std::cout << "\nPunti nel plot: DATA=" << nD
              << ", MC=" << nM << ", DATA/MC=" << nR << "\n";

    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gD->SetLineColor(kRed+1);  gD->SetMarkerColor(kRed+1);
    gD->SetMarkerStyle(20);    gD->SetMarkerSize(1.1); gD->SetLineWidth(2);
    gM->SetLineColor(kBlue+1); gM->SetMarkerColor(kBlue+1);
    gM->SetMarkerStyle(24);    gM->SetMarkerSize(1.1); gM->SetLineWidth(2);
    gR->SetLineColor(kBlack);  gR->SetMarkerColor(kBlack);
    gR->SetMarkerStyle(21);    gR->SetMarkerSize(1.0); gR->SetLineWidth(2);

    TCanvas* c=new TCanvas("c_HET_pHe","HET efficiency p+He",1000,950);
    TPad* pad1=new TPad("pad1_HET","efficiencies",0.,0.30,1.,1.);
    pad1->SetLeftMargin(0.10); pad1->SetRightMargin(0.03);
    pad1->SetBottomMargin(0.0); pad1->SetTopMargin(0.04);
    pad1->SetLogx(); pad1->SetGridx(); pad1->SetGridy();
    pad1->SetTickx(); pad1->SetTicky(); pad1->Draw(); pad1->cd();
    TH1F* frame1=pad1->DrawFrame(xmin,0.05,xmax,1.06);
    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0.);
    frame1->GetYaxis()->SetLabelSize(0.034);
    frame1->GetYaxis()->SetNdivisions(505);
    gM->Draw("P SAME"); gD->Draw("P SAME");
    TLatex y1; y1.SetTextAngle(90); y1.SetTextFont(42);
    y1.SetTextSize(0.037); y1.DrawLatexNDC(0.035,0.55,"HET efficiency");
    TLegend* leg=new TLegend(0.19,0.56,0.44,0.76);
    leg->SetLineColor(kBlack); leg->SetFillColor(0);
    leg->SetTextSize(0.037); leg->SetHeader("#scale[1.1]{p + He}","C");
    leg->AddEntry(gD,"Flight data","lep");
    leg->AddEntry(gM,"MC (70-30)","lep"); leg->Draw();

    c->cd();
    TPad* pad2=new TPad("pad2_HET","DATA / MC",0.,0.,1.,0.30);
    pad2->SetLeftMargin(0.10); pad2->SetRightMargin(0.03);
    pad2->SetTopMargin(0.0); pad2->SetBottomMargin(0.30);
    pad2->SetLogx(); pad2->SetGridx(); pad2->SetGridy();
    pad2->SetTickx(); pad2->SetTicky(); pad2->Draw(); pad2->cd();
    TH1F* frame2=pad2->DrawFrame(xmin,yRatioMin,xmax,yRatioMax);
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
    TLine ref95(xmin,0.965,xmax,0.965);
    ref95.SetLineColor(kRed); ref95.SetLineStyle(2);
    ref95.SetLineWidth(3); ref95.DrawClone("SAME");
    gR->Draw("P SAME");
    TLatex y2; y2.SetTextAngle(90); y2.SetTextFont(42);
    y2.SetTextSize(0.076); y2.DrawLatexNDC(0.035,0.45,"DATA / MC");
    c->cd(); c->Update();

    const TString outBase="HET_pHe_DATA_MC_70-30_3bins_merge10_100TeV";
    c->SaveAs(outBase+".png");
    c->SaveAs(outBase+".pdf");
    //c->SaveAs(outBase+".eps");
    TFile* fout=TFile::Open(outBase+".root","RECREATE");
    if (fout && !fout->IsZombie()) {
        fout->cd();
        // Copie indipendenti per conservare i conteggi originali nel file di output.
        for (auto item : {std::pair<const char*,TH1D*>("hDATA_Unb",dTot),
                          std::pair<const char*,TH1D*>("hDATA_And",dPass),
                          std::pair<const char*,TH1D*>("hMC_Unb",mTot),
                          std::pair<const char*,TH1D*>("hMC_And",mPass)}) {
            TH1D* tmp=(TH1D*)item.second->Clone(item.first);
            tmp->Write();
            delete tmp;
        }
        hNeffMC->Write();
        effD->Write(); effM->Write();
        gD->Write(); gM->Write(); gR->Write(); c->Write();
        fout->Close(); delete fout;
    } else {
        std::cerr << "WARNING: impossibile creare file ROOT di output.\n";
        if (fout) { fout->Close(); delete fout; }
    }
    fData->Close(); fMC->Close();
    delete fData; delete fMC;
    std::cout << "Output: " << outBase << ".{png,pdf,eps,root}\n";
}

