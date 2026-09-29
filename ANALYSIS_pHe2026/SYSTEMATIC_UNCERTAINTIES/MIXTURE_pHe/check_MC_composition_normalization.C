// Diagnostica: verifica se il flusso SBPL alto deriva da una diversa
// normalizzazione globale delle matrici MC.
// Esecuzione: root -l -b -q check_MC_composition_normalization.C

#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

void check_MC_composition_normalization()
{
    gStyle->SetOptStat(0);

    // =========================================================
    // NOMI FILE: MODIFICARE QUI (due ROOT contenenti le matrici).
    // =========================================================
    const char *originalRoot = "../../ROOT_FILES/PHe_MC_p_He_5PeV_5binperdecade_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_STKvert0e7_24sett26.root";
    const char *sbplRoot = "ROOT_FILES/PHe_MC_p_He_5PeV_5binperdecade_3sLow_6Up_PSDprogr_STKch450_comb_vert0e7_SBPLmix.root";

    TFile *f0 = TFile::Open(originalRoot, "READ");
    TFile *f1 = TFile::Open(sbplRoot, "READ");
    if (!f0 || f0->IsZombie() || !f1 || f1->IsZombie()) {
        cerr << "ERRORE: impossibile leggere uno dei ROOT di ingresso." << endl;
        if (f0) { f0->Close(); delete f0; }
        if (f1) { f1->Close(); delete f1; }
        return;
    }

    TH2 *all0 = dynamic_cast<TH2 *>(f0->Get("h2Ntrig_wgt_all"));
    TH2 *end0 = dynamic_cast<TH2 *>(f0->Get("h2Ntrig_wgt"));
    TH2 *all1 = dynamic_cast<TH2 *>(f1->Get("h2Ntrig_wgt_all"));
    TH2 *end1 = dynamic_cast<TH2 *>(f1->Get("h2Ntrig_wgt"));

    if (!all0 || !end0 || !all1 || !end1) {
        cerr << "ERRORE: mancano h2Ntrig_wgt_all o h2Ntrig_wgt." << endl;
        f0->Close(); f1->Close(); delete f0; delete f1;
        return;
    }

    const int nBins = all0->GetNbinsX();
    if (all1->GetNbinsX() != nBins || end0->GetNbinsX() != nBins ||
        end1->GetNbinsX() != nBins ||
        all0->GetNbinsY() != all1->GetNbinsY()) {
        cerr << "ERRORE: binning delle matrici incompatibile." << endl;
        f0->Close(); f1->Close(); delete f0; delete f1;
        return;
    }
    for (int i = 1; i <= nBins; ++i) {
        const double a = all0->GetXaxis()->GetBinLowEdge(i);
        const double b = all1->GetXaxis()->GetBinLowEdge(i);
        if (std::fabs(a - b) > 1e-5 * a) {
            cerr << "ERRORE: bin di energia vera diversi al bin " << i << endl;
            f0->Close(); f1->Close(); delete f0; delete f1;
            return;
        }
    }

    // Proiezione sull'energia VERA, sommando SOLO i bin ricostruiti interni
    // (i.e. senza underflow/overflow). Trattamento uguale per entrambi.
    TH1D *projAll0 = all0->ProjectionX("projAll0_cmp", 1, all0->GetNbinsY(), "e");
    TH1D *projAll1 = all1->ProjectionX("projAll1_cmp", 1, all1->GetNbinsY(), "e");
    TH1D *projEnd0 = end0->ProjectionX("projEnd0_cmp", 1, end0->GetNbinsY(), "e");
    TH1D *projEnd1 = end1->ProjectionX("projEnd1_cmp", 1, end1->GetNbinsY(), "e");
    const TH1D *fractionP = dynamic_cast<TH1D *>(f1->Get("h_ab_SBPL_p_used"));
    const TH1D *fractionHe = dynamic_cast<TH1D *>(f1->Get("h_ab_SBPL_He_used"));

    vector<double> x, ratioAll, ratioFinal, ratioEfficiency;
    gSystem->mkdir("TXT_FILES", kTRUE);
    ofstream out("TXT_FILES/pHe_MC_composition_normalization_check.dat");
    if (!out) {
        cerr << "ERRORE: impossibile aprire il file di diagnostica in output." << endl;
        f0->Close(); f1->Close(); delete f0; delete f1;
        return;
    }
    // final/all e' soltanto la frazione di selezione all'interno del range
    // ricostruito salvato negli istogrammi: NON e' l'accettanza assoluta.
    out << "# EtrueGeV originalAll mixAll mixOverOriginalAll "
           "originalFinal mixFinal mixOverOriginalFinal "
           "effOriginal effMix effMixOverOriginal fp fHe\n";
    out << setprecision(12);

    cout << "\n# Etrue[GeV]  mix/orig(all)  mix/orig(final) "
            " eff_mix/eff_orig  f_p  f_He\n";
    for (int i = 1; i <= nBins; ++i) {
        const double eLo = all0->GetXaxis()->GetBinLowEdge(i);
        const double eHi = all0->GetXaxis()->GetBinUpEdge(i);
        const double E = sqrt(eLo * eHi);
        const double a0 = projAll0->GetBinContent(i);
        const double a1 = projAll1->GetBinContent(i);
        const double z0 = projEnd0->GetBinContent(i);
        const double z1 = projEnd1->GetBinContent(i);
        const double fP = fractionP ? fractionP->GetBinContent(i) : -1.0;
        const double fHe = fractionHe ? fractionHe->GetBinContent(i) : -1.0;

        const double rAll = a0 > 0. ? a1 / a0 : -1.;
        const double rFinal = z0 > 0. ? z1 / z0 : -1.;
        const double eff0 = a0 > 0. ? z0 / a0 : -1.;
        const double eff1 = a1 > 0. ? z1 / a1 : -1.;
        const double rEff = eff0 > 0. ? eff1 / eff0 : -1.;

        out << E << ' ' << a0 << ' ' << a1 << ' ' << rAll << ' '
            << z0 << ' ' << z1 << ' ' << rFinal << ' '
            << eff0 << ' ' << eff1 << ' ' << rEff << ' '
            << fP << ' ' << fHe << '\n';

        if (rAll <= 0. || rFinal <= 0. || rEff <= 0.) continue;
        x.push_back(E);
        ratioAll.push_back(rAll);
        ratioFinal.push_back(rFinal);
        ratioEfficiency.push_back(rEff);

        cout << E << "    " << rAll << "    " << rFinal << "    "
             << rEff << "    " << fP << "    " << fHe << '\n';
    }
    out.close();

    if (!x.empty()) {
        TCanvas *c = new TCanvas("c_mix_diagnostics", "MC mix: normalization", 1050, 870);
        TPad *upper = new TPad("pad_norm_up", "", 0., 0.43, 1., 1.);
        TPad *lower = new TPad("pad_norm_dn", "", 0., 0., 1., 0.43);
        upper->SetBottomMargin(0.03); upper->SetLeftMargin(0.13);
        lower->SetTopMargin(0.03); lower->SetBottomMargin(0.19); lower->SetLeftMargin(0.13);
        upper->Draw(); lower->Draw();

        double maxR = 1.0;
        for (size_t i = 0; i < x.size(); ++i) {
            maxR = maxR > ratioAll[i] ? maxR : ratioAll[i];
            maxR = maxR > ratioFinal[i] ? maxR : ratioFinal[i];
        }
        upper->cd(); upper->SetLogx(); upper->SetGrid();
        TH1F *frameUp = upper->DrawFrame(x.front()/1.2, 0.0, x.back()*1.2, 1.25*maxR);
        frameUp->SetTitle("; ;(MC SBPL) / (MC originale)");
        frameUp->GetXaxis()->SetLabelSize(0);
        TGraph *grAll = new TGraph((int)x.size(), x.data(), ratioAll.data());
        TGraph *grFinal = new TGraph((int)x.size(), x.data(), ratioFinal.data());
        grAll->SetMarkerStyle(20); grAll->SetMarkerColor(kBlue+1); grAll->SetLineColor(kBlue+1);
        grFinal->SetMarkerStyle(21); grFinal->SetMarkerColor(kRed+1); grFinal->SetLineColor(kRed+1);
        grAll->Draw("PL SAME"); grFinal->Draw("PL SAME");
        TLine *unityUp = new TLine(x.front()/1.2, 1., x.back()*1.2, 1.);
        unityUp->SetLineStyle(2); unityUp->Draw();
        TLegend *leg = new TLegend(0.55, 0.72, 0.92, 0.88);
        leg->AddEntry(grAll, "h2Ntrig_wgt_all", "pl");
        leg->AddEntry(grFinal, "h2Ntrig_wgt", "pl"); leg->Draw();

        lower->cd(); lower->SetLogx(); lower->SetGrid();
        double maxEff = 1.0;
        double minEff = 1.0;
        for (size_t i = 0; i < ratioEfficiency.size(); ++i) {
            if (ratioEfficiency[i] > maxEff) maxEff = ratioEfficiency[i];
            if (ratioEfficiency[i] < minEff) minEff = ratioEfficiency[i];
        }
        double span = (maxEff - minEff) > 0.05 ? (maxEff - minEff) : 0.05;
        TH1F *frameDn = lower->DrawFrame(x.front()/1.2,
             minEff - 0.25*span, x.back()*1.2, maxEff + 0.25*span);
        frameDn->SetTitle(";MC true energy [GeV];(final/all) SBPL / original");
        frameDn->GetXaxis()->SetTitleSize(0.07);
        frameDn->GetXaxis()->SetLabelSize(0.06);
        frameDn->GetYaxis()->SetTitleSize(0.065);
        frameDn->GetYaxis()->SetLabelSize(0.055);
        frameDn->GetYaxis()->SetTitleOffset(0.9);
        TGraph *grEff = new TGraph((int)x.size(), x.data(), ratioEfficiency.data());
        grEff->SetMarkerStyle(20); grEff->SetMarkerColor(kGreen+2);
        grEff->SetLineColor(kGreen+2);
        grEff->Draw("PL SAME");
        TLine *unityDn = new TLine(x.front()/1.2, 1., x.back()*1.2, 1.);
        unityDn->SetLineStyle(2); unityDn->Draw();

        gSystem->mkdir("PLOTS", kTRUE);
        c->SaveAs("PLOTS/pHe_MC_composition_normalization_check.pdf");
        c->SaveAs("PLOTS/pHe_MC_composition_normalization_check.png");
    }

    cout << "\nDiagnostica salvata in "
            "TXT_FILES/pHe_MC_composition_normalization_check.dat\n";
    f0->Close(); f1->Close(); delete f0; delete f1;
}

