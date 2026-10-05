
// plot_pHe_systematics.C
//
// ROOT macro per disegnare le incertezze relative del flusso p+He
// a partire da sys_summary.txt.
//
// Uso:
//   root -l
//   .x plot_pHe_systematics.C
//
// oppure:
//   root -l -q 'plot_pHe_systematics.C("sys_summary.txt")'
//
// Nota importante:
//   - tutte le Sys_* nel file sono gia' FRAZIONI relative;
//   - Stat_err e' invece l'errore ASSOLUTO sul flusso,
//     quindi viene convertito in Stat_err / Flux.

#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TColor.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

struct SystBin {
    double Emin;
    double Emax;
    double Emean;
    double Flux;
    double Stat_err;
    double Sys_HET;
    double Sys_STKtrack;
    double Sys_ChargeSel;
    double Sys_unfolding;
    double Sys_pHe_mix;
    double Sys_ana;
    double Sys_had;
    double Sys_tot;
};

void SetHistStyle(TH1D *h, Color_t color, Style_t style = 1, Width_t width = 2)
{
    h->SetLineColor(color);
    h->SetLineStyle(style);
    h->SetLineWidth(width);
    h->SetFillStyle(0);
}

void plot_pHe_systematics(const char *filename = "sys_summary_10TeV_EPOSLHC.txt")
{
    // ------------------------------------------------------------
    // 1. Lettura del file
    // ------------------------------------------------------------
    std::ifstream fin(filename);

    if (!fin.is_open()) {
        std::cerr << "ERRORE: non riesco ad aprire il file: "
                  << filename << std::endl;
        return;
    }

    std::vector<SystBin> bins;

    std::string line;
    while (std::getline(fin, line)) {

        // Salta righe vuote
        if (line.empty())
            continue;

        // Salta commenti/header
        std::size_t first = line.find_first_not_of(" \t");
        if (first == std::string::npos)
            continue;

        if (line[first] == '#')
            continue;

        std::istringstream iss(line);

        SystBin b;

        if (!(iss >> b.Emin
                  >> b.Emax
                  >> b.Emean
                  >> b.Flux
                  >> b.Stat_err
                  >> b.Sys_HET
                  >> b.Sys_STKtrack
                  >> b.Sys_ChargeSel
                  >> b.Sys_unfolding
                  >> b.Sys_pHe_mix
                  >> b.Sys_ana
                  >> b.Sys_had
                  >> b.Sys_tot)) {

            std::cerr << "ATTENZIONE: riga non riconosciuta, la salto:\n"
                      << line << std::endl;
            continue;
        }

        bins.push_back(b);
    }

    fin.close();

    if (bins.empty()) {
        std::cerr << "ERRORE: nessun bin letto dal file." << std::endl;
        return;
    }

    const int nBins = static_cast<int>(bins.size());

    // ------------------------------------------------------------
    // 2. Bordi energetici dei bin
    //
    //    TH1D a binning variabile permette di ottenere direttamente
    //    l'andamento "a scalini" come nella figura di riferimento.
    // ------------------------------------------------------------
    std::vector<double> edges(nBins + 1);

    edges[0] = bins[0].Emin;
    for (int i = 0; i < nBins; ++i)
        edges[i + 1] = bins[i].Emax;

    // ------------------------------------------------------------
    // 3. Istogrammi
    // ------------------------------------------------------------
    TH1D *hHET       = new TH1D("hHET",       "", nBins, edges.data());
    TH1D *hSTK       = new TH1D("hSTK",       "", nBins, edges.data());
    TH1D *hCharge    = new TH1D("hCharge",    "", nBins, edges.data());
    TH1D *hUnfold    = new TH1D("hUnfold",    "", nBins, edges.data());
    TH1D *hPHeMix    = new TH1D("hPHeMix",    "", nBins, edges.data());
    TH1D *hStat      = new TH1D("hStat",      "", nBins, edges.data());
    TH1D *hAna       = new TH1D("hAna",       "", nBins, edges.data());
    TH1D *hHad       = new TH1D("hHad",       "", nBins, edges.data());
    TH1D *hTot       = new TH1D("hTot",       "", nBins, edges.data());

    double maxDiffAna = 0.0;
    double maxDiffTot = 0.0;

    for (int i = 0; i < nBins; ++i) {

        const int ibin = i + 1;
        const SystBin &b = bins[i];

        // Le sistematiche sono gia' frazioni relative.
        hHET    ->SetBinContent(ibin, b.Sys_HET);
        hSTK    ->SetBinContent(ibin, b.Sys_STKtrack);
        hCharge ->SetBinContent(ibin, b.Sys_ChargeSel);
        hUnfold ->SetBinContent(ibin, b.Sys_unfolding);
        hPHeMix ->SetBinContent(ibin, b.Sys_pHe_mix);
        hAna    ->SetBinContent(ibin, b.Sys_ana);
        hHad    ->SetBinContent(ibin, b.Sys_had);
        hTot    ->SetBinContent(ibin, b.Sys_tot);

        // Stat_err e' assoluto sul flusso -> errore statistico relativo.
        const double statRelative = (b.Flux != 0.0) ? b.Stat_err / b.Flux : 0.0;

        hStat->SetBinContent(ibin, statRelative);

        // --------------------------------------------------------
        // Check delle somme in quadratura
        // --------------------------------------------------------
        const double anaCheck = std::sqrt(
              b.Sys_HET       * b.Sys_HET
            + b.Sys_STKtrack  * b.Sys_STKtrack
            + b.Sys_ChargeSel * b.Sys_ChargeSel
            + b.Sys_unfolding * b.Sys_unfolding
            + b.Sys_pHe_mix   * b.Sys_pHe_mix
        );

        const double totCheck = std::sqrt( b.Sys_ana * b.Sys_ana + b.Sys_had * b.Sys_had );

        maxDiffAna = std::max(maxDiffAna, std::fabs(anaCheck - b.Sys_ana));
        maxDiffTot = std::max(maxDiffTot, std::fabs(totCheck - b.Sys_tot));
    }

    std::cout << "\nNumero di bin letti: " << nBins << std::endl;
    //std::cout << "Max |Sys_ana - quadrature check| = " << maxDiffAna << std::endl;
    //std::cout << "Max |Sys_tot - quadrature check| = " << maxDiffTot << std::endl;

    // ------------------------------------------------------------
    // 4. Stile delle curve
    // ------------------------------------------------------------
    SetHistStyle(hHET,    kBlue,     1, 2);
    SetHistStyle(hSTK,    kMagenta,  1, 2);
    SetHistStyle(hCharge, kGreen,    1, 2);
    SetHistStyle(hUnfold, kOrange,   1, 2);
    SetHistStyle(hPHeMix, kCyan+1,   1, 2);

    SetHistStyle(hStat,   kRed,      1, 2);

    // Totale delle sole sistematiche di analisi: dash-dot
    SetHistStyle(hAna,    kBlack,   10, 2);

    // Hadronic: grigio
    SetHistStyle(hHad,    kGray+1,   1, 2);

    // Totale finale: nero, un po' piu' spesso
    SetHistStyle(hTot,    kBlack,    1, 3);

    // ------------------------------------------------------------
    // 5. Canvas
    // ------------------------------------------------------------
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_pHe_systematics","p+He relative uncertainties",800,700);

    c->SetLogx();
    c->SetLogy();

    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13);
    c->SetTopMargin(0.06);

    c->SetTicks(1, 1);

    // ------------------------------------------------------------
    // 6. Assi
    //
    // hTot viene disegnato per primo e quindi definisce il frame.
    // ------------------------------------------------------------
    hTot->SetTitle("");

    hTot->GetXaxis()->SetTitle("Kinetic energy (GeV)");
    hTot->GetYaxis()->SetTitle("Relative uncertainty");

    hTot->GetXaxis()->SetTitleSize(0.040);
    hTot->GetYaxis()->SetTitleSize(0.040);

    hTot->GetXaxis()->SetLabelSize(0.040);
    hTot->GetYaxis()->SetLabelSize(0.040);

    hTot->GetXaxis()->SetTitleOffset(1.20);
    hTot->GetYaxis()->SetTitleOffset(1.30);

    hTot->GetXaxis()->SetMoreLogLabels(false);
    hTot->GetXaxis()->SetNoExponent(false);

    // Modifica questi due numeri se vuoi cambiare il range verticale.
    hTot->SetMinimum(1.0e-3);
    hTot->SetMaximum(2.0);

    // ------------------------------------------------------------
    // 7. Disegno
    // ------------------------------------------------------------
    hTot->Draw("HIST");

    hHET    ->Draw("HIST SAME");
    hSTK    ->Draw("HIST SAME");
    hCharge ->Draw("HIST SAME");
    hUnfold ->Draw("HIST SAME");
    hPHeMix ->Draw("HIST SAME");

    hStat ->Draw("HIST SAME");

    hAna  ->Draw("HIST SAME");
    hHad  ->Draw("HIST SAME");

    // Ridisegnato per tenerlo visivamente in primo piano.
    hTot  ->Draw("HIST SAME");

    // ------------------------------------------------------------
    // 8. Legenda a due colonne, nello stile della figura allegata
    // ------------------------------------------------------------
    TLegend *leg = new TLegend(0.17, 0.71, 0.91, 0.92);

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    leg->SetNColumns(2);
    leg->SetColumnSeparation(0.05);

    leg->AddEntry(hHET,    "HE trigger",                    "l");
    leg->AddEntry(hStat,   "Statistical",                   "l");

    leg->AddEntry(hSTK,    "STK track",                     "l");
    leg->AddEntry(hAna,    "Total systematics analysis",    "l");

    leg->AddEntry(hCharge, "Charge selection",              "l");
    leg->AddEntry(hHad,    "Hadronic",                      "l");

    leg->AddEntry(hUnfold, "Unfolding",                     "l");
    leg->AddEntry(hTot,    "Total systematics ana. #oplus had.", "l");

    leg->AddEntry(hPHeMix, "p+He mixture",                  "l");

    leg->Draw();

    // ------------------------------------------------------------
    // 9. Label p+He
    // ------------------------------------------------------------
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(62);
    latex.SetTextSize(0.050);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.93, 0.18, "p+He");

    c->RedrawAxis();

    // ------------------------------------------------------------
    // 10. Output
    // ------------------------------------------------------------
    c->SaveAs("pHe_relative_uncertainties_summary_10TeV_EPOSLHC.png");
    c->SaveAs("pHe_relative_uncertainties_summary_10TeV_EPOSLHC.pdf");

    std::cout << "\nSalvati:" << std::endl;
    std::cout << "  pHe_relative_uncertainties_summary_10TeV_EPOSLHC.png" << std::endl;
    std::cout << "  pHe_relative_uncertainties_summary_10TeV_EPOSLHC.pdf" << std::endl;
}
