#include <TCanvas.h>
#include <TColor.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

void draw_pHe_fluxes_perYear()
{
    // Cartella contenente i file .dat (cambiala se necessario).
    const std::string directory = "TXT_FILES/";
    const std::string suffix =
        "_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC.dat";

    const int primoAnno = 2016;
    const int ultimoAnno = 2025;
    const int scartaInizio = 3;
    const int scartaFine = 5;

    // Flusso rappresentato: E^(-2.6) * Phi(E).
    // Per usare invece E^(+2.6) * Phi(E), imposta exponentE = +2.6.
    const double exponentE = 2.6;

    // Palette categoriale: 10 colori distinti e leggibili su fondo bianco.
    const char* hex[10] = {
        "#4E79A7", // 2016 blu
        "#F28E2B", // 2017 arancio
        "#59A14F", // 2018 verde
        "#E15759", // 2019 rosso
        "#9467BD", // 2020 viola
        "#1E9EAE", // 2021 turchese
        "#9C755F", // 2022 marrone
        "#B58900", // 2023 ocra
        "#CC6677", // 2024 rosa scuro
        "#333333"  // 2025 antracite
    };
    const int markers[10] = {20, 21, 22, 23, 29, 33, 34, 24, 25, 26};

    struct Punto {
        double E, flux, errFlux, errE, Emin, Emax;
    };

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(3);

    std::vector<std::pair<int, TGraphErrors*>> grafici;
    double xmin = std::numeric_limits<double>::infinity();
    double xmax = 0.0;
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = 0.0;

    for (int anno = primoAnno; anno <= ultimoAnno; ++anno) {
        const std::string path = directory + "flux_spectrum_pHe_"
                               + std::to_string(anno) + suffix;
        std::ifstream file(path);
        if (!file.is_open()) {
            std::cerr << "ATTENZIONE: file non trovato: " << path << '\n';
            continue;
        }

        std::vector<Punto> dati;
        std::string riga;
        while (std::getline(file, riga)) {
            if (riga.find_first_not_of(" \t\r\n") == std::string::npos) continue;
            std::istringstream ss(riga);
            Punto p;
            if (ss >> p.E >> p.flux >> p.errFlux >> p.errE >> p.Emin >> p.Emax)
                dati.push_back(p);
        }

        if (dati.size() <= static_cast<size_t>(scartaInizio + scartaFine)) {
            std::cerr << "ATTENZIONE: troppo pochi bin nel file " << path << '\n';
            continue;
        }

        auto* g = new TGraphErrors();
        g->SetName(("gFlux_" + std::to_string(anno)).c_str());
        const int idx = anno - primoAnno;
        const int col = TColor::GetColor(hex[idx]);
        g->SetLineColor(col);
        g->SetMarkerColor(col);
        g->SetLineWidth(2);
        g->SetMarkerStyle(markers[idx]);
        g->SetMarkerSize(1.0);

        // Indici mantenuti: 3 ... N-6 (estremi inclusi).
        for (size_t i = scartaInizio; i < dati.size() - scartaFine; ++i) {
            const Punto& p = dati[i];
            if (p.E <= 0 || p.flux <= 0) continue; // incompatibili con assi log

            // Scala identica per il flusso e l'errore statistico.
            const double factor = std::pow(p.E, exponentE);
            const double fluxScaled = p.flux * factor;
            const double errFluxScaled = std::max(0.0, p.errFlux) * factor;

            const int n = g->GetN();
            g->SetPoint(n, p.E, fluxScaled);
            g->SetPointError(n, std::max(0.0, p.errE), errFluxScaled);

            xmin = std::min(xmin, p.Emin);
            xmax = std::max(xmax, p.Emax);
            ymin = std::min(ymin, fluxScaled);
            ymax = std::max(ymax, fluxScaled);
        }

        if (g->GetN() == 0) {
            delete g;
            std::cerr << "ATTENZIONE: nessun punto positivo per " << anno << '\n';
            continue;
        }
        grafici.emplace_back(anno, g);
        std::cout << anno << ": " << g->GetN() << " punti plottati\n";
    }

    if (grafici.empty()) {
        std::cerr << "ERRORE: nessun flusso da plottare.\n";
        return;
    }

    auto* c = new TCanvas("c_flux", "p+He flux 2016-2025", 1200, 800);
    c->SetLogx();
    //c->SetLogy();
    c->SetTicks(1, 1);
    c->SetLeftMargin(0.13);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.07);

    auto* frame = c->DrawFrame(xmin * 0.85, 3.5e3, xmax * 1.15, 16.5e3);
    //auto* frame = c->DrawFrame(xmin * 0.85, 7.2e3, xmax * 1.15, 14.8e3);
    frame->SetTitle("p+He flux by year (2016-2025)");
    frame->GetXaxis()->SetTitle("Kinetic energy [GeV]");
    std::ostringstream yTitle;
    yTitle << "E^{2.6} Flux [m^{-2} s^{-1} sr^{-1} (GeV)^{1.6}]";
    frame->GetYaxis()->SetTitle(yTitle.str().c_str());
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleSize(0.048);
    frame->GetYaxis()->SetTitleSize(0.048);
    frame->GetXaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetXaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetTitleOffset(1.25);

    // Disegna prima le linee e poi punti+barre statistiche per evitare
    // che le linee di anni successivi coprano i marker precedenti.
    for (const auto& item : grafici) item.second->Draw("LX SAME");
    for (const auto& item : grafici) item.second->Draw("PZ SAME");

    auto* legenda = new TLegend(0.59, 0.72, 0.94, 0.91);
    legenda->SetNColumns(2);
    legenda->SetBorderSize(0);
    legenda->SetFillStyle(0);
    legenda->SetTextFont(42);
    legenda->SetTextSize(0.033);
    for (const auto& item : grafici) {
        const std::string anno = std::to_string(item.first);
        legenda->AddEntry(item.second, anno.c_str(), "lp");
    }
    legenda->Draw();

    c->RedrawAxis();
    c->Update();
    c->SaveAs("PLOTS/fluxes_pHe_perYear_2016_2025.pdf");
    c->SaveAs("PLOTS/fluxes_pHe_perYear_2016_2025.png");
}
