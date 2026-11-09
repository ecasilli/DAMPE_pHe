
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TColor.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

struct FluxData {
    string path;
    string label;
    int color;

    vector<double> E;
    vector<double> flux;
    vector<double> err;
};

void plot_FLUKA_GEANT4_comparison()
{
    gStyle->SetOptStat(0);

    const double index = 2.6;
    const int removeFirst = 3;
    const int removeLast = 5;

    vector<FluxData> data = {
        {
            "../../../FLUKA_ANALYSIS/TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_FLUKA_XScorr_5PeV.dat",
            "FLUKA",
            kRed+1
        },
        {
            "../../TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC_XScorr.dat",
            "GEANT4",
            kBlue+1
        }
    };

    // GEANT4 come riferimento
    const int refIndex = 1;

    // --------------------------------------------------
    // LETTURA FILE
    // --------------------------------------------------

    for (auto &d : data) {

        ifstream file(d.path);

        if (!file.is_open()) {
            cerr << "Impossibile aprire: " << d.path << endl;
            return;
        }

        string line;
        double E, F, err, dummy, low, high;

        while (getline(file, line)) {

            istringstream iss(line);

            if (!(iss >> E >> F >> err >> dummy >> low >> high))
                continue;

            d.E.push_back(E);
            d.flux.push_back(F);
            d.err.push_back(err);
        }

        file.close();

        int n = d.E.size();

        if (n <= removeFirst + removeLast) {
            cerr << "Troppi punti esclusi per " << d.path << endl;
            return;
        }

        // Escludi primi 3 e ultimi 5 bin
        auto trim = [&](vector<double> &v) {
            v.erase(v.end() - removeLast, v.end());
            v.erase(v.begin(), v.begin() + removeFirst);
        };

        trim(d.E);
        trim(d.flux);
        trim(d.err);
    }

    const FluxData &ref = data[refIndex];
    const int n = ref.E.size();

    // --------------------------------------------------
    // VERIFICA BIN ENERGETICI
    // --------------------------------------------------

    for (const auto &d : data) {

        if ((int)d.E.size() != n) {
            cerr << "Numero di bin differente: " << d.path << endl;
            return;
        }

        for (int i = 0; i < n; i++) {

            if (abs(d.E[i] - ref.E[i]) > 1e-6 * abs(ref.E[i])) {
                cerr << "Bin energetici non coincidenti in " << d.path << " al bin " << i << endl;
                return;
            }
        }
    }

    // --------------------------------------------------
    // CREAZIONE GRAFICI
    // --------------------------------------------------

    vector<TGraphErrors*> fluxGraphs;
    vector<TGraphErrors*> ratioGraphs;
    vector<double> envelope(n, 0.0);

    double minFlux = 1e100;
    double maxFlux = 0.0;
    double maxRatio = 0.0;

    for (size_t j = 0; j < data.size(); j++) {

        const auto &d = data[j];

        vector<double> y(n), ey(n), ex(n, 0.0);

        for (int i = 0; i < n; i++) {

            double factor = pow(d.E[i], index);

            y[i] = d.flux[i] * factor;
            ey[i] = d.err[i] * factor;

            minFlux = min(minFlux, y[i] - ey[i]);
            maxFlux = max(maxFlux, y[i] + ey[i]);
        }

        auto gr = new TGraphErrors(
            n, d.E.data(), y.data(), ex.data(), ey.data()
        );

        gr->SetMarkerStyle(20 + j);
        gr->SetMarkerSize(1.15);
        gr->SetMarkerColor(d.color);
        gr->SetLineColor(d.color);

        fluxGraphs.push_back(gr);

        // Non disegnare il rapporto del riferimento
        if ((int)j == refIndex) {
            ratioGraphs.push_back(nullptr);
            continue;
        }

        vector<double> rx, ry, rey, rex;

        for (int i = 0; i < n; i++) {

            double Fref = ref.flux[i];
            double Eref = ref.err[i];

            double F = d.flux[i];
            double Eerr = d.err[i];

            if (Fref == 0.0) continue;

            // Differenza percentuale:
            // (FLUKA - GEANT4) / GEANT4 * 100
            double ratio = 100.0 * (F - Fref) / Fref;
            envelope[i] = max(envelope[i], abs(ratio));

            // Propagazione errori non correlati
            double ratioErr = 100.0 * sqrt(pow(Eerr / Fref, 2) + pow(F * Eref / (Fref * Fref), 2));

            rx.push_back(d.E[i]);
            ry.push_back(ratio);
            rey.push_back(ratioErr);
            rex.push_back(0.0);

            maxRatio = max(maxRatio, abs(ratio) + ratioErr);
        }

        auto grRatio = new TGraphErrors(rx.size(), rx.data(), ry.data(), rex.data(), rey.data());

        grRatio->SetMarkerStyle(20);
        grRatio->SetMarkerSize(1.15);
        grRatio->SetMarkerColor(kBlack);
        grRatio->SetLineColor(kBlack);
        grRatio->SetLineWidth(2);

        ratioGraphs.push_back(grRatio);
    }

    // --------------------------------------------------
    // CANVAS
    // --------------------------------------------------

    auto c = new TCanvas("c_FLUKA_GEANT4","pHe FLUKA vs GEANT4", 1000, 950);

    auto pad1 = new TPad("pad1", "", 0, 0.3, 1, 1);
    auto pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);

    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.04);

    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.23);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.04);

    pad1->Draw();
    pad2->Draw();

    double xmin = ref.E.front() / 1.2;
    double xmax = ref.E.back() * 1.2;

    // --------------------------------------------------
    // PAD SUPERIORE: FLUSSI
    // --------------------------------------------------

    pad1->cd();
    pad1->SetLogx();
    pad1->SetGrid();

    auto frame1 = pad1->DrawFrame(xmin, minFlux * 0.83, xmax, maxFlux * 1.35);

    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0);

    frame1->GetYaxis()->SetTitle("Flux #times E^{2.6} [m^{-2} s^{-1} sr^{-1} GeV^{1.6}]");

    frame1->GetYaxis()->SetTitleOffset(1.5);

    for (auto gr : fluxGraphs)
        gr->Draw("P SAME");

    auto legend = new TLegend(0.18, 0.68, 0.43, 0.87);

    legend->SetTextSize(0.034);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);

    for (size_t j = 0; j < data.size(); j++)
        legend->AddEntry(fluxGraphs[j], data[j].label.c_str(), "p");

    legend->Draw();

    // --------------------------------------------------
    // PAD INFERIORE: DIFFERENZA PERCENTUALE
    // --------------------------------------------------

    pad2->cd();
    pad2->SetLogx();
    pad2->SetGrid();

    // Range automatico, almeno +/- 7%
    double ylim = max(7.0, maxRatio * 1.2);
    //double ylim = 22.;

    auto frame2 = pad2->DrawFrame(xmin, -ylim, xmax, ylim);

    frame2->SetTitle("");

    frame2->GetXaxis()->SetTitle("Energy [GeV]");
    frame2->GetYaxis()->SetTitle("(FLUKA - G4) / G4 (%)");

    frame2->GetXaxis()->SetTitleSize(0.10);
    frame2->GetYaxis()->SetTitleSize(0.085);

    frame2->GetXaxis()->SetLabelSize(0.09);
    frame2->GetYaxis()->SetLabelSize(0.08);

    frame2->GetYaxis()->SetTitleOffset(0.65);
    frame2->GetXaxis()->SetTitleOffset(1.0);

    frame2->GetYaxis()->SetNdivisions(505);

    auto line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kGray+2);
    line->SetLineStyle(2);
    line->Draw();

    for (auto gr : ratioGraphs)
        if (gr) gr->Draw("PL SAME");

    // --------------------------------------------------
    // OUTPUT
    // --------------------------------------------------

    gSystem->mkdir("PLOTS", true);

    c->SaveAs("PLOTS/pHe_FLUKA_vs_GEANT4_XScorr_5PeV.pdf");
    c->SaveAs("PLOTS/pHe_FLUKA_vs_GEANT4_XScorr_5PeV.png");

    // --------------------------------------------------
    // SALVATAGGIO ENVELOPE FLUKA vs GEANT4
    // --------------------------------------------------

    gSystem->mkdir("TXT_FILES", true);

    ofstream output("TXT_FILES/pHe_FLUKA_GEANT4_envelope_XScorr_5PeV.dat");

    output << "# E_GeV max_deviation_percent\n";
    double maxDeviation = 0.0;
    double energyMaxDeviation = 0.0;

    for (int i = 0; i < n; i++) {

        output << ref.E[i] << " " << envelope[i] << "\n";

        if (envelope[i] > maxDeviation) {
            maxDeviation = envelope[i];
            energyMaxDeviation = ref.E[i];
        }
    }

    output.close();

    cout << "Envelope salvato in TXT_FILES/pHe_FLUKA_GEANT4_envelope_XScorr_5PeV.dat" << endl;

    cout << "Massima differenza percentuale FLUKA-GEANT4: " << maxDeviation << " %" << endl;

    cout << "Energia della massima differenza: " << energyMaxDeviation << " GeV" << endl;

    cout << "Confronto FLUKA vs GEANT4 completato: " << n << " bin energetici utilizzati." << endl;
}

