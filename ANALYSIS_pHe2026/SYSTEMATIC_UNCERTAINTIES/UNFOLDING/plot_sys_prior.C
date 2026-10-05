
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
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

struct FluxData {
    string suffix;
    string label;
    int color;

    vector<double> E;
    vector<double> flux;
    vector<double> err;
};

void plot_sys_prior()
{
    gStyle->SetOptStat(0);

    const double index = 2.6;
    const int removeFirst = 3;
    const int removeLast = 5;

    const string prefix =
        "TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_PSDprogr_3sLow_6sUp_STKch450_comb_vert0e7_10TeV_EPOSLHC_E2e7_";

    vector<FluxData> data = {
        {"priorE2e5_10iter_smooth.dat", "prior E^{ -2.5}", kBlue},
        {"priorE2e6_10iter_smooth.dat", "prior E^{ -2.6}", kGreen+2},
        {"priorE2e7_10iter_smooth.dat", "prior E^{ -2.7}", kBlack},
        {"priorE2e8_10iter_smooth.dat", "prior E^{ -2.8}", kOrange+1},
        {"priorE2e9_10iter_smooth.dat", "prior E^{ -2.9}", kRed}
    };

    // Riferimento: 10 iterazioni con smoothing
    const int refIndex = 2;

    // --------------------------------------------------
    // LETTURA FILE
    // --------------------------------------------------

    for (auto &d : data) {

        ifstream file(prefix + d.suffix);

        if (!file.is_open()) {
            cerr << "Impossibile aprire: "
                 << prefix + d.suffix << endl;
            return;
        }

        double E, F, err, dummy, low, high;

        while (file >> E >> F >> err >> dummy >> low >> high) {
            d.E.push_back(E);
            d.flux.push_back(F);
            d.err.push_back(err);
        }

        file.close();

        int n = d.E.size();

        if (n <= removeFirst + removeLast) {
            cerr << "Troppi punti esclusi per " << d.suffix << endl;
            return;
        }

        // Escludi ultimi 5 e primi 3 bin
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

    // Verifica corrispondenza dei bin energetici
    for (const auto &d : data) {
        if ((int)d.E.size() != n) {
            cerr << "Numero di bin differente: " << d.suffix << endl;
            return;
        }

        for (int i = 0; i < n; i++) {
            if (abs(d.E[i] - ref.E[i]) > 1e-6 * ref.E[i]) {
                cerr << "Bin energetici non coincidenti in " << d.suffix << endl;
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

        auto gr = new TGraphErrors(n, d.E.data(), y.data(), ex.data(), ey.data());

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

            // Differenza percentuale
            double ratio = 100.0 * (Fref - F) / Fref;

            // Propagazione degli errori (non correlati)
            double ratioErr = 100.0 * sqrt( pow(Eerr / Fref, 2) + pow(F * Eref / (Fref * Fref), 2) );

            rx.push_back(d.E[i]);
            ry.push_back(ratio);
            rey.push_back(ratioErr);
            rex.push_back(0.0);

            envelope[i] = max(envelope[i], abs(ratio));
            maxRatio = max(maxRatio, abs(ratio) + ratioErr);
        }

        auto grRatio = new TGraphErrors( rx.size(), rx.data(), ry.data(), rex.data(), rey.data() );

        grRatio->SetMarkerStyle(20 + j);
        grRatio->SetMarkerSize(1.15);
        grRatio->SetMarkerColor(d.color);
        grRatio->SetLineColor(d.color);
        grRatio->SetLineWidth(2);

        ratioGraphs.push_back(grRatio);
    }

    // --------------------------------------------------
    // CANVAS
    // --------------------------------------------------

    auto c = new TCanvas("c", "pHe unfolding", 1000, 950);

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

    auto frame1 = pad1->DrawFrame( xmin, minFlux * 0.83, xmax, maxFlux * 1.35 );

    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0);

    frame1->GetYaxis()->SetTitle( "Flux #times E^{2.6} [m^{-2} s^{-1} sr^{-1} GeV^{1.6}]" );
    frame1->GetYaxis()->SetTitleOffset(1.5);

    for (auto gr : fluxGraphs)
        gr->Draw("P SAME");

    auto legend = new TLegend(0.18, 0.63, 0.4, 0.87);

    legend->SetHeader("p + He", "C");
    legend->SetTextSize(0.034);
    //legend->SetBorderSize(0);

    for (size_t j = 0; j < data.size(); j++)
        legend->AddEntry(fluxGraphs[j], data[j].label.c_str(), "p");

    legend->Draw();

    // --------------------------------------------------
    // PAD INFERIORE: DIFFERENZE PERCENTUALI
    // --------------------------------------------------

    pad2->cd();
    pad2->SetLogx();
    pad2->SetGrid();

    //double ylim = max(5.0, maxRatio * 1.2);

    auto frame2 = pad2->DrawFrame(xmin, -7, xmax, 7);

    frame2->SetTitle("");

    frame2->GetXaxis()->SetTitle("Energy [GeV]");
    frame2->GetYaxis()->SetTitle("Difference (%)");

    frame2->GetXaxis()->SetTitleSize(0.10);
    frame2->GetYaxis()->SetTitleSize(0.09);

    frame2->GetXaxis()->SetLabelSize(0.09);
    frame2->GetYaxis()->SetLabelSize(0.08);

    frame2->GetYaxis()->SetTitleOffset(0.55);
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

    c->SaveAs("PLOTS/pHe_prior_comparison_10TeV_EPOSLHC_sameW_diffP.pdf");
    c->SaveAs("PLOTS/pHe_prior_comparison_10TeV_EPOSLHC_sameW_diffP.png");

    ofstream output("TXT_FILES/pHe_prior_envelope_10TeV_EPOSLHC_sameW_diffP.dat");

    output << "# E_GeV max_deviation_percent\n";

    for (int i = 0; i < n; i++)
        output << ref.E[i] << " " << envelope[i] << "\n";

    output.close();

    cout << "Confronto completato: " << n << " bin energetici utilizzati." << endl;
}
