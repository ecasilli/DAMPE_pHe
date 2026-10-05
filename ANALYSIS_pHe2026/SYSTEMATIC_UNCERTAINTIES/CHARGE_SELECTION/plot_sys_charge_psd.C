#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>
#include <TColor.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

struct FluxDataPSD {
    string filename;
    string label;
    int color;
    int marker;
    vector<double> E, flux, err;
};

bool ReadFluxPSD(FluxDataPSD &d, int removeFirst, int removeLast)
{
    ifstream file(d.filename);
    if (!file) {
        cerr << "ERROR: file non trovato: " << d.filename << endl;
        return false;
    }

    string line;
    int lineNumber = 0;
    while (getline(file, line)) {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r\n") == string::npos ||
            line.find_first_not_of(" \t\r\n") == line.find('#')) continue;

        istringstream ss(line);
        double e, f, ef;
        if (!(ss >> e >> f >> ef)) {
            cerr << "ERROR: riga " << lineNumber << " non leggibile in "
                 << d.filename << endl;
            return false;
        }
        if (!isfinite(e) || !isfinite(f) || !isfinite(ef) || e <= 0 || ef < 0) {
            cerr << "ERROR: valori non validi alla riga " << lineNumber
                 << " in " << d.filename << endl;
            return false;
        }
        d.E.push_back(e);
        d.flux.push_back(f);
        d.err.push_back(ef);
    }

    const int originalN = static_cast<int>(d.E.size());
    if (originalN <= removeFirst + removeLast) {
        cerr << "ERROR: troppo pochi bin in " << d.filename << endl;
        return false;
    }
    const int last = originalN - removeLast;
    auto trim = [removeFirst, last](vector<double> &v) {
        v = vector<double>(v.begin() + removeFirst, v.begin() + last);
    };
    trim(d.E);
    trim(d.flux);
    trim(d.err);
    cout << d.label << ": " << originalN << " bin letti, "
         << d.E.size() << " utilizzati" << endl;
    return true;
}

void plot_sys_charge_psd()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    const double index = 2.6;
    const int removeFirst = 3;
    const int removeLast = 5;
    const int refIndex = 0; // nominale

    // MODIFICA SOLTANTO QUESTI TRE PERCORSI con i nomi reali dei tuoi .dat.
    // I nomi tight/loose qui sotto sono un esempio, NON vengono generati
    // automaticamente con questi nomi dalle macro charge_data/charge_mc.
    const string prefix =
        "TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_";

    vector<FluxDataPSD> data = {
        {prefix + "nominal.dat", "Nominal", kBlack, 20, {}, {}, {}},
        {prefix + "psd_tight05.dat", "PSD tight (-5%)", kBlue+1, 21, {}, {}, {}},
        {prefix + "psd_loose05.dat", "PSD loose (+5%)", kRed+1, 22, {}, {}, {}}
    };

    for (auto &d : data)
        if (!ReadFluxPSD(d, removeFirst, removeLast)) return;

    const FluxDataPSD &ref = data[refIndex];
    const int n = static_cast<int>(ref.E.size());
    for (const auto &d : data) {
        if (static_cast<int>(d.E.size()) != n) {
            cerr << "ERROR: numero di bin diverso per " << d.label << endl;
            return;
        }
        for (int i = 0; i < n; ++i) {
            if (abs(d.E[i] - ref.E[i]) > 1e-6 * max(1., abs(ref.E[i]))) {
                cerr << "ERROR: bin energetici non coincidenti: " << d.label << " al bin " << i << endl;
                return;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        if (ref.flux[i] <= 0) {
            cerr << "ERROR: flusso nominale non positivo a E=" << ref.E[i] << endl;
            return;
        }
    }

    vector<TGraphErrors*> fluxGraphs;
    vector<TGraphErrors*> ratioGraphs(data.size(), nullptr);
    vector<double> envelope(n, 0.0), deltaTight(n, 0.0), deltaLoose(n, 0.0);
    vector<string> maxSource(n, "");

    double yMin = numeric_limits<double>::max();
    double yMax = -numeric_limits<double>::max();
    double maxDeviation = 0.0;

    for (size_t j = 0; j < data.size(); ++j) {
        const auto &d = data[j];
        vector<double> y(n), ey(n), ex(n, 0.0);
        for (int i = 0; i < n; ++i) {
            const double factor = pow(d.E[i], index);
            y[i]  = d.flux[i] * factor;
            ey[i] = d.err[i] * factor;
            yMin = min(yMin, y[i] - ey[i]);
            yMax = max(yMax, y[i] + ey[i]);
        }
        auto gr = new TGraphErrors(n, d.E.data(), y.data(), ex.data(), ey.data());
        gr->SetLineColor(d.color);
        gr->SetMarkerColor(d.color);
        gr->SetMarkerStyle(d.marker);
        gr->SetMarkerSize(1.1);
        fluxGraphs.push_back(gr);

        if (static_cast<int>(j) == refIndex) continue;

        vector<double> ratio(n), ratioErr(n);
        for (int i = 0; i < n; ++i) {
            const double Fref = ref.flux[i];
            const double F = d.flux[i];
            ratio[i] = 100.0 * (Fref - F) / Fref;

            // Formula solo indicativa: trascura la covarianza con il nominale.
            ratioErr[i] = 100.0 * sqrt(pow(d.err[i]/Fref, 2) + pow(F * ref.err[i] / (Fref * Fref), 2));

            if (j == 1) deltaTight[i] = ratio[i];
            if (j == 2) deltaLoose[i] = ratio[i];
            if (abs(ratio[i]) > envelope[i]) {
                envelope[i] = abs(ratio[i]);
                maxSource[i] = d.label;
            }
            maxDeviation = max(maxDeviation, abs(ratio[i]));
        }
        auto grRatio = new TGraphErrors(n, d.E.data(), ratio.data(), ex.data(), ratioErr.data());
        grRatio->SetLineColor(d.color);
        grRatio->SetMarkerColor(d.color);
        grRatio->SetMarkerStyle(d.marker);
        grRatio->SetMarkerSize(1.1);
        grRatio->SetLineWidth(2);
        ratioGraphs[j] = grRatio;
    }

    auto c = new TCanvas("cChargePSD", "PSD charge selection", 1000, 950);
    auto pad1 = new TPad("padChargePSD1", "", 0, 0.3, 1, 1);
    auto pad2 = new TPad("padChargePSD2", "", 0, 0, 1, 0.3);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.04);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.23);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.04);
    pad1->Draw();
    pad2->Draw();

    const double xmin = ref.E.front()/1.2;
    const double xmax = ref.E.back()*1.2;

    pad1->cd();
    pad1->SetLogx();
    pad1->SetGrid();
    const double lower = (yMin > 0.0) ? yMin * 0.83 : yMin * 1.1;
    auto frame1 = pad1->DrawFrame(xmin, lower, xmax, yMax*1.35);
    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0);
    frame1->GetYaxis()->SetTitle("Flux #times E^{2.6} [m^{-2} s^{-1} sr^{-1} GeV^{1.6}]");
    frame1->GetYaxis()->SetTitleOffset(1.5);
    for (auto gr : fluxGraphs) gr->Draw("P SAME");

    auto legend = new TLegend(0.17, 0.69, 0.45, 0.88);
    legend->SetHeader("p + He", "C");
    legend->SetTextSize(0.034);
    for (size_t j = 0; j < data.size(); ++j)
        legend->AddEntry(fluxGraphs[j], data[j].label.c_str(), "p");
    legend->Draw();

    pad2->cd();
    pad2->SetLogx();
    pad2->SetGrid();
    // >= 7% per mantenere la scala del grafico delle iterazioni,
    // ma si espande da sola se le differenze sono maggiori.
    const double ylim = max(7.0, 1.25*maxDeviation);
    auto frame2 = pad2->DrawFrame(xmin, -ylim, xmax, ylim);
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

    auto zero = new TLine(xmin, 0, xmax, 0);
    zero->SetLineColor(kGray+2);
    zero->SetLineStyle(2);
    zero->Draw();
    for (auto gr : ratioGraphs)
        if (gr) gr->Draw("PL SAME");

    // Le due curve mostrano esattamente l'envelope esportato nel .dat.
    vector<double> negativeEnvelope(n);
    for (int i = 0; i < n; ++i) negativeEnvelope[i] = -envelope[i];
    auto grUpper = new TGraph(n, ref.E.data(), envelope.data());
    auto grLower = new TGraph(n, ref.E.data(), negativeEnvelope.data());
    for (auto gr : {grUpper, grLower}) {
        gr->SetLineColor(kMagenta+2);
        gr->SetLineStyle(2);
        gr->SetLineWidth(2);
        gr->Draw("L SAME");
    }

    gSystem->mkdir("PLOTS", true);
    gSystem->mkdir("TXT_FILES", true);
    c->SaveAs("PLOTS/pHe_PSD_charge_comparison_10TeV_EPOSLHC.pdf");
    c->SaveAs("PLOTS/pHe_PSD_charge_comparison_10TeV_EPOSLHC.png");

    ofstream envelopeFile("TXT_FILES/pHe_PSD_charge_envelope_10TeV_EPOSLHC.dat");
    ofstream differencesFile("TXT_FILES/pHe_PSD_charge_differences_10TeV_EPOSLHC.dat");
    if (!envelopeFile || !differencesFile) {
        cerr << "ERROR: impossibile scrivere gli output .dat" << endl;
        return;
    }
    envelopeFile << "# E_GeV max_deviation_percent\n";
    differencesFile << "# E_GeV delta_tight_percent delta_loose_percent envelope_percent\n";
    envelopeFile << setprecision(12);
    differencesFile << setprecision(12);
    for (int i = 0; i < n; ++i) {
        envelopeFile << ref.E[i] << " " << envelope[i] << "\n";
        differencesFile << ref.E[i] << " " << deltaTight[i] << " "
                        << deltaLoose[i] << " " << envelope[i] << "\n";
        cout << "E=" << ref.E[i] << " GeV: tight=" << deltaTight[i]
             << "%, loose=" << deltaLoose[i] << "%, envelope="
             << envelope[i] << "% (" << maxSource[i] << ")\n";
    }
    cout << "Confronto completato: " << n << " bin" << endl;
}

