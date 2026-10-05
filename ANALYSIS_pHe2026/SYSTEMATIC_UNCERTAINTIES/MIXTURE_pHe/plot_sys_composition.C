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
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

struct CompositionFluxData {
    string name;
    vector<double> E, flux, err, E_low, E_high;
};

// Formato atteso per ogni riga:
// E_GeV  flux  flux_stat_error  dummy  E_low  E_high
// Sono accettate righe commentate che iniziano con #.
bool ReadCompositionFlux(const string &filename, CompositionFluxData &d)
{
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "ERRORE: impossibile aprire " << filename << endl;
        return false;
    }

    string line;
    int lineNumber = 0;
    while (getline(in, line)) {
        ++lineNumber;
        const size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos || line[first] == '#') continue;

        istringstream ss(line);
        double e, f, ef, dummy, low, high;
        if (!(ss >> e >> f >> ef >> dummy >> low >> high)) {
            cerr << "ERRORE: formato non valido in " << filename
                 << " alla riga " << lineNumber << endl;
            return false;
        }
        if (!isfinite(e) || !isfinite(f) || !isfinite(ef) ||
            !isfinite(low) || !isfinite(high) ||
            e <= 0. || f < 0. || ef < 0. ||
            low <= 0. || high <= low) {
            cerr << "ERRORE: dati non validi in " << filename << " alla riga " << lineNumber << endl;
            return false;
        }

        d.E.push_back(e);
        d.flux.push_back(f);
        d.err.push_back(ef);
        d.E_low.push_back(low);
        d.E_high.push_back(high);
    }

    if (d.E.empty()) {
        cerr << "ERRORE: nessun dato in " << filename << endl;
        return false;
    }
    return true;
}

void TrimCompositionFlux(CompositionFluxData &d, int first, int last)
{
    auto trim = [first, last](vector<double> &v) {
        v.erase(v.end() - last, v.end());
        v.erase(v.begin(), v.begin() + first);
    };
    trim(d.E);
    trim(d.flux);
    trim(d.err);
    trim(d.E_low);
    trim(d.E_high);
}

bool SameEnergy(double a, double b)
{
    return fabs(a - b) <= 1e-5 * max(fabs(a), fabs(b));
}

void plot_sys_composition()
{
    gStyle->SetOptStat(0);

    // ============================================================
    // MODIFICA QUI I NOMI DEI TUOI DUE FILE DI FLUSSO (.dat).
    // Il nome _SBPLmix_ratioGen_ sotto e' indicativo: va sostituito con
    // quello realmente prodotto dal tuo script di unfolding.
    // ============================================================
    const string originalFile = "../../TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC.dat";
    const string fitFile = "TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6Up_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_SBPLmix_ratioGen.dat";

    const int removeFirst = 3;
    const int removeLast = 5;

    const double spectralIndex = 2.6;

    CompositionFluxData original, fit;
    original.name = "Originale";
    fit.name = "Composizione SBPL p+He";

    if (!ReadCompositionFlux(originalFile, original) ||
        !ReadCompositionFlux(fitFile, fit)) return;

    if (removeFirst < 0 || removeLast < 0 ||
        (int)original.E.size() <= removeFirst + removeLast ||
        (int)fit.E.size() <= removeFirst + removeLast) {
        cerr << "ERRORE: numero di bin insufficiente per il trimming." << endl;
        return;
    }

    TrimCompositionFlux(original, removeFirst, removeLast);
    TrimCompositionFlux(fit, removeFirst, removeLast);

    const int n = static_cast<int>(original.E.size());
    if ((int)fit.E.size() != n) {
        cerr << "ERRORE: numero di bin differente tra i due spettri." << endl;
        return;
    }

    for (int i = 0; i < n; ++i) {
        if (!SameEnergy(original.E[i], fit.E[i]) ||
            !SameEnergy(original.E_low[i], fit.E_low[i]) ||
            !SameEnergy(original.E_high[i], fit.E_high[i])) {
            cerr << "ERRORE: energie o estremi di bin differenti al bin " << i << endl;
            return;
        }
    }

    vector<double> yOriginal(n), eyOriginal(n), yFit(n), eyFit(n);
    vector<double> delta(n), absDelta(n), ex(n, 0.0), eyZero(n, 0.0);

    double minY = numeric_limits<double>::max();
    double maxY = -numeric_limits<double>::max();
    double maxAbsDelta = 0.0;

    for (int i = 0; i < n; ++i) {
        const double eFactor = pow(original.E[i], spectralIndex);
        yOriginal[i]  = original.flux[i] * eFactor;
        eyOriginal[i] = original.err[i] * eFactor;
        yFit[i]       = fit.flux[i] * eFactor;
        eyFit[i]      = fit.err[i] * eFactor;

        minY = min(minY, min(yOriginal[i] - eyOriginal[i], yFit[i] - eyFit[i]));
        maxY = max(maxY, max(yOriginal[i] + eyOriginal[i], yFit[i] + eyFit[i]));

        // Differenza rispetto all'originale: positiva se il fit aumenta il flusso.
        delta[i] = 100.0 * (fit.flux[i] - original.flux[i]) / original.flux[i];
        absDelta[i] = fabs(delta[i]);
        maxAbsDelta = max(maxAbsDelta, absDelta[i]);
    }

    auto grOriginal = new TGraphErrors(n, original.E.data(), yOriginal.data(), ex.data(), eyOriginal.data());
    auto grFit = new TGraphErrors(n, fit.E.data(), yFit.data(), ex.data(), eyFit.data());
    // Niente errori statistici nel rapporto: original e fit sono correlati.
    auto grDelta = new TGraphErrors(n, original.E.data(), delta.data(), ex.data(), eyZero.data());

    grOriginal->SetMarkerStyle(20);
    grOriginal->SetMarkerSize(1.10);
    grOriginal->SetMarkerColor(kBlack);
    grOriginal->SetLineColor(kBlack);

    grFit->SetMarkerStyle(24);
    grFit->SetMarkerSize(1.15);
    grFit->SetMarkerColor(kRed + 1);
    grFit->SetLineColor(kRed + 1);

    grDelta->SetMarkerStyle(20);
    grDelta->SetMarkerSize(1.05);
    grDelta->SetMarkerColor(kRed + 1);
    grDelta->SetLineColor(kRed + 1);
    grDelta->SetLineWidth(2);

    auto c = new TCanvas("c_composition", "p+He composition systematic", 1000, 950);
    auto pad1 = new TPad("pad_composition_top", "", 0, 0.30, 1, 1);
    auto pad2 = new TPad("pad_composition_bottom", "", 0, 0, 1, 0.30);

    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.04);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.23);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.04);
    pad1->Draw();
    pad2->Draw();

    const double xmin = original.E.front() / 1.2;
    const double xmax = original.E.back() * 1.2;

    // Pannello superiore: flussi pesati con E^2.6.
    pad1->cd();
    pad1->SetLogx();
    pad1->SetGrid();

    const double spanY = maxY - minY;
    const double paddingY = spanY > 0.0 ? 0.15 * spanY : 0.15 * fabs(maxY);
    //auto frame1 = pad1->DrawFrame(xmin, minY - paddingY, xmax, maxY + paddingY);
    auto frame1 = pad1->DrawFrame(xmin, 5.8e3, xmax, 13.6e3);
    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0);
    frame1->GetYaxis()->SetTitle("Flux #times E^{2.6} [m^{-2} s^{-1} sr^{-1} GeV^{1.6}]");
    frame1->GetYaxis()->SetTitleOffset(1.5);

    grOriginal->Draw("P SAME");
    grFit->Draw("P SAME");

    auto legend = new TLegend(0.17, 0.70, 0.52, 0.88);
    //legend->SetHeader("p + He", "C");
    legend->SetTextSize(0.033);
    //legend->SetBorderSize(0);
    legend->AddEntry(grOriginal, "p+He flux 50-50", "p");
    legend->AddEntry(grFit, "DAMPE composition model", "p");
    legend->Draw();

    // Pannello inferiore: differenza percentuale con segno.
    pad2->cd();
    pad2->SetLogx();
    pad2->SetGrid();
    const double ylim = max(2.0, 1.25 * maxAbsDelta);
    //auto frame2 = pad2->DrawFrame(xmin, -ylim, xmax, ylim);
    auto frame2 = pad2->DrawFrame(xmin, -5.4, xmax, 5.4);
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

    auto lineZero = new TLine(xmin, 0.0, xmax, 0.0);
    lineZero->SetLineColor(kGray + 2);
    lineZero->SetLineStyle(2);
    lineZero->Draw();
    grDelta->Draw("PL SAME");

    gSystem->mkdir("PLOTS", true);
    gSystem->mkdir("TXT_FILES", true);
    c->SaveAs("PLOTS/pHe_composition_SBPLmix_ratioGen_comparison_10TeV_EPOSLHC.pdf");
    c->SaveAs("PLOTS/pHe_composition_SBPLmix_ratioGen_comparison_10TeV_EPOSLHC.png");

    ofstream out("TXT_FILES/pHe_composition_SBPLmix_ratioGen_difference_10TeV_EPOSLHC.dat");
    if (!out.is_open()) {
        cerr << "ERRORE: non riesco a creare il file delle differenze." << endl;
        return;
    }
    out << "# E_GeV  delta_percent_signed  abs_delta_percent  E_low_GeV  E_high_GeV\n";
    out << "# delta = 100 * (flux_SBPLmix_ratioGen - flux_original) / flux_original\n";
    out << setprecision(12);
    for (int i = 0; i < n; ++i)
        out << original.E[i] << " " << delta[i] << " " << absDelta[i]
            << " " << original.E_low[i] << " " << original.E_high[i] << "\n";
    out.close();

    cout << "Confronto completato: " << n << " bin utilizzati." << endl;
    cout << "Massima differenza assoluta: " << maxAbsDelta << " %" << endl;
    cout << "Output: PLOTS/pHe_composition_SBPLmix_ratioGen_comparison_10TeV_EPOSLHC.{pdf,png}" << endl;
    cout << "        TXT_FILES/pHe_composition_SBPLmix_ratioGen_difference_10TeV_EPOSLHC.dat" << endl;
}

