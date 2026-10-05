
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TAxis.h>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>

using namespace std;

struct FluxData {
    string suffix;
    string label;
    int color;

    vector<double> E;
    vector<double> flux;
    vector<double> err;
};

void plot_sys_weight_MC()
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // ------------------------------------------
    // CONFIGURAZIONE
    // ------------------------------------------

    const double index = 2.6;
    const int removeFirst = 3;
    const int removeLast = 5;
    const int refIndex = 2;  // E2e7 nominale

    const string prefix = "TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_PSDprogr_3sLow_6sUp_STKch450_comb_vert0e7_10TeV_EPOSLHC_";

    vector<FluxData> data = {
        {"E2e5_priorE2e5_10iter_smooth.dat", "E^{ -2e5}", kBlue},
        {"E2e6_priorE2e6_10iter_smooth.dat", "E^{ -2e6}", kGreen+2},
        {"E2e7_priorE2e7_10iter_smooth.dat", "E^{ -2e7} (ref.)", kBlack},
        {"E2e8_priorE2e8_10iter_smooth.dat", "E^{ -2e8}", kRed},
        {"E2e9_priorE2e9_10iter_smooth.dat", "E^{ -2e9}", kOrange+1}
    };

    // ------------------------------------------
    // LETTURA FILE .DAT
    // ------------------------------------------

    for (auto &d : data) {

        ifstream file(prefix + d.suffix);

        if (!file.is_open()) {
            cerr << "Errore apertura: " << prefix + d.suffix << endl;
            return;
        }

        string line;

        while (getline(file, line)) {

            istringstream iss(line);

            double E, F, err;

            // Legge solo le prime tre colonne
            // Ignora righe vuote o commenti
            if (!(iss >> E >> F >> err))
                continue;

            if (!isfinite(E) || !isfinite(F) || !isfinite(err) || E <= 0 || err < 0) {
                cerr << "Dati non validi in " << d.suffix << endl;
                return;
            }

            d.E.push_back(E);
            d.flux.push_back(F);
            d.err.push_back(err);
        }

        file.close();

        int n = d.E.size();

        if (n <= removeFirst + removeLast) {
            cerr << "Troppi punti esclusi in " << d.suffix << endl;
            return;
        }

        // Esclude primi 3 e ultimi 5 punti
        auto trim = [&](vector<double> &v) {
            v.erase(v.end()-removeLast, v.end());
            v.erase(v.begin(), v.begin()+removeFirst);
        };

        trim(d.E);
        trim(d.flux);
        trim(d.err);

        cout << d.label << ": " << d.E.size() << " bin utilizzati" << endl;
    }

    // ------------------------------------------
    // CONTROLLO BIN
    // ------------------------------------------

    const FluxData &ref = data[refIndex];
    const int n = ref.E.size();

    for (const auto &d : data) {

        if ((int)d.E.size() != n) {
            cerr << "Numero di bin differente: " << d.suffix << endl;
            return;
        }

        for (int i = 0; i < n; i++) {

            double tol = 1e-6 * max(1.0, abs(ref.E[i]));

            if (abs(d.E[i]-ref.E[i]) > tol) {
                cerr << "Bin energetici non coincidenti: " << d.suffix << endl;
                return;
            }
        }
    }

    // Il riferimento deve essere non nullo
    for (int i = 0; i < n; i++) {
        if (ref.flux[i] == 0.) {
            cerr << "Flusso nominale nullo a E = " << ref.E[i] << endl;
            return;
        }
    }

    // ------------------------------------------
    // CREAZIONE GRAFICI
    // ------------------------------------------

    vector<TGraphErrors*> fluxGraphs;
    vector<TGraphErrors*> ratioGraphs;

    vector<double> envelope(n, 0.0);

    double minFlux = numeric_limits<double>::max();
    double maxFlux = 0.;
    double maxRatio = 0.;

    for (size_t j = 0; j < data.size(); j++) {

        const auto &d = data[j];

        vector<double> y(n), ey(n), ex(n, 0.);

        for (int i = 0; i < n; i++) {

            double factor = pow(d.E[i], index);

            y[i]  = d.flux[i] * factor;
            ey[i] = d.err[i] * factor;

            minFlux = min(minFlux, y[i]-ey[i]);
            maxFlux = max(maxFlux, y[i]+ey[i]);
        }

        auto gr = new TGraphErrors(
            n, d.E.data(), y.data(),
            ex.data(), ey.data()
        );

        gr->SetMarkerStyle(20+j);
        gr->SetMarkerSize(1.15);
        gr->SetMarkerColor(d.color);
        gr->SetLineColor(d.color);

        fluxGraphs.push_back(gr);

        // Non disegna il rapporto del riferimento
        if ((int)j == refIndex) {
            ratioGraphs.push_back(nullptr);
            continue;
        }

        vector<double> rx(n), ry(n);
        vector<double> rex(n, 0.), rey(n);

        for (int i = 0; i < n; i++) {

            double Fref = ref.flux[i];
            double F = d.flux[i];

            double ErrRef = ref.err[i];
            double Err = d.err[i];

            double ratio = 100. * (Fref-F)/Fref;

            // Errori trattati come non correlati
            double ratioErr = 100. * sqrt( pow(Err/Fref, 2) + pow(F*ErrRef/(Fref*Fref), 2) );

            rx[i] = d.E[i];
            ry[i] = ratio;
            rey[i] = ratioErr;

            envelope[i] = max(envelope[i], abs(ratio));

            maxRatio = max(maxRatio, abs(ratio)+ratioErr);
        }

        auto grRatio = new TGraphErrors( n, rx.data(), ry.data(), rex.data(), rey.data() );

        grRatio->SetMarkerStyle(20+j);
        grRatio->SetMarkerSize(1.0);
        grRatio->SetMarkerColor(d.color);
        grRatio->SetLineColor(d.color);
        grRatio->SetLineWidth(2);

        ratioGraphs.push_back(grRatio);
    }


    cout << "\n========== CHECK ENVELOPE ==========\n";

    for (int i = 0; i < n; i++) {

        double maxFromGraph = 0.0;
        string maxLabel = "";

        cout << "\nE = " << ref.E[i] << " GeV\n";

        for (size_t j = 0; j < data.size(); j++) {

            if ((int)j == refIndex) continue;

            double x, y;

            ratioGraphs[j]->GetPoint(i, x, y);

            cout << "  " << data[j].label
                 << " : " << y << " %"
                 << "  (E graph = " << x << ")\n";

            if (abs(y) > maxFromGraph) {
                maxFromGraph = abs(y);
                maxLabel = data[j].label;
            }
        }

        cout << "  MAX FROM GRAPH = "
             << maxFromGraph << " % ("
             << maxLabel << ")\n";

        cout << "  ENVELOPE = "
             << envelope[i] << " %\n";

        if (abs(maxFromGraph - envelope[i]) > 1e-8) {
            cerr << "  ATTENZIONE: DISCREPANZA!\n";
        }
    }

    cout << "\n====================================\n";


    // ------------------------------------------
    // CANVAS
    // ------------------------------------------

    auto c = new TCanvas("cWeight", "MC weight comparison", 1000, 950);

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

    double xmin = ref.E.front()/1.2;
    double xmax = ref.E.back()*1.2;

    // ------------------------------------------
    // PAD SUPERIORE: FLUSSI
    // ------------------------------------------

    pad1->cd();
    pad1->SetLogx();
    pad1->SetGrid();

    auto frame1 = pad1->DrawFrame(xmin, 6.2e3, xmax,14.5e3);

    frame1->SetTitle("");
    frame1->GetXaxis()->SetLabelSize(0);

    frame1->GetYaxis()->SetTitle("Flux #times E^{2.6} [m^{-2} s^{-1} sr^{-1} GeV^{1.6}]");

    frame1->GetYaxis()->SetTitleOffset(1.5);

    for (auto gr : fluxGraphs)
        gr->Draw("P SAME");

    auto legend = new TLegend( 0.24, 0.52, 0.46, 0.86 );

    legend->SetHeader("p + He ", "C");
    legend->SetTextSize(0.035);
    //legend->SetBorderSize(0);

    for (size_t j = 0; j < data.size(); j++)
        legend->AddEntry(fluxGraphs[j], data[j].label.c_str(), "p");

    legend->Draw();

    // ------------------------------------------
    // PAD INFERIORE: DIFFERENZE %
    // ------------------------------------------

    pad2->cd();
    pad2->SetLogx();
    pad2->SetGrid();

    //double ylim = max(5.0, maxRatio*1.2);

    auto frame2 = pad2->DrawFrame(xmin, -4.5, xmax, 4.5);

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

    // ------------------------------------------
    // SALVATAGGIO
    // ------------------------------------------

    gSystem->mkdir("PLOTS", true);

    c->SaveAs("PLOTS/pHe_MC_weight_comparison_10TeV_EPOSLHC_diffP.png");
    c->SaveAs("PLOTS/pHe_MC_weight_comparison_10TeV_EPOSLHC_diffP.pdf");
    c->SaveAs("PLOTS/pHe_MC_weight_comparison_10TeV_EPOSLHC_diffP.eps");

    ofstream output(
        "TXT_FILES/pHe_MC_weight_envelope_10TeV_EPOSLHC_diffP.dat"
    );

    output << "# E_GeV max_deviation_percent\n";

    for (int i = 0; i < n; i++)
        output << ref.E[i] << " " << envelope[i] << "\n";

    output.close();

    cout << "Confronto completato!" << endl;
    cout << "Bin utilizzati: " << n << endl;
}
