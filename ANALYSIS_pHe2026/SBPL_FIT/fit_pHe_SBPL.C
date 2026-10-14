// ROOT macro: PL vs SBPL fit of DAMPE p+He data with piecewise nuisance parameters.
// Run e.g.:
//   root -l -q 'fit_pHe_SBPL.C("pHe_flux.dat",2,2.0e4,1.0e6,5.0,true)'
//
// Convention:
//   Phi_SBPL = Phi0 (E/1 TeV)^(-gamma) [1 + (E/Eb)^s]^(DeltaGamma/s)
// so DeltaGamma > 0 corresponds to a HARDENING 

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TH1.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TString.h>

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/DistFunc.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace PHeFit {

struct Bin {
    double Emin = 0.0;
    double Emax = 0.0;
    double E    = 0.0;   // central/mean energy used in the fit
    double flux = 0.0;
    double stat = 0.0;   // absolute statistical error

    // The following systematic columns are RELATIVE/FRACTIONAL in the input file.
    double sysAna = 0.0;
    double sysHad = 0.0;
    double sysTot = 0.0;
};

struct Context {
    std::vector<Bin> fitBins;
    std::vector<double> regionEdges;
    std::vector<int> regionOfBin;
    std::vector<double> sigmaNuis; // Gaussian prior widths of nuisance parameters

    int m = 2;
    double E0 = 1.0e3;       // 1 TeV in GeV
    double fitEmin = 2.0e4;  // GeV
    double fitEmax = 1.0e6;  // GeV
    double sFixed  = 5.0;
};

struct FitResult {
    bool ok = false;
    int status = -999;
    double chi2 = std::numeric_limits<double>::infinity();
    int ndf = -1;
    std::vector<double> p;
    std::vector<double> e;

    // Diagnostica della matrice di covarianza
    int covStatus = -1;
    double edm = 0.0;
    std::vector<double> rhoEb;
};

Context gCtx;
int gNLowStat = 4;

// -----------------------------------------------------------------------------
// Models
// -----------------------------------------------------------------------------
double PL(double E, double phi0, double gamma)
{
    return phi0 * std::pow(E / gCtx.E0, -gamma);
}

double SBPL(double E, double phi0, double gamma,
            double Eb, double dGamma, double s)
{
    return phi0 * std::pow(E / gCtx.E0, -gamma)
         * std::pow(1.0 + std::pow(E / Eb, s), dGamma / s);
}

// -----------------------------------------------------------------------------
// Chi2 functions.
// PL parameters:   [Phi0, gamma, w0, ..., w_(m-1)]
// SBPL parameters: [Phi0, gamma, Eb, DeltaGamma, w0, ..., w_(m-1)]
// -----------------------------------------------------------------------------
double Chi2PL(const double *p)
{
    double chi2 = 0.0;

    for (size_t i = 0; i < gCtx.fitBins.size(); ++i) {
        const auto &b = gCtx.fitBins[i];
        const int j = gCtx.regionOfBin[i];
        const double w = p[2 + j];
        const double pred = PL(b.E, p[0], p[1]) * w;
        const double pull = (pred - b.flux) / b.stat;
        chi2 += pull * pull;
    }

    for (int j = 0; j < gCtx.m; ++j) {
        const double w = p[2 + j];
        const double pullNuis = (1.0 - w) / gCtx.sigmaNuis[j];
        chi2 += pullNuis * pullNuis;
    }

    return chi2;
}

double Chi2SBPL(const double *p)
{
    double chi2 = 0.0;

    for (size_t i = 0; i < gCtx.fitBins.size(); ++i) {
        const auto &b = gCtx.fitBins[i];
        const int j = gCtx.regionOfBin[i];
        const double w = p[4 + j];
        const double pred = SBPL(b.E, p[0], p[1], p[2], p[3], gCtx.sFixed) * w;
        const double pull = (pred - b.flux) / b.stat;
        chi2 += pull * pull;
    }

    for (int j = 0; j < gCtx.m; ++j) {
        const double w = p[4 + j];
        const double pullNuis = (1.0 - w) / gCtx.sigmaNuis[j];
        chi2 += pullNuis * pullNuis;
    }

    return chi2;
}

double Chi2PLowStat(const double *p)
{
    double chi2 = 0.0;

    for (int i=0; i < gNLowStat; ++i) {
        const auto &b = gCtx.fitBins[i];
        const double pred = PL(b.E, p[0], p[1]);
        const double pull = (pred - b.flux)/b.stat;
        chi2 += pull * pull;
    }
    return chi2;
}

// -----------------------------------------------------------------------------
// Input
// -----------------------------------------------------------------------------
std::vector<Bin> LoadData(const char *filename)
{
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "ERROR: cannot open " << filename << std::endl;
        return {};
    }

    std::vector<Bin> bins;
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream ss(line);
        Bin b;
        double sysHET, sysSTK, sysCharge, sysUnfold, sysMix, sysSat;

        if (!(ss >> b.Emin >> b.Emax >> b.E >> b.flux >> b.stat
                 >> sysHET >> sysSTK >> sysCharge >> sysUnfold >> sysMix >> sysSat 
                 >> b.sysAna >> b.sysHad >> b.sysTot)) {
            continue;
        }

        bins.push_back(b);
    }

    std::cout << "Loaded " << bins.size() << " bins from " << filename << std::endl;
    return bins;
}

int RegionIndex(double E)
{
    // Edges: [e0,e1), [e1,e2), ... last region includes upper boundary.
    for (int j = 0; j < gCtx.m - 1; ++j) {
        if (E < gCtx.regionEdges[j + 1]) return j;
    }
    return gCtx.m - 1;
}

bool BuildContext(const std::vector<Bin> &allBins,
                  int m, double fitEmin, double fitEmax, double sFixed)
{
    gCtx = Context{};
    gCtx.m = m;
    gCtx.fitEmin = fitEmin;
    gCtx.fitEmax = fitEmax;
    gCtx.sFixed = sFixed;

    if (m < 1 || fitEmin <= 0.0 || fitEmax <= fitEmin) {
        std::cerr << "ERROR: invalid fit configuration." << std::endl;
        return false;
    }

    // DAMPE-style logarithmic division of the fit range into m nuisance regions.
    gCtx.regionEdges.resize(m + 1);
    const double l0 = std::log(fitEmin);
    const double l1 = std::log(fitEmax);
    for (int j = 0; j <= m; ++j) {
        gCtx.regionEdges[j] = std::exp(l0 + (l1 - l0) * double(j) / double(m));
    }
    // TEST: confine tra i due nuisance a 250 TeV
    //if (m == 2) {
        //gCtx.regionEdges[1] = 2.5e5;
    //    gCtx.regionEdges[1] = 1.0e5;
    //}

    for (const auto &b : allBins) {
        // We apply the fit-range cut to the central energy, as in the quoted method.
        if (b.E >= fitEmin && b.E <= fitEmax) {
            gCtx.fitBins.push_back(b);
            gCtx.regionOfBin.push_back(RegionIndex(b.E));
        }
    }

    if (gCtx.fitBins.empty()) {
        std::cerr << "ERROR: no bins in fit range." << std::endl;
        return false;
    }

    // Since the supplied Sys_ana and Sys_had are already relative fractions, the
    // per-bin relative total is simply sqrt(Sys_ana^2 + Sys_had^2).
    gCtx.sigmaNuis.assign(m, 0.0);
    std::vector<int> counts(m, 0);

    for (size_t i = 0; i < gCtx.fitBins.size(); ++i) {
        const auto &b = gCtx.fitBins[i];
        const int j = gCtx.regionOfBin[i];
        const double relSys = std::hypot(b.sysAna, b.sysHad);
        gCtx.sigmaNuis[j] += relSys;
        counts[j]++;
    }

    for (int j = 0; j < m; ++j) {
        if (counts[j] == 0) {
            std::cerr << "ERROR: nuisance region " << j
                      << " contains zero fit bins. Reduce m or change fit range."
                      << std::endl;
            return false;
        }
        gCtx.sigmaNuis[j] /= counts[j];
    }

    std::cout << "\nFit range = [" << fitEmin / 1e3 << ", "
              << fitEmax / 1e3 << "] TeV (cut on E_mean)" << std::endl;
    std::cout << "N fit bins = " << gCtx.fitBins.size()
              << ", nuisance parameters m = " << m
              << ", s = " << sFixed << " (fixed)\n" << std::endl;

    for (int j = 0; j < m; ++j) {
        std::cout << " nuisance " << j
                  << ": E in [" << gCtx.regionEdges[j] / 1e3
                  << ", " << gCtx.regionEdges[j + 1] / 1e3 << "] TeV"
                  << ", sigma_w = " << gCtx.sigmaNuis[j]
                  << ", Nbins = " << counts[j] << std::endl;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Minimization
// -----------------------------------------------------------------------------
std::unique_ptr<ROOT::Math::Minimizer> MakeMinimizer()
{
    std::unique_ptr<ROOT::Math::Minimizer> min(ROOT::Math::Factory::CreateMinimizer("Minuit","Migrad"));

    if (!min) {
        std::cerr << "ERROR: impossibile creare Minuit/Migrad" << std::endl;
        return nullptr;
    }

    min->SetMaxFunctionCalls(1000000);
    min->SetMaxIterations(100000);
    min->SetTolerance(1e-7);
    min->SetPrintLevel(0);
    min->SetErrorDef(1.0);
    min->SetStrategy(2);
    return min;
}

FitResult FitPL()
{
    const int npar = 2 + gCtx.m;
    ROOT::Math::Functor f(&Chi2PL, npar);
    auto min = MakeMinimizer();
    min->SetFunction(f);

    const double gamma0 = 2.8;
    const auto &b0 = gCtx.fitBins.front();
    const double phi00 = b0.flux * std::pow(b0.E / gCtx.E0, gamma0);

    min->SetLimitedVariable(0, "Phi0", phi00, std::max(phi00 * 1e-3, 1e-10), 1e-8, 1e-1);
    min->SetLimitedVariable(1, "gamma", gamma0, 1e-3, 1.5, 4.0);

    //for (int j = 0; j < gCtx.m; ++j)
    //    min->SetLimitedVariable(2 + j, Form("w%d", j), 1.0, 1e-3, 0.3, 1.7);
    // only statistical errors:
    for (int j = 0; j < gCtx.m; ++j)
        min->SetFixedVariable(2 + j, Form("w%d", j), 1.0);

    const bool ok = min->Minimize();
    min->Hesse();

    FitResult r;
    r.ok = ok;
    r.status = min->Status();
    r.chi2 = min->MinValue();
    //r.ndf = int(gCtx.fitBins.size()) - npar;
    r.ndf = int(gCtx.fitBins.size()) - 2; // only stat
    r.p.assign(min->X(), min->X() + npar);
    r.e.assign(min->Errors(), min->Errors() + npar);
    return r;
}

FitResult RunSBPLOnce(const FitResult &pl, double ebSeed, double dgSeed)
{
    const int npar = 4 + gCtx.m;
    ROOT::Math::Functor f(&Chi2SBPL, npar);
    auto min = MakeMinimizer();
    min->SetFunction(f);

    const double phiSeed = pl.p.empty() ? 3e-4 : pl.p[0];
    const double gamSeed = pl.p.empty() ? 2.8  : pl.p[1];

    min->SetLimitedVariable(0, "Phi0", phiSeed, std::max(phiSeed * 1e-3, 1e-10), 1e-8, 1e-1);
    min->SetLimitedVariable(1, "gamma", gamSeed, 1e-3, 1.5, 4.0);
    min->SetLimitedVariable(2, "Eb", ebSeed, std::max(ebSeed * 1e-3, 1.0),
                            gCtx.fitEmin, gCtx.fitEmax);
    min->SetLimitedVariable(3, "DeltaGamma", dgSeed, 1e-3, -1.5, 1.5);

    //for (int j = 0; j < gCtx.m; ++j) {
    //    const double wSeed = (int(pl.p.size()) > 2 + j) ? pl.p[2 + j] : 1.0;
    //    min->SetLimitedVariable(4 + j, Form("w%d", j), wSeed, 1e-3, 0.3, 1.7);
    //}
    // only statistical errors:
    for (int j = 0; j < gCtx.m; ++j) {
        min->SetFixedVariable(4 + j, Form("w%d", j), 1.0);
    }

    const bool ok = min->Minimize();
    min->Hesse();

    FitResult r;
    r.ok = ok;
    r.status = min->Status();
    r.chi2 = min->MinValue();
    //r.ndf = int(gCtx.fitBins.size()) - npar;
    r.ndf = int(gCtx.fitBins.size()) - 4; // only stat
    r.p.assign(min->X(), min->X() + npar);
    r.e.assign(min->Errors(), min->Errors() + npar);

    // Stato della matrice di covarianza
    r.covStatus = min->CovMatrixStatus();
    r.edm = min->Edm();

    // Correlazioni tra Eb e gli altri parametri
    r.rhoEb.resize(
        npar,
        std::numeric_limits<double>::quiet_NaN()
    );

    const double varEb = min->CovMatrix(2, 2);

    if (varEb > 0.0) {
        for (int k = 0; k < npar; ++k) {
            const double varK = min->CovMatrix(k, k);

            if (varK > 0.0) {
                r.rhoEb[k] =
                    min->CovMatrix(2, k) /
                    std::sqrt(varEb * varK);
            }
        }
    }

    return r;
}

FitResult FitSBPL(const FitResult &pl)
{
    FitResult best;

    // Multi-start in log(Eb) and DeltaGamma to reduce the chance of landing
    // in a local/edge minimum. Both positive and negative DeltaGamma seeds
    // are tried; the fitted sign tells us whether the preferred feature is
    // a hardening (>0) or softening (<0).
    const std::vector<double> logFractions = {0.20, 0.40, 0.60, 0.80};
    const std::vector<double> dgSeeds = {0.15, 0.35, -0.20};

    const double l0 = std::log(gCtx.fitEmin);
    const double l1 = std::log(gCtx.fitEmax);

    for (double f : logFractions) {
        const double ebSeed = std::exp(l0 + f * (l1 - l0));
        for (double dgSeed : dgSeeds) {
            FitResult r = RunSBPLOnce(pl, ebSeed, dgSeed);
            if (std::isfinite(r.chi2) && r.chi2 < best.chi2)
                best = r;
        }
    }

    return best;
}

void CheckEbMinos(const FitResult &best)
{
    const int npar = 4 + gCtx.m;

    if (!best.ok || int(best.p.size()) != npar) {
        std::cout << "Invalid SBPL fit for MINOS\n";
        return;
    }

    ROOT::Math::Functor f(&Chi2SBPL, npar);

    auto min = MakeMinimizer();
    if (!min) return;

    min->SetFunction(f);
    // Ripartiamo dalla soluzione SBPL migliore
    min->SetLimitedVariable(0, "Phi0", best.p[0], std::max(best.p[0]*1e-3, 1e-10), 1e-8, 1e-1);
    min->SetLimitedVariable(1, "gamma", best.p[1], 1e-3, 1.5, 4.0);
    min->SetLimitedVariable(2, "Eb", best.p[2], std::max(best.p[2]*1e-3, 1.0), gCtx.fitEmin, gCtx.fitEmax);
    min->SetLimitedVariable(3, "DeltaGamma", best.p[3], 1e-3, -1.5, 1.5);

    //for (int j = 0; j < gCtx.m; ++j) {
    //    min->SetLimitedVariable(4+j, Form("w%d", j), best.p[4+j], 1e-3, 0.3, 1.7);
    //}
    for (int j = 0; j < gCtx.m; ++j) {
        min->SetFixedVariable(4 + j, Form("w%d", j), 1.0);
    }

    // Rifacciamo il fit partendo dal minimo migliore
    const bool ok = min->Minimize();
    min->Hesse();

    const double chi2Refit = min->MinValue();

    std::cout << "\n========== EB MINOS ==========\n";
    std::cout << "Refit chi2 = " << chi2Refit << "\n";
    std::cout << "Original chi2 = " << best.chi2 << "\n";

    if (!ok ||
        !std::isfinite(chi2Refit) ||
        std::abs(chi2Refit - best.chi2) > 1e-4) {
        std::cout << "WARNING: refit does not reproduce " << "the original minimum.\n";
        return;
    }

    // Correlazione globale di Eb con gli altri parametri
    std::cout << "GlobalCC(Eb) = "
              << min->GlobalCC(2) << "\n";

    // Errori asimmetrici MINOS
    double errLow = 0.0;
    double errUp  = 0.0;

    const bool minosOK =
        min->GetMinosError(2, errLow, errUp);

    std::cout << "MINOS ok = " << minosOK << "\n";

    if (minosOK) {
        std::cout << "Eb MINOS = " << min->X()[2]/1e3 << " -" << -errLow/1e3 << " +" << errUp/1e3 << " TeV\n";
    } else {
        std::cout << "MINOS failed or reached a limit." << std::endl;
    }
}

void ProfileEb(const FitResult &best)
{
    const int npar = 4 + gCtx.m;
    ROOT::Math::Functor f(&Chi2SBPL, npar);

    const int nPoints = 81;

    const double Emin = gCtx.fitEmin * 1.001;
    const double Emax = gCtx.fitEmax * 0.999;

    auto *g = new TGraph();
    std::ofstream out("profile_Eb_SBPL_StatOnly.dat");

    out << "# Eb_TeV  DeltaChi2\n";

    int nGood = 0;

    for (int i = 0; i < nPoints; ++i) {

        double fraction = double(i) / (nPoints - 1);
        double Eb = std::exp(std::log(Emin) + fraction * (std::log(Emax) - std::log(Emin)));
        auto min = MakeMinimizer();
        if (!min) return;

        min->SetFunction(f);

        // Normalizzazione
        min->SetLimitedVariable(0, "Phi0", best.p[0], std::max(best.p[0]*1e-3, 1e-10), 1e-8, 1e-1);
        // Indice spettrale
        min->SetLimitedVariable(1, "gamma", best.p[1], 1e-3, 1.5, 4.0);
        // Eb fissato durante questo fit
        min->SetFixedVariable(2, "Eb", Eb);
        // DeltaGamma libero
        min->SetLimitedVariable(3, "DeltaGamma", best.p[3], 1e-3, -1.5, 1.5);
        // Nuisance liberi
        //for (int j = 0; j < gCtx.m; ++j) {
        //    min->SetLimitedVariable(4+j, Form("w%d", j), best.p[4+j], 1e-3, 0.3, 1.7);
        //}
        for (int j = 0; j < gCtx.m; ++j) {
            min->SetFixedVariable(4 + j, Form("w%d", j), 1.0);
        }

        const bool ok = min->Minimize();

        if (!ok || min->Status() != 0 ||
            !std::isfinite(min->MinValue())) {
            std::cout << "Profile failed at Eb = " << Eb/1e3 << " TeV\n";
            continue;
        }

        double dChi2 = min->MinValue() - best.chi2;

        if (dChi2 < -1e-4) {
            std::cout << "WARNING: better minimum at Eb = " << Eb/1e3 << " TeV\n";
        }

        if (dChi2 < 0 && dChi2 > -1e-4)
            dChi2 = 0.0;

        g->SetPoint(nGood++, Eb/1e3, dChi2);
        out << Eb/1e3 << " " << dChi2 << "\n";
    }

    out.close();
    std::cout << "\nProfile points: " << nGood << "/" << nPoints << "\n";

    if (nGood < 2) return;

    auto *c = new TCanvas("c_profile_Eb", "Eb profile", 850, 600);

    c->SetLogx();
    c->SetGridy();

    g->SetTitle(";E_{b} (TeV);#Delta#chi^{2}");
    g->SetLineColor(kBlue+1);
    g->SetLineWidth(3);
    g->Draw("AL");
    g->GetYaxis()->SetRangeUser(0, 6);

    // Delta chi2 = 1
    auto *l1 = new TLine(Emin/1e3, 1, Emax/1e3, 1);
    l1->SetLineColor(kRed);
    l1->SetLineStyle(2);
    l1->Draw("SAME");

    // Delta chi2 = 4
    auto *l2 = new TLine(Emin/1e3, 4, Emax/1e3, 4);
    l2->SetLineColor(kGray+2);
    l2->SetLineStyle(3);
    l2->Draw("SAME");

    c->Modified();
    c->Update();

    c->SaveAs("profile_Eb_SBPL_StatOnly.pdf");
    c->SaveAs("profile_Eb_SBPL_StatOnly.png");
}

void PrintFit(const char *name, const FitResult &r, bool sbpl)
{
    std::cout << "\n========== " << name << " ==========" << std::endl;
    std::cout << "Minuit ok = " << r.ok << ", status = " << r.status << std::endl;
    std::cout << std::setprecision(8);
    std::cout << "chi2 / ndf = " << r.chi2 << " / " << r.ndf << std::endl;
    std::cout << "Phi0  = " << r.p[0] << " +/- " << r.e[0] << std::endl;
    std::cout << "gamma = " << r.p[1] << " +/- " << r.e[1] << std::endl;

    if (sbpl) {
        std::cout << "Eb    = " << r.p[2] / 1e3 << " +/- " << r.e[2] / 1e3 << " TeV" << std::endl;
        std::cout << "DeltaGamma = " << r.p[3] << " +/- " << r.e[3] << std::endl;
        //std::cout << "gamma_high = gamma - DeltaGamma = " << r.p[1] - r.p[3] << std::endl;
        std::cout << "s = " << gCtx.sFixed << " (fixed)" << std::endl;
        for (int j = 0; j < gCtx.m; ++j)
            std::cout << "w" << j << " = " << r.p[4 + j] << " +/- " << r.e[4 + j] << std::endl;
        std::cout << "\n--- SBPL DIAGNOSTICS ---" << std::endl;
        std::cout << "CovMatrixStatus = " << r.covStatus << std::endl;
        std::cout << "EDM = " << r.edm << std::endl;
        std::cout << "rho(Eb, Phi0) = " << r.rhoEb[0] << std::endl;
        std::cout << "rho(Eb, gamma) = " << r.rhoEb[1] << std::endl;
        std::cout << "rho(Eb, DeltaGamma) = " << r.rhoEb[3] << std::endl;
        for (int j = 0; j < gCtx.m; ++j) {
            std::cout << "rho(Eb, w" << j << ") = "
                      << r.rhoEb[4 + j] << std::endl;
        }
    } else {
        for (int j = 0; j < gCtx.m; ++j)
            std::cout << "w" << j << " = " << r.p[2 + j] << " +/- " << r.e[2 + j] << std::endl;
    }
}

// -----------------------------------------------------------------------------
// Plot
// -----------------------------------------------------------------------------
void DrawPlot(const std::vector<Bin> &bins,
              const FitResult &pl, const FitResult &sbpl,
              bool drawPL, bool drawSBPL,
              const char *outTag = "")
{
    const double plotPower = 2.6;

    const int n = int(bins.size());
    auto *gData = new TGraphErrors(n);
    auto *gAna  = new TGraphAsymmErrors(n+2);
    auto *gTot  = new TGraphAsymmErrors(n+2);
    //auto *upper = new TGraph(n+2);
    //auto *lower = new TGraph(n+2);

    for (int i = 0; i < n; ++i) {
        const auto &b = bins[i];
        const double weight = std::pow(b.E, plotPower);
        const double y = b.flux * weight;
        const double eyStat = b.stat * weight;
        const double relTot = std::hypot(b.sysAna, b.sysHad);

        gData->SetPoint(i, b.E, y);
        gData->SetPointError(i, 0.0, eyStat);

        const int j = i + 1;

        gAna->SetPoint(j, b.E, y);
        gAna->SetPointError(j, b.E - b.Emin, b.Emax - b.E, y * b.sysAna, y * b.sysAna);

        gTot->SetPoint(j, b.E, y);
        gTot->SetPointError(j, b.E - b.Emin, b.Emax - b.E, y * relTot, y * relTot);

        //upper->SetPoint(i, b.E, y + y*relTot);
        //lower->SetPoint(i, b.E, y - y*relTot);
    }

    auto extendBand = [&](TGraphAsymmErrors *g) {
        double x, y;

        // Primo bin: estendi fino a Emin
        g->GetPoint(1, x, y);
        double eyLow  = g->GetErrorYlow(1);
        double eyHigh = g->GetErrorYhigh(1);

        g->SetPoint(0, bins.front().Emin, y);
        g->SetPointError(0, 0, 0, eyLow, eyHigh);

        // Ultimo bin: estendi fino a Emax
        g->GetPoint(n, x, y);
        eyLow  = g->GetErrorYlow(n);
        eyHigh = g->GetErrorYhigh(n);

        g->SetPoint(n+1, bins.back().Emax, y);
        g->SetPointError(n+1, 0, 0, eyLow, eyHigh);
    };

    extendBand(gTot);
    extendBand(gAna);

    auto *upper = new TGraph();
    auto *lower = new TGraph();

    for (int i = 0; i < gTot->GetN(); ++i) {
        double x, y;
        gTot->GetPoint(i, x, y);

        const double yup = y + gTot->GetErrorYhigh(i);
        const double ylow = y - gTot->GetErrorYlow(i);

        upper->SetPoint(i, x, yup);
        lower->SetPoint(i, x, ylow);
    }

    upper->SetLineColor(17);
    upper->SetLineWidth(1);

    lower->SetLineColor(17);
    lower->SetLineWidth(1);


    gTot->SetFillColorAlpha(18, 0.90);
    gTot->SetLineColor(17);
    gAna->SetFillColorAlpha(16, 0.95);
    gAna->SetLineColor(16);

    gData->SetMarkerStyle(20);
    gData->SetMarkerSize(1.3);
    gData->SetLineColor(kRed+1);
    gData->SetLineWidth(3);
    gData->SetMarkerColor(kRed+1);

    upper->SetLineColor(17);
    upper->SetLineWidth(1);
    lower->SetLineColor(17);
    lower->SetLineWidth(1);

    const int nCurve = 500;
    auto *gSBPL = new TGraph(nCurve);
    auto *gPL = new TGraph(nCurve);

    const double l0 = std::log(gCtx.fitEmin);
    const double l1 = std::log(gCtx.fitEmax);
    for (int i = 0; i < nCurve; ++i) {
        const double E = std::exp(l0 + (l1 - l0) * double(i) / double(nCurve - 1));
        const double weight = std::pow(E, plotPower);

        const double ySBPL = SBPL(E, sbpl.p[0], sbpl.p[1], sbpl.p[2], sbpl.p[3], gCtx.sFixed) * weight;
        const double yPL   = PL(E, pl.p[0], pl.p[1]) * weight;

        // Plot the smooth physical model, NOT the piecewise nuisance multiplier.
        gSBPL->SetPoint(i, E, ySBPL);
        gPL->SetPoint(i, E, yPL);
    }

    gSBPL->SetLineWidth(3);
    gSBPL->SetLineColor(kBlue + 1);
    gPL->SetLineWidth(2);
    gPL->SetLineStyle(2);
    gPL->SetLineColor(kRed + 1);

    //auto *c = new TCanvas("c_pHe_paperMethod", "p+He SBPL paper method", 950, 500);
    auto *c = new TCanvas(Form("c_pHe_paperMethod%s", outTag), "p+He SBPL paper method", 950, 480);
    c->SetTopMargin(0.03);
    c->SetRightMargin(0.04);
    c->SetBottomMargin(0.14);
    c->SetLeftMargin(0.12);
    c->SetTicks(1,1);
    c->SetLogx();

    TH1F *frame = nullptr;
    if (plotPower == 2.7)
        frame = c->DrawFrame(3e1, 12e3, 1.5e6, 41e3);
    else
        frame = c->DrawFrame(2.5e1, 5.5e3, 1.5e6, 19e3);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("Kinetic energy (GeV)");
    if (plotPower == 2.7)
        frame->GetYaxis()->SetTitle("Flux #times E^{2.7} (m^{-2} sr^{-1} s^{-1} GeV^{1.7})");
    else
        frame->GetYaxis()->SetTitle("Flux #times E^{2.6} (m^{-2} sr^{-1} s^{-1} GeV^{1.6})");

    frame->GetXaxis()->SetLabelSize(0.05);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetTitleOffset(1.20);
    frame->GetYaxis()->SetLabelSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();

    gTot->Draw("E3 SAME");
    gAna->Draw("E3 SAME");
    upper->Draw("L SAME");
    lower->Draw("L SAME");
    gData->Draw("P SAME");
    if (drawSBPL) gSBPL->Draw("L SAME");
    if (drawPL) gPL->Draw("L SAME");

    auto *leg = new TLegend(0.16, 0.67, 0.42, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.05);
    leg->AddEntry(gData, "DAMPE", "P");
    leg->AddEntry(gAna, "syst. uncert. (analysis only)", "f");
    leg->AddEntry(gTot, "syst. uncert. (analysis + had. model)", "f");
    if (drawSBPL) leg->AddEntry(gSBPL, "SBPL", "l");
    if (drawPL) leg->AddEntry(gPL, "PL", "l");
    leg->Draw();

    TLatex pHe_tex;
    pHe_tex.SetNDC();
    pHe_tex.SetTextFont(62);
    pHe_tex.SetTextSize(0.06);
    pHe_tex.SetTextAlign(31);
    pHe_tex.DrawLatex(0.88, 0.88, "p+He");

    // Optional diagnostic: nuisance-region boundaries.
    // Uncomment if you want to see where each w_j acts.
    /*
    for (int j = 1; j < gCtx.m; ++j) {
        auto *line = new TLine(gCtx.regionEdges[j], 5e3,
                               gCtx.regionEdges[j], 17.5e3);
        line->SetLineStyle(3);
        line->SetLineColor(kGray + 2);
        line->Draw("SAME");
    }
    */

    c->Modified();
    c->Update();

    //TString outname = Form("fit_pHe_SBPL_10TeV_EPOSLHC_%dnuisance%s", gCtx.m, drawPL ? "_wPL" : "");
    TString outname = Form("fit_pHe_SBPL_10TeV_EPOSLHC_%dnuisance%s%s_statOnly", gCtx.m, drawPL ? "_wPL" : "", outTag);
    c->SaveAs(outname + ".pdf");
    c->SaveAs(outname + ".png");
}


} // namespace PHeFit

// =============================================================================
// Main entry point
// =============================================================================
void fit_pHe_SBPL(const char *filename = "pHe_flux_fit_10TeV_EPOSLHC_wSatCorr.dat",
                  int N_NUIS = 2,
                  double FIT_EMIN = 2.5e4,
                  double FIT_EMAX = 1.0e6,
                  double S_FIXED = 5.0,
                  bool DRAW_PL_TOO = true)
{
    using namespace PHeFit;

    const auto bins = LoadData(filename);
    if (bins.empty()) return;

    if (!BuildContext(bins, N_NUIS, FIT_EMIN, FIT_EMAX, S_FIXED)) return;

    const FitResult pl = FitPL();
    const FitResult sbpl = FitSBPL(pl);

    //std::cout << "\n========== FIT BIN REGIONS ==========\n";
    //std::cout << "E (TeV)" << "\tRegion 158" << "\tRegion 250\n";

    std::cout << "\n======= NUISANCE REGIONS =======:\n";
    for (size_t i = 0; i < gCtx.fitBins.size(); ++i) {
        std::cout << gCtx.fitBins[i].E / 1e3
                  << " TeV -> region "
                  << gCtx.regionOfBin[i]
                  << std::endl;
    }

    const double edge158 = std::sqrt(2.5e4 * 1.0e6);
    const double edge250 = 2.5e5;

    for (const auto &b : gCtx.fitBins) {

        int r158 = (b.E < edge158) ? 0 : 1;
        int r250 = (b.E < edge250) ? 0 : 1;

        std::cout << b.E / 1e3 << "\t\t" << r158 << "\t\t" << r250;

        if (r158 != r250)
            std::cout << "  <-- CHANGED REGION";
        std::cout << std::endl;
    }

    PrintFit("PL", pl, false);
    PrintFit("SBPL", sbpl, true);

    CheckEbMinos(sbpl);
    ProfileEb(sbpl);

    // With s fixed, SBPL has two additional shape parameters relative to PL:
    // Eb and DeltaGamma. This is the DAMPE-style Delta-chi2 comparison with Delta dof = 2.
    const double deltaChi2 = pl.chi2 - sbpl.chi2;
    const int deltaDof = 2;

    double pvalue = 1.0;
    double Z = 0.0;
    if (deltaChi2 > 0.0) {
        pvalue = ROOT::Math::chisquared_cdf_c(deltaChi2, deltaDof);
        // One-sided Gaussian-equivalent significance: P(N(0,1) > Z) = pvalue.
        Z = ROOT::Math::normal_quantile_c(pvalue, 1.0);
    }

    std::cout << "\n========== MODEL COMPARISON ==========" << std::endl;
    std::cout << "Delta chi2 = chi2_PL - chi2_SBPL = " << deltaChi2 << std::endl;
    std::cout << "Delta dof  = " << deltaDof << " (s fixed)" << std::endl;
    std::cout << "p-value    = " << pvalue << std::endl;
    std::cout << "Gaussian-equivalent Z = " << Z << " sigma" << std::endl;

    if (sbpl.p.size() >= 4) {
        if (sbpl.p[3] > 0.0)
            std::cout << "\nPreferred SBPL feature has DeltaGamma > 0: HARDENING." << std::endl;
        else
            std::cout << "\nPreferred SBPL feature has DeltaGamma < 0: SOFTENING." << std::endl;
    }

    if (sbpl.ndf <= 1) {
        std::cout << "WARNING: SBPL has only " << sbpl.ndf
                  << " residual degree(s) of freedom. Treat this configuration as a robustness test."
                  << std::endl;
    }

    DrawPlot(bins, pl, sbpl, false, true, "");


}

