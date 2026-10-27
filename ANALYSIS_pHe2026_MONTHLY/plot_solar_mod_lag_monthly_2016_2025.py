#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analisi p+He: adattamento del codice originale di Irene (all-particle).
Python 2.7.14 + PyROOT + numpy + scipy. NON usa matplotlib/pandas/pyunfold.

Mantiene i canvas originali: c, cy, cym, cdy0-2, cpippo*, ccorr0-2,
ctime, cpol0-2, c3d, c_asb_relvar, c_abs_relvar3.
Mantiene i due fit sinusoidali, le correlazioni e la scansione temporale.

Modifiche rispetto a Irene:
  - nome dei ROOT p+He e directory corretta;
  - istogramma 'unfold_cts' al posto di 'unf_cts';
  - esattamente 3 bin: 4, 5, 6 (etichette 50, 80, 125 GeV);
  - 2016-2025 (asse fino alla tacca 2026) e ctime con legenda a sinistra;
  - rimosse f-string e chiamate incompatibili con Python 2.7;
  - scale dei grafici e inizializzazione dei fit adeguate al flusso p+He.

I valori dei ROOT 'flux' sono GIA' flussi. NON si divide per la livetime.
La livetime non viene usata in questa versione, come nel codice di Irene.
Le medie e gli errori semestrali seguono ESATTAMENTE la sua ricetta:
  media = mean(6 mesi); err = mean(err mensili) + std(6 flussi)/sqrt(6).
"""
from __future__ import print_function

import os
import json
import math
import array as ary
import ROOT
import numpy as np
from scipy import stats

from ROOT import TCanvas, TFile, gStyle, gPad, TF1, TLatex

# ---------------------------------------------------------------------------
# CONFIGURAZIONE (percorsi assoluti: lo script gira da qualsiasi directory)
# ---------------------------------------------------------------------------
BASE_DIR = "./"
ROOT_DIR = os.path.join(BASE_DIR, "ROOT_FILES")
ROOT_TEMPLATE = ("unfold_results_pHe_%04d_%02d_3sigmaLow_6sigmaUp_"
                 "PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC.root")
SSN_NAME = "observed-solar-cycle-indices_untilJune2026.json"

FIRST_YEAR = 2016
Nyrs = 10                  # 2016, ..., 2025; tacca di destra = 2026
Nperiods = 2 * Nyrs        # 20 semestri
Nmonths = 12 * Nyrs        # 120 mesi
ENERGY_BINS = [4, 5, 6]  # ROOT conta i bin da 1
ENERGY_GEV = [50, 80, 125]
NENERGY = len(ENERGY_BINS)

# Salva ctime in PDF/PNG e, opzionalmente, tutti gli altri canvas originali.
OUTPUT_DIR = os.path.join(BASE_DIR, "PLOTS")
SAVE_CTIME = True
SAVE_ALL_ORIGINAL_CANVASES = False  # Irene non salvava automaticamente i canvas
BATCH = True                       # True funziona anche senza display X11

# Per ROOT e per la visualizzazione del flusso come nell'originale
alpha = 2.6
gamma = 2.6
Amc = 2.0 * math.pi * math.pi  # dato originale; NON applicato al flusso gia' pronto

# ---------------------------------------------------------------------------
# FUNZIONI DEL CODICE DI IRENE
# ---------------------------------------------------------------------------
def median_energy(e_min, e_max):
    return (((e_max**(1-alpha) - e_min**(1-alpha)) /
             ((1-alpha) * (e_max-e_min)))**(-1.0/alpha))


def sin_fit(x, par):
    time = x[0]
    arg = 2.0 * math.pi * (time - par[3]) / par[2]
    return par[0] + par[1] * math.sin(arg)


def double_sin_fit(x, par):
    time = x[0]
    arg1 = 2.0 * math.pi * (time - par[3]) / par[2]
    arg2 = 2.0 * math.pi * (time - par[6]) / par[5]
    return par[0] + par[1] * math.sin(arg1) + par[4] * math.sin(arg2)


def array_d(values):
    return ary.array('d', [float(v) for v in values])


def check_file(path):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise IOError("ROOT mancante/non leggibile: " + path)
    return f


def get_hist(f, name, path):
    h = f.Get(name)
    if not h:
        raise KeyError("Istogramma '%s' assente in %s" % (name, path))
    h.SetDirectory(0)
    if h.GetNbinsX() < max(ENERGY_BINS):
        raise ValueError("Istogramma '%s' con troppo pochi bin: %s" % (name, path))
    return h


def find_ssn_json():
    candidates = [os.path.join(BASE_DIR, SSN_NAME),
                  os.path.join(BASE_DIR, "TXT_FILES", SSN_NAME),
                  os.path.join(os.getcwd(), SSN_NAME)]
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise IOError("JSON SSN non trovato. Copialo in %s o TXT_FILES.\nPercorsi provati: %s" %
                  (BASE_DIR, ", ".join(candidates)))


def read_ssn_monthly(path):
    # Equivalente del read_json + indice mensile nel codice Irene.
    with open(path, 'r') as handle:
        rows = json.load(handle)
    if isinstance(rows, dict):
        # Formato opzionale contenente una lista nel campo 'data'.
        rows = rows.get('data', rows.get('results', []))
    ssn = {}
    for r in rows:
        date = str(r.get('time-tag', ''))[:7]
        if len(date) != 7 or date[4] != '-':
            continue
        year, month = int(date[:4]), int(date[5:7])
        try:
            val = float(r['ssn'])
        except (ValueError, TypeError, KeyError):
            continue
        if not (val == val) or val < 0.0:  # NOAA: -1 = missing
            continue
        ssn[(year, month)] = val
    return ssn


def ssn_for_halfyear(ssn, year, half):
    months = range(1, 7) if half == 1 else range(7, 13)
    vals = []
    for month in months:
        if (year, month) not in ssn:
            raise ValueError("Manca SSN %04d-%02d (semestre %iH%i)" %
                             (year, month, year, half))
        vals.append(ssn[(year, month)])
    return np.array(vals, dtype=float)


def try_ssn_for_halfyear(ssn, year, half):
    months = range(1, 7) if half == 1 else range(7, 13)
    if any((year, month) not in ssn for month in months):
        return None
    return np.array([ssn[(year, month)] for month in months], dtype=float)


def sem_from_index(index):
    index = int(index)
    return FIRST_YEAR + (index // 2), 1 + (index % 2)


def safe_corr(x, y):
    if len(x) < 3 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return 0.0, 1.0
    # Esattamente il coefficiente di ROOT: TGraph::GetCorrelationFactor()
    g = ROOT.TGraph(len(x), array_d(x), array_d(y))
    r = float(g.GetCorrelationFactor())
    r = max(-1.0, min(1.0, r))
    if abs(r) >= 1.0:
        p = 0.0
    else:
        t = abs(r) * math.sqrt((len(x) - 2.0) / (1.0 - r*r))
        p = 2.0 * (1.0 - stats.t.cdf(t, df=len(x)-2))
    return r, p


def y_range(values, fraction=0.10):
    low = min(values)
    high = max(values)
    delta = high - low
    if delta == 0.0:
        delta = max(abs(low) * 0.05, 1.e-12)
    return low - fraction*delta, high + fraction*delta


def save_canvas(canvas):
    if not os.path.isdir(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    canvas.Modified()
    canvas.Update()
    canvas.SaveAs(os.path.join(OUTPUT_DIR, canvas.GetName() + ".png"))
    canvas.SaveAs(os.path.join(OUTPUT_DIR, canvas.GetName() + ".pdf"))


# Conserva gli oggetti Python, affinche' PyROOT non li distrugga.
keep_objects = []
original_canvases = []


def keep(obj):
    keep_objects.append(obj)
    return obj


def new_canvas(name, title, width=900, height=700):
    canvas = keep(ROOT.TCanvas(name, title, width, height))
    original_canvases.append(canvas)
    return canvas


def draw_line(x, ymin, ymax, color, style=9, width=2):
    line = keep(ROOT.TLine(x, ymin, x, ymax))
    line.SetLineColor(color)
    line.SetLineStyle(style)
    line.SetLineWidth(width)
    line.Draw("same")
    return line


# ---------------------------------------------------------------------------
# ANALISI PRINCIPALE
# ---------------------------------------------------------------------------
def main():
    ROOT.gROOT.SetBatch(BATCH)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)

    period_times = np.array([year + (0.25 if half == 1 else 0.75)
                             for year in range(FIRST_YEAR, FIRST_YEAR + Nyrs)
                             for half in (1, 2)], dtype=float)
    period_labels = ["%iH%i" % (year, half)
                     for year in range(FIRST_YEAR, FIRST_YEAR + Nyrs)
                     for half in (1, 2)]
    print("Periodi:", period_labels)

    hother_unf_yrs_pow = []
    hother_unf_yrs = []
    hother_unf_monthly = []
    flux_PTstatistics = np.full(NENERGY, 1.e9)

    # Come Irene: 6 ROOT mensili -> 1 istogramma medio semestrale.
    for iyr in range(FIRST_YEAR, FIRST_YEAR + Nyrs):
        for half in (1, 2):
            months = range(1, 7) if half == 1 else range(7, 13)
            monthly_histograms = {"flux_pow": [], "flux": []}
            for month in months:
                path = os.path.join(ROOT_DIR, ROOT_TEMPLATE % (iyr, month))
                infile = check_file(path)
                for hist_name in ("flux_pow", "flux"):
                    hist = get_hist(infile, hist_name, path)
                    monthly_histograms[hist_name].append(hist)
                    if hist_name == "flux":
                        hother_unf_monthly.append(hist)
                hist_cts = get_hist(infile, "unfold_cts", path)
                for ipt, bin_number in enumerate(ENERGY_BINS):
                    v = hist_cts.GetBinContent(bin_number)
                    if v < flux_PTstatistics[ipt]:
                        flux_PTstatistics[ipt] = v
                infile.Close()
                print("Loaded", "%i %02i" % (iyr, month), path)

            label = "%iH%i" % (iyr, half)
            for hist_name, collection in (("flux_pow", hother_unf_yrs_pow),
                                          ("flux", hother_unf_yrs)):
                hists = monthly_histograms[hist_name]
                result = hists[0].Clone("%s_%s_mean" % (hist_name, label))
                result.SetDirectory(0)
                for h in hists[1:]:
                    if h.GetNbinsX() != result.GetNbinsX():
                        raise ValueError("Bin incoerenti: %s %s" % (hist_name, label))
                for ibin in range(result.GetNbinsX() + 2):
                    values = np.array([h.GetBinContent(ibin) for h in hists], dtype=float)
                    errors = np.array([h.GetBinError(ibin) for h in hists], dtype=float)
                    result.SetBinContent(ibin, float(np.mean(values)))
                    # IDENTICA prescrizione di Irene (non sostituita con quadratura).
                    result.SetBinError(ibin,
                        float(np.mean(errors) + np.std(values)/math.sqrt(6.0)))
                collection.append(result)
            print("Averaged six monthly histograms:", label)

    # SSN mensili e medie esatte dei 6 mesi (senza pandas/f-string).
    ssn_data = read_ssn_monthly(find_ssn_json())
    ssn_halfyear_months = np.array([
        ssn_for_halfyear(ssn_data, iyr, half)
        for iyr in range(FIRST_YEAR, FIRST_YEAR + Nyrs)
        for half in (1, 2)], dtype=float)
    xssn = array_d(np.mean(ssn_halfyear_months, axis=1))
    exssn = array_d(np.std(ssn_halfyear_months, axis=1) / math.sqrt(6.))
    ssn_yrs_array = ssn_halfyear_months
    dmp_ssn = [ssn_data[(year, month)]
               for year in range(FIRST_YEAR, FIRST_YEAR + Nyrs)
               for month in range(1, 13)]

    # -----------------------------------------------------------------------
    # c: spettri p+He dei 20 semestri (come Irene)
    # -----------------------------------------------------------------------
    c = new_canvas("c", "c")
    c.cd()
    c.SetGrid()
    c.SetLogx()
    c.SetLogy()
    for idx, label in enumerate(period_labels):
        h = hother_unf_yrs_pow[idx]
        h.SetMarkerStyle(8)
        h.SetMarkerColor(51 + (idx % 50))
        h.SetLineColor(51 + (idx % 50))
        h.SetTitle("p+He: %s" % label)
        h.Draw("PEHIST" if idx == 0 else "samePEHIST")
    c.BuildLegend()

    # -----------------------------------------------------------------------
    # cy, cym, cdy: stessi grafici di Irene, ora nei bin 4, 5, 6
    # -----------------------------------------------------------------------
    cy = new_canvas("cy", "cy")
    cy.SetGrid()
    hy = keep(ROOT.TH1D("hy", "hy; year (6-month bins); third pt",
                        Nperiods, FIRST_YEAR, FIRST_YEAR + Nyrs))
    hym = keep(ROOT.TH1D("hym", "Mean of 50-125 GeV flux; year (6-month bins); 3 pts mean",
                         Nperiods, FIRST_YEAR, FIRST_YEAR + Nyrs))
    y_pow = ary.array('d', [0.] * Nperiods)
    ey_pow = ary.array('d', [0.] * Nperiods)
    yy = np.zeros((NENERGY, Nperiods))
    eyy = np.zeros((NENERGY, Nperiods))
    dene_max = np.zeros(NENERGY)
    edene_max = np.zeros(NENERGY)
    energies_2d = []
    fluxes_2d = []
    years_2d = []
    mean_spotN = 0.0

    for idx, timestamp in enumerate(period_times):
        pt3_pow = hother_unf_yrs_pow[idx].GetBinContent(ENERGY_BINS[2])
        y_pow[idx] = pt3_pow
        ey_pow[idx] = hother_unf_yrs_pow[idx].GetBinError(ENERGY_BINS[2])
        mean_spotN += xssn[idx]
        mean = 0.0
        emean = 0.0
        for ipt, bin_index in enumerate(ENERGY_BINS):
            mean += hother_unf_yrs_pow[idx].GetBinContent(bin_index)
            emean += hother_unf_yrs_pow[idx].GetBinError(bin_index)
            val_flu = hother_unf_yrs[idx].GetBinContent(bin_index)
            err_flu = hother_unf_yrs[idx].GetBinError(bin_index)
            yy[ipt][idx] = val_flu
            eyy[ipt][idx] = err_flu
            ene_center = hother_unf_yrs[idx].GetBinCenter(bin_index)
            if ene_center <= 0.0:
                raise ValueError("Energia non positiva nel bin ROOT %d" % bin_index)
            energies_2d.append(math.log10(ene_center))
            fluxes_2d.append(val_flu)
            years_2d.append(timestamp)
            dene_max[ipt] = ene_center
            edene_max[ipt] = (ene_center -
                             hother_unf_yrs[idx].GetBinLowEdge(bin_index))
        mean /= float(NENERGY)
        emean /= float(NENERGY)
        hy.SetBinContent(idx+1, pt3_pow)
        hy.SetBinError(idx+1, ey_pow[idx])
        hym.SetBinContent(idx+1, mean)
        hym.SetBinError(idx+1, emean)

    mean_spotN /= float(Nperiods)
    mean_flux = np.mean(yy, axis=1)
    print("Mean p+He flux:", mean_flux)
    hym.SetLineColor(910)
    hym.SetLineWidth(6)
    hym.Draw()
    hy.SetLineWidth(4)
    hy.SetLineColor(394)
    hy.SetTitle("E = %d GeV" % ENERGY_GEV[2])
    hy.Draw("same")
    cy.BuildLegend()

    hym_spot = keep(ROOT.TGraphErrors(Nperiods, xssn,
                                     array_d(yy[2]), exssn, array_d(eyy[2])))
    cym = new_canvas("cym", "cym")
    cym.SetGrid()
    cym.SetLeftMargin(0.13)
    cym.SetRightMargin(0.03)
    cym.SetTopMargin(0.05)
    cym.SetBottomMargin(0.09)
    hym_spot.SetMarkerStyle(8)
    hym_spot.GetXaxis().SetTitle("Mean SunSpot Number over 6 months")
    hym_spot.GetYaxis().SetTitle(
        "#Phi_{p+He}(E = %i GeV) [GeV^{-1} (m^{2} s sr)^{-1}]" % ENERGY_GEV[2])
    hym_spot.GetXaxis().SetTitleSize(0.04)
    hym_spot.GetXaxis().SetTitleOffset(0.9)
    hym_spot.GetYaxis().SetTitleSize(0.04)
    hym_spot.GetYaxis().SetTitleOffset(1.3)
    hym_spot.Draw("AP")
    for ipt in range(Nperiods):
        latex = keep(ROOT.TLatex(hym_spot.GetX()[ipt],
                                 hym_spot.GetY()[ipt], period_labels[ipt]))
        latex.SetTextSize(0.03)
        latex.SetTextColor(633)
        latex.Draw("same")

    dx = ary.array('d', [0.] * Nperiods)
    edx = ary.array('d', [0.] * Nperiods)
    dy = np.zeros((NENERGY, Nperiods))
    edy = np.zeros((NENERGY, Nperiods))
    dy_max = np.zeros(NENERGY)
    edy_max = np.zeros(NENERGY)
    print("Mean spot number:", mean_spotN)
    for i in range(Nperiods):
        dx[i] = (xssn[i] - mean_spotN) / mean_spotN
        edx[i] = 2.0 * exssn[i] / mean_spotN
        for ipt in range(NENERGY):
            if mean_flux[ipt] == 0.0:
                raise ValueError("Flusso medio zero per %i GeV" % ENERGY_GEV[ipt])
            dy[ipt][i] = (yy[ipt][i] - mean_flux[ipt]) / mean_flux[ipt]
            edy[ipt][i] = 2.0 * eyy[ipt][i] / mean_flux[ipt]
            if dy_max[ipt] < abs(dy[ipt][i]):
                dy_max[ipt] = abs(dy[ipt][i])
                edy_max[ipt] = edy[ipt][i]
    print("Max absolute relative variations:")
    for ipt in range(NENERGY):
        print("at %d GeV: %.6g +- %.6g" %
              (ENERGY_GEV[ipt], dy_max[ipt], edy_max[ipt]))

    for iene in range(NENERGY):
        cdy = new_canvas("cdy%i" % iene, "cdy%i" % iene, 1100, 1100)
        cdy.SetGrid()
        cdy.SetLeftMargin(0.13)
        cdy.SetRightMargin(0.03)
        cdy.SetTopMargin(0.05)
        cdy.SetBottomMargin(0.09)
        hdym2 = keep(ROOT.TGraphErrors(Nperiods, dx, array_d(dy[iene]),
                                      edx, array_d(edy[iene])))
        hdym = keep(ROOT.TGraph(Nperiods, dx, array_d(dy[iene])))
        hdym.SetTitle("")
        hdym.GetXaxis().SetTitle("#Delta(SunSpot N.) / <SunSpot N.>")
        hdym.GetYaxis().SetTitle(
            "#Delta(#Phi_{%i GeV}) / <#Phi_{%i GeV}>" %
            (ENERGY_GEV[iene], ENERGY_GEV[iene]))
        hdym.GetXaxis().SetTitleSize(0.04)
        hdym.GetXaxis().SetTitleOffset(0.9)
        hdym.GetYaxis().SetTitleSize(0.04)
        hdym.GetYaxis().SetTitleOffset(1.5)
        ymin, ymax = y_range(dy[iene], 0.20)
        hdym.GetYaxis().SetRangeUser(ymin, ymax)
        hdym.GetXaxis().SetLimits(-2.5, 2.5)
        hdym.SetLineStyle(7)
        hdym.SetLineWidth(1)
        hdym.Draw("AL")
        corr_root, pval = safe_corr(dx, dy[iene])
        print("E=%i GeV: Correlation coefficient (ROOT)=%g; p-value=%g" %
              (ENERGY_GEV[iene], corr_root, pval))
        sigma_x = np.std(dx)
        slope = corr_root * np.std(dy[iene]) / sigma_x if sigma_x else 0.0
        intercept = np.mean(dy[iene]) - slope*np.mean(dx)
        corr_line = keep(ROOT.TF1("corr_line_%d" % iene,
                                  "[0] + [1]*x", min(dx), max(dx)))
        corr_line.SetParameters(intercept, slope)
        corr_line.SetLineColor(632)
        corr_line.SetLineWidth(2)
        corr_line.SetLineStyle(9)
        corr_line.Draw("same")
        hdym2.SetMarkerStyle(8)
        hdym2.SetMarkerColor(633)
        hdym2.SetMarkerSize(1.2)
        hdym2.SetLineColor(633)
        hdym2.SetLineWidth(2)
        for ipt in range(Nperiods):
            latex = keep(ROOT.TLatex(hdym2.GetX()[ipt],
                                     hdym2.GetY()[ipt], period_labels[ipt]))
            latex.SetTextSize(0.03)
            latex.SetTextColor(1)
            latex.Draw("same")
        hdym2.Draw("P SAME")
        legenda = keep(ROOT.TLegend(0.6, 0.7, 0.88, 0.88))
        legenda.AddEntry(hdym2,
                         "Rel. variations in SunSpot N. and #Phi_{%i GeV}" %
                         ENERGY_GEV[iene], "Pe")
        legenda.Draw()

    # -----------------------------------------------------------------------
    # Lag scan: l'SSN traslato di +1 semestre e' POSTERIORE al flusso
    # Identica logica senza ricircolo dei punti ai bordi temporali.
    # -----------------------------------------------------------------------
    shift_values = np.arange(-10, 11, dtype=int)
    shifts = shift_values.astype(float) * 0.5
    Nshifts = len(shift_values)
    corr_facs = np.full((NENERGY, Nshifts), np.nan)
    pval_corr_facs = np.full((NENERGY, Nshifts), np.nan)

    for indx_shift, shift in enumerate(shift_values):
        flux_indices = []
        shifted_ssn = []
        shifted_ssn_err = []
        for iflux in range(Nperiods):
            # Permette shift all'indietro e fino a giugno 2026, se presente.
            total_half = FIRST_YEAR * 2 + iflux + int(shift)
            shifted_year = total_half // 2
            shifted_half = 1 + (total_half % 2)
            months_ssn = try_ssn_for_halfyear(ssn_data, shifted_year, shifted_half)
            if months_ssn is None:
                continue
            flux_indices.append(iflux)
            shifted_ssn.append(float(np.mean(months_ssn)))
            shifted_ssn_err.append(float(np.std(months_ssn) / math.sqrt(6.)))
        if len(flux_indices) < 3:
            print("AVVISO: shift %+0.1f anni saltato (<3 coppie)" % shifts[indx_shift])
            continue

        shifted_ssn = np.array(shifted_ssn, dtype=float)
        shifted_ssn_err = np.array(shifted_ssn_err, dtype=float)
        n_pairs = len(flux_indices)
        mean_shifted = np.mean(shifted_ssn)
        shift_dx = (shifted_ssn - mean_shifted) / mean_shifted
        eshift_dx = 2.0 * shifted_ssn_err / mean_shifted
        shift_dy = dy[:, flux_indices].copy()
        eshift_dy = edy[:, flux_indices].copy()
        pippo_xyr_ssn = period_times[flux_indices] + shift * 0.5
        pippo_xyr_flu = period_times[flux_indices].copy()
        zeros = ary.array('d', [0.] * n_pairs)

        cpippo = new_canvas("cpippo_%+d" % shift, "cpippo_%+d" % shift)
        cpippo.Divide(1, 4, 1.e-3, 1.e-3)
        for ipad in range(1, 5):
            pad = cpippo.cd(ipad)
            pad.SetGrid()
            pad.SetTopMargin(0.07)
            pad.SetBottomMargin(0.18)
            pad.SetLeftMargin(0.15)
            pad.SetRightMargin(0.05)
        cpippo.cd(1)
        grpippo_ssn = keep(ROOT.TGraphErrors(n_pairs, array_d(pippo_xyr_ssn),
                                             array_d(shift_dx), zeros,
                                             array_d(eshift_dx)))
        grpippo_ssn.SetMarkerStyle(8)
        grpippo_ssn.SetMarkerColor(9)
        grpippo_ssn.SetMarkerSize(1.3)
        grpippo_ssn.SetLineColor(9)
        grpippo_ssn.SetLineWidth(2)
        grpippo_ssn.Draw("APL")

        for iene in range(NENERGY):
            cpippo.cd(iene+2)
            grpippo_flux = keep(ROOT.TGraphErrors(n_pairs,
                                    array_d(pippo_xyr_flu), array_d(shift_dy[iene]),
                                    zeros, array_d(eshift_dy[iene])))
            grpippo_flux.SetMarkerStyle(8)
            grpippo_flux.SetMarkerColor(633)
            grpippo_flux.SetMarkerSize(1.3)
            grpippo_flux.SetLineColor(633)
            grpippo_flux.SetLineWidth(2)
            grpippo_flux.Draw("APL")
            r, pval = safe_corr(shift_dx, shift_dy[iene])
            corr_facs[iene][indx_shift] = r
            pval_corr_facs[iene][indx_shift] = pval
            print("SSN shift %+0.1f year, E=%i GeV: r=%g p=%g" %
                  (shifts[indx_shift], ENERGY_GEV[iene], r, pval))

    for iene in range(NENERGY):
        ccorr = new_canvas("ccorr%i" % iene, "ccorr%i" % iene)
        ccorr.SetGrid()
        good = [idx for idx in range(Nshifts)
                if np.isfinite(corr_facs[iene, idx])]
        # Come Irene, il p-value e' riportato su y nel TGraphErrors.
        # NON e' un'incertezza statistica su Pearson r: grafico mantenuto fedele.
        gr_cor = keep(ROOT.TGraphErrors(len(good),
                        array_d([shifts[k] for k in good]),
                        array_d([corr_facs[iene, k] for k in good]),
                        ary.array('d', [0.]*len(good)),
                        array_d([pval_corr_facs[iene, k] for k in good])))
        gr_cor.SetTitle("")
        gr_cor.GetXaxis().SetTitle("SSN shift [years] (0.5-year steps)")
        gr_cor.GetYaxis().SetTitle("Corr. Factor at %i GeV" % ENERGY_GEV[iene])
        gr_cor.GetXaxis().SetTitleSize(0.04)
        gr_cor.GetXaxis().SetTitleOffset(0.9)
        gr_cor.GetYaxis().SetTitleSize(0.04)
        gr_cor.GetYaxis().SetTitleOffset(1.5)
        gr_cor.SetLineWidth(1)
        gr_cor.SetMarkerStyle(23)
        gr_cor.SetMarkerSize(1.8)
        gr_cor.SetMarkerColor(2)
        gr_cor.Draw("APL")

    # -----------------------------------------------------------------------
    # ctime: COME FIGURA, 4 pannelli (SSN e 50/80/125 GeV)
    # -----------------------------------------------------------------------
    ctime = new_canvas("ctime", "ctime - DAMPE p+He", 1100, 1400)
    ctime.Divide(1, 4, 1.e-3, 1.e-3)
    for ipad in range(1, 5):
        pad = ctime.cd(ipad)
        pad.SetGrid()
        pad.SetTopMargin(0.07)
        pad.SetBottomMargin(0.18 if ipad != 1 else 0.20)
        pad.SetLeftMargin(0.15)
        pad.SetRightMargin(0.05)

    x_moth = ary.array('d', [0.] * Nmonths)
    ex_moth = ary.array('d', [0.] * Nmonths)
    x_halfyear = ary.array('d', [0.] * Nperiods)
    ex_halfyear = ary.array('d', [0.] * Nperiods)
    for idx in range(Nmonths):
        x_moth[idx] = FIRST_YEAR + (idx + 0.5) / 12.0
    for idx in range(Nperiods):
        x_halfyear[idx] = FIRST_YEAR + (idx + 0.5) / 2.0

    ctime.cd(1)
    gr_month_ssn = keep(ROOT.TGraphErrors(Nmonths, x_moth, array_d(dmp_ssn),
                                         ex_moth, ex_moth))
    gr_month_ssn.SetTitle("")
    gr_month_ssn.SetMarkerStyle(8)
    gr_month_ssn.SetMarkerSize(0.7)
    gr_month_ssn.SetMarkerColor(30)
    gr_month_ssn.SetLineColor(29)
    gr_month_ssn.GetYaxis().SetTitle("SunSpot Number")
    gr_month_ssn.GetXaxis().SetTitle("year")
    gr_month_ssn.GetXaxis().SetLabelSize(0.075)
    gr_month_ssn.GetXaxis().SetTitleSize(0.08)
    gr_month_ssn.GetXaxis().SetTitleOffset(0.9)
    gr_month_ssn.GetXaxis().SetLabelOffset(0.018)
    gr_month_ssn.GetYaxis().SetLabelSize(0.075)
    gr_month_ssn.GetYaxis().SetTitleSize(0.10)
    gr_month_ssn.GetYaxis().SetTitleOffset(0.7)
    gr_month_ssn.GetYaxis().SetNdivisions(6, True)
    gr_month_ssn.GetXaxis().SetLimits(2016.0, 2026.0)
    gr_month_ssn.GetXaxis().SetNdivisions(510)
    gr_month_ssn.GetYaxis().SetRangeUser(0., 240.)
    gr_month_ssn.Draw("APL")

    end_24cycle = 2019.96
    peak_25cycle = 2024.62
    draw_line(end_24cycle, 0., 240., 12, 9, 1)
    draw_line(peak_25cycle, 0., 240., 12, 9, 1)
    latex25 = keep(ROOT.TLatex(2021.5, 205, "Solar cycle 25"))
    latex25.SetTextSize(0.09)
    latex25.SetTextColor(9)
    latex25.Draw()
    latex24 = keep(ROOT.TLatex(2017.1, 205, "Solar cycle 24"))
    latex24.SetTextSize(0.09)
    latex24.SetTextColor(9)
    latex24.Draw()

    gr_halfyear_ssn = keep(ROOT.TGraphErrors(Nperiods, x_halfyear,
                                            xssn, ex_halfyear, exssn))
    gr_halfyear_ssn.SetTitle("")
    gr_halfyear_ssn.SetMarkerStyle(8)
    gr_halfyear_ssn.SetMarkerColor(9)
    gr_halfyear_ssn.SetMarkerSize(1.3)
    gr_halfyear_ssn.SetLineColor(9)
    gr_halfyear_ssn.SetLineWidth(2)
    gr_halfyear_ssn.Draw("P SAME")

    sin_ssn = keep(ROOT.TF1("sin_ssn", sin_fit,
                            FIRST_YEAR, FIRST_YEAR + float(Nyrs), 4))
    sin_ssn.SetParameter(0, float(np.mean(xssn)))
    sin_ssn.SetParameter(1, float((max(xssn) - min(xssn))/2.))
    sin_ssn.FixParameter(2, 11.)
    sin_ssn.SetParameter(3, 2017.21)
    sin_ssn.SetParLimits(3, 2017., 2017.5)
    sin_ssn.SetLineColor(871)
    sin_ssn.SetLineWidth(2)
    gr_halfyear_ssn.Fit(sin_ssn, "Q0")
    ssn_min_sinX = x_halfyear[int(np.argmin(xssn))]
    ssn_max_X = x_halfyear[int(np.argmax(xssn))]
    draw_line(ssn_min_sinX, 0., 240., 871)
    draw_line(ssn_max_X, 0., 240., 871)

    sin_fluxes = []
    flu_max_sinX = np.zeros(NENERGY)
    for ipt in range(NENERGY):
        ctime.cd(2 + ipt)
        # Curve mensili: stesse grandezze 'flux' di Irene; stile della figura.
        monthly_yvals = array_d([h.GetBinContent(ENERGY_BINS[ipt])
                                 for h in hother_unf_monthly])
        monthly_yerrs = array_d([h.GetBinError(ENERGY_BINS[ipt])
                                 for h in hother_unf_monthly])
        gr_monthly_flux = keep(ROOT.TGraphErrors(
            Nmonths, x_moth, monthly_yvals, ex_moth, monthly_yerrs))
        gr_monthly_flux.SetTitle("")
        gr_monthly_flux.SetMarkerStyle(8)
        gr_monthly_flux.SetMarkerSize(0.7)
        gr_monthly_flux.SetMarkerColor(633)
        gr_monthly_flux.SetLineColor(ROOT.kGray + 2)
        gr_monthly_flux.SetLineStyle(7)
        gr_monthly_flux.SetLineWidth(1)
        gr_monthly_flux.GetXaxis().SetTitle("year")
        gr_monthly_flux.GetYaxis().SetTitle("at E = %i GeV" % ENERGY_GEV[ipt])
        gr_monthly_flux.GetXaxis().SetTitleSize(0.08)
        gr_monthly_flux.GetXaxis().SetTitleOffset(0.9)
        gr_monthly_flux.GetXaxis().SetLabelSize(0.075)
        gr_monthly_flux.GetXaxis().SetLabelOffset(0.015)
        gr_monthly_flux.GetYaxis().SetTitleSize(0.085)
        gr_monthly_flux.GetYaxis().SetTitleOffset(0.65)
        gr_monthly_flux.GetYaxis().SetLabelSize(0.075)
        gr_monthly_flux.GetYaxis().SetNdivisions(6, True)
        gr_monthly_flux.GetXaxis().SetLimits(2016.0, 2026.0)
        gr_monthly_flux.GetXaxis().SetNdivisions(510)
        # Range dinamico per il flusso p+He, diverso da all-particle.
        ymin, ymax = y_range([monthly_yvals[i] - monthly_yerrs[i]
                              for i in range(Nmonths)] +
                             [monthly_yvals[i] + monthly_yerrs[i]
                              for i in range(Nmonths)], 0.10)
        gr_monthly_flux.GetYaxis().SetRangeUser(ymin, ymax)
        gr_monthly_flux.Draw("APL")

        gr_halfyear_flux = keep(ROOT.TGraph(Nperiods, x_halfyear,
                                           array_d(yy[ipt])))
        gr_halfyear_flux.SetLineStyle(7)
        gr_halfyear_flux.SetLineWidth(1)
        gr_halfyear_flux.SetTitle("")
        gr_halfyear_flux2 = keep(ROOT.TGraphErrors(Nperiods, x_halfyear,
                             array_d(yy[ipt]), ex_halfyear, array_d(eyy[ipt])))
        gr_halfyear_flux2.SetTitle("")
        gr_halfyear_flux2.SetMarkerStyle(21)  # quadrato pieno
        gr_halfyear_flux2.SetMarkerColor(ROOT.kBlue + 1)
        gr_halfyear_flux2.SetMarkerSize(1.1)
        gr_halfyear_flux2.SetLineColor(ROOT.kBlue + 1)
        gr_halfyear_flux2.SetLineWidth(2)
        gr_halfyear_flux2.Draw("P SAME")

        # Stessa funzione e stessi intervalli di periodicita' di Irene.
        # Amplitudini inizializzate con la scala reale di p+He.
        sin_flu = keep(ROOT.TF1("sin_flu%i" % ipt, double_sin_fit,
                                2019., FIRST_YEAR+float(Nyrs), 7))
        avg_flu = float(np.mean(yy[ipt]))
        amp_flu = max((max(yy[ipt]) - min(yy[ipt]))/2., abs(avg_flu)*0.001)
        sin_flu.SetParameter(0, avg_flu)
        sin_flu.SetParameter(1, amp_flu)
        sin_flu.FixParameter(2, 11.)
        sin_flu.SetParameter(3, 2020.)
        sin_flu.SetParLimits(3, 2016., 2038.)
        sin_flu.SetParameter(4, amp_flu / 3.0)
        sin_flu.SetParLimits(4, amp_flu/100., amp_flu * 10.)
        sin_flu.SetParameter(5, 5.)
        sin_flu.SetParLimits(5, 4., 6.)
        sin_flu.SetParameter(6, 2020.)
        sin_flu.SetLineColor(632)
        sin_flu.SetLineWidth(2)
        gr_halfyear_flux2.Fit(sin_flu, "Q0")
        sin_flu.SetRange(2016., 2026.)
        sin_fluxes.append(sin_flu)
        flu_max_sinX[ipt] = x_halfyear[int(np.argmax(yy[ipt]))]
        print("Flux max at %d GeV, semester center: %.2f" %
              (ENERGY_GEV[ipt], flu_max_sinX[ipt]))
        # La riga verticale resta, ma la legenda NON copre il 2023-2024.
        draw_line(flu_max_sinX[ipt], ymin, ymax, 632)

        # Modifica di layout richiesta: legenda a sinistra (2019 circa).
        legend = keep(ROOT.TLegend(0.30, 0.26, 0.65, 0.42))
        legend.AddEntry(gr_halfyear_flux2,
                        "DAMPE p+He at %i GeV (Stat. err.)" % ENERGY_GEV[ipt],
                        "Pe")
        legend.AddEntry(0, "Min. obs. statistics %i" %
                        int(flux_PTstatistics[ipt]), "")
        legend.SetFillStyle(1001)
        legend.SetFillColor(ROOT.kWhite)
        legend.SetBorderSize(1)
        legend.SetTextSize(0.055)
        legend.Draw()

    # Titolo verticale comune, come nell'originale.
    ctime.cd()
    pad_title = keep(ROOT.TPad("pad_title", "", 0.005, 0.01, 0.2, 0.8))
    pad_title.SetFillStyle(4000)
    pad_title.Draw()
    pad_title.cd()
    latexytitle = keep(ROOT.TLatex(0.14, 0.32,"#bf{#Phi(E) [GeV^{-1} (m^{2} s sr)^{-1}]}"))
    latexytitle.SetTextAngle(90)
    latexytitle.SetTextSize(0.16)
    latexytitle.SetTextColor(1)
    latexytitle.Draw()
    pad_title.Update()

    # -----------------------------------------------------------------------
    # cpol: stessi 3 grafici polari dell'originale, ora a 50/80/125 GeV
    # -----------------------------------------------------------------------
    hlight0 = ary.array('d', [0.] * 12)
    Dtheta_12yr = array_d([iyr * (math.pi/5.5) for iyr in range(12)])
    eDtheta_12yr = ary.array('d', [0.] * 12)
    Dtheta_yr = array_d([(v - FIRST_YEAR) * (2.*math.pi/11.)
                         for v in period_times])
    eDtheta_yr = ary.array('d', [0.] * Nperiods)
    pol_axis_up_orig = [0.03, 0.025, 0.015]
    fac_reduce_ssn_var = [1.5e-2, 1.e-2, 6.e-3]

    gr_highlight0 = keep(ROOT.TGraphPolar(12, Dtheta_12yr, hlight0,
                                         eDtheta_12yr, eDtheta_12yr))
    gr_highlight0.SetTitle("")
    gr_highlight0.SetLineColor(1)
    gr_highlight0.SetLineWidth(2)
    for ipt in range(NENERGY):
        cpol = new_canvas("cpol%i" % ipt, "cpol%i" % ipt, 1100, 1100)
        cpol.SetGrid()
        cpol.SetTopMargin(0.05)
        cpol.SetBottomMargin(0.05)
        cpol.SetLeftMargin(0.05)
        cpol.SetRightMargin(0.05)
        rmax = max(pol_axis_up_orig[ipt],
                   1.3 * max(abs(v) for v in dy[ipt]))
        gr_polDy = keep(ROOT.TGraphPolar(Nperiods, Dtheta_yr,
                         array_d(dy[ipt]), eDtheta_yr, array_d(edy[ipt])))
        gr_polDy.SetTitle("")
        gr_polDy2 = keep(ROOT.TGraphPolar(Nperiods))
        gr_polDy2.SetTitle("")
        gr_polDy2.SetMarkerStyle(8)
        gr_polDy2.SetMarkerColor(633)
        gr_polDy2.SetMarkerSize(2.)
        gr_polDy2.SetLineColor(633)
        gr_polDy2.SetLineWidth(3)
        gr_polDy.SetLineColor(633)
        gr_polDy.SetLineStyle(9)
        gr_polDy.SetLineWidth(3)
        gr_polDy.Draw("APOL")

        dx_dim = array_d([dx[k] * fac_reduce_ssn_var[ipt]
                          for k in range(Nperiods)])
        edx_dim = array_d([edx[k] * fac_reduce_ssn_var[ipt]
                           for k in range(Nperiods)])
        gr_polDx = keep(ROOT.TGraphPolar(Nperiods, Dtheta_yr, dx_dim,
                         eDtheta_yr, edx_dim))
        gr_polDx.SetTitle("")
        gr_polDx2 = keep(ROOT.TGraphPolar(Nperiods))
        gr_polDx2.SetTitle("")
        gr_polDx2.SetMarkerStyle(8)
        gr_polDx2.SetMarkerColor(9)
        gr_polDx2.SetMarkerSize(2.)
        gr_polDx2.SetLineColor(9)
        gr_polDx2.SetLineWidth(3)
        gr_polDx.SetLineColor(9)
        gr_polDx.SetLineWidth(3)
        gr_polDx.SetLineStyle(9)
        cpol.Update()
        polgram = gr_polDy.GetPolargram()
        if polgram:
            polgram.SetTextSize(0.01)
            polgram.SetPolarLabelSize(0.001)
            polgram.SetNdivPolar(11)
            for iyear in range(11):
                polgram.SetPolarLabel(iyear, "%d" % (2016+iyear))
            polgram.SetPolarLabelSize(0.001)
            polgram.SetPolarOffset(0.015)
            polgram.SetToRadian()
            polgram.SetNdivRadial(210)
            polgram.SetRadialLabelSize(0.)
            polgram.SetRangeRadial(-rmax, rmax)
        gr_polDy2.Draw("PE")
        gr_polDx.Draw("LPE")
        gr_polDx2.Draw("PE")
        gr_highlight0.Draw("C")

        # Marker polari delle date di minimo SSN e massimo flusso.
        zero = ary.array('d', [0.])
        rerr = ary.array('d', [rmax])
        theta_ssn = array_d([(math.pi / 5.5) * (ssn_min_sinX - FIRST_YEAR)])
        gr_ssn_min = keep(ROOT.TGraphPolar(1, theta_ssn, zero, zero, rerr))
        gr_ssn_min.SetTitle("")
        gr_ssn_min.SetLineColor(871)
        gr_ssn_min.SetLineWidth(3)
        gr_ssn_min.SetLineStyle(9)
        gr_ssn_min.SetMarkerSize(0)
        gr_ssn_min.Draw("PE")
        theta_flux = array_d([(math.pi / 5.5) * (flu_max_sinX[ipt] - FIRST_YEAR)])
        gr_flu_max = keep(ROOT.TGraphPolar(1, theta_flux, zero, zero, rerr))
        gr_flu_max.SetTitle("")
        gr_flu_max.SetLineColor(632)
        gr_flu_max.SetLineWidth(3)
        gr_flu_max.SetLineStyle(9)
        gr_flu_max.SetMarkerSize(0)
        gr_flu_max.Draw("PE")

        lege = keep(ROOT.TLegend(0.1, 0.96, 0.9, 1.))
        lege.SetTextAlign(12)
        lege.SetNColumns(2)
        lege.AddEntry(gr_polDy2,
                      "#Delta(#Phi_{%i GeV}) / <#Phi_{%i GeV}>" %
                      (ENERGY_GEV[ipt], ENERGY_GEV[ipt]), "pe")
        lege.AddEntry(gr_polDx2,
                      "(#Delta(SunSpot N.) / <SunSpot N.>) x (%.0f 10^{-3})" %
                      (fac_reduce_ssn_var[ipt] * 1.e3), "pe")
        lege.Draw("same")
        cpol.Update()

    # -----------------------------------------------------------------------
    # c3d: E x year x variazione relativa del flusso (stesso ROOT TGraph2D)
    # -----------------------------------------------------------------------
    c3d = new_canvas("c3d", "c3d")
    c3d.SetRightMargin(0.1)
    c3d.SetLeftMargin(0.1)
    flux_array = np.array(fluxes_2d).reshape(Nperiods, NENERGY)
    flux_mean_alongYrs = np.mean(flux_array, axis=0)
    relvar_fluxes = (flux_array - flux_mean_alongYrs) / flux_mean_alongYrs
    gr3D = keep(ROOT.TGraph2D(len(energies_2d), array_d(energies_2d),
                              array_d(years_2d), array_d(relvar_fluxes.flatten())))
    gr3D.SetTitle(";log_{10}(E [GeV]);Year;#Delta(#Phi) / <#Phi>")
    gr3D.GetXaxis().SetTitleOffset(2)
    gr3D.GetYaxis().SetLimits(FIRST_YEAR, FIRST_YEAR+Nyrs)
    gr3D.GetYaxis().SetTitleOffset(2)
    gr3D.GetZaxis().SetTitleOffset(2)
    gr3D.Draw("surf1z")
    c3d.Update()

    # -----------------------------------------------------------------------
    # c_asb_relvar e c_abs_relvar3: max delle deviazioni ai 3 bin p+He
    # -----------------------------------------------------------------------
    c_asb_relvar = new_canvas("c_asb_relvar", "c_asb_relvar")
    c_asb_relvar.SetGrid()
    c_asb_relvar.SetLogx()
    c_asb_relvar.SetRightMargin(0.1)
    c_asb_relvar.SetLeftMargin(0.1)
    gr_abs_relvar = keep(ROOT.TGraphErrors(NENERGY, array_d(dene_max),
                         array_d(dy_max), array_d(edene_max), array_d(edy_max)))
    gr_abs_relvar.SetTitle(";E [GeV];max ( | #Delta#Phi_{E} / <#Phi_{E}> | )")
    gr_abs_relvar.GetYaxis().SetRangeUser(0.,
                                         max(0.02, 1.2*max(dy_max+edy_max)))
    gr_abs_relvar.GetXaxis().SetLimits(3.e1, 2.e2)
    gr_abs_relvar.GetXaxis().SetTitleOffset(1.2)
    gr_abs_relvar.SetMarkerStyle(8)
    gr_abs_relvar.SetMarkerColor(633)
    gr_abs_relvar.SetMarkerSize(2.)
    gr_abs_relvar.SetLineColor(633)
    gr_abs_relvar.SetLineWidth(3)
    gr_abs_relvar.Draw("APE")

    c_abs_relvar3 = new_canvas("c_abs_relvar3", "c_abs_relvar3")
    c_abs_relvar3.SetGrid()
    c_abs_relvar3.SetLogx()
    c_abs_relvar3.SetRightMargin(0.1)
    c_abs_relvar3.SetLeftMargin(0.1)
    gr_abs_relvar3 = keep(ROOT.TGraphErrors(NENERGY, array_d(dene_max),
                         array_d(dy_max*100.), array_d(edene_max),
                         array_d(edy_max*100.)))
    gr_abs_relvar3.SetTitle(";E [GeV];max ( | #Delta#Phi_{E} / <#Phi_{E}> | ) %")
    gr_abs_relvar3.GetYaxis().SetRangeUser(0.,
                             max(2.0, 120.*max(dy_max+edy_max)))
    gr_abs_relvar3.GetXaxis().SetLimits(3.e1, 2.e2)
    gr_abs_relvar3.GetXaxis().SetTitleOffset(1.2)
    gr_abs_relvar3.GetXaxis().SetNdivisions(110, True)
    gr_abs_relvar3.GetXaxis().SetMoreLogLabels()
    gr_abs_relvar3.GetXaxis().SetNoExponent()
    gr_abs_relvar3.SetMarkerStyle(8)
    gr_abs_relvar3.SetMarkerColor(633)
    gr_abs_relvar3.SetMarkerSize(2.)
    gr_abs_relvar3.SetLineColor(633)
    gr_abs_relvar3.SetLineWidth(3)
    gr_abs_relvar3.Draw("APE")

    # Nessun plot/CSV aggiuntivo di Pearson, ne' grafici nuovi.
    # Conserva le figure esistenti; di default si salva soltanto ctime.
    if SAVE_ALL_ORIGINAL_CANVASES:
        for cv in original_canvases:
            save_canvas(cv)
    elif SAVE_CTIME:
        save_canvas(ctime)
    print("Analisi p+He terminata. Output:", OUTPUT_DIR)
    if not BATCH:
        raw_input("Press enter..")


if __name__ == "__main__":
    main()
