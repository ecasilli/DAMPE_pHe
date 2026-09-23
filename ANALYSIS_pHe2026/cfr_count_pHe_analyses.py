import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# ==========================================================
# FILE INPUT
# ==========================================================

file_pHe = "TXT_FILES/DAMPE_p+He_120M_paperDraft.dat"
file_paul_He = "GENEVA_pHe_FILES/Helium_files_20260922/Counts-and-background.txt"
file_andrii_p = "GENEVA_pHe_FILES/COUNTS_PROTON_SEP2026_flux_noescale_ekin_p_2026.root"


# ==========================================================
# 1. LETTURA DEI MIEI CONTEGGI p+He
# ==========================================================

# Colonne:
# Emean_GeV Emin_GeV Emax_GeV Counts Flux FluxErrStat

data_pHe = np.loadtxt(file_pHe, comments="#")

Emean = data_pHe[:, 0]
Emin = data_pHe[:, 1]
Emax = data_pHe[:, 2]

counts_pHe = data_pHe[:, 3]

# ==========================================================
# 2. LETTURA DEI CONTEGGI DI PAUL (He)
# ==========================================================

# Colonne:
# E_BGO_min E_BGO_max E_BGO_center Counts Background (%)

dati_paul_He = np.loadtxt(file_paul_He, skiprows=1)

Emin_paul_He = dati_paul_He[:, 0]
Emax_paul_He = dati_paul_He[:, 1]
Ecenter_paul_He = dati_paul_He[:, 2]

counts_paul_He = dati_paul_He[:, 3]

background_paul_He = dati_paul_He[:, 4]

# Conteggi dopo sottrazione del background
counts_paul_He_bgsub = counts_paul_He * (1.0 - background_paul_He / 100.0)

# ==========================================================
# 3. LETTURA DEI CONTEGGI DI ANDRII (p)
# ==========================================================

f = ROOT.TFile.Open(file_andrii_p, "READ")

if not f or f.IsZombie():
    raise RuntimeError("ERROR: cannot open ROOT file")

h_data = f.Get("data_pregion")
h_rslt = f.Get("rslt_pregion")

if not h_data or not h_rslt:
    raise RuntimeError("ERROR: cannot find histograms")

# ==========================================================
# 4. ALLINEAMENTO DEI BIN ENERGETICI
# ==========================================================

# I bin dei tre file sono praticamente coincidenti,
# ma i bordi degli istogrammi ROOT e quelli dei TXT
# hanno precisioni numeriche differenti.
#
# Per questo motivo confrontiamo i bordi con una
# tolleranza relativa.

def same_bin(emin1, emax1, emin2, emax2, tol=1e-3):

    return (abs(emin1 - emin2) / emin1 < tol and abs(emax1 - emax2) / emax1 < tol)

# Vettori finali dei conteggi, ordinati secondo
# i bin energetici della mia analisi p+He.

x = []

y_pHe = []

y_paul_He = []
y_andrii_p = []

y_paul_He_bg = []
y_andrii_p_bg = []

err_pHe = []
err_paul_He = []
err_andrii_p = []

err_paul_He_bg = []
err_andrii_p_bg = []


for i in range(len(Emean)):

    e_min = Emin[i]
    e_max = Emax[i]

    # ------------------------------------------------------
    # Cerco il bin corrispondente nei dati di Paul
    # ------------------------------------------------------

    bin_paul_He = None

    for j in range(len(counts_paul_He)):

        if same_bin(e_min, e_max, Emin_paul_He[j], Emax_paul_He[j]):
            bin_paul_He = j
            break
    if bin_paul_He is None:
        raise RuntimeError("Bin not found: {} - {} GeV".format(e_min, e_max))

    # ------------------------------------------------------
    # Cerco il bin corrispondente negli istogrammi di Andrii
    # ------------------------------------------------------

    bin_andrii_p = None

    for j in range(1, h_data.GetNbinsX() + 1):

        # Conversione MeV -> GeV
        e_min_andrii_p = (h_data.GetXaxis().GetBinLowEdge(j) / 1000.0)
        e_max_andrii_p = (h_data.GetXaxis().GetBinUpEdge(j) / 1000.0)

        if same_bin(e_min, e_max, e_min_andrii_p, e_max_andrii_p):
            bin_andrii_p = j
            break

    if bin_andrii_p is None:
        raise RuntimeError("Bin not found: {} - {} GeV".format(e_min, e_max))

    # ------------------------------------------------------
    # Riempio i vettori dei conteggi
    # ------------------------------------------------------

    x.append(Emean[i])
    y_pHe.append(counts_pHe[i])
    y_paul_He.append(counts_paul_He[bin_paul_He])
    y_andrii_p.append(h_data.GetBinContent(bin_andrii_p))
    y_paul_He_bg.append(counts_paul_He_bgsub[bin_paul_He])
    y_andrii_p_bg.append(h_rslt.GetBinContent(bin_andrii_p))

    err_pHe.append(np.sqrt(counts_pHe[i]))
    err_paul_He.append(np.sqrt(counts_paul_He[bin_paul_He]))
    err_andrii_p.append(h_data.GetBinError(bin_andrii_p))

    err_paul_He_bg.append(np.sqrt(counts_paul_He[bin_paul_He]))

    # Andrii: utilizzo gli errori memorizzati nell'istogramma
    # corretto, se presenti
    if h_rslt.GetSumw2N() > 0:
        err_andrii_p_bg.append(h_rslt.GetBinError(bin_andrii_p))
    else:
        # Approssimazione in assenza di errori propagati
        # nell'istogramma corretto
        err_andrii_p_bg.append(np.sqrt(h_data.GetBinContent(bin_andrii_p)))

# ==========================================================
# 5. SOMMA DEI CONTEGGI DELLE ANALISI SINGOLE
# ==========================================================

y_sum = []
y_sum_bg = []

err_sum = []
err_sum_bg = []

# Residui
y_residual = []
y_residual_bg = []

# Errori sui residui
err_residual = []
err_residual_bg = []

for i in range(len(x)):

    # Somma dei conteggi totali
    n_sum = y_paul_He[i] + y_andrii_p[i]
    y_sum.append(n_sum)
    sigma_sum = np.sqrt(err_paul_He[i]**2 + err_andrii_p[i]**2)
    err_sum.append(sigma_sum)

    residual = y_pHe[i] - n_sum
    y_residual.append(residual)
    # Errore sul residuo (correlazioni trascurate)
    sigma_residual = np.sqrt(err_pHe[i]**2 + sigma_sum**2)
    err_residual.append(sigma_residual)

    # Somma dei conteggi dopo sottrazione del background
    n_sum_bg = (y_paul_He_bg[i] + y_andrii_p_bg[i])
    y_sum_bg.append(n_sum_bg)
    sigma_sum_bg = np.sqrt(err_paul_He_bg[i]**2 + err_andrii_p_bg[i]**2)
    err_sum_bg.append(sigma_sum_bg)

    residual_bg = y_pHe[i] - n_sum_bg
    y_residual_bg.append(residual_bg)
    # Errore sul residuo BG-subtracted
    sigma_residual_bg = np.sqrt(err_pHe[i]**2 + sigma_sum_bg**2)
    err_residual_bg.append(sigma_residual_bg)


# ==========================================================
# STAMPA DEI RESIDUI
# ==========================================================

print("\n========================================================")
print("               CONFRONTO DEI CONTEGGI")
print("========================================================")

print(
    "{:>12} {:>15} {:>15} {:>15} {:>15}".format(
        "E [GeV]",
        "N_pHe",
        "N_p+N_He",
        "Residual",
        "Residual [%]"
    )
)

for i in range(len(x)):

    # Residuo assoluto
    residual = y_residual[i]

    # Residuo relativo
    if y_pHe[i] > 0:
        residual_rel = 100.0 * residual / y_pHe[i]
    else:
        residual_rel = float("nan")

    print(
        "{:12.2f} {:15.2f} {:15.2f} {:15.2f} {:15.2f}".format(
            x[i],
            y_pHe[i],
            y_sum[i],
            residual,
            residual_rel
        )
    )


print("\n========================================================")
print("        CONFRONTO CON SOTTRAZIONE DEL BACKGROUND")
print("========================================================")

print(
    "{:>12} {:>15} {:>15} {:>15} {:>15}".format(
        "E [GeV]",
        "N_pHe",
        "N_p+N_He BG",
        "Residual",
        "Residual [%]"
    )
)

for i in range(len(x)):

    # Residuo assoluto
    residual_bg = y_residual_bg[i]

    # Residuo relativo
    if y_pHe[i] > 0:
        residual_bg_rel = 100.0 * residual_bg / y_pHe[i]
    else:
        residual_bg_rel = float("nan")

    print(
        "{:12.2f} {:15.2f} {:15.2f} {:15.2f} {:15.2f}".format(
            x[i],
            y_pHe[i],
            y_sum_bg[i],
            residual_bg,
            residual_bg_rel
        )
    )

# ==========================================================
# RESIDUI RELATIVI IN PERCENTUALE
# ==========================================================

y_residual_percent = []
y_residual_bg_percent = []

err_residual_percent = []
err_residual_bg_percent = []

for i in range(len(x)):

    N_pHe = y_pHe[i]

    N_sum = y_sum[i]
    N_sum_bg = y_sum_bg[i]

    sigma_pHe = err_pHe[i]

    sigma_sum = err_sum[i]
    sigma_sum_bg = err_sum_bg[i]

    if N_pHe <= 0:

        y_residual_percent.append(float("nan"))
        y_residual_bg_percent.append(float("nan"))

        err_residual_percent.append(float("nan"))
        err_residual_bg_percent.append(float("nan"))

        continue

    # ------------------------------------------------------
    # RESIDUI PERCENTUALI DEI CONTEGGI TOTALI
    # ------------------------------------------------------

    residual_percent = 100.0 * (N_pHe - N_sum) / N_pHe
    y_residual_percent.append(residual_percent)
    # Propagazione degli errori assumendo indipendenza
    sigma_residual_percent = 100.0 * np.sqrt((sigma_sum / N_pHe)**2 + (N_sum * sigma_pHe / N_pHe**2)**2)
    err_residual_percent.append(sigma_residual_percent)

    # ------------------------------------------------------
    # RESIDUI PERCENTUALI DOPO SOTTRAZIONE DEL BACKGROUND
    # ------------------------------------------------------

    residual_bg_percent = 100.0 * (N_pHe - N_sum_bg) / N_pHe
    y_residual_bg_percent.append(residual_bg_percent)
    sigma_residual_bg_percent = 100.0 * np.sqrt((sigma_sum_bg / N_pHe)**2 + (N_sum_bg * sigma_pHe / N_pHe**2)**2)
    err_residual_bg_percent.append(sigma_residual_bg_percent)

# ==========================================================
# 6. FUNZIONE PER CREARE I GRAFICI ROOT
# ==========================================================

def make_graph(name,counts_pHe,counts_He,counts_p,counts_sum,errors_pHe,errors_He,errors_p,errors_sum,residuals,errors_residual,bg_subtracted=False):

    c = ROOT.TCanvas(name, name, 1000, 900)

    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)

    pad_top = ROOT.TPad(name + "_top", "", 0.0, 0.30, 1.0, 1.0)

    pad_top.SetLeftMargin(0.13)
    pad_top.SetRightMargin(0.05)

    pad_top.SetTopMargin(0.06)
    pad_top.SetBottomMargin(0.03)

    pad_top.SetLogx()
    pad_top.SetLogy()

    pad_top.SetGridx()
    pad_top.SetGridy()

    pad_top.Draw()

    pad_bottom = ROOT.TPad(name + "_bottom", "", 0.0, 0.0, 1.0, 0.30)

    pad_bottom.SetLeftMargin(0.13)
    pad_bottom.SetRightMargin(0.05)

    pad_bottom.SetTopMargin(0.03)
    pad_bottom.SetBottomMargin(0.34)

    pad_bottom.SetLogx()

    pad_bottom.SetGridx()
    pad_bottom.SetGridy()

    pad_bottom.Draw()

    # ------------------------------------------------------
    # Funzione per costruire un TGraph
    # ------------------------------------------------------

    def make_root_tgraph(y, ey, color, marker, positive_only=True):

        points = []

        for energy, counts, error in zip(x, y, ey):
            if positive_only and counts <= 0:
                continue
            points.append((energy, counts, error))

        gr = ROOT.TGraphErrors(len(points))

        for i, (energy, counts, error) in enumerate(points):
            gr.SetPoint(i, energy, counts)
            gr.SetPointError( i, 0.0, error)

        gr.SetMarkerStyle(marker)

        gr.SetMarkerColor(color)
        gr.SetLineColor(color)

        gr.SetMarkerSize(1.1)

        gr.SetLineWidth(2)

        return gr

    # ------------------------------------------------------
    # Grafici dei conteggi
    # ------------------------------------------------------

    gr_pHe = make_root_tgraph(counts_pHe, errors_pHe, ROOT.kBlack, 20)
    gr_paul_He = make_root_tgraph( counts_He, errors_He, ROOT.kOrange + 7, 21)
    gr_andrii_p = make_root_tgraph( counts_p, errors_p, ROOT.kAzure + 2, 22)

    # ------------------------------------------------------
    # Grafico della somma p + He
    # ------------------------------------------------------

    gr_sum = make_root_tgraph(counts_sum, errors_sum, ROOT.kRed + 1, 1)
    gr_sum.SetLineWidth(3)

    # ======================================================
    # DISEGNO DEL PANNELLO SUPERIORE
    # ======================================================

    pad_top.cd()

    all_counts = [n for serie in ( counts_pHe, counts_He, counts_p, counts_sum) for n in serie if n > 0]

    ymin = min(all_counts) / 2.0
    ymax = max(all_counts) * 2.0
    xmin = min(Emin) / 1.1
    xmax = max(Emax) * 1.1

    frame_top = pad_top.DrawFrame(xmin,ymin,xmax,ymax)

    frame_top.SetTitle("")
    frame_top.GetYaxis().SetTitle("Counts")

    frame_top.GetYaxis().SetTitleSize(0.055)
    frame_top.GetYaxis().SetLabelSize(0.045)
    frame_top.GetYaxis().SetTitleOffset(1.05)

    # Nascondiamo le etichette dell'asse X superiore,
    # che saranno presenti nel pannello inferiore.

    frame_top.GetXaxis().SetLabelSize(0)
    frame_top.GetXaxis().SetTitleSize(0)

    # ------------------------------------------------------
    # Disegno della somma come linea continua
    # ------------------------------------------------------

    gr_sum.Draw("L SAME")

    # ------------------------------------------------------
    # Disegno dei conteggi come punti con barre di errore
    # ------------------------------------------------------

    gr_pHe.Draw("P E1 SAME")
    gr_paul_He.Draw("P E1 SAME")
    gr_andrii_p.Draw("P E1 SAME")

    # ======================================================
    # LEGENDA
    # ======================================================

    leg = ROOT.TLegend(0.48,0.68,0.93,0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.AddEntry(gr_pHe, "p+He analysis", "p")
    leg.AddEntry(gr_andrii_p, "p" #+ (" (BG subtracted)" if bg_subtracted else "")
        , "p")
    leg.AddEntry(gr_paul_He, "He" #+ (" (BG subtracted)" if bg_subtracted else "")
        , "p")
    leg.AddEntry(gr_sum, "N_{p} + N_{He}" # + (" (BG subtracted)" if bg_subtracted else "")
        , "l")
    leg.Draw()

    # ======================================================
    # PANNELLO INFERIORE: RESIDUI
    # ======================================================

    pad_bottom.cd()

    # ------------------------------------------------------
    # Costruzione del TGraphErrors dei residui
    # ------------------------------------------------------

    gr_residual = make_root_tgraph(residuals,errors_residual,ROOT.kBlack, 20, positive_only=False)
    
    '''
    print("\n====================================")
    print("RESIDUALS: {}".format(name))
    print("====================================")

    for i in range(gr_residual.GetN()):

        energy = ROOT.Double(0.0)
        residual = ROOT.Double(0.0)

        gr_residual.GetPoint(i, energy, residual)

        print(
            "E = {:.2f} GeV, Residual = {:+.4f}".format(
                float(energy),
                float(residual)
            )
        )
    '''
    # ------------------------------------------------------
    # Range dell'asse Y
    # ------------------------------------------------------

    # Consideriamo anche gli errori per evitare
    # che le barre escano dal pannello.

    #residual_min = min(r - e for r, e in zip(residuals, errors_residual))
    #residual_max = max(r + e for r, e in zip(residuals, errors_residual))
    residual_min = -20
    residual_max = 40

    # Includiamo sempre lo zero nel range verticale.

    residual_min = min(residual_min, 0.0)
    residual_max = max(residual_max, 0.0)
    residual_range = residual_max - residual_min

    # Protezione nel caso in cui tutti i residui
    # siano uguali.

    if residual_range == 0:
        residual_range = 1.0

    #residual_min -= 0.12 * residual_range
    #residual_max += 0.12 * residual_range

    # ------------------------------------------------------
    # Frame dei residui
    # ------------------------------------------------------

    frame_bottom = pad_bottom.DrawFrame(xmin, residual_min, xmax, residual_max)
    frame_bottom.SetTitle("")

    # ------------------------------------------------------
    # Asse X
    # ------------------------------------------------------

    frame_bottom.GetXaxis().SetTitle("Energy [GeV]")
    frame_bottom.GetXaxis().SetTitleSize(0.12)
    frame_bottom.GetXaxis().SetLabelSize(0.10)
    frame_bottom.GetXaxis().SetTitleOffset(1.1)

    # ------------------------------------------------------
    # Asse Y
    # ------------------------------------------------------

    frame_bottom.GetYaxis().SetTitle("Residuals [%]") #("N_{p+He} - (N_{p}+N_{He})")
    frame_bottom.GetYaxis().SetTitleSize(0.095)
    frame_bottom.GetYaxis().SetLabelSize(0.085)
    frame_bottom.GetYaxis().SetTitleOffset(0.55)
    frame_bottom.GetYaxis().SetNdivisions(505)

    # ------------------------------------------------------
    # Linea orizzontale a zero
    # ------------------------------------------------------

    zero_line = ROOT.TLine( xmin, 0.0, xmax, 0.0)

    zero_line.SetLineColor(ROOT.kRed + 1)
    zero_line.SetLineStyle(2)
    zero_line.SetLineWidth(2)

    zero_line.Draw("SAME")

    # ------------------------------------------------------
    # Disegno dei residui con le barre di errore
    # ------------------------------------------------------

    gr_residual.Draw("P SAME")

    # ======================================================
    # AGGIORNAMENTO DEL CANVAS
    # ======================================================

    c.cd()
    c.Modified()
    c.Update()

    # ======================================================
    # SALVATAGGIO
    # ======================================================

    c.SaveAs("PLOTS/" + name + ".png")
    c.SaveAs("PLOTS/" + name + ".pdf")

    # Manteniamo i riferimenti agli oggetti Python.
    # Questo evita che PyROOT possa eliminarli mentre
    # il canvas li sta ancora utilizzando.
    c._keepalive = (pad_top,pad_bottom,frame_top,frame_bottom,gr_pHe,gr_paul_He,gr_andrii_p,gr_sum,gr_residual,zero_line,leg)

    return c


# ==========================================================
# 7. GRAFICO DEI CONTEGGI TOTALI
# ==========================================================

c_tot = make_graph(
    "counts_comparison_total_residuals_percent",

    y_pHe,
    y_paul_He,
    y_andrii_p,
    y_sum,

    err_pHe,
    err_paul_He,
    err_andrii_p,
    err_sum,

    y_residual_percent,
    err_residual_percent,

    bg_subtracted=False
)

# ==========================================================
# 8. GRAFICO DEI CONTEGGI DOPO SOTTRAZIONE DEL BACKGROUND
# ==========================================================

c_bg = make_graph(
    "counts_comparison_bgsub_residuals_percent",

    y_pHe,
    y_paul_He_bg,
    y_andrii_p_bg,
    y_sum_bg,

    err_pHe,
    err_paul_He_bg,
    err_andrii_p_bg,
    err_sum_bg,

    y_residual_bg_percent,
    err_residual_bg_percent,

    bg_subtracted=True
)

# ==========================================================
# 9. SALVATAGGIO DEI GRAFICI IN UN FILE ROOT
# ==========================================================

fout = ROOT.TFile("ROOT_FILES/counts_comparison_residuals.root", "RECREATE")

c_tot.Write()
c_bg.Write()

fout.Close()
f.Close()
