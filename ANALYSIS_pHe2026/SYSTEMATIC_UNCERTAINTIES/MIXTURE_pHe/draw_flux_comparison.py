 
import sys
from ROOT import gStyle, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex, TLegend, TLine, TBox, TCanvas, gPad, kBlack, kGray, kRed, kBlue, kGreen, kAzure, kOrange, kMagenta, TPad, TGaxis
import array as ary
import numpy as np
import math

def make_flux_graph_DAMPE2024(filename, color, marker, size, alpha):
    Emean    = np.loadtxt(filename, skiprows=0, usecols=(0,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=0, usecols=(1,), unpack=True)
    Flux_stat= np.loadtxt(filename, skiprows=0, usecols=(2,), unpack=True)
    Ene_err  = np.loadtxt(filename, skiprows=0, usecols=(3,), unpack=True)
    Flux_sysA= np.loadtxt(filename, skiprows=0, usecols=(4,), unpack=True)
    Flux_sysH= np.loadtxt(filename, skiprows=0, usecols=(5,), unpack=True)

    Flux    = (Flux_2) * Emean**alpha
    Flux_err= (Flux_stat*Flux)/100.
    syst    = (Flux_sysA*Flux) 
    syst_had= (Flux_sysH*Flux)

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_err, Flux_err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    gr_in  = TGraphErrors(len(Emean), Emean, Flux, null, syst) #TGraphAsymmErrors(len(Emean)) #+2)
    gr_out = TGraphErrors(len(Emean), Emean, Flux, null, syst_had) #TGraphAsymmErrors(len(Emean)) #+2)
    gr_in.SetLineColor(0)
    gr_in.SetMarkerColor(color)
    gr_in.SetMarkerStyle(marker)
    gr_in.SetMarkerSize(size)
    gr_in.SetFillColor(17)
    gr_in.SetFillStyle(1001)
    #gr_in.SetPoint(0, Elow[0], Flux[0])
    #gr_in.SetPointError(0, 0, 0, syst[0], syst[0])
    
    gr_out.SetLineColor(0)
    gr_out.SetMarkerColor(color)
    gr_out.SetMarkerStyle(marker)
    gr_out.SetMarkerSize(size)
    gr_out.SetFillColor(18)
    gr_out.SetFillStyle(1001)
    #gr_out.SetPoint(0, Elow[0], Flux[0])
    #gr_out.SetPointError(0, 0, 0, syst_had[0], syst_had[0])

    return gr, gr_in, gr_out

def make_flux_graph_DAMPE2026(filename, color, marker, size, alpha):
    Emean    = np.loadtxt(filename, skiprows=0, usecols=(0,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=0, usecols=(1,), unpack=True)
    Flux_stat= np.loadtxt(filename, skiprows=0, usecols=(2,), unpack=True)
    Ene_err  = np.loadtxt(filename, skiprows=0, usecols=(3,), unpack=True)

    Flux    = (Flux_2) * Emean**alpha
    Flux_err= (Flux_stat) * Emean**alpha

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_err, Flux_err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr 

def make_flux_graph_pHe(filenameP, filenameHe, color, marker, size, alpha):
    # Qty   <E>  Elo  Eup   y   ystat_lo  ystat_up  ysyst_lo  ysyst_up  yerrtot_lo 
    # NEW: emin    emax      ene      flux      stat         ana         had        pow
    EmeanP       = np.loadtxt(filenameP, skiprows=1, usecols=(2,), unpack=True)
    Flux_2P      = np.loadtxt(filenameP, skiprows=1, usecols=(3,), unpack=True)
    Flux_stat_loP= np.loadtxt(filenameP, skiprows=1, usecols=(4,), unpack=True)
    Flux_stat_upP= np.loadtxt(filenameP, skiprows=1, usecols=(4,), unpack=True)

    FluxP    = (Flux_2P) * EmeanP**alpha
    Flux_err_loP= (Flux_stat_loP) * EmeanP**alpha
    Flux_err_upP= (Flux_stat_upP) * EmeanP**alpha


    # NEW: #E_min (Gev) E_max (GeV)  E_c (GeV)    Flux         Stat_err     Sys_err_ana  Sys_err_had  Sys_err_tot
    #      #doerrorbands
    EmeanHe       = np.loadtxt(filenameHe, skiprows=2, usecols=(2,), unpack=True)
    Flux_2He      = np.loadtxt(filenameHe, skiprows=2, usecols=(3,), unpack=True)
    Flux_stat_loHe= np.loadtxt(filenameHe, skiprows=2, usecols=(4,), unpack=True)
    Flux_stat_upHe= np.loadtxt(filenameHe, skiprows=2, usecols=(4,), unpack=True)

    FluxHe    = (Flux_2He) * EmeanHe**alpha
    Flux_err_loHe= (Flux_stat_loHe) * EmeanHe**alpha
    Flux_err_upHe= (Flux_stat_upHe) * EmeanHe**alpha

    null = np.zeros(len(EmeanP))

    assert np.allclose(EmeanP, EmeanHe, rtol=1e-3)

    Flux_sum = FluxP + FluxHe
    Stat_sum_lo = np.sqrt(Flux_err_loP**2 + Flux_err_loHe**2)
    Stat_sum_up = np.sqrt(Flux_err_upP**2 + Flux_err_upHe**2)

    Flux_stat_pHe_lo = np.sqrt( Flux_stat_loP**2 + Flux_stat_loHe**2)
    Flux_stat_pHe_up = np.sqrt( Flux_stat_upP**2 + Flux_stat_upHe**2)
    Flux_err_pHe_lo = (Flux_stat_pHe_lo) * EmeanP**alpha
    Flux_err_pHe_up = (Flux_stat_pHe_up) * EmeanP**alpha


    #gr = TGraphAsymmErrors(len(EmeanP), EmeanP, Flux_sum, null, null, Stat_sum_lo, Stat_sum_up)
    gr = TGraphAsymmErrors(len(EmeanP), EmeanP, Flux_sum, null, null, Flux_err_pHe_lo, Flux_err_pHe_up)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr 

def make_flux_graph_pHe_sum(filename, color, marker, size, alpha):
    Emean    = np.loadtxt(filename, skiprows=1, usecols=(2,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=1, usecols=(3,), unpack=True)
    Flux_stat= np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)

    Flux    = (Flux_2) * Emean**alpha
    Flux_err= (Flux_stat) * Emean**alpha

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_err, Flux_err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr 

def make_flux_graph_from_ROOT(filename, histname, color, marker, size, alpha):
    from ROOT import TFile
    
    f = TFile.Open(filename)
    h = f.FindObjectAny(histname)
    h.SetDirectory(0)
    #f.Close()

    n = h.GetNbinsX()
    Emean, Flux, Flux_err = [], [], []

    for i in range(1, n+1):
        xlow = h.GetBinLowEdge(i)
        xup  = xlow + h.GetBinWidth(i)
        x    = np.sqrt(xlow * xup)   # centro geometrico, corretto per bin log
        y    = h.GetBinContent(i)
        e    = h.GetBinError(i)
        if y <= 0:
            continue
        Emean.append(x)
        Flux.append(y) #* x**alpha)
        Flux_err.append(e) #* x**alpha)

    Emean    = np.array(Emean, dtype='d')
    Flux     = np.array(Flux,  dtype='d')
    Flux_err = np.array(Flux_err, dtype='d')
    null     = np.zeros(len(Emean), dtype='d')

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_err, Flux_err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)
    return gr


if __name__ == '__main__':

    file_DAMPE2026 = '../../TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_24sett26_wPHe_kernel_5bin_PLOT.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 20, 1.3, 2.6)

    file_DAMPE2026_SBPL = 'TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6Up_PSDprogr_STKch450_comb_vert0e7_SBPLmix_ratioGen_PLOT.dat'
    gr_DAMPE2026_SBPL = make_flux_graph_DAMPE2026(file_DAMPE2026_SBPL, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026_SBPL_EPOSLHC = 'TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6Up_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_SBPLmix_ratioGen_PLOT.dat'
    gr_DAMPE2026_SBPL_EPOSLHC = make_flux_graph_DAMPE2026(file_DAMPE2026_SBPL_EPOSLHC, kBlue+1, 20, 1.3, 2.6)

    file_DAMPE2026_EPOSLHC = '../../TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC_PLOT.dat'
    gr_DAMPE2026_EPOSLHC = make_flux_graph_DAMPE2026(file_DAMPE2026_EPOSLHC, kBlue+1, 24, 1.3, 2.6)

    filename_Geneva_sum = '../../GENEVA_pHe_FILES/Proton_plus_Helium_20260922.txt'
    gr_DAMPE2026_pHe_Geneva_sum = make_flux_graph_pHe_sum(filename_Geneva_sum, kGreen+2, 21, 1.3, 2.6)



    cc = TCanvas("cc", "Flux", 1050, 750)
    cc.SetLeftMargin(0.13)
    cc.SetRightMargin(0.04)
    cc.SetTopMargin(0.05)
    cc.SetBottomMargin(0.12)
    cc.SetTicks(1,1)
    cc.SetLogx()

    #frame = cc.DrawFrame(1e1, 1e3, 2e7, 20e3)
    frame = cc.DrawFrame(1e1, 3.5e3, 5e6, 17.e3)

    frame.GetXaxis().SetTitle("Kinetic energy [GeV]")
    frame.GetYaxis().SetTitle("E^{2.6} Flux [m^{-2} s^{-1} sr^{-1} (GeV)^{1.6}]")

    frame.GetXaxis().SetLabelSize(0.035)
    frame.GetYaxis().SetLabelSize(0.035)

    frame.GetXaxis().SetTitleSize(0.035)
    frame.GetYaxis().SetTitleSize(0.035)

    frame.GetXaxis().SetTitleOffset(1.4)
    frame.GetYaxis().SetTitleOffset(1.7)

    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    
    gr_DAMPE2026_pHe_Geneva_sum.Draw("P SAME")
    
    gr_DAMPE2026.Draw("P SAME")
    gr_DAMPE2026_SBPL.Draw("P SAME")
    gr_DAMPE2026_EPOSLHC.Draw("P SAME")
    gr_DAMPE2026_SBPL_EPOSLHC.Draw("P SAME")

    # ------------------- LEGEND

    leg = TLegend(0.17, 0.67, 0.42, 0.91)  # x1,y1,x2,y2 in NDC pad1
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.023)

    leg.AddEntry(gr_DAMPE2026_pHe_Geneva_sum, "#Phi_{p} + #Phi_{He} DAMPE (2026) ", "PE")
    leg.AddEntry(gr_DAMPE2026, "p+He DAMPE flux 50-50 ", "PE")
    leg.AddEntry(gr_DAMPE2026_SBPL, "p+He DAMPE SBPL composition model ", "PE")
    leg.AddEntry(gr_DAMPE2026_EPOSLHC, "p+He DAMPE 50-50 (10-100TeV EPOSLHC_FTFP) ", "PE")
    leg.AddEntry(gr_DAMPE2026_SBPL_EPOSLHC, "p+He DAMPE SBPL composition model (10-100TeV EPOSLHC_FTFP) ", "PE")
    leg.Draw()


    cc.Update()

    cc.SaveAs('PLOTS/flux_pHe_mix_Geneva_comparison_10TeV_EPOSLHC.pdf')
    cc.SaveAs('PLOTS/flux_pHe_mix_Geneva_comparison_10TeV_EPOSLHC.png')

    raw_input("Press enter..")



