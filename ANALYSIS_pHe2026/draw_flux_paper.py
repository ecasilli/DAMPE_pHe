 
import sys
from ROOT import gStyle, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex, TLegend, TMarker, TCanvas, gPad, kBlack, kGray, kRed, kBlue, kGreen, kAzure, kOrange, kMagenta, TPad, TGaxis
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

    gr_in  = TGraphErrors(len(Emean), Emean, Flux, null, syst) 
    gr_out = TGraphErrors(len(Emean), Emean, Flux, null, syst_had) 
    gr_in.SetLineColor(0)
    gr_in.SetMarkerColor(color)
    gr_in.SetMarkerStyle(marker)
    gr_in.SetMarkerSize(size)
    gr_in.SetFillColor(17)
    gr_in.SetFillStyle(1001)
    
    gr_out.SetLineColor(0)
    gr_out.SetMarkerColor(color)
    gr_out.SetMarkerStyle(marker)
    gr_out.SetMarkerSize(size)
    gr_out.SetFillColor(18)
    gr_out.SetFillStyle(1001)

    return gr, gr_in, gr_out

def make_flux_graph_DAMPE2026(filename, color, marker, size, alpha):
    Emean    = np.loadtxt(filename, skiprows=1, usecols=(0,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)
    Flux_stat= np.loadtxt(filename, skiprows=1, usecols=(5,), unpack=True)

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

def make_flux_graph_LHAASO(filename, color, marker, size, alpha):
    logEmin    = np.loadtxt(filename, skiprows=1, usecols=(0,), unpack=True)
    logEmax    = np.loadtxt(filename, skiprows=1, usecols=(1,), unpack=True)
    Flux_2   = np.loadtxt(filename, skiprows=1, usecols=(2,), unpack=True)
    Flux_stat= np.loadtxt(filename, skiprows=1, usecols=(3,), unpack=True)
    Sys_err  = np.loadtxt(filename, skiprows=1, usecols=(4,), unpack=True)

    
    Emean = 10**( (logEmin+logEmax)/2. + 6 )
    Flux    = (Flux_2)    * Emean**alpha * 1e-6
    Flux_err= (Flux_stat) * Emean**alpha * 1e-6
    Flux_sys= (Sys_err) * Emean**alpha * 1e-6

    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_err, Flux_err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    gr_in  = TGraphErrors(len(Emean), Emean, Flux, null, Flux_sys) 
    gr_in.SetLineColor(0)
    gr_in.SetMarkerColor(color)
    gr_in.SetMarkerStyle(marker)
    gr_in.SetMarkerSize(size)
    gr_in.SetFillColor(kBlue-10)
    gr_in.SetFillStyle(1001)

    return gr, gr_in 

def make_flux_graph_KASCADE_SIBYLL(filename, color, marker, size, alpha):
    Emean = np.loadtxt(filename, usecols=(0,), unpack=True)
    FluxP = np.loadtxt(filename, usecols=(1,), unpack=True)
    ErrUpP   = np.loadtxt(filename, usecols=(2,), unpack=True)
    ErrDownP = np.loadtxt(filename, usecols=(3,), unpack=True)
    FluxHe = np.loadtxt(filename, usecols=(4,), unpack=True)
    ErrUpHe   = np.loadtxt(filename, usecols=(5,), unpack=True)
    ErrDownHe = np.loadtxt(filename, usecols=(6,), unpack=True)

    Flux = (FluxP + FluxHe) * Emean**alpha
    ErrUp = np.sqrt(ErrUpP**2 + ErrUpHe**2) * Emean**alpha
    ErrDown = np.sqrt(ErrDownP**2 + ErrDownHe**2) * Emean**alpha
    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, ErrDown, ErrUp)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def make_flux_graph_ARGO(filename, color, marker, size, alpha):
    logE_TeV      = np.loadtxt(filename, usecols=(0,), unpack=True)
    Flux_in       = np.loadtxt(filename, usecols=(1,), unpack=True)
    Stat_up       = np.loadtxt(filename, usecols=(2,), unpack=True)
    Stat_down     = np.loadtxt(filename, usecols=(3,), unpack=True)
    Sys_up_perc   = np.loadtxt(filename, usecols=(4,), unpack=True)
    Sys_down_perc = np.loadtxt(filename, usecols=(5,), unpack=True)
    Cont          = np.loadtxt(filename, usecols=(6,), unpack=True)

    Emean = (10**logE_TeV) * 1e3
    Flux = (Flux_in / 1e12) * Emean**alpha
    Stat_up_E   = (Stat_up   / 1e12) * Emean**alpha
    Stat_down_E = (Stat_down / 1e12) * Emean**alpha
    Sys_up = Flux + (Sys_up_perc * Flux)/100.0
    Sys_down = ( Flux - (Sys_down_perc * Flux)/100.0 - (Cont * Flux)/100.0 )

    exl = np.full(len(Emean), 0.01)
    exh = np.full(len(Emean), 0.01)

    gr_stat = TGraphAsymmErrors(len(Emean), Emean, Flux, exl, exh, Stat_down_E, Stat_up_E)
    gr_sys_up = TGraph(len(Emean), Emean, Sys_up)
    gr_sys_down = TGraph(len(Emean), Emean, Sys_down)

    for gr in [gr_stat, gr_sys_up, gr_sys_down]:
        gr.SetMarkerStyle(marker)
        gr.SetMarkerSize(size)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)

    return gr_stat, gr_sys_up, gr_sys_down


def make_flux_graph_HAWC(filename, color, marker, size, alpha):
    Emean       = np.loadtxt(filename, usecols=(0,), unpack=True)
    Flux        = np.loadtxt(filename, usecols=(1,), unpack=True)
    StatErr     = np.loadtxt(filename, usecols=(2,), unpack=True)
    Sys_up_in   = np.loadtxt(filename, usecols=(3,), unpack=True)
    Sys_down_in = np.loadtxt(filename, usecols=(4,), unpack=True)

    Flux_E = Flux * Emean**alpha
    Stat_E = StatErr * Emean**alpha
    Sys_up   = Flux_E + (Sys_up_in * Emean**alpha)
    Sys_down = Flux_E - (Sys_down_in * Emean**alpha)
    Ene_err = Emean * 1e-5

    gr_stat = TGraphErrors( len(Emean), Emean, Flux_E, Ene_err, Stat_E)
    gr_sys_up = TGraph(len(Emean), Emean, Sys_up)
    gr_sys_down = TGraph(len(Emean), Emean, Sys_down)

    for gr in [gr_stat, gr_sys_up, gr_sys_down]:
        gr.SetMarkerStyle(marker)
        gr.SetMarkerSize(size)
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
    gr_sys_up.SetLineStyle(9)
    gr_sys_down.SetLineStyle(9)

    return gr_stat, gr_sys_up, gr_sys_down


def make_flux_graph_KASCADE_QGSJet(filename, color, marker, size, alpha):
    Emean = np.loadtxt(filename, usecols=(0,), unpack=True)
    FluxP = np.loadtxt(filename, usecols=(1,), unpack=True)
    ErrP  = np.loadtxt(filename, usecols=(3,), unpack=True)
    FluxHe = np.loadtxt(filename, usecols=(4,), unpack=True)
    ErrHe  = np.loadtxt(filename, usecols=(6,), unpack=True)

    Flux = (FluxP + FluxHe) * Emean**alpha
    Err = np.sqrt(ErrP**2 + ErrHe**2) * Emean**alpha
    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Err, Err)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def make_flux_graph_KASCADE_SIBYLL(filename, color, marker, size, alpha):
    Emean = np.loadtxt(filename, usecols=(0,), unpack=True)
    FluxP = np.loadtxt(filename, usecols=(1,), unpack=True)
    ErrUpP   = np.loadtxt(filename, usecols=(2,), unpack=True)
    ErrDownP = np.loadtxt(filename, usecols=(3,), unpack=True)
    FluxHe = np.loadtxt(filename, usecols=(4,), unpack=True)
    ErrUpHe   = np.loadtxt(filename, usecols=(5,), unpack=True)
    ErrDownHe = np.loadtxt(filename, usecols=(6,), unpack=True)

    Flux = (FluxP + FluxHe) * Emean**alpha
    ErrUp = np.sqrt(ErrUpP**2 + ErrUpHe**2) * Emean**alpha
    ErrDown = np.sqrt(ErrDownP**2 + ErrDownHe**2) * Emean**alpha
    null = np.zeros(len(Emean))

    gr = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, ErrDown, ErrUp)
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)

    return gr

def draw_legend_marker(graph, x, y):

    marker = TMarker(x, y, graph.GetMarkerStyle())

    marker.SetNDC()
    marker.SetMarkerColor(graph.GetMarkerColor())
    marker.SetMarkerSize(graph.GetMarkerSize())

    marker.Draw()

    return marker


if __name__ == '__main__':

    file_DAMPE2024 = 'TXT_FILES/DAMPE_p+He_72M26.dat'
    gr_DAMPE2024, gr_DAMPE2024_sys, gr_DAMPE2024_sys_had = make_flux_graph_DAMPE2024(file_DAMPE2024, kRed+1, 24, 1.3, 2.6)

    #file_DAMPE2026 = 'TXT_FILES/DAMPE_p+He_120M_paperDraft.dat'
    file_DAMPE2026 = 'TXT_FILES/DAMPE_p+He_120M_paperDraft_24sett26_3sigmaLow_6e5sigmaUp_450adc.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 20, 1.4, 2.6)

    filename_GenevaP  = 'TXT_FILES/PROTON_SEP2026_ANDRII_flux_noescale_ekin_p_2026.txt'
    filename_GenevaHe = 'TXT_FILES/HELIUM_SEP2026_PAUL_Helium_Paul_Geant4.txt'
    gr_DAMPE2026_pHe_Geneva = make_flux_graph_pHe(filename_GenevaP, filename_GenevaHe, kGreen+2, 21, 1.3, 2.6)

    filename_Geneva_sum = 'GENEVA_pHe_FILES/Proton_plus_Helium_20260922.txt'
    gr_DAMPE2026_pHe_Geneva_sum = make_flux_graph_pHe_sum(filename_Geneva_sum, kGreen+2, 21, 1.3, 2.6)

    file_LHAASO_EPOSLHC = 'TXT_FILES/light_component_LHAASO_EPOSLHC.dat'
    gr_LHAASO_EPOSLHC, gr_LHAASO_EPOSLHC_sys = make_flux_graph_LHAASO(file_LHAASO_EPOSLHC, kGray+2, 22, 1.3, 2.6)

    file_LHAASO_QGSJET = 'TXT_FILES/light_component_LHAASO_QGSJET.dat'
    gr_LHAASO_QGSJET, gr_LHAASO_QGSJET_sys = make_flux_graph_LHAASO(file_LHAASO_QGSJET, kGray+2, 26, 1.3, 2.6)

    file_LHAASO_SIBYLL = 'TXT_FILES/light_component_LHAASO_SIBYLL.dat'
    gr_LHAASO_SIBYLL, gr_LHAASO_SIBYLL_sys = make_flux_graph_LHAASO(file_LHAASO_SIBYLL, kGray+2, 26, 1.3, 2.6)

    file_KASCADE_QGSJet = '../pHe2024/SPECTRUM/Indirect/KASCADE_QGSjet01.txt'
    gr_KASCADE_QGSJet = make_flux_graph_KASCADE_QGSJet(file_KASCADE_QGSJet, kGray+1, 25, 1.3, 2.6)

    file_KASCADE_SIBYLL = '../pHe2024/SPECTRUM/Indirect/KASCADE_SIBYLL.txt'
    gr_KASCADE_SIBYLL = make_flux_graph_KASCADE_SIBYLL(file_KASCADE_SIBYLL, kGray+1, 21, 1.3, 2.6)

    file_ARGO = '../pHe2024/SPECTRUM/Indirect/ARGO2New.txt'
    gr_ARGO, gr_ARGO_sys_up, gr_ARGO_sys_down = make_flux_graph_ARGO(file_ARGO, kGray+1, 27, 1.4, 2.6)

    file_HAWC = '../pHe2024/SPECTRUM/Indirect/HAWC2022.txt'
    gr_HAWC, gr_HAWC_sys_up, gr_HAWC_sys_down = make_flux_graph_HAWC(file_HAWC, kGray+1, 30, 1.4, 2.6)

    cc = TCanvas("cc", "Flux", 1200, 800)

    cc.SetTopMargin(0.02)
    cc.SetRightMargin(0.04)
    cc.SetBottomMargin(0.13)
    cc.SetLeftMargin(0.13)

    cc.SetTicks(1,1)
    cc.SetLogx()
    cc.SetLogy()

    frame = cc.DrawFrame(1e1, 1.5e3, 5e7, 20e3)

    frame.GetXaxis().SetTitle("Kinetic energy (GeV)")
    frame.GetYaxis().SetTitle("Flux #times E^{2.6} (m^{-2} sr^{-1} s^{-1} GeV^{1.6})")

    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetTitleSize(0.050)
    frame.GetXaxis().SetTitleOffset(1.29)

    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.29)

    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    
    gr_DAMPE2024_sys_had.Draw("E3 SAME")
    gr_DAMPE2024_sys.Draw("E3 SAME")
    #gr_LHAASO_EPOSLHC_sys.Draw("E3 SAME")
    gr_DAMPE2024.Draw("P SAME")

    gr_DAMPE2026_pHe_Geneva_sum.Draw("P SAME")
    
    gr_DAMPE2026.Draw("P SAME")
    
    gr_LHAASO_QGSJET.Draw("P SAME")
    gr_LHAASO_EPOSLHC.Draw("P SAME")
    #gr_LHAASO_SIBYLL.Draw("P SAME")

    gr_ARGO.Draw("P SAME")
    gr_ARGO_sys_up.Draw("L SAME")
    gr_ARGO_sys_down.Draw("L SAME")

    gr_HAWC.Draw("P SAME")
    gr_HAWC_sys_up.Draw("L SAME")
    gr_HAWC_sys_down.Draw("L SAME")

    gr_KASCADE_QGSJet.Draw("P SAME")
    gr_KASCADE_SIBYLL.Draw("P SAME")

    label = TLatex()
    label.SetNDC()
    label.SetTextFont(62)  
    label.SetTextSize(0.045)
    label.SetTextAlign(31)  
    label.DrawLatex(0.88, 0.88, "p+He")

    # ------------------- LEGEND

    #leg = TLegend(0.17, 0.69, 0.42, 0.9)  
    leg = TLegend(0.20,0.19,0.43,0.66)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.032)
    #leg.SetHeader("p+He ");
    #leg.SetNColumns(2)

    leg.AddEntry(gr_DAMPE2024, "DAMPE (2024)", "P")
    leg.AddEntry(gr_DAMPE2024_sys,"ana. error (2024)","f")
    leg.AddEntry(gr_DAMPE2024_sys_had,"ana. #oplus had. error (2024)","f")

    leg.AddEntry(gr_DAMPE2026_pHe_Geneva_sum, "DAMPE #Phi_{p} + #Phi_{He}", "P")
    leg.AddEntry(gr_DAMPE2026, "DAMPE ", "P")
    #leg.Draw()

    
    #leg1 = TLegend(0.5, 0.7, 0.73, 0.87)  
    #leg1.SetBorderSize(0)
    #leg1.SetFillStyle(0)
    #leg1.SetTextSize(0.023)


    #leg.AddEntry(gr_LHAASO_QGSJET, "LHAASO (QGSJET-II-04)", "P")
    #leg.AddEntry(gr_LHAASO_EPOSLHC, "LHAASO (EPOS-LHC)", "PE")
    #leg1.AddEntry(gr_LHAASO_SIBYLL, "LHAASO (SIBYLL 2.3d)", "PE")
    #leg1.AddEntry(gr_LHAASO_EPOSLHC_sys, "sys. error LHAASO (EPOS-LHC, PRL 2026)", "f")
    #leg1.Draw()

    leg.AddEntry(gr_ARGO,"ARGO (YBJ+WFCT)","p")
    leg.AddEntry(gr_HAWC,"HAWC (QGSJET-II-04)","p")

    leg.AddEntry(gr_KASCADE_QGSJet,"KASCADE (QGSJET-01), (SIBYLL-2.1)","P")
    leg.AddEntry(gr_LHAASO_QGSJET,"LHAASO (QGSJET-II-04), (EPOS-LHC)","P")
    leg.Draw()

    marker_kascade = draw_legend_marker(gr_KASCADE_SIBYLL, 0.615, 0.268)
    marker_lhaaso  = draw_legend_marker(gr_LHAASO_EPOSLHC, 0.615, 0.217)


    cc.Update()

    cc.SaveAs('PLOTS/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6e5sigmaUp_450adc.pdf')
    cc.SaveAs('PLOTS/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6e5sigmaUp_450adc.png')

    cc.SaveAs('ROOT_FILES/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6e5sigmaUp_450adc.root')

    raw_input("Press enter..")



