 
import sys
from ROOT import gStyle, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex, TLegend, TCanvas, gPad, kBlack, kGray, kRed, kBlue, kGreen, kAzure, kOrange, kMagenta
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
    EmeanP       = np.loadtxt(filenameP, skiprows=2, usecols=(1,), unpack=True)
    Flux_2P      = np.loadtxt(filenameP, skiprows=2, usecols=(4,), unpack=True)
    Flux_stat_loP= np.loadtxt(filenameP, skiprows=2, usecols=(5,), unpack=True)
    Flux_stat_upP= np.loadtxt(filenameP, skiprows=2, usecols=(6,), unpack=True)

    FluxP    = (Flux_2P) * EmeanP**alpha
    Flux_err_loP= (Flux_stat_loP) * EmeanP**alpha
    Flux_err_upP= (Flux_stat_upP) * EmeanP**alpha


    EmeanHe       = np.loadtxt(filenameHe, skiprows=2, usecols=(1,), unpack=True)
    Flux_2He      = np.loadtxt(filenameHe, skiprows=2, usecols=(4,), unpack=True)
    Flux_stat_loHe= np.loadtxt(filenameHe, skiprows=2, usecols=(5,), unpack=True)
    Flux_stat_upHe= np.loadtxt(filenameHe, skiprows=2, usecols=(6,), unpack=True)

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


if __name__ == '__main__':

    file_DAMPE2024 = 'TXT_FILES/DAMPE_p+He_72M26.dat'
    gr_DAMPE2024, gr_DAMPE2024_sys, gr_DAMPE2024_sys_had = make_flux_graph_DAMPE2024(file_DAMPE2024, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026 = 'TXT_FILES/DAMPE_p+He_120M_paperDraft.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 20, 1.33, 2.6)

    filename_GenevaP  = 'TXT_FILES/DAMPE_p_2026_pHePaperDraft.txt'
    filename_GenevaHe = 'TXT_FILES/DAMPE_He_2026_pHePaperDraft.txt'
    gr_DAMPE2026_pHe_Geneva = make_flux_graph_pHe(filename_GenevaP, filename_GenevaHe, kGreen+1, 21, 1.3, 2.6)

    file_LHAASO_EPOSLHC = 'TXT_FILES/light_component_LHAASO_EPOSLHC.dat'
    gr_LHAASO_EPOSLHC, gr_LHAASO_EPOSLHC_sys = make_flux_graph_LHAASO(file_LHAASO_EPOSLHC, kBlue+1, 25, 1.3, 2.6)

    file_LHAASO_QGSJET = 'TXT_FILES/light_component_LHAASO_QGSJET.dat'
    gr_LHAASO_QGSJET, gr_LHAASO_QGSJET_sys = make_flux_graph_LHAASO(file_LHAASO_QGSJET, kBlue+1, 24, 1.3, 2.6)

    file_LHAASO_SIBYLL = 'TXT_FILES/light_component_LHAASO_SIBYLL.dat'
    gr_LHAASO_SIBYLL, gr_LHAASO_SIBYLL_sys = make_flux_graph_LHAASO(file_LHAASO_SIBYLL, kBlue+1, 26, 1.3, 2.6)

    cc = TCanvas("", "", 1200, 800)

    cc.SetTopMargin(0.02)
    cc.SetRightMargin(0.04)
    cc.SetBottomMargin(0.13)
    cc.SetLeftMargin(0.13)

    cc.SetTicks(1,1)
    cc.SetLogx()

    frame = cc.DrawFrame(1e1, 1e3, 2e7, 21.5e3)

    frame.GetXaxis().SetTitle("Kinetic energy (GeV)")
    frame.GetYaxis().SetTitle("Flux #times E^{2.6} (m^{-2} sr^{-1} s^{-1} GeV^{1.6})")

    frame.GetXaxis().SetLabelSize(0.045)
    frame.GetXaxis().SetTitleSize(0.050)
    frame.GetXaxis().SetTitleOffset(1.29)

    frame.GetYaxis().SetLabelSize(0.045)
    frame.GetYaxis().SetTitleSize(0.050)
    frame.GetYaxis().SetTitleOffset(1.39)

    frame.GetXaxis().CenterTitle()
    frame.GetYaxis().CenterTitle()
    
    gr_DAMPE2024_sys_had.Draw("E3 SAME")
    gr_DAMPE2024_sys.Draw("E3 SAME")
    gr_LHAASO_EPOSLHC_sys.Draw("E3 SAME")
    gr_DAMPE2024.Draw("P SAME")

    gr_DAMPE2026_pHe_Geneva.Draw("P SAME")
    
    gr_DAMPE2026.Draw("P SAME")
    
    gr_LHAASO_QGSJET.Draw("PEZ SAME")
    gr_LHAASO_EPOSLHC.Draw("P SAME")
    gr_LHAASO_SIBYLL.Draw("P SAME")

    label = TLatex()
    label.SetNDC()
    label.SetTextFont(62)  
    label.SetTextSize(0.045)
    label.SetTextAlign(31)  
    label.DrawLatex(0.9, 0.9, "p+He")

    # ------------------- LEGEND

    leg = TLegend(0.17, 0.69, 0.42, 0.9)  
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.023)
    #leg.SetHeader("p+He ");
    #leg.SetNColumns(2)
    leg.AddEntry(gr_DAMPE2024, "p+He DAMPE (PRL 2024)", "PE")
    leg.AddEntry(gr_DAMPE2024_sys,"ana. error (PRL 2024)","f")
    leg.AddEntry(gr_DAMPE2024_sys_had,"ana. #oplus had. error (PRL 2024)","f")

    leg.AddEntry(gr_DAMPE2026_pHe_Geneva, "#Phi_{p} + #Phi_{He} DAMPE (2026) ", "PE")
    leg.AddEntry(gr_DAMPE2026, "p+He DAMPE (this work - preliminary) ", "PE")
    leg.Draw()

    
    leg1 = TLegend(0.5, 0.7, 0.73, 0.87)  
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextSize(0.023)


    leg1.AddEntry(gr_LHAASO_QGSJET, "light comp. LHAASO (QGSJET-II-04, PRL 2026)", "PE")
    leg1.AddEntry(gr_LHAASO_EPOSLHC, "light comp. LHAASO (EPOS-LHC, PRL 2026)", "PE")
    leg1.AddEntry(gr_LHAASO_SIBYLL, "light comp. LHAASO (SIBYLL 2.3d, PRL 2026)", "PE")
    leg1.AddEntry(gr_LHAASO_EPOSLHC_sys, "sys. error LHAASO (EPOS-LHC, PRL 2026)", "f")
    leg1.Draw()


    cc.Update()

    cc.SaveAs('PLOTS/flux_pHe_update2026_wPHe_5bins.pdf')

    raw_input("Press enter..")



