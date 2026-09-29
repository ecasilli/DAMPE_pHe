 
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

def draw_legend_marker(graph, x, y):

    marker = TMarker(x, y, graph.GetMarkerStyle())

    marker.SetNDC()
    marker.SetMarkerColor(graph.GetMarkerColor())
    marker.SetMarkerSize(graph.GetMarkerSize())

    marker.Draw()

    return marker


if __name__ == '__main__':

    file_DAMPE2026 = 'TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_24sett26_wPHe_kernel_5bin.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 24, 1.4, 2.6)

    file_DAMPE2026_nominal = 'TXT_FILES/flux_spectrum_pHe_2026_Orb120Month_3sLow_6sUp_PSDprogr_STKch450_comb_vert0e7_nominal.dat'
    gr_DAMPE2026_nominal = make_flux_graph_DAMPE2026(file_DAMPE2026_nominal, kRed+1, 20, 1.4, 2.6)

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
    
    gr_DAMPE2026.Draw("P SAME")
    gr_DAMPE2026_nominal.Draw("P SAME")

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

    leg.AddEntry(gr_DAMPE2026, "DAMPE ", "P")
    #leg.Draw()

    cc.Update()

    #cc.SaveAs('PLOTS/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6sigmaUp_450adc.pdf')
    #cc.SaveAs('PLOTS/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6sigmaUp_450adc.png')
    #cc.SaveAs('ROOT_FILES/flux_pHe_update2026_wPHe_24sett26_3sigmaLow_6sigmaUp_450adc.root')

    raw_input("Press enter..")



