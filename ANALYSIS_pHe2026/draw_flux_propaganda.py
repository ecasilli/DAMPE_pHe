 
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

    '''
    for i in range(len(Emean)):
        gr_in.SetPoint(i+1, Emean[i], Flux[i])
        gr_in.SetPointError(i+1, 0, 0, syst[i], syst[i])
        gr_out.SetPoint(i+1, Emean[i], Flux[i])
        gr_out.SetPointError(i+1, 0, 0, syst_had[i], syst_had[i])
    gr_in.SetPoint(len(Emean)+1, Eup[-1], Flux[-1])
    gr_in.SetPointError(len(Emean)+1, 0, 0, syst[-1], syst[-1])
    gr_out.SetPoint(len(Emean)+1, Eup[-1], Flux[-1])
    gr_out.SetPointError(len(Emean)+1, 0, 0, syst_had[-1], syst_had[-1])
    '''

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

    gr_in  = TGraphErrors(len(Emean), Emean, Flux, null, Flux_sys) #TGraphAsymmErrors(len(Emean)) #+2)
    gr_in.SetLineColor(0)
    gr_in.SetMarkerColor(color)
    gr_in.SetMarkerStyle(marker)
    gr_in.SetMarkerSize(size)
    gr_in.SetFillColor(38)
    gr_in.SetFillStyle(3002)

    return gr, gr_in 


if __name__ == '__main__':

    file_DAMPE2024 = 'TXT_FILES/DAMPE_p+He_72M26.dat'
    gr_DAMPE2024, gr_DAMPE2024_sys, gr_DAMPE2024_sys_had = make_flux_graph_DAMPE2024(file_DAMPE2024, kRed+1, 20, 1.3, 2.6)

    file_DAMPE2026 = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth_PLOT.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026_all = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.dat'
    gr_DAMPE2026_all = make_flux_graph_DAMPE2026(file_DAMPE2026_all, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026_first6years = 'TXT_FILES/flux_spectrum_pHe_2026_72months_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.dat'
    gr_DAMPE2026_first6years = make_flux_graph_DAMPE2026(file_DAMPE2026_first6years, kOrange+1, 20, 1.3, 2.6)

    file_DAMPE2026_last4years = 'TXT_FILES/flux_spectrum_pHe_2026_rem48months_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.dat'
    gr_DAMPE2026_last4years = make_flux_graph_DAMPE2026(file_DAMPE2026_last4years, kGreen+1, 20, 1.3, 2.6)

    file_DAMPE2026_COR = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth_CORR_PLOT.dat'
    gr_DAMPE2026_COR = make_flux_graph_DAMPE2026(file_DAMPE2026_COR, kGreen, 24, 1.3, 2.6)
    
    file_LHAASO_EPOSLHC = 'TXT_FILES/light_component_LHAASO_EPOSLHC.dat'
    gr_LHAASO_EPOSLHC, gr_LHAASO_EPOSLHC_sys = make_flux_graph_LHAASO(file_LHAASO_EPOSLHC, kBlue+1, 25, 1.3, 2.6)

    file_LHAASO_QGSJET = 'TXT_FILES/light_component_LHAASO_QGSJET.dat'
    gr_LHAASO_QGSJET, gr_LHAASO_QGSJET_sys = make_flux_graph_LHAASO(file_LHAASO_QGSJET, kBlue+1, 24, 1.3, 2.6)

    file_LHAASO_SIBYLL = 'TXT_FILES/light_component_LHAASO_SIBYLL.dat'
    gr_LHAASO_SIBYLL, gr_LHAASO_SIBYLL_sys = make_flux_graph_LHAASO(file_LHAASO_SIBYLL, kBlue+1, 26, 1.3, 2.6)

    cc = TCanvas("cc", "Flux", 1050, 750)
    cc.SetLeftMargin(0.13)
    cc.SetRightMargin(0.04)
    cc.SetTopMargin(0.05)
    cc.SetBottomMargin(0.12)
    cc.SetTicks(1,1)
    cc.SetLogx()

    frame = cc.DrawFrame(1e1, 1e3, 2e7, 20e3)

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
    
    gr_DAMPE2024_sys_had.Draw("E3 SAME")
    gr_DAMPE2024_sys.Draw("E3 SAME")
    #gr_LHAASO_EPOSLHC_sys.Draw("E3 SAME")
    gr_DAMPE2024.Draw("P SAME")
    #gr_DAMPE2026_COR.Draw("P SAME")
    gr_DAMPE2026.Draw("P SAME")
    gr_DAMPE2026_first6years.Draw("P SAME")
    gr_DAMPE2026_last4years.Draw("P SAME")
    gr_DAMPE2026_all.Draw("P SAME")
    #gr_LHAASO_QGSJET.Draw("PEZ SAME")
    #gr_LHAASO_EPOSLHC.Draw("P SAME")
    #gr_LHAASO_SIBYLL.Draw("P SAME")


    # ------------------- LEGEND

    leg = TLegend(0.17, 0.72, 0.42, 0.9)  # x1,y1,x2,y2 in NDC pad1
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.023)
    #leg.SetHeader("p+He ");
    #leg.SetNColumns(2)
    leg.AddEntry(gr_DAMPE2024, "p+He DAMPE (PRL 2024)", "P")
    leg.AddEntry(gr_DAMPE2024_sys,"ana. error (PRL 2024)","f")
    leg.AddEntry(gr_DAMPE2024_sys_had,"ana. #oplus had. error (PRL 2024)","f")
    leg.AddEntry(gr_DAMPE2026, "p+He DAMPE (this work 2026 - in progress) ", "P")
    leg.Draw()

    '''
    leg1 = TLegend(0.55, 0.72, 0.78, 0.9)  # x1,y1,x2,y2 in NDC pad1
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextSize(0.023)


    leg1.AddEntry(gr_LHAASO_QGSJET, "light comp. LHAASO (QGSJET-II-04, PRL 2026)", "P")
    leg1.AddEntry(gr_LHAASO_EPOSLHC, "light comp. LHAASO (EPOS-LHC, PRL 2026)", "P")
    leg1.AddEntry(gr_LHAASO_SIBYLL, "light comp. LHAASO (SIBYLL 2.3d, PRL 2026)", "P")
    leg1.AddEntry(gr_LHAASO_EPOSLHC_sys, "sys. error LHAASO (EPOS-LHC, PRL 2026)", "f")
    leg1.Draw()
    '''

    leg1 = TLegend(0.55, 0.79, 0.78, 0.9)  # x1,y1,x2,y2 in NDC pad1
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextSize(0.023)
    leg1.AddEntry(gr_DAMPE2026_first6years, "p+He DAMPE (this work 2016 - 2021)", "P")
    leg1.AddEntry(gr_DAMPE2026_last4years, "p+He DAMPE (this work 2022 - 2025)", "P")
    leg1.Draw()


    cc.Update()

    cc.SaveAs('PLOTS/flux_pHe_update2026_all_cfr6-4years.pdf')
    cc.SaveAs('PLOTS/flux_pHe_update2026_all_cfr6-4years.png')

    raw_input("Press enter..")



