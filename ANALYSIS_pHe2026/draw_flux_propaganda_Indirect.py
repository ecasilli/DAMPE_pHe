 
import sys
from ROOT import gStyle, TGraph, TGraphErrors, TGraphAsymmErrors, TLatex, TLegend, TLine, TBox, TCanvas, gPad, kBlack, kViolet, kTeal, kGray, kRed, kBlue, kGreen, kAzure, kOrange, kMagenta, TPad, TGaxis
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
    gr_in.SetFillColor(kBlue-10)
    gr_in.SetFillStyle(1001)

    return gr, gr_in 

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

def make_flux_graph_EASTOP(color=28, marker=30, size=2.2, alpha=2.6):
    Emean = np.array([80e3], dtype='double')
    Flux_in = np.array([1.80e-09], dtype='double')
    Err_in = np.array([0.44e-09], dtype='double')

    Flux = Flux_in * Emean**alpha
    Err  = Err_in  * Emean**alpha

    Ene_err = np.array([0.01], dtype='double')

    gr = TGraphErrors(1, Emean, Flux, Ene_err, Err)
    gr.SetMarkerStyle(marker)
    gr.SetMarkerSize(size)
    gr.SetMarkerColor(color)
    gr.SetLineColor(color)

    return gr


if __name__ == '__main__':

    file_DAMPE2024 = 'TXT_FILES/DAMPE_p+He_72M26.dat'
    gr_DAMPE2024, gr_DAMPE2024_sys, gr_DAMPE2024_sys_had = make_flux_graph_DAMPE2024(file_DAMPE2024, kRed+1, 20, 1.3, 2.6)

    file_DAMPE2026 = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth_PLOT.dat'
    gr_DAMPE2026 = make_flux_graph_DAMPE2026(file_DAMPE2026, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026_all = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.dat'
    gr_DAMPE2026_all = make_flux_graph_DAMPE2026(file_DAMPE2026_all, kRed+1, 24, 1.3, 2.6)

    file_DAMPE2026_COR = 'TXT_FILES/flux_spectrum_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth_CORR_PLOT.dat'
    gr_DAMPE2026_COR = make_flux_graph_DAMPE2026(file_DAMPE2026_COR, kGreen, 24, 1.3, 2.6)
    
    file_LHAASO_EPOSLHC = 'TXT_FILES/light_component_LHAASO_EPOSLHC.dat'
    gr_LHAASO_EPOSLHC, gr_LHAASO_EPOSLHC_sys = make_flux_graph_LHAASO(file_LHAASO_EPOSLHC, kBlue+1, 25, 1.3, 2.6)

    file_LHAASO_QGSJET = 'TXT_FILES/light_component_LHAASO_QGSJET.dat'
    gr_LHAASO_QGSJET, gr_LHAASO_QGSJET_sys = make_flux_graph_LHAASO(file_LHAASO_QGSJET, kBlue+1, 24, 1.3, 2.6)

    file_LHAASO_SIBYLL = 'TXT_FILES/light_component_LHAASO_SIBYLL.dat'
    gr_LHAASO_SIBYLL, gr_LHAASO_SIBYLL_sys = make_flux_graph_LHAASO(file_LHAASO_SIBYLL, kBlue+1, 26, 1.3, 2.6)

    file_ARGO = '../pHe2024/SPECTRUM/Indirect/ARGO2New.txt'
    gr_ARGO, gr_ARGO_sys_up, gr_ARGO_sys_down = make_flux_graph_ARGO(file_ARGO, 42, 26, 1.4, 2.6)

    file_HAWC = '../pHe2024/SPECTRUM/Indirect/HAWC2022.txt'
    gr_HAWC, gr_HAWC_sys_up, gr_HAWC_sys_down = make_flux_graph_HAWC(file_HAWC, kGray+2, 32, 1.4, 2.6)

    file_KASCADE_QGSJet = '../pHe2024/SPECTRUM/Indirect/KASCADE_QGSjet01.txt'
    gr_KASCADE_QGSJet = make_flux_graph_KASCADE_QGSJet(file_KASCADE_QGSJet, kGreen+2, 27, 2.0, 2.6)

    file_KASCADE_SIBYLL = '../pHe2024/SPECTRUM/Indirect/KASCADE_SIBYLL.txt'
    gr_KASCADE_SIBYLL = make_flux_graph_KASCADE_SIBYLL(file_KASCADE_SIBYLL, kGreen+2, 28, 2.0, 2.6)

    gr_EAS_TOP = make_flux_graph_EASTOP()


    cc = TCanvas("cc", "Flux", 1050, 750)
    cc.SetLeftMargin(0.13)
    cc.SetRightMargin(0.04)
    cc.SetTopMargin(0.05)
    cc.SetBottomMargin(0.12)
    cc.SetTicks(1,1)
    cc.SetLogx()
    cc.SetLogy()

    frame = cc.DrawFrame(1e1, 5e2, 5e7, 20e3)

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
    
    gr_LHAASO_EPOSLHC_sys.Draw("E3 SAME")
    gr_DAMPE2024_sys_had.Draw("E3 SAME")
    gr_DAMPE2024_sys.Draw("E3 SAME")

    gr_EAS_TOP.Draw("P SAME")

    gr_ARGO.Draw("P SAME")
    gr_ARGO_sys_up.Draw("L SAME")
    gr_ARGO_sys_down.Draw("L SAME")

    gr_HAWC.Draw("P SAME")
    gr_HAWC_sys_up.Draw("L SAME")
    gr_HAWC_sys_down.Draw("L SAME")

    gr_KASCADE_QGSJet.Draw("P SAME")
    gr_KASCADE_SIBYLL.Draw("P SAME")

    gr_LHAASO_QGSJET.Draw("PEZ SAME")
    gr_LHAASO_EPOSLHC.Draw("P SAME")
    gr_LHAASO_SIBYLL.Draw("P SAME")

    gr_DAMPE2024.Draw("P SAME")
    #gr_DAMPE2026_COR.Draw("P SAME")
    gr_DAMPE2026.Draw("P SAME")
    gr_DAMPE2026_all.Draw("P SAME")


    # ------------------- LEGEND

    leg = TLegend(0.20,0.20,0.44,0.65)  
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.024)
    leg.SetHeader("p+He indirect measurements")
    leg.SetFillColor(0)


    leg.AddEntry(gr_EAS_TOP,"EAS-TOP+MACRO (2004)","ep")
    leg.AddEntry(gr_ARGO,"ARGO YBJ+WFCT (2015)","ep")
    leg.AddEntry(gr_HAWC,"HAWC (2022)","ep")
    leg.AddEntry(gr_KASCADE_QGSJet,"KASCADE QGSjet01 (2005)","ep")
    leg.AddEntry(gr_KASCADE_SIBYLL,"KASCADE SIBYLL-2.1 (2005)","ep")
    leg.AddEntry(gr_DAMPE2024,"DAMPE (PRL 2024)","ep")
    leg.AddEntry(gr_DAMPE2026,"DAMPE (this work 2026 - in progress)","ep")
    
    leg.AddEntry(gr_LHAASO_QGSJET, "LHAASO (QGSJET-II-04, PRL 2026)", "EP")
    leg.AddEntry(gr_LHAASO_EPOSLHC,"LHAASO (EPOS-LHC, PRL 2026)", "EP")
    leg.AddEntry(gr_LHAASO_SIBYLL, "LHAASO (SIBYLL 2.3d, PRL 2026)", "EP")
    #leg.AddEntry(gr_LHAASO_EPOSLHC_sys, "sys. error LHAASO (EPOS-LHC, PRL 2026)", "f")
    leg.Draw()


    cc.Update()

    cc.SaveAs('PLOTS/flux_pHe_update2026_cfrLHAASO_indirect_all.pdf')
    cc.SaveAs('PLOTS/flux_pHe_update2026_cfrLHAASO_indirect_all.png')

    raw_input("Press enter..")



