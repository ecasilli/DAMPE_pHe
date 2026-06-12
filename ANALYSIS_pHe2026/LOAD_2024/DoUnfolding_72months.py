# DoUnfolding.py performs the unfolding of the proton spectrum
# input root files provide by Antonio Surdo
# author M. Nicola Mazziotta (mazziotta@ba.infn.it), Nov 2016
# Plese cite: arXiv:0912.1236
#             arXiv:0912.3695
#             NIM A362 (1995) 487 

import sys
from ROOT import *
import array as ary
import numpy as np
import math

#from FastBayes import *
from FastBayes_Iteration5 import *


if __name__ == '__main__':

    nrbt = 8 # to rebin TRUE energy axis
    nrbo = 8 # to rebin observed energy axis
    #nrbt = 12 # to rebin TRUE energy axis
    #nrbo = 12 # to rebin observed energy axis

    days = 2183 # days of data collection, 72 months - 9 days of bad data acquisition 
    #days = 1818 # days of data collection, 60 months - 9 days of bad data acquisition
    #days = 1339.# days of data collection, 44 months.
    #days = 1178.# days of data collection, 39 months.
    alpha = 2.7
    DataFileOutput='DAMPE_p+He_120M.dat'
    ResultNameLegend='DAMPE_p+He_This Work'



    # fsm is the root file with MC used for the IRF 

    #fsm = "PHe_MC_p-1PeV_He-200TeV_HeSigma6_flux5bins.root"
#    fsm = "PHe_MC_p-1PeV_He-200TeV_HeSigma6_44Months.root"
    fsm = "PHe_MC_p_He-1PeV_v6r0p10p12p13.root"
    fsm_h1="PHe_MC_p_He-1PeV_h1Ngen_48binsPerDecade.root"

    #Orb_file= "PHe_skim_Orb39Month_He6_5bin.root"
    #Orb_file= "PHe_skim_Orb39Month_He6.root"
    Orb_file= "PHe_skim_Orb120Month_He6_48binperdecadeTEST.root"
    BGO_data = "h1SelBGO_orb6"
    frsm = TFile(fsm)
    frsm_h1 = TFile(fsm_h1)
    

    # get MC generated distribution, Ng(Etrue), i.e. in each bin there is the number of MC events used to evluate the detector response function
    # e.g. in case of input spectrum 1/E, in each bin the number of events is N_job*Evt_job*(log10(E_i+1) - log10(E_i))/(log10(E_max) - log10(E_min))
    # where: N_job=number of simulated jobs, Evt_job=number of true events for each job, E_min=min energy used for the MC, E_max=max energy used in MC
    # E_i+1 and E_i are the upper and lower values of the i-energy bin   
    h1ng = frsm_h1.Get("h1p1")

    

    #  MC acceptance used for the simulation in unit of m^2 sr
    Amc =  2.*TMath.Pi()*6. # to be checked if the simulation has been performed over 4 pi then a factor 2 should be added !
   # Amc = 2.*TMath.Pi()*TMath.Pi()
    # get MC migration matrix for the triggered/detected events (i.e. events after cuts), as a function of true MC energy and observed (reconstructed) energy
    # Nd(Etrue | Eobs)
    #h2nd0 = frsm.Get("h2NtrigCor")
    h2nd0 = frsm.Get("h2NtrigCor6") 

    h1ng.Rebin(nrbt) # rebin true energy axis
    h2nd0.RebinX(nrbt) # rebin true energy axis
    h2nd0.RebinY(nrbo) # rebin observed energy axis

    nbt = h2nd0.GetNbinsX()
    nbo = h2nd0.GetNbinsY()
    print "Number of true energy bin = ", nbt
    print "Number of observed energy bin = ", nbo

    Et = ary.array("d",range(nbt+1))
    Eo = ary.array("d",range(nbo+1))
    for i in range(nbt+1):
        Et[i] = h2nd0.GetXaxis().GetBinLowEdge(i+1)
    print Et
    for j in range(nbo+1):
        Eo[j] = h2nd0.GetYaxis().GetBinLowEdge(j+1)
    print Eo
    
    h2nd = TH2D("h2nd","Nd(Eo | Et) - Migration matrix; MC true energy (GeV); Reconstruced energy (GeV)", nbt, Et, nbo, Eo)
    h2sm = TH2D("h2sm","P(Eo | Et) - Smearing matrix; MC true energy (GeV); Reconstruced energy (GeV)", nbt, Et, nbo, Eo)
    h1ef = TH1D("h1ef","Detector efficiency; MC true energy (GeV); Efficiency", nbt, Et)
    h1ac = TH1D("h1ac","Detector Acceptance; MC true energy (GeV); Acceptance (m^{2} sr)", nbt, Et)
    h1ef.Sumw2()
    h1ac.Sumw2()
    for i in range(nbt):
        cg = h1ng.GetBinContent(i+1)  # if the bins of h1ng and h2nd are the same, please use this one!. 
	#cg = h1ng.GetBinContent(i+1+4) # the bin of h1ng starts at 0.1,meanwhile, h2nd starts at 1.
        cdt = 0.
        for j in range(nbo):
            cd = h2nd0.GetBinContent(i+1, j+1)
	    if(j<5.):
            	h2nd.SetBinContent(i+1, j+1, cd*1.0)	#1.21 = HET eff; 0.95 = Track eff #09FEB2018_NEWCUT02_HET*TRACK=1.1685
		if(cg>0.):
			h2sm.SetBinContent(i+1, j+1, (cd*1.0)/cg)
		cdt += cd*1.0
	    if(j>4.):
	   	h2nd.SetBinContent(i+1, j+1, cd*1.0)
		if(cg>0.):
			h2sm.SetBinContent(i+1, j+1, (cd*1.0)/cg)
		cdt += cd*1.0

        eff = 0.
        eef = 0.
        if(cg>0.):
            eff = cdt/cg
            # eff should be <= 1, if not there is something wrong
            if(eff<=1.):
                eef = TMath.Sqrt(eff*(1.-eff)/cg)
            else:
                print "Warning efficiency > 1 in the bin ", i 
        # print i, cdt, cg, eff
        h1ef.SetBinContent(i+1, eff)
        h1ef.SetBinError(i+1, eef)
        acc = eff*Amc
        eac = eef*Amc
        h1ac.SetBinContent(i+1, acc)
        h1ac.SetBinError(i+1, eac)


    #sys.exit()

    # Evaluate live time and acceptance
    
    #days = 16.
    #days = 2.

    #TotTime = days*86400*0.7667 #days*seconds per day*18.4 hours per day (which is the time left when I exclude dead time etc..)
    TotTime = 242576599.4 # 10 years
    #TotTime = 143806419.9 # 6 years
    #TotTime = 143747867 #David - 132485.76 (2 giorni * 86400 * 0.7667)
    print "Live Time = ",str(TotTime)," seconds"
    TotTime = TotTime*Amc
    print "Live Time * A_MC = ",str(TotTime)," [cm^2 s sr]"

    h1cts0 = TH1F("h1cts0","h1cts0", nbo, Eo)

    h1cts = h1cts0.Clone("h1cts")
    
    h1cts.SetNameTitle("h1cts", "Observed counts; Reconstructed Energy (GeV); Number of Events")

    h1cts.GetXaxis().SetLabelSize(0.04)
    h1cts.GetYaxis().SetLabelSize(0.04)
    h1cts.GetXaxis().SetTitleSize(0.04)
    h1cts.GetYaxis().SetTitleSize(0.04)

    h1cts_ams02 = h1cts.Clone("h1cts_ams02")
    h1cts_ams02.SetNameTitle("h1cts_ams02", "Observed counts; Reconstructed Energy (GeV); Number of Events")
    h1cts_ams02.SetLineColor(2)
    h1cts_ams02.SetMarkerColor(2)
    h1cts_ams02.SetMarkerStyle(20)

#'''
#
#    data = 'AMS02_2015_Proton_prl.dat'  #Proton
#    Emin,Emax = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
#    Emean = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
#    Flux  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
#    Flux0  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
#    Flux_low  = np.loadtxt(data,skiprows=2,usecols=(8,),unpack=True)
#    Flux_up   = np.loadtxt(data,skiprows=2,usecols=(9,),unpack=True)
#'''
    
     
    data = 'AMS02_2015_Helium_prl.dat'  #Helium
    Emin,Emax = np.loadtxt(data,skiprows=2,usecols=(1,2),unpack=True)
    Emean = np.loadtxt(data,skiprows=2,usecols=(0,),unpack=True)
    Flux  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
    Flux0  = np.loadtxt(data,skiprows=2,usecols=(3,),unpack=True)
    Flux_low  = np.loadtxt(data,skiprows=2,usecols=(8,),unpack=True)
    Flux_up   = np.loadtxt(data,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(Emean))
    grprams02prl = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_low, Flux_up)
    grprams02prl.SetLineColor(2)
    grprams02prl.SetMarkerColor(2)
    grprams02prl.SetMarkerStyle(20)

    Emin0 = ary.array("f",range(len(Emean)))
    Emax0 = ary.array("f",range(len(Emean)))
    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(Emean)):
        e0 = Emin[i]
        e1 = Emax[i]
        Emin0[i] = e0
        Emax0[i] = e1
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        Emean[i] = em
        #Flux[i] *= TMath.Power(Emean[i], alpha)
        #Flux_low[i] *= TMath.Power(Emean[i], alpha)
        #Flux_up[i] *= TMath.Power(Emean[i], alpha)
        Flux[i] *= TMath.Power(em, alpha)
        Flux_low[i] *= TMath.Power(em, alpha)
        Flux_up[i] *= TMath.Power(em, alpha)

    grprams02pow = TGraphAsymmErrors(len(Emean), Emean, Flux, null, null, Flux_low, Flux_up)
    grprams02pow.SetLineColor(2)
    grprams02pow.SetMarkerColor(2)
    grprams02pow.SetMarkerStyle(20)

    Fluxams = TH1F("Fluxams","AMS_He+P intensity; MC true energy (GeV); J(E) (GeV^{-1} m^{-2} s^{-1} sr^{-1}", nbt, Et)
    if(alpha==0.):
	ytit = "J(E) (GeV s^{-1} m^{-2} sr^{-1})"
    if(alpha==2.):
	ytit = "E^{2} J(E) (GeV s^{-1} m^{-2} sr^{-1})"
    #if(alpha==2.6):
    if(alpha==2.7):
	ytit = "E^{2.7} J(E) (GeV^{1.7} s^{-1} m^{-2} sr^{-1})"
    if(alpha==2.6):
	ytit = "E^{2.6} J(E) (GeV^{1.7} s^{-1} m^{-2} sr^{-1})"
    Fluxamspow = TH1F("Fluxamspow","AMS_He+P intensity; MC true energy (GeV);"+ytit, nbt, Et)

    fpl = TF1("fpl","[0]*pow(x/100., -[1])", 500., 6000.)
    fpl.SetNpx(100000)
    grprams02prl.Fit("fpl","R")
    p0 = fpl.GetParameter(0)
    p1 = fpl.GetParameter(1)
    #fpl1 = TF1("fpl1","[0]*pow(x/100., -[1])", 350., 50000.)
    fpl1 = TF1("fpl1","[0]*pow(x/100., -[1])", 500., 100000.)
    fpl1.SetParameters(p0, p1)
    fpl1.SetNpx(100000)

    #fpl2 = TF1("fpl2","[0]*pow(x/100., -[1])*pow(x, [2])", 350., 50000.)
    fpl2 = TF1("fpl2","[0]*pow(x/100., -[1])*pow(x, [2])", 500., 100000.)
    fpl2.SetParameters(p0, p1, alpha)
    fpl2.SetNpx(100000)

 
    data_p= 'AMS02_2015_Proton_prl.dat'  #Proton
    Emin_p,Emax_p = np.loadtxt(data_p,skiprows=2,usecols=(1,2),unpack=True)
    Emean_p = np.loadtxt(data_p,skiprows=2,usecols=(0,),unpack=True)
    Flux_p  = np.loadtxt(data_p,skiprows=2,usecols=(3,),unpack=True)
    Flux0_p  = np.loadtxt(data_p,skiprows=2,usecols=(3,),unpack=True)
    Flux_low_p  = np.loadtxt(data_p,skiprows=2,usecols=(8,),unpack=True)
    Flux_up_p   = np.loadtxt(data_p,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(Emean_p))
    grprams02prl_p = TGraphAsymmErrors(len(Emean_p), Emean_p, Flux_p, null, null, Flux_low_p, Flux_up_p)
    grprams02prl_p.SetLineColor(3)
    grprams02prl_p.SetMarkerColor(3)
    grprams02prl_p.SetMarkerStyle(21)

    Emin0_p = ary.array("f",range(len(Emean_p)))
    Emax0_p = ary.array("f",range(len(Emean_p)))
    #alpha = 2.7
   # alpha = 2.6
    for i in range(len(Emean_p)):
        e0 = Emin_p[i]
        e1 = Emax_p[i]
        Emin0_p[i] = e0
        Emax0_p[i] = e1
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        Emean_p[i] = em
        #Flux[i] *= TMath.Power(Emean[i], alpha)
        #Flux_low[i] *= TMath.Power(Emean[i], alpha)
        #Flux_up[i] *= TMath.Power(Emean[i], alpha)
        Flux_p[i] *= TMath.Power(em, alpha)
        Flux_low_p[i] *= TMath.Power(em, alpha)
        Flux_up_p[i] *= TMath.Power(em, alpha)

    grprams02pow_p = TGraphAsymmErrors(len(Emean_p), Emean_p, Flux_p, null, null, Flux_low_p, Flux_up_p)
    grprams02pow_p.SetLineColor(3)
    grprams02pow_p.SetMarkerColor(3)
    grprams02pow_p.SetMarkerStyle(21)

    Fluxams_p = TH1F("Fluxams_p","AMS_P intensity; MC true energy (GeV); J(E) (GeV^{-1} m^{-2} s^{-1} sr^{-1}", nbt, Et)
    if(alpha==0.):
	ytit = "J(E) (GeV s^{-1} m^{-2} sr^{-1})"
    if(alpha==2.):
	ytit = "E^{2} J(E) (GeV s^{-1} m^{-2} sr^{-1})"
    if(alpha==2.7):
	ytit = "E^{2.7} J(E) (GeV^{1.7} s^{-1} m^{-2} sr^{-1})"
    if(alpha==2.6):
	ytit = "E^{2.6} J(E) (GeV^{1.7} s^{-1} m^{-2} sr^{-1})"
   
    Fluxamspow_p = TH1F("Fluxamspow_p","AMS_P intensity; MC true energy (GeV);"+ytit, nbt, Et)

    fpl_p = TF1("fpl_p","[0]*pow(x/100., -[1])", 350., 2000.)
    fpl_p.SetNpx(100000)
    grprams02prl_p.Fit("fpl_p","R")
    p0_p = fpl_p.GetParameter(0)
    p1_p = fpl_p.GetParameter(1)
    #fpl1 = TF1("fpl1","[0]*pow(x/100., -[1])", 350., 50000.)
    fpl1_p = TF1("fpl1_p","[0]*pow(x/100., -[1])", 350., 100000.)
    fpl1_p.SetParameters(p0_p, p1_p)
    fpl1_p.SetNpx(100000)

    #fpl2 = TF1("fpl2","[0]*pow(x/100., -[1])*pow(x, [2])", 350., 50000.)
    fpl2_p = TF1("fpl2_p","[0]*pow(x/100., -[1])*pow(x, [2])", 350., 100000.)
    fpl2_p.SetParameters(p0_p, p1_p, alpha)
    fpl2_p.SetNpx(100000)





    for j in range(nbo):
        cp = 0.
        for i in range(nbt):
            e0 = h2sm.GetXaxis().GetBinLowEdge(i+1)
            e1 = h2sm.GetXaxis().GetBinLowEdge(i+2)
            y0 = 0.
            y1 = 0.
            y0_p = 0.
            y1_p = 0.
            yy=0.
	    yy_p=0.   
            if(e0>=Emin[0] and e1<500.):
                y0 = grprams02prl.Eval(e0, 0, "S")
                y1 = grprams02prl.Eval(e1, 0, "S")
                #y1 = Flux0[ie0] + (Flux0[ie1] - Flux0[ie1])/(e11-e10)*(e1-e10)
                yy = 0.
                nk = 100
                for k in range(nk):
                    eem = e0+(e1-e0)*float(k)/float(nk)
                    yy += grprams02prl.Eval(eem, 0, "S")*(e1-e0)/float(nk)
                yy0 = yy/(e1-e0)
            else:
                y0 = p0*pow(e0/100., -p1)
                y1 = p0*pow(e1/100., -p1)
                yy = fpl1.Integral(e0, e1)
                yy0 = yy/(e1-e0)
            if(e0>=Emin[0] and e1<350.):
                y0_p = grprams02prl_p.Eval(e0, 0, "S")
                y1_p = grprams02prl_p.Eval(e1, 0, "S")
                #y1 = Flux0[ie0] + (Flux0[ie1] - Flux0[ie1])/(e11-e10)*(e1-e10)
                yy_p = 0.
                nk = 100
                for k in range(nk):
                    eem = e0+(e1-e0)*float(k)/float(nk)
                    yy_p += grprams02prl_p.Eval(eem, 0, "S")*(e1-e0)/float(nk)
                yy0_p = yy_p/(e1-e0)
            else:
                y0_p = p0_p*pow(e0/100., -p1_p)
                y1_p = p0_p*pow(e1/100., -p1_p)
                yy_p = fpl1_p.Integral(e0, e1)
                yy0_p = yy_p/(e1-e0)
            
            if(j==0):
                Fluxams.SetBinContent(i+1, yy0+yy0_p)
                em = TMath.Sqrt(e0*e1)
                em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
                em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
                Fluxamspow.SetBinContent(i+1, (yy0+yy0_p)*pow(em, alpha))
                counts = h1ef.GetBinContent(i+1)*TotTime*(yy+yy_p)
                errf = 0.
                if(counts>0.):
                    errf =  (yy+yy_p)/TMath.Sqrt(counts)
                Fluxams.SetBinError(i+1, errf)
                Fluxamspow.SetBinError(i+1, errf*pow(em, alpha))
            cp += h2sm.GetBinContent(i+1, j+1)*(yy+yy_p)
        cp *= TotTime
        h1cts_ams02.SetBinContent(j+1, cp)
        h1cts.SetBinContent(j+1, cp)
        
    for i in range(nbt):
        e0 = h2sm.GetXaxis().GetBinLowEdge(i+1)
        e1 = h2sm.GetXaxis().GetBinLowEdge(i+2)
        counts = h1cts.GetBinContent(i+1)
        flux = Fluxams.GetBinContent(i+1)
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        errf = 0.
        if(counts>0.):
            errf =  0.5*(y0+y1)/TMath.Sqrt(counts)
        errf = 0.
        Fluxams.SetBinError(i+1, errf)
        Fluxamspow.SetBinError(i+1, errf*pow(em, alpha))

    # read real data



    # read real data


    fdat = TFile(Orb_file)
    h1ctsdata = fdat.Get(BGO_data)
    #print "Data = ", h1ctsdata.GetEntries()
    print "Data = ", h1ctsdata.GetSumOfWeights()
    h1ctsdata.Rebin(nrbo)

    # Unfolding !!!

    maxeff =  1./h1ef.GetMaximum()
    #maxeff = 1.
    print "Max eff = ", 1./maxeff


    # set regularization
    creg = "none"
    #creg = "PowerLaw"
    #creg = "SmoothArray"
    ireg = 0
    if (creg == "none"):
        ireg = 0
    if (creg == "PowerLaw"):
        ireg = 1
    if (creg == "ExpCutoff"):
        ireg = 2
    if (creg == "BrokenPowerLaw"):
        ireg = 3
    if (creg == "SmoothArray"):
        ireg = 4

    # error bars want to include in unfolding histograms

    cerr = "statistical"
    #cerr = "overall"
    isys = 0
    if (cerr == "statistical"):
        isys = 0
    if (cerr == "systematic"):
        isys = 1
    if (cerr == "overall"):
        isys = 2

    # Initialize vectors for the unfolding

    noc = nbt
    noe = nbo
    
    nc = ary.array('f',range(noc))
    enc = ary.array('f',range(noc))
    enc1 = ary.array('f',range(noc))
    nc_mc = ary.array('f',range(noc))
    
    ne = ary.array('f',range(noe))
    neeq = ary.array('f',range(noe))

    pec = [[0.0 for col in range(noc)] for raw in range(noe)]
    fec = [[0.0 for col in range(noc)] for raw in range(noe)]

    Vc0 = [[0.0 for col in range(noc)]for raw in range(noc)]
    Vc1 = [[0.0 for col in range(noc)]for raw in range(noc)]
     

    noc11 = 0.    
    for i in range(noc):
        
        # mc_mc[i] = number of MC events in the i-th cause bin
        #nc_mc[i] = h1ng.GetBinContent(i+1+4) #if the bin numbers of h1ng and h2nd are different with 4(after rebin(8)),please use this one.
	nc_mc[i] = h1ng.GetBinContent(i+1) #if the bins of h1ng and h2nd are the same, please use this one!. 

        for j in range(noe):
            
            # pec[j][i] = number of events in the j-th effect bin from the i-th cause bin
            xj = h2nd.GetBinContent(i+1,j+1)*maxeff
            ej = h1ctsdata.GetBinLowEdge(j+1)
            #if (ej<90.):
            #    xj=0.            
            #    h2sm.SetBinContent(i+1, j+1, 0.)
            if ( nc_mc[i] > 0. ):
                pec[j][i] = xj / nc_mc[i]
                
            if (j == 0):
                fec[j][i] = pec[j][i]
            else: 
                fec[j][i] = fec[j-1][i] + pec[j][i]
                

        #SEtrue.SetBinContent(i+1,fec[noe-1][i])
        if(fec[noe-1][i]>0):
            noc11 = noc11 + 1
                
##     # Set small probabilities to zero 
##     for i in range(noc):
##         if (fec[noe-1][i] < 5.e-5):
##             for j in range(noe):
##                 pec[j][i] = 0.
                
    
  
    print "AMS Total observed Events = ",str(h1cts.Integral())

    # Evaluate equivalent number of events and total number of events
    net = 0
    for j in range (noe):
        ne[j] = 0.
#                                         Usa gli "eventi equivalenti" di AMS per l'unfolding:
        xj = h1cts.GetBinContent(j+1)
        ej = h1ctsdata.GetBinLowEdge(j+1)
#                                         Usa i Dati per l'unfolding:
        xj = h1ctsdata.GetBinContent(j+1)
        #if (xj < 0.0 or ej<90.):
        #    xj=0.            
        #    h1ctsdata.SetBinContent(j+1, 0.)
        ne[j] = xj
        neeq[j] = TMath.Sqrt(xj)
        neeq[j] = neeq[j] * neeq[j]
        net = net + ne[j]

    print "net = ", str(net) 

    # Build vector of fractional systematic errors on energies
                
      



    # Define histograms
    
    Flux = TH1F("Flux","Unfolded intensity; MC true energy (GeV); J(E) (GeV^{-1} m^{-2} s^{-1} sr^{-1}", noc, Et)
    Fluxpow = TH1F("Fluxpow","Unfolded intensity; MC true energy (GeV);"+ytit, noc, Et)
    CovMatrix0 = TH2F("CovMatrix0", "Covariance Matrix - data contribution", noc, Et, noc, Et)
    CovMatrix1 = TH2F("CovMatrix1", "Covariance Matrix - MC contribution", noc, Et, noc, Et)

    unfsmooth = TH1F("unfsmooth", "Test histogram", noc, Et)

    print noc
    print Et

    # Define power-law function to initialize nc
    fl = TF1("fl","[0]*TMath::Power(x, -[1])", Et[0], Et[noc])
    gamma = 2.7
    p0 = (gamma-1.) * TMath.Power(Et[0], gamma-1.)
    fl.SetParameters(p0, gamma)

    for i in range(noc):
        enc[i] = 0.
        enc1[i] = 0.
        # starting hypotehsis: uniform distribution
        # nc[i] = net / noc11
        # starting hypothesis: power law distribution
        E1 = Flux.GetBinLowEdge(i+1)
        E2 = Flux.GetBinLowEdge(i+2)
        nc[i] = fl.Integral(E1, E2) 
        #print i, nc[i]
  
    ier = -1
    iregul = ireg
    errcal = 2
    ibkg = 0

    # Perform the unfolding 
    ier = fastbayes(noe, ne, neeq, noc, nc, pec, iter, enc, enc1, nc_mc, iregul, errcal, unfsmooth, ibkg, Vc0, Vc1, maxeff)
    print "fastbayes output error = ",str(ier)

    # Now build the histograms (and the associated errors)

    fobj=open(DataFileOutput,'w')
    fobj.write('#DAMPE_results\n')
    fobj.write('#   <E>  Elo  Eup   y   ystat_lo  ystat_up  ysyst_lo  ysyst_up  yerrtot_lo  yerrtot_up \n')

    print("noc is:",noc)
    # Evaluate unfolded distribution, flux and SED       
    for i in range(noc):
        DEne = Et[i+1] - Et[i]
        xi = nc[i] / (DEne*TotTime)
        if (isys == 0):
            xj = enc[i]/(DEne*TotTime)
        if (isys == 1):
            xj = enc1[i]/(DEne*TotTime)
        if (isys == 2):
            xj = TMath.Sqrt(enc[i]*enc[i]+ enc1[i]*enc1[i])/(DEne*TotTime)
        Flux.SetBinContent(i+1,xi)
        Flux.SetBinError(i+1,xj)
        fobj.write(str(TMath.Sqrt(Et[i]*Et[i+1]))+'\t'+str(Et[i])+'\t'+str(Et[i+1])+'\t'+str(xi)+'\t'+str(xj)+'\t'+str(xj)+'\t'+'0.00'+'\t'+'0.00'+'\t'+str(xj)+'\t'+str(xj)+'\n')
    fobj.close() 

    h1cntUnfod = TH1F("h1cntUnfod","Counts after Unfolding", noc, Et)
    for i in range (noc):
         
          DEne = Et[i+1] - Et[i]
          print ("FIRST: Energy should be in",Et[i],"and ",Et[i+1])
          a =(Flux.GetBinContent(i+1)*(DEne*TotTime)*(h1ef.GetBinContent(i+1))) 
          print("In bin-%d:%dGeV-%dGeV Count is:%d" %(i+1,Flux.GetBinLowEdge(i+1),Flux.GetBinLowEdge(i+2),a))
          h1cntUnfod.SetBinContent(i+1,a*pow(DEne,2.6))
    # Fill histograms with covariance matrix contents
    #exit(0)
    for i in range (noc):
        E1 = Flux.GetBinLowEdge(i+1)
        E2 = Flux.GetBinLowEdge(i+2)
        em = TMath.Sqrt(E1*E2)
        em = math.pow(0.5*(math.pow(E2 ,-alpha+1)+math.pow(E1 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(E2 ,-alpha+1.)-math.pow(E1 ,-alpha+1.))/((E2-E1)*(-alpha+1.)), 1./(-alpha))
        y = Flux.GetBinContent(i+1)*pow(em, alpha)
        Fluxpow.SetBinContent(i+1, y)
        y = Flux.GetBinError(i+1)*pow(em, alpha)
        Fluxpow.SetBinError(i+1, y)
        for j in range (noc):
            CovMatrix0.SetBinContent(i+1, j+1, Vc0[i][j]/(TotTime*TotTime))
            CovMatrix1.SetBinContent(i+1, j+1, Vc1[i][j]/(TotTime*TotTime))

 

    gStyle.SetOptStat(0)
    c0 = TCanvas("c0", "Efficency")
    gPad.SetLogx()
    h1ef.SetLineColor(0)
    h1ef.SetMarkerColor(2)
    h1ef.SetMarkerStyle(20)
    h1ef.Draw("p")
   # histeff= frsm.Get("h1Eff_w");
    #histeff.Rebin(nrbt);
    #histeff.SetMarkerStyle(21)
    #histeff.SetMarkerColor(1)

    #histeff.Draw("sames")
    c00 = TCanvas("c00", "Counts_After Unfolding")
    gPad.SetLogx()
    gPad.SetLogy()

    
    h1cntUnfod.SetLineColor(1)
    h1cntUnfod.SetMarkerColor(1)
    h1cntUnfod.SetMarkerStyle(20)
    h1cntUnfod.GetXaxis().SetTitle("Energy (GeV)")
    h1cntUnfod.GetXaxis().CenterTitle()
    h1cntUnfod.GetYaxis().SetTitle("Counts After Unfolding")
    h1cntUnfod.Draw("p")

    c1 = TCanvas("c1", "Acceptance")
    gPad.SetLogx()
    h1ac.SetLineColor(1)
    h1ac.SetMarkerColor(1)
    h1ac.SetMarkerStyle(20)
    h1ac.Draw("p")


    #h2sm.Scale(Amc)
    c2 = TCanvas("c2", "Smearing Matrix")
    c2.SetRightMargin(0.15)
    gPad.SetLogx()
    gPad.SetLogy()
    gPad.SetLogz()
    h2sm.Draw("colz")

    c3 = TCanvas("c3", "Observed counts")
    gPad.SetLogx()
    gPad.SetLogy()
    h1cts.SetLineColor(1)
    h1cts.SetMarkerColor(1)
    h1cts.SetMarkerStyle(20)
    #h1cts_ams02.Draw("e")
    #h1cts.Draw("esame")
    h1ctsdata.GetXaxis().SetTitle("Energy (GeV)")
    h1ctsdata.GetXaxis().CenterTitle()
    h1ctsdata.GetYaxis().SetTitle("Counts Before Unfolding")
    h1ctsdata.SetLineColor(1)
    h1ctsdata.SetMarkerColor(1)
    h1ctsdata.SetMarkerStyle(20)
    h1ctsdata.Draw("esame")
   
    for ii in range(h1ctsdata.GetNbinsX()):

    	DEne = Et[i+1] - Et[i]
        print ("SECOND: Energy should be in",Et[i],"and ",Et[i+1])
        a =(h1ctsdata.GetBinContent(i+1)*(DEne*TotTime)*(h1ef.GetBinContent(i+1))) 
        print("In bin-%d:%dGeV-%dGeV Count is:%d" %(i+1,h1ctsdata.GetBinLowEdge(i+1),h1ctsdata.GetBinLowEdge(i+2),a))
        #  h1cntUnfod.SetBinContent(i+1,a)

    c4 = TCanvas("c4", "Proton Intensity")
    gPad.SetLogx()
    gPad.SetLogy()

    #frame2 = gPad.DrawFrame(0.100, 1.e-8, 5.e4, 2.01e3)
    frame2 = gPad.DrawFrame(0.100, 1.e-10, 2.e5, 2.01e3)
    frame2.SetTitle('')
    frame2.GetXaxis().CenterTitle()
    frame2.GetYaxis().CenterTitle()
    frame2.GetYaxis().SetTitleOffset(1.2)
    frame2.GetXaxis().SetTitle("Energy (GeV)")
    frame2.GetYaxis().SetTitle("J(E) (GeV^{-1} s^{-1} m^{-1} sr^{-1})")
    frame2.GetXaxis().SetLabelOffset(0.001)
    frame2.GetXaxis().SetLabelOffset(0.0001)

    grprams02prl.Draw("psame")
    grprams02prl_p.Draw("psame")
    
    fpl1.Draw("same")
    fpl1_p.Draw("same")
    Flux.SetLineColor(1)
    Flux.SetMarkerColor(1)
    Flux.SetMarkerStyle(20)
    Flux.Draw("lpsame")
    Fluxams.SetLineColor(4)
    Fluxams.Draw("same")

    '''
    

    dataC3He = 'CREAM-III_Helium0.dat'  #Helium
    EminC3He,EmaxC3He = np.loadtxt(dataC3He,skiprows=2,usecols=(1,2),unpack=True)
    EmeanC3He = np.loadtxt(dataC3He,skiprows=2,usecols=(0,),unpack=True)
    FluxC3He  = np.loadtxt(dataC3He,skiprows=2,usecols=(3,),unpack=True)
    Flux0C3He  = np.loadtxt(dataC3He,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowC3He  = np.loadtxt(dataC3He,skiprows=2,usecols=(8,),unpack=True)
    Flux_upC3He   = np.loadtxt(dataC3He,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanC3He))
    grC3He= TGraphAsymmErrors(len(EmeanC3He), EmeanC3He, FluxC3He, null, null, Flux_lowC3He, Flux_upC3He)
    grC3He.SetLineColor(2)
    grC3He.SetMarkerColor(2)
    grC3He.SetMarkerStyle(20)

    dataC3P = 'CREAM-III_Proton0.dat'  #proton
    EminC3P,EmaxC3P = np.loadtxt(dataC3P,skiprows=2,usecols=(1,2),unpack=True)
    EmeanC3P = np.loadtxt(dataC3P,skiprows=2,usecols=(0,),unpack=True)
    FluxC3P  = np.loadtxt(dataC3P,skiprows=2,usecols=(3,),unpack=True)
    Flux0C3P  = np.loadtxt(dataC3P,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowC3P  = np.loadtxt(dataC3P,skiprows=2,usecols=(8,),unpack=True)
    Flux_upC3P   = np.loadtxt(dataC3P,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanC3P))
    grC3P= TGraphAsymmErrors(len(EmeanC3P), EmeanC3P, FluxC3P, null, null, Flux_lowC3P, Flux_upC3P)
    grC3P.SetLineColor(3)
    grC3P.SetMarkerColor(3)
    grC3P.SetMarkerStyle(21)

    EmeanC3PHe = np.zeros(len(EmeanC3P))

    FluxC3PHe =  np.zeros(len(EmeanC3P))
    FluxC3_lowPHe =  np.zeros(len(EmeanC3P))
    FluxC3_upPHe =  np.zeros(len(EmeanC3P))

    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanC3He)):
        e0 = EminC3P[i]
        e1 = EmaxC3P[i]
	FluxC3PHe[i]= FluxC3P[i]+FluxC3He[i]
	FluxC3_lowPHe[i] = math.sqrt(Flux_lowC3P[i]**2+Flux_lowC3He[i]**2)
	FluxC3_upPHe[i] = math.sqrt(Flux_upC3P[i]**2+Flux_upC3He[i]**2)
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanC3PHe[i] = em
        #Flux[i] *= TMath.Power(Emean[i], alpha)
        #Flux_low[i] *= TMath.Power(Emean[i], alpha)
        #Flux_up[i] *= TMath.Power(Emean[i], alpha)
        FluxC3PHe[i] *= TMath.Power(em, alpha)
        FluxC3_lowPHe[i] *= TMath.Power(em, alpha)
        FluxC3_upPHe[i] *= TMath.Power(em, alpha)

    null = np.zeros(len(EmeanC3PHe))
    print("len EmeanC3PHe:", len(EmeanC3PHe))
    grC3PHepow = TGraphAsymmErrors(len(EmeanC3PHe), EmeanC3PHe, FluxC3PHe, null, null, FluxC3_lowPHe, FluxC3_upPHe)
    grC3PHepow.SetLineColor(4)
    grC3PHepow.SetMarkerColor(4)
    grC3PHepow.SetMarkerStyle(22)


    dataNuHe = 'Nucleon_He.dat'  #Helium
    EminNuHe,EmaxNuHe = np.loadtxt(dataNuHe,skiprows=2,usecols=(1,2),unpack=True)
    EmeanNuHe = np.loadtxt(dataNuHe,skiprows=2,usecols=(0,),unpack=True)
    FluxNuHe  = np.loadtxt(dataNuHe,skiprows=2,usecols=(3,),unpack=True)
    Flux0NuHe  = np.loadtxt(dataNuHe,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowNuHe  = np.loadtxt(dataNuHe,skiprows=2,usecols=(8,),unpack=True)
    Flux_upNuHe   = np.loadtxt(dataNuHe,skiprows=2,usecols=(9,),unpack=True)
    print EminNuHe
    
    null = np.zeros(len(EmeanNuHe))
    grNuHe= TGraphAsymmErrors(len(EmeanNuHe), EmeanNuHe, FluxNuHe, null, null, Flux_lowNuHe, Flux_upNuHe)
    #grA2He.SetLineColor(2)
    #grA2He.SetMarkerColor(2)
    #grA2He.SetMarkerStyle(20)

    dataNuP = 'Nucleon_p.dat'  # proton
    EminNuP,EmaxNuP = np.loadtxt(dataNuP,skiprows=2,usecols=(1,2),unpack=True)
    EmeanNuP = np.loadtxt(dataNuP,skiprows=2,usecols=(0,),unpack=True)
    FluxNuP  = np.loadtxt(dataNuP,skiprows=2,usecols=(3,),unpack=True)
    Flux0NuP  = np.loadtxt(dataNuP,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowNuP  = np.loadtxt(dataNuP,skiprows=2,usecols=(8,),unpack=True)
    Flux_upNuP   = np.loadtxt(dataNuP,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanNuP))
    grNuP= TGraphAsymmErrors(len(EmeanNuP), EmeanNuP, FluxNuP, null, null, Flux_lowNuP, Flux_upNuP)
    #grA2P.SetLineColor(3)
    #grA2P.SetMarkerColor(3)
    #grA2P.SetMarkerStyle(21)

    EmeanNuPHe = np.zeros(len(EmeanNuP))

    FluxNuPHe =  np.zeros(len(EmeanNuP))
    FluxNu_lowPHe =  np.zeros(len(EmeanNuP))
    FluxNu_upPHe =  np.zeros(len(EmeanNuP))

    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanNuHe)):
        e0 = EminNuP[i]
        e1 = EmaxNuP[i]
        #print "eo:",e0," e1: ", e1
	FluxNuPHe[i]= FluxNuP[i]+FluxNuHe[i]
	FluxNu_lowPHe[i] = math.sqrt(Flux_lowNuP[i]**2+Flux_lowNuHe[i]**2)
	FluxNu_upPHe[i] = math.sqrt(Flux_upNuP[i]**2+Flux_upNuHe[i]**2)
        #print i, e0, e1
       # em1  = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em1  = EmeanNuP[i]
       # em1 = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanNuPHe[i] = em1 
        #Flux[i] *= TMath.Power(Emean[i], alpha)
        #Flux_low[i] *= TMath.Power(Emean[i], alpha)
        #Flux_up[i] *= TMath.Power(Emean[i], alpha)
        FluxNuPHe[i] *= TMath.Power(em1, alpha)
        FluxNu_lowPHe[i] *= TMath.Power(em1, alpha)
        FluxNu_upPHe[i] *= TMath.Power(em1, alpha)

    null = np.zeros(len(EmeanNuPHe))
    grNuPHepow = TGraphAsymmErrors(len(EmeanNuPHe), EmeanNuPHe, FluxNuPHe, null, null, FluxNu_lowPHe, FluxNu_upPHe)
    grNuPHepow.SetLineColor(2)
    grNuPHepow.SetMarkerColor(2)
    grNuPHepow.SetMarkerStyle(23)


    dataHAWC = 'HAWC_AllParticles.dat' #DAMPE
    EminHAWC,EmaxHAWC = np.loadtxt(dataHAWC,skiprows=2,usecols=(1,2),unpack=True)
    EmeanHAWC = np.loadtxt(dataHAWC,skiprows=2,usecols=(0,),unpack=True)
    FluxHAWC  = np.loadtxt(dataHAWC,skiprows=2,usecols=(3,),unpack=True)
    Flux0HAWC  = np.loadtxt(dataHAWC,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowHAWC  = np.loadtxt(dataHAWC,skiprows=2,usecols=(8,),unpack=True)
    Flux_upHAWC  = np.loadtxt(dataHAWC,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanHAWC))
    grHAWC= TGraphAsymmErrors(len(EmeanHAWC), EmeanHAWC, FluxHAWC, null, null, Flux_lowHAWC, Flux_upHAWC)
    grHAWC.SetLineColor(5)
    grHAWC.SetMarkerColor(5)
    grHAWC.SetMarkerStyle(25)



    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanHAWC)):
        e0 = EminHAWC[i]
        e1 = EmaxHAWC[i]
	
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanHAWC[i] = em

        FluxHAWC[i] *= TMath.Power(em, alpha)
        Flux_lowHAWC[i] *= TMath.Power(em, alpha)
        Flux_upHAWC[i] *= TMath.Power(em, alpha)

    null = np.zeros(len(EmeanHAWC))
    print("len EmeanDMP:", len(EmeanHAWC))
    grHAWCpow = TGraphAsymmErrors(len(EmeanHAWC), EmeanHAWC, FluxHAWC, null, null, Flux_lowHAWC, Flux_upHAWC)
    grHAWCpow.SetLineColor(3)
    grHAWCpow.SetMarkerColor(3)
    grHAWCpow.SetMarkerStyle(23)
    grHAWCpow.SetMarkerSize(1.5)
      



      
    dataA2He = 'ATIC2_Helium0.dat'  #Helium
    EminA2He,EmaxA2He = np.loadtxt(dataA2He,skiprows=2,usecols=(1,2),unpack=True)
    EmeanA2He = np.loadtxt(dataA2He,skiprows=2,usecols=(0,),unpack=True)
    FluxA2He  = np.loadtxt(dataA2He,skiprows=2,usecols=(3,),unpack=True)
    Flux0A2He  = np.loadtxt(dataA2He,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowA2He  = np.loadtxt(dataA2He,skiprows=2,usecols=(8,),unpack=True)
    Flux_upA2He   = np.loadtxt(dataA2He,skiprows=2,usecols=(9,),unpack=True)
   # print EminA2He
   
    
    null = np.zeros(len(EmeanA2He))
    grA2He= TGraphAsymmErrors(len(EmeanA2He), EmeanA2He, FluxA2He, null, null, Flux_lowA2He, Flux_upA2He)
    #grA2He.SetLineColor(2)
    #grA2He.SetMarkerColor(2)
    #grA2He.SetMarkerStyle(20)

    dataA2P = 'ATIC2_Proton0.dat'  #Helium
    EminA2P,EmaxA2P = np.loadtxt(dataA2P,skiprows=2,usecols=(1,2),unpack=True)
    EmeanA2P = np.loadtxt(dataA2P,skiprows=2,usecols=(0,),unpack=True)
    FluxA2P  = np.loadtxt(dataA2P,skiprows=2,usecols=(3,),unpack=True)
    Flux0A2P  = np.loadtxt(dataA2P,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowA2P  = np.loadtxt(dataA2P,skiprows=2,usecols=(8,),unpack=True)
    Flux_upA2P   = np.loadtxt(dataA2P,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanA2P))
    grA2P= TGraphAsymmErrors(len(EmeanA2P), EmeanA2P, FluxA2P, null, null, Flux_lowA2P, Flux_upA2P)
    #grA2P.SetLineColor(3)
    #grA2P.SetMarkerColor(3)
    #grA2P.SetMarkerStyle(21)

    EmeanA2PHe = np.zeros(len(EmeanA2P))

    FluxA2PHe =  np.zeros(len(EmeanA2P))
    FluxA2_lowPHe =  np.zeros(len(EmeanA2P))
    FluxA2_upPHe =  np.zeros(len(EmeanA2P))

    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanA2He)):
        e0 = EminA2P[i]
        e1 = EmaxA2P[i]
        #print "eo:",e0," e1: ", e1
	FluxA2PHe[i]= FluxA2P[i]+FluxA2He[i]
	FluxA2_lowPHe[i] = math.sqrt(Flux_lowA2P[i]**2+Flux_lowA2He[i]**2)
	FluxA2_upPHe[i] = math.sqrt(Flux_upA2P[i]**2+Flux_upA2He[i]**2)
        #print i, e0, e1
       # em1  = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em1  = EmeanA2P[i]
       # em1 = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanA2PHe[i] = em1
        #Flux[i] *= TMath.Power(Emean[i], alpha)
        #Flux_low[i] *= TMath.Power(Emean[i], alpha)
        #Flux_up[i] *= TMath.Power(Emean[i], alpha)
        FluxA2PHe[i] *= TMath.Power(em1, alpha)
        FluxA2_lowPHe[i] *= TMath.Power(em1, alpha)
        FluxA2_upPHe[i] *= TMath.Power(em1, alpha)

    null = np.zeros(len(EmeanA2PHe))
    grA2PHepow = TGraphAsymmErrors(len(EmeanA2PHe), EmeanA2PHe, FluxA2PHe, null, null, FluxA2_lowPHe, FluxA2_upPHe)
    grA2PHepow.SetLineColor(6)
    grA2PHepow.SetMarkerColor(6)
    grA2PHepow.SetMarkerStyle(23)

  

   
    dataARGOG4 = 'ARGO_PHe_G4.dat' #DAMPE
    EminARGOG4,EmaxARGOG4 = np.loadtxt(dataARGOG4,skiprows=2,usecols=(1,2),unpack=True)
    EmeanARGOG4 = np.loadtxt(dataARGOG4,skiprows=2,usecols=(0,),unpack=True)
    FluxARGOG4  = np.loadtxt(dataARGOG4,skiprows=2,usecols=(3,),unpack=True)
    Flux0ARGOG4  = np.loadtxt(dataARGOG4,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowARGOG4  = np.loadtxt(dataARGOG4,skiprows=2,usecols=(8,),unpack=True)
    Flux_upARGOG4  = np.loadtxt(dataARGOG4,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanARGOG4))
    grARGOG4= TGraphAsymmErrors(len(EmeanARGOG4), EmeanARGOG4, FluxARGOG4, null, null, Flux_lowARGOG4, Flux_upARGOG4)
    grARGOG4.SetLineColor(5)
    grARGOG4.SetMarkerColor(5)
    grARGOG4.SetMarkerStyle(25)



    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanARGOG4)):
        e0 = EminARGOG4[i]
        e1 = EmaxARGOG4[i]
	
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanARGOG4[i] = em

        FluxARGOG4[i] *= TMath.Power(em, alpha)
        Flux_lowARGOG4[i] *= TMath.Power(em, alpha)
        Flux_upARGOG4[i] *= TMath.Power(em, alpha)

    null = np.zeros(len(EmeanARGOG4))
    print("len EmeanDMP:", len(EmeanARGOG4))
    grARGOG4pow = TGraphAsymmErrors(len(EmeanARGOG4), EmeanARGOG4, FluxARGOG4, null, null, Flux_lowARGOG4, Flux_upARGOG4)
    grARGOG4pow.SetLineColor(2)
    grARGOG4pow.SetMarkerColor(2)
    grARGOG4pow.SetMarkerStyle(22)
    grARGOG4pow.SetMarkerSize(1.5)
    '''
    """   


     
    dataPChuan = 'DAMPE_Proton.dat' #DAMPE
    EminPChuan,EmaxPChuan = np.loadtxt(dataPChuan,skiprows=2,usecols=(1,2),unpack=True)
    EmeanPChuan = np.loadtxt(dataPChuan,skiprows=2,usecols=(0,),unpack=True)
    FluxPChuan  = np.loadtxt(dataPChuan,skiprows=2,usecols=(3,),unpack=True)
    Flux0PChuan  = np.loadtxt(dataPChuan,skiprows=2,usecols=(3,),unpack=True)
    Flux_lowPChuan  = np.loadtxt(dataPChuan,skiprows=2,usecols=(8,),unpack=True)
    Flux_upPChuan  = np.loadtxt(dataPChuan,skiprows=2,usecols=(9,),unpack=True)

    null = np.zeros(len(EmeanPChuan))
    grPChuan= TGraphAsymmErrors(len(EmeanPChuan), EmeanPChuan, FluxPChuan, null, null, Flux_lowPChuan, Flux_upPChuan)
    grPChuan.SetLineColor(5)
    grPChuan.SetMarkerColor(5)
    grPChuan.SetMarkerStyle(25)



    #alpha = 2.6
    #alpha = 2.7
    for i in range(len(EmeanPChuan)):
        e0 = EminPChuan[i]
        e1 = EmaxPChuan[i]
	
        #print i, e0, e1
        em = math.pow(0.5*(math.pow(e1 ,-alpha+1)+math.pow(e0 ,-alpha+1)), -1./(alpha-1))
        em = math.pow( (math.pow(e1 ,-alpha+1.)-math.pow(e0 ,-alpha+1.))/((e1-e0)*(-alpha+1.)), 1./(-alpha))
        EmeanPChuan[i] = em

        FluxPChuan[i] *= TMath.Power(em, alpha)
        Flux_lowPChuan[i] *= TMath.Power(em, alpha)
        Flux_upPChuan[i] *= TMath.Power(em, alpha)

    null = np.zeros(len(EmeanPChuan))
    print("len EmeanDMP:", len(EmeanPChuan))
    grPChuanpow = TGraphAsymmErrors(len(EmeanPChuan), EmeanPChuan, FluxPChuan, null, null, Flux_lowPChuan, Flux_upPChuan)
    grPChuanpow.SetLineColor(6)
    grPChuanpow.SetMarkerColor(6)
    grPChuanpow.SetMarkerStyle(25)
    grPChuanpow.SetMarkerSize(1.5)
    """   



    #c5 = TCanvas("c5", "Proton Intensity pow")
    c5 = TCanvas("c5", "Proton Intensity pow (Data 6.0 (20 months))")
    gPad.SetLogx()
    gPad.SetLogy()

    #frame3 = gPad.DrawFrame(0.100, 4.999e2, 5.e4, 5.001e4)
    frame3 = gPad.DrawFrame(10.000, 4.e3, 1.e6, 6.e4)
    frame3.SetTitle('')
    frame3.GetXaxis().CenterTitle()
    frame3.GetYaxis().CenterTitle()
    frame3.GetYaxis().SetTitleOffset(1.2)
    frame3.GetXaxis().SetTitle("Energy (GeV)")
    frame3.GetYaxis().SetTitle("J(E) (GeV^{1.7} s^{-1} m^{-1} sr^{-1})")
    frame3.GetYaxis().SetTitle(ytit)
    frame3.GetXaxis().SetLabelOffset(0.001)
    frame3.GetXaxis().SetLabelOffset(0.0001)
    



   # grprams02pow.Draw("psame")
   # grprams02pow_p.Draw("psame")
    #grA2PHepow.Draw("psame")
    #grNuPHepow.Draw("psame")
    #grC3PHepow.Draw("psame")
   # grARGOG4pow.Draw("psame")
   # grHAWCpow.Draw("psame")
    #grPChuanpow.Draw("psame")




    fpl2.SetLineColor(2)
    fpl2_p.SetLineColor(3)
   # fpl2.Draw("sames")  #AMS proton HE line
   # fpl2_p.Draw("sames")
    #fpl2.Draw("same")
    Fluxamspow.SetMarkerColor(4)
    Fluxamspow.SetLineColor(4)
    Fluxamspow.SetMarkerStyle(24)
    #Fluxamspow.Draw("histsame")
   # Fluxamspow.Draw("sames")
    Fluxpow.SetLineColor(4)
    Fluxpow.SetMarkerColor(4)
    Fluxpow.SetMarkerStyle(20)
    Fluxpow.SetMarkerSize(1.1)
    Fluxpow.Draw("psame")
    
    legend =TLegend(0.1,0.7,0.48,0.9)
    #legend.AddEntry(grprams02pow_p,"AMS Proton","lep")
   # legend.AddEntry(grprams02pow,"AMS Helium","lep")
    #legend.AddEntry(grA2PHepow,"ATIC02","lep")
    #legend.AddEntry(grNuPHepow,"Nucleon(KLEM)","lep")
    #legend.AddEntry(grC3PHepow,"CREAM-III","lep")
   # legend.AddEntry(grARGOG4pow,"ARGO 2015","lep")
    #legend.AddEntry(grHAWCpow,"HAWC All-Particles 2017","lep")
   # legend.AddEntry(Fluxamspow,"AMS-2","l")
    #legend.AddEntry(grPChuanpow,"DAMPE-Proton","lep")
    legend.AddEntry(Fluxpow,ResultNameLegend,"lep")

    legend.Draw()


    raw_input("Press enter..")

