# DoUnfolding.py performs the unfolding of the proton spectrum
# author M. Nicola Mazziotta (mazziotta@ba.infn.it), Nov 2016
# Plese cite: arXiv:0912.1236
#             arXiv:0912.3695
#             NIM A362 (1995) 487 

import sys
from ROOT import *
import array as ary
import numpy as np
import math

class PowerLaw:
    def __call__(self, x, par):
        if (par[1]>0.0):
            y = par[0] * TMath.Power(x[0], -par[1])
        else:
            y = par[0]
        return y

class Linear:
    def __call__(self, x, par):
        y = par[0] - x[0]*par[1]
        return y

class PowerLawExpCut:
    def __call__(self, x, par):
        if (par[1]>0.0 and par[2]>0.0):
            y = par[0] * pow(x[0], -par[1]) * TMath.Exp(-x[0]/par[2])
        else:
            y = par[0]
        return y

class LinearExpCut:
    def __call__(self, x, par):
        y = par[0] - par[1]*x[0] - par[2]*TMath.Power(10.0, x[0]) 
        return y

class BrokenPowerLaw:
    def __call__(self, x, par):
        if (par[1]>0.0 and par[2]>0.0 and par[3]>0.0):
            if (x[0] < par[2]):
                y = par[0] * TMath.Power(x[0]/par[2], -par[1])
            else:
                y = par[0] * TMath.Power(x[0]/par[2], -par[3])
        else:
            y = par[0]
        return y

class LinearBrokenPowerLaw:
    def __call__(self, x, par):
        if (x[0] < par[3]):            
            y = par[0] - par[1]*(x[0]-par[3])
        else:
            y = par[0] - par[2]*(x[0]-par[3])
        return y

# Definition of fastbayes function

def fastbayes(noe, ne, neeq, noc, nc, pec, iter, enc, enc1, nc_mc, iregul, errcal, unfsmooth, ibkg, Vc0, Vc1, maxeff):

    #   SOME NOTATIONS
    # 
    #   noe = number of effect bins
    #
    #   ne[i] = number of entries in the i-th effect bin
    #
    #   neeq[i] = number of equivalent entries in the i-th effect bin
    #
    #   noc = number of cause bins
    #
    #   nc[j] = number of reconstructed events in the j-th cause bin
    #
    #   pec[i][j] = conditional probability P(E_i|C_j) (estimated by MC)
    #
    #   iter = iteration number
    #
    #   enc[j] = error on the events in the j-th cause bin -
    #            contribution due to the statistical errors on the observed events ne[i] 
    #
    #   enc1[j] = error on the events in the j-th cause bin -
    #             contribution due to the statistical errors on the MC used to build the
    #             smearing matrix
    #
    #   Vc0[i][j]   : contribution to the covariance matrix of nc[i] due
    #                 to the observed events and assuming the constraint 
    #                 to the total number of observed events,
    #                 and taking into account also the 
    #                 overall normalization uncertainty due to the
    #                 total number of observed events
    #
    #   Vc1[i][j]   : contribution to the covariance matrix of nc[i] due
    #                 to limited statistic used to evaluate the smearing matrix  
    #
    #
    #   Vc2[i][j]   : contribution to the covariance matrix of nc[i] due
    #                 to systematics  
    #
    #   unfsmooth = histogram needed to perform regularization (defined as the observed spectrum histogram) 
    #
    #   iregul = regularization parameter
    #
    #
    #   nc_mc[i] = number of MC events generated in the i-th cause bin 
    #
    #   pce[i][j]   : P(C_i|E_j) : conditional probabilities of the
    #                 i-th cell of the true values to come from j-th 
    #                 cell of the measured values
    #
    #   pc[i]      :  Input: initial probability of the i-th true value cell
    #                 Output : final probability of the i-th true value cell
    #
    #   nc[i]     : estimated number of events in the i-th cell of the 
    #               true values
    #
    #   eff[i]    : efficiency of detecting the cause Ci in any of the possible effects
    #
    #   nc1[i]    : temporary estimated number of events
    #
    #   M_unf[i][j]: unfolding matrix
    #
    #   ibkg = bkg flag (set to a value !=0 if bkg is included in the last cause cell)
    #
    #   errcal = error calculation flag (set = 0 if no error calculation is required)
    #

    # Logarithmic smoothing histogram

    ymin = TMath.Log10(unfsmooth.GetBinLowEdge(1))
    ymax = TMath.Log10(unfsmooth.GetBinLowEdge(noc)+unfsmooth.GetBinWidth(noc))
    unfsmoothlog = TH1F("unfsmoothlog", "unfsmoothlog", noc, ymin, ymax)
    unfsmooth2 = unfsmooth.Clone("unfsmooth2")

    # Definitions and initializations
    
    pc  = ary.array('f',range(noc))
    eff = ary.array('f',range(noc))
    nc1 = ary.array('f',range(noc))

    pce = [[0.0 for col in range(noe)]for raw in range(noc)]
    M_unf = [[0.0 for col in range(noe)]for raw in range(noc)]
    
    # Additional variables used for speedup of Vc1 calculation
    
    nc_inv_mc = ary.array('f',range(noc))
    npec_inv = [[0.0 for col in range(noc)]for raw in range(noe)]
    M_tmp = [[0.0 for col in range(noe)]for raw in range(noe)]

    
    xfit = ary.array('f',range(noc))
    yfit = ary.array('f',range(noc))
    coeff = ary.array('f',range(10))
    coeff1 = ary.array('f',range(10))
    pc1 = ary.array('d',range(noc))
  
    ier = 0
    ntimes = 1

    # New variables for increasing the speed of calculations

    temp = noe*noc
    jevec = ary.array('i',range(temp))
    icvec = ary.array('i',range(temp))

    print "\n BAYESIAN UNFOLDING \n"
    if (iregul != 0):
        if (iregul==1):
            print "REGULARIZATION METHOD CHOSEN: POWER LAW FIT"
        elif (iregul==2):
            print "REGULARIZATION METHOD CHOSEN: POWER + EXPONENTIAL CUT-OFF LAW FIT" 
        elif (iregul==3):
            print "REGULARIZATION METHOD CHOSEN: BROKEN POWER LAW FIT" 
        elif (iregul==4):
            print "REGULARIZATION METHOD CHOSEN: TH1::SMOOTHARRAY" 
        elif(iregul <1 or iregul > 4):
            print "REGULARIZATION METHOD NOT AVAILABLE: NO REGULARIZATION WILL BE APPLIED!"
            iregul = 0
            
    else:
        print "NO REGULARIZATION CHOSEN \n"        


    # Consider only the bins where pec[j][i] > 0
    # jevec[k] = effect index
    # icvec[k] = cause index
    # In this loop the efficiencies are also initialized
    
    nptot = 0
    for i in range(noc):
        eff[i] = 0.
        for j in range(noe):
            if (pec[j][i]>0):
                jevec[nptot] = j
                icvec[nptot] = i
                nptot = nptot +1
                
    print " ===> noe * noc = " + str(noe*noc) + " ===> nptot = " + str(nptot)

    # Count the total observed events and number of bins with at least an observed event

    Nobs = 0
    NNe = 0
    for j in range(noe):
        Nobs = Nobs + ne[j]
        if (ne[j] != 0):
            NNe = NNe + 1
            
    print " Number of observed events = " + str(Nobs)
    print " Number of bins with observed events entries > 0 = " + str(NNe)
    if (Nobs <= 0):
        print " No observed events!"
        ier = 1
        return ier


    # Set here the delta_chi2 to stop the iterations

    if (NNe <= 0):
        NNe = noe
    if (NNe > 10):
        delta_chi2 = 0.01 * (NNe-1)
    else:
        delta_chi2 = 0.01 * (NNe-1)
        
    print "Delta chi square = " + str(delta_chi2)


    # Evaluate here the efficiencies 

    for k in range(nptot):
        eff[icvec[k]] = eff[icvec[k]] + pec[jevec[k]][icvec[k]]
        

    # First iteration    
        
    iter = 0
    chi2 = 1.0e6

    # Continue iterations up to a maximum of 1000

    #while ( (iter<999) and (chi2>delta_chi2)):
    while ( (iter<5) and (chi2>delta_chi2)):

        # Starting

        # Count total number of causes
        nct = 0
        for i in range(noc):
            nct = nct + nc[i]

        # Evaluate probabilities of each cause
        for i in range(noc):
            pc[i] = nc[i] / nct

        # Evaluate total probability and total efficiency
        # Re-evaluate probabilities and number of causes
    
        eff_0 = 0.
        ptot = 0.
        for i in range(noc):
            ptot = ptot + pc[i]
            eff_0 = eff_0 + eff[i] * pc[i]
            nc[i] = pc[i] * Nobs

        if (ptot < 0.999 or ptot > 1.001): 
            print " Un-normalized initial probabilities: ptot = " + str(ptot)
            ier = 8
            return
        else:
            eff_0 = eff_0/ptot
            if (eff_0 < 0.0):
                print "Initial efficiency = " + str(eff_0)
                ier = 9
                return

        for i in range(noc):
            nc[i] = pc[i]*Nobs/eff_0


        # Iterations
        # print "iter: " + str(iter) + " pc[0]= " + str(pc[0])

        # Regularization of the cause distribution

        if (iregul == 1 and iter>0):
            # Power law regularization
            lfirst = 0
            for i in range(noc):
                unfsmooth.SetBinContent(i+1, pc[i]/unfsmooth.GetBinWidth(i+1) )
                xfit[i] = unfsmooth.GetBinCenter(i+1)
                if (pc[i]>0):
                    unfsmoothlog.SetBinContent(i+1, TMath.Log10(pc[i]/unfsmooth.GetBinWidth(i+1)) )
                    ilast = i
                    if(lfirst == 0):
                        ifirst = i
                        lfirst = 1
                    
            x1 = unfsmooth.GetBinLowEdge(ifirst+1)
            y1 = unfsmoothlog.GetBinLowEdge(ifirst+1)
            if (ibkg==0):
                # ilast = noc+1
                ilast += 1
                nlast = noc
            #else:
            #    ilast = noc
            #    nlast = noc-1
                
            x2 = unfsmooth.GetBinLowEdge(ilast+1)
            y2 = unfsmoothlog.GetBinLowEdge(ilast+1)
            
            #print ifirst, ilast, x1, x2, unfsmooth.GetBinContent(ifirst+1), unfsmooth.GetBinContent(ilast), pc[ifirst], pc[ilast-1]
            f1 = TF1("f1", PowerLaw(), x1, x2, 2)
            g1 = TF1("g1", Linear(), y1, y2, 2)
            gamma = 2.7
            p0 = (gamma-1.) * TMath.Power(x1, gamma-1.)
            p0 = unfsmooth.GetBinContent(ifirst+1)
            p0 = 0.
            #print p0, gamma
            f1.SetParameters(p0, gamma)
            lp0 = TMath.Log10(p0)
            lp0 = unfsmoothlog.GetBinContent(ifirst+1)
            lp0 = 0.
            g1.SetParameters(lp0, gamma-1)

            #print "f1 parameters ", f1.GetParameter(0),  f1.GetParameter(1)
            #print "g1 parameters ", g1.GetParameter(0),  g1.GetParameter(1)

            unfsmooth.Fit("f1","W0NQ")
            unfsmoothlog.Fit("g1","W0NQ")
            
            for i in range(2):
                coeff[i] = f1.GetParameter(i)

            #print coeff[0], coeff[1]
            coeff1[0] = TMath.Power(10.0, g1.GetParameter(0))
            coeff1[1] = g1.GetParameter(1)            
            #print coeff[0], coeff[1]


            unfsmooth.Reset()
            unfsmoothlog.Reset()
            
            ptot = 0.
            for i in range (ifirst, nlast):
                yfit[i] = coeff[0]*TMath.Power(xfit[i], -coeff[1]) * unfsmooth.GetBinWidth(i+1)
                ptot = ptot + yfit[i]
                
            if (ibkg != 0):
                ptot = ptot + pc[noc-1]

            for i in range(ifirst, nlast):
                pc[i] = yfit[i]/ptot
                

        elif (iregul == 2 and iter > 0):
            # Power law + exponential cutoff regularization
            for i in range(noc):
                unfsmooth.SetBinContent(i+1, pc[i]/unfsmooth.GetBinWidth(i+1) )
                xfit[i] = unfsmooth.GetBinCenter(i+1)
                if (pc[i]>0):
                    unfsmoothlog.SetBinContent(i+1, TMath.Log10(pc[i]/unfsmooth.GetBinWidth(i+1)) )
                    
            x1 = unfsmooth.GetBinLowEdge(1)
            y1 = unfsmoothlog.GetBinLowEdge(1)
            if (ibkg==0):
                ilast = noc+1
                nlast = noc
            else:
                ilast = noc
                nlast = noc-1
                
            x2 = unfsmooth.GetBinLowEdge(ilast)
            y2 = unfsmoothlog.GetBinLowEdge(ilast)
                
            f1 = TF1("f1", PowerLawExpCut(), x1, x2, 3)
            g1 = TF1("g1", LinearExpCut(), y1, y2, 3)
 
            unfsmooth.Fit("f1","W0NQ")
            unfsmoothlog.Fit("g1","W0NQ")
            
            for i in range(3):
                coeff[i] = f1.GetParameter(i)

            coeff1[0] = TMath.Power(10.0, g1.GetParameter(0))
            coeff1[1] = g1.GetParameter(1)            
            coeff1[2] = 1.0/ (g1.GetParameter(2) * TMath.Log(10.0))             

            unfsmooth.Reset()
            unfsmoothlog.Reset()
            
            ptot = 0.
            for i in range (nlast):
                yfit[i] = coeff1[0]*TMath.Power(xfit[i],-coeff1[1])*TMath.Exp(-xfit[i]/coeff1[2]) * unfsmooth.GetBinWidth(i+1)
                ptot = ptot + yfit[i]
                
            if (ibkg != 0):
                ptot = ptot + pc[noc-1]

            for i in range(nlast):
                pc[i] = yfit[i]/ptot

        elif (iregul == 3 and iter> 0):
            # Broken Power law 
            for i in range(noc):
                unfsmooth.SetBinContent(i+1, pc[i]/unfsmooth.GetBinWidth(i+1) )
                xfit[i] = unfsmooth.GetBinCenter(i+1)
                if (pc[i]>0):
                    unfsmoothlog.SetBinContent(i+1, TMath.Log10(pc[i]/unfsmooth.GetBinWidth(i+1)) )
                    
            x1 = unfsmooth.GetBinLowEdge(1)
            y1 = unfsmoothlog.GetBinLowEdge(1)
            if (ibkg==0):
                ilast = noc+1
                nlast = noc
            else:
                ilast = noc
                nlast = noc-1
                
            x2 = unfsmooth.GetBinLowEdge(ilast)
            y2 = unfsmoothlog.GetBinLowEdge(ilast)
                
            f1 = TF1("f1", BrokenPowerLaw(), x1, x2, 4)
            g1 = TF1("g1", LinearBrokenPowerLaw(), y1, y2, 4)
 
            unfsmooth.Fit("f1","W0NQ")
            unfsmoothlog.Fit("g1","W0NQ")
            
            for i in range(4):
                coeff[i] = f1.GetParameter(i)

            coeff1[0] = TMath.Power(10.0, g1.GetParameter(0))
            coeff1[1] = g1.GetParameter(1)            
            coeff1[2] = TMath.Power(10.0, g1.GetParameter(3))             
            coeff1[3] = g1.GetParameter(2)            

            unfsmooth.Reset()
            unfsmoothlog.Reset()
            
            ptot = 0.
            for i in range (nlast):
                if (xfit[i] < coeff[2]):                    
                    yfit[i] = coeff[0]*TMath.Power(xfit[i],-coeff[1]) * unfsmooth.GetBinWidth(i+1)
                else:
                    yfit[i] = coeff[0]*TMath.Power(xfit[i],-coeff[3]) * unfsmooth.GetBinWidth(i+1)                   
                ptot = ptot + yfit[i]
                
            if (ibkg != 0):
                ptot = ptot + pc[noc-1]

            for i in range(nlast):
                pc[i] = yfit[i]/ptot

        elif (iregul == 4):
            # SMOOTHARRAY regularization
            ifirst = -1
            ilast1 = noc
            if (ibkg==1):
                ilast = noc-3
            i1 = 0
            for i in range (ilast1):
                unfsmooth.SetBinContent(i+1,pc[i]/unfsmooth.GetBinWidth(i+1))
                if (pc[i]> 0.):
                    if (ifirst ==-1):
                        ifirst = i
                    pc1[i1] = pc[i] / unfsmooth.GetBinWidth(i+1)
                    i1 = i1+1

                if (iter==0):
                    if ( (nc[i]>0.01) or  (eff[i]>0)):
                        ilast = i
                else:
                    if (nc[i]>0.01):
                        ilast = i

            nbins = ilast - ifirst + 1
                
            TH1.SmoothArray(nbins, pc1, ntimes)

            ptot = 0.
            for i in range(nbins):               
                ptot += pc1[i]*unfsmooth.GetBinWidth(i+1)

            for i in range(nbins):
                pc[i+ifirst] = pc1[i]*unfsmooth.GetBinWidth(i+1)/ptot

            for i in range(noc):
                unfsmooth2.SetBinContent(i+1,pc[i]/unfsmooth2.GetBinWidth(i+1))
 

        # Here starts the inversion

        # Application of Bayes'theorem to invert probabilities:
        # the pce[i][j] are obtained from the pec[i][j]

        for j in range(noe):
            denom = 0.
            for i in range(noc):
                denom += pec[j][i]*pc[i]
                
            if(denom > 1.E-9):
                denom = 1./denom
            else:
                denom = 0.

            for i in range(noc):
                pce[i][j] = pec[j][i]*pc[i]*denom


        # Evaluate unfolded distribution nc1[i] and true number of events Ntrue

        for i in range(noc):
            nc1[i] = 0.

        Ntrue = 0.
        for k in range(nptot):
            M_unf[icvec[k]][jevec[k]] = 0.
            if (eff[icvec[k]]>0):
                M_unf[icvec[k]][jevec[k]] = pce[icvec[k]][jevec[k]] / eff[icvec[k]]

            nc1[icvec[k]] = nc1[icvec[k]] + M_unf[icvec[k]][jevec[k]]*ne[jevec[k]]
            Ntrue = Ntrue + M_unf[icvec[k]][jevec[k]]*ne[jevec[k]]


        # Compare the initial and the final distribution

        chi2 = 0.    
        for i in range(noc):
            if ((nc[i] + nc1[i])>1.0):
                chi2 = chi2 + TMath.Power( (nc1[i]-nc[i]), 2) / (nc[i]+nc1[i]) 
            else:
                chi2 = chi2 + TMath.Power( (nc1[i]-nc[i]), 2) 

        # Increment iteration counter

    	print " iter = ",str(iter)
    	print " chi2 = ",str(chi2)

        iter = iter+1

        # Replace old cause distribution nc[i] with new one nc1[i] and evaluate new cause probabilities pc[i]

        for i in range(noc):
            nc[i] = 0
            if( nc1[i]>0 and TMath.IsNaN(nc1[i])==0):
                nc[i] = nc1[i]
                pc[i] = nc1[i]/Ntrue
                    
    ## end of iterations ##
    ## rescaling wtih maxeff
    Ntrue = 0.
    for i in range(noc):
        nc[i] *= maxeff
        eff[i] *= 1./maxeff
        Ntrue += nc[i]
        for j in range(noe):
            M_unf[i][j] *= maxeff
            pec[j][i]  *= 1./maxeff
 

    for i in range(noc):
        pc[i] = nc[i]/Ntrue

    eff_true = Nobs / Ntrue
    
    print " ## Bayes unfolding terminated  ## " 
    print "    Number of iterations      = ",str(iter)
    print "    Observed number of events = ",str(Nobs)
    print "    Unfolded number of events = ",str(Ntrue)
    print "    Initial efficiency        = ",str(eff_0)
    print "    Final efficiency          = ",str(eff_true)
    print "    Final chi2                = ",str(chi2)
    
    if (errcal == 0):
        print "No error calculation "
        return ier

    # Start here error calculation

    # Vc0 = covariance matrix from uncertainties on observed distributions
    # Vc1 = covariance matrix from uncertainties on MC simulation used to evaluate the pec[i][j]
  
    for i in range(noc):
        for j in range(noc):
            Vc0[i][j] = 0.
            Vc1[i][j] = 0.
            # adding the normalization uncertainty to the contribution from the observed data       
            Vc0[i][j] = nc[i]*nc[j]/Nobs
          

    # errors from observed number of events 

    # Multinomial approximation (not correct)

##    for k in range (noc):
##        for l in range (k,noc):
##            for i in range (noe):
##                for j in range (noe):
##                    if (i==j):
##                        Vc0[k][l] = Vc0[k][l] + M_unf[k][i]*M_unf[l][j]*neeq[i]*(1. - neeq[j]/Ntrue)
##                    else:
##                        Vc0[k][l] = Vc0[k][l] + M_unf[k][i]*M_unf[l][j]*neeq[i]*neeq[j]/Ntrue            

    # Poisson approximation

    for k in range(noc):
        for l in range(noc):
            for i in range(noe):
                j = i
                Vc0[k][l] = Vc0[k][l] + M_unf[k][i]*M_unf[l][j]*neeq[i]

                        
    # Take diagonal elements of Vc0 as errors enc[i]

    for i in range(noc):
        enc[i] = 0.
        if (Vc0[i][i] >= 0.):
            enc[i] = TMath.Sqrt(Vc0[i][i])
        else:
            print "Warning: coavariance matrix Vc0 has negative value on the diagonal!"
            ier = 2
    
    if (errcal == 1):
        print "No error calculation from smearing function"
        return ier
  

    # errors from smearing matrix     
  
    # calculate first some temporary quantities to speed up
    # processing inside the main 4-fold nested do-loop
    
    for k in range(noc):
        nc_inv_mc[k] = 0.
        if (nc_mc[k] !=0):
            nc_inv_mc[k] = 1./nc_mc[k]

        for i in range(noe):
            npec_inv[i][k] = 0.
            if (pec[i][k] != 0):
                npec_inv[i][k] = nc_inv_mc[k] / pec[i][k]
      
    for i in range(noe):
        for j in range (noe):
            M_tmp[i][j] = 0.
          
    for i in range(nptot):
        for j in range (nptot):
            if (icvec[i] == icvec[j]):
                denom = - nc_inv_mc[icvec[i]]
                if (jevec[i] == jevec[j]):
                    denom = denom + npec_inv[jevec[i]][icvec[i]]
                    M_tmp[jevec[i]][jevec[j]] = M_tmp[jevec[i]][jevec[j]] + M_unf[icvec[i]][jevec[i]]*M_unf[icvec[j]][jevec[j]]*eff[icvec[i]]*eff[icvec[j]]*denom


    for i1 in range(nptot):
        for i2 in range(nptot):
            k = icvec[i1]
            i = jevec[i1]
            l = icvec[i2]
            j = jevec[i2]
            denom = M_unf[l][i]*nc_inv_mc[l] + M_unf[k][j]*nc_inv_mc[k] + M_tmp[i][j]
            if (k==l):
                denom = denom - nc_inv_mc[k] / eff[k]
                if (i==j):
                    denom = denom + npec_inv[i][k]
                  
            if (i==j):
                denom = denom - M_unf[l][i]*eff[l]*npec_inv[i][l]-M_unf[k][i]*eff[k]*npec_inv[i][k] 
      
            Vc1[k][l] = Vc1[k][l] + M_unf[k][i]*M_unf[l][j]*denom*ne[i]*ne[j]
    
  

    for i in range(noc):
        enc1[i] = 0.
        if ( (Vc1[i][i] ) >= 0. ):
            enc1[i] = TMath.Sqrt(Vc1[i][i])
        else: 
            print  "Warning: covariance matrix Vc1 has negative value on the diagonal!" 
            ier = 3

            
    return ier


# end definition of fastbayes function



