"""
@brief: dirty skimmer, starting from Stephan's skim -> use Ions & Co
@vale: used for skimming of protons/ions

"""

from sys import argv
import numpy as np
import math
#from xrd_utils import isfile
from os.path import isdir, abspath
from os import remove, mkdir, system as os_call
from glob import glob
from argparse import ArgumentParser
from itertools import product as iproduct
from array import array

# ML tools allocation
from dampeapi import DampeApi
from dampe_bgo_sat_api import DampeBgoApi

from ROOT import *

from datetime import datetime
start=datetime.now()

#mc_pdg_particle_id = 100010 # Neon
e_min = 100. # GeV
#delta_theta_max = 25. # Degrees

def mc_particle_info(mc_particle_name):
    particles = {"p":[2212, 1],"He":[100002, 2], "C":[100006, 6], "N":[100007, 7], "O":[100008, 8], "F":[100009, 9], "Ne":[100010, 10], "Na":[100011, 11], "Mg":[100012, 12], "Al":[100013, 13], "Si":[100014, 14], "P":[100015, 15], "S":[100016, 16], "Cl":[100017, 17], "Ar":[100018, 18], "K":[100019, 19], "Ca":[100020, 20], "Fe":[100026, 26]}
    try:
        particles[mc_particle_name]
    except:
        print "Error! Wrong particle name given in arguments, check among those available\n"
        pass

    return particles[mc_particle_name]

def spectral_weight(E, gamma, edges):
    # weight for a poer-law spectrum normalized in [Emin, Emax] from edges
    bin_idx = np.searchsorted(edges, E, side='right') - 1
    e1 = edges[bin_idx]
    e2 = edges[bin_idx + 1]
    norm = (gamma-1.) * math.log(e2/e1) / (e1**(1.-gamma) - e2**(1.-gamma))
    return norm * E**(1.-gamma)

def main(args=None):

    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="run dirty skimmer (for e/gamma)")
    parser.add_argument("-i","--input",dest='input', help='the input files to be processed ')
    parser.add_argument("-v","--verbose", action='store_true', default=False, dest='verbose', help='run in high verbosity mode')
    parser.add_argument("-q","--quiet", action='store_true', default=False, dest='quiet', help='suppress a lot of output, quiet mode')
    #parser.add_argument("-f","--fileMode",default=False, action="store_true", dest='FileMode', help='use this flag to pass textfile with file names to tool')
    parser.add_argument("-o","--diagnostics_output",default="diagnostics.root", type=str, dest='diag_output', help='name of diagnostics file to be written')
    parser.add_argument("-O","--output",default="skimmed_data/", type=str, dest='output', help='name of output directory')
    parser.add_argument("-s","--skimFile",default="skim_out.root", type=str, dest='skimFile', help='skim files for roofit')
    parser.add_argument("-e","--text_output_missing",default="livetime_missing.txt", type=str, dest='txt_output_missing', help='name of text output file with missing runs')
    parser.add_argument("-mc",default=False, action="store_true", dest='mc', help='use this flag for MC data')
    parser.add_argument("-pid", "--mc_particle_name", default="p", type=str, dest="mc_particle_name", help="particle simulated in the MC file you want to skim, e.g.[p, He, C, N, O, F, Ne, na, Mg, Al, Si, Fe]")
    parser.add_argument("-data",default=False, action="store_true", dest='data', help='use data and output .root ')
    parser.add_argument("-skim",default=False, action="store_true", dest='skim', help='active for skim')
    opts = parser.parse_args(args)


    verbosity = 2002 # ignore ROOT warnings
    if opts.verbose:
        verbosity = 1000
    if opts.quiet:
        verbosity = 3003

    mc_pdg_particle_id = mc_particle_info(opts.mc_particle_name)[0]
    mc_particle_Z = mc_particle_info(opts.mc_particle_name)[1]
    print 'mc_pdg_particle_id ', mc_pdg_particle_id, ', mc_particle_Z ', mc_particle_Z

    inp = opts.input
    print '====> reading files in  ... ', inp
    stream_name = "ions"

    oout_text_miss = opts.txt_output_missing
    oout_text_skim = opts.skimFile
    ######### log bins defintion #######
    nEDBINs_6bins = 36 #6 bins per 6 decadi
    EDmin0=10
    EDmax0=10000000
    EDmin_6bins = []
    EDmax_6bins = []
    EDEdge_6bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_6bins
    for iedge in range(0, nEDBINs_6bins):
        EDEdge_6bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_6bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_6bins.append(EDmax_6bins[nEDBINs_6bins-1])
    Edges_6bins= array('d',EDEdge_6bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_3bins = 18 # 3 bins per 6 decadi
    EDmin_3bins = []
    EDmax_3bins = []
    EDEdge_3bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_3bins
    for iedge in range(0, nEDBINs_3bins):
        EDEdge_3bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_3bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_3bins.append(EDmax_3bins[nEDBINs_3bins-1])
    Edges_3bins= array('d',EDEdge_3bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_7bins = 42 # 7 bins per 6 decadi
    EDmin_7bins = []
    EDmax_7bins = []
    EDEdge_7bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_7bins
    for iedge in range(0, nEDBINs_7bins):
        EDEdge_7bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_7bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_7bins.append(EDmax_7bins[nEDBINs_7bins-1])
    Edges_7bins= array('d',EDEdge_7bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_5bins = 30 # 5 bins per 6 decadi
    EDmin_5bins = []
    EDmax_5bins = []
    EDEdge_5bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_5bins
    for iedge in range(0, nEDBINs_5bins):
        EDEdge_5bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_5bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_5bins.append(EDmax_5bins[nEDBINs_5bins-1])
    Edges_5bins= array('d',EDEdge_5bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_4bins = 24 # 4 bins per 6 decadi
    EDmin_4bins = []
    EDmax_4bins = []
    EDEdge_4bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_4bins
    for iedge in range(0, nEDBINs_4bins):
        EDEdge_4bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_4bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_4bins.append(EDmax_4bins[nEDBINs_4bins-1])
    Edges_4bins= array('d',EDEdge_4bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_8bins = 48 # 8 bins per 6 decadi
    EDmin_8bins = []
    EDmax_8bins = []
    EDEdge_8bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_8bins
    for iedge in range(0, nEDBINs_8bins):
        EDEdge_8bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_8bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_8bins.append(EDmax_8bins[nEDBINs_8bins-1])
    Edges_8bins= array('d',EDEdge_8bins) # this make a bound array for TH1F

    ######### log bins defintion #######
    nEDBINs_120bins = 720 # 120 bins per 6 decadi
    EDmin_120bins = []
    EDmax_120bins = []
    EDEdge_120bins = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs_120bins
    for iedge in range(0, nEDBINs_120bins):
        EDEdge_120bins.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax_120bins.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        
    EDEdge_120bins.append(EDmax_120bins[nEDBINs_120bins-1])
    Edges_120bins= array('d',EDEdge_120bins) # this make a bound array for TH1F
    
    BGOz = [58.5,87.5,116.5,145.5,174.5,203.5,232.5,261.5,290.5,319.5,348.5,377.5,406.5,435.5]

    BGOz_X = BGOz[1::2]  
    BGOz_Y  =BGOz[0::2]   

    spectral_indices = { 
        "E2e4": 2.4, 
        "E2e5": 2.5, 
        "E2e6": 2.6, 
        "E2e7": 2.7, 
        "E2e8": 2.8, 
        "E2e9": 2.9, 
        "E3e0": 3.0, 
    }
    energy_binnings = { 
        "3bins":  (nEDBINs_3bins,  Edges_3bins),
        "4bins":  (nEDBINs_4bins,  Edges_4bins),
        "5bins":  (nEDBINs_5bins,  Edges_5bins),
        "6bins":  (nEDBINs_6bins,  Edges_6bins),
        "7bins":  (nEDBINs_7bins,  Edges_7bins),
        "8bins":  (nEDBINs_8bins,  Edges_8bins),
        "120bins":(nEDBINs_120bins,Edges_120bins),
    }
    
    
    pevspace_api = DampeApi()
    #check_api = DampeApi()
    bgo_api = DampeBgoApi()
    
    pevspace_api_quench = DampeApi()
    bgo_api_quench = DampeBgoApi()
    
    gSystem.Load("libDmpEvent.so")
    gSystem.Load("libDmpEventFilter.so")
    
    gSystem.Load("libDmpKernel.so")
    gSystem.Load("libDmpService.so")

    gROOT.SetBatch(True)
    gROOT.ProcessLine( "gErrorIgnoreLevel = %i;"%verbosity)
    import DMPSW
    
    ############## create skim file ############
    root_dir = abspath(opts.output)
    print " skimming files "
    print '* * will create output here: %s'%root_dir
    if not isdir(root_dir): mkdir(root_dir)
    
    file_skim_out = root_dir+'/'+oout_text_skim
    print file_skim_out
    tf_skim = TFile(file_skim_out,"RECREATE") 
    newtree = TTree("newtree","test")
    
    nlayer_psd = 2
    nlayer_bgo = 14
 
    #### PSD
    fPSD_EnergyY = array('d', 2*[ -999.] )
    fPSD_EnergyX = array('d', 2*[ -999.] )
    newtree.Branch('PSD_EnergyY', fPSD_EnergyY, 'PSD_EnergyY[2]/D')
    newtree.Branch('PSD_EnergyX', fPSD_EnergyX, 'PSD_EnergyX[2]/D')
    
    fPSD_EnergyY_corr = array('d', 2*[ -999.] )
    fPSD_EnergyX_corr = array('d', 2*[ -999.] )
    newtree.Branch('PSD_EnergyY_corr', fPSD_EnergyY_corr, 'PSD_EnergyY_corr[2]/D')
    newtree.Branch('PSD_EnergyX_corr', fPSD_EnergyX_corr, 'PSD_EnergyX_corr[2]/D')

    fPSD_psdX_projY = array('d', 2*[ -999. ] )
    fPSD_psdY_projX = array('d', 2*[ -999. ] )
    newtree.Branch("PSD_psdX_projY",fPSD_psdX_projY,"PSD_psdX_projY[2]/D")
    newtree.Branch("PSD_psdY_projX",fPSD_psdY_projX,"PSD_psdY_projX[2]/D")
    
    fPSD_psdX_position  = array('d', 2*[ -999. ] )
    fPSD_psdY_position  = array('d', 2*[ -999. ] )
    newtree.Branch("PSD_psdX_position",fPSD_psdX_position,"PSD_psdX_position[2]/D")
    newtree.Branch("PSD_psdY_position",fPSD_psdY_position,"PSD_psdY_position[2]/D")

    fPSD_psdX_pathlength = array('d', 2*[ -999.] )
    fPSD_psdY_pathlength = array('d', 2*[ -999.] )
    newtree.Branch("PSD_psdX_pathlength",fPSD_psdX_pathlength,"PSD_psdX_pathlength[2]/D")
    newtree.Branch("PSD_psdY_pathlength",fPSD_psdY_pathlength,"PSD_psdY_pathlength[2]/D")
    
    fPSD_psdX_Emax_layer     = array('d', [ -999.] )
    fPSD_psdY_Emax_layer     = array('d', [ -999.] )
    newtree.Branch("PSD_psdX_Emax_layer",fPSD_psdX_Emax_layer,"PSD_psdX_Emax_layer/D")
    newtree.Branch("PSD_psdY_Emax_layer",fPSD_psdY_Emax_layer,"PSD_psdY_Emax_layer/D")
    
    ### MARGHERITA&ANTONIO - PSD_Emax (Oct 23, 2017)
    fPSD_EmaxY = array('d', 2*[ -999.] )
    fPSD_EmaxX = array('d', 2*[ -999.] )
    newtree.Branch('PSD_EmaxY', fPSD_EmaxY, 'PSD_EmaxY[2]/D')
    newtree.Branch('PSD_EmaxX', fPSD_EmaxX, 'PSD_EmaxX[2]/D')

    ### MARGHERITA&ANTONIO - PSD_Etracks (Oct 13, 2017)
    fPSD_cpsdmax    = array('d', 4*[-999.])
    fPSD_cpsdstk    = array('d', 4*[-999.])
    newtree.Branch("PSD_cpsdmax",fPSD_cpsdmax,"PSD_cpsdmax[4]/D")
    newtree.Branch("PSD_cpsdstk",fPSD_cpsdstk,"PSD_cpsdstk[4]/D")

    fPSD_Etrack01   = array('d', [-999.])       
    fPSD_Etrack23   = array('d', [-999.])   
    newtree.Branch("PSD_Etrack01",fPSD_Etrack01,"PSD_Etrack01/D")
    newtree.Branch("PSD_Etrack23",fPSD_Etrack23,"PSD_Etrack23/D")

    ### GSSI - PSD Charge corrected after trunk update 8989
    fpsdchargeX0 = array('d', [-999.])
    fpsdchargeX1 = array('d', [-999.])
    newtree.Branch("PSD_ChargeX0", fpsdchargeX0, "PSD_ChargeX0/D")
    newtree.Branch("PSD_ChargeX1", fpsdchargeX1, "PSD_ChargeX1/D")

    fpsdchargeY0 = array('d', [-999.])
    fpsdchargeY1 = array('d', [-999.])
    newtree.Branch("PSD_ChargeY0", fpsdchargeY0, "PSD_ChargeY0/D")
    newtree.Branch("PSD_ChargeY1", fpsdchargeY1, "PSD_ChargeY1/D")

    fPSD_Global_Charge = array('d', [-999.])
    newtree.Branch("PSD_Global_Charge", fPSD_Global_Charge, "PSD_Global_Charge/D")

    fPSD_Global_Charge_PathAverage = array('d', [-999.])
    newtree.Branch("PSD_Global_Charge_PathAverage", fPSD_Global_Charge_PathAverage, "PSD_Global_Charge_PathAverage/D")

    fPSD_PathWeighted_Charge = array('d', [-999.])
    newtree.Branch("PSD_PathWeighted_Charge", fPSD_PathWeighted_Charge, "PSD_PathWeighted_Charge/D")

    fPSD_nhits = array('i', [ 0 ] )
    newtree.Branch("PSD_nhits",fPSD_nhits,"PSD_nhits/I")
    
    #### STK
    fSTK_chargeX = array('d', 6*[ -999. ] )
    fSTK_chargeY = array('d', 6*[ -999. ] )
    newtree.Branch("STK_chargeX",fSTK_chargeX,"STK_chargeX[6]/D")
    newtree.Branch("STK_chargeY",fSTK_chargeY,"STK_chargeY[6]/D")

    fSTK_chargeX_etaCorr = array('d', 6*[ -999. ] )
    fSTK_chargeY_etaCorr = array('d', 6*[ -999. ] )
    newtree.Branch("STK_chargeX_etaCorr",fSTK_chargeX_etaCorr,"STK_chargeX_etaCorr[6]/D")
    newtree.Branch("STK_chargeY_etaCorr",fSTK_chargeY_etaCorr,"STK_chargeY_etaCorr[6]/D")

    fSTK_ntracks = array('i', [ -999 ] )
    fSTK_nclusters = array('i', [ -999 ] )
    fSTK_theta_correction = array('d', [ -999. ] )
    newtree.Branch("STK_ntracks",fSTK_ntracks,"STK_ntracks/I")
    newtree.Branch("STK_nclusters",fSTK_nclusters,"STK_nclusters/I")
    newtree.Branch("STK_theta_correction",fSTK_theta_correction,"STK_theta_correction/D")

    fSTK_trackIP = array('d', 3*[ -999. ] )
    fSTK_trackDirection = array('d', 2*[ -999. ] ) 
    fSTK_vertexPrediction = array('d', [ -999. ] )
    newtree.Branch("STK_trackIP",fSTK_trackIP,"STK_trackIP[3]/D")
    newtree.Branch("STK_trackDirection",fSTK_trackDirection,"STK_trackDirection[2]/D")
    newtree.Branch("STK_vertexPrediction",fSTK_vertexPrediction,"STK_vertexPrediction/D")

    #### BGO
    fBGO_EnergyG = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG",fBGO_EnergyG,"BGO_EnergyG/D")

    fBGO_EnergyG_SatCorrp = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorrp",fBGO_EnergyG_SatCorrp,"BGO_EnergyG_SatCorrp/D")
    fBGO_EnergyG_SatCorrHe = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorrHe",fBGO_EnergyG_SatCorrHe,"BGO_EnergyG_SatCorrHe/D")
    fBGO_EnergyG_SatCorrBCNO = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorrBCNO",fBGO_EnergyG_SatCorrBCNO,"BGO_EnergyG_SatCorrBCNO/D")
    fBGO_EnergyG_QuenchSatCorrBCNO = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorrBCNO",fBGO_EnergyG_QuenchSatCorrBCNO,"BGO_EnergyG_QuenchSatCorrBCNO/D")
    fBGO_EnergyG_SatCorrNeMgSi = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorrNeMgSi",fBGO_EnergyG_SatCorrNeMgSi,"BGO_EnergyG_SatCorrNeMgSi/D")
    fBGO_EnergyG_QuenchSatCorr = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorr",fBGO_EnergyG_QuenchSatCorr,"BGO_EnergyG_QuenchSatCorr/D")
    fBGO_EnergyG_SatCorrFe = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorrFe",fBGO_EnergyG_SatCorrFe,"BGO_EnergyG_SatCorrFe/D")

    fBGO_EnergyG_SatCorr_ML = array('d', 29*[ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorr_ML",fBGO_EnergyG_SatCorr_ML,"BGO_EnergyG_SatCorr_ML[29]/D")
    fBGO_EnergyG_QuenchSatCorr_ML = array('d', 29*[ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorr_ML",fBGO_EnergyG_QuenchSatCorr_ML,"BGO_EnergyG_QuenchSatCorr_ML[29]/D")

    fBGO_EnergyG_SatCorr_ML_He = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorr_ML_He",fBGO_EnergyG_SatCorr_ML_He,"BGO_EnergyG_SatCorr_ML_He/D")
    fBGO_EnergyG_SatCorr_ML_Fe = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorr_ML_Fe",fBGO_EnergyG_SatCorr_ML_Fe,"BGO_EnergyG_SatCorr_ML_Fe/D")
    fBGO_EnergyG_SatCorr_ML_ions = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_SatCorr_ML_ions",fBGO_EnergyG_SatCorr_ML_ions,"BGO_EnergyG_SatCorr_ML_ions/D")
    
    fBGO_EnergyG_QuenchSatCorr_ML_He = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorr_ML_He",fBGO_EnergyG_QuenchSatCorr_ML_He,"BGO_EnergyG_QuenchSatCorr_ML_He/D")
    fBGO_EnergyG_QuenchSatCorr_ML_Fe = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorr_ML_Fe",fBGO_EnergyG_QuenchSatCorr_ML_Fe,"BGO_EnergyG_QuenchSatCorr_ML_Fe/D")
    fBGO_EnergyG_QuenchSatCorr_ML_ions = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG_QuenchSatCorr_ML_ions",fBGO_EnergyG_QuenchSatCorr_ML_ions,"BGO_EnergyG_QuenchSatCorr_ML_ions/D")
    
    fBGO_EnergyG_Quench = array('d', [-999.])
    newtree.Branch("BGO_EnergyG_Quench",fBGO_EnergyG_Quench,"BGO_EnergyG_Quench/D")

    fBGO_HET = array('i', [ 0 ] )
    fBGO_LET = array('i', [ 0 ] )
    fBGO_UBT = array('i', [ 0 ] )
    fBGO_UBE = array('i', [ 0 ] )
    newtree.Branch("BGO_HET",fBGO_HET,"BGO_HET/I")
    newtree.Branch("BGO_LET",fBGO_LET,"BGO_LET/I")
    newtree.Branch("BGO_UBT",fBGO_UBT,"BGO_UBT/I")
    newtree.Branch("BGO_UBE",fBGO_UBE,"BGO_UBE/I")

    fBGO_EneLay = array('d', 14*[ -999.] )
    newtree.Branch("BGO_EneLay",fBGO_EneLay,"BGO_EneLay[14]/D")
    fBGO_nhits_layers = array('i', 14*[ -999 ] )
    newtree.Branch("BGO_nhits_layers",fBGO_nhits_layers,"BGO_nhits_layers[14]/I")

    fBGO_nbars = array('i', [ -999 ] )
    newtree.Branch("BGO_nbars",fBGO_nbars,"BGO_nbars/I")

    fBGO_Xtr = array('d',[-999.])
    newtree.Branch("BGO_xtr",fBGO_Xtr,"BGO_xtr/D")

    ### MARGHERITA&ANTONIO - COORDINATES MAX ENERGY IN BGO LAYERS and BGO/STK (Oct 03, 2017)
    fBGO_cbgomax = array('d',14*[-999.])
    newtree.Branch("BGO_cbgomax",fBGO_cbgomax,"BGO_cbgomax[14]/D")

    fBGO_cbgostk = array('d',14*[-999.])
    newtree.Branch("BGO_cbgostk",fBGO_cbgostk,"BGO_cbgostk[14]/D")

    ### MARGHERITA&ANTONIO - TIME AND TIMEMS (Oct 13, 2017)
    ftime = array('i', [ 0 ] )
    ftimems = array('i', [ 0 ] )
    newtree.Branch("time",ftime,"time/I")
    newtree.Branch("timems",ftimems,"timems/I")

    ### 
    flong = array('d', [-9999.])
    flat  = array('d', [-9999.])
    newtree.Branch("long",flong,"long/D")
    newtree.Branch("lat", flat, "lat/D")

    fSTKtrack_to_PSD_topY = array('d', [-999.] )
    fSTKtrack_to_PSD_topX = array('d', [-999.] )
    fBGOtrack_to_PSD_topY = array('d', [-999.] )
    newtree.Branch("STKtrack_to_PSD_topY",fSTKtrack_to_PSD_topY,"STKtrack_to_PSD_topY/D")
    newtree.Branch("STKtrack_to_PSD_topX",fSTKtrack_to_PSD_topX,"STKtrack_to_PSD_topX/D")
    newtree.Branch("BGOtrack_to_PSD_topY",fBGOtrack_to_PSD_topY,"BGOtrack_to_PSD_topY/D")

    fSTK_theta = array('d', [-999.])
    fBGO_theta = array('d', [-999.])
    newtree.Branch("STK_theta",fSTK_theta,"STK_theta/D")
    newtree.Branch("BGO_theta",fBGO_theta,"BGO_theta/D")

    fSTK_phi = array('d', [-999.])
    fBGO_phi = array('d', [-999.])
    newtree.Branch("STK_phi", fSTK_phi, "STK_phi/D")
    newtree.Branch("BGO_phi", fBGO_phi, "BGO_phi/D")

    fBGO_slopeXZ_analy = array('d', [-999.])
    fBGO_slopeYZ_analy = array('d', [-999.])
    fBGO_interceptXZ_analy = array('d', [-999.])
    fBGO_interceptYZ_analy = array('d', [-999.])
    newtree.Branch("BGO_slopeXZ_analy", fBGO_slopeXZ_analy, "BGO_slopeXZ_analy/D")
    newtree.Branch("BGO_slopeYZ_analy", fBGO_slopeYZ_analy, "BGO_slopeYZ_analy/D")
    newtree.Branch("BGO_interceptXZ_analy", fBGO_interceptXZ_analy, "BGO_interceptXZ_analy/D")
    newtree.Branch("BGO_interceptYZ_analy", fBGO_interceptYZ_analy, "BGO_interceptYZ_analy/D")

    fBGO_theta_analy = array('d', [-999.])
    fBGO_phi_analy = array('d', [-999.])
    newtree.Branch("BGO_theta_analy", fBGO_theta_analy, "BGO_theta_analy/D")
    newtree.Branch("BGO_phi_analy", fBGO_phi_analy, "BGO_phi_analy/D")

    fBGO_slopeX_ML = array('d', [-999.])
    fBGO_slopeY_ML = array('d', [-999.])
    fBGO_interceptX_ML = array('d', [-999.])
    fBGO_interceptY_ML = array('d', [-999.])
    newtree.Branch("BGO_slopeX_ML", fBGO_slopeX_ML, "BGO_slopeX_ML/D")
    newtree.Branch("BGO_slopeY_ML", fBGO_slopeY_ML, "BGO_slopeY_ML/D")
    newtree.Branch("BGO_interceptX_ML", fBGO_interceptX_ML, "BGO_interceptX_ML/D")
    newtree.Branch("BGO_interceptY_ML", fBGO_interceptY_ML, "BGO_interceptY_ML/D")

    fSTK_slopeX_ML = array('d', [-999.])
    fSTK_slopeY_ML = array('d', [-999.])
    fSTK_interceptX_ML = array('d', [-999.])
    fSTK_interceptY_ML = array('d', [-999.])
    fSTK_interceptZ_ML = array('d', [-999.])
    newtree.Branch("STK_slopeX_ML", fSTK_slopeX_ML, "STK_slopeX_ML/D")
    newtree.Branch("STK_slopeY_ML", fSTK_slopeY_ML, "STK_slopeY_ML/D")
    newtree.Branch("STK_interceptX_ML", fSTK_interceptX_ML, "STK_interceptX_ML/D")
    newtree.Branch("STK_interceptY_ML", fSTK_interceptY_ML, "STK_interceptY_ML/D")
    newtree.Branch("STK_interceptZ_ML", fSTK_interceptZ_ML, "STK_interceptZ_ML/D")

    #### MC 
    fMC_EnergyT = array('d', [ -999. ] )
    fMC_theta  = array('d', [ -999. ] )
    fMC_phi = array('d', [-999.] )
    newtree.Branch("MC_EnergyT",fMC_EnergyT,"MC_EnergyT/D")
    newtree.Branch("MC_theta",fMC_theta,"MC_theta/D")
    newtree.Branch("MC_phi",fMC_phi,"MC_phi/D")
    fMC_stopZ = array('d', [ -99999. ] )
    fMC_stopX = array('d', [ -99999. ] )
    fMC_stopY = array('d', [ -99999. ] )
    newtree.Branch("MC_stopZ",fMC_stopZ,"MC_stopZ/D")
    newtree.Branch("MC_stopX",fMC_stopX,"MC_stopX/D")
    newtree.Branch("MC_stopY",fMC_stopY,"MC_stopY/D")

    ###########################################################################

    dmpch = DmpChain("CollectionTree")
    metach = TChain("RunMetadataTree")

    files = [f.replace("\n","") for f in open(inp,'r').readlines()]    
   
    for ifile, f in enumerate(files):
        DMPSW.IOSvc.Set("InData/Read" if ifile == 0 else "InData/ReadMore",f)
        if True:#isfile(f):
            dmpch.Add(f)
            metach.Add(f)
            print ifile , f
        else:
            with open(oout_text_miss,'ar') as file:
                file.write(f+'\n')
                continue
    
    if not opts.mc:
        DMPSW.IOSvc.Set("OutData/NoOutput", "True")
        DMPSW.IOSvc.Initialize()
        pFilter = DmpFilterOrbit("EventHeader")
        pFilter.ActiveMe()
    
    
    runmetadata = DmpRunMetadata()                      # [optional] - only if you use STK and run on the data (needed for the alignment)
    metach.SetBranchAddress("RunMetadata", runmetadata)     # [optional]   ......


    BGOzTop = 46.
    BGOzBot = 448.
    bgo_acceptance_range = 280.

    ############# histo defintion #########
    h_terrestrial_lat_vs_long =  TH2F("h_terrestrial_lat_vs_long","h_terrestrial_lat_vs_long",360,0,360,180,-90,90)
    ncuts = 12 # first 5 cuts are standard, then analysis cuts
    cut_string = ['All',   'etot > 20 GeV ',  'BGO Acceptance', 'EratioLay','BGOLateral','SAA','HET','NTracks','BGO-STK','PSD Fiducial','STK-PSD','NLayers']
    h_bgo_projection_before = TH2F("h_bgo_projection_before","h_bgo_projection_before",1000,-500,500,1000,-500,500)
    h_bgo_projection_after  = TH2F("h_bgo_projection_after", "h_bgo_projection_after",1000,-500,500,1000,-500,500)
    h_bgo_projectionTop_before = TH2F("h_bgo_projectionTop_before","h_bgo_projectionTop_before",1000,-500,500,1000,-500,500)
    h_bgo_projectionTop_after  = TH2F("h_bgo_projectionTop_after", "h_bgo_projectionTop_after",1000,-500,500,1000,-500,500)
    h_bgo_EratioLay_tot_before  = TH1F("h_bgo_EratioLay_tot_before","h_bgo_EratioLay_tot_before",100,0,1)
    h_bgo_EratioLay_goo_before  = TH1F("h_bgo_EratioLay_goo_before","h_bgo_EratioLay_goo_before",100,0,1)
    h_bgo_EratioLay_bad_before  = TH1F("h_bgo_EratioLay_bad_before","h_bgo_EratioLay_bad_before",100,0,1)
    h_bgo_EratioLay_tot_after  = TH1F("h_bgo_EratioLay_tot_after","h_bgo_EratioLay_tot_after",100,0,1)
    h_bgo_lateralpos_after  = TH1F("h_bgo_lateralpos_after","h_bgo_lateralpos_after",25,-1.5,23.5)
    h_bgo_lateralpos_before = TH1F("h_bgo_lateralpos_before","h_bgo_lateralpos_before",25,-1.5,23.5)
    h_stk_ntracks_before = TH1F("h_stk_ntracks_before","h_stk_ntracks_before",1001,-0.5,1000.5)
    h_stk_ntracks_after = TH1F("h_stk_ntracks_after","h_stk_ntracks_after",1001,-0.5,1000.5)
    h_stk_XY_check = TH1F("h_stk_XY_check","1 = cluster X; 2 = cluster Y; 3 = cluster X and Y",5,0.5,5.5)
    h_stk_IP_XY = TH2F("h_stk_IP_XY","h_stk_IP_XY",1000,-500,500,1000,-500,500)
    h_stk_resX_top = TH1F("h_stk_resX_top","h_stk_resX_top",200,-1000,1000)
    h_stk_resY_top = TH1F("h_stk_resY_top","h_stk_resY_top",200,-1000,1000)
    h_stk_chi2_norm = TH1F("h_stk_chi2_norm","h_stk_chi2_norm",1000,0,1000)
    h_stk_resX = TH1F("h_stk_resX","h_stk_resX",200,-1000,1000)
    h_stk_resY = TH1F("h_stk_resY","h_stk_resY",200,-1000,1000)
    h_stk_resX_top_deltaTheta = TH1F("h_stk_resX_top_deltaTheta","h_stk_resX_top_deltaTheta",200,-1000,1000)
    h_stk_resY_top_deltaTheta = TH1F("h_stk_resY_top_deltaTheta","h_stk_resY_top_deltaTheta",200,-1000,1000)
    h_stk_resX_deltaTheta = TH1F("h_stk_resX_deltaTheta","h_stk_resX_deltaTheta",200,-1000,1000)
    h_stk_resY_deltaTheta = TH1F("h_stk_resY_deltaTheta","h_stk_resY_deltaTheta",200,-1000,1000)
    h_stk_resX_top_resBot = TH1F("h_stk_resX_top_resBot","h_stk_resX_top_resBot",200,-1000,1000)
    h_stk_resY_top_resBot = TH1F("h_stk_resY_top_resBot","h_stk_resY_top_resBot",200,-1000,1000)
    h_stk_resX_top_selection = TH1F("h_stk_resX_top_selection","h_stk_resX_top_selection",200,-1000,1000)
    h_stk_resY_top_selection = TH1F("h_stk_resY_top_selection","h_stk_resY_top_selection",200,-1000,1000)
    h_stk_resX_selection = TH1F("h_stk_resX_selection","h_stk_resX_selection",200,-1000,1000)
    h_stk_resY_selection = TH1F("h_stk_resY_selection","h_stk_resY_selection",200,-1000,1000)
    h_stk_clusterX = TH1F("h_stk_clusterX","h_stk_clusterX",10000,0,10000)
    h_stk_clusterY = TH1F("h_stk_clusterY","h_stk_clusterY",10000,0,10000)
    h_psd_chargeX_size = TH1F("h_psd_chargeX_size","h_psd_chargeX_size",10,0,10)
    h_psd_chargeY_size = TH1F("h_psd_chargeY_size","h_psd_chargeY_size",10,0,10)
    h_psd_chargeX = TH1F("h_psd_chargeX","h_psd_chargeX",300,0,30)
    h_psd_chargeY = TH1F("h_psd_chargeY","h_psd_chargeY",300,0,30)
    h_psd_chargeX_corr = TH1F("h_psd_chargeX_corr","h_psd_chargeX_corr",300,0,30)
    h_psd_chargeY_corr = TH1F("h_psd_chargeY_corr","h_psd_chargeY_corr",300,0,30)
    h_psd_energy_hit = TH1F("h_psd_energy_hit","h_psd_energy_hit",1000,0,100)
    h_psd_chargeX_stk_selection    = TH1F("h_psd_chargeX_stk_selection","h_psd_chargeX_stk_selection",300,0,30)
    h_psd_chargeY_stk_selection    = TH1F("h_psd_chargeY_stk_selection","h_psd_chargeY_stk_selection",300,0,30)
    h_psd_chargeYvsE_stk_selection = TH2F("h_psd_chargeYvsE_stk_selection","h_psd_chargeYvsE_stk_selection",6000,0,6,300,0,30)
    h_bgo_theta = TH1F("h_bgo_theta","h_bgo_theta",110,-10,100)
    h_bgo_stk_deltaTheta = TH1F ("h_bgo_stk_deltaTheta","h_bgo_stk_deltaTheta",400,-100,100)
    h_bgo_stk_deltaTheta_selection = TH1F ("h_bgo_stk_deltaTheta_selection","h_bgo_stk_deltaTheta_selection",400,-100,100)
    h_bgo_stk_deltaTheta_tracksel = TH1F ("h_bgo_stk_deltaTheta_tracksel","h_bgo_stk_deltaTheta_tracksel",400,-100,100)
    h_psd_GlobalID = TH1F("h_psd_GlobalID","h_psd_GlobalID",7000,0,7000)
    h_psd_corr_ID512  = TH2F("h_psd_corr_ID512","h_psd_corr_ID512",1000,-500,500,200,0,2.)
    h_psd_corr_ID5504 = TH2F("h_psd_corr_ID5504","h_psd_corr_ID5504",1000,-500,500,200,0,2.)
    h_stk_to_psd = TH2F("h_stk_to_psd","h_stk_to_psd",1000,-1000,1000,1000,-1000,1000)
    h_bgo_firedLayer_before = TH1F("h_bgo_firedLayer_before","h_bgo_firedLayer_before",15,0,15)
    h_bgo_ene4lay_before = TH1F("h_bgo_ene4lay_before","h_bgo_ene4lay_before",120,0,1.2)
    h_bgo_ene4layLast_before = TH1F("h_bgo_ene4layLast_before","h_bgo_ene4layLast_before",120,0,1.2)
    h_bgo_ene4Lay_top_vs_bot_before = TH2F("h_bgo_ene4Lay_top_vs_bot_before","h_bgo_ene4Lay_top_vs_bot_before",120,0,1.2,120,0,1.2)
    h_bgo_ene4lay_after = TH1F("h_bgo_ene4lay_after","h_bgo_ene4lay_after",120,0,1.2)
    h_bgo_ene4layLast_after = TH1F("h_bgo_ene4layLast_after","h_bgo_ene4layLast_after",120,0,1.2)
    h_bgo_ene4Lay_top_vs_bot_after = TH2F("h_bgo_ene4Lay_top_vs_bot_after","h_bgo_ene4Lay_top_vs_bot_after",120,0,1.2,120,0,1.2)
    h_psd_chargeY_before = TH1F("h_psd_chargeY_before","h_psd_chargeY PSDX or PSDY",1200,0,30)
    h_psd_chargeX_before = TH1F("h_psd_chargeX_before","h_psd_chargeX PSDX or PSDY",1200,0,30)
    h_psd_chargeY_after  = TH1F("h_psd_chargeY_after","h_psd_chargeY PSDX and PSDY",1200,0,30)
    h_psd_chargeX_after  = TH1F("h_psd_chargeX_after","h_psd_chargeX PSDX and PSDY",1200,0,30)
    h_psd_resX = []
    h_psd_resY = []
    h_psd_resX_top = TH1F("h_psd_resX_top","h_psd_resX_top",400,-20,20) 
    h_psd_resY_top = TH1F("h_psd_resY_top","h_psd_resY_top",400,-20,20) 
    h_psd_resX_bot = TH1F("h_psd_resX_bot","h_psd_resX_bot",400,-20,20) 
    h_psd_resY_bot = TH1F("h_psd_resY_bot","h_psd_resY_bot",400,-20,20) 

    for iNpsd in xrange(0,2):
        h_psd_resX_title = "h_psd_resX_" + str(iNpsd) 
        h_psd_resY_title = "h_psd_resY_" + str(iNpsd) 
        h_psdX_histo = TH1F(h_psd_resX_title,h_psd_resX_title,2000,-1000,1000)
        h_psdY_histo = TH1F(h_psd_resY_title,h_psd_resY_title,2000,-1000,1000)
        h_psd_resX.append(h_psdX_histo)
        h_psd_resY.append(h_psdY_histo)


    has_quenching = False
    if opts.mc:
        for branch in dmpch.GetListOfBranches():
            if (branch.GetName() == "DmpEvtBgoQuenchRec"):
                has_quenching = True
                print "BGO Quenching: ON"
                break
                
        if has_quenching == False:
            print "BGO Quenching: OFF"
        
    if opts.mc:
        h_mc_stopZ = TH1F("h_mc_stopZ","h_mc_stopZ",2000,-500,500)
        h_mc_stopX_vs_Y = TH2F("h_mc_stopX_vs_Y","h_mc_stopX_vs_Y",2000,-500,500,2000,-500,500)
        ### check the good definition 
        h_stopZ_good = TH1F("h_stopZ_good","h_stopZ_good",2000,-500,500)
        h_stopX_good = TH1F("h_stopX_good","h_stopX_good",1000,-500,500)
        h_stopY_good = TH1F("h_stopY_good","h_stopY_good",1000,-500,500)
        h_bgo_projection_good_Top = TH2F("h_bgo_projection_good_Top","h_bgo_projection_good_Top",2000,-500,500,2000,-500,500)
        h_bgo_projection_good_Bot = TH2F("h_bgo_projection_good_Bot","h_bgo_projection_good_Bot",2000,-500,500,2000,-500,500)
        h_bgo_deltaTheta_bgo = TH1F ("h_bgo_deltaTheta_bgo","h_bgo_deltaTheta_bgo",200,-5,5)
        h_bgo_deltaTheta_stk = TH1F ("h_bgo_deltaTheta_stk","h_bgo_deltaTheta_stk",200,-5,5)


    if opts.mc:
        h_energy_good = []
        h_energy_bad = []
        h_energy_matrix = []
        h_energy_truth = {}
        h_energy_truth_weight = {}
        h_energy_weight = {}
        for ihisto_mc in range(0,ncuts): 
            h_title_good  = 'h_energy_good_cut_' + str(ihisto_mc)
            h_ene_cut_good = TH1F(h_title_good, h_title_good , nEDBINs_120bins , Edges_120bins)
            h_ene_cut_good.Sumw2()
            h_energy_good.append(h_ene_cut_good)
            h_title_bad = 'h_energy_bad_cut_' + str(ihisto_mc)
            h_ene_cut_bad = TH1F(h_title_bad, h_title_bad , nEDBINs_120bins , Edges_120bins)
            h_ene_cut_bad.Sumw2()
            h_energy_bad.append(h_ene_cut_bad)
            h_energy_matrix_title = "h_energy_matrix_cut_" + str(ihisto_mc)
            h_energy_matrix_histo = TH2F(h_energy_matrix_title,h_energy_matrix_title, nEDBINs_120bins, Edges_120bins, nEDBINs_120bins, Edges_120bins)
            h_energy_matrix.append(h_energy_matrix_histo)
            h_energy_truth_weight[ihisto_mc] = {}
            h_energy_weight[ihisto_mc] = {}
            h_energy_truth[ihisto_mc] = {}
            for bname, (nbins, edges) in energy_binnings.items():
                name_t = "h_energy_truth_{}_cut_{}".format(bname, ihisto_mc)
                h_t = TH1F(name_t, name_t, nbins, edges)
                h_t.Sumw2()
                h_energy_truth[ihisto_mc][bname]= h_t
                h_energy_truth_weight[ihisto_mc][bname] = {}
                h_energy_weight[ihisto_mc][bname] = {}
                for tag, gamma in spectral_indices.items():
                    name_tw = "h_energy_truth_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc)
                    h_tw = TH1F(name_tw, name_tw, nbins, edges)
                    h_tw.Sumw2()
                    h_energy_truth_weight[ihisto_mc][bname][tag] = h_tw
                    name = "h_energy_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc)
                    h = TH1F(name, name, nbins, edges)
                    h.Sumw2()
                    h_energy_weight[ihisto_mc][bname][tag] = h

    h_energy_NoTracks = TH1F("h_energy_NoTracks","h_energy_NoTracks",nEDBINs_120bins , Edges_120bins)
    h_energy_trigger_check = TH1F("h_energy_trigger_check","h_energy_trigger_check",nEDBINs_120bins , Edges_120bins) 
    h_energy = []
    
    for ihisto in range(0,ncuts):
        # print i 
        h_title = 'h_energy_cut_' + str(ihisto)
        h_ene_cut = TH1F(h_title, h_title , nEDBINs_120bins , Edges_120bins)
        h_energy.append(h_ene_cut)

    h_stk_cluster_XvsY = []
    for iplane in range(0,6):
        h_name = 'stk_cluster_XvsY_plane_'+str(iplane)
        h_title = 'cluster X vs Y - plane '+str(iplane)
        h_stk_cluster_plot = TH2F(h_name,h_title,1000,-500,500,1000,-500,500)
        h_stk_cluster_XvsY.append(h_stk_cluster_plot)


    h_bgo_ene_lay = []                                                                                                                               
    for ilay in range(0,14):                                                                                                                         
        h_title = 'h_ngo_ene_lay_'+str(ilay)                                                                                                         
        h_bgo_ene_layer = TH1F(h_title,h_title,1000,0,100)                                                                                           
        h_bgo_ene_lay.append(h_bgo_ene_layer)


    ###########################################################################

    nevents = dmpch.GetEntries()
    print ' Number of events: ',nevents
    
    ntrack_newsel = 0
 
    dmpch.SetOutputDir(root_dir,stream_name) #
    if opts.skim:
        print ' skim option activeted'


    # PeVSpace DAMPE Api Setup
    wrong_id = 0
    entry_ntupla = 0

    for iev in xrange(0,nevents):  ## in the skim the event 0 is not taken into account
        if(iev % 10000 == 0): print 'Processing event: ',iev
        if opts.mc:
            DmpVSvc.gPsdECor.SetMCflag(1)

        pev=dmpch.GetDmpEvent(iev)
        metach.GetEntry(iev)

        good_event = False

        stopZ_parent  = -1000
        stopX_parent  = -1000
        stopY_parent  = -1000
        startX_parent = -1000
        startY_parent = -1000
        startZ_parent = -1000
           
        bgorec = pev.pEvtBgoRec()
        bgoquenchrec = pev.pEvtBgoQuenchRec()

        stkclusters = pev.GetStkSiClusterCollection()
        stkladders = pev.GetStkLadderAdcCollection()

        # Initialization of pevspace api on first event
        if opts.mc and iev == 0:
            pevspace_api.Initialize(False, bgorec, stkclusters, stkladders, config="ions")
            #check_api.Initialize(False, bgorec, stkclusters, stkladders, config="electronfull2")
            bgo_api.Initialize(bgorec)
            bgo_api_quench.Initialize(bgoquenchrec)
        elif not opts.mc and iev == 0:
            #pevspace_api.Initialize(True, bgorec, stkclusters, runmetadata, config="ions")
            pevspace_api.Initialize(True, bgorec, stkclusters, None, runmetadata, config="ions")
            #check_api.Initialize(True, bgorec, stkclusters, None, runmetadata, config="electronfull2")
            bgo_api.Initialize(bgorec)
        # Rebinding of pevspace api in successive events
        elif opts.mc and iev > 0:
            pevspace_api.Rebind(bgorec, stkclusters, stkladders)
            #check_api.Rebind(bgorec, stkclusters, stkladders)
            bgo_api.BindBgoRec(bgorec)
            bgo_api_quench.BindBgoRec(bgoquenchrec)
        elif not opts.mc and iev > 0:
            pevspace_api.Rebind(bgorec, stkclusters, 0, runmetadata)
            #check_api.Rebind(bgorec, stkclusters, 0, runmetadata)
            bgo_api.BindBgoRec(bgorec)
        
        pevspace_api.Predict(bgodirection=True, stkvertex=True, stktrack=True)
        pevspace_api.ObtainTrackHits()

        #check_api.Predict(bgodirection=True, stkvertex=True, stktrack=True)
        #check_api.ObtainTrackHits()

        ## for MC nuclei
        if opts.mc:
            parent_particleID = int(pev.pEvtSimuPrimaries().pvpart_pdg/10000.)
#           print '----------- parent: ',parent_particleID
            if(int(parent_particleID) != mc_pdg_particle_id):
                print "*************WRONG PARTICLE ID, particle ID ", pev.pEvtSimuPrimaries().pvpart_pdg
                wrong_id += 1
                continue
            
            for it in xrange(0, pev.NSimuTrajectory()):
                if (pev.pSimuTrajectory(0).parentID == 0):
                    stopZ_parent  =  pev.pSimuTrajectory(0).stop_z
                    stopX_parent  =  pev.pSimuTrajectory(0).stop_x
                    stopY_parent  =  pev.pSimuTrajectory(0).stop_y
                    startX_parent = pev.pSimuTrajectory(0).start_x
                    startY_parent = pev.pSimuTrajectory(0).start_y
                    startZ_parent = pev.pSimuTrajectory(0).start_z
                else:
                    print " ERROR NO PARENT!!! "
                    break
            
            h_mc_stopX_vs_Y.Fill(stopX_parent,stopY_parent)
            h_mc_stopZ.Fill(stopZ_parent)
            def_good_truth_stopXY = False
            if(stopZ_parent  < 0.):
                def_good_truth_stopXY = abs(stopX_parent) < 400. and abs(stopY_parent )< 400.
            else:
                def_good_truth_stopXY = abs(stopX_parent) < 300. and abs(stopY_parent )< 300.
                                        
            def_good_truth_stopZ = (stopZ_parent <= 350. and stopZ_parent >= -280.)
            if startZ_parent != stopZ_parent:
                interceptXZ_truth = ((startZ_parent*stopX_parent) - (stopZ_parent*startX_parent) )/(startZ_parent-stopZ_parent)
                interceptYZ_truth = ((startZ_parent*stopY_parent) - (stopZ_parent*startY_parent) )/(startZ_parent-stopZ_parent)
                slope_X = (  startX_parent - stopX_parent)/(  startZ_parent - stopZ_parent)
                slope_Y = (  startY_parent - stopY_parent)/(  startZ_parent - stopZ_parent)
                x_projection_truth_onBGO_Top = interceptXZ_truth  + BGOzTop  * slope_X 
                y_projection_truth_onBGO_Top = interceptYZ_truth  + BGOzTop  * slope_Y 
                x_projection_truth_onBGO_Bot = interceptXZ_truth  + BGOzBot  * slope_X 
                y_projection_truth_onBGO_Bot = interceptYZ_truth  + BGOzBot  * slope_Y 
                def_good_truth_bgoAcc = abs(x_projection_truth_onBGO_Top) < bgo_acceptance_range and  abs(y_projection_truth_onBGO_Top) < bgo_acceptance_range and abs(x_projection_truth_onBGO_Bot) < bgo_acceptance_range and  abs(y_projection_truth_onBGO_Bot) < bgo_acceptance_range
            else: 
                def_good_truth_bgoAcc = False

            if(def_good_truth_stopXY==True and def_good_truth_stopZ == True and  def_good_truth_bgoAcc == True):
                good_event = True
                h_stopZ_good.Fill(stopZ_parent)
                h_stopX_good.Fill(stopX_parent)
                h_stopY_good.Fill(stopY_parent)
                h_bgo_projection_good_Top.Fill(x_projection_truth_onBGO_Top,y_projection_truth_onBGO_Top)
                h_bgo_projection_good_Bot.Fill(x_projection_truth_onBGO_Bot,y_projection_truth_onBGO_Bot)
                if(abs(stopX_parent) > 400. or abs(stopY_parent) > 400.):
                    print iev, stopX_parent, stopY_parent


        # Get energy corrected for saturation and quenching
        etot=pev.pEvtBgoRec().GetTotalEnergy()/1000. 
        etot_quench=-1
        etot_satcorrZmc=pev.pEvtBgoRec().GetDy2SatEcor(mc_particle_Z)/1000.
        fBGO_EnergyG_SatCorr_ML[mc_particle_Z]=etot
        fBGO_EnergyG_SatCorr_ML_ions[0]=etot
        fBGO_EnergyG_QuenchSatCorr_ML_ions[0]=etot
        etot_quench_satcorrZmc=-1 
        if has_quenching:
            try:
                etot_quench=pev.pEvtBgoQuenchRec().GetTotalEnergy()/1000.
                etot_quench_satcorrZmc=pev.pEvtBgoQuenchRec().GetDy2SatEcor(mc_particle_Z)/1000.
                fBGO_EnergyG_QuenchSatCorr_ML[mc_particle_Z]=etot
                fBGO_EnergyG_QuenchSatCorr_ML_ions[0]=etot_quench
            except:
                etot_quench=etot 
                etot_quench_satcorrZmc=etot_satcorrZmc
                fBGO_EnergyG_QuenchSatCorr_ML_ions[0]=etot 
        else:
            etot_quench=etot 
            etot_quench_satcorrZmc=etot_satcorrZmc
            fBGO_EnergyG_QuenchSatCorr_ML_ions[0]=etot 

        #saturation doesn't occur below 1TeV
        if pev.pEvtBgoRec().GetTotalEnergy() > 1e6:
            bgo_api.Predict(mc_particle_Z)
            if bgo_api.IsSaturated():
                fBGO_EnergyG_SatCorr_ML[mc_particle_Z] = bgo_api.GetReconstructedBGOE()/1000.

        if pev.pEvtBgoRec().GetTotalEnergy() > 1e6:
            bgo_api.Predict(mc_particle_Z, method = "ions") #Z is a dummy value
            if bgo_api.IsSaturated():
                fBGO_EnergyG_SatCorr_ML_ions[0]=bgo_api.GetReconstructedBGOE()/1000.
                fBGO_EnergyG_QuenchSatCorr_ML_ions[0]=bgo_api.GetReconstructedBGOE()/1000.

        ##################################
        etot_truth = -9.
        if opts.mc:
            etot_truth = pev.pEvtSimuPrimaries().pvpart_ekin/1000. 
        
        if not opts.mc:
            longitude = pev.pEvtAttitude().lon_geo 
            latitude  = pev.pEvtAttitude().lat_geo
            
        else:
            longitude= 0
            latitude = 0

        flong[0] = longitude
        flat[0] = latitude
            
        h_energy[0].Fill(etot)
        
        if opts.mc:
            h_energy_matrix[0].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[0][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[0][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[0][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[0].Fill(etot)
            else:
                h_energy_bad[0].Fill(etot)
                
        ####### Energy cut
        if(etot < e_min): continue
        h_energy[1].Fill(etot)
        if opts.mc:
            h_energy_matrix[1].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[1][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[1][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[1][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[1].Fill(etot)
            else:
                h_energy_bad[1].Fill(etot)

        ###### BGO acceptance at the 3rd last layer of the BGO
        x_projection_fromBGO_to_BGO =  pevspace_api.GetDirectionBGOInteceptX() + BGOzBot* pevspace_api.GetDirectionBGOSlopeX()
        y_projection_fromBGO_to_BGO =  pevspace_api.GetDirectionBGOInteceptY() + BGOzBot* pevspace_api.GetDirectionBGOSlopeY()

        x_projection_fromBGO_to_BGOTop =  pevspace_api.GetDirectionBGOInteceptX() + BGOzTop  * pevspace_api.GetDirectionBGOSlopeX()
        y_projection_fromBGO_to_BGOTop =  pevspace_api.GetDirectionBGOInteceptY() + BGOzTop  * pevspace_api.GetDirectionBGOSlopeY()
        
        h_bgo_projection_before.Fill(x_projection_fromBGO_to_BGO,y_projection_fromBGO_to_BGO)
        h_bgo_projectionTop_before.Fill(x_projection_fromBGO_to_BGOTop,y_projection_fromBGO_to_BGOTop)
        
        bgo_acceptance = abs(x_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(x_projection_fromBGO_to_BGOTop) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGOTop) < bgo_acceptance_range

        if( ( pevspace_api.GetDirectionBGOSlopeX() == 0 and pevspace_api.GetDirectionBGOInteceptX() == 0 ) or ( pevspace_api.GetDirectionBGOSlopeY() == 0  and pevspace_api.GetDirectionBGOInteceptY() == 0 ) ): 
            #print " Event with all 0 ", iev
            continue

        fBGO_interceptX_ML[0] = pevspace_api.GetDirectionBGOInteceptX()
        fBGO_interceptY_ML[0] = pevspace_api.GetDirectionBGOInteceptY() 
        fBGO_slopeX_ML[0] = pevspace_api.GetDirectionBGOSlopeX()
        fBGO_slopeY_ML[0] = pevspace_api.GetDirectionBGOSlopeY()

        if not (bgo_acceptance): continue
        h_energy[2].Fill(etot)
        if opts.mc:
            h_energy_matrix[2].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[2][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[2][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[2][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[2].Fill(etot)
            else:
                h_energy_bad[2].Fill(etot)

        h_bgo_projection_after.Fill(x_projection_fromBGO_to_BGO,y_projection_fromBGO_to_BGO)
        h_bgo_projectionTop_after.Fill(x_projection_fromBGO_to_BGOTop,y_projection_fromBGO_to_BGOTop)
        
        
        ######## MaxEneLay #########
        v_bgolayer  = np.array([pev.pEvtBgoRec().GetELayer(ibgo) for ibgo in range(14)])
        frac_max_lay = np.max(v_bgolayer)/1000./etot
        
        h_bgo_EratioLay_tot_before.Fill(frac_max_lay)
        if(good_event == True): 
            h_bgo_EratioLay_goo_before.Fill(frac_max_lay)
        else:
            h_bgo_EratioLay_bad_before.Fill(frac_max_lay)

        if(frac_max_lay > 0.35):   continue      
        h_bgo_EratioLay_tot_after.Fill(frac_max_lay)

        h_energy[3].Fill(etot)
        if opts.mc:
            h_energy_matrix[3].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[3][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[3][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[3][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[3].Fill(etot)
            else:
                h_energy_bad[3].Fill(etot)
                

        ######### BGO LATERAL CUT 
        pos_max_lateral = False
        lateral_pos_vec = []
        for ilay in xrange(1,4):
            v_bgolayer_bars  = np.array([pev.pEvtBgoRec().GetEdepPos(ilay,ibar) for ibar in xrange(0,23)])
            lateral_pos = np.argmax(v_bgolayer_bars)
            ene_lay = np.max(v_bgolayer_bars)

            if(ene_lay <= 0.):
                lateral_pos = -1
            h_bgo_lateralpos_before.Fill(lateral_pos)
            lateral_pos_vec.append(lateral_pos)
  

        ############# Test Xin method #######
        #  print 'nbgo hits ',pev.NEvtBgoHits()
        barNumberMaxEBarLayer = []
        MaxEBarLayer = []
        for ia in xrange(0,3):
            barNumberMaxEBarLayer.append(-1)
            MaxEBarLayer.append(0)
            
        for ihit in xrange(0,pev.NEvtBgoHits()):
            hitE = pev.pEvtBgoHits().fEnergy[ihit]
            layID = pev.pEvtBgoHits().GetLayerID(ihit)
            if(layID == 1 or layID == 2 or layID ==3):
                if(hitE > MaxEBarLayer[layID -1]):
                    barNumberMaxEBarLayer[layID -1] = (pev.pEvtBgoHits().fGlobalBarID[ihit] >> 6) & 0x1f
                    MaxEBarLayer[layID -1] = hitE
        
        for jj in xrange(0,len(barNumberMaxEBarLayer)):
            if(barNumberMaxEBarLayer[jj] <= 0 or barNumberMaxEBarLayer[jj] == 21):
                pos_max_lateral = True ### bug solved!

        if(pos_max_lateral == True):
            continue
    

        for ilay in xrange(0,3):
            h_bgo_lateralpos_after.Fill(lateral_pos_vec[ilay])
           
        h_energy[4].Fill(etot)
        if opts.mc:
            h_energy_matrix[4].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[4][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[4][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[4][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[4].Fill(etot)
            else:
                h_energy_bad[4].Fill(etot)

        if  opts.skim:
            dmpch.SaveCurrentEvent(stream_name)

        ################################## AS SKIMMED FILES UP TO HERE ############################
        # SAA only for data:
        if not opts.mc:
            inSAA = pFilter.IsInSAA(pev.pEvtHeader().GetSecond())
            #inSAA = False
            if (inSAA): continue 
            h_energy[5].Fill(etot)
            h_terrestrial_lat_vs_long.Fill(longitude,latitude)
        else:
            h_energy[5].Fill(etot)
            h_energy_matrix[5].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[5][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[5][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[5][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[5].Fill(etot)
            else:
                h_energy_bad[5].Fill(etot)
                
        
        if(pev.pEvtHeader().GeneratedTrigger(3) == True or pev.pEvtHeader().GeneratedTrigger(0) == True ):
            h_energy_trigger_check.Fill(etot)
        
        h_energy[6].Fill(etot)
        
                
                
        # TRACKING
        fSTK_vertexPrediction[0] = pevspace_api.GetVertexPrediction()

        theta_bgo = np.arctan(np.sqrt(pevspace_api.GetDirectionBGOSlopeX()**2 + pevspace_api.GetDirectionBGOSlopeY()**2))
        theta_stk = np.arctan(np.sqrt(pevspace_api.GetDirectionSTKSlopeX()**2 + pevspace_api.GetDirectionSTKSlopeY()**2))
        track_correction = np.cos(theta_stk)
        theta_bgo_deg = theta_bgo*180/np.pi 
        theta_stk_deg = theta_stk*180/np.pi 

        phi_bgo = math.atan(pevspace_api.GetDirectionBGOSlopeY() / pevspace_api.GetDirectionBGOSlopeX())
        phi_stk = math.atan(pevspace_api.GetDirectionSTKSlopeY() / pevspace_api.GetDirectionSTKSlopeX())
        phi_bgo_deg = phi_bgo*180/np.pi 
        phi_stk_deg = phi_stk*180/np.pi

        fBGO_slopeX_ML[0] = pevspace_api.GetDirectionBGOSlopeX()
        fBGO_slopeY_ML[0] = pevspace_api.GetDirectionBGOSlopeY()
        fBGO_slopeXZ_analy[0] = pev.pEvtBgoRec().GetSlopeXZ()
        fBGO_slopeYZ_analy[0] = pev.pEvtBgoRec().GetSlopeYZ()
        fBGO_theta_analy[0] = (np.arctan(np.sqrt(pev.pEvtBgoRec().GetSlopeXZ()**2 + pev.pEvtBgoRec().GetSlopeYZ()**2)))*180/np.pi
        if pev.pEvtBgoRec().GetSlopeXZ() > 0:
            fBGO_phi_analy[0] = (math.atan(pev.pEvtBgoRec().GetSlopeYZ() / pev.pEvtBgoRec().GetSlopeXZ()))*180/np.pi
        fBGO_interceptXZ_analy[0] = pev.pEvtBgoRec().GetInterceptXZ()
        fBGO_interceptYZ_analy[0] = pev.pEvtBgoRec().GetInterceptYZ()

        fSTK_slopeX_ML[0] = pevspace_api.GetDirectionSTKSlopeX()
        fSTK_slopeY_ML[0] = pevspace_api.GetDirectionSTKSlopeY()
        fSTK_interceptX_ML[0] = pevspace_api.GetTrackHitImpactX()
        fSTK_interceptY_ML[0] = pevspace_api.GetTrackHitImpactY()
        fSTK_interceptZ_ML[0] = pevspace_api.GetTrackHitImpactZ()

        hitClusterSum = [0,0]

        etaCorr = DmpStkClusterEtaCorrection()
        cluster_chargeX =  array('d', 6*[ -999.] )
        cluster_chargeY =  array('d', 6*[ -999.] )
        cluster_chargeX_etacorr =  array('d', 6*[ -999.] )
        cluster_chargeY_etacorr =  array('d', 6*[ -999.] )


        for j in range(12):
            hitclustersignal = pevspace_api.GetTrackHitSignal(j)
            hitclusterid = pevspace_api.GetTrackHitID(j)     # ID of hit cluster in the [stkclusters] collection
            if hitclusterid >= 0:             #  - cluster/hit for this layer
                hitClusterSum[j%2] += 1
                cluster = stkclusters.ConstructedAt(hitclusterid)
                if j%2 == 0:
                    cluster_chargeX[j/2] = hitclustersignal
                    cluster_chargeX_etacorr[j/2] = etaCorr.getStkEtaCorrEnergy(cluster, track_correction)*track_correction
                elif j%2 == 1:
                    cluster_chargeY[j/2] = hitclustersignal
                    cluster_chargeY_etacorr[j/2] = etaCorr.getStkEtaCorrEnergy(cluster, track_correction)*track_correction

        
        # No good tracks found
        if ( hitClusterSum[0] == 0 or hitClusterSum[1] == 0): continue
                                
        h_energy[7].Fill(etot)
        if opts.mc:
            h_energy_matrix[7].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[7][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[7][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[7][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[7].Fill(etot)
            else:
                h_energy_bad[7].Fill(etot)


        #if abs(theta_bgo - theta_stk) > delta_theta_max: continue
        
        h_energy[8].Fill(etot)

        #    print 'theta bgo', theta_bgo
        deltaTheta_rec_sel = theta_bgo_deg - theta_stk_deg
        h_bgo_stk_deltaTheta_tracksel.Fill(deltaTheta_rec_sel)

        theta_truth  = -9.
        phi_truth = -9.
        if opts.mc:
            theta_truth  = math.atan(np.sqrt( (pev.pEvtSimuPrimaries().pvpart_px/pev.pEvtSimuPrimaries().pvpart_pz*pev.pEvtSimuPrimaries().pvpart_px/pev.pEvtSimuPrimaries().pvpart_pz)+(pev.pEvtSimuPrimaries().pvpart_py/pev.pEvtSimuPrimaries().pvpart_pz*pev.pEvtSimuPrimaries().pvpart_py/pev.pEvtSimuPrimaries().pvpart_pz)))*180./math.pi    
            phi_truth = math.atan(pev.pEvtSimuPrimaries().pvpart_py / pev.pEvtSimuPrimaries().pvpart_px)*180./math.pi
            #print "truth ", theta_truth, " stk ", theta_track, "bgo ", theta_bgo
            delta_theta_bgo = theta_truth - theta_bgo
            delta_theta_stk = theta_truth - theta_stk
            h_bgo_deltaTheta_bgo.Fill(delta_theta_bgo)
            h_bgo_deltaTheta_stk.Fill(delta_theta_stk)


        #### projection on PSD for charge selection
        # cluster Charge X and Y
        h_stk_clusterX.Fill(cluster_chargeX[0])
        h_stk_clusterY.Fill(cluster_chargeY[0])        


        ########## PSD FIDUCIAL VOLUME CUT
        track_slope = TVector3(pevspace_api.GetDirectionSTKSlopeX(), pevspace_api.GetDirectionSTKSlopeY(), 1)
        track_impact_point = TVector3(pevspace_api.GetTrackHitImpactX(), pevspace_api.GetTrackHitImpactY(), pevspace_api.GetTrackHitImpactZ())

        psd_YZ_top = -324.7
        psd_XZ_top = -298.5
        stk_to_psd_topY = (track_slope.Y()*(psd_YZ_top - track_impact_point.Z()) + track_impact_point.Y())
        stk_to_psd_topX = (track_slope.X()*(psd_XZ_top - track_impact_point.Z()) + track_impact_point.X())

        #if(abs(stk_to_psd_topX) > 400.): continue
        #if(abs(stk_to_psd_topY) > 400.): continue
        fSTKtrack_to_PSD_topY[0] = stk_to_psd_topY # ELISABETTA: STK track projection on top of PSD
        fSTKtrack_to_PSD_topX[0] = stk_to_psd_topX # ELISABETTA: STK track projection on top of PSD
        
        
        h_energy[9].Fill(etot)


        ######### STK - PSD MATCH #############

        # psd Charge X and Y        
        psd_vec_chargeX     = [[]for _ in range(2)]
        psd_vec_chargeX0    = [[]for _ in range(2)]
        psd_vec_chargeX1    = [[]for _ in range(2)]
        psd_vec_gidX        = [[]for _ in range(2)]
        psd_vec_pathlengthX = [[]for _ in range(2)]
        psd_vec_positionX   = [[]for _ in range(2)]
        psd_bar_id_x        = [[]for _ in range(2)]

        psd_vec_chargeY     = [[]for _ in range(2)]
        psd_vec_chargeY0    = [[]for _ in range(2)]
        psd_vec_chargeY1    = [[]for _ in range(2)]
        psd_vec_gidY        = [[]for _ in range(2)]
        psd_vec_pathlengthY = [[]for _ in range(2)]
        psd_vec_positionY   = [[]for _ in range(2)]
        psd_bar_id_y        = [[]for _ in range(2)]

        for ipsd in xrange(0,pev.NEvtPsdHits()):
            if(pev.pEvtPsdHits().IsHitMeasuringX(ipsd)):
                crossingX = False
                lengthX = [-99999.,-99999.]
                array_lengthX = array('d',lengthX)
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -298.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -298.51): #-298.5 ):
                    npsdX = 0
    
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -284.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -284.51):
                    npsdX = 1
        
                if opts.mc:
                    crossingX = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_slope, track_impact_point, array_lengthX)
                else:
                    crossingX = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_slope,track_impact_point, array_lengthX)

                if(crossingX == True):
                    psd_vec_chargeX[npsdX].append(pev.pEvtPsdHits().fEnergy[ipsd]) 
                    psd_vec_chargeX0[npsdX].append(pev.pEvtPsdHits().fES0[ipsd]) 
                    psd_vec_chargeX1[npsdX].append(pev.pEvtPsdHits().fES1[ipsd])

                    psd_vec_gidX[npsdX].append(pev.pEvtPsdHits().fGlobalBarID[ipsd]) 
                    psd_bar_id_x[npsdX].append(pev.pEvtPsdHits().GetHittedBarNumber()) 
                    psd_vec_pathlengthX[npsdX].append(array_lengthX[1]) 
                    psd_vec_positionX[npsdX].append(pev.pEvtPsdHits().GetHitX(ipsd))
                    

            elif(pev.pEvtPsdHits().IsHitMeasuringY(ipsd)):
                crossingY = False
                lengthY = [-99999.,-99999.]
                array_lengthY = array('d',lengthY)
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -324.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -324.71):
                    npsdY = 0 
        
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -310.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -310.71):
                    npsdY = 1
        
                    
                if opts.mc:
                    crossingY = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd], track_slope, track_impact_point, array_lengthY)
                else:
                    crossingY = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_slope, track_impact_point, array_lengthY)

                if (crossingY == True):
                    psd_vec_chargeY[npsdY].append(pev.pEvtPsdHits().fEnergy[ipsd]) 
                    psd_vec_chargeY0[npsdY].append(pev.pEvtPsdHits().fES0[ipsd])
                    psd_vec_chargeY1[npsdY].append(pev.pEvtPsdHits().fES1[ipsd])
                    psd_vec_gidY[npsdY].append(pev.pEvtPsdHits().fGlobalBarID[ipsd]) 
                    psd_bar_id_y[npsdY].append(pev.pEvtPsdHits().GetHittedBarNumber()) 
                    psd_vec_pathlengthY[npsdY].append(array_lengthY[1]) 
                    psd_vec_positionY[npsdY].append(pev.pEvtPsdHits().GetHitY(ipsd))

        check_psd = False
        for ipsd in xrange(0,2):
            if(len(psd_vec_chargeY[ipsd])>1 or len(psd_vec_chargeX[ipsd])>1) :
                check_psd = True
    
        if(check_psd == True):
            print '******* ERROR!!!!!!!!!!!!! *****'
            print 'Event with more than 1 psd bar per layer', iev
            print psd_vec_chargeX, psd_vec_chargeY
            print psd_vec_pathlengthX, psd_vec_pathlengthY
            print psd_vec_positionX, psd_vec_positionY
            #print "track selected", track_ID, track_sel.getImpactPoint().x(), track_sel.getImpactPoint().y(), track_sel.getImpactPoint().z()


        psdchargeX = [-999.,-999.]
        psdchargeY = [-999.,-999.]
        psdchargeX_corr = [-999.,-999.]
        psdchargeX0_corr = [-999.,-999.]
        psdchargeX1_corr = [-999.,-999.]
        psdchargeX01_corr = [-999.,-999.]
        psdchargeY_corr = [-999.,-999.]
        psdchargeY0_corr = [-999., -999.]
        psdchargeY1_corr = [-999.,-999.]
        psdchargeY01_corr = [-999.,-999.]
        psdvec = []
        psdvec_pathlength = []
        psdchargeX_proj = [-999.,-999.]
        psdchargeY_proj = [-999.,-999.]
        psdX_pathlength  = [-999.,-999.]
        psdY_pathlength  = [-999.,-999.]
        psdX_position = [-999.,-999.]
        psdY_position = [-999.,-999.]
        PsdEC_tmpX = 0.
        PsdEC_tmpY = 0.
        
        for ipsd in xrange(0,2):
             
            if(len(psd_vec_chargeY[ipsd]) > 0):

                pos_max_len = np.argmax(psd_vec_pathlengthY[ipsd])
                #print "********* pos min layer", pos_max_len
                lengthY = [-99999.,-99999.]
                array_lengthY = array('d',lengthY)
                test_pos = False 
                if opts.mc:
                    
                    test_pos = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidY[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthY)
                else:
                    test_pos = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidY[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthY)
                    
                PsdEC_tmpY = -1.
                if(test_pos == True):
                    PsdEC_tmpY = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidY[ipsd][pos_max_len], array_lengthY[0]/10.)
                else:
                    print "****** ERROR! it should be true!!!"

                AttCorY0 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidY[ipsd][pos_max_len], 0, array_lengthY[0])    
                AttCorY1 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidY[ipsd][pos_max_len], 1, array_lengthY[0])

                lay = DmpPsdBase.GetLayerID(psd_vec_gidY[ipsd][pos_max_len])
                bar = DmpPsdBase.GetStripID(psd_vec_gidY[ipsd][pos_max_len])                                                              

                psdchargeY[ipsd] = psd_vec_chargeY[ipsd][pos_max_len]
                psdchargeY_corr[ipsd] = psd_vec_chargeY[ipsd][pos_max_len]*PsdEC_tmpY
                if(len(psd_vec_chargeY0[ipsd]) > 0):
                    #Charge After Attenuation Correction
                    psdchargeY0_corr[ipsd] = np.sqrt(psd_vec_chargeY0[ipsd][pos_max_len]/AttCorY0*10./array_lengthY[1])
                    #Quenching & Equalization Correction (only for FD)
                    if opts.data:
                            psdchargeY0_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 0, psdchargeY0_corr[ipsd])
                else:
                    psdchargeY0_corr[ipsd] = -1.

                if(len(psd_vec_chargeY1[ipsd]) > 0):
                    psdchargeY1_corr[ipsd] = np.sqrt(psd_vec_chargeY1[ipsd][pos_max_len]/AttCorY1*10./array_lengthY[1])
                    if opts.data: 
                        psdchargeY1_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 1, psdchargeY1_corr[ipsd])
                else:
                    psdchargeY1_corr[ipsd] = -1.
                psdchargeY01_corr[ipsd] = (psdchargeY0_corr[ipsd] + psdchargeY1_corr[ipsd])/2.0
                #Position Correction (only for FD)
                if opts.data:
                    psdchargeY01_corr[ipsd] = DmpVSvc.gPsdECor.GetPosCorrectedCharge(lay, bar, array_lengthY[0], psdchargeY01_corr[ipsd])
                psdchargeY_proj[ipsd] = array_lengthY[0]
                psdY_pathlength[ipsd] = array_lengthY[1]
                psdY_position[ipsd] =  psd_vec_positionY[ipsd][pos_max_len]


                h_psd_GlobalID.Fill(psd_vec_gidY[ipsd][pos_max_len])
                if(psdchargeY[ipsd] > 0.):
                    h_psd_chargeY.Fill(np.sqrt(psdchargeY[ipsd]/2.))
                if(psdchargeY_corr[ipsd] > 0):
                    h_psd_chargeY_corr.Fill(np.sqrt(psdchargeY_corr[ipsd]/2.))
                if(psd_vec_gidY[ipsd][pos_max_len] == 512 and array_lengthY[0]!= -99999.):
                    h_psd_corr_ID512.Fill(array_lengthY[0],PsdEC_tmpY)

            if(len(psd_vec_chargeX[ipsd]) > 0 ):
                pos_max_len = np.argmax(psd_vec_pathlengthX[ipsd])
                lengthX = [-99999.,-99999.]
                array_lengthX = array('d',lengthX)
                test_posX = False
                if opts.mc:
                    test_posX = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidX[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthX)
                else:
                    test_posX = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidX[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthX)
                PsdEC_tmpX = -1.
                if(test_posX == True):
                    PsdEC_tmpX = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidX[ipsd][pos_max_len],array_lengthX[0]/10.)
                
                AttCorX0 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidX[ipsd][pos_max_len], 0, array_lengthX[0])    
                AttCorX1 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidX[ipsd][pos_max_len], 1, array_lengthX[0]) 
                lay = DmpPsdBase.GetLayerID(psd_vec_gidX[ipsd][pos_max_len])
                bar = DmpPsdBase.GetStripID(psd_vec_gidX[ipsd][pos_max_len])                                                              


                psdchargeX[ipsd] = psd_vec_chargeX[ipsd][pos_max_len]
                psdchargeX_corr[ipsd] = psd_vec_chargeX[ipsd][pos_max_len]*PsdEC_tmpX
                if(len(psd_vec_chargeX0[ipsd]) > 0):
                    #Charge After Attenuation Correction
                    psdchargeX0_corr[ipsd] = np.sqrt(psd_vec_chargeX0[ipsd][pos_max_len]/AttCorX0*10./array_lengthX[1])
                    #Quenching & Equalization Correction (only for FD)
                    if opts.data: 
                        psdchargeX0_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 0, psdchargeX0_corr[ipsd])
                else:
                    psdchargeX0_corr[ipsd] = -1.

                if(len(psd_vec_chargeX1[ipsd]) > 0):
                    psdchargeX1_corr[ipsd] = np.sqrt(psd_vec_chargeX1[ipsd][pos_max_len]/AttCorX1*10./array_lengthX[1])
                    if opts.data: 
                        psdchargeX1_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 1, psdchargeX1_corr[ipsd])
                else:
                    psdchargeX1_corr[ipsd] = -1.
                psdchargeX01_corr[ipsd] = (psdchargeX0_corr[ipsd] + psdchargeX1_corr[ipsd])/2.0
                #Position Correction (only for FD)
                if opts.data: 
                    psdchargeX01_corr[ipsd] = DmpVSvc.gPsdECor.GetPosCorrectedCharge(lay, bar, array_lengthX[0], psdchargeX01_corr[ipsd])
                psdchargeX_proj[ipsd] = array_lengthX[0] 
                psdX_pathlength[ipsd] = array_lengthX[1]
                psdX_position[ipsd] =  psd_vec_positionX[ipsd][pos_max_len]


                h_psd_GlobalID.Fill(psd_vec_gidX[ipsd][pos_max_len])
                if(psdchargeX[ipsd] > 0.):
                    h_psd_chargeX.Fill(np.sqrt(psdchargeX[ipsd]/2.))          
                if(psdchargeX_corr[ipsd] > 0.):
                    h_psd_chargeX_corr.Fill(np.sqrt(psdchargeX_corr[ipsd]/2.))    
                if(psd_vec_gidX[ipsd][pos_max_len] == 5504 and array_lengthX[0] != -99999.):
                    h_psd_corr_ID5504.Fill(array_lengthX[0],PsdEC_tmpX)


        for ii in xrange(2):
            if(psdchargeY01_corr[ii] > 0.): 
                psdvec.append(psdchargeY01_corr[ii])
                psdvec_pathlength.append(psdY_pathlength[ii])
        for ii in xrange(2):
            if(psdchargeX01_corr[ii] > 0.): 
                psdvec.append(psdchargeX01_corr[ii])
                psdvec_pathlength.append(psdX_pathlength[ii])

        h_bgo_theta.Fill(theta_bgo)
        

        ########## contine from here with the ntupla construction        
        sum_len_psd_vec_chargeX = len(psd_vec_chargeX[0]) + len(psd_vec_chargeX[1])
        sum_len_psd_vec_chargeY = len(psd_vec_chargeY[0]) + len(psd_vec_chargeY[1])

        for ipsd in xrange(0,2):
            if(psdchargeY_corr[ipsd] > 0.):
                h_psd_chargeY_before.Fill(np.sqrt(psdchargeY_corr[ipsd]/2.))
            if(psdchargeX_corr[ipsd] > 0.):
                h_psd_chargeX_before.Fill(np.sqrt(psdchargeX_corr[ipsd]/2.))
        
        ### new cuts 
        if(sum_len_psd_vec_chargeY == 0): continue
        if(sum_len_psd_vec_chargeX == 0): continue

        for ipsd in xrange(0,2):
            if(psdchargeY_corr[ipsd] > 0.):
                h_psd_chargeY_after.Fill(np.sqrt(psdchargeY_corr[ipsd]/2.))
            if(psdchargeX_corr[ipsd] > 0.):
                h_psd_chargeX_after.Fill(np.sqrt(psdchargeX_corr[ipsd]/2.))


        h_energy[10].Fill(etot)
        if opts.mc:
            h_energy_matrix[10].Fill(etot_truth,etot)
            for bname, (nbins, edges) in energy_binnings.items():
                h_energy_truth[10][bname].Fill(etot_truth)
                for tag, gamma in spectral_indices.items():
                    weight = spectral_weight(etot_truth, gamma, edges)
                    h_energy_weight[10][bname][tag].Fill(etot, weight)
                    h_energy_truth_weight[10][bname][tag].Fill(etot_truth, weight)
            if(good_event == True):
                h_energy_good[10].Fill(etot)
            else:
                h_energy_bad[10].Fill(etot)


        ###### HITS ON ALL THE BGO LAYERS 
        ###### nlayers > 14
        fired_layer = np.count_nonzero(v_bgolayer)
        h_bgo_firedLayer_before.Fill(fired_layer)

        h_energy[11].Fill(etot)

        Ene4Layers = sum(v_bgolayer[:4])/1000/etot
        Ene4LayersLast = sum(v_bgolayer[-4:])/1000/etot

        Ene2Layers = sum(v_bgolayer[:2])/1000/etot
        Ene2LayersLast = sum(v_bgolayer[-2:])/1000/etot

        ####### Fill of the NTUPLE 
        vec_psdX_max = []
        vec_psdY_max = []
        for ipsd in xrange(0,pev.NEvtPsdHits()):
            if(pev.pEvtPsdHits().GetLayerID(ipsd) == 0):
                vec_psdY_max.append(pev.pEvtPsdHits().fEnergy[ipsd])
            elif (pev.pEvtPsdHits().GetLayerID(ipsd) == 1):
                vec_psdX_max.append(pev.pEvtPsdHits().fEnergy[ipsd])


        #### PSD 
        zpsd = [-324.7,-310.7,-298.5,-284.5]                                                    #M e A 17-10-17
                                        
        fpsdchargeX0[0] = psdchargeX01_corr[0]
        fpsdchargeX1[0] = psdchargeX01_corr[1]
        fpsdchargeY0[0] = psdchargeY01_corr[0]
        fpsdchargeY1[0] = psdchargeY01_corr[1]

        globalcharge = -999.
        globalaverage = -999.
        #print 'Event: ', iev, psdvec
        if(len(psdvec)==0):
            print 'There is no charge value!'
        elif(len(psdvec)==1):
            globalcharge = psdvec[0]
            globalaverage = psdvec[0]
            #print 'Only one value used for charge calculation: ', globalcharge
        elif(len(psdvec)>1 and len(psdvec)<5):
            PSDdiff = True
            for iQ in xrange(1,len(psdvec)):
                PSDdiff = PSDdiff and (abs(psdvec[iQ-1] - psdvec[iQ]) < 1) # consider only data before fragmentation (diff = 1)
                if((iQ == 1) and (not PSDdiff)):
                    psd_sum = psdvec[iQ-1]
                    psd_sum_weight = psdvec[iQ-1]*psdvec_pathlength[iQ-1]
                    psd_nlay = 1.
                    psd_nlayweight = psdvec_pathlength[iQ-1]
                elif((iQ == 1) and (PSDdiff)):
                    psd_sum = psdvec[iQ-1] + psdvec[iQ]
                    psd_sum_weight = (psdvec[iQ-1]*psdvec_pathlength[iQ-1]) + (psdvec[iQ]*psdvec_pathlength[iQ]) 
                    psd_nlay = 2.
                    psd_nlayweight = (psdvec_pathlength[iQ-1] + psdvec_pathlength[iQ])
                elif((iQ !=1) and (PSDdiff)):
                    psd_sum = psd_sum + psdvec[iQ]
                    psd_sum_weight = psd_sum_weight + (psdvec[iQ] * psdvec_pathlength[iQ])
                    psd_nlay = psd_nlay + 1.
                    psd_nlayweight = psd_nlayweight + psdvec_pathlength[iQ] 
                else:
                    break
            globalcharge = psd_sum/psd_nlay
            globalaverage = psd_sum_weight/psd_nlayweight
            #print 'The global charge value is ', psd_sum, ' divided by ', psd_nlay, ' that gives: ', globalcharge
        else:
            print 'Houston, we have a problem!'
        fPSD_Global_Charge[0] = globalcharge
        fPSD_Global_Charge_PathAverage[0] = globalaverage

        average_charge = ( (((TMath.Sign(1.,psdchargeY01_corr[0])+1.)/2.)*psdchargeY01_corr[0]*psdY_pathlength[0]) + (((TMath.Sign(1.,psdchargeY01_corr[1])+1.)/2.)*psdchargeY01_corr[1]*psdY_pathlength[1]) + (((TMath.Sign(1.,psdchargeX01_corr[0])+1.)/2.)*psdchargeX01_corr[0]*psdX_pathlength[0]) + (((TMath.Sign(1.,psdchargeX01_corr[1])+1.)/2.)*psdchargeX01_corr[1]*psdX_pathlength[1]) ) / ( (((TMath.Sign(1.,psdchargeY01_corr[0])+1.)/2.)*psdY_pathlength[0]) + (((TMath.Sign(1.,psdchargeY01_corr[1])+1.)/2.)*psdY_pathlength[1]) + (((TMath.Sign(1.,psdchargeX01_corr[0])+1.)/2.)*psdX_pathlength[0]) + (((TMath.Sign(1.,psdchargeX01_corr[1])+1.)/2.)*psdX_pathlength[1]) )
        fPSD_PathWeighted_Charge[0] = average_charge


        for ipsd in xrange(0,nlayer_psd):
            fPSD_EnergyY[ipsd] = psdchargeY[ipsd]
            fPSD_EnergyX[ipsd] = psdchargeX[ipsd]
        
            if (ipsd < 2):                                                          #M e A 17-10-17
                fPSD_cpsdstk[ipsd]=track_impact_point.y()+track_slope.y()*(zpsd[ipsd]-track_impact_point.z())      #M e A 17-10-17
            if (ipsd > 1):                                                                                                              #M e A 17-10-17
                fPSD_cpsdstk[ipsd]=track_impact_point.x()+track_slope.x()*(zpsd[ipsd]-track_impact_point.z())      #M e A 17-10-17   
 
           
            fPSD_EnergyY_corr[ipsd] = psdchargeY_corr[ipsd]
            fPSD_EnergyX_corr[ipsd] = psdchargeX_corr[ipsd]
            fPSD_psdX_projY[ipsd] = psdchargeX_proj[ipsd]
            fPSD_psdY_projX[ipsd] = psdchargeY_proj[ipsd]
            fPSD_psdX_pathlength[ipsd] = psdX_pathlength[ipsd]
            fPSD_psdY_pathlength[ipsd] = psdY_pathlength[ipsd] 
            fPSD_psdX_position[ipsd] = psdX_position[ipsd]
            fPSD_psdY_position[ipsd] = psdY_position[ipsd]  


        fPSD_nhits[0] = pev.NEvtPsdHits()
        fPSD_psdX_Emax_layer[0] = np.max(vec_psdX_max)
        fPSD_psdY_Emax_layer[0] = np.max(vec_psdY_max)   
    
        ######################### MARGHERITA&ANTONIO - PSD Emax cpsdmax (Oct 23, 2017)
        fPSD_EmaxY[0] = 0.0
        fPSD_EmaxY[1] = 0.0
        fPSD_EmaxX[0] = 0.0
        fPSD_EmaxX[1] = 0.0

        for ipsd in xrange(0,pev.NEvtPsdHits()):
            if(pev.pEvtPsdHits().GetHitZ(ipsd)== -324.7):
                if(fPSD_EmaxY[0] < pev.pEvtPsdHits().fEnergy[ipsd]):
                    fPSD_EmaxY[0] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[0] = pev.pEvtPsdHits().GetHitY(ipsd)
            if(pev.pEvtPsdHits().GetHitZ(ipsd)== -310.7):
                if(fPSD_EmaxY[1] < pev.pEvtPsdHits().fEnergy[ipsd]):
                    fPSD_EmaxY[1] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[1] = pev.pEvtPsdHits().GetHitY(ipsd)
            if(pev.pEvtPsdHits().GetHitZ(ipsd)== -298.5):
                if(fPSD_EmaxX[0] < pev.pEvtPsdHits().fEnergy[ipsd]):
                    fPSD_EmaxX[0] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[2] = pev.pEvtPsdHits().GetHitX(ipsd)+0.4
            if(pev.pEvtPsdHits().GetHitZ(ipsd)== -284.5):               
                if(fPSD_EmaxX[1] < pev.pEvtPsdHits().fEnergy[ipsd]):
                    fPSD_EmaxX[1] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[3] = pev.pEvtPsdHits().GetHitX(ipsd)+0.6
                #print " coor psd x1 ", pev.pEvtPsdHits().GetHitX(ipsd)
                #print " coor psd z ", pev.pEvtPsdHits().GetHitZ(ipsd)      

        ### MARGHERITA&ANTONIO - PSD Etracks (Oct 13, 2017)
        fPSD_Etrack01[0] = (((TMath.Sign(1,fPSD_EnergyY_corr[0])+1)/2)*fPSD_EnergyY_corr[0]+((TMath.Sign(1,fPSD_EnergyY_corr[1])+1)/2)*fPSD_EnergyY_corr[1])*10./(((TMath.Sign(1,fPSD_psdY_pathlength[0])+1)/2)*fPSD_psdY_pathlength[0]+((TMath.Sign(1,fPSD_psdY_pathlength[1])+1)/2)*fPSD_psdY_pathlength[1])
        fPSD_Etrack23[0] = (((TMath.Sign(1,fPSD_EnergyX_corr[0])+1)/2)*fPSD_EnergyX_corr[0]+((TMath.Sign(1,fPSD_EnergyX_corr[1])+1)/2)*fPSD_EnergyX_corr[1])*10./(((TMath.Sign(1,fPSD_psdX_pathlength[0])+1)/2)*fPSD_psdX_pathlength[0]+((TMath.Sign(1,fPSD_psdX_pathlength[1])+1)/2)*fPSD_psdX_pathlength[1])

  
        
        #### STK
        for iplane in xrange(0,6):
            fSTK_chargeX[iplane] = cluster_chargeX[iplane]
            fSTK_chargeY[iplane] = cluster_chargeY[iplane]
            fSTK_chargeX_etaCorr[iplane] = cluster_chargeX_etacorr[iplane]
            fSTK_chargeY_etaCorr[iplane] = cluster_chargeY_etacorr[iplane]
        
        fSTK_ntracks[0] = pev.NStkKalmanTrack()
        fSTK_nclusters[0] =  pev.NStkSiCluster()
        fSTK_theta_correction[0] = track_correction
        fSTK_trackIP[0] = track_impact_point.x()
        fSTK_trackIP[1] = track_impact_point.y()
        fSTK_trackIP[2] = track_impact_point.z() #track_sel.getImpactPoint().x()   MARGHERITA&ANTONIO (we changed x->z)
        fSTK_trackDirection[0] = track_slope.x()
        fSTK_trackDirection[1] = track_slope.y()

                      
        for ibgo in xrange(0,14):
            fBGO_EneLay[ibgo] = pev.pEvtBgoRec().GetELayer(ibgo)/1000.
            fBGO_nhits_layers[ibgo] = pev.pEvtBgoRec().GetLayerHits()[ibgo]

        EBarLayers   = [[]for _ in range(14)]
        PosBarLayers = [[]for _ in range(14)]
        PosXBarLayers = [[]for _ in range(14)]
        PosYBarLayers = [[]for _ in range(14)]
        rmsLayers = [[]for _ in range(14)]
        eCoreMaxLayer = [[]for _ in range(14)]
        eCoreCoordMaxLayer = [[]for _ in range(14)]        
        eCoreLayer = [[]for _ in range(14)]
        eCoreCoordLayer = [[]for _ in range(14)]        
        eCoreCoorAll = [[]for _ in range(14)]        
        rmsLayer = [[0.]for _ in range(14)]        
        sumRMS = 0.

        #      print " Event number ", iev
        for ihit in xrange(0,pev.NEvtBgoHits()):
            layID = pev.pEvtBgoHits().GetLayerID(ihit)
            hitE = pev.pEvtBgoHits().fEnergy[ihit]
            barNumber = (pev.pEvtBgoHits().fGlobalBarID[ihit] >> 6) & 0x1f
            EBarLayers[layID].append(hitE)
            PosBarLayers[layID].append(barNumber)
            PosXBarLayers[layID].append(pev.pEvtBgoHits().GetHitX(ihit))
            PosYBarLayers[layID].append(pev.pEvtBgoHits().GetHitY(ihit))


        zbgo = array('d', 14*[-999.] )
        for ibgo in xrange(0,14):
            zbgo[ibgo]=58.5+ibgo*29.0           #MARGHERITA&ANTONIO- (Oct 03, 2017)
            maxEneLay = 0
            maxEneLaypos = -9
            maxPosLay = -999
            maxPosX   = -999.
            maxPosY   = -999.
            if len(EBarLayers[ibgo])!= 0 : 
                maxEneLay = np.max(EBarLayers[ibgo])
                maxEneLaypos = np.argmax(EBarLayers[ibgo])
                maxPosLay = PosBarLayers[ibgo][maxEneLaypos]
                maxPosX   = PosXBarLayers[ibgo][maxEneLaypos]
                maxPosY   = PosYBarLayers[ibgo][maxEneLaypos] 
                if(ibgo%2 == 0):                   #MARGHERITA&ANTONIO- (Oct 03, 2017)
                    fBGO_cbgomax[ibgo]=maxPosY     #MARGHERITA&ANTONIO- (Oct 03, 2017)
                if(ibgo%2 == 1):                   #MARGHERITA&ANTONIO- (Oct 03, 2017) 
                    fBGO_cbgomax[ibgo]=maxPosX     #MARGHERITA&ANTONIO- (Oct 03, 2017)

            if(ibgo%2==0):                                                                                                                  #MARGHERITA&ANTONIO- (Oct 03, 2017)
                fBGO_cbgostk[ibgo]=track_impact_point.y()+track_slope.y()*(zbgo[ibgo]-track_impact_point.z())      #MARGHERITA&ANTONIO- (Oct 03, 2017)
            if(ibgo%2==1):                                                                                                                  #MARGHERITA&ANTONIO- (Oct 03, 2017)
                fBGO_cbgostk[ibgo]=track_impact_point.x()+track_slope.x()*(zbgo[ibgo]-track_impact_point.z())          #MARGHERITA&ANTONIO- (Oct 03, 2017)  

            if(maxEneLay > 0.):
                barpos = maxPosX  if  ibgo%2 else  maxPosY
                                
                eCoreMaxLayer[ibgo] = maxEneLay
                eCoreCoordMaxLayer[ibgo] = maxEneLay*barpos

                eCoreLayer[ibgo] = maxEneLay 
                eCoreCoordLayer[ibgo] = maxEneLay*barpos
 

                bar_less1 = -9
                bar_plus1 = -9
                if(PosBarLayers[ibgo][maxEneLaypos] > 0 and PosBarLayers[ibgo][maxEneLaypos] < 21):
                    if((PosBarLayers[ibgo][maxEneLaypos] - 1) in PosBarLayers[ibgo]):
                        bar_less1 = PosBarLayers[ibgo].index(PosBarLayers[ibgo][maxEneLaypos] - 1)
                        if((PosBarLayers[ibgo][maxEneLaypos] + 1) in PosBarLayers[ibgo]):
                            bar_plus1 = PosBarLayers[ibgo].index(PosBarLayers[ibgo][maxEneLaypos] + 1)

                if(bar_less1 != -9):
                    eCoreLayer[ibgo] += EBarLayers[ibgo][bar_less1]
                    eCoreCoordLayer[ibgo] += EBarLayers[ibgo][bar_less1]*PosXBarLayers[ibgo][bar_less1] if ibgo%2 else EBarLayers[ibgo][bar_less1]*PosYBarLayers[ibgo][bar_less1]
                    
                if(bar_plus1 != -9):
                    eCoreLayer[ibgo] += EBarLayers[ibgo][bar_plus1]
                    eCoreCoordLayer[ibgo] += EBarLayers[ibgo][bar_plus1]*PosXBarLayers[ibgo][bar_plus1] if ibgo%2 else EBarLayers[ibgo][bar_plus1]*PosYBarLayers[ibgo][bar_plus1]
                
                
                #print "OK up to here..."
                eCoreCoordLayer[ibgo] /= eCoreLayer[ibgo]
                # ok print "lay ", ibgo,   eCoreCoordLayer[ibgo]
                
                for ibar in xrange(0,len(EBarLayers[ibgo])):
                    rmsLayer[ibgo] += EBarLayers[ibgo][ibar]*(PosXBarLayers[ibgo][ibar]-eCoreCoordLayer[ibgo])*(PosXBarLayers[ibgo][ibar]-eCoreCoordLayer[ibgo]) if ibgo%2 else EBarLayers[ibgo][ibar]*(PosYBarLayers[ibgo][ibar]-eCoreCoordLayer[ibgo])*(PosYBarLayers[ibgo][ibar]-eCoreCoordLayer[ibgo])
                    #print " rms ibar ", PosXBarLayers[ibgo][ibar], PosYBarLayers[ibgo][ibar] , " lay ", ibgo, rmsLayer[ibgo] 
                
                #print rmsLayer[ibgo], len(EBarLayers[ibgo])
            if(maxEneLay > 0.):
                rmsLayer[ibgo]  = np.sqrt(rmsLayer[ibgo]/sum(EBarLayers[ibgo]))
                
            else :
                rmsLayer[ibgo] = 0
            sumRMS += rmsLayer[ibgo]

        Xtr = sumRMS*sumRMS*sumRMS*sumRMS*(sum(EBarLayers[13])/pev.pEvtBgoRec().GetTotalEnergy())/8000000.

        zpsd = [-324.7,-310.7,-298.5,-284.5]
        y_projection_fromBGO_to_PSD_up =  pevspace_api.GetDirectionBGOInteceptY() + zpsd[0]* pevspace_api.GetDirectionBGOSlopeY() #BGO-track projection on PSD YZ view of "high" layer
        
        fBGOtrack_to_PSD_topY[0] = y_projection_fromBGO_to_PSD_up
        fBGO_theta[0] = theta_bgo_deg
        fSTK_theta[0] = theta_stk_deg
        fBGO_phi[0] = phi_bgo_deg
        fSTK_phi[0] = phi_stk_deg

        #### BGO
        fBGO_EnergyG[0] = etot
        #fBGO_EnergyG_SatCorrp[0] = etot_satcorrp
        #fBGO_EnergyG_SatCorrHe[0] = etot_satcorrhe
        #fBGO_EnergyG_SatCorrBCNO[0] = etot_satcorrbcno
        #fBGO_EnergyG_SatCorrNeMgSi[0] = etot_satcorrnemgsi
        #fBGO_EnergyG_SatCorrFe[0] = etot_satcorrfe
        fBGO_EnergyG_Quench[0] = etot_quench
        fBGO_EnergyG_QuenchSatCorr[0] = etot_quench_satcorrZmc
        
        fBGO_HET[0] = pev.pEvtHeader().GeneratedTrigger(3)
        fBGO_LET[0] = pev.pEvtHeader().GeneratedTrigger(4)
        fBGO_UBT[0] = pev.pEvtHeader().GeneratedTrigger(0)
        fBGO_UBE[0] = pev.pEvtHeader().EnabledTrigger(0)
        fBGO_nbars[0] = pev.NEvtBgoHits()
        fBGO_Xtr[0] = Xtr

        ### MARGHERITA&ANTONIO - TIME AND TIMEMS (Oct 13, 2017)
        ftime[0] = pev.pEvtHeader().GetSecond()
        ftimems[0] = pev.pEvtHeader().GetMillisecond()


        #### MC 
        fMC_EnergyT[0] = etot_truth
        fMC_theta[0] = theta_truth
        fMC_phi[0] = phi_truth
        fMC_stopZ[0] = stopZ_parent
        fMC_stopX[0] = stopX_parent
        fMC_stopY[0] = stopY_parent
    
        entry_ntupla += 1

        if(psdchargeX_corr[0] == -100 and psdchargeX_corr[1] == -100):
            print " -100 find!" 
            print sum_len_psd_vec_chargeY, ' ', sum_len_psd_vec_chargeX

        for ibgo in range(14):
            h_bgo_ene_lay[ibgo].Fill((pev.pEvtBgoRec().GetELayer(ibgo))/1000.)

        newtree.Fill()

    dmpch.Terminate()

    
    print "All statistic:   ", h_energy[0].GetEntries()
    for i in range( 1, ncuts):
        percent = h_energy[i].GetEntries()/h_energy[0].GetEntries()*100 
        print cut_string[i], "  n events (all) :  " , h_energy[i].GetEntries() ,  "percent: ",    percent
        if(opts.mc):
            percent_good = h_energy_good[i].GetEntries()/h_energy_good[0].GetEntries()*100 
            print cut_string[i], "  n events (good):  " , h_energy_good[i].GetEntries() ,  "percent: ",    percent_good
            percent_bad = h_energy_bad[i].GetEntries()/h_energy_bad[0].GetEntries()*100 
            print cut_string[i], "  n events (bad) :  " , h_energy_bad[i].GetEntries() ,  "percent: ",    percent_bad

    print "Wrong ID for MC ",wrong_id
    print "Events not passint the new BGO-STK match ", ntrack_newsel
    tf_skim.Write()
    tf_skim.Close()

   
    
    if opts.data:
        oout = root_dir+'/'+opts.diag_output
        tf = TFile(oout, "RECREATE")
        h_terrestrial_lat_vs_long.Write()
        h_bgo_projection_before.Write()
        h_bgo_projection_after.Write()
        h_bgo_projectionTop_before.Write()
        h_bgo_projectionTop_after.Write()
        h_bgo_EratioLay_tot_before.Write()
        h_bgo_EratioLay_tot_after.Write()
        h_bgo_lateralpos_after.Write()
        h_bgo_lateralpos_before.Write()
        h_stk_ntracks_before.Write()
        h_stk_ntracks_after.Write()
        h_stk_XY_check.Write()
        h_stk_IP_XY.Write()
        h_bgo_theta.Write()
        h_bgo_stk_deltaTheta.Write()
        h_bgo_stk_deltaTheta_selection.Write()
        h_bgo_stk_deltaTheta_tracksel.Write()
        h_stk_chi2_norm.Write()
        h_stk_resX_top.Write()
        h_stk_resY_top.Write()
        h_psd_resX[0].Write()
        h_psd_resY[0].Write()
        h_psd_resX[1].Write()
        h_psd_resY[1].Write()
        h_psd_resX_top.Write()
        h_psd_resY_top.Write()
        h_psd_resX_bot.Write()
        h_psd_resY_bot.Write()
        h_stk_resX_top_deltaTheta.Write()
        h_stk_resY_top_deltaTheta.Write()
        h_stk_resX_deltaTheta.Write()
        h_stk_resY_deltaTheta.Write()
        h_stk_resX_top_resBot.Write()
        h_stk_resY_top_resBot.Write()
        h_stk_resX_top_selection.Write()
        h_stk_resY_top_selection.Write()
        h_stk_resX_selection.Write()
        h_stk_resY_selection.Write()      
        h_psd_chargeX.Write()
        h_psd_chargeY.Write()      
        h_psd_chargeX_corr.Write()
        h_psd_chargeY_corr.Write()
        h_psd_GlobalID.Write()
        h_psd_corr_ID512.Write()
        h_psd_corr_ID5504.Write()
        h_stk_to_psd.Write()
        h_psd_energy_hit.Write()
        h_bgo_firedLayer_before.Write()
        h_bgo_ene4lay_before.Write()
        h_bgo_ene4layLast_before.Write()
        h_bgo_ene4Lay_top_vs_bot_before.Write()
        h_bgo_ene4lay_after.Write()
        h_bgo_ene4layLast_after.Write()
        h_bgo_ene4Lay_top_vs_bot_after.Write()        
        h_psd_chargeY_before.Write()
        h_psd_chargeX_before.Write()
        h_psd_chargeY_after.Write()
        h_psd_chargeX_after.Write()

        for i in range(0, ncuts):
            h_energy[i].Write()

        for iplane in range(0,6):
            h_stk_cluster_XvsY[iplane].Write()

        for ibgo in range(14):
            h_bgo_ene_lay[ibgo].Write()


        tf.Write()
        tf.Close()   



if __name__ == "__main__":
    main()
    print datetime.now()-start
        


