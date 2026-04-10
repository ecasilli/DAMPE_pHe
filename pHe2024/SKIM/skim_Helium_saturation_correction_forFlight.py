"""
@brief: dirty skimmer, starting from Stephan's skim -> use Ions & Co
@vale: used for skimming of protons/ions

"""


from sys import argv
import numpy as np
import math
from xrd_utils import isfile
from os.path import isdir, abspath
from os import remove, mkdir, system as os_call
from glob import glob
from argparse import ArgumentParser
from itertools import product as iproduct
from array import array



from datetime import datetime
start=datetime.now()


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
    parser.add_argument("-data",default=False, action="store_true", dest='data', help='use data and output .root ')
    parser.add_argument("-skim",default=False, action="store_true", dest='skim', help='active for skim')
    opts = parser.parse_args(args)

    verbosity = 2002 # ignore ROOT warnings
    if opts.verbose:
        verbosity = 1000
    if opts.quiet:
        verbosity = 3003

    inp = opts.input
    print '====> reading files in  ... ', inp
    stream_name = "ions"

    oout_text_miss = opts.txt_output_missing
    oout_text_skim = opts.skimFile
    ######### log bins defintion #######
    nEDBINs = 160
    EDmin0=10
    EDmax0=1000000
    EDmin = []
    EDmax = []
    EDEdge = [] 
    EDstepX=np.log10(EDmax0/EDmin0)/nEDBINs
    for iedge in range(0, nEDBINs):
        EDEdge.append(EDmin0*pow(10,iedge*EDstepX))
        EDmax.append(EDmin0*pow(10,(iedge+1)*EDstepX))
        #        print EDmin[i], EDmax[i]
        
    EDEdge.append(EDmax[nEDBINs-1])
    Edges= array('d',EDEdge) # this make a bound array for TH1F
    
    BGOz = [58.5,87.5,116.5,145.5,174.5,203.5,232.5,261.5,290.5,319.5,348.5,377.5,406.5,435.5]

    BGOz_X = BGOz[1::2]  
    BGOz_Y  =BGOz[0::2]   
    
    #from ROOT import *
    #import ROOT
    from ROOT import TClonesArray, TFile, TTree, gSystem, gROOT, AddressOf, gInterpreter
    from ROOT import TH2F, TH1F, TMath, TGraphAsymmErrors
    
    gSystem.Load("libDmpEvent.so")
    gSystem.Load("libDmpEventFilter.so")
    
    gSystem.Load("libDmpKernel.so")
    gSystem.Load("libDmpService.so")

#    
    gROOT.SetBatch(True)
    gROOT.ProcessLine( "gErrorIgnoreLevel = %i;"%verbosity)
    from ROOT import DmpChain, DmpEvent, DmpFilterOrbit, DmpPsdBase, DmpCore, DmpVAlg
    from ROOT import DmpSvcPsdEposCor, DmpVSvc   #DmpRecPsdManager
    import DMPSW
#    import 

#    gInterpreter.AddIncludePath('/cvmfs/dampe.cern.ch/centos7/opt/DMPSW/trunk/include/')
#    gInterpreter.ProcessLine('#include "/storage/gpfs_data/dampe/users/falemann/AlgDy2SatCor/include/DmpAlgDy2SatCor.h"')
#    gSystem.Load('/storage/gpfs_data/dampe/users/falemann/AlgDy2SatCor/libDmpAlgDy2SatCor.so')
    ############## create skim file ############
#    if opts.skimFile:
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

    fPSD_cpsdmax	= array('d', 4*[-999.])
    fPSD_cpsdstk	= array('d', 4*[-999.])
    newtree.Branch("PSD_cpsdmax",fPSD_cpsdmax,"PSD_cpsdmax[4]/D")
    newtree.Branch("PSD_cpsdstk",fPSD_cpsdstk,"PSD_cpsdstk[4]/D")


    fPSD_Etrack01	= array('d', [-999.])		
    fPSD_Etrack23	= array('d', [-999.])	
    newtree.Branch("PSD_Etrack01",fPSD_Etrack01,"PSD_Etrack01/D")
    newtree.Branch("PSD_Etrack23",fPSD_Etrack23,"PSD_Etrack23/D")

    
    fPSD_nhits = array('i', [ 0 ] )
    newtree.Branch("PSD_nhits",fPSD_nhits,"PSD_nhits/I")
    
    #### STK

    fSTK_chargeX = array('d', [ -999. ] )
    fSTK_chargeY = array('d', [ -999. ] )
    newtree.Branch("STK_chargeX",fSTK_chargeX,"STK_chargeX/D")
    newtree.Branch("STK_chargeY",fSTK_chargeY,"STK_chargeY/D")
    
    fSTK_ntracks = array('i', [ -999 ] )
    fSTK_nclusters = array('i', [ -999 ] )
    fSTK_theta_correction = array('d', [ -999. ] )
    newtree.Branch("STK_ntracks",fSTK_ntracks,"STK_ntracks/I")
    newtree.Branch("STK_nclusters",fSTK_nclusters,"STK_nclusters/I")
    newtree.Branch("STK_theta_correction",fSTK_theta_correction,"STK_theta_correction/D")

    fSTK_trackIP = array('d', 3*[ -999. ] )
    fSTK_trackDirection = array('d', 2*[ -999. ] ) 
    newtree.Branch("STK_trackIP",fSTK_trackIP,"STK_trackIP[3]/D")
    newtree.Branch("STK_trackDirection",fSTK_trackDirection,"STK_trackDirection[2]/D")

    #### BGO

    fBGO_EnergyG = array('d', [ -999. ] )
    newtree.Branch("BGO_EnergyG",fBGO_EnergyG,"BGO_EnergyG/D")

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

#############################################################################




    #### MC 
    fMC_EnergyT = array('d', [ -999. ] )
    fMC_theta  = array('d', [ -999. ] )
    newtree.Branch("MC_EnergyT",fMC_EnergyT,"MC_EnergyT/D")
    newtree.Branch("MC_theta",fMC_theta,"MC_theta/D")
    fMC_stopZ = array('d', [ -99999. ] )
    fMC_stopX = array('d', [ -99999. ] )
    fMC_stopY = array('d', [ -99999. ] )
    newtree.Branch("MC_stopZ",fMC_stopZ,"MC_stopZ/D")
    newtree.Branch("MC_stopX",fMC_stopX,"MC_stopX/D")
    newtree.Branch("MC_stopY",fMC_stopY,"MC_stopY/D")
       
    
#    newtree.Branch("",f,"/F")
#    newtree.Branch("",f,"/F")
#    newtree.Branch("",f,"/F")
#    newtree.Branch("",f,"/F")


    dmpch = DmpChain("CollectionTree")
    #    dmpch.SetBranchesToDrop("DmpEvtBgoRec_V0");
    files = [f.replace("\n","") for f in open(inp,'r').readlines()]
    
   
    for ifile, f in enumerate(files):
        DMPSW.IOSvc.Set("InData/Read" if ifile == 0 else "InData/ReadMore",f)
        if isfile(f):
            dmpch.Add(f)
            print ifile , f
        else:
            with open(oout_text_miss,'ar') as file:
                file.write(f+'\n')
                continue

        #dmpch.Add(f)
        #print ifile , f
            

#    dmpch.Print()

    if not opts.mc:
        DMPSW.IOSvc.Set("OutData/NoOutput", "True")
        DMPSW.IOSvc.Initialize()
        pFilter = DmpFilterOrbit("EventHeader")
        pFilter.ActiveMe()
    


#    BGOz = [58.5,87.5,116.5,145.5,174.5,203.5,232.5,261.5,290.5,319.5,348.5,377.5,406.5,435.5]


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
#    h_bgo_LRMS_before = TH1F("h_bgo_LRMS_before","h_bgo_LRMS_before",2000,0,200)

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







#    h_bgo_TRMS_vs_LRMS =  TH2F("h_bgo_TRMS_vs_LRMS","h_bgo_TRMS_vs_LRMS",1000,0,1000,1000,0,1000)
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
        #h_energy_truth_box =  TH1F("h_energy_truth_box","h_energy_truth_box", nEDBINs , Edges)
        h_energy_good = []
        h_energy_bad = []
        h_energy_matrix = []
        h_energy_truth = []
        h_energy_truth_weight = []
        h_energy_weight = []
        for ihisto_mc in range(0,ncuts):
   
            # print i 
            h_title_good  = 'h_energy_good_cut_' + str(ihisto_mc)
            h_ene_cut_good = TH1F(h_title_good, h_title_good , nEDBINs , Edges)
            h_energy_good.append(h_ene_cut_good)
            h_title_bad = 'h_energy_bad_cut_' + str(ihisto_mc)
            h_ene_cut_bad = TH1F(h_title_bad, h_title_bad , nEDBINs , Edges)
            h_energy_bad.append(h_ene_cut_bad)
            
            h_energy_matrix_title = "h_energy_matrix_cut_" + str(ihisto_mc)
            h_energy_matrix_histo = TH2F(h_energy_matrix_title,h_energy_matrix_title, nEDBINs , Edges,nEDBINs , Edges)
            h_energy_matrix.append(h_energy_matrix_histo)

            h_energy_truth_title = "h_energy_truth_cut_" + str(ihisto_mc)
            h_energy_truth_histo = TH1F(h_energy_truth_title,h_energy_truth_title, nEDBINs , Edges)
            h_energy_truth.append(h_energy_truth_histo)

            h_energy_truth_weight_title = "h_energy_truth_weight_cut_" + str(ihisto_mc)
            h_energy_truth_weight_histo = TH1F(h_energy_truth_weight_title,h_energy_truth_weight_title, nEDBINs , Edges)
            h_energy_truth_weight_histo.Sumw2()
            h_energy_truth_weight.append(h_energy_truth_weight_histo)

            
            h_energy_weight_title = "h_energy_weight_cut_" + str(ihisto_mc)
            h_energy_weight_histo = TH1F(h_energy_weight_title,h_energy_weight_title, nEDBINs , Edges)
            h_energy_weight_histo.Sumw2()
            h_energy_weight.append(h_energy_weight_histo)
            
            
            #print '*******', len(h_energy_good),len(h_energy_bad)


    
    h_energy_NoTracks = TH1F("h_energy_NoTracks","h_energy_NoTracks",nEDBINs , Edges)
    h_energy_trigger_check = TH1F("h_energy_trigger_check","h_energy_trigger_check",nEDBINs , Edges) 
    h_energy = []
    
    for ihisto in range(0,ncuts):
        # print i 
        h_title = 'h_energy_cut_' + str(ihisto)
        h_ene_cut = TH1F(h_title, h_title , nEDBINs , Edges)
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



    nevents = dmpch.GetEntries()
    #nevents = 5000
    print ' Number of events: ',nevents
    
    ntrack_newsel = 0

 
    dmpch.SetOutputDir(root_dir,stream_name) #
    if opts.skim:
        print ' skim option activeted'



    wrong_id = 0
    entry_ntupla = 0
#    for iev in xrange(0,2000):  ## in the skim the event 0 is not taken into account
    for iev in xrange(0,nevents):  ## in the skim the event 0 is not taken into account
        if opts.mc:
            DmpVSvc.gPsdECor.SetMCflag(1)

        pev=dmpch.GetDmpEvent(iev)
        good_event = False

      #  if (opts.mc and  pev.pEvtSimuPrimaries().pvpart_ekin/1000. < 10000.): continue

#        ## for Litium
#        if opts.mc:
#            parent_particleID = int(pev.pEvtSimuPrimaries().pvpart_pdg/10000.)
##            print '----------- parent: ',parent_particleID
#            if(int(parent_particleID) != 100003):
#                print "*************WRONG PARTICLE ID, particle ID ", pev.pEvtSimuPrimaries().pvpart_pdg
#                wrong_id += 1
#                continue
#            ######## energy check range
#            if(pev.pEvtSimuPrimaries().pvpart_ekin/1000. < 10.):
#                wrong_id += 1
#                continue
#            if(pev.pEvtSimuPrimaries().pvpart_ekin/1000. > 1000.):
#                wrong_id += 1
#                continue
            
 


        stopZ_parent  = -1000
        stopX_parent  = -1000
        stopY_parent  = -1000
        startX_parent = -1000
        startY_parent = -1000
        startZ_parent = -1000
           


       
        ## for helium
        if opts.mc:
            parent_particleID = int(pev.pEvtSimuPrimaries().pvpart_pdg/10000.)
#            print '----------- parent: ',parent_particleID
            if(int(parent_particleID) != 100002):
                print "*************WRONG PARTICLE ID, particle ID ", pev.pEvtSimuPrimaries().pvpart_pdg
                wrong_id += 1
                continue
            ######## energy check range
       #     if(pev.pEvtSimuPrimaries().pvpart_ekin/1000. < 10000.):
        #        wrong_id += 1
         #       continue
          #  if(pev.pEvtSimuPrimaries().pvpart_ekin/1000. > 100000.):
           #     wrong_id += 1
            #    continue
 

        # for protons
#        if opts.mc:
#            parent_particleID = pev.pEvtSimuPrimaries().pvpart_pdg
##            print '----------- parent: ',parent_particleID
#            if(int(parent_particleID) != 2212):
#                print "*************WRONG PARTICLE ID, particle ID ", parent_particleID
#                continue

            
            
            
            #print 'ntrajectories ' ,pev.NSimuTrajectory()
            for it in xrange(0, pev.NSimuTrajectory()):
                if (pev.pSimuTrajectory(0).parentID == 0):
                    stopZ_parent  =  pev.pSimuTrajectory(0).stop_z
                    stopX_parent  =  pev.pSimuTrajectory(0).stop_x
                    stopY_parent  =  pev.pSimuTrajectory(0).stop_y
                    startX_parent = pev.pSimuTrajectory(0).start_x;
                    startY_parent = pev.pSimuTrajectory(0).start_y;
                    startZ_parent = pev.pSimuTrajectory(0).start_z;
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
            interceptXZ_truth = ((startZ_parent*stopX_parent) - (stopZ_parent*startX_parent) )/(startZ_parent-stopZ_parent)
            interceptYZ_truth = ((startZ_parent*stopY_parent) - (stopZ_parent*startY_parent) )/(startZ_parent-stopZ_parent)
            slope_X = (  startX_parent - stopX_parent)/(  startZ_parent - stopZ_parent)
            slope_Y = (  startY_parent - stopY_parent)/(  startZ_parent - stopZ_parent)
            x_projection_truth_onBGO_Top = interceptXZ_truth  + BGOzTop  * slope_X 
            y_projection_truth_onBGO_Top = interceptYZ_truth  + BGOzTop  * slope_Y 
            x_projection_truth_onBGO_Bot = interceptXZ_truth  + BGOzBot  * slope_X 
            y_projection_truth_onBGO_Bot = interceptYZ_truth  + BGOzBot  * slope_Y 
                                                              

            def_good_truth_bgoAcc = abs(x_projection_truth_onBGO_Top) < bgo_acceptance_range and  abs(y_projection_truth_onBGO_Top) < bgo_acceptance_range and abs(x_projection_truth_onBGO_Bot) < bgo_acceptance_range and  abs(y_projection_truth_onBGO_Bot) < bgo_acceptance_range
            


            if(def_good_truth_stopXY==True and def_good_truth_stopZ == True and  def_good_truth_bgoAcc == True):
                good_event = True
		h_stopZ_good.Fill(stopZ_parent)
                h_stopX_good.Fill(stopX_parent)
                h_stopY_good.Fill(stopY_parent)
                h_bgo_projection_good_Top.Fill(x_projection_truth_onBGO_Top,y_projection_truth_onBGO_Top)
                h_bgo_projection_good_Bot.Fill(x_projection_truth_onBGO_Bot,y_projection_truth_onBGO_Bot)
                if(abs(stopX_parent) > 400. or abs(stopY_parent) > 400.):
                    print iev, stopX_parent, stopY_parent


#	from ROOT import DmpAlgDy2SatCor                    
        
        #etot=pev.pEvtBgoRec().GetTotalEnergy()/1000. 
	#etot = DmpAlgDy2SatCor.GetDy2Cor()
	
	from ROOT import DmpEvtBgoRec
	#chargevalue = 2
	etot = pev.pEvtBgoRec().GetDy2SatEcor(2)/1000.
#	if(etot < 0.1): continue
#	print etot
##################################





        etot_truth = -9.
        if opts.mc:
            etot_truth = pev.pEvtSimuPrimaries().pvpart_ekin/1000. 
       
        
        if not opts.mc:
            longitude = 0 #pev.pEvtAttitude().lon_geo -> not saved in Andrii skim
            latitude  = 0 #pev.pEvtAttitude().lat_geo
            
        else:
            longitude= 0
            latitude = 0
            
        h_energy[0].Fill(etot)
        
        if opts.mc:
            h_energy_truth[0].Fill(etot_truth)
            h_energy_matrix[0].Fill(etot_truth,etot)
            h_energy_weight[0].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[0].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[0].Fill(etot)
            else:
                h_energy_bad[0].Fill(etot)
           
#            print " detector flag ", pev.pEvtSimuHeader().GetDetectorFlag()
#            if(pev.pEvtSimuHeader().GetDetectorFlag() == True):
#                h_energy_truth_box.Fill(etot_truth)
                
        ####### Energy cut
        if(etot < 20.): continue
      #  if(etot >= 100.): continue
        h_energy[1].Fill(etot)
        if opts.mc:
            h_energy_truth[1].Fill(etot_truth)
            h_energy_matrix[1].Fill(etot_truth,etot)
            h_energy_weight[1].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[1].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[1].Fill(etot)
            else:
               h_energy_bad[1].Fill(etot)

     
        ###### BGO acceptance at the 3rd last layer of the BGO
        x_projection_fromBGO_to_BGO =  pev.pEvtBgoRec().GetInterceptXZ() + BGOzBot* pev.pEvtBgoRec().GetSlopeXZ()
        y_projection_fromBGO_to_BGO =  pev.pEvtBgoRec().GetInterceptYZ() + BGOzBot* pev.pEvtBgoRec().GetSlopeYZ()

        x_projection_fromBGO_to_BGOTop =  pev.pEvtBgoRec().GetInterceptXZ() +BGOzTop  * pev.pEvtBgoRec().GetSlopeXZ()
        y_projection_fromBGO_to_BGOTop =  pev.pEvtBgoRec().GetInterceptYZ() +BGOzTop  * pev.pEvtBgoRec().GetSlopeYZ()
        
 
        h_bgo_projection_before.Fill(x_projection_fromBGO_to_BGO,y_projection_fromBGO_to_BGO)
        h_bgo_projectionTop_before.Fill(x_projection_fromBGO_to_BGOTop,y_projection_fromBGO_to_BGOTop)
        
        bgo_acceptance = abs(x_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(x_projection_fromBGO_to_BGOTop) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGOTop) < bgo_acceptance_range

#        print 'Bottom projection', x_projection_fromBGO_to_BGO, y_projection_fromBGO_to_BGO
#        print '  Top  projection', x_projection_fromBGO_to_BGOTop, y_projection_fromBGO_to_BGOTop
#        print '  Bool flag      ', bgo_acceptance

        if( (pev.pEvtBgoRec().GetSlopeXZ()==0 and pev.pEvtBgoRec().GetInterceptXZ()==0) or (pev.pEvtBgoRec().GetSlopeYZ()==0 and pev.pEvtBgoRec().GetInterceptYZ()==0)): 
            #print " Event with all 0 ", iev
            continue


        if not (bgo_acceptance): continue
        h_energy[2].Fill(etot)
        if opts.mc:
            h_energy_truth[2].Fill(etot_truth)
            h_energy_matrix[2].Fill(etot_truth,etot)
            h_energy_weight[2].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[2].Fill(etot_truth,math.pow(etot_truth,-1.7))
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
#        #            if(good_event == True): print iev#        

        h_bgo_EratioLay_tot_after.Fill(frac_max_lay)


        h_energy[3].Fill(etot)
        if opts.mc:
            h_energy_truth[3].Fill(etot_truth)
            h_energy_matrix[3].Fill(etot_truth,etot)
            h_energy_weight[3].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[3].Fill(etot_truth,math.pow(etot_truth,-1.7))
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
        #print barNumberMaxEBarLayer
        #print MaxEBarLayer
        
        for jj in xrange(0,len(barNumberMaxEBarLayer)):
            if(barNumberMaxEBarLayer[jj] <= 0 or barNumberMaxEBarLayer[jj] == 21):
                pos_max_lateral = True ### bug solved!

        if(pos_max_lateral == True):
            #  print "entry that should not be here! ", iev
             continue
    

        for ilay in xrange(0,3):
            h_bgo_lateralpos_after.Fill(lateral_pos_vec[ilay])
           
        h_energy[4].Fill(etot)
        if opts.mc:
            h_energy_truth[4].Fill(etot_truth)
            h_energy_matrix[4].Fill(etot_truth,etot)
            h_energy_weight[4].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[4].Fill(etot_truth,math.pow(etot_truth,-1.7))
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
            h_energy_truth[5].Fill(etot_truth)
            h_energy_matrix[5].Fill(etot_truth,etot)
            h_energy_weight[5].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[5].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[5].Fill(etot)
            else:
                h_energy_bad[5].Fill(etot)
                
            
        
        
        
        if(pev.pEvtHeader().GeneratedTrigger(3) == True or pev.pEvtHeader().GeneratedTrigger(0) == True ):
            h_energy_trigger_check.Fill(etot)

        # HET trigger or UB or LET
            #if(pev.pEvtHeader().GeneratedTrigger(3) == False and pev.pEvtHeader().GeneratedTrigger(0) == False and pev.pEvtHeader().GeneratedTrigger(4) == False): continue
        
        
        h_energy[6].Fill(etot)
        if opts.mc:
            h_energy_truth[6].Fill(etot_truth)
            h_energy_matrix[6].Fill(etot_truth,etot)
            h_energy_weight[6].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[6].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[6].Fill(etot)
            else:
                h_energy_bad[6].Fill(etot)
                

        # ntrcks != 0 cut
        
        ntracks = pev.NStkKalmanTrack()
        if(ntracks < 0): 
            print "###ERROR!!!!! on ntracks ",ntracks
            break
        h_stk_ntracks_before.Fill(ntracks)
        if(ntracks == 0): 
            #print "0 tracks, entry:  ",iev ," energy: ",etot 
            h_energy_NoTracks.Fill(etot)
            continue
        
        h_stk_ntracks_after.Fill(ntracks)
        h_energy[7].Fill(etot)
        if opts.mc:
            h_energy_truth[7].Fill(etot_truth)
            h_energy_matrix[7].Fill(etot_truth,etot)
            h_energy_weight[7].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[7].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[7].Fill(etot)
            else:
                h_energy_bad[7].Fill(etot)



        # BGO-STK-PSD match
        is_match_with_PSDY = False
        is_match_with_PSDX = False
        is_match_with_STKY = False
        is_match_with_STKX = False
        
        res_x_min = 1000.
        res_y_min = 1000.
        trackID_X = -9
        trackID_Y = -9
        
        tgZ = math.atan(np.sqrt( (pev.pEvtBgoRec().GetSlopeXZ()*pev.pEvtBgoRec().GetSlopeXZ()) + (pev.pEvtBgoRec().GetSlopeYZ()*pev.pEvtBgoRec().GetSlopeYZ()) ) );
        theta_bgo = tgZ*180./math.pi 

        vec_track_IDX = []
        vec_track_IDY = []

        vec_track_resX = []
        vec_track_resY = []

#        passed_chi2 = False

#        if(iev == 54):
#            print 'ntracks ', ntracks

#        print 'start new event', iev
        for itrack in xrange(0,ntracks):
            track = pev.pStkKalmanTrack(itrack)
            chi2_norm = track.getChi2()/(track.getNhitX()+track.getNhitY()-4)
            has_plane0_clusterX = False
            has_plane0_clusterY = False
            h_stk_chi2_norm.Fill(chi2_norm)


#            if(iev == 54): 
#                print "chi2", itrack, chi2_norm
            

            if(chi2_norm > 25.): continue

            for iclu in xrange(0,track.GetNPoints()):
                clux = track.pClusterX(iclu)
                cluy = track.pClusterY(iclu)
                if (clux and clux.getPlane() == 0):
                    has_plane0_clusterX = True
                if (cluy and cluy.getPlane() == 0):
                    has_plane0_clusterY = True
                    
                # check plot for the dead region of STK
                if(clux and cluy):
                    if(clux.getPlane() == cluy.getPlane()):
                        h_stk_cluster_XvsY[clux.getPlane()].Fill(clux.GetX(),cluy.GetY())

            if(has_plane0_clusterX and not has_plane0_clusterY):
                h_stk_XY_check.Fill(1)
            if(has_plane0_clusterY and not has_plane0_clusterX):
                h_stk_XY_check.Fill(2)
            if(has_plane0_clusterY and  has_plane0_clusterX):
                h_stk_XY_check.Fill(3)


#            if(iev == 54): 
#                print has_plane0_clusterY, has_plane0_clusterX
            if(has_plane0_clusterY == False and  has_plane0_clusterX == False):  continue

            ######loop on track
            
            theta_track =math.acos(track.getDirection().CosTheta())*180./math.pi
        #    print 'theta bgo', theta_bgo
            deltaTheta_rec = theta_bgo - theta_track
            h_bgo_stk_deltaTheta.Fill(deltaTheta_rec)
                        
            # print track.getImpactPoint().x(), track.getImpactPoint().y(), track.getImpactPoint().z() 
            h_stk_IP_XY.Fill(track.getImpactPoint().x(),track.getImpactPoint().y())

            # charge distribution => flat skim for RooFit plots
            # residual at the TRACK IMPACT POINT
            res_stk_bgo_x_top = track.getImpactPoint().x() - (pev.pEvtBgoRec().GetInterceptXZ() + track.getImpactPoint().z() * pev.pEvtBgoRec().GetSlopeXZ())
            res_stk_bgo_y_top = track.getImpactPoint().y() - (pev.pEvtBgoRec().GetInterceptYZ() + track.getImpactPoint().z() * pev.pEvtBgoRec().GetSlopeYZ())



            h_stk_resX_top.Fill(res_stk_bgo_x_top)
            h_stk_resY_top.Fill(res_stk_bgo_y_top)
            
            # residual at the TOP OF BGO
            track_projection_x = (track.getDirection().x()*(BGOzTop - track.getImpactPoint().z()) + track.getImpactPoint().x())
            track_projection_y = (track.getDirection().y()*(BGOzTop - track.getImpactPoint().z()) + track.getImpactPoint().y())
            
                
            res_stk_bgo_x = x_projection_fromBGO_to_BGOTop - track_projection_x
            res_stk_bgo_y = y_projection_fromBGO_to_BGOTop - track_projection_y


                
            h_stk_resX.Fill(res_stk_bgo_x)
            h_stk_resY.Fill(res_stk_bgo_y)

#            if(iev == 151 or iev == 888):
#                print 'Ev. ' , iev, '  delta theta ',deltaTheta_rec
            
            if(abs(deltaTheta_rec) > 25.): continue
            
            h_stk_resX_top_deltaTheta.Fill(res_stk_bgo_x_top)
            h_stk_resY_top_deltaTheta.Fill(res_stk_bgo_y_top)
            h_stk_resX_deltaTheta.Fill(res_stk_bgo_x)
            h_stk_resY_deltaTheta.Fill(res_stk_bgo_y)



#            if(abs(res_stk_bgo_x) > 200.): continue
#            if(abs(res_stk_bgo_y) > 200.): continue
            h_stk_resX_top_resBot.Fill(res_stk_bgo_x_top)
            h_stk_resY_top_resBot.Fill(res_stk_bgo_y_top)


            
            if(abs(res_stk_bgo_x) < 40. and abs(res_stk_bgo_y)< 40. and abs(res_stk_bgo_x_top)< 40. and abs(res_stk_bgo_y_top)<40.):
                h_bgo_stk_deltaTheta_selection.Fill(deltaTheta_rec)
                #                if(abs(deltaTheta_rec) > 20.):
                #                    print "delta Theta ", deltaTheta_rec, " event ", iev 
                

            if(abs(deltaTheta_rec) < 5):
                h_stk_resX_top_selection.Fill(res_stk_bgo_x_top)
                h_stk_resY_top_selection.Fill(res_stk_bgo_y_top)
                h_stk_resX_selection.Fill(res_stk_bgo_x)
                h_stk_resY_selection.Fill(res_stk_bgo_y)

                
            if(abs(res_stk_bgo_x_top) < 200. and abs(res_stk_bgo_x) < 60.): ##
                is_match_with_STKX = True
                vec_track_IDX.append(itrack)
                vec_track_resX.append(res_x_min)

                if( res_x_min > abs(res_stk_bgo_x_top)):
                    res_x_min = abs(res_stk_bgo_x_top)
                    trackID_X = itrack
                
            if(abs(res_stk_bgo_y_top) < 200. and abs(res_stk_bgo_y) < 60.): ## 
                is_match_with_STKY = True
                vec_track_IDY.append(itrack)
                vec_track_resY.append(res_y_min)

                if( res_y_min > abs(res_stk_bgo_y_top)):
                    res_y_min = abs(res_stk_bgo_y_top)
                    trackID_Y = itrack
                    
            
#            if(iev == 39):
#                print itrack, res_stk_bgo_x , res_stk_bgo_y, chi2_norm
#                print vec_track_IDX, vec_track_IDY


        
##### deltaTheta < 25 degree and residual to the top of the BGO < 60 mm

#            if(iev == 9948):
#                print iev, 'residuals ', res_stk_bgo_y_top, res_stk_bgo_x_top
#            print vec_track_IDY, vec_track_IDX
#            for ipsd in xrange(0,pev.NEvtPsdHits()):
#                test = np.sqrt(pev.pEvtPsdHits().fEnergy[ipsd]/2.)
#                print "charge ",test, "pos", pev.pEvtPsdHits().GetHitX(ipsd),  pev.pEvtPsdHits().GetHitY(ipsd), pev.pEvtPsdHits().GetHitZ(ipsd) 

    


        if(trackID_X == -9): continue
        if(trackID_Y == -9): continue
        
       
        track_ID = -9
#        if(iev == 54):
#            print trackID_X, trackID_Y
        ######## both tracks are needed for PSD charge identification! 
        if(trackID_X == -9 or trackID_Y == -9): 
            print "ERROR!!!!! Wrong track id!!!! both at -9? break"
            break

        if(trackID_X == trackID_Y):
            track_ID = trackID_X
        else:
            
            trackX = pev.pStkKalmanTrack(trackID_X)
            trackY = pev.pStkKalmanTrack(trackID_Y)
            chi2KX = trackX.getChi2() /(trackX.getNhitX()+trackX.getNhitY()-4)
	    chi2KY = trackY.getChi2() /(trackY.getNhitX()+trackY.getNhitY()-4)
            npointX = trackX.GetNPoints()
            npointY = trackY.GetNPoints()
            
            
            
          #  print 'chi2     ', chi2KX,  chi2KY
          #  print 'n points ', npointX, npointY
            if(npointX == npointY or abs(npointX - npointY) == 1):
                if(chi2KX < chi2KY):
                    #print trackID_X, trackID_Y
                    #print  "X", vec_track_IDX , "Y",  vec_track_IDY
                    if trackID_X in vec_track_IDY:
                        track_ID = trackID_X
                    elif trackID_Y in vec_track_IDX:
                            track_ID = trackID_Y
                    else:
                        common_id = list(set(vec_track_IDX).intersection(vec_track_IDY))
                        #    print common_id
                        res_min = 1000.
                        for ids in xrange(0, len(common_id)):
                            #print common_id[ids]
                            pos_X = vec_track_IDX.index(common_id[ids])
                            pos_Y = vec_track_IDY.index(common_id[ids])
                            #print "position in vec ", pos_X , pos_Y
                            res_tot = abs(vec_track_resX[pos_X]) + abs(vec_track_resY[pos_Y])
                            #print " res tot ",  res_tot
                            if(res_min > res_tot):
                                res_min = res_tot
                                #print "min ", res_tot
                                track_ID = common_id[ids]
                        #print "res X ", vec_track_resX, "res Y ",vec_track_resY
                    #print "results", track_ID
                    #print ""
                else:
                    if trackID_Y in vec_track_IDX:
                        track_ID = trackID_Y
                    elif trackID_X in vec_track_IDY:
                            track_ID = trackID_X
                    else:
                        common_id = list(set(vec_track_IDX).intersection(vec_track_IDY))
                        res_min = 1000.
                        for ids in xrange(0, len(common_id)):
                            pos_X = vec_track_IDX.index(common_id[ids])
                            pos_Y = vec_track_IDY.index(common_id[ids])
                            res_tot = abs(vec_track_resX[pos_X]) + abs(vec_track_resY[pos_Y])
                            if(res_min > res_tot):
                                res_min = res_tot
                                track_ID = common_id[ids]


            else:
                if(npointX > npointY):
                    if trackID_X in vec_track_IDY:
                        track_ID = trackID_X
                    elif trackID_Y in vec_track_IDX:
                            track_ID = trackID_Y
                    else:
                        common_id = list(set(vec_track_IDX).intersection(vec_track_IDY))
                        res_min = 1000.
                        for ids in xrange(0, len(common_id)):
                            pos_X = vec_track_IDX.index(common_id[ids])
                            pos_Y = vec_track_IDY.index(common_id[ids])
                            res_tot = abs(vec_track_resX[pos_X]) + abs(vec_track_resY[pos_Y])
                            if(res_min > res_tot):
                                res_min = res_tot
                                track_ID = common_id[ids]

                    
                else:
                    if trackID_Y in vec_track_IDX:
                        track_ID = trackID_Y
                    elif trackID_X in vec_track_IDY:
                            track_ID = trackID_X
                    else:
                        common_id = list(set(vec_track_IDX).intersection(vec_track_IDY))
                        res_min = 1000.
                        for ids in xrange(0, len(common_id)):
                            pos_X = vec_track_IDX.index(common_id[ids])
                            pos_Y = vec_track_IDY.index(common_id[ids])
                            res_tot = abs(vec_track_resX[pos_X]) + abs(vec_track_resY[pos_Y])
                            if(res_min > res_tot):
                                res_min = res_tot
                                track_ID = common_id[ids]
                    


                                
        if(track_ID == -9): 
            #print " event not passing the match of track id ", iev, "delta theta rec ",deltaTheta_rec
            ntrack_newsel += 1 
            continue

        h_energy[8].Fill(etot)
        if opts.mc:
            h_energy_truth[8].Fill(etot_truth)
            h_energy_matrix[8].Fill(etot_truth,etot)
            h_energy_weight[8].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[8].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[8].Fill(etot)
            else:
                h_energy_bad[8].Fill(etot)



#        print 'selected track ID', track_ID, '\n'
        track_sel = pev.pStkKalmanTrack(track_ID)
        
        

        theta_track_sel =math.acos(track_sel.getDirection().CosTheta())*180./math.pi
        #    print 'theta bgo', theta_bgo
        deltaTheta_rec_sel = theta_bgo - theta_track_sel
        h_bgo_stk_deltaTheta_tracksel.Fill(deltaTheta_rec_sel)

        theta_truth  = -9.
        if opts.mc:
            theta_truth  = math.atan(np.sqrt( (pev.pEvtSimuPrimaries().pvpart_px/pev.pEvtSimuPrimaries().pvpart_pz*pev.pEvtSimuPrimaries().pvpart_px/pev.pEvtSimuPrimaries().pvpart_pz)+(pev.pEvtSimuPrimaries().pvpart_py/pev.pEvtSimuPrimaries().pvpart_pz*pev.pEvtSimuPrimaries().pvpart_py/pev.pEvtSimuPrimaries().pvpart_pz)))*180./math.pi    
            #print "truth ", theta_truth, " stk ", theta_track, "bgo ", theta_bgo
            delta_theta_bgo = theta_truth - theta_bgo
            delta_theta_stk = theta_truth - theta_track
            h_bgo_deltaTheta_bgo.Fill(delta_theta_bgo)
            h_bgo_deltaTheta_stk.Fill(delta_theta_stk)
        



        #### projection on PSD for charge selection
        # cluster Charge X and Y
        track_correction = track_sel.getDirection().CosTheta()


        cluster_chargeX = -1000.
        cluster_chargeY = -1000.

        for iclu in xrange(0,track_sel.GetNPoints()):
            clux = track_sel.pClusterX(iclu)
            cluy = track_sel.pClusterY(iclu)
            if (clux and clux.getPlane() == 0):
                cluster_chargeX = clux.getEnergy()*track_correction
                #cluster_cor_energy = DmpStkTrackHelper.CorrEnergy(pev.pStkKalmanTrack(track_ID), track_sel.pClusterX(iclu))
                #print " stk correction ", cluster_chargeX , cluster_cor_energy


            if (cluy and cluy.getPlane() == 0):
                cluster_chargeY = cluy.getEnergy()*track_correction


        h_stk_clusterX.Fill(cluster_chargeX)
        h_stk_clusterY.Fill(cluster_chargeY)        



        ########## PSD FIDUCIAL VOLUME CUT
        #double PSD_ZY[NLAYERSPSD]  ={ -324.7, -310.7};
        #double PSD_ZX[NLAYERSPSD]  ={ -298.5,-284.5};
        

        psd_YZ_top = -324.7
        psd_XZ_top = -298.5
        stk_to_psd_topY = (track_sel.getDirection().y()*(psd_YZ_top - track_sel.getImpactPoint().z()) + track_sel.getImpactPoint().y())
        stk_to_psd_topX = (track_sel.getDirection().x()*(psd_XZ_top - track_sel.getImpactPoint().z()) + track_sel.getImpactPoint().x())
        

#        if(iev == 54):
#            print stk_to_psd_topX, stk_to_psd_topY



        if(abs(stk_to_psd_topX) > 400.): continue
        if(abs(stk_to_psd_topY) > 400.): continue
        
        
        h_energy[9].Fill(etot)
        if opts.mc:
            h_energy_truth[9].Fill(etot_truth)
            h_energy_matrix[9].Fill(etot_truth,etot)
            h_energy_weight[9].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[9].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[9].Fill(etot)
            else:
                h_energy_bad[9].Fill(etot)


                
        ######### STK - PSD MATCH #############

        # psd Charge X and Y

        
        psd_vec_chargeX     = [[]for _ in range(2)]
        psd_vec_gidX        = [[]for _ in range(2)]
        psd_vec_pathlengthX = [[]for _ in range(2)]
        psd_vec_positionX   = [[]for _ in range(2)]

        psd_vec_chargeY     = [[]for _ in range(2)]
        psd_vec_gidY        = [[]for _ in range(2)]
        psd_vec_pathlengthY = [[]for _ in range(2)]
        psd_vec_positionY   = [[]for _ in range(2)]

#        print "****Event ", iev, " second ", pev.pEvtHeader().GetSecond(), " milli second ", pev.pEvtHeader().GetMillisecond()
#        print " track ID ", track_ID
#        print " pathlengh X ", psd_vec_pathlengthX
#        print " theta track ", theta_track

        ######### STK - PSD MATCH using the new tool #############
        ################# QUI LE MODIFICHE ################### ADB


	
#	offset = [0.0,0.0,0.4,0.6]
        for ipsd in xrange(0,pev.NEvtPsdHits()):
            if(pev.pEvtPsdHits().IsHitMeasuringX(ipsd)):
                crossingX = False
                lenghtX = [-99999.,-99999.]
                array_lenghtX = array('d',lenghtX)
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -298.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -298.51): #-298.5 ):_
                    npsdX = 0
	        if(pev.pEvtPsdHits().GetHitZ(ipsd) < -284.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -284.51):
                    npsdX = 1
		
                
                if opts.mc:
                    crossingX = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtX)
                else:
                    crossingX = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtX)

#                if (iev == 54):
#                    print 'crossing X ',crossingX, pev.pEvtPsdHits().GetHitX(ipsd), np.sqrt(pev.pEvtPsdHits().fEnergy[ipsd]/2.)

                #        print "is crossing X bar? ", crossingX, " glob bar ID ", pev.pEvtPsdHits().fGlobalBarID[ipsd]
                if(crossingX == True):
                    psd_vec_chargeX[npsdX].append(pev.pEvtPsdHits().fEnergy[ipsd]) 
                    psd_vec_gidX[npsdX].append(pev.pEvtPsdHits().fGlobalBarID[ipsd]) 
                    psd_vec_pathlengthX[npsdX].append(array_lenghtX[1]) 
#                    print " X bar ",npsdX, array_lenghtX[1], pev.pEvtPsdHits().GetHitZ(ipsd), " glob bar ID", pev.pEvtPsdHits().fGlobalBarID[ipsd]
                    psd_vec_positionX[npsdX].append(pev.pEvtPsdHits().GetHitX(ipsd))
                    
                    

            elif(pev.pEvtPsdHits().IsHitMeasuringY(ipsd)):
                crossingY = False
                lenghtY = [-99999.,-99999.]
                array_lenghtY = array('d',lenghtY)
                
                if(pev.pEvtPsdHits().GetHitZ(ipsd) < -324.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -324.71):
                    npsdY = 0 
		if(pev.pEvtPsdHits().GetHitZ(ipsd) < -310.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -310.71):
                    npsdY = 1
		
                    
                if opts.mc:
                    crossingY = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtY)
                else:
                    crossingY = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtY)
#                if (iev == 54):
#                    print 'crossing Y ',crossingY, pev.pEvtPsdHits().GetHitY(ipsd), np.sqrt(pev.pEvtPsdHits().fEnergy[ipsd]/2.)

                if (crossingY == True):
                    psd_vec_chargeY[npsdY].append(pev.pEvtPsdHits().fEnergy[ipsd]) 
                    psd_vec_gidY[npsdY].append(pev.pEvtPsdHits().fGlobalBarID[ipsd]) 
                    psd_vec_pathlengthY[npsdY].append(array_lenghtY[1]) 
                    psd_vec_positionY[npsdY].append(pev.pEvtPsdHits().GetHitY(ipsd))
                 #   print " Y bar ",npsdY, array_lenghtY[1], pev.pEvtPsdHits().GetHitZ(ipsd), " glob bar ID", pev.pEvtPsdHits().fGlobalBarID[ipsd]

#        print " psd vec posX  ", psd_vec_positionX, " psd vec posY " ,psd_vec_positionY
        
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
            print "track selected", track_ID, track_sel.getImpactPoint().x(), track_sel.getImpactPoint().y(), track_sel.getImpactPoint().z()

        
#        if(iev == 17):
#            print psd_vec_chargeX, psd_vec_chargeX
#            print psd_vec_pathlengthX, psd_vec_pathlengthY
#            print psd_vec_positionX, psd_vec_positionY
            

#        print " after loop the psd vec ", psd_vec_pathlengthX
#        print ""


        psdchargeX = [-999.,-999.]
        psdchargeY = [-999.,-999.]

        psdchargeX_corr = [-999.,-999.]
        psdchargeY_corr = [-999.,-999.]
        
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
                lenghtY = [-99999.,-99999.]
                array_lenghtY = array('d',lenghtY)
                test_pos = False 
                if opts.mc:
                    test_pos = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidY[ipsd][pos_max_len],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtY)
                else:
                    test_pos = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidY[ipsd][pos_max_len],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtY)
                 
                    
                PsdEC_tmpY = -1.
                if(test_pos == True):
                    PsdEC_tmpY = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidY[ipsd][pos_max_len], array_lenghtY[0]/10.)
                else:
                    print "****** ERROR! it should be true!!!"

                psdchargeY[ipsd] = psd_vec_chargeY[ipsd][pos_max_len]
                psdchargeY_corr[ipsd] = psd_vec_chargeY[ipsd][pos_max_len]*PsdEC_tmpY
                psdchargeY_proj[ipsd] = array_lenghtY[0]
                psdY_pathlength[ipsd] = array_lenghtY[1]
                psdY_position[ipsd] =  psd_vec_positionY[ipsd][pos_max_len]


                h_psd_GlobalID.Fill(psd_vec_gidY[ipsd][pos_max_len])
                if(psdchargeY[ipsd] > 0.):
                    h_psd_chargeY.Fill(np.sqrt(psdchargeY[ipsd]/2.))
                if(psdchargeY_corr[ipsd] > 0):
                    h_psd_chargeY_corr.Fill(np.sqrt(psdchargeY_corr[ipsd]/2.))
                if(psd_vec_gidY[ipsd][pos_max_len] == 512 and array_lenghtY[0]!= -99999.):
                    h_psd_corr_ID512.Fill(array_lenghtY[0],PsdEC_tmpY)

            if(len(psd_vec_chargeX[ipsd]) > 0 ):
                pos_max_len = np.argmax(psd_vec_pathlengthX[ipsd])
                lenghtX = [-99999.,-99999.]
                array_lenghtX = array('d',lenghtX)
                test_posX = False
                if opts.mc:
                    test_posX = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidX[ipsd][pos_max_len],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtX)
                else:
                    test_posX = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidX[ipsd][pos_max_len],track_sel.getDirection(),track_sel.getImpactPoint(), array_lenghtX)
                PsdEC_tmpX = -1.
                if(test_posX == True):
                    PsdEC_tmpX = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidX[ipsd][pos_max_len],array_lenghtX[0]/10.)
                

                psdchargeX[ipsd] = psd_vec_chargeX[ipsd][pos_max_len]
                psdchargeX_corr[ipsd] = psd_vec_chargeX[ipsd][pos_max_len]*PsdEC_tmpX
                psdchargeX_proj[ipsd] = array_lenghtX[0] 
                psdX_pathlength[ipsd] = array_lenghtX[1]
                psdX_position[ipsd] =  psd_vec_positionX[ipsd][pos_max_len]

                h_psd_GlobalID.Fill(psd_vec_gidX[ipsd][pos_max_len])
                if(psdchargeX[ipsd] > 0.):
                    h_psd_chargeX.Fill(np.sqrt(psdchargeX[ipsd]/2.))          
                if(psdchargeX_corr[ipsd] > 0.):
                    h_psd_chargeX_corr.Fill(np.sqrt(psdchargeX_corr[ipsd]/2.))    
                if(psd_vec_gidX[ipsd][pos_max_len] == 5504 and array_lenghtX[0] != -99999.):
                    h_psd_corr_ID5504.Fill(array_lenghtX[0],PsdEC_tmpX)
                    
 

        h_bgo_theta.Fill(theta_bgo)
        

        ########## contine from here with the ntupla construction

        
        sum_len_psd_vec_chargeX = len(psd_vec_chargeX[0]) + len(psd_vec_chargeX[1])
        sum_len_psd_vec_chargeY = len(psd_vec_chargeY[0]) + len(psd_vec_chargeY[1])

        for ipsd in xrange(0,2):
            if(psdchargeY_corr[ipsd] > 0.):
                h_psd_chargeY_before.Fill(np.sqrt(psdchargeY_corr[ipsd]/2.))
            if(psdchargeX_corr[ipsd] > 0.):
                h_psd_chargeX_before.Fill(np.sqrt(psdchargeX_corr[ipsd]/2.))


#        if(iev == 54):
#            print "psd charge ",psd_vec_chargeX, psd_vec_chargeY
        
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
            h_energy_truth[10].Fill(etot_truth)
            h_energy_matrix[10].Fill(etot_truth,etot)
            h_energy_weight[10].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[10].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[10].Fill(etot)
            else:
                h_energy_bad[10].Fill(etot)




        ###### HITS ON ALL THE BGO LAYERS 
        ###### nlayers > 14
        fired_layer = np.count_nonzero(v_bgolayer)
        h_bgo_firedLayer_before.Fill(fired_layer)

#        if(fired_layer != 14): continue

        h_energy[11].Fill(etot)
        if opts.mc:
            h_energy_truth[11].Fill(etot_truth)
            h_energy_matrix[11].Fill(etot_truth,etot)
            h_energy_weight[11].Fill(etot,math.pow(etot_truth,-1.7))
            h_energy_truth_weight[11].Fill(etot_truth,math.pow(etot_truth,-1.7))
            if(good_event == True):
                h_energy_good[11].Fill(etot)
            else:
                h_energy_bad[11].Fill(etot)


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
#       offset = [0.0,0.0,0.4,0.6]
	zpsd = [-324.7,-310.7,-298.5,-284.5]													#M e A 17-10-17
	                                    
	
        for ipsd in xrange(0,nlayer_psd):
            fPSD_EnergyY[ipsd] = psdchargeY[ipsd]
            fPSD_EnergyX[ipsd] = psdchargeX[ipsd]
	    
	    if (ipsd < 2):															#M e A 17-10-17
 	    	fPSD_cpsdstk[ipsd]=track_sel.getImpactPoint().y()+track_sel.getDirection().y()*(zpsd[ipsd]-track_sel.getImpactPoint().z())  	#M e A 17-10-17
            if (ipsd > 1):                                                                                                   			#M e A 17-10-17
            	fPSD_cpsdstk[ipsd]=track_sel.getImpactPoint().x()+track_sel.getDirection().x()*(zpsd[ipsd]-track_sel.getImpactPoint().z())  	#M e A 17-10-17   
 
           
	    fPSD_EnergyY_corr[ipsd] = psdchargeY_corr[ipsd]
            fPSD_EnergyX_corr[ipsd] = psdchargeX_corr[ipsd]
            fPSD_psdX_projY[ipsd] = psdchargeX_proj[ipsd]
            fPSD_psdY_projX[ipsd] = psdchargeY_proj[ipsd]
            fPSD_psdX_pathlength[ipsd] = psdX_pathlength[ipsd]
            fPSD_psdY_pathlength[ipsd] = psdY_pathlength[ipsd] 
            fPSD_psdX_position[ipsd] = psdX_position[ipsd]
            fPSD_psdY_position[ipsd] = psdY_position[ipsd]
	    
#	    if(ipsd<2):
#			 fPSD_cpsdmax[ipsd] = pev.pEvtPsdHits().GetHitY(ipsd) 
#	    if(ipsd>1):
#			 fPSD_cpsdmax[ipsd] = pev.pEvtPsdHits().GetHitX(ipsd)+offset[ipsd]	   


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



       
#        print  psdchargeY[0], psdchargeY[1], psdchargeX[0], psdchargeX[1]
#        print  'max X ',fPSD_psdX_Emax_layer, ' max Y ', fPSD_psdY_Emax_layer
#        print "npsd hits ", fPSD_nhits[0]
        #             
        #           
        
        #### STK

        fSTK_chargeX[0] = cluster_chargeX
        fSTK_chargeY[0] = cluster_chargeY
        fSTK_ntracks[0] = pev.NStkKalmanTrack()
        fSTK_nclusters[0] =  pev.NStkSiCluster()
        fSTK_theta_correction[0] = track_correction
        fSTK_trackIP[0] = track_sel.getImpactPoint().x()
        fSTK_trackIP[1] = track_sel.getImpactPoint().y()
        fSTK_trackIP[2] = track_sel.getImpactPoint().z() #track_sel.getImpactPoint().x()   MARGHERITA&ANTONIO (we changed x->z)
        fSTK_trackDirection[0] = track_sel.getDirection().x()
        fSTK_trackDirection[1] = track_sel.getDirection().y()

        


             
        for ibgo in xrange(0,14):
            
            fBGO_EneLay[ibgo] = pev.pEvtBgoRec().GetELayer(ibgo)/1000.
           # if(CorLyrE[ibgo]>fBGO_EneLay[ibgo]): fBGO_EneLay[ibgo] = CorLyrE[ibgo] #saturation correction
            
            fBGO_nhits_layers[ibgo] = pev.pEvtBgoRec().GetLayerHits()[ibgo]

            #print "Raw Lyr: ", pev.pEvtBgoRec().GetELayer(ibgo)/1000., "Corr Lyr:", CorLyrE[ibgo], "Lyr:", ibgo
	
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

 #       test_meno9 = False

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
            zbgo[ibgo]=58.5+ibgo*29.0                               #MARGHERITA&ANTONIO- (Oct 03, 2017)
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
            	if(ibgo%2 == 0):                                    #MARGHERITA&ANTONIO- (Oct 03, 2017)
			fBGO_cbgomax[ibgo]=maxPosY                  #MARGHERITA&ANTONIO- (Oct 03, 2017)
            	if(ibgo%2 == 1):                                    #MARGHERITA&ANTONIO- (Oct 03, 2017)	
                	fBGO_cbgomax[ibgo]=maxPosX                  #MARGHERITA&ANTONIO- (Oct 03, 2017)

            if(ibgo%2==0):                                                                                                                  #MARGHERITA&ANTONIO- (Oct 03, 2017)
		fBGO_cbgostk[ibgo]=track_sel.getImpactPoint().y()+track_sel.getDirection().y()*(zbgo[ibgo]-track_sel.getImpactPoint().z())  	#MARGHERITA&ANTONIO- (Oct 03, 2017)
	    if(ibgo%2==1):                                                                                                                  #MARGHERITA&ANTONIO- (Oct 03, 2017)
                fBGO_cbgostk[ibgo]=track_sel.getImpactPoint().x()+track_sel.getDirection().x()*(zbgo[ibgo]-track_sel.getImpactPoint().z())  	    #MARGHERITA&ANTONIO- (Oct 03, 2017)  
     
     

            




#            if(maxEneLaypos == -9):
#                test_meno9 = True
#                print "Entry ",iev
#                print maxEneLay, maxEneLaypos , maxPosLay, maxPosX, maxPosY

            if(maxEneLay > 0.):
                barpos = maxPosX  if  ibgo%2 else  maxPosY
                                
                eCoreMaxLayer[ibgo] = maxEneLay
                eCoreCoordMaxLayer[ibgo] = maxEneLay*barpos

                
                

                eCoreLayer[ibgo] = maxEneLay 
                eCoreCoordLayer[ibgo] = maxEneLay*barpos
 

               
                

                bar_less1 = -9
                bar_plus1 = -9
 #               print "Bar number in layer ", PosBarLayers[ibgo][maxEneLaypos]
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
                


               #print " ecore coord ", eCoreCoordLayer[ibgo]
                    
                    #print bar_less1, bar_plus1
                
                
                
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
#            print  ibgo, " rms layer ",rmsLayer[ibgo], sum(EBarLayers[ibgo])
#        print " sumRMS ",sumRMS

        Xtr = sumRMS*sumRMS*sumRMS*sumRMS*(sum(EBarLayers[13])/pev.pEvtBgoRec().GetTotalEnergy())/8000000.
#        if(test_meno9 == True):
#            print iev, Xtr

                #### BGO
        fBGO_EnergyG[0] = etot
        
        fBGO_HET[0] = pev.pEvtHeader().GeneratedTrigger(3)
        fBGO_LET[0] = pev.pEvtHeader().GeneratedTrigger(4)
        fBGO_UBT[0] = pev.pEvtHeader().GeneratedTrigger(0)
        fBGO_UBE[0] = pev.pEvtHeader().EnabledTrigger(0)
        fBGO_nbars[0] = pev.NEvtBgoHits()
        fBGO_Xtr[0] = Xtr

#        print "Entry ", iev, " Xtr ", Xtr


### MARGHERITA&ANTONIO - TIME AND TIMEMS (Oct 13, 2017)


	ftime[0] = pev.pEvtHeader().GetSecond()
	ftimems[0] = pev.pEvtHeader().GetMillisecond()





        #### MC 
            

        fMC_EnergyT[0] = etot_truth
        fMC_theta[0] = theta_truth
        fMC_stopZ[0] = stopZ_parent
        fMC_stopX[0] = stopX_parent
        fMC_stopY[0] = stopY_parent
    

        entry_ntupla += 1

        #if(iev>=40 and iev <=60):
#        if(etot > 500.):
#            print "event passing ", iev


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

        
#        h_psd_resX.Write()
#        h_psd_resY.Write()

#      h_stk_clusterX.Write()
#      h_stk_clusterY.Write()        
      
#      h_psd_chargeX_size.Write()
#      h_psd_chargeY_size.Write()
      
        h_psd_chargeX.Write()
        h_psd_chargeY.Write()
      
        h_psd_chargeX_corr.Write()
        h_psd_chargeY_corr.Write()
      
        

  #    h_energy_NoTracks.Write()

  #    h_energy_trigger_check.Write()

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
        

#        h_bgo_LRMS_before.Write()
        
        h_psd_chargeY_before.Write()
        h_psd_chargeX_before.Write()
        h_psd_chargeY_after.Write()
        h_psd_chargeX_after.Write()

        for i in range(0, ncuts):
            h_energy[i].Write()
      

      #h_psd_chargeX_stk_selection.Write()   
      #h_psd_chargeY_stk_selection.Write()
      #h_psd_chargeYvsE_stk_selection.Write()

        for iplane in range(0,6):
            h_stk_cluster_XvsY[iplane].Write()

        for ibgo in range(14):
            h_bgo_ene_lay[ibgo].Write()


        tf.Write()
        tf.Close()
########## comment up to here ################
#     if opts.skimFile:

 


            #   
#   
            

         



if __name__ == "__main__":
    main()
    print datetime.now()-start
        


