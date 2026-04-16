"""
@brief: dirty skimmer, starting from Stephan's skim -> use Ions & Co
@vale: used for skimming of protons/ions

Refactor goals:
- one code path that works for both MC and flight data
- SAA cut only on flight data
- quenching handled only for MC and only if the branch exists
- cut number from documented "DAMPE p+He spectrum update"
"""

from sys import argv
import numpy as np
import math
from os.path import isdir, abspath
from os import mkdir
from argparse import ArgumentParser
from array import array
from datetime import datetime

# ML tools allocation
from dampeapi import DampeApi
from dampe_bgo_sat_api import DampeBgoApi

from ROOT import *

start = datetime.now()

e_min = 20.  # GeV


def mc_particle_info(mc_particle_name):
    particles = {
        "p": [2212, 1],
        "He": [100002, 2],
        "C": [100006, 6],
        "N": [100007, 7],
        "O": [100008, 8],
        "F": [100009, 9],
        "Ne": [100010, 10],
        "Na": [100011, 11],
        "Mg": [100012, 12],
        "Al": [100013, 13],
        "Si": [100014, 14],
        "P": [100015, 15],
        "S": [100016, 16],
        "Cl": [100017, 17],
        "Ar": [100018, 18],
        "K": [100019, 19],
        "Ca": [100020, 20],
        "Fe": [100026, 26],
    }
    if mc_particle_name not in particles:
        raise ValueError(
            "Wrong particle name given in arguments. Available: {}".format(
                ", ".join(sorted(particles.keys()))
            )
        )
    return particles[mc_particle_name]


def spectral_weight(E, gamma, edges):
    # weight for a power-law spectrum normalized in [Emin, Emax] from edges
    bin_idx = np.searchsorted(edges, E, side='right') - 1
    if bin_idx < 0:
        bin_idx = 0
    if bin_idx >= len(edges) - 1:
        bin_idx = len(edges) - 2
    e1 = edges[bin_idx]
    e2 = edges[bin_idx + 1]
    norm = (gamma - 1.) * math.log(e2 / e1) / (e1 ** (1. - gamma) - e2 ** (1. - gamma))
    return norm * E ** (1. - gamma)


def build_log_edges(edmin, edmax, nbins):
    edges = []
    step = np.log10(edmax / edmin) / nbins
    for iedge in range(nbins):
        edges.append(edmin * pow(10, iedge * step))
    edges.append(edmin * pow(10, nbins * step))
    return array('d', edges)


def make_branch(tree, name, typecode, length=None, init_value=None):
    if length is None:
        buf = array(typecode, [init_value])
        leaf = "%s/%s" % (name, typecode.upper())
    else:
        buf = array(typecode, [init_value] * length)
        leaf = "%s[%d]/%s" % (name, length, typecode.upper())
    tree.Branch(name, buf, leaf)
    return buf


def reset_array(buf, value):
    for i in xrange(len(buf)):
        buf[i] = value


def has_branch(chain, branch_name):
    branches = chain.GetListOfBranches()
    if not branches:
        return False
    for branch in branches:
        if branch.GetName() == branch_name:
            return True
    return False


def fill_energy_cut(idx, etot, etot_truth, good_event, is_mc,
                    h_energy, h_energy_matrix,
                    h_energy_truth, h_energy_weight, h_energy_truth_weight,
                    energy_binnings, spectral_indices,
                    h_energy_good, h_energy_bad):
    h_energy[idx].Fill(etot)

    if not is_mc:
        return

    h_energy_matrix[idx].Fill(etot_truth, etot)
    for bname, (nbins, edges) in energy_binnings.items():
        h_energy_truth[idx][bname].Fill(etot_truth)
        for tag, gamma in spectral_indices.items():
            weight = spectral_weight(etot_truth, gamma, edges)
            h_energy_weight[idx][bname][tag].Fill(etot, weight)
            h_energy_truth_weight[idx][bname][tag].Fill(etot_truth, weight)

    if good_event:
        h_energy_good[idx].Fill(etot)
    else:
        h_energy_bad[idx].Fill(etot)


def initialize_or_rebind_apis(iev, is_mc, bgorec, stkclusters,
                              stkladders, runmetadata, pevspace_api,
                              bgo_api):
    if is_mc and iev == 0:
        pevspace_api.Initialize(False, bgorec, stkclusters, stkladders, config="ions")
        bgo_api.Initialize(bgorec)
        #if has_quenching and bgoquenchrec:
            #bgo_api_quench.Initialize(bgoquenchrec)
    elif (not is_mc) and iev == 0:
        pevspace_api.Initialize(True, bgorec, stkclusters, None, runmetadata, config="ions")
        bgo_api.Initialize(bgorec)
    elif is_mc and iev > 0:
        pevspace_api.Rebind(bgorec, stkclusters, stkladders)
        bgo_api.BindBgoRec(bgorec)
        #if has_quenching and bgoquenchrec:
            #bgo_api_quench.BindBgoRec(bgoquenchrec)
    else:
        pevspace_api.Rebind(bgorec, stkclusters, 0, runmetadata)
        bgo_api.BindBgoRec(bgorec)


def compute_bgo_energies(pev, bgo_api, is_mc, has_quenching, mc_particle_Z):
    raw_total_mev = pev.pEvtBgoRec().GetTotalEnergy()
    etot = raw_total_mev / 1000.

    etot_satcorr_ions = etot
    etot_satcorr_ions_v3 = etot
    if raw_total_mev > 1e6:
        bgo_api.Predict(mc_particle_Z, method="ions")
        if bgo_api.IsSaturated():
            etot_satcorr_ions = bgo_api.GetReconstructedBGOE() / 1000.

        bgo_api.Predict(mc_particle_Z, method="ions_v3")
        if bgo_api.IsSaturated():
            etot_satcorr_ions_v3 = bgo_api.GetReconstructedBGOE() / 1000.

    # Quenching branch: only for MC and only if available
    etot_quench = etot
    etot_quench_satcorr_ions = etot
    etot_quench_satcorr_ions_v3 = etot
    if is_mc and has_quenching:
        try:
            quench_rec = pev.pEvtBgoQuenchRec()
            etot_quench = quench_rec.GetTotalEnergy() / 1000.
            etot_quench_satcorr_ions = etot_quench
            etot_quench_satcorr_ions_v3 = etot_quench
            if quench_rec.GetTotalEnergy() > 1e6:
                bgo_api.Predict(mc_particle_Z, method="ions") # Z is a dummy value here. N.B. without quenching!
                if bgo_api.IsSaturated():
                    etot_quench_satcorr_ions = bgo_api.GetReconstructedBGOE() / 1000.

                bgo_api.Predict(mc_particle_Z, method="ions_v3") # Z is a dummy value here. New ML SATcorr (0.4.0-beta.1) N.B. without quenching!
                if bgo_api.IsSaturated():
                    etot_quench_satcorr_ions_v3 = bgo_api.GetReconstructedBGOE() / 1000.
        
        except Exception:
            etot_quench = etot
            etot_quench_satcorr_ions = etot_satcorr_ions
            etot_quench_satcorr_ions_v3 = etot_satcorr_ions_v3

    return (
        etot,
        etot_satcorr_ions,
        etot_satcorr_ions_v3,
        etot_quench,
        etot_quench_satcorr_ions,
        etot_quench_satcorr_ions_v3,
    )


def main(args=None):
    parser = ArgumentParser(usage="Usage: %(prog)s [options]", description="run dirty skimmer (for protons/ions)")
    parser.add_argument("-i", "--input", dest='input', help='the input files to be processed')
    parser.add_argument("-v", "--verbose", action='store_true', default=False, dest='verbose', help='run in high verbosity mode')
    parser.add_argument("-q", "--quiet", action='store_true', default=False, dest='quiet', help='suppress a lot of output, quiet mode')
    parser.add_argument("-o", "--diagnostics_output", default="diagnostics.root", type=str, dest='diag_output', help='name of diagnostics file to be written')
    parser.add_argument("-O", "--output", default="skimmed_data/", type=str, dest='output', help='name of output directory')
    parser.add_argument("-s", "--skimFile", default="skim_out.root", type=str, dest='skimFile', help='skim files for roofit')
    parser.add_argument("-e", "--text_output_missing", default="livetime_missing.txt", type=str, dest='txt_output_missing', help='name of text output file with missing runs')
    parser.add_argument("-mc", default=False, action="store_true", dest='mc', help='use this flag for MC data')
    parser.add_argument("-pid", "--mc_particle_name", default="p", type=str, dest="mc_particle_name", help="particle simulated in the MC file you want to skim")
    parser.add_argument("-data", default=False, action="store_true", dest='data', help='use data and output .root')
    parser.add_argument("-skim", default=False, action="store_true", dest='skim', help='active for skim')
    opts = parser.parse_args(args)

    is_mc = opts.mc
    is_data = opts.data
    if is_mc == is_data:
        raise RuntimeError("Specify exactly one mode: -mc or -data")

    verbosity = 2002
    if opts.verbose:
        verbosity = 1000
    if opts.quiet:
        verbosity = 3003

    mc_pdg_particle_id = None
    mc_particle_Z = None
    if is_mc:
        mc_pdg_particle_id, mc_particle_Z = mc_particle_info(opts.mc_particle_name)
        print 'mc_pdg_particle_id ', mc_pdg_particle_id, ', mc_particle_Z ', mc_particle_Z

    inp = opts.input
    print '====> reading files in  ... ', inp
    stream_name = "ions"

    oout_text_miss = opts.txt_output_missing
    oout_text_skim = opts.skimFile

    nEDBINs_6bins = 36
    nEDBINs_3bins = 18
    nEDBINs_7bins = 42
    nEDBINs_5bins = 30
    nEDBINs_4bins = 24
    nEDBINs_8bins = 48
    nEDBINs_120bins = 720
    EDmin0 = 10
    EDmax0 = 10000000

    Edges_3bins = build_log_edges(EDmin0, EDmax0, nEDBINs_3bins)
    Edges_4bins = build_log_edges(EDmin0, EDmax0, nEDBINs_4bins)
    Edges_5bins = build_log_edges(EDmin0, EDmax0, nEDBINs_5bins)
    Edges_6bins = build_log_edges(EDmin0, EDmax0, nEDBINs_6bins)
    Edges_7bins = build_log_edges(EDmin0, EDmax0, nEDBINs_7bins)
    Edges_8bins = build_log_edges(EDmin0, EDmax0, nEDBINs_8bins)
    Edges_120bins = build_log_edges(EDmin0, EDmax0, nEDBINs_120bins)

    BGOz = [58.5, 87.5, 116.5, 145.5, 174.5, 203.5, 232.5, 261.5, 290.5, 319.5, 348.5, 377.5, 406.5, 435.5]
    BGOz_X = BGOz[1::2]
    BGOz_Y = BGOz[0::2]

    spectral_indices = {
        "E2e3": 2.3,
        "E2e4": 2.4,
        "E2e5": 2.5,
        "E2e6": 2.6,
        "E2e7": 2.7,
        "E2e8": 2.8,
        "E2e9": 2.9,
        "E3e0": 3.0,
    }
    energy_binnings = {
        "3bins": (nEDBINs_3bins, Edges_3bins),
        "4bins": (nEDBINs_4bins, Edges_4bins),
        "5bins": (nEDBINs_5bins, Edges_5bins),
        "6bins": (nEDBINs_6bins, Edges_6bins),
        "7bins": (nEDBINs_7bins, Edges_7bins),
        "8bins": (nEDBINs_8bins, Edges_8bins),
        "120bins": (nEDBINs_120bins, Edges_120bins),
    }

    pevspace_api = DampeApi()
    bgo_api = DampeBgoApi()
    #bgo_api_quench = DampeBgoApi()

    gSystem.Load("libDmpEvent.so")
    gSystem.Load("libDmpEventFilter.so")
    gSystem.Load("libDmpKernel.so")
    gSystem.Load("libDmpService.so")

    gROOT.SetBatch(True)
    gROOT.ProcessLine("gErrorIgnoreLevel = %i;" % verbosity)
    import DMPSW

    root_dir = abspath(opts.output)
    print " skimming files "
    print '* * will create output here: %s' % root_dir
    if not isdir(root_dir):
        mkdir(root_dir)

    file_skim_out = root_dir + '/' + oout_text_skim
    print file_skim_out
    tf_skim = TFile(file_skim_out, "RECREATE")
    newtree = TTree("newtree", "test")

    nlayer_psd = 2
    nlayer_bgo = 14

    # -----------------------------
    # Branches
    # -----------------------------
    fPSD_EnergyY = make_branch(newtree, 'PSD_EnergyY', 'd', 2, -999.)
    fPSD_EnergyX = make_branch(newtree, 'PSD_EnergyX', 'd', 2, -999.)
    fPSD_EnergyY_corr = make_branch(newtree, 'PSD_EnergyY_corr', 'd', 2, -999.)
    fPSD_EnergyX_corr = make_branch(newtree, 'PSD_EnergyX_corr', 'd', 2, -999.)
    fPSD_psdX_projY = make_branch(newtree, 'PSD_psdX_projY', 'd', 2, -999.)
    fPSD_psdY_projX = make_branch(newtree, 'PSD_psdY_projX', 'd', 2, -999.)
    fPSD_psdX_position = make_branch(newtree, 'PSD_psdX_position', 'd', 2, -999.)
    fPSD_psdY_position = make_branch(newtree, 'PSD_psdY_position', 'd', 2, -999.)
    fPSD_psdX_pathlength = make_branch(newtree, 'PSD_psdX_pathlength', 'd', 2, -999.)
    fPSD_psdY_pathlength = make_branch(newtree, 'PSD_psdY_pathlength', 'd', 2, -999.)
    fPSD_psdX_Emax_layer = make_branch(newtree, 'PSD_psdX_Emax_layer', 'd', None, -999.)
    fPSD_psdY_Emax_layer = make_branch(newtree, 'PSD_psdY_Emax_layer', 'd', None, -999.)
    fPSD_EmaxY = make_branch(newtree, 'PSD_EmaxY', 'd', 2, -999.)
    fPSD_EmaxX = make_branch(newtree, 'PSD_EmaxX', 'd', 2, -999.)
    fPSD_cpsdmax = make_branch(newtree, 'PSD_cpsdmax', 'd', 4, -999.)
    fPSD_cpsdstk = make_branch(newtree, 'PSD_cpsdstk', 'd', 4, -999.)
    fPSD_Etrack01 = make_branch(newtree, 'PSD_Etrack01', 'd', None, -999.)
    fPSD_Etrack23 = make_branch(newtree, 'PSD_Etrack23', 'd', None, -999.)
    fpsdchargeX0 = make_branch(newtree, 'PSD_ChargeX0', 'd', None, -999.)
    fpsdchargeX1 = make_branch(newtree, 'PSD_ChargeX1', 'd', None, -999.)
    fpsdchargeY0 = make_branch(newtree, 'PSD_ChargeY0', 'd', None, -999.)
    fpsdchargeY1 = make_branch(newtree, 'PSD_ChargeY1', 'd', None, -999.)
    fPSD_Global_Charge = make_branch(newtree, 'PSD_Global_Charge', 'd', None, -999.)
    fPSD_Global_Charge_PathAverage = make_branch(newtree, 'PSD_Global_Charge_PathAverage', 'd', None, -999.)
    fPSD_PathWeighted_Charge = make_branch(newtree, 'PSD_PathWeighted_Charge', 'd', None, -999.)
    fPSD_nhits = make_branch(newtree, 'PSD_nhits', 'i', None, 0)

    fSTK_chargeX = make_branch(newtree, 'STK_chargeX', 'd', 6, -999.)
    fSTK_chargeY = make_branch(newtree, 'STK_chargeY', 'd', 6, -999.)
    fSTK_chargeX_etaCorr = make_branch(newtree, 'STK_chargeX_etaCorr', 'd', 6, -999.)
    fSTK_chargeY_etaCorr = make_branch(newtree, 'STK_chargeY_etaCorr', 'd', 6, -999.)
    fSTK_ntrack = make_branch(newtree, 'STK_ntrack', 'i', None, -999)
    fSTK_nclusters = make_branch(newtree, 'STK_nclusters', 'i', None, -999)
    fSTK_theta_correction = make_branch(newtree, 'STK_theta_correction', 'd', None, -999.)
    fSTK_trackIP = make_branch(newtree, 'STK_trackIP', 'd', 3, -999.)
    fSTK_trackDirection = make_branch(newtree, 'STK_trackDirection', 'd', 2, -999.)
    fSTK_vertexPrediction = make_branch(newtree, 'STK_vertexPrediction', 'd', None, -999.)

    fBGO_EnergyG = make_branch(newtree, 'BGO_EnergyG', 'd', None, -999.)
    fBGO_EnergyG_SatCorr_ML_ions = make_branch(newtree, 'BGO_EnergyG_SatCorr_ML_ions', 'd', None, -999.)
    fBGO_EnergyG_SatCorr_ML_ions_v3 = make_branch(newtree, 'BGO_EnergyG_SatCorr_ML_ions_v3', 'd', None, -999.)
    
    fBGO_EnergyG_Quench = None
    fBGO_EnergyG_QuenchSatCorr_ML_ions = None
    fBGO_EnergyG_QuenchSatCorr_ML_ions_v3 = None
    if is_mc:
        fBGO_EnergyG_Quench = make_branch(newtree, 'BGO_EnergyG_Quench', 'd', None, -999.)
        fBGO_EnergyG_QuenchSatCorr_ML_ions = make_branch(newtree, 'BGO_EnergyG_QuenchSatCorr_ML_ions', 'd', None, -999.)
        fBGO_EnergyG_QuenchSatCorr_ML_ions_v3 = make_branch(newtree, 'BGO_EnergyG_QuenchSatCorr_ML_ions_v3', 'd', None, -999.)

    fBGO_HET = make_branch(newtree, 'BGO_HET', 'i', None, 0)
    fBGO_LET = make_branch(newtree, 'BGO_LET', 'i', None, 0)
    fBGO_UBT = make_branch(newtree, 'BGO_UBT', 'i', None, 0)
    fBGO_UBE = make_branch(newtree, 'BGO_UBE', 'i', None, 0)
    fBGO_EneLay = make_branch(newtree, 'BGO_EneLay', 'd', 14, -999.)
    fBGO_nhits_layers = make_branch(newtree, 'BGO_nhits_layers', 'i', 14, -999)
    fBGO_nbars = make_branch(newtree, 'BGO_nbars', 'i', None, -999)
    fBGO_Xtr = make_branch(newtree, 'BGO_xtr', 'd', None, -999.)
    fBGO_cbgomax = make_branch(newtree, 'BGO_cbgomax', 'd', 14, -999.)
    fBGO_cbgostk = make_branch(newtree, 'BGO_cbgostk', 'd', 14, -999.)

    ftime = make_branch(newtree, 'time', 'i', None, 0)
    ftimems = make_branch(newtree, 'timems', 'i', None, 0)
    flong = make_branch(newtree, 'long', 'd', None, -9999.)
    flat = make_branch(newtree, 'lat', 'd', None, -9999.)

    fSTKtrack_to_PSD_topY = make_branch(newtree, 'STKtrack_to_PSD_topY', 'd', None, -999.)
    fSTKtrack_to_PSD_topX = make_branch(newtree, 'STKtrack_to_PSD_topX', 'd', None, -999.)
    fBGOtrack_to_PSD_topY = make_branch(newtree, 'BGOtrack_to_PSD_topY', 'd', None, -999.)
    fSTK_theta = make_branch(newtree, 'STK_theta', 'd', None, -999.)
    fBGO_theta = make_branch(newtree, 'BGO_theta', 'd', None, -999.)
    fSTK_phi = make_branch(newtree, 'STK_phi', 'd', None, -999.)
    fBGO_phi = make_branch(newtree, 'BGO_phi', 'd', None, -999.)

    fBGO_slopeXZ_analy = make_branch(newtree, 'BGO_slopeXZ_analy', 'd', None, -999.)
    fBGO_slopeYZ_analy = make_branch(newtree, 'BGO_slopeYZ_analy', 'd', None, -999.)
    fBGO_interceptXZ_analy = make_branch(newtree, 'BGO_interceptXZ_analy', 'd', None, -999.)
    fBGO_interceptYZ_analy = make_branch(newtree, 'BGO_interceptYZ_analy', 'd', None, -999.)
    fBGO_theta_analy = make_branch(newtree, 'BGO_theta_analy', 'd', None, -999.)
    fBGO_phi_analy = make_branch(newtree, 'BGO_phi_analy', 'd', None, -999.)

    fBGO_slopeX_ML = make_branch(newtree, 'BGO_slopeX_ML', 'd', None, -999.)
    fBGO_slopeY_ML = make_branch(newtree, 'BGO_slopeY_ML', 'd', None, -999.)
    fBGO_interceptX_ML = make_branch(newtree, 'BGO_interceptX_ML', 'd', None, -999.)
    fBGO_interceptY_ML = make_branch(newtree, 'BGO_interceptY_ML', 'd', None, -999.)

    fSTK_slopeX_ML = make_branch(newtree, 'STK_slopeX_ML', 'd', None, -999.)
    fSTK_slopeY_ML = make_branch(newtree, 'STK_slopeY_ML', 'd', None, -999.)
    fSTK_interceptX_ML = make_branch(newtree, 'STK_interceptX_ML', 'd', None, -999.)
    fSTK_interceptY_ML = make_branch(newtree, 'STK_interceptY_ML', 'd', None, -999.)
    fSTK_interceptZ_ML = make_branch(newtree, 'STK_interceptZ_ML', 'd', None, -999.)

    fMC_EnergyT = fMC_theta = fMC_phi = None
    fMC_stopZ = fMC_stopX = fMC_stopY = None
    if is_mc:
        fMC_EnergyT = make_branch(newtree, 'MC_EnergyT', 'd', None, -999.)
        fMC_theta = make_branch(newtree, 'MC_theta', 'd', None, -999.)
        fMC_phi = make_branch(newtree, 'MC_phi', 'd', None, -999.)
        fMC_stopZ = make_branch(newtree, 'MC_stopZ', 'd', None, -99999.)
        fMC_stopX = make_branch(newtree, 'MC_stopX', 'd', None, -99999.)
        fMC_stopY = make_branch(newtree, 'MC_stopY', 'd', None, -99999.)

    # -----------------------------
    # Input chains
    # -----------------------------
    dmpch = DmpChain("CollectionTree")
    metach = TChain("RunMetadataTree")

    files = [f.replace("\n", "") for f in open(inp, 'r').readlines()]
    for ifile, f in enumerate(files):
        DMPSW.IOSvc.Set("InData/Read" if ifile == 0 else "InData/ReadMore", f)
        dmpch.Add(f)
        metach.Add(f)
        print ifile, f

    has_quenching = is_mc and has_branch(dmpch, "DmpEvtBgoQuenchRec")
    print "BGO Quenching:", "ON" if has_quenching else "OFF"

    pFilter = None
    if not is_mc:
        DMPSW.IOSvc.Set("OutData/NoOutput", "True")
        DMPSW.IOSvc.Initialize()
        pFilter = DmpFilterOrbit("EventHeader")
        pFilter.ActiveMe()

    runmetadata = DmpRunMetadata()
    metach.SetBranchAddress("RunMetadata", runmetadata)

    BGOzTop = 46.
    BGOzBot = 448.
    bgo_acceptance_range = 280.

    # -----------------------------
    # Histograms
    # -----------------------------
    h_terrestrial_lat_vs_long = TH2F("h_terrestrial_lat_vs_long", "h_terrestrial_lat_vs_long", 360, 0, 360, 180, -90, 90)
    ncuts = 12
    cut_string = ['All', 'etot > 20 GeV ', 'BGO Acceptance', 'EratioLay', 'BGOLateral', 'SAA', 'HET', 'NTracks', 'BGO-STK', 'PSD Fiducial', 'STK-PSD', 'NLayers']
    h_bgo_projection_before = TH2F("h_bgo_projection_before", "h_bgo_projection_before", 1000, -500, 500, 1000, -500, 500)
    h_bgo_projection_after = TH2F("h_bgo_projection_after", "h_bgo_projection_after", 1000, -500, 500, 1000, -500, 500)
    h_bgo_projectionTop_before = TH2F("h_bgo_projectionTop_before", "h_bgo_projectionTop_before", 1000, -500, 500, 1000, -500, 500)
    h_bgo_projectionTop_after = TH2F("h_bgo_projectionTop_after", "h_bgo_projectionTop_after", 1000, -500, 500, 1000, -500, 500)
    h_bgo_EratioLay_tot_before = TH1F("h_bgo_EratioLay_tot_before", "h_bgo_EratioLay_tot_before", 100, 0, 1)
    h_bgo_EratioLay_goo_before = TH1F("h_bgo_EratioLay_goo_before", "h_bgo_EratioLay_goo_before", 100, 0, 1)
    h_bgo_EratioLay_bad_before = TH1F("h_bgo_EratioLay_bad_before", "h_bgo_EratioLay_bad_before", 100, 0, 1)
    h_bgo_EratioLay_tot_after = TH1F("h_bgo_EratioLay_tot_after", "h_bgo_EratioLay_tot_after", 100, 0, 1)
    h_bgo_lateralpos_after = TH1F("h_bgo_lateralpos_after", "h_bgo_lateralpos_after", 25, -1.5, 23.5)
    h_bgo_lateralpos_before = TH1F("h_bgo_lateralpos_before", "h_bgo_lateralpos_before", 25, -1.5, 23.5)
    h_stk_ntracks_before = TH1F("h_stk_ntracks_before", "h_stk_ntracks_before", 1001, -0.5, 1000.5)
    h_stk_ntracks_after = TH1F("h_stk_ntracks_after", "h_stk_ntracks_after", 1001, -0.5, 1000.5)
    h_stk_XY_check = TH1F("h_stk_XY_check", "1 = cluster X; 2 = cluster Y; 3 = cluster X and Y", 5, 0.5, 5.5)
    h_stk_IP_XY = TH2F("h_stk_IP_XY", "h_stk_IP_XY", 1000, -500, 500, 1000, -500, 500)
    h_stk_resX_top = TH1F("h_stk_resX_top", "h_stk_resX_top", 200, -1000, 1000)
    h_stk_resY_top = TH1F("h_stk_resY_top", "h_stk_resY_top", 200, -1000, 1000)
    h_stk_chi2_norm = TH1F("h_stk_chi2_norm", "h_stk_chi2_norm", 1000, 0, 1000)
    h_stk_resX = TH1F("h_stk_resX", "h_stk_resX", 200, -1000, 1000)
    h_stk_resY = TH1F("h_stk_resY", "h_stk_resY", 200, -1000, 1000)
    h_stk_resX_top_deltaTheta = TH1F("h_stk_resX_top_deltaTheta", "h_stk_resX_top_deltaTheta", 200, -1000, 1000)
    h_stk_resY_top_deltaTheta = TH1F("h_stk_resY_top_deltaTheta", "h_stk_resY_top_deltaTheta", 200, -1000, 1000)
    h_stk_resX_deltaTheta = TH1F("h_stk_resX_deltaTheta", "h_stk_resX_deltaTheta", 200, -1000, 1000)
    h_stk_resY_deltaTheta = TH1F("h_stk_resY_deltaTheta", "h_stk_resY_deltaTheta", 200, -1000, 1000)
    h_stk_resX_top_resBot = TH1F("h_stk_resX_top_resBot", "h_stk_resX_top_resBot", 200, -1000, 1000)
    h_stk_resY_top_resBot = TH1F("h_stk_resY_top_resBot", "h_stk_resY_top_resBot", 200, -1000, 1000)
    h_stk_resX_top_selection = TH1F("h_stk_resX_top_selection", "h_stk_resX_top_selection", 200, -1000, 1000)
    h_stk_resY_top_selection = TH1F("h_stk_resY_top_selection", "h_stk_resY_top_selection", 200, -1000, 1000)
    h_stk_resX_selection = TH1F("h_stk_resX_selection", "h_stk_resX_selection", 200, -1000, 1000)
    h_stk_resY_selection = TH1F("h_stk_resY_selection", "h_stk_resY_selection", 200, -1000, 1000)
    h_stk_clusterX = TH1F("h_stk_clusterX", "h_stk_clusterX", 10000, 0, 10000)
    h_stk_clusterY = TH1F("h_stk_clusterY", "h_stk_clusterY", 10000, 0, 10000)
    h_psd_chargeX_size = TH1F("h_psd_chargeX_size", "h_psd_chargeX_size", 10, 0, 10)
    h_psd_chargeY_size = TH1F("h_psd_chargeY_size", "h_psd_chargeY_size", 10, 0, 10)
    h_psd_chargeX = TH1F("h_psd_chargeX", "h_psd_chargeX", 300, 0, 30)
    h_psd_chargeY = TH1F("h_psd_chargeY", "h_psd_chargeY", 300, 0, 30)
    h_psd_chargeX_corr = TH1F("h_psd_chargeX_corr", "h_psd_chargeX_corr", 300, 0, 30)
    h_psd_chargeY_corr = TH1F("h_psd_chargeY_corr", "h_psd_chargeY_corr", 300, 0, 30)
    h_psd_energy_hit = TH1F("h_psd_energy_hit", "h_psd_energy_hit", 1000, 0, 100)
    h_psd_chargeX_stk_selection = TH1F("h_psd_chargeX_stk_selection", "h_psd_chargeX_stk_selection", 300, 0, 30)
    h_psd_chargeY_stk_selection = TH1F("h_psd_chargeY_stk_selection", "h_psd_chargeY_stk_selection", 300, 0, 30)
    h_psd_chargeYvsE_stk_selection = TH2F("h_psd_chargeYvsE_stk_selection", "h_psd_chargeYvsE_stk_selection", 6000, 0, 6, 300, 0, 30)
    h_bgo_theta = TH1F("h_bgo_theta", "h_bgo_theta", 110, -10, 100)
    h_bgo_stk_deltaTheta = TH1F("h_bgo_stk_deltaTheta", "h_bgo_stk_deltaTheta", 400, -100, 100)
    h_bgo_stk_deltaTheta_selection = TH1F("h_bgo_stk_deltaTheta_selection", "h_bgo_stk_deltaTheta_selection", 400, -100, 100)
    h_bgo_stk_deltaTheta_tracksel = TH1F("h_bgo_stk_deltaTheta_tracksel", "h_bgo_stk_deltaTheta_tracksel", 400, -100, 100)
    h_psd_GlobalID = TH1F("h_psd_GlobalID", "h_psd_GlobalID", 7000, 0, 7000)
    h_psd_corr_ID512 = TH2F("h_psd_corr_ID512", "h_psd_corr_ID512", 1000, -500, 500, 200, 0, 2.)
    h_psd_corr_ID5504 = TH2F("h_psd_corr_ID5504", "h_psd_corr_ID5504", 1000, -500, 500, 200, 0, 2.)
    h_stk_to_psd = TH2F("h_stk_to_psd", "h_stk_to_psd", 1000, -1000, 1000, 1000, -1000, 1000)
    h_bgo_firedLayer_before = TH1F("h_bgo_firedLayer_before", "h_bgo_firedLayer_before", 15, 0, 15)
    h_bgo_ene4lay_before = TH1F("h_bgo_ene4lay_before", "h_bgo_ene4lay_before", 120, 0, 1.2)
    h_bgo_ene4layLast_before = TH1F("h_bgo_ene4layLast_before", "h_bgo_ene4layLast_before", 120, 0, 1.2)
    h_bgo_ene4Lay_top_vs_bot_before = TH2F("h_bgo_ene4Lay_top_vs_bot_before", "h_bgo_ene4Lay_top_vs_bot_before", 120, 0, 1.2, 120, 0, 1.2)
    h_bgo_ene4lay_after = TH1F("h_bgo_ene4lay_after", "h_bgo_ene4lay_after", 120, 0, 1.2)
    h_bgo_ene4layLast_after = TH1F("h_bgo_ene4layLast_after", "h_bgo_ene4layLast_after", 120, 0, 1.2)
    h_bgo_ene4Lay_top_vs_bot_after = TH2F("h_bgo_ene4Lay_top_vs_bot_after", "h_bgo_ene4Lay_top_vs_bot_after", 120, 0, 1.2, 120, 0, 1.2)
    h_psd_chargeY_before = TH1F("h_psd_chargeY_before", "h_psd_chargeY PSDX or PSDY", 1200, 0, 30)
    h_psd_chargeX_before = TH1F("h_psd_chargeX_before", "h_psd_chargeX PSDX or PSDY", 1200, 0, 30)
    h_psd_chargeY_after = TH1F("h_psd_chargeY_after", "h_psd_chargeY PSDX and PSDY", 1200, 0, 30)
    h_psd_chargeX_after = TH1F("h_psd_chargeX_after", "h_psd_chargeX PSDX and PSDY", 1200, 0, 30)
    h_psd_resX = []
    h_psd_resY = []
    h_psd_resX_top = TH1F("h_psd_resX_top", "h_psd_resX_top", 400, -20, 20)
    h_psd_resY_top = TH1F("h_psd_resY_top", "h_psd_resY_top", 400, -20, 20)
    h_psd_resX_bot = TH1F("h_psd_resX_bot", "h_psd_resX_bot", 400, -20, 20)
    h_psd_resY_bot = TH1F("h_psd_resY_bot", "h_psd_resY_bot", 400, -20, 20)

    for iNpsd in xrange(0, 2):
        h_psd_resX.append(TH1F("h_psd_resX_%d" % iNpsd, "h_psd_resX_%d" % iNpsd, 2000, -1000, 1000))
        h_psd_resY.append(TH1F("h_psd_resY_%d" % iNpsd, "h_psd_resY_%d" % iNpsd, 2000, -1000, 1000))

    if is_mc:
        h_mc_stopZ = TH1F("h_mc_stopZ", "h_mc_stopZ", 2000, -500, 500)
        h_mc_stopX_vs_Y = TH2F("h_mc_stopX_vs_Y", "h_mc_stopX_vs_Y", 2000, -500, 500, 2000, -500, 500)
        h_stopZ_good = TH1F("h_stopZ_good", "h_stopZ_good", 2000, -500, 500)
        h_stopX_good = TH1F("h_stopX_good", "h_stopX_good", 1000, -500, 500)
        h_stopY_good = TH1F("h_stopY_good", "h_stopY_good", 1000, -500, 500)
        h_bgo_projection_good_Top = TH2F("h_bgo_projection_good_Top", "h_bgo_projection_good_Top", 2000, -500, 500, 2000, -500, 500)
        h_bgo_projection_good_Bot = TH2F("h_bgo_projection_good_Bot", "h_bgo_projection_good_Bot", 2000, -500, 500, 2000, -500, 500)
        h_bgo_deltaTheta_bgo = TH1F("h_bgo_deltaTheta_bgo", "h_bgo_deltaTheta_bgo", 200, -5, 5)
        h_bgo_deltaTheta_stk = TH1F("h_bgo_deltaTheta_stk", "h_bgo_deltaTheta_stk", 200, -5, 5)

        h_energy_good = []
        h_energy_bad = []
        h_energy_matrix = []
        h_energy_truth = {}
        h_energy_truth_weight = {}
        h_energy_weight = {}
        for ihisto_mc in range(0, ncuts):
            h_ene_cut_good = TH1F('h_energy_good_cut_' + str(ihisto_mc), 'h_energy_good_cut_' + str(ihisto_mc), nEDBINs_120bins, Edges_120bins)
            h_ene_cut_good.Sumw2()
            h_energy_good.append(h_ene_cut_good)

            h_ene_cut_bad = TH1F('h_energy_bad_cut_' + str(ihisto_mc), 'h_energy_bad_cut_' + str(ihisto_mc), nEDBINs_120bins, Edges_120bins)
            h_ene_cut_bad.Sumw2()
            h_energy_bad.append(h_ene_cut_bad)

            h_energy_matrix_histo = TH2F('h_energy_matrix_cut_' + str(ihisto_mc), 'h_energy_matrix_cut_' + str(ihisto_mc), nEDBINs_120bins, Edges_120bins, nEDBINs_120bins, Edges_120bins)
            h_energy_matrix.append(h_energy_matrix_histo)

            h_energy_truth[ihisto_mc] = {}
            h_energy_truth_weight[ihisto_mc] = {}
            h_energy_weight[ihisto_mc] = {}
            for bname, (nbins, edges) in energy_binnings.items():
                h_t = TH1F("h_energy_truth_{}_cut_{}".format(bname, ihisto_mc), "h_energy_truth_{}_cut_{}".format(bname, ihisto_mc), nbins, edges)
                h_t.Sumw2()
                h_energy_truth[ihisto_mc][bname] = h_t
                h_energy_truth_weight[ihisto_mc][bname] = {}
                h_energy_weight[ihisto_mc][bname] = {}
                for tag, gamma in spectral_indices.items():
                    h_tw = TH1F("h_energy_truth_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc), "h_energy_truth_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc), nbins, edges)
                    h_tw.Sumw2()
                    h_energy_truth_weight[ihisto_mc][bname][tag] = h_tw
                    h = TH1F("h_energy_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc), "h_energy_{}_weight_{}_cut_{}".format(bname, tag, ihisto_mc), nbins, edges)
                    h.Sumw2()
                    h_energy_weight[ihisto_mc][bname][tag] = h
    else:
        h_energy_good = None
        h_energy_bad = None
        h_energy_matrix = None
        h_energy_truth = None
        h_energy_truth_weight = None
        h_energy_weight = None

    h_energy_NoTracks = TH1F("h_energy_NoTracks", "h_energy_NoTracks", nEDBINs_120bins, Edges_120bins)
    h_energy_trigger_check = TH1F("h_energy_trigger_check", "h_energy_trigger_check", nEDBINs_120bins, Edges_120bins)
    h_energy = []
    for ihisto in range(0, ncuts):
        h_energy.append(TH1F('h_energy_cut_' + str(ihisto), 'h_energy_cut_' + str(ihisto), nEDBINs_120bins, Edges_120bins))

    h_stk_cluster_XvsY = []
    for iplane in range(0, 6):
        h_stk_cluster_XvsY.append(TH2F('stk_cluster_XvsY_plane_' + str(iplane), 'cluster X vs Y - plane ' + str(iplane), 1000, -500, 500, 1000, -500, 500))

    h_bgo_ene_lay = []
    for ilay in range(0, 14):
        h_bgo_ene_lay.append(TH1F('h_ngo_ene_lay_' + str(ilay), 'h_ngo_ene_lay_' + str(ilay), 1000, 0, 100))

    # -----------------------------
    # Event loop
    # -----------------------------
    nevents = dmpch.GetEntries()
    print ' Number of events: ', nevents

    ntrack_newsel = 0
    dmpch.SetOutputDir(root_dir, stream_name)
    if opts.skim:
        print ' skim option activeted'

    wrong_id = 0
    entry_ntupla = 0

    for iev in xrange(0, nevents):
        if (iev % 10000 == 0):
            print 'Processing event: ', iev

        if is_mc:
            DmpVSvc.gPsdECor.SetMCflag(1)

        pev = dmpch.GetDmpEvent(iev)
        metach.GetEntry(iev)

        # reset event-dependent scalar branches that may be absent in data/MC modes
        fBGO_EnergyG_SatCorr_ML_ions[0] = -999.
        fBGO_EnergyG_SatCorr_ML_ions_v3[0] = -999.
        if is_mc:
            fBGO_EnergyG_Quench[0] = -999.
            fBGO_EnergyG_QuenchSatCorr_ML_ions[0] = -999.
            fBGO_EnergyG_QuenchSatCorr_ML_ions_v3[0] = -999.

        good_event = False

        stopZ_parent = -1000
        stopX_parent = -1000
        stopY_parent = -1000
        startX_parent = -1000
        startY_parent = -1000
        startZ_parent = -1000

        # ------- ENERGY COMPUTATION AND SKIM CUTS 

        bgorec = pev.pEvtBgoRec()
        #bgoquenchrec = pev.pEvtBgoQuenchRec() if has_quenching else None
        stkclusters = pev.GetStkSiClusterCollection()
        stkladders = pev.GetStkLadderAdcCollection()

        initialize_or_rebind_apis(
            iev, is_mc, bgorec, stkclusters, stkladders,
            runmetadata, pevspace_api, bgo_api
        )

        pevspace_api.Predict(bgodirection=True, stkvertex=True, stktrack=True)
        pevspace_api.ObtainTrackHits()

        if is_mc:
            parent_particleID = int(pev.pEvtSimuPrimaries().pvpart_pdg / 10000.)
            if int(parent_particleID) != mc_pdg_particle_id:
                print "*************WRONG PARTICLE ID, particle ID ", pev.pEvtSimuPrimaries().pvpart_pdg
                wrong_id += 1
                continue

            for it in xrange(0, pev.NSimuTrajectory()):
                if pev.pSimuTrajectory(0).parentID == 0:
                    stopZ_parent = pev.pSimuTrajectory(0).stop_z
                    stopX_parent = pev.pSimuTrajectory(0).stop_x
                    stopY_parent = pev.pSimuTrajectory(0).stop_y
                    startX_parent = pev.pSimuTrajectory(0).start_x
                    startY_parent = pev.pSimuTrajectory(0).start_y
                    startZ_parent = pev.pSimuTrajectory(0).start_z
                else:
                    print " ERROR NO PARENT!!! "
                    break

            h_mc_stopX_vs_Y.Fill(stopX_parent, stopY_parent)
            h_mc_stopZ.Fill(stopZ_parent)
            if stopZ_parent < 0.:
                def_good_truth_stopXY = abs(stopX_parent) < 400. and abs(stopY_parent) < 400.
            else:
                def_good_truth_stopXY = abs(stopX_parent) < 300. and abs(stopY_parent) < 300.

            def_good_truth_stopZ = (stopZ_parent <= 350. and stopZ_parent >= -280.)
            if startZ_parent != stopZ_parent:
                interceptXZ_truth = ((startZ_parent * stopX_parent) - (stopZ_parent * startX_parent)) / (startZ_parent - stopZ_parent)
                interceptYZ_truth = ((startZ_parent * stopY_parent) - (stopZ_parent * startY_parent)) / (startZ_parent - stopZ_parent)
                slope_X = (startX_parent - stopX_parent) / (startZ_parent - stopZ_parent)
                slope_Y = (startY_parent - stopY_parent) / (startZ_parent - stopZ_parent)
                x_projection_truth_onBGO_Top = interceptXZ_truth + BGOzTop * slope_X
                y_projection_truth_onBGO_Top = interceptYZ_truth + BGOzTop * slope_Y
                x_projection_truth_onBGO_Bot = interceptXZ_truth + BGOzBot * slope_X
                y_projection_truth_onBGO_Bot = interceptYZ_truth + BGOzBot * slope_Y
                def_good_truth_bgoAcc = abs(x_projection_truth_onBGO_Top) < bgo_acceptance_range and abs(y_projection_truth_onBGO_Top) < bgo_acceptance_range and abs(x_projection_truth_onBGO_Bot) < bgo_acceptance_range and abs(y_projection_truth_onBGO_Bot) < bgo_acceptance_range
            else:
                def_good_truth_bgoAcc = False

            if def_good_truth_stopXY and def_good_truth_stopZ and def_good_truth_bgoAcc:
                good_event = True
                h_stopZ_good.Fill(stopZ_parent)
                h_stopX_good.Fill(stopX_parent)
                h_stopY_good.Fill(stopY_parent)
                h_bgo_projection_good_Top.Fill(x_projection_truth_onBGO_Top, y_projection_truth_onBGO_Top)
                h_bgo_projection_good_Bot.Fill(x_projection_truth_onBGO_Bot, y_projection_truth_onBGO_Bot)
                if abs(stopX_parent) > 400. or abs(stopY_parent) > 400.:
                    print iev, stopX_parent, stopY_parent

        (etot,
         etot_satcorr_ions,
         etot_satcorr_ions_v3,
         etot_quench,
         etot_quench_satcorr_ions,
         etot_quench_satcorr_ions_v3) = compute_bgo_energies(
            pev, bgo_api, is_mc, has_quenching, mc_particle_Z
        )

        etot_truth = -9.
        if is_mc:
            etot_truth = pev.pEvtSimuPrimaries().pvpart_ekin / 1000.

        if is_data:
            longitude = pev.pEvtAttitude().lon_geo
            latitude = pev.pEvtAttitude().lat_geo
        else:
            longitude = 0
            latitude = 0

        flong[0] = longitude
        flat[0] = latitude

        fill_energy_cut(0, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        # cut 1
        if etot < e_min:
            continue

        fill_energy_cut(1, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        x_projection_fromBGO_to_BGO = pevspace_api.GetDirectionBGOInteceptX() + BGOzBot * pevspace_api.GetDirectionBGOSlopeX()
        y_projection_fromBGO_to_BGO = pevspace_api.GetDirectionBGOInteceptY() + BGOzBot * pevspace_api.GetDirectionBGOSlopeY()
        x_projection_fromBGO_to_BGOTop = pevspace_api.GetDirectionBGOInteceptX() + BGOzTop * pevspace_api.GetDirectionBGOSlopeX()
        y_projection_fromBGO_to_BGOTop = pevspace_api.GetDirectionBGOInteceptY() + BGOzTop * pevspace_api.GetDirectionBGOSlopeY()

        h_bgo_projection_before.Fill(x_projection_fromBGO_to_BGO, y_projection_fromBGO_to_BGO)
        h_bgo_projectionTop_before.Fill(x_projection_fromBGO_to_BGOTop, y_projection_fromBGO_to_BGOTop)

        bgo_acceptance = abs(x_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGO) < bgo_acceptance_range and abs(x_projection_fromBGO_to_BGOTop) < bgo_acceptance_range and abs(y_projection_fromBGO_to_BGOTop) < bgo_acceptance_range

        # cut 2
        if ((pevspace_api.GetDirectionBGOSlopeX() == 0 and pevspace_api.GetDirectionBGOInteceptX() == 0) or
                (pevspace_api.GetDirectionBGOSlopeY() == 0 and pevspace_api.GetDirectionBGOInteceptY() == 0)):
            continue

        fBGO_interceptX_ML[0] = pevspace_api.GetDirectionBGOInteceptX()
        fBGO_interceptY_ML[0] = pevspace_api.GetDirectionBGOInteceptY()
        fBGO_slopeX_ML[0] = pevspace_api.GetDirectionBGOSlopeX()
        fBGO_slopeY_ML[0] = pevspace_api.GetDirectionBGOSlopeY()

        # cut 3
        if not bgo_acceptance:
            continue

        fill_energy_cut(2, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        h_bgo_projection_after.Fill(x_projection_fromBGO_to_BGO, y_projection_fromBGO_to_BGO)
        h_bgo_projectionTop_after.Fill(x_projection_fromBGO_to_BGOTop, y_projection_fromBGO_to_BGOTop)

        v_bgolayer = np.array([pev.pEvtBgoRec().GetELayer(ibgo) for ibgo in range(14)])
        frac_max_lay = np.max(v_bgolayer) / 1000. / etot

        h_bgo_EratioLay_tot_before.Fill(frac_max_lay)
        if good_event:
            h_bgo_EratioLay_goo_before.Fill(frac_max_lay)
        else:
            h_bgo_EratioLay_bad_before.Fill(frac_max_lay)

        # cut 4
        if frac_max_lay > 0.35:
            continue
        h_bgo_EratioLay_tot_after.Fill(frac_max_lay)

        fill_energy_cut(3, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        pos_max_lateral = False
        lateral_pos_vec = []
        for ilay in xrange(1, 4):
            v_bgolayer_bars = np.array([pev.pEvtBgoRec().GetEdepPos(ilay, ibar) for ibar in xrange(0, 23)])
            lateral_pos = np.argmax(v_bgolayer_bars)
            ene_lay = np.max(v_bgolayer_bars)
            if ene_lay <= 0.:
                lateral_pos = -1
            h_bgo_lateralpos_before.Fill(lateral_pos)
            lateral_pos_vec.append(lateral_pos)

        barNumberMaxEBarLayer = [-1, -1, -1]
        MaxEBarLayer = [0, 0, 0]
        for ihit in xrange(0, pev.NEvtBgoHits()):
            hitE = pev.pEvtBgoHits().fEnergy[ihit]
            layID = pev.pEvtBgoHits().GetLayerID(ihit)
            if layID == 1 or layID == 2 or layID == 3:
                if hitE > MaxEBarLayer[layID - 1]:
                    barNumberMaxEBarLayer[layID - 1] = (pev.pEvtBgoHits().fGlobalBarID[ihit] >> 6) & 0x1f
                    MaxEBarLayer[layID - 1] = hitE

        for jj in xrange(0, len(barNumberMaxEBarLayer)):
            if barNumberMaxEBarLayer[jj] <= 0 or barNumberMaxEBarLayer[jj] == 21:
                pos_max_lateral = True

        # cut 5
        if pos_max_lateral:
            continue

        for ilay in xrange(0, 3):
            h_bgo_lateralpos_after.Fill(lateral_pos_vec[ilay])

        fill_energy_cut(4, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        if opts.skim:
            dmpch.SaveCurrentEvent(stream_name)

        # ------- SAA EXCLUSION (FLIGHT DATA)

        # cut 6
        if is_data:
            inSAA = pFilter.IsInSAA(pev.pEvtHeader().GetSecond())
            if inSAA:
                continue
            h_energy[5].Fill(etot)
            h_terrestrial_lat_vs_long.Fill(longitude, latitude)
        else:
            fill_energy_cut(5, etot, etot_truth, good_event, is_mc,
                            h_energy, h_energy_matrix,
                            h_energy_truth, h_energy_weight, h_energy_truth_weight,
                            energy_binnings, spectral_indices,
                            h_energy_good, h_energy_bad)

        if pev.pEvtHeader().GeneratedTrigger(3) == True or pev.pEvtHeader().GeneratedTrigger(0) == True:
            h_energy_trigger_check.Fill(etot)

        h_energy[6].Fill(etot)

        # ------- TRACKING 

        # cut 7 --> ntrack: save a flag to eventually add this cut! (even with ML tracking)
        fSTK_ntrack[0] = pev.NStkKalmanTrack()
        fSTK_nclusters[0] = pev.NStkSiCluster()

        fSTK_vertexPrediction[0] = pevspace_api.GetVertexPrediction()

        theta_bgo = np.arctan(np.sqrt(pevspace_api.GetDirectionBGOSlopeX() ** 2 + pevspace_api.GetDirectionBGOSlopeY() ** 2))
        theta_stk = np.arctan(np.sqrt(pevspace_api.GetDirectionSTKSlopeX() ** 2 + pevspace_api.GetDirectionSTKSlopeY() ** 2))
        track_correction = np.cos(theta_stk)
        theta_bgo_deg = theta_bgo * 180 / np.pi
        theta_stk_deg = theta_stk * 180 / np.pi

        phi_bgo = math.atan(pevspace_api.GetDirectionBGOSlopeY() / pevspace_api.GetDirectionBGOSlopeX())
        phi_stk = math.atan(pevspace_api.GetDirectionSTKSlopeY() / pevspace_api.GetDirectionSTKSlopeX())
        phi_bgo_deg = phi_bgo * 180 / np.pi
        phi_stk_deg = phi_stk * 180 / np.pi

        fBGO_slopeX_ML[0] = pevspace_api.GetDirectionBGOSlopeX()
        fBGO_slopeY_ML[0] = pevspace_api.GetDirectionBGOSlopeY()
        fBGO_slopeXZ_analy[0] = pev.pEvtBgoRec().GetSlopeXZ()
        fBGO_slopeYZ_analy[0] = pev.pEvtBgoRec().GetSlopeYZ()
        fBGO_theta_analy[0] = np.arctan(np.sqrt(pev.pEvtBgoRec().GetSlopeXZ() ** 2 + pev.pEvtBgoRec().GetSlopeYZ() ** 2)) * 180 / np.pi
        if pev.pEvtBgoRec().GetSlopeXZ() > 0:
            fBGO_phi_analy[0] = math.atan(pev.pEvtBgoRec().GetSlopeYZ() / pev.pEvtBgoRec().GetSlopeXZ()) * 180 / np.pi
        fBGO_interceptXZ_analy[0] = pev.pEvtBgoRec().GetInterceptXZ()
        fBGO_interceptYZ_analy[0] = pev.pEvtBgoRec().GetInterceptYZ()

        fSTK_slopeX_ML[0] = pevspace_api.GetDirectionSTKSlopeX()
        fSTK_slopeY_ML[0] = pevspace_api.GetDirectionSTKSlopeY()
        fSTK_interceptX_ML[0] = pevspace_api.GetTrackHitImpactX()
        fSTK_interceptY_ML[0] = pevspace_api.GetTrackHitImpactY()
        fSTK_interceptZ_ML[0] = pevspace_api.GetTrackHitImpactZ()

        hitClusterSum = [0, 0]
        etaCorr = DmpStkClusterEtaCorrection()
        cluster_chargeX = array('d', 6 * [-999.])
        cluster_chargeY = array('d', 6 * [-999.])
        cluster_chargeX_etacorr = array('d', 6 * [-999.])
        cluster_chargeY_etacorr = array('d', 6 * [-999.])

        for j in range(12):
            hitclustersignal = pevspace_api.GetTrackHitSignal(j)
            hitclusterid = pevspace_api.GetTrackHitID(j)
            if hitclusterid >= 0:
                hitClusterSum[j % 2] += 1
                cluster = stkclusters.ConstructedAt(hitclusterid)
                if j % 2 == 0:
                    cluster_chargeX[j / 2] = hitclustersignal
                    cluster_chargeX_etacorr[j / 2] = etaCorr.getStkEtaCorrEnergy(cluster, track_correction) * track_correction
                else:
                    cluster_chargeY[j / 2] = hitclustersignal
                    cluster_chargeY_etacorr[j / 2] = etaCorr.getStkEtaCorrEnergy(cluster, track_correction) * track_correction

        # cut 9 
        if hitClusterSum[0] == 0 or hitClusterSum[1] == 0:
            continue

        fill_energy_cut(7, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        h_energy[8].Fill(etot)
        deltaTheta_rec_sel = theta_bgo_deg - theta_stk_deg
        h_bgo_stk_deltaTheta_tracksel.Fill(deltaTheta_rec_sel)

        theta_truth = -9.
        phi_truth = -9.
        if is_mc:
            theta_truth = math.atan(np.sqrt((pev.pEvtSimuPrimaries().pvpart_px / pev.pEvtSimuPrimaries().pvpart_pz * pev.pEvtSimuPrimaries().pvpart_px / pev.pEvtSimuPrimaries().pvpart_pz) + (pev.pEvtSimuPrimaries().pvpart_py / pev.pEvtSimuPrimaries().pvpart_pz * pev.pEvtSimuPrimaries().pvpart_py / pev.pEvtSimuPrimaries().pvpart_pz))) * 180. / math.pi
            phi_truth = math.atan(pev.pEvtSimuPrimaries().pvpart_py / pev.pEvtSimuPrimaries().pvpart_px) * 180. / math.pi
            h_bgo_deltaTheta_bgo.Fill(theta_truth - theta_bgo)
            h_bgo_deltaTheta_stk.Fill(theta_truth - theta_stk)

        h_stk_clusterX.Fill(cluster_chargeX[0])
        h_stk_clusterY.Fill(cluster_chargeY[0])

        track_slope = TVector3(pevspace_api.GetDirectionSTKSlopeX(), pevspace_api.GetDirectionSTKSlopeY(), 1)
        track_impact_point = TVector3(pevspace_api.GetTrackHitImpactX(), pevspace_api.GetTrackHitImpactY(), pevspace_api.GetTrackHitImpactZ())

        psd_YZ_top = -324.7
        psd_XZ_top = -298.5
        stk_to_psd_topY = (track_slope.Y() * (psd_YZ_top - track_impact_point.Z()) + track_impact_point.Y())
        stk_to_psd_topX = (track_slope.X() * (psd_XZ_top - track_impact_point.Z()) + track_impact_point.X())

        fSTKtrack_to_PSD_topY[0] = stk_to_psd_topY
        fSTKtrack_to_PSD_topX[0] = stk_to_psd_topX

        h_energy[9].Fill(etot)

        psd_vec_chargeX = [[] for _ in range(2)]
        psd_vec_chargeX0 = [[] for _ in range(2)]
        psd_vec_chargeX1 = [[] for _ in range(2)]
        psd_vec_gidX = [[] for _ in range(2)]
        psd_vec_pathlengthX = [[] for _ in range(2)]
        psd_vec_positionX = [[] for _ in range(2)]
        psd_bar_id_x = [[] for _ in range(2)]

        psd_vec_chargeY = [[] for _ in range(2)]
        psd_vec_chargeY0 = [[] for _ in range(2)]
        psd_vec_chargeY1 = [[] for _ in range(2)]
        psd_vec_gidY = [[] for _ in range(2)]
        psd_vec_pathlengthY = [[] for _ in range(2)]
        psd_vec_positionY = [[] for _ in range(2)]
        psd_bar_id_y = [[] for _ in range(2)]

        for ipsd in xrange(0, pev.NEvtPsdHits()):
            if pev.pEvtPsdHits().IsHitMeasuringX(ipsd):
                lengthX = [-99999., -99999.]
                array_lengthX = array('d', lengthX)
                if pev.pEvtPsdHits().GetHitZ(ipsd) < -298.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -298.51:
                    npsdX = 0
                if pev.pEvtPsdHits().GetHitZ(ipsd) < -284.49 and pev.pEvtPsdHits().GetHitZ(ipsd) > -284.51:
                    npsdX = 1

                if is_mc:
                    crossingX = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd], track_slope, track_impact_point, array_lengthX)
                else:
                    crossingX = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd], track_slope, track_impact_point, array_lengthX)

                if crossingX == True:
                    psd_vec_chargeX[npsdX].append(pev.pEvtPsdHits().fEnergy[ipsd])
                    psd_vec_chargeX0[npsdX].append(pev.pEvtPsdHits().fES0[ipsd])
                    psd_vec_chargeX1[npsdX].append(pev.pEvtPsdHits().fES1[ipsd])
                    psd_vec_gidX[npsdX].append(pev.pEvtPsdHits().fGlobalBarID[ipsd])
                    psd_bar_id_x[npsdX].append(pev.pEvtPsdHits().GetHittedBarNumber())
                    psd_vec_pathlengthX[npsdX].append(array_lengthX[1])
                    psd_vec_positionX[npsdX].append(pev.pEvtPsdHits().GetHitX(ipsd))

            elif pev.pEvtPsdHits().IsHitMeasuringY(ipsd):
                lengthY = [-99999., -99999.]
                array_lengthY = array('d', lengthY)
                if pev.pEvtPsdHits().GetHitZ(ipsd) < -324.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -324.71:
                    npsdY = 0
                if pev.pEvtPsdHits().GetHitZ(ipsd) < -310.69 and pev.pEvtPsdHits().GetHitZ(ipsd) > -310.71:
                    npsdY = 1

                if is_mc:
                    crossingY = DmpVSvc.gPsdECor.GetPathLPMC(pev.pEvtPsdHits().fGlobalBarID[ipsd], track_slope, track_impact_point, array_lengthY)
                else:
                    crossingY = DmpVSvc.gPsdECor.GetPathLengthPosition(pev.pEvtPsdHits().fGlobalBarID[ipsd], track_slope, track_impact_point, array_lengthY)

                if crossingY == True:
                    psd_vec_chargeY[npsdY].append(pev.pEvtPsdHits().fEnergy[ipsd])
                    psd_vec_chargeY0[npsdY].append(pev.pEvtPsdHits().fES0[ipsd])
                    psd_vec_chargeY1[npsdY].append(pev.pEvtPsdHits().fES1[ipsd])
                    psd_vec_gidY[npsdY].append(pev.pEvtPsdHits().fGlobalBarID[ipsd])
                    psd_bar_id_y[npsdY].append(pev.pEvtPsdHits().GetHittedBarNumber())
                    psd_vec_pathlengthY[npsdY].append(array_lengthY[1])
                    psd_vec_positionY[npsdY].append(pev.pEvtPsdHits().GetHitY(ipsd))

        check_psd = False
        for ipsd in xrange(0, 2):
            if len(psd_vec_chargeY[ipsd]) > 1 or len(psd_vec_chargeX[ipsd]) > 1:
                check_psd = True
        if check_psd:
            print '******* ERROR!!!!!!!!!!!!! *****'
            print 'Event with more than 1 psd bar per layer', iev
            print psd_vec_chargeX, psd_vec_chargeY
            print psd_vec_pathlengthX, psd_vec_pathlengthY
            print psd_vec_positionX, psd_vec_positionY

        psdchargeX = [-999., -999.]
        psdchargeY = [-999., -999.]
        psdchargeX_corr = [-999., -999.]
        psdchargeX0_corr = [-999., -999.]
        psdchargeX1_corr = [-999., -999.]
        psdchargeX01_corr = [-999., -999.]
        psdchargeY_corr = [-999., -999.]
        psdchargeY0_corr = [-999., -999.]
        psdchargeY1_corr = [-999., -999.]
        psdchargeY01_corr = [-999., -999.]
        psdvec = []
        psdvec_pathlength = []
        psdchargeX_proj = [-999., -999.]
        psdchargeY_proj = [-999., -999.]
        psdX_pathlength = [-999., -999.]
        psdY_pathlength = [-999., -999.]
        psdX_position = [-999., -999.]
        psdY_position = [-999., -999.]

        for ipsd in xrange(0, 2):
            if len(psd_vec_chargeY[ipsd]) > 0:
                pos_max_len = np.argmax(psd_vec_pathlengthY[ipsd])
                lengthY = [-99999., -99999.]
                array_lengthY = array('d', lengthY)
                if is_mc:
                    test_pos = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidY[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthY)
                else:
                    test_pos = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidY[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthY)

                if test_pos == True:
                    PsdEC_tmpY = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidY[ipsd][pos_max_len], array_lengthY[0] / 10.)
                else:
                    print "****** ERROR! it should be true!!!"
                    PsdEC_tmpY = -1.

                AttCorY0 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidY[ipsd][pos_max_len], 0, array_lengthY[0])
                AttCorY1 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidY[ipsd][pos_max_len], 1, array_lengthY[0])
                lay = DmpPsdBase.GetLayerID(psd_vec_gidY[ipsd][pos_max_len])
                bar = DmpPsdBase.GetStripID(psd_vec_gidY[ipsd][pos_max_len])

                psdchargeY[ipsd] = psd_vec_chargeY[ipsd][pos_max_len]
                psdchargeY_corr[ipsd] = psd_vec_chargeY[ipsd][pos_max_len] * PsdEC_tmpY
                if len(psd_vec_chargeY0[ipsd]) > 0:
                    psdchargeY0_corr[ipsd] = np.sqrt(psd_vec_chargeY0[ipsd][pos_max_len] / AttCorY0 * 10. / array_lengthY[1])
                    if is_data:
                        psdchargeY0_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 0, psdchargeY0_corr[ipsd])
                else:
                    psdchargeY0_corr[ipsd] = -1.

                if len(psd_vec_chargeY1[ipsd]) > 0:
                    psdchargeY1_corr[ipsd] = np.sqrt(psd_vec_chargeY1[ipsd][pos_max_len] / AttCorY1 * 10. / array_lengthY[1])
                    if is_data:
                        psdchargeY1_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 1, psdchargeY1_corr[ipsd])
                else:
                    psdchargeY1_corr[ipsd] = -1.

                psdchargeY01_corr[ipsd] = (psdchargeY0_corr[ipsd] + psdchargeY1_corr[ipsd]) / 2.0
                if is_data:
                    psdchargeY01_corr[ipsd] = DmpVSvc.gPsdECor.GetPosCorrectedCharge(lay, bar, array_lengthY[0], psdchargeY01_corr[ipsd])
                psdchargeY_proj[ipsd] = array_lengthY[0]
                psdY_pathlength[ipsd] = array_lengthY[1]
                psdY_position[ipsd] = psd_vec_positionY[ipsd][pos_max_len]

                h_psd_GlobalID.Fill(psd_vec_gidY[ipsd][pos_max_len])
                if psdchargeY[ipsd] > 0.:
                    h_psd_chargeY.Fill(np.sqrt(psdchargeY[ipsd] / 2.))
                if psdchargeY_corr[ipsd] > 0:
                    h_psd_chargeY_corr.Fill(np.sqrt(psdchargeY_corr[ipsd] / 2.))
                if psd_vec_gidY[ipsd][pos_max_len] == 512 and array_lengthY[0] != -99999.:
                    h_psd_corr_ID512.Fill(array_lengthY[0], PsdEC_tmpY)

            if len(psd_vec_chargeX[ipsd]) > 0:
                pos_max_len = np.argmax(psd_vec_pathlengthX[ipsd])
                lengthX = [-99999., -99999.]
                array_lengthX = array('d', lengthX)
                if is_mc:
                    test_posX = DmpVSvc.gPsdECor.GetPathLPMC(psd_vec_gidX[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthX)
                else:
                    test_posX = DmpVSvc.gPsdECor.GetPathLengthPosition(psd_vec_gidX[ipsd][pos_max_len], track_slope, track_impact_point, array_lengthX)

                if test_posX == True:
                    PsdEC_tmpX = DmpVSvc.gPsdECor.GetPsdECor(psd_vec_gidX[ipsd][pos_max_len], array_lengthX[0] / 10.)
                else:
                    PsdEC_tmpX = -1.

                AttCorX0 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidX[ipsd][pos_max_len], 0, array_lengthX[0])
                AttCorX1 = DmpVSvc.gPsdECor.GetPsdMipAttESp3(psd_vec_gidX[ipsd][pos_max_len], 1, array_lengthX[0])
                lay = DmpPsdBase.GetLayerID(psd_vec_gidX[ipsd][pos_max_len])
                bar = DmpPsdBase.GetStripID(psd_vec_gidX[ipsd][pos_max_len])

                psdchargeX[ipsd] = psd_vec_chargeX[ipsd][pos_max_len]
                psdchargeX_corr[ipsd] = psd_vec_chargeX[ipsd][pos_max_len] * PsdEC_tmpX
                if len(psd_vec_chargeX0[ipsd]) > 0:
                    psdchargeX0_corr[ipsd] = np.sqrt(psd_vec_chargeX0[ipsd][pos_max_len] / AttCorX0 * 10. / array_lengthX[1])
                    if is_data:
                        psdchargeX0_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 0, psdchargeX0_corr[ipsd])
                else:
                    psdchargeX0_corr[ipsd] = -1.

                if len(psd_vec_chargeX1[ipsd]) > 0:
                    psdchargeX1_corr[ipsd] = np.sqrt(psd_vec_chargeX1[ipsd][pos_max_len] / AttCorX1 * 10. / array_lengthX[1])
                    if is_data:
                        psdchargeX1_corr[ipsd] = DmpVSvc.gPsdECor.GetQuenchEqalCalFunTwoEnd(lay, bar, 1, psdchargeX1_corr[ipsd])
                else:
                    psdchargeX1_corr[ipsd] = -1.

                psdchargeX01_corr[ipsd] = (psdchargeX0_corr[ipsd] + psdchargeX1_corr[ipsd]) / 2.0
                if is_data:
                    psdchargeX01_corr[ipsd] = DmpVSvc.gPsdECor.GetPosCorrectedCharge(lay, bar, array_lengthX[0], psdchargeX01_corr[ipsd])
                psdchargeX_proj[ipsd] = array_lengthX[0]
                psdX_pathlength[ipsd] = array_lengthX[1]
                psdX_position[ipsd] = psd_vec_positionX[ipsd][pos_max_len]

                h_psd_GlobalID.Fill(psd_vec_gidX[ipsd][pos_max_len])
                if psdchargeX[ipsd] > 0.:
                    h_psd_chargeX.Fill(np.sqrt(psdchargeX[ipsd] / 2.))
                if psdchargeX_corr[ipsd] > 0.:
                    h_psd_chargeX_corr.Fill(np.sqrt(psdchargeX_corr[ipsd] / 2.))
                if psd_vec_gidX[ipsd][pos_max_len] == 5504 and array_lengthX[0] != -99999.:
                    h_psd_corr_ID5504.Fill(array_lengthX[0], PsdEC_tmpX)

        for ii in xrange(2):
            if psdchargeY01_corr[ii] > 0.:
                psdvec.append(psdchargeY01_corr[ii])
                psdvec_pathlength.append(psdY_pathlength[ii])
        for ii in xrange(2):
            if psdchargeX01_corr[ii] > 0.:
                psdvec.append(psdchargeX01_corr[ii])
                psdvec_pathlength.append(psdX_pathlength[ii])

        h_bgo_theta.Fill(theta_bgo)

        sum_len_psd_vec_chargeX = len(psd_vec_chargeX[0]) + len(psd_vec_chargeX[1])
        sum_len_psd_vec_chargeY = len(psd_vec_chargeY[0]) + len(psd_vec_chargeY[1])

        for ipsd in xrange(0, 2):
            if psdchargeY_corr[ipsd] > 0.:
                h_psd_chargeY_before.Fill(np.sqrt(psdchargeY_corr[ipsd] / 2.))
            if psdchargeX_corr[ipsd] > 0.:
                h_psd_chargeX_before.Fill(np.sqrt(psdchargeX_corr[ipsd] / 2.))

        # cut 13 
        if sum_len_psd_vec_chargeY == 0:
            continue
        if sum_len_psd_vec_chargeX == 0:
            continue

        for ipsd in xrange(0, 2):
            if psdchargeY_corr[ipsd] > 0.:
                h_psd_chargeY_after.Fill(np.sqrt(psdchargeY_corr[ipsd] / 2.))
            if psdchargeX_corr[ipsd] > 0.:
                h_psd_chargeX_after.Fill(np.sqrt(psdchargeX_corr[ipsd] / 2.))

        fill_energy_cut(10, etot, etot_truth, good_event, is_mc,
                        h_energy, h_energy_matrix,
                        h_energy_truth, h_energy_weight, h_energy_truth_weight,
                        energy_binnings, spectral_indices,
                        h_energy_good, h_energy_bad)

        fired_layer = np.count_nonzero(v_bgolayer)
        h_bgo_firedLayer_before.Fill(fired_layer)
        h_energy[11].Fill(etot)

        Ene4Layers = sum(v_bgolayer[:4]) / 1000 / etot
        Ene4LayersLast = sum(v_bgolayer[-4:]) / 1000 / etot
        Ene2Layers = sum(v_bgolayer[:2]) / 1000 / etot
        Ene2LayersLast = sum(v_bgolayer[-2:]) / 1000 / etot

        vec_psdX_max = []
        vec_psdY_max = []
        for ipsd in xrange(0, pev.NEvtPsdHits()):
            if pev.pEvtPsdHits().GetLayerID(ipsd) == 0:
                vec_psdY_max.append(pev.pEvtPsdHits().fEnergy[ipsd])
            elif pev.pEvtPsdHits().GetLayerID(ipsd) == 1:
                vec_psdX_max.append(pev.pEvtPsdHits().fEnergy[ipsd])

        zpsd = [-324.7, -310.7, -298.5, -284.5]

        fpsdchargeX0[0] = psdchargeX01_corr[0]
        fpsdchargeX1[0] = psdchargeX01_corr[1]
        fpsdchargeY0[0] = psdchargeY01_corr[0]
        fpsdchargeY1[0] = psdchargeY01_corr[1]

        globalcharge = -999.
        globalaverage = -999.
        if len(psdvec) == 0:
            print 'There is no charge value!'
        elif len(psdvec) == 1:
            globalcharge = psdvec[0]
            globalaverage = psdvec[0]
        elif len(psdvec) > 1 and len(psdvec) < 5:
            PSDdiff = True
            for iQ in xrange(1, len(psdvec)):
                PSDdiff = PSDdiff and (abs(psdvec[iQ - 1] - psdvec[iQ]) < 1)
                if (iQ == 1) and (not PSDdiff):
                    psd_sum = psdvec[iQ - 1]
                    psd_sum_weight = psdvec[iQ - 1] * psdvec_pathlength[iQ - 1]
                    psd_nlay = 1.
                    psd_nlayweight = psdvec_pathlength[iQ - 1]
                elif (iQ == 1) and PSDdiff:
                    psd_sum = psdvec[iQ - 1] + psdvec[iQ]
                    psd_sum_weight = (psdvec[iQ - 1] * psdvec_pathlength[iQ - 1]) + (psdvec[iQ] * psdvec_pathlength[iQ])
                    psd_nlay = 2.
                    psd_nlayweight = (psdvec_pathlength[iQ - 1] + psdvec_pathlength[iQ])
                elif (iQ != 1) and PSDdiff:
                    psd_sum = psd_sum + psdvec[iQ]
                    psd_sum_weight = psd_sum_weight + (psdvec[iQ] * psdvec_pathlength[iQ])
                    psd_nlay = psd_nlay + 1.
                    psd_nlayweight = psd_nlayweight + psdvec_pathlength[iQ]
                else:
                    break
            globalcharge = psd_sum / psd_nlay
            globalaverage = psd_sum_weight / psd_nlayweight
        else:
            print 'Houston, we have a problem!'

        fPSD_Global_Charge[0] = globalcharge
        fPSD_Global_Charge_PathAverage[0] = globalaverage

        den_average_charge = (
            (((TMath.Sign(1., psdchargeY01_corr[0]) + 1.) / 2.) * psdY_pathlength[0]) +
            (((TMath.Sign(1., psdchargeY01_corr[1]) + 1.) / 2.) * psdY_pathlength[1]) +
            (((TMath.Sign(1., psdchargeX01_corr[0]) + 1.) / 2.) * psdX_pathlength[0]) +
            (((TMath.Sign(1., psdchargeX01_corr[1]) + 1.) / 2.) * psdX_pathlength[1])
        )

        if den_average_charge > 0:
            average_charge = (
                ((((TMath.Sign(1., psdchargeY01_corr[0]) + 1.) / 2.) * psdchargeY01_corr[0] * psdY_pathlength[0]) +
                 (((TMath.Sign(1., psdchargeY01_corr[1]) + 1.) / 2.) * psdchargeY01_corr[1] * psdY_pathlength[1]) +
                 (((TMath.Sign(1., psdchargeX01_corr[0]) + 1.) / 2.) * psdchargeX01_corr[0] * psdX_pathlength[0]) +
                 (((TMath.Sign(1., psdchargeX01_corr[1]) + 1.) / 2.) * psdchargeX01_corr[1] * psdX_pathlength[1]))
                / den_average_charge
            )
        else:
            average_charge = -999.

        fPSD_PathWeighted_Charge[0] = average_charge

        for ipsd in xrange(0, nlayer_psd):
            fPSD_EnergyY[ipsd] = psdchargeY[ipsd]
            fPSD_EnergyX[ipsd] = psdchargeX[ipsd]
            if ipsd < 2:
                fPSD_cpsdstk[ipsd] = track_impact_point.y() + track_slope.y() * (zpsd[ipsd] - track_impact_point.z())
            if ipsd > 1:
                fPSD_cpsdstk[ipsd] = track_impact_point.x() + track_slope.x() * (zpsd[ipsd] - track_impact_point.z())
            fPSD_EnergyY_corr[ipsd] = psdchargeY_corr[ipsd]
            fPSD_EnergyX_corr[ipsd] = psdchargeX_corr[ipsd]
            fPSD_psdX_projY[ipsd] = psdchargeX_proj[ipsd]
            fPSD_psdY_projX[ipsd] = psdchargeY_proj[ipsd]
            fPSD_psdX_pathlength[ipsd] = psdX_pathlength[ipsd]
            fPSD_psdY_pathlength[ipsd] = psdY_pathlength[ipsd]
            fPSD_psdX_position[ipsd] = psdX_position[ipsd]
            fPSD_psdY_position[ipsd] = psdY_position[ipsd]

        fPSD_nhits[0] = pev.NEvtPsdHits()
        fPSD_psdX_Emax_layer[0] = np.max(vec_psdX_max) if len(vec_psdX_max) else -999.
        fPSD_psdY_Emax_layer[0] = np.max(vec_psdY_max) if len(vec_psdY_max) else -999.

        fPSD_EmaxY[0] = 0.0
        fPSD_EmaxY[1] = 0.0
        fPSD_EmaxX[0] = 0.0
        fPSD_EmaxX[1] = 0.0
        for ipsd in xrange(0, pev.NEvtPsdHits()):
            if pev.pEvtPsdHits().GetHitZ(ipsd) == -324.7:
                if fPSD_EmaxY[0] < pev.pEvtPsdHits().fEnergy[ipsd]:
                    fPSD_EmaxY[0] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[0] = pev.pEvtPsdHits().GetHitY(ipsd)
            if pev.pEvtPsdHits().GetHitZ(ipsd) == -310.7:
                if fPSD_EmaxY[1] < pev.pEvtPsdHits().fEnergy[ipsd]:
                    fPSD_EmaxY[1] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[1] = pev.pEvtPsdHits().GetHitY(ipsd)
            if pev.pEvtPsdHits().GetHitZ(ipsd) == -298.5:
                if fPSD_EmaxX[0] < pev.pEvtPsdHits().fEnergy[ipsd]:
                    fPSD_EmaxX[0] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[2] = pev.pEvtPsdHits().GetHitX(ipsd) + 0.4
            if pev.pEvtPsdHits().GetHitZ(ipsd) == -284.5:
                if fPSD_EmaxX[1] < pev.pEvtPsdHits().fEnergy[ipsd]:
                    fPSD_EmaxX[1] = pev.pEvtPsdHits().fEnergy[ipsd]
                    fPSD_cpsdmax[3] = pev.pEvtPsdHits().GetHitX(ipsd) + 0.6

        den_psd_y = (
            (((TMath.Sign(1, fPSD_psdY_pathlength[0]) + 1) / 2) * fPSD_psdY_pathlength[0]) +
            (((TMath.Sign(1, fPSD_psdY_pathlength[1]) + 1) / 2) * fPSD_psdY_pathlength[1])
        )

        if den_psd_y > 0:
            fPSD_Etrack01[0] = (
                ((((TMath.Sign(1, fPSD_EnergyY_corr[0]) + 1) / 2) * fPSD_EnergyY_corr[0]) +
                 (((TMath.Sign(1, fPSD_EnergyY_corr[1]) + 1) / 2) * fPSD_EnergyY_corr[1]))
                * 10. / den_psd_y
            )
        else:
            fPSD_Etrack01[0] = -999.

        den_psd_x = (
            (((TMath.Sign(1, fPSD_psdX_pathlength[0]) + 1) / 2) * fPSD_psdX_pathlength[0]) +
            (((TMath.Sign(1, fPSD_psdX_pathlength[1]) + 1) / 2) * fPSD_psdX_pathlength[1])
        )

        if den_psd_x > 0:
            fPSD_Etrack23[0] = (
                ((((TMath.Sign(1, fPSD_EnergyX_corr[0]) + 1) / 2) * fPSD_EnergyX_corr[0]) +
                 (((TMath.Sign(1, fPSD_EnergyX_corr[1]) + 1) / 2) * fPSD_EnergyX_corr[1]))
                * 10. / den_psd_x
            )
        else:
            fPSD_Etrack23[0] = -999.

        for iplane in xrange(0, 6):
            fSTK_chargeX[iplane] = cluster_chargeX[iplane]
            fSTK_chargeY[iplane] = cluster_chargeY[iplane]
            fSTK_chargeX_etaCorr[iplane] = cluster_chargeX_etacorr[iplane]
            fSTK_chargeY_etaCorr[iplane] = cluster_chargeY_etacorr[iplane]

        fSTK_theta_correction[0] = track_correction
        fSTK_trackIP[0] = track_impact_point.x()
        fSTK_trackIP[1] = track_impact_point.y()
        fSTK_trackIP[2] = track_impact_point.z()
        fSTK_trackDirection[0] = track_slope.x()
        fSTK_trackDirection[1] = track_slope.y()

        for ibgo in xrange(0, 14):
            fBGO_EneLay[ibgo] = pev.pEvtBgoRec().GetELayer(ibgo) / 1000.
            fBGO_nhits_layers[ibgo] = pev.pEvtBgoRec().GetLayerHits()[ibgo]

        EBarLayers = [[] for _ in range(14)]
        PosBarLayers = [[] for _ in range(14)]
        PosXBarLayers = [[] for _ in range(14)]
        PosYBarLayers = [[] for _ in range(14)]
        rmsLayer = [[0.] for _ in range(14)]
        sumRMS = 0.

        for ihit in xrange(0, pev.NEvtBgoHits()):
            layID = pev.pEvtBgoHits().GetLayerID(ihit)
            hitE = pev.pEvtBgoHits().fEnergy[ihit]
            barNumber = (pev.pEvtBgoHits().fGlobalBarID[ihit] >> 6) & 0x1f
            EBarLayers[layID].append(hitE)
            PosBarLayers[layID].append(barNumber)
            PosXBarLayers[layID].append(pev.pEvtBgoHits().GetHitX(ihit))
            PosYBarLayers[layID].append(pev.pEvtBgoHits().GetHitY(ihit))

        zbgo = array('d', 14 * [-999.])
        for ibgo in xrange(0, 14):
            zbgo[ibgo] = 58.5 + ibgo * 29.0
            maxEneLay = 0
            maxEneLaypos = -9
            maxPosLay = -999
            maxPosX = -999.
            maxPosY = -999.
            if len(EBarLayers[ibgo]) != 0:
                maxEneLay = np.max(EBarLayers[ibgo])
                maxEneLaypos = np.argmax(EBarLayers[ibgo])
                maxPosLay = PosBarLayers[ibgo][maxEneLaypos]
                maxPosX = PosXBarLayers[ibgo][maxEneLaypos]
                maxPosY = PosYBarLayers[ibgo][maxEneLaypos]
                if ibgo % 2 == 0:
                    fBGO_cbgomax[ibgo] = maxPosY
                if ibgo % 2 == 1:
                    fBGO_cbgomax[ibgo] = maxPosX

            if ibgo % 2 == 0:
                fBGO_cbgostk[ibgo] = track_impact_point.y() + track_slope.y() * (zbgo[ibgo] - track_impact_point.z())
            if ibgo % 2 == 1:
                fBGO_cbgostk[ibgo] = track_impact_point.x() + track_slope.x() * (zbgo[ibgo] - track_impact_point.z())

            if maxEneLay > 0.:
                bar_less1 = -9
                bar_plus1 = -9
                if PosBarLayers[ibgo][maxEneLaypos] > 0 and PosBarLayers[ibgo][maxEneLaypos] < 21:
                    if (PosBarLayers[ibgo][maxEneLaypos] - 1) in PosBarLayers[ibgo]:
                        bar_less1 = PosBarLayers[ibgo].index(PosBarLayers[ibgo][maxEneLaypos] - 1)
                        if (PosBarLayers[ibgo][maxEneLaypos] + 1) in PosBarLayers[ibgo]:
                            bar_plus1 = PosBarLayers[ibgo].index(PosBarLayers[ibgo][maxEneLaypos] + 1)

                eCoreLayer = maxEneLay
                eCoreCoordLayer = maxEneLay * (maxPosX if ibgo % 2 else maxPosY)
                if bar_less1 != -9:
                    eCoreLayer += EBarLayers[ibgo][bar_less1]
                    eCoreCoordLayer += EBarLayers[ibgo][bar_less1] * (PosXBarLayers[ibgo][bar_less1] if ibgo % 2 else PosYBarLayers[ibgo][bar_less1])
                if bar_plus1 != -9:
                    eCoreLayer += EBarLayers[ibgo][bar_plus1]
                    eCoreCoordLayer += EBarLayers[ibgo][bar_plus1] * (PosXBarLayers[ibgo][bar_plus1] if ibgo % 2 else PosYBarLayers[ibgo][bar_plus1])
                eCoreCoordLayer /= eCoreLayer

                for ibar in xrange(0, len(EBarLayers[ibgo])):
                    delta = (PosXBarLayers[ibgo][ibar] - eCoreCoordLayer) if ibgo % 2 else (PosYBarLayers[ibgo][ibar] - eCoreCoordLayer)
                    rmsLayer[ibgo] += EBarLayers[ibgo][ibar] * delta * delta

            if maxEneLay > 0.:
                rmsLayer[ibgo] = np.sqrt(rmsLayer[ibgo] / sum(EBarLayers[ibgo]))
            else:
                rmsLayer[ibgo] = 0
            sumRMS += rmsLayer[ibgo]

        Xtr = sumRMS * sumRMS * sumRMS * sumRMS * (sum(EBarLayers[13]) / pev.pEvtBgoRec().GetTotalEnergy()) / 8000000.

        zpsd = [-324.7, -310.7, -298.5, -284.5]
        y_projection_fromBGO_to_PSD_up = pevspace_api.GetDirectionBGOInteceptY() + zpsd[0] * pevspace_api.GetDirectionBGOSlopeY()

        fBGOtrack_to_PSD_topY[0] = y_projection_fromBGO_to_PSD_up
        fBGO_theta[0] = theta_bgo_deg
        fSTK_theta[0] = theta_stk_deg
        fBGO_phi[0] = phi_bgo_deg
        fSTK_phi[0] = phi_stk_deg

        fBGO_EnergyG[0] = etot
        fBGO_EnergyG_SatCorr_ML_ions[0] = etot_satcorr_ions
        fBGO_EnergyG_SatCorr_ML_ions_v3[0] = etot_satcorr_ions_v3

        if is_mc:
            fBGO_EnergyG_Quench[0] = etot_quench
            fBGO_EnergyG_QuenchSatCorr_ML_ions[0] = etot_quench_satcorr_ions
            fBGO_EnergyG_QuenchSatCorr_ML_ions_v3[0] = etot_quench_satcorr_ions_v3

        fBGO_HET[0] = pev.pEvtHeader().GeneratedTrigger(3)
        fBGO_LET[0] = pev.pEvtHeader().GeneratedTrigger(4)
        fBGO_UBT[0] = pev.pEvtHeader().GeneratedTrigger(0)
        fBGO_UBE[0] = pev.pEvtHeader().EnabledTrigger(0)
        fBGO_nbars[0] = pev.NEvtBgoHits()
        fBGO_Xtr[0] = Xtr

        ftime[0] = pev.pEvtHeader().GetSecond()
        ftimems[0] = pev.pEvtHeader().GetMillisecond()

        if is_mc:
            fMC_EnergyT[0] = etot_truth
            fMC_theta[0] = theta_truth
            fMC_phi[0] = phi_truth
            fMC_stopZ[0] = stopZ_parent
            fMC_stopX[0] = stopX_parent
            fMC_stopY[0] = stopY_parent

        entry_ntupla += 1

        if psdchargeX_corr[0] == -100 and psdchargeX_corr[1] == -100:
            print " -100 find!"
            print sum_len_psd_vec_chargeY, ' ', sum_len_psd_vec_chargeX

        for ibgo in range(14):
            h_bgo_ene_lay[ibgo].Fill((pev.pEvtBgoRec().GetELayer(ibgo)) / 1000.)

        newtree.Fill()

    dmpch.Terminate()

    print "All statistic:   ", h_energy[0].GetEntries()
    for i in range(1, ncuts):
        percent = h_energy[i].GetEntries() / h_energy[0].GetEntries() * 100 if h_energy[0].GetEntries() else 0.
        print cut_string[i], "  n events (all) :  ", h_energy[i].GetEntries(), "percent: ", percent
        if is_mc:
            percent_good = h_energy_good[i].GetEntries() / h_energy_good[0].GetEntries() * 100 if h_energy_good[0].GetEntries() else 0.
            print cut_string[i], "  n events (good):  ", h_energy_good[i].GetEntries(), "percent: ", percent_good
            percent_bad = h_energy_bad[i].GetEntries() / h_energy_bad[0].GetEntries() * 100 if h_energy_bad[0].GetEntries() else 0.
            print cut_string[i], "  n events (bad) :  ", h_energy_bad[i].GetEntries(), "percent: ", percent_bad

    print "Wrong ID for MC ", wrong_id
    print "Events not passint the new BGO-STK match ", ntrack_newsel
    tf_skim.Write()
    tf_skim.Close()

    if opts.data:
        oout = root_dir + '/' + opts.diag_output
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
        for iplane in range(0, 6):
            h_stk_cluster_XvsY[iplane].Write()
        for ibgo in range(14):
            h_bgo_ene_lay[ibgo].Write()

        tf.Write()
        tf.Close()


if __name__ == "__main__":
    main()
    print datetime.now() - start
