import numpy as np
import matplotlib.pyplot as plt

gamma = 2.6
###### DAMPE DATA ######

# dampe header: TotalEnergy Flux Stat_err EnergyErr Sys_err SysHad 
dampe_2024 = "../pHe2024/SPECTRUM/Direct/DAMPE_p+He_72M26.dat" ### the flux have to be multiplied E^2.7; energy GeV
def load_dampe2024(filename):
    data = np.loadtxt(filename, skiprows=0)
    energy = data[:, 0]
    flux = data[:, 1]
    stat_err = data[:, 2]
    energy_err = data[:, 3]
    sys_err = data[:, 4]
    sys_had = data[:, 5]

    flux_E26 = flux * energy**gamma
    stat_err_E26 = (stat_err * flux_E26)/100.
    sys_err_E26 = sys_err * flux_E26
    sys_had_E26 = sys_had * flux_E26

    sys_err_up = flux_E26 + sys_err_E26
    sys_err_down = flux_E26 - sys_err_E26
    had_err_up = flux_E26 + sys_had_E26
    had_err_down = flux_E26 - sys_had_E26

    return energy, flux_E26, stat_err_E26, sys_err_E26, sys_had_E26, sys_err_up, sys_err_down, had_err_up, had_err_down

dampe_2026 = "flux_spectrum_pHe_2026_MLionsv3_smooth.dat"
def load_dampe2026(filename):
    data = np.loadtxt(filename, skiprows=2)
    energy = data[:, 0]
    flux = data[:, 1]
    stat_err = data[:, 2]
    energy_err = data[:, 3]

    flux_E26 = flux * energy**gamma
    stat_err_E26 = stat_err * energy**gamma

    return energy, flux_E26, stat_err_E26

####### INDIRECT EXPERIMENT DATA ######

argo = "../pHe2024/SPECTRUM/Indirect/ARGO2New.txt"
# Energy_ARGOTeV Flux_ARGO_E265 Stat_err_up_ARGO Stat_err_down_ARGO Sys_err_up_ARGO Sys_err_down_ARGO Cont_ARGO
def load_argo(filename):
    data = np.loadtxt(filename, skiprows=0)
    energy = data[:, 0]
    flux = data[:, 1]
    stat_err_up = data[:, 2]
    stat_err_down = data[:, 3]
    sys_err_up = data[:, 4]
    sys_err_down = data[:, 5]
    cont = data[:, 6]

    energy_GeV = 10**energy * 1000
    flux_E26 = flux/1e12 * energy_GeV**gamma
    stat_err_down_E26 = stat_err_down/1e12 * energy_GeV**gamma
    stat_err_up_E26 = stat_err_up/1e12 * energy_GeV**gamma
    sys_err_up_E26 = flux_E26 * sys_err_up/100.0
    #sys_err_down_E26 = flux_E26 * sys_err_down/100.0
    cont = cont * flux_E26/100.0 
    sys_err_down_E26 = flux_E26 - cont - (flux_E26 * sys_err_down/100.0)

    return energy_GeV, flux_E26, stat_err_down_E26, stat_err_up_E26, sys_err_up_E26, sys_err_down_E26

energy_argo, flux_argo, stat_err_down_argo, stat_err_up_argo, sys_err_up_argo, sys_err_down_argo = load_argo(argo)

hawc = "../pHe2024/SPECTRUM/Indirect/HAWC2022.txt"
# Energy_HAWC Flux_HAWC Stat_err_HAWC Sys_err_up_HAWC Sys_err_down_HAWC
def load_hawc(filename):
    data = np.loadtxt(filename, skiprows=0)
    energy = data[:, 0]
    flux = data[:, 1]
    stat_err = data[:, 2]
    sys_err_up = data[:, 3]
    sys_err_down = data[:, 4]

    flux_E26 = flux * energy**gamma
    stat_err_E26 = stat_err * energy**gamma
    sys_err_down_E26 = flux_E26 - (sys_err_down * energy**gamma)
    sys_err_up_E26 = flux_E26 + (sys_err_up * energy**gamma)

    return energy, flux_E26, stat_err_E26, sys_err_up_E26, sys_err_down_E26

energy_hawc, flux_hawc, stat_err_hawc, sys_err_up_hawc, sys_err_down_hawc = load_hawc(hawc)

hawc_2025 = "../pHe2024/SPECTRUM/Indirect/HAWC_ICRC25.dat"
#E(GeV)	Flux * E^2,6	Flux err up	Flux err donw	Sys Up	Sys Down
energy_hawc_2025, flux_hawc_2025, sys_err_up_hawc_2025, sys_err_down_hawc_2025 = np.loadtxt(hawc_2025, skiprows=1, unpack=True, usecols=(0,1,4,5))
energy_hawc_2025 = 10**energy_hawc_2025

kascade_qgsjet01 = "../pHe2024/SPECTRUM/Indirect/KASCADE_QGSjet01.txt"
# Energy_KASCADE_QGSJet_GeV Flux_KASCADE_Protons Err_down_KASCADE_Protons Err_up_KASCADE_Protons Flux_KASCADE_Helium Err_down_KASCADE_Helium Err_up_KASCADE_Helium
kascade_sibyll = "../pHe2024/SPECTRUM/Indirect/KASCADE_SIBYLL.txt"
# Energy_KASCADE_SIBYLL_GeV Flux_KASCADE_SIBYLL_Protons Err_up_KASCADE_SIBYLL_Protons Err_down_KASCADE_SIBYLL_Protons Flux_KASCADE_SIBYLL_Helium Err_up_KASCADE_SIBYLL_Helium Err_down_KASCADE_SIBYLL_Helium 

def load_kascade(filename):
    data = np.loadtxt(filename, skiprows=0)
    energy = data[:, 0]
    flux_proton = data[:, 1]
    err_down_proton = data[:, 2]
    err_up_proton = data[:, 3]
    flux_helium = data[:, 4]
    err_down_helium = data[:, 5]
    err_up_helium = data[:, 6]

    flux_total = (flux_proton + flux_helium) * energy**gamma
    err_up_E26 = np.sqrt((err_up_proton * energy**gamma)**2 + (err_up_helium * energy**gamma)**2)
    err_down_E26 = np.sqrt((err_down_proton * energy**gamma)**2 + (err_down_helium * energy**gamma)**2)

    return energy, flux_total, err_down_E26, err_up_E26

energy_kascade_qgsjet01, flux_kascade_qgsjet01, err_down_kascade_qgsjet01, err_up_kascade_qgsjet01 = load_kascade(kascade_qgsjet01)
energy_kascade_sibyll, flux_kascade_sibyll, err_down_kascade_sibyll, err_up_kascade_sibyll = load_kascade(kascade_sibyll)

ice_top = "../pHe2024/SPECTRUM/Indirect/IceTop_ICRC25.dat"
#E GeV	Flux E^2.6	Flux Err Up	Flux Err Down
energy_icetop, flux_icetop, err_up_icetop, err_down_icetop = np.loadtxt(ice_top, skiprows=1, unpack=True, usecols=(0,1,2,3))
energy_icetop = 10**energy_icetop

'''
lhaaso = "../pHe2024/SPECTRUM/Indirect/LHAASO_2026.dat"
#E (PeV)	Flux * E^2.75	Flux Err Up	
energy_lhaaso, flux_lhaaso, err_up_lhaaso = np.loadtxt(lhaaso, skiprows=1, unpack=True, usecols=(0,1,2))
energy_lhaaso = 10**energy_lhaaso * 1e6
flux_lhaaso_E26 = flux_lhaaso * (energy_lhaaso**gamma / energy_lhaaso**2.75)  # fix the normalization
err_lhaaso_E26 = err_up_lhaaso * (energy_lhaaso**gamma / energy_lhaaso**2.75)
'''


#### PLOT 
plt.figure(figsize=(10, 6))

def plot_flux(energy, flux, stat_err, sys_err_up, sys_err_down, had_err_up, had_err_down, label, color):
    plt.errorbar(energy, flux, yerr=stat_err, fmt='.', label=label + ' (stat)', markersize=5, color=color)
    plt.fill_between(energy, sys_err_down, sys_err_up, alpha=0.3, label=label + ' (sys)', color=color)
    plt.fill_between(energy, had_err_down, had_err_up, alpha=0.2, label=label + ' (total)', color=color)

#### DAMPE
energy_dampe, flux_dampe, stat_err_dampe, sys_err_dampe, sys_had_dampe, sys_err_up_dampe, sys_err_down_dampe, had_err_up_dampe, had_err_down_dampe = load_dampe2024(dampe_2024)
energy_dampe2026, flux_dampe2026, stat_err_dampe2026 = load_dampe2026(dampe_2026)

plot_flux(energy_dampe, flux_dampe, stat_err_dampe, sys_err_up_dampe, sys_err_down_dampe, had_err_up_dampe, had_err_down_dampe, 'DAMPE 2024', 'blue')
plt.errorbar(energy_dampe2026[:-4], flux_dampe2026[:-4], yerr=stat_err_dampe2026[:-4], fmt='o', label='DAMPE 2026', markersize=5)

#### INDIRECT EXPERIMENTS
plt.errorbar(energy_argo, flux_argo, yerr=[stat_err_down_argo, stat_err_up_argo], fmt='s', label='ARGO', markersize=5, alpha=0.7)

#plt.errorbar(energy_hawc, flux_hawc, yerr=stat_err_hawc, fmt='^', label='HAWC', markersize=5, alpha=0.7)
plt.errorbar(energy_hawc_2025, flux_hawc_2025, yerr=[sys_err_down_hawc_2025, sys_err_up_hawc_2025], fmt='v', label='HAWC', markersize=5, alpha=0.7)

plt.errorbar(energy_kascade_qgsjet01, flux_kascade_qgsjet01, yerr=[err_down_kascade_qgsjet01, err_up_kascade_qgsjet01], fmt='d', label='KASCADE QGSJet01', markersize=5, alpha=0.7)
plt.errorbar(energy_kascade_sibyll, flux_kascade_sibyll, yerr=[err_down_kascade_sibyll, err_up_kascade_sibyll], fmt='x', label='KASCADE SIBYLL', markersize=5, alpha=0.7)

plt.errorbar(energy_icetop, flux_icetop, yerr=[err_down_icetop, err_up_icetop], fmt='P', label='IceTop', markersize=5, alpha=0.7)

#plt.errorbar(energy_lhaaso, flux_lhaaso_E26, yerr=err_lhaaso_E26, fmt='X', label='LHAASO 2026', markersize=5, alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlim(10, 5*1e8)
#plt.ylim(1e4, 1e5)

plt.xlabel('Energy (GeV)', fontsize=14)
plt.ylabel(r'$E^{2.6} \times$ Flux (m$^{-2}$ s$^{-1}$ sr$^{-1}$ GeV$^{1.6}$)', fontsize=14)


plt.legend(fontsize=14)
plt.grid(True, which='both', ls='--', lw=0.5)


#### save the plot

#plt.savefig('indirect_flux_comparison.png', dpi=300)

plt.show()