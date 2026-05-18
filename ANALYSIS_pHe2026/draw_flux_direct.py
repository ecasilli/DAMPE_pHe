import numpy as np
import matplotlib.pyplot as plt


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

    flux_E27 = flux * energy**2.7
    stat_err_E27 = (stat_err * flux_E27)/100.
    sys_err_E27 = sys_err * flux_E27
    sys_had_E27 = sys_had * flux_E27

    sys_err_up = flux_E27 + sys_err_E27
    sys_err_down = flux_E27 - sys_err_E27
    had_err_up = flux_E27 + sys_had_E27
    had_err_down = flux_E27 - sys_had_E27

    return energy, flux_E27, stat_err_E27, sys_err_E27, sys_had_E27, sys_err_up, sys_err_down, had_err_up, had_err_down

dampe_2026 = "flux_spectrum_pHe_2026_MLionsv3_smooth.dat"
### flux_spectrum_pHe_2026_MLionsv3.dat no smoothing
dampe_v1 = "flux_spectrum_pHe_2026.dat"

def load_dampe2026(filename):
    data = np.loadtxt(filename, skiprows=2)
    energy = data[:, 0]
    flux = data[:, 1]
    stat_err = data[:, 2]
    energy_err = data[:, 3]

    flux_E27 = flux * energy**2.7
    stat_err_E27 = stat_err * energy**2.7

    return energy, flux_E27, stat_err_E27

####### DIRECT EXPERIMENT DATA ######

# header: Energy FluxProton ErrProton FluxHelium ErrHelium  ErrEnergy 
### the flux have to be multiplied E^2.7
### separated data for proton and helium, sum them together to get the total flux
atic = "../pHe2024/SPECTRUM/Direct/ATIC2_ProtonStat_HeliumStat.dat"
cream = "../pHe2024/SPECTRUM/Direct/CREAM_Proton_Helium.dat"
nucleon = "../pHe2024/SPECTRUM/Direct/Nucleon_p_He.dat"

def load_data(filename):
    data = np.loadtxt(filename, skiprows=0)
    energy = data[:, 0]
    flux_proton = data[:, 1]
    err_proton = data[:, 2]
    flux_helium = data[:, 3]
    err_helium = data[:, 4]
    err_energy = data[:, 5]

    flux_pHe = (flux_proton + flux_helium)* energy**2.7
    err_flux_pHe = np.sqrt(err_proton**2 + err_helium**2) * energy**2.7

    return energy, flux_pHe, err_flux_pHe



plt.figure(figsize=(10, 6))

# Load DAMPE data
energy_dampe, flux_dampe, stat_err_dampe, sys_err_dampe, sys_had_dampe, sys_err_up_dampe, sys_err_down_dampe, had_err_up_dampe, had_err_down_dampe = load_dampe2024(dampe_2024)
energy_dampe2026, flux_dampe2026, stat_err_dampe2026 = load_dampe2026(dampe_2026)
energy_v1, flux_v1, stat_err_v1 = load_dampe2026(dampe_v1)

# Load direct experiment datas
energy_atic, flux_atic, err_atic = load_data(atic)
energy_cream, flux_cream, err_cream = load_data(cream)
energy_nucleon, flux_nucleon, err_nucleon = load_data(nucleon)

def plot_flux(energy, flux, stat_err, sys_err_up, sys_err_down, had_err_up, had_err_down, label, color):
    plt.errorbar(energy, flux, yerr=stat_err, fmt='.', label=label + ' (stat)', markersize=5, color=color)
    plt.fill_between(energy, sys_err_down, sys_err_up, alpha=0.3, label=label + ' (sys)', color=color)
    plt.fill_between(energy, had_err_down, had_err_up, alpha=0.2, label=label + ' (total)', color=color)

plot_flux(energy_dampe, flux_dampe, stat_err_dampe, sys_err_up_dampe, sys_err_down_dampe, had_err_up_dampe, had_err_down_dampe, 'DAMPE 2024', 'blue')
plt.errorbar(energy_atic, flux_atic, yerr=err_atic, fmt='s', label='ATIC', markersize=3, alpha=0.7)
plt.errorbar(energy_cream, flux_cream, yerr=err_cream, fmt='^', label='CREAM', markersize=3, alpha=0.7)
plt.errorbar(energy_nucleon, flux_nucleon, yerr=err_nucleon, fmt='d', label='Nucleon', markersize=3, alpha=0.7)

#plt.errorbar(energy_v1[:-4], flux_v1[:-4], yerr=stat_err_v1[:-4], fmt='s', label='DAMPE v1', markersize=6)
plt.errorbar(energy_dampe2026[:-4], flux_dampe2026[:-4], yerr=stat_err_dampe2026[:-4], fmt='o', label='DAMPE 2026', markersize=5)

plt.xscale('log')
#plt.yscale('log')
plt.xlim(10, 1e7)

plt.xlabel('Energy (GeV)', fontsize=14)
plt.ylabel(r'$E^{2.7} \times$ Flux (m$^{-2}$ s$^{-1}$ sr$^{-1}$ GeV$^{1.7}$)', fontsize=14)


plt.legend(fontsize=14)
plt.grid(True, which='both', ls='--', lw=0.5)


#### save the plot
#plt.savefig('dampe_flux_comparison_v1-v3.png', dpi=300)
#plt.savefig('direct_flux_comparison.png', dpi=300)

plt.show()
