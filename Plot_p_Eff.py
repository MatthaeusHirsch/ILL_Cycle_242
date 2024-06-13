import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat



def t_to_vs(countsumT, chopper):
    v_corr_countsum = []
    v_corr_countsum_uncertainty = []
    tp = np.linspace(0, chopper[3]-chopper[2], len(countsumT))
    vp = 0.631/tp
    vp = np.asarray(vp)

    for i in range(len(countsumT)-1):
        #v_corr_countsum.append(countsumT[i]/(vp[i]-vp[i+1]))
        v_corr_countsum.append(countsumT[i])
        ##### rausgenommen
        v_corr_countsum_uncertainty.append(np.sqrt(np.absolute(countsumT[i]))/(vp[i]-vp[i+1]))
    v_corr_countsum = np.asarray(v_corr_countsum)
    v_corr_countsum_uncertainty = np.asarray(v_corr_countsum_uncertainty)
    return vp, v_corr_countsum, v_corr_countsum_uncertainty



def tRate_to_vs(countsumT, uncertainty_countsumT, bintimes):
    v_Spect = np.multiply(countsumT , np.power(bintimes, 2))* 1/0.631
    uncertainty_v_Spect = np.multiply(uncertainty_countsumT , np.power(bintimes, 2)) * 1/0.631
    return v_Spect, uncertainty_v_Spect



def getTSpectrum(time, chopper):
    size = 100
    #Add bins
    countsum = np.zeros(size)
    for i in range(len(chopper)-1):
        n, bins = np.histogram(time, bins=size, range=(chopper[i], chopper[i+1]))

        countsum = countsum + n
    return [countsum, np.sqrt(countsum)]


def read_files(pfad):
    data = np.loadtxt(pfad, skiprows=1, delimiter=';', usecols=(1,2,3,4))
    ch = data[:,0]
    time = data[:,1]/1e12
    ph = data[:,2]
    
    #Dets
    chopper = []
    td = []

    for i in range(len(ch)):
        if ch[i] == 1:
            td.append(time[i])
        elif ch[i] == 0:
            if ph[i] > 7000:
                chopper.append(time[i])

    chopper = np.asarray(chopper)+0.02182 #Chopper offset 21.82ms
    td = np.asarray(td)

    return (td, chopper)


def substractBackground(tSpectrum):
    fitting_spectrum = [tSpectrum[i] for i in range(2*len(tSpectrum)//3, len(tSpectrum))]
    tSpectrum = np.asarray(tSpectrum)
    background = np.mean(fitting_spectrum)
    corrected_t_Spectrum = tSpectrum - background
    return corrected_t_Spectrum


def getSpectrum(run):
    
    bins = 100

    td = []
    tc = []

    if noOfFiles > 1 :
        for i in range(noOfFiles):
            pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}_{i}.CSV'
            result = read_files(pfad)
            td.extend(result[0])
            tc.extend(result[1])
    else: 
        pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'
        result = read_files(pfad)
        td.extend(result[0])
        tc.extend(result[1])

    td = np.asarray(td)
    tc = np.asarray(tc)

    
    tSpectrum_with_BG = getTSpectrum(td, tc)
    tSpectrum = substractBackground(tSpectrum_with_BG[0])

    tRate = np.divide(tSpectrum, max(tc))*(1/ ((tc[4]-tc[3]) / bins))
    uncertainty_tRate = np.divide(np.sqrt(np.absolute(tSpectrum)), max(tc))*(1/ ((tc[4]-tc[3]) / bins))


    
    bin_mid = np.add(np.linspace(0,0.4,bins+1) , (0.4/(2*bins)))[:-1]

    vRateSpectrum , uncertainty_vRateSpectrum = tRate_to_vs(tRate , uncertainty_tRate, bin_mid)


    return np.divide(0.631 , bin_mid), vRateSpectrum, uncertainty_vRateSpectrum


def calc_PolEff_i(X_i, X_i2, X_0):
    p = unp.sqrt(unp.sqrt((X_i * X_i2/X_0)**2))

    return p

def plot_runs_with_pol_eff(runs, xlim, lineplot = False):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed
    integrals = []
    i = 1
    Spectra = []
    for run in runs:
        curr_run = getSpectrum(run)
        Spectra.append(curr_run)
        color = colors[i % len(colors)]  # Select a color from the list
        plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
        if lineplot == True:
            plt.plot(curr_run[0], curr_run[1], color=color)

        integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))
        i += 1
    Spectra = np.asarray(Spectra)
    
    
    # Convert Spectra values and uncertainties to ufloat arrays
    Spectra_values = [unp.uarray(Spectra[i][1], Spectra[i][2]) for i in range(len(Spectra))]

    # Calculate X12, X32, X31, X34 using uncertainties
    X12 = (Spectra_values[0] - Spectra_values[1]) / (Spectra_values[0] + Spectra_values[1])
    X32 = (Spectra_values[2] - Spectra_values[3]) / (Spectra_values[2] + Spectra_values[3])
    X31 = (Spectra_values[4] - Spectra_values[5]) / (Spectra_values[4] + Spectra_values[5])
    X34 = (Spectra_values[6] - Spectra_values[7]) / (Spectra_values[6] + Spectra_values[7])

    # Assuming calc_PolEff_i function works correctly with ufloat inputs
    upper_bound = 10
    lower_bound = -64

    p1 = calc_PolEff_i(X12, X31, X32)[upper_bound:lower_bound]
    p2 = calc_PolEff_i(X12, X32, X31)[upper_bound:lower_bound]
    p3 = calc_PolEff_i(X31, X32, X12)[upper_bound:lower_bound]
    p4 = X34[upper_bound:lower_bound] / p3
    


    plt.errorbar(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p1),unp.std_devs(p1), fmt='.', label='Polarizing-Eff. P1')
    if lineplot == True:
        plt.plot(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p1))

    plt.errorbar(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p2), unp.std_devs(p2), fmt='.', label='Polarizing-Eff. P2')
    if lineplot == True:
        plt.plot(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p2))

    plt.errorbar(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p3), unp.std_devs(p3), fmt='.', label='Polarizing-Eff. P3')
    if lineplot == True:
        plt.plot(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p3))

    plt.errorbar(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p4), unp.std_devs(p4), fmt='.', label='Polarizing-Eff. P4')
    if lineplot == True:
        plt.plot(curr_run[0][upper_bound:lower_bound], unp.nominal_values(p4))

    plt.grid()
    plt.legend()
    plt.xlabel('velocity in [ms]')
    plt.ylabel('UCN count rate [(m/s)$^{-1}$ * s$^{-1}§] ')
    plt.xlim(xlim)
    plt.ylim(0,7)

    
    # Save the figure with the runs in the filename
    if lineplot == True:
        plt.savefig(f'pf2data/Plots/PolEffs_and_V-Spectrum_with_lines_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    else: 
        plt.savefig(f'pf2data/Plots/PolEffs_and_V-Spectrum_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    

    plt.show()









#run the code

noOfFiles = 1
plot_runs_with_pol_eff([164, 165, 168, 169, 172, 173, 174, 175], [0,25], lineplot=True) # runs to plot and velocity limits to plot in; Order: first _0 then _1: T12, T32 , T31 , T34
