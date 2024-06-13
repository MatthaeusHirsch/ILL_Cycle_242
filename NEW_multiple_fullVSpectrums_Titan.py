import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import unumpy as unp


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



def tRate_to_vs(TRate, uncertainty_TRate, bintimes):
    ufloats = unp.uarray(TRate, uncertainty_TRate)
    v_Spect = np.multiply(ufloats , np.power(bintimes, 2))* 1/0.631
 
    return unp.nominal_values(v_Spect), unp.std_devs(v_Spect)

def tRate_to_vs_Titan(TRate, uncertainty_TRate, bintimes):
    ufloats = unp.uarray(TRate, uncertainty_TRate)
    v_Spect = np.multiply(ufloats , np.power(bintimes, 2))* 1/0.5
    print('Titan used!')
    return unp.nominal_values(v_Spect), unp.std_devs(v_Spect)



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


    return np.divide(0.631 , bin_mid), vRateSpectrum, uncertainty_vRateSpectrum, tRate, uncertainty_tRate, bin_mid


def getSpectrumTitan(run):
    
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
    vRateSpectrum , uncertainty_vRateSpectrum = tRate_to_vs_Titan(tRate , uncertainty_tRate, bin_mid)


    return np.divide(0.5 , bin_mid), vRateSpectrum, uncertainty_vRateSpectrum, tRate, uncertainty_tRate, bin_mid


def flip_Array(x):
    x_resorted = []
    for i in range(len(x)):
        x_resorted.append(x[len(x)-1-i])
    return np.asarray(x_resorted)

def getfraction(S1, S2):


    x1 = flip_Array(S1[0])
    x2 = flip_Array(S2[0])
    print(x2)
    y1 = flip_Array(S1[1])
    y2 = flip_Array(S2[1])
    print(y1/y2)

    
    # Create a common set of x-values (you can adjust the range and density as needed)
    x_common = np.linspace(3.5, 25, 100)
    # Interpolate y-values for both functions at these common x-values
    y1_interp = np.interp(x_common, x1, y1)
    y2_interp = np.interp(x_common, x2, y2)

    frac = y1_interp / y2_interp

    return x_common, frac


def plot_runs_vSpect_with_frac(runs, xlim = [0,25], lineplot = False):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed
    integrals = []
    i = 0
    dataRuns = []
    for run in runs:
        if run >= 190:
            curr_run = getSpectrumTitan(run)
            dataRuns.append(curr_run)
            color = colors[i % len(colors)]  # Select a color from the list
            plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
            if lineplot == True:
                plt.plot(curr_run[0], curr_run[1], color=color)

            integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))
            #print(curr_run[0][18:-62])
            
        else:
            if run == 185:
                curr_run = getSpectrum(run)
                dataRuns.append(curr_run)
                color = colors[i % len(colors)]  # Select a color from the list
                plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
                if lineplot == True:
                    plt.plot(curr_run[0], curr_run[1], color=color)

                integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))
                #print(curr_run[0][18:-62])
                
            else:                                                                                                                                                                                                                                                                                                                                                               
                curr_run = getSpectrum(run)
                dataRuns.append(curr_run)
                color = colors[i % len(colors)]  # Select a color from the list
                plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
                if lineplot == True:
                    plt.plot(curr_run[0], curr_run[1], color=color)

                integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))    
        print(f'run{run} is in pos {i}')
        i += 1


    dataRuns = np.asarray(dataRuns)
    frac = getfraction(dataRuns[3], dataRuns[2])
    print(frac[0], frac[1])
    plt.plot(frac[0], frac[1], label = 'Fraction')

    plt.grid()
    plt.legend()
    plt.xlabel('velocity in [ms]')
    plt.ylabel('UCN count rate [(m/s)$^{-1}$]')
    plt.xlim(xlim)
    plt.ylim(0,1)


    
    # Save the figure with the runs in the filename
    if lineplot == True:
        plt.savefig(f'pf2data/Plots/correct_V-Spectrum_with_lines_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    else: 
        plt.savefig(f'pf2data/Plots/correct_V-Spectrum_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    
    print(integrals)
    plt.show()




def plot_runs_vSpect(runs, xlim = [0,25], lineplot = False):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed
    integrals = []
    i = 0
    for run in runs:
        if run >= 190:
            curr_run = getSpectrumTitan(run)
            color = colors[i % len(colors)]  # Select a color from the list
            
            if run == 193:
                plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
                if lineplot == True:
                    plt.plot(curr_run[0], curr_run[1], color=color)
            else:
                plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
                if lineplot == True:
                    plt.plot(curr_run[0], curr_run[1], color=color)

            integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))
            print(curr_run[0][18:-62])
            i += 1
        else:
            curr_run = getSpectrum(run)
            color = colors[i % len(colors)]  # Select a color from the list
            plt.errorbar(curr_run[0], curr_run[1], curr_run[2], fmt='.', label=f'run{run}', color=color)
            if lineplot == True:
                plt.plot(curr_run[0], curr_run[1], color=color)

            integrals.append(-np.trapz(curr_run[1][18:-62], curr_run[0][18:-62]))
            print(curr_run[0][18:-62])
            i += 1
        print(f'run{run} is in pos {i}')


    plt.grid()
    plt.legend()
    plt.xlabel('velocity in [ms]')
    plt.ylabel('UCN count rate [(m/s)$^{-1}$]')
    plt.xlim(xlim)


    
    # Save the figure with the runs in the filename
    if lineplot == True:
        plt.savefig(f'pf2data/Plots/correct_V-Spectrum_with_lines_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    else: 
        plt.savefig(f'pf2data/Plots/correct_V-Spectrum_runs{"_".join(map(str, runs))}__xlim_{xlim[0]}-{xlim[1]}mps.png', dpi=300)
    
    print(integrals)
    plt.show()


def plot_runs_tSpect(runs, lineplot = False):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed
    integrals = []
    i = 0
    dataRuns = []
    for run in runs:
        curr_run = getSpectrum(run)
        dataRuns.append(curr_run[1])
        color = colors[i % len(colors)]  # Select a color from the list
        plt.errorbar(curr_run[5], curr_run[3], curr_run[4], fmt='.', label=f'run{run}', color=color)
        if lineplot == True:
            plt.plot(curr_run[5], curr_run[3], color=color)
        print(curr_run[5])
        integrals.append(np.trapz(curr_run[3], curr_run[5]))
        i += 1
    dataRuns = np.asarray(dataRuns)



    plt.grid()
    plt.legend()
    plt.xlabel('velocity in [ms]')
    plt.ylabel('UCN count rate [(s)$^{-2}$]')



    
    # Save the figure with the runs in the filename
    if lineplot == True:
        plt.savefig(f'pf2data/Plots/correct_T-Spectrum_with_lines_runs{"_".join(map(str, runs))}.png', dpi=300)
    else: 
        plt.savefig(f'pf2data/Plots/correct_T-Spectrum_runs{"_".join(map(str, runs))}.png', dpi=300)
    
    print(integrals)
    plt.show()






#run the code

noOfFiles = 1
plot_runs_vSpect_with_frac([186, 192, 185, 193], lineplot=True) # runs to plot and velocity limits to plot in
