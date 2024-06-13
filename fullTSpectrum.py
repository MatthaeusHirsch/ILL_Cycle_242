import numpy as np
import matplotlib.pyplot as plt

def t_to_vs(countsumT, chopper):
    v_corr_countsum = []
    v_corr_countsum_uncertainty = []
    tp = np.linspace(0, chopper[3]-chopper[2], len(countsumT))
    vp = 0.631/tp
    vp = np.asarray(vp)

    for i in range(len(countsumT)-1):
        v_corr_countsum.append(countsumT[i]/(vp[i]-vp[i+1]))
        v_corr_countsum_uncertainty.append(np.sqrt(countsumT[i])/(vp[i]-vp[i+1]))
    v_corr_countsum = np.asarray(v_corr_countsum)
    v_corr_countsum_uncertainty = np.asarray(v_corr_countsum_uncertainty)
    return vp, v_corr_countsum, v_corr_countsum_uncertainty


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
            if ph[i] > 6000:
                chopper.append(time[i])

    chopper = np.asarray(chopper)+0.02182 #Chopper offset 21.82ms
    td = np.asarray(td)


    return (td, chopper)



def substractBackground(tSpectrum):
    fitting_spectrum = [tSpectrum[i] for i in range(2*len(tSpectrum)//3, len(tSpectrum))]
    tSpectrum = np.asarray(tSpectrum)
    background = np.mean(fitting_spectrum)
    corrected_t_Spectrum = tSpectrum-background
    return corrected_t_Spectrum





def plot(x, y , yErr, xlim_low = 0, xlim_high = 25):
    plt.errorbar(x, y, yErr, fmt='.', color='blue')
    plt.grid()
    plt.xlabel('time in [s]')
    plt.ylabel('UCN count rate [s$^{-1}$]')
    plt.xlim([xlim_low,xlim_high])
    plt.savefig( f'pf2data/Plots/T-Spectrum_run{run}.png', dpi = 300 )
    plt.show()



def main(xlim_low = 0, xlim_high = 25):
    

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

    xt = np.linspace(0,0.4,100)



    plot(xt, tSpectrum/max(tc), np.sqrt(np.absolute(tSpectrum))/max(tc), xlim_low, xlim_high)


#run the code
run = 185
noOfFiles = 1
main(xlim_low = 0, xlim_high = 0.4)