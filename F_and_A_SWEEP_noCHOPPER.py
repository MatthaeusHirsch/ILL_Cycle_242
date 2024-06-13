import numpy as np
import matplotlib.pyplot as plt


#Calculate v from M in T
def calc_ts(M):
    Emax = (210+M*60)*1e-9
    Emin = (210-60*M)*1e-9

    vmax = np.sqrt(2*Emax/(939e6/3e8**2))
    vmin = np.sqrt(2*Emin/(939e6/3e8**2))

    tmax = 0.631/vmin
    tmin = 0.631/vmax
    return (tmax , tmin)



def v_to_vs(countsum, M):
    v_corr_countsum = []
    tmax , tmin = calc_ts(M)
    tp = np.linspace(tmin, tmax, len(countsum))
    vp = 0.631/tp
    vp = np.asarray(vp)
    for i in range(len(countsum)-1):
        v_corr_countsum.append(countsum[i]/(vp[i]-vp[i+1]))

    v_corr_countsum = np.asarray(v_corr_countsum)

    return v_corr_countsum

def getCountsOLD(M, data):
    tmax , tmin = calc_ts(M)
        
    det = 1
    size = 12

    ch = data[:,0]
    time = data[:,1]/1e12
    ph = data[:,2]


    #Dets
    d = []
    chopper = []
    td = []

    for i in range(len(ch)):
        if ch[i] == det:
            if time[i]<1000:
                d.append(ph[i])
                td.append(time[i])
        elif ch[i] == 0:
            if ph[i] > 0:
                if time[i]<1000:
                    chopper.append(time[i])

    chopper = np.asarray(chopper)+0.02182 #Chopper offset 21.82ms

    #Add bins
    countsum = np.zeros(size)

    for i in range(len(chopper)-1):
        n, bins = np.histogram(td, bins=size, range=(chopper[i]+tmin, chopper[i]+tmax))

        countsum = countsum + n

    return sum(v_to_vs(countsum, 1.2))



def getCountsNoChopper(time):
    return len(time)


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

"""
def splitTimes(time, chopper, intervall):
    hilf1 = []
    hilf2 = []
    for i in range(len(time)):
        if time[i] < (len(hilf2)+1)*intervall:
            hilf1.append(time[i])
        else:
            hilf2.append(hilf1)
    
    hilf2.append(hilf1)

    return hilf2
"""

def splitTimes(timestamps1, timestamps2, interval, noOfInterval):
    result1 = []
    result2 = []

    start = 0
    end = interval
    EndEnd = noOfInterval*interval
    
    while start <= EndEnd:
        subarray1 = timestamps1[(timestamps1 >= start) & (timestamps1 < end)]
        subarray2 = timestamps2[(timestamps2 >= start) & (timestamps2 < end)]
        
        result1.append(subarray1)
        result2.append(subarray2)
        
        start = end
        end = start + interval

    return result1, result2





def main():

    tdGes = []
    tcGes = []

    for i in range(noOfFiles):
        pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}_{i}.CSV'
        result = read_files(pfad)
        tdGes.extend(result[0])
        tcGes.extend(result[1])

    tdGes = np.asarray(tdGes)
    tcGes = np.asarray(tcGes)

    #print(tcGes)

    tdSplit, tcSplit = splitTimes(tdGes, tcGes, Minterval_in_s, noOfFrequencies * noOfAmplitudes)
    
    N_ROI = [] #Counts inside region of interest in the i-th time-interval
    for j in range(noOfAmplitudes):
        N_ROI_local = []
        for i in range(noOfFrequencies):
            N_ROI_local.append(getCountsNoChopper(tdSplit[j*noOfFrequencies+i])) #N_ROI_local = Array with j th Amplitude and the frequency Sweep
        N_ROI.append(N_ROI_local)  

    N_ROI = np.asarray(N_ROI)     #N_ROI[i][j] = Number of Neutrons of Interest measured with ith Amplitude jth Frequencs
    print(np.shape(N_ROI))
    print(N_ROI)
    x = np.linspace(fMin, fMax, noOfFrequencies)



    

    for i in WantedAmplitudes:
        plt.errorbar(x, N_ROI[int(i/Amplitude_Steps - Starting_Amplitude/Amplitude_Steps)], np.sqrt(N_ROI[int(i/Amplitude_Steps - Starting_Amplitude/Amplitude_Steps)]), fmt='.', label = f'Amplitude = {i}mV pp')

    plt.grid()

    plt.xlabel('Frequency $f$ [Hz]')
    plt.ylabel('UCN-counts $N$')
    plt.legend()
    plt.savefig( f'pf2data/Plots/V_and_A_Sweep_run{run}_Amplitudes_{WantedAmplitudes}.png', dpi = 300 )
    plt.show()



run = 177
WantedAmplitudes = [85,95,105]
noOfFiles = 3

Minterval_in_s = 60


noOfAmplitudes = 3
Amplitude_Steps = 10
Starting_Amplitude = 85

fMin = 10000
fMax = 20000
noOfFrequencies = 11

main()

"""
data = np.loadtxt('pf2data/DAQ/run_17/FILTERED/DataF_run_17.CSV', skiprows=1, delimiter=';', usecols=(1,2,3,4))

ch = data[:,0]
time = data[:,1]/1e12
ph = data[:,2]


print(len(splitTimes(time, ch, 50)))


#Dets
d = []
chopper = []
td = []

for i in range(len(ch)):
    if ch[i] == det:
        if time[i]<1000:
            d.append(ph[i])
            td.append(time[i])
    elif ch[i] == 0:
        if ph[i] > 0:
            if time[i]<1000:
                chopper.append(time[i])

d = np.asarray(d)
chopper = np.asarray(chopper)+0.02182 #Chopper offset 21.82ms
td = np.asarray(td)

#print(tmin, tmax, chopper[0], chopper[1])

#np.savez('spectrumTest.npz', array1=d, array2=chopper, array3=td)


#Add bins
countsum = getCounts(1.2,data)


t0 = np.linspace(chopper[0], chopper[1], size)
tp = [0]

for i in range(1, len(t0)):
    tp.append(t0[i]-t0[0])
tp = np.asarray(tp)


#print(len(vp[1:-1]), len(v_corr_countsum[1:]))
#print(v_corr_countsum)

#Plot
#plt.scatter(0.631/tp, countsum*(tp**2/0.631)/sum(countsum), color='blue', marker='.')
#plt.scatter(tp, countsum/sum(countsum), color='blue', marker='.')
tp = np.linspace(tmin, tmax, size)
vp = 0.631/tp
vp = np.asarray(vp)
plt.errorbar(vp[1:-1], countsum[:-1],np.sqrt(countsum[:-1]), fmt='.', color = 'blue')

#integr = np.trapz(countsum/sum(countsum), x=np.linspace(5.5,7,20))


plt.xlabel('UCN-velocity $v_{UCN}$ [m/s]')
plt.ylabel('Normalised velocity-distribution [a.u.]')

plt.grid('both')

#plt.xlim([0,25])

plt.show()
"""