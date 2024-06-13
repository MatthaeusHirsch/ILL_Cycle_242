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



def getCounts(M, time, chopper):

    return [sum(time), np.sqrt(sum(time))]


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


def plot(x, y):
    plt.errorbar(x, y[:,0], y[:,1], fmt='.', color='blue')

    plt.grid()

    plt.xlabel('Frequency $f$ [Hz]')
    plt.ylabel('UCN-counts inside the ROI $N$')

    plt.show()


def main():
    run = 24
    noOfFiles = 4
    interval = 50
    noOfInterval = 40
    M = 1.2
    fMin = 10000
    fMax = 30000

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

    tdSplit, tcSplit = splitTimes(tdGes, tcGes, interval, noOfInterval)
    
    N_ROI = [] #Counts inside region of interest in the i-th time-interval
    for i in range(len(tdSplit)):
        N_ROI.append(getCounts(M, tdSplit[i], tcSplit[i]))

    N_ROI = np.asarray(N_ROI)

    x = np.linspace(fMin, fMax, len(N_ROI))

    plot(x, N_ROI)





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