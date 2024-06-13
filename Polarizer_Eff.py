import numpy as np
import matplotlib.pyplot as plt
import uncertainties as u
from uncertainties import umath as um




det = 1

#Calculate v from M in T
def calc_ts(M):
    Emax = (210+M*60)*1e-9
    Emin = (210-60*M)*1e-9

    vmax = np.sqrt(2*Emax/(939e6/3e8**2))
    vmin = np.sqrt(2*Emin/(939e6/3e8**2))

    tmax = 0.631/vmin
    tmin = 0.631/vmax
    return (tmax , tmin)


def getCounts(M, time, chopper):
    tmax , tmin = calc_ts(M)
    size = 2
    #Add bins
    countsum = np.zeros(size)
    for i in range(len(chopper)-1):
        n, bins = np.histogram(time, bins=size, range=(chopper[i]+tmin, chopper[i]+tmax))

        countsum = countsum + n
    return [sum(countsum), np.sqrt(sum(countsum))]


def getCountsbytimes(time, chopper, tmin, tmax):
    size = 100
    #Add bins
    countsum = np.zeros(size)
    for i in range(len(chopper)-1):
        n, bins = np.histogram(time, bins=size, range=(chopper[i]+tmin, chopper[i]+tmax))

        countsum = countsum + n
    return [sum(countsum), np.sqrt(sum(countsum))]


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

    chopper = np.asarray(chopper)+0.021 #Chopper offset 21.82ms
    td = np.asarray(td)

    return (td, chopper)








def getRates_with_M(run, M):
    #with automated time analysis
    pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'
    data = read_files(pfad)
    
    res = getCounts(M, data[0], data[1])
    return u.ufloat(res[0]/data[0].max(), np.sqrt(res[0])/data[0].max())
    
def getRates_with_velocity(run, vmin , vmax):
    #with manual time analysis
    pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'
    data = read_files(pfad)
    tmin = 0.631/vmax
    tmax = 0.631/vmin
    res = getCountsbytimes(data[0], data[1], tmin , tmax)
    #res = getCountsbytimes(data[0], data[1], 0,0.4)
    return u.ufloat(res[0]/data[0].max(), np.sqrt(res[0])/data[0].max()) / max(data[1])

def getRates_full_Spectrum(run, M):
    pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'
    data = read_files(pfad)
    time = data[0]
    return u.ufloat(len(time)/time.max(), np.sqrt(len(time))/time.max())
    




def calc_PolEff_i(X_i , X_i2 , X_0):
    return um.sqrt(np.divide(X_i*X_i2 , X_0))



M = 2
T12_0 = getRates_with_M(164, M)
T12_1 = getRates_with_M(165, M)
print(f"T12_0 = {T12_0}")
print(f"T12_1 = {T12_1}")

T32_0 = getRates_with_M(168, M)
T32_1 = getRates_with_M(169, M)
print(f"T32_0 = {T32_0}")
print(f"T32_1 = {T32_1}")

T31_0 = getRates_with_M(172, M)
T31_1 = getRates_with_M(173, M)
print(f"T31_0 = {T31_0}")
print(f"T31_1 = {T31_1}")

T34_0 = getRates_with_M(174, M)
T34_1 = getRates_with_M(175, M)
print(f"T34_0 = {T34_0}")
print(f"T34_1 = {T34_1}")

X12 = np.divide(T12_0-T12_1 , T12_0+T12_1)
X32 = np.divide(T32_0-T32_1 , T32_0+T32_1)
X31 = np.divide(T31_0-T31_1 , T31_0+T31_1)
X34 = np.divide(T34_0-T34_1 , T34_0+T34_1)

p1 = calc_PolEff_i(X12 , X31 , X32)
p2 = calc_PolEff_i(X12 , X32 , X31)
p3 = calc_PolEff_i(X31 , X32 , X12)
p4 = np.divide(X34 , p3)

print("Druuuummmrooooolllll \n")
print(f"Polarizer Eff. 1 = {p1}")
print("\n")
print(f"Polarizer Eff. 2 = {p2}")
print("\n")
print(f"Polarizer Eff. 3 = {p3}")
print("\n")
print(f"Polarizer Eff. 4 = {p4}")



