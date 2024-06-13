import numpy as np
import matplotlib.pyplot as plt
import uncertainties as u
from uncertainties import umath as um



M = 2
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


def calcEff(run):


    pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'

    data = read_files(pfad)
    



    #OLD: without time analysis
    '''
    time = data[0]
    return u.ufloat(len(time)/time.max(), np.sqrt(len(time))/time.max())
    '''

    #with automated time analysis
    res = getCounts(M, data[0], data[1])
    return u.ufloat(res[0]/data[0].max(), np.sqrt(res[0])/data[0].max())
    

    #with manual time analysis
    '''
    res = getCountsbytimes(data[0], data[1], 0.08320803680467945 ,0.13156345796695304)
    #res = getCountsbytimes(data[0], data[1], 0,0.4)
    return u.ufloat(res[0]/data[0].max(), np.sqrt(res[0])/data[0].max()) / max(data[1])
    '''


def calc_sfEff_F2(OO,II,OI,IO):
    return(0.5*(1+((II-OI) / (OO-IO))))
def calc_sfEff_F1(OO,II,OI,IO):
    return(0.5*(1+((II-IO) / (OO-OI))))
def calc_sqrtsfEff_F1(OO,II,OI,IO):
    print(type(0.5*(1+((II-IO) / (OO-OI)))))
    return(um.sqrt(0.5*(1+((II-IO) / (OO-OI)))))




''' #Daten vom 3.6.
II = calcEff(34)
OO = calcEff(31)
IO = calcEff(33)
OI = calcEff(32)
'''

#Daten vom 4.6.
II = calcEff(166)
OO = calcEff(164)
IO = calcEff(165)
OI = calcEff(167)


print(f"II = {II}     OO = {OO}     IO = {IO}     OI = {OI}")
print("\n")
print(f"Spinflip-Eff.IO = {calc_sfEff_F1(OO,II,OI,IO)}")
print("\n")
print(f"Spinflip-Eff.OI = {calc_sfEff_F2(OO,II,OI,IO)}")
print("\n")
print(f"sqrt = {calc_sqrtsfEff_F1(OO,II,OI,IO)}")



