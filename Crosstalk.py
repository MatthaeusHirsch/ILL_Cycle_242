import numpy as np
import matplotlib.pyplot as plt
import uncertainties as u



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



def getCounts(run):

    td = []
    tc = []

    if noOfFiles > 1 :
        for i in range(noOfFiles):
            pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}_{i}.CSV'
            result = read_files(pfad)
            td.extend(result[0])
    else: 
        pfad = f'pf2data/DAQ/run_{run}/FILTERED/DataF_run_{run}.CSV'
        result = read_files(pfad)
        td.extend(result[0])


    td = np.asarray(td)


    return len(td)/max(td), np.sqrt(len(td))/max(td)


def plot_runs(runs, xlim, lineplot = False):
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # Add more colors if needed

    O = getCounts(180)
    plt.errorbar(0, O[0], O[1] , fmt='x', label = 'N000' , markersize = 10)
    I = getCounts(181)
    plt.errorbar(0, I[0], I[1] , fmt='x', label = 'N010' , markersize = 10)
    II = getCounts(182)
    plt.errorbar(0, II[0], II[1] , fmt='x', label = 'N001' , markersize = 10)
    III = getCounts(183)
    plt.errorbar(0, III[0], III[1] , fmt='x', label = 'N011' , markersize = 10)
    plt.grid()
    plt.legend()
    plt.ylabel('UCN count rate [s$^{-1}$]')
    plt.xlim(xlim)
    plt.ylim(O[0]*0.99,O[0]*1.01)
    plt.xticks([])

    print(O)
    print(I)
    print(II)
    print(III)

    
    # Save the figure with the runs in the filename
    
    plt.savefig(f'pf2data/Plots/Crosstalk_runs{"_".join(map(str, runs))}_plotted_on_1percent_deviation_around_N000.png', dpi=300)
    
    plt.show()









#run the code

noOfFiles = 1
plot_runs([180, 181, 182, 183], [-1,1], lineplot=True) # runs to plot and velocity limits to plot in
