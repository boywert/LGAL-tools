from mass_fn import *
from globalconf import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct
import add_observations
sys.path.append("../python/")
import read_lgal_advance as read_lgal
import timeit
import random
rank = "0"
SEC_PER_YEAR = 3600*24*365.25
Msun2kg = 1.989e30
h_mass = 1.6737237e-27 #kg
ranki = str(random.randint(0,1000000))
os.system("mkdir -p ../tmp/"+ranki)
pylab.rc('text', usetex=True)
pylab.rc('lines', linewidth=2)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8

def plot_tau():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        tau = numpy.loadtxt(tau_folder+"/"+model_names[i]+".tau")
        if tau[0,0] == "1":
            tau[:,0] = tau[:,0] + 1.
        ax.plot(tau[:,0],tau[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])

    leg = ax.legend(loc=4, handlelength = 10,ncol=1, fancybox=True, prop={'size':12})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"$\tau_e$")
    ax.set_xlim([6,18])
    fig.savefig("taue.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)
    
def main():
    plot_tau()

if __name__=="__main__":
    main()
