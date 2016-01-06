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

def plot_xfrac():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        xfrac = numpy.loadtxt(tau_folder+"/"+model_names[i]+"_"+str(model_fesc[i])+".log")
        ax.plot(xfrac[:,0],xfrac[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])

    leg = ax.legend(loc="lower left", handlelength = 7,ncol=1, fancybox=True, prop={'size':12})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"$\langle x^{\mathrm{v}}_{\mathrm{HII}}\rangle$")
    ax.set_xlim([6,20])
    #ax.set_xscale("log")
    ax.set_ylim([0.0,1.0])
    fig.savefig("xfrachist.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)
    
def main():
    plot_xfrac()

if __name__=="__main__":
    main()
