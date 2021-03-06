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
    x = numpy.arange(0.,20.,0.1)
    y = tau_e*numpy.ones(len(x))
    yu = y+delta_tau_e
    yl = y-delta_tau_e
    ax.plot(x,y,"k-",linewidth=1)
    ax.fill_between(x, yu, yl, alpha=0.25, facecolor="#606060",linewidth=0)
    ax.text(7.5, 0.069, 'PLANCK 2015', fontsize=12)
    for i in range(len(model_names)):
        index = model_names[i]
        tau = numpy.loadtxt(tau_folder+"/"+model_names[i]+".tau")
        print tau[0,0]
        if tau[0,0] == 1:
            tau[:,0] = tau[:,0] - 1.
        ax.plot(tau[:,0],tau[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])

    leg = ax.legend(loc=4, handlelength = 10,ncol=1, fancybox=True, prop={'size':12})
    ax.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax.yaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax.set_axisbelow(True)
    leg.get_frame().set_alpha(0)
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"redshift")
    ax.set_ylabel(r"$\tau_e$")
    ax.set_xlim([6,18])
    ax.set_ylim([0.02,0.085])
    fig.savefig("taue.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)
    
def main():
    plot_tau()

if __name__=="__main__":
    main()
