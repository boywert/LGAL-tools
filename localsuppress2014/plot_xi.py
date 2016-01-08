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
folder = "xi/"
def plot_size(ax,m,z,pos):    

    for i in range(len(model_names)):
        index = model_names[i]
        data = numpy.loadtxt(folder+"/StellarMass_"+index+"_"+m+"_"+z+".txt")
        ax.plot(data[:,0],data[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])


    ax.set_xscale("log")
    #ax.set_xlim([1e-1,100])
    #ax.set_ylim([1e-2,1])
    ax.set_yscale("log")
    ax.set_xlabel(r"$\xi(\mathrm{r})$")
    if pos%2 == 1:
        ax.set_ylabel(r"$\mathrm{r~[h^{-1}Mpc]}$")
    if pos == "r":
        leg = ax.legend(loc="lower right", handlelength = 7,ncol=1, fancybox=True, prop={'size':11})
        leg.get_frame().set_linewidth(0)

    ax.text(0.95, 0.9, z+m,
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=18)
    if pos != "l":
        ax.yaxis.set_ticklabels([])
        labels = ["",r"$10^{0}$",r"$10^{1}$",r"$10^2$",r"$10^3$"]
        ax.xaxis.set_ticklabels(labels)

    
def main():
    fig = plt.figure()
    plt.subplots_adjust(wspace = 0)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    fig.canvas.draw()
    plot_size(ax1, "6.0","6.00", 1)
    plot_size(ax2, "6.0","9.03", 2)
    plot_size(ax3, "7.0","6.00", 3)
    plot_size(ax4, "7.0","9.03", 4)
    fig.savefig("xi_.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
