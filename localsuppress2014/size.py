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

def plot_size(ax,folder,pos):    
    print folder
    for i in range(len(model_names)):
        index = model_names[i]
        print index
        data = numpy.loadtxt(folder+"/"+index)
        x = data[:,0]
        print x
        dx = numpy.zeros(len(x))
        for j in range(len(x)-1):
            dx[j] = x[j] - x[j+1]
        dx[len(x)-1] = x[len(x)-1]
        ax.plot(x,x*data[:,3]/numpy.sum(data[:,3])/dx,color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])


    ax.set_xscale("log")
    ax.set_xlim([1e-1,30])
    ax.set_ylim([1e-2,1])
    ax.set_yscale("log")

    ax.set_ylabel(r"$\mathrm{R ~dP(R)/dR}$")
    if pos == "b":
        ax.set_xlabel(r"$\mathrm{R/Mpc}$")
        
    if pos == "t":
        leg = ax.legend(loc="lower left", handlelength = 5,ncol=1, fancybox=True, prop={'size':14})
        leg.get_frame().set_linewidth(0)

    ax.text(0.05, 0.85, r'$\langle x^{\mathrm{m}}_{\mathrm{HII}}\rangle = %s$' %(folder),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=18)
    if pos != "b":
        ax.xaxis.set_ticklabels([])
        labels = ["",r"$0.1$",r"$1.0$",r"$10^2$",r"$10^3$"]
        ax.yaxis.set_ticklabels(labels)

    
def main():
    fig = plt.figure(figsize=(8, 18))
    plt.subplots_adjust(hspace = 0)
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    fig.canvas.draw()
    plot_size(ax1,"0.3","t")
    plot_size(ax2,"0.5","m")
    plot_size(ax3,"0.7","b")
    fig.savefig("sizedist.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
