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

    for i in range(len(model_names)):
        index = model_names[i]
        data = numpy.loadtxt(folder+"/"+index)
        x = data[:,0]
        dx = numpy.zeros(len(x))
        for i in range(len(x)-1):
            dx = x[i] - x[i+1]
        dx[len(x)-1] = x[len(x)-1]
        ax.plot(x,data[:,3]/numpy.sum(data[:,3])/,color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])

    if pos != "l":
        ax.yaxis.set_ticklabels([])
        #labels = ["",r"$8.5$",r"$9.0$",r"$9.5$",r"$10.0$",r"$10.5$",r"$11.0$",r"$11.5$",r"$12.0$"]
        #ax.xaxis.set_ticklabels(labels)
    ax.set_xscale("log")
    #ax.set_yscale("log")
    ax.set_xlabel(r"$\log_{10}(M_{\mathrm{200c}}/\mathrm{M_\odot})$")
    if pos == "l":
        ax.set_ylabel(r"$\log_{10}(m_*/\mathrm{M_\odot})$")
        leg = ax.legend(loc="upper right", handlelength = 6,ncol=1, fancybox=True, prop={'size':10})
        leg.get_frame().set_linewidth(0)

    ax.text(0.1, 0.9, r'$\langle x^{\mathrm{m}}_{\mathrm{HII}}\rangle = %s$' %(folder),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=14)


    
def main():
    fig = plt.figure(figsize=(24, 6))
    plt.subplots_adjust(wspace = 0)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    fig.canvas.draw()
    plot_size(ax1,"0.3","l")
    plot_size(ax2,"0.5","c")
    plot_size(ax3,"0.7","r")
    fig.savefig("sizedist.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
