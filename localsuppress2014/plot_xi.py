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
pylab.rc('lines', linewidth=1.5)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8
folder = "xi/"
def plot_size(ax,m,z,pos):    
    for i in range(len(model_names)):
        index = model_names[i]
        data = numpy.loadtxt(folder+"/StellarMass_"+index+"_"+m+"_"+z+".txt")
        ax.plot(data[:,0],data[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])


    ax.set_xscale("log")
    ax.set_xlim([1e-1,10])
    ax.set_ylim([1e-1,1e3])
    ax.set_yscale("log")
    
    
    if pos == 1:
        labels = ["","", r"$1$",r"$10$",r"$100$",r"$1000$"]
        ax.yaxis.set_ticklabels(labels)
    if pos == 3:
        labels = ["",r"$0.1$",r"$1$",r"$10$",r"$100$",r"$1000$"]
        ax.yaxis.set_ticklabels(labels)
    if pos == 3:
        labels = ["",r"$0.1$",r"$1$",r"$10$"]
        ax.xaxis.set_ticklabels(labels)
    if pos == 4:
        labels = ["","",r"$1$",r"$10$"]
        ax.xaxis.set_ticklabels(labels)
    if pos == 1:
        leg = ax.legend(loc="upper right", handlelength = 7,ncol=1, fancybox=True, prop={'size':13})
        leg.get_frame().set_linewidth(0)
    if (pos%2==1):
        ax.set_ylabel(r"$\xi(\mathrm{r})$")
    if (pos== 5) | (pos==6):
        ax.set_xlabel(r"$r[h^{-1}\mathrm{Mpc}]$")
    if (pos%2 == 0):
        ax.yaxis.set_ticklabels([])
    if (pos <= 4 ):
        ax.xaxis.set_ticklabels([])
    ax.text(0.95, 0.5, r"$10^{%2.0f} < m_*/\mathrm{M_\odot} < 10^{%2.0f}$" % (float(m),float(m)+1.),
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=18)
    ax.text(0.1, 0.1, "z = %d" % (int(float(z)+0.5)),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=18)

    
def main():
    fig = plt.figure(figsize=(16, 18))
    plt.subplots_adjust(wspace = 0)
    plt.subplots_adjust(hspace = 0)
    ax1 = fig.add_subplot(321)
    ax2 = fig.add_subplot(322)
    ax3 = fig.add_subplot(323)
    ax4 = fig.add_subplot(324)
    ax5 = fig.add_subplot(325)
    ax6 = fig.add_subplot(326)
    fig.canvas.draw()
    plot_size(ax1, "4.0","6.00", 1)
    plot_size(ax2, "4.0","9.03", 2)
    plot_size(ax3, "6.0","6.00", 3)
    plot_size(ax4, "6.0","9.03", 4)
    plot_size(ax5, "7.0","6.00", 5)
    plot_size(ax6, "7.0","9.03", 6)
    fig.savefig("xi_.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
