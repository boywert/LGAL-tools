from mass_fn import *
from globalconf import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
from matplotlib import gridspec
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
    fig = plt.figure(figsize=(8, 6)) 
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 3, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    plt.subplots_adjust(hspace = 0)
    for i in range(len(model_names)):
        index = model_names[i]
        xfrac = numpy.loadtxt(tau_folder+"/"+model_names[i]+".log")
        if i== 0:
            ref = xfrac[:,2]
        ax0.plot(xfrac[:,0], xfrac[:,2]/xfrac[:,1],color=model_plot_colors[i],linestyle=model_plot_patterns[i])
        ax1.plot(xfrac[:,0],xfrac[:,2],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
        ax2.plot(xfrac[:,0], (xfrac[:,2]-ref)/ref*100,color=model_plot_colors[i],linestyle=model_plot_patterns[i])
    leg = ax1.legend(loc="upper right", handlelength = 7,ncol=1, fancybox=True, prop={'size':12})
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_alpha(0)
    ax2.set_xlabel(r"redshift")
    ax0.set_ylabel(r"$\langle x^{\mathrm{m}}_{\mathrm{HII}}\rangle/\langle x^{\mathrm{v}}_{\mathrm{HII}}\rangle$")
    ax1.set_ylabel(r"$\langle x^{\mathrm{m}}_{\mathrm{HII}}\rangle$")
    ax2.set_ylabel(r"$\% \rm relative~residual$")
    ax1.set_xlim([6,20])
    ax0.set_xlim([6,20])
    ax2.set_xlim([6,20])
    ax2.set_ylim([-8,4])
    ax0.set_axisbelow(True)
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    ax0.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax0.yaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax1.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax1.yaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax2.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax2.yaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax0.xaxis.set_ticklabels([])
    ax1.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticklabels([r"$-8$",r"$-6$",r"$-4$",r"$-2$",r"$0$",r"$2$",""])
    #ax.set_xscale("log")
    ax1.set_ylim([0.0,1.0])
    ax0.set_ylim([1.0,3.0])
    fig.savefig("xfrachist.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)
    
def main():
    plot_xfrac()

if __name__=="__main__":
    main()
