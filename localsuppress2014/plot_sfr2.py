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
rank = "0"
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False    
    filter['Sfr'] = True
    filter['BulgeMass'] = True
    filter['DiskMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

#filter model
filter_tmp = []
dt_tmp = []
model_names_tmp = []
struct_file_tmp = []
model_labels_tmp = []
model_paths_tmp = []
for i in range(len(use_model)):
    if use_model[i]:
        filter_tmp.append(filter[i])
        dt_tmp.append(dt[i])
        model_names_tmp.append(model_names[i])
        struct_file_tmp.append(struct_file[i])
        model_labels_tmp.append(model_labels[i])
        model_paths_tmp.append(model_paths[i])

filter = filter_tmp
dt = dt_tmp
model_names = model_names_tmp
struct_file = struct_file_tmp
model_labels = model_labels_tmp
model_paths = model_paths_tmp       



pylab.rc('text', usetex=True)
pylab.rc('lines', linewidth=2)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8
zlist = open(zlistfile,"r").readlines()


def plot_smf_z8(ax):
    z = "7.96"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
   
    sfr_x = {}
    sfr_y = {}
   
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        (sfr_x[index],sfr_y[index]) = sfr_density_fn(gal[index],mass_min=10**-2,mass_max=10.**2,nbins=20)

    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])

    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_xlim([-2,2])
    ax.set_ylim([1e-5,10])
    labels = ["",r"$-1.5$",r"$-1.0$",r"$-0.5$",r"$0$",r"$0.5$",r"$1.0$",r"$1.5$",r"$2.0$"]
    ax.xaxis.set_ticklabels(labels)
    #ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")     
    ax.text(0.9, 0.9, 'z = 8',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=15)
    ax.yaxis.set_ticklabels([])
    
def plot_smf_z7(ax):
    z = "6.98"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
   
    smf_x = {}
    smf_y = {}
   
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        (smf_x[index],smf_y[index]) = stellar_mass_fn(gal[index],mass_min=1.e4,mass_max=1e11,nbins=40)

    add_observations.add_obs_smf_z7("observations/SMF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(smf_x[index],smf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    #leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    #leg.get_frame().set_linewidth(0)
    #ax.set_xlabel(r"$\mathrm{\log_{10}[m_*/M_\odot]}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_xlim([1.e-2,1e2])
    ax.set_xlim([1e-5,10])
    ax.set_yscale("log")     
    ax.text(0.9, 0.9, 'z = 7',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=15)

    
def plot_smf_z6(ax):
    z = "6.00"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    sfr_x = {}
    sfr_y = {}

    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
            (sfr_x[index],sfr_y[index]) = sfr_density_fn(gal[index],mass_min=10**-2,mass_max=10.**2,nbins=20)
            
    add_observations.add_obs_sfr_z6("observations/SFR/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='lower left', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    ax.set_xlim([-2,2])
    ax.set_ylim([1e-5,10])
    ax.text(0.9, 0.9, 'z = 6',
            verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontsize=15)


    
def main():
    fig = plt.figure(figsize=(16, 6))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    #ax3 = fig.add_subplot(3,1,3)
    plt.subplots_adjust(wspace = 0)
    plot_smf_z6(ax1)
    plot_smf_z7(ax2)
    #plot_smf_z8(ax2)
    fig.savefig("SFR67.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)

if __name__=="__main__":
    main()
