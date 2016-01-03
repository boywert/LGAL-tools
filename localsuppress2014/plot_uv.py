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
    filter['MagDust'] = True
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


def plot_uv_z8(ax):
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
    luvlf_x = {}
    luvlf_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-13,nbins=24)

    
    add_observations.add_obs_uv_z8("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(luvlf_x[index],luvlf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label="")
    leg = ax.legend(loc=4, handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1500 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    ax.set_xlim([-20.5,-16.])
    ax.set_ylim([1e-5,1e-1])
    ax.text(0.1, 0.9, 'z = 8',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=15)



def plot_uv_z7(ax):
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
   
    luvlf_x = {}
    luvlf_y = {}
   
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-15,nbins=20)
     
    # UVLF
    add_observations.add_obs_uv_z7("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(luvlf_x[index],luvlf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label="")
    leg = ax.legend(loc=4, handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    #ax.set_xlabel(r"M1500 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    ax.set_xlim([-21,-16])
    ax.text(0.1, 0.9, 'z = 7',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=15)

    
def plot_uv_z6(ax):
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

    luvlf_x = {}
    luvlf_y = {}

    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-15,nbins=20)
         
    # UVLF

    add_observations.add_obs_uv_z6("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(luvlf_x[index],luvlf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc=4, handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    #ax.set_xlabel(r"M1500 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    ax.set_xlim([-22.,-15.5])
    ax.set_ylim([1e-6,1e-1])
    ax.text(0.1, 0.9, 'z = 6',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=15)


    
def main():
    fig = plt.figure(figsize=(8, 18))
    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2)
    ax3 = fig.add_subplot(3,1,3)
    plt.subplots_adjust(hspace = .4)
    plot_uv_z6(ax1)
    plot_uv_z7(ax2)
    plot_uv_z8(ax3)
    fig.savefig("UVLF678.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)

if __name__=="__main__":
    main()
