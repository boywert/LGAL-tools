from mass_fn import *
from globalconf import *
import globalconf as models
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

def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+ranki)
    os.system("cp "+structfile+" ../tmp/"+ranki+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+ranki+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['Mvir'] = True
    filter['HaloM_Crit200'] = True
    filter['Sfr'] = True
    filter['Type'] = True
    filter['CumulativeSFR'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

def setfilter(models):
    dt = []
    filter = []
    for i in range(len(models.struct_file)):
        (f,t) = loadfilter(models.struct_file[i])
        filter.append(f)
        dt.append(t)
    return dt,filter


# filter_tmp = []
# dt_tmp = []
# model_names_tmp = []
# struct_file_tmp = []
# model_labels_tmp = []
# model_paths_tmp = []
# for i in range(len(use_model)):
#     if use_model[i]:
#         filter_tmp.append(filter[i])
#         dt_tmp.append(dt[i])
#         model_names_tmp.append(model_names[i])
#         struct_file_tmp.append(struct_file[i])
#         model_labels_tmp.append(model_labels[i])
#         model_paths_tmp.append(model_paths[i])

# filter = filter_tmp
# dt = dt_tmp
# model_names = model_names_tmp
# struct_file = struct_file_tmp
# model_labels = models.model_labels_tmp
# models.model_paths = models.model_paths_tmp
def plot_z(z,ax,pos):    
    dt,filter = setfilter(models)
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    cmass_x={}
    cmass_y={}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (cmass_x[index],cmass_y[index]) = integrated_stellar_mass_fn(gal[index],mass_min=1e3,mass_max=1e12,nbins=30)
         
    # UVLF
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(cmass_x[index],cmass_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc=4, handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    #ax.set_xlabel(r"M1500 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    #ax.set_xlim([-22.,-15.5])
    #ax.set_ylim([1e-6,1e-1])
    ax.text(0.1, 0.9, 'z = '+long(z),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=15)

    
def main():
    zlist = open(zlistfile).readlines()
    zi = zlist[long(sys.argv[1])].strip()
    fig = plt.figure(figsize=(8, 18))
    plt.subplots_adjust(hspace = 0.1)
    ax1 = fig.add_subplot(311)
    plot_z("6.00",ax1,"t")
    ax2 = fig.add_subplot(312)
    plot_z("6.98",ax2,"m")
    ax2 = fig.add_subplot(312)
    plot_z("7.96",ax2,"b")
    fig.savefig("Baryons_"+zi+".pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)

if __name__=="__main__":
    main()
