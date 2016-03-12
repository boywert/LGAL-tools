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
#import add_observations
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
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
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
#zlist = open(zlistfile,"r").readlines()


def plot_smf():
    z = "0.00"
    file_prefix = "SA_z"+z
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
        print index
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        (smf_x[index],smf_y[index]) = stellar_mass_fn(gal[index],mass_min=1.e2,mass_max=1e9,nbins=40)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(smf_x[index],smf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    ax.set_xlabel(r"$\mathrm{\log_{10}[m_*/M_\odot]}$")
    ax.set_ylim([1.e-5,1e1])
    ax.set_xlim([2,9])
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")     
    fig.savefig("SMF.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)    
    
def main():
    plot_smf()


if __name__=="__main__":
    main()
