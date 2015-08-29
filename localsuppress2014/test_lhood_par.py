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
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

firstfile = 0
lastfile = 127
zlistfile = "/scratch/01937/cs390/data//snap_z.txt"
config = {}

h0 = 0.7
gadgetmass = 1.e10

model_names = []
struct_file = []
model_labels = []
model_paths = []
use_model = []

for i in range(20):
    model_names.append("input_"+str(i))
    struct_file.append("/scratch/01937/cs390/test_run/python/LGalaxyStruct.py")
    model_labels.append("input_"+str(i))
    model_paths.append("/scratch/01937/cs390/data/outputs/"+str(i)+"/")
    use_model.append(True)

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

zlist = open(zlistfile,"r").readlines()
offset = 18.0
    
def plot_uv_z6():
    z = "6.00"
    file_prefix = "SA_z"+z
    firstfile = 0
    lastfile = 127
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
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        (sfr_x[index],sfr_y[index]) =  sfr_density_fn(gal[index],mass_min=10.**-0.5,mass_max=1.e3,nbins=10)
 
        

    # SFR
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_sfr_z6("observations/SFR/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("mix_sfr_z6.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)

def main():
     plot_uv_z6()

if __name__=="__main__":
    main()
