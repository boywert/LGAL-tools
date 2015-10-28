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
SEC_PER_YEAR = 3600*24*365.25
Msun2kg = 1.989e30
h_mass = 1.6737237e-27 #kg
hubble_h = 0.7
boxsize = 47.0
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['NPhotReion'] = True
    filter['HaloM_Crit200'] = True
    filter['Sfr'] = True
    filter['CumulativeSFR'] = True
    filter["Mvir"] = True
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


def get_uv(z):
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sum_photon = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        gal[index]["NPhotReion"] = gal[index]["NPhotReion"] + numpy.log10(11.6e6*SEC_PER_YEAR)
        
        rangen = (7.5,11.5)
        bins = 40                
        sum_photon[index] = numpy.sum((numpy.float64(10)**gal[index]["NPhotReion"].astype(numpy.float64)-1.),dtype=numpy.float64)
        print sum_photon[index]
        
    del(gal)
    return sum_photon


def main_plot():
    zlist = open(zlistfile,"r").readlines()
    uv_gamma = {}
    z_plot = {}
    for i in range(len(model_names)):
        index = model_names[i]
        uv_gamma[index] = numpy.zeros(len(zlist),dtype=numpy.float64)
        z_plot[index] = numpy.zeros(len(zlist),dtype=numpy.float64)
    for iz in range(len(zlist)):
        z = zlist[iz]
        uv = get_uv(z.strip())
        for i in range(len(model_names)):
            index = model_names[i]
            z_plot[index][iz] = float(z.strip())
            uv_gamma[index][iz] = uv[index]
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    print uv_gamma['oka_infall']
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(z_plot[index],uv_gamma[index]/(boxsize/hubble_h/143.)/1e70,model_plot_patterns[i],label=model_labels[i])
    
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$N_{photon} Ilian plot$")
    ax.set_yscale("log")
    fig.savefig("gammavsz.pdf",bbox_inches='tight',pad_inches=0)
def main():
    main_plot()

if __name__=="__main__":
    main()
