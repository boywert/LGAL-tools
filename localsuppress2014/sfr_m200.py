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
    filter['Halo_M200Crit'] = True
    filter['Sfr'] = True
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

zlist = open(zlistfile,"r").readlines()


def plot_uv_z8():
    z = "6.00"
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    sum_SFR = {}
    sub_SFR_sq = {}
    N = {}
    mean_SFR = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        sum_SFR[index] = numpy.histogram(gal[index]["Halo_M200Crit"],bins=50,weigh=gal[index]["Sfr"])
        sum_SFR_sq[index] = numpy.histogram(gal[index]["Halo_M200Crit"],bins=50,weigh=gal[index]["Sfr"]**2)
        N[index] = numpy.histogram(gal[index]["Halo_M200Crit"],bins=50)
        mean_SFR[index] = sum_SFR[:,0]/N[:,0]
        m200c[index] = []
        for i in range(len(sub_SFR[:,0])):
            m200c[index].append(0.5*(sub_SFR[i,1]+sub_SFR[i+1,1]))
        del(gal[index])
        del(nTreeGals[index])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(m200c[index],mean_SFR[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$M_{200c}[M_\odot]$")
    ax.set_ylabel(r"$\mathrm{SFR [M_\odot/year]}$")
    ax.set_yscale("log")
    fig.savefig("SFRvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    
def main():
    #plot_uv_z6()
    #plot_uv_z7()
    plot_uv_z8()

if __name__=="__main__":
    main()
