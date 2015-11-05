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
sys.path.append("../clustering/")
import read_lgal_advance as read_lgal
import xi
import timeit
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+str(rank))
    os.system("cp "+structfile+" ../tmp/"+str(rank)+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+str(rank)+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['Pos'] = True
    filter['Mag'] = True
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


def plot_xi(z):
    file_prefix = "SA_z"+z
    if rank == 0:
        try:
            gal
        except NameError:
            gal = {}
            nTrees = {}
            nGals = {}
            nTreeGals = {}
            for i in range(len(model_names)):
                index = model_names[i]
                if not index in gal:
                    (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
                gal[index] = gal[index][numpy.where(gal[index]["Mag"][:,5]<-15.)]
    else:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
        (nTrees[index],nGals[index],nTreeGals[index],gal[index])  = (None,None,None,None)
        
    mpi.bcast(gal,root=0)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],mean_logphoton[index],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{NPHOT}$")
    # ax.set_yscale("log")
    # fig.savefig("NPHOTvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
def main():
    plot_xi("6.00")
    #plot_z("7.96")

if __name__=="__main__":
    main()
