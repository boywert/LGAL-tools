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
import xi as CF
import timeit
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def loadfilter(structfile):
    os.system("mkdir -p ../tmp/"+str(rank))
    sys.path.insert(0,"../tmp/"+str(rank))
    os.system("cp "+structfile+" ../tmp/"+str(rank)+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+str(rank)+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['Type'] = True
    filter['Pos'] = True
    filter['Mag'] = True
    filter['Mvir'] = True
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
    xi = {}
        
    gal = {}
    nTrees = {}
    nGals = {}
    nTreeGals = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            if rank == 0:
                (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
    comm.Barrier()
    slot = "Mag"
    mlist = numpy.arange(-18,0,0.5)
    for m in mlist:
        mag  = m
        if rank == 0:
            print "mag",mag
        mag1 = mag+0.5
        for i in range(len(model_names)):
            index = model_names[i]
            if rank == 0:
                #data = gal[index][numpy.where((numpy.log10(gal[index][slot]*1e10)>mag) & (numpy.log10(gal[index][slot]*1e10)<mag1))]["Pos"]
                data = gal[index][numpy.where((gal[index]["Mag"][:,5]>mag) & (gal[index]["Mag"][:,5]<mag1))]["Pos"]
            else:
                data = None
            data = comm.bcast(data,root=0)
            (r,xi[index]) = CF.calNN(data,47.0)            
        if rank == 0:
            print "plotting figure"
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for i in range(len(model_names)):
                index = model_names[i]
                print "adding",model_labels[i]
                ax.plot(r[1:],xi[index][1:]-1.,model_plot_patterns[i],label=model_labels[i])
                leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
            leg.get_frame().set_linewidth(0)
            ax.set_xlabel(r"$r[h^{-1}Mpc]$")
            ax.set_ylabel(r"$\xi(r)$")
            ax.set_yscale("log")
            ax.set_xscale("log")
            print "saving fig",slot+"_"+str(abs(mag))+"_xi"+str(z)+".pdf"
            fig.savefig(slot+"_"+str(abs(mag))+"_xi"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
            print "done"
def main():
    plot_xi("6.00")
    plot_xi("7.96")

if __name__=="__main__":
    main()
