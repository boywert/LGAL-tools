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
    filter['StellarMass'] = True
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
z3list = open(z3listfile,"r").readlines()

class xfrac:
    grid = 0
    data = 0
    
def read_xfrac(filename):
    f = open(filename,"rb")
    output = xfrac()
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    output.grid = numpy.fromfile(f,numpy.int32,3)
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    output.data = numpy.fromfile(f,numpy.float32,output.grid[0]**3).reshape(( output.grid[0], output.grid[1], output. grid[2]))
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    return output

def plot_xi(snap):
    z = zlist[long(snap)].strip()
    z3 = z3list[long(snap)].strip()
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
                sys.stdout.flush()
                # if "infall" in index:
                #     gal[index]["MagDust"][:,5] = gal[index]["MagDust"][:,5] - 0.87
    comm.Barrier()
    slot = "StellarMass"
    m_i = 6.
    m_f = 8.
    dm = 1.0
    mlist = numpy.arange(m_i,m_f,dm)
    for m in mlist:
        mag  = m
        if rank == 0:
            print "mag",mag
        mag1 = mag+dm
        for i in range(len(model_names)):
            index = model_names[i]
            if rank == 0:
                xfilename = model_xfrac_path[i]+"xfrac3d_"+ z3 +".bin"
                #Xfrac3d = read_xfrac(xfilename)
                data = gal[index][numpy.where((numpy.log10(gal[index][slot]*1e10/hubble_h)>mag) & (numpy.log10(gal[index][slot]*1e10/hubble_h)<mag1))]["Pos"]
                xmask = numpy.ones(len(data),dtype=numpy.float32)
                # for iii in range(len(data)):
                #     iix =long(data[iii]['Pos'][0]/(sim_boxsize/Xfrac3d.grid[0]))%Xfrac3d.grid[0]
                #     iiy =long(data[iii]['Pos'][1]/(sim_boxsize/Xfrac3d.grid[1]))%Xfrac3d.grid[1]
                #     iiz =long(data[iii]['Pos'][2]/(sim_boxsize/Xfrac3d.grid[2]))%Xfrac3d.grid[2]
                #     iblock = iix+iiy*Xfrac3d.grid[0]+iiz*Xfrac3d.grid[0]*Xfrac3d.grid[1]
                #     xmask[iii] = Xfrac3d.data[iblock]
                data = data[xmask > 0.99]
            else:
                data = None
            data = comm.bcast(data,root=0)
            (r,xi[index]) = CF.calNN(data,sim_boxsize)            
        if rank == 0:
            print "plotting figure"
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for i in range(len(model_names)):
                index = model_names[i]
                print "adding",model_labels[i]
                ff = open(slot+"_"+model_names[i]+"_"+str(mag)+"_"+str(z)+".txt","w")
                print "creating",slot+"_"+model_names[i]+"_"+str(mag)+"_"+str(z)+".txt"
                for ii in range(len(r)):
                    print >> ff, r[ii],xi[index][ii]-1.
                ff.close()
                
                ax.plot(r[1:],xi[index][1:]-1.,model_plot_patterns[i],label=model_labels[i])
                leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
            leg.get_frame().set_linewidth(0)
            ax.set_xlabel(r"$\mathrm{r}[h^{-1}\mathrm{Mpc}]$")
            ax.set_ylabel(r"$\xi(\mathrm{r})$")
            ax.set_yscale("log")
            ax.set_xscale("log")
            ax.text(0.1, 0.1, r'M1500 $\in$ $[%2.0f,%2.0f]$' % (mag,mag1),
                verticalalignment='bottom', horizontalalignment='left',
                transform=ax.transAxes, fontsize=14)
            #print "saving fig",slot+"_"+str(abs(mag))+"_xi"+str(z)+".pdf"
            #fig.savefig(slot+"_"+str(abs(mag))+"_xi"+str(z)+".pdf",bbox_inches='tight',pad_inches=0.01)
            print "done"
            plt.close(fig)
def main():
    isnap = sys.argv[1]
    plot_xi(isnap)

if __name__=="__main__":
    main()
