import matplotlib.cm as cm
from mass_fn import *
from globalconf import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("touch LGalaxyStruct.py")
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
    filter['StellarMass'] = True
    filter['Rvir'] = True
    filter['Mvir'] = True
    filter['Type'] = True
    filter['Pos'] = True
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
    count = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    
   
    for i in range(len(model_names)):
        index = model_names[i]
        print index
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile[i],lastfile[i],filter[i],dt[i],1)
        #gal[index] = gal[index][numpy.where(gal[index]['Pos'][:,2] >51.75275192)]
        #gal[index] = gal[index][numpy.where(gal[index]['Pos'][:,2] <51.80275192)]
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        sc = ax.scatter(gal[index]['Pos'][:,0],gal[index]['Pos'][:,1],c=numpy.log10(gal[index]['StellarMass']*1e10),cmap='Reds')
        #cbar = fig.colorbar(sc)
        ax.set_xlim([47,54])
        ax.set_ylim([43,49])
        cbar = plt.colorbar.ColorbarBase(sc, cmap='Reds',
                       norm=plt.colors.Normalize(vmin=4, vmax=10))
        # leg = ax.legend(loc='upper left', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
        # leg.get_frame().set_linewidth(0)
        fig.savefig("xy_plane_"+str(i)+".png",bbox_inches='tight',pad_inches=0.1)
        plt.close(fig)   


    
    
def main():
    plot_smf()


if __name__=="__main__":
    main()
    
