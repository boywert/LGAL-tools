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
        rangen = (0,1.)
        bins = 100
        step = (rangen[1]-rangen[0])/bins
        count = numpy.arange(rangen[0],rangen[1],bins,dtype=numpy.int64)
        firstgal = numpy.where(gal[index]["Type"] == 0)[0]
        for ii in range(len(firstgal)-1):
            for j in range(firstgal[ii+1]-firstgal[ii]+1):
                this_gal = firstgal[ii]+j
                distance = numpy.sqrt((gal[index][this_gal]['Pos'][0] - gal[index][firstgal[ii]]['Pos'][0])**2.+(gal[index][this_gal]['Pos'][1] - gal[index][firstgal[ii]]['Pos'][1])**2.+(gal[index][this_gal]['Pos'][2] - gal[index][firstgal[ii]]['Pos'][2])**2.)/(1.+float(z))
                if (distance < gal[index][firstgal[ii]]['Rvir']):
                    slot = int(distance/gal[index][firstgal[ii]]['Rvir']/step)
                    print distance/gal[index][firstgal[ii]]['Rvir'],slot

    # fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     pgal = numpy.where(gal[index]["StellarMass"]>0.)
    #     ax.scatter(numpy.log10(1e10*gal[index]["Mvir"][pgal]),numpy.log10(1e10*gal[index]["StellarMass"][pgal]))
    # ax.set_ylabel(r"$\mathrm{\log_{10}[h^{-1}M_*/M_\odot]}$")
    # ax.set_xlabel(r"$\mathrm{\log_{10}[h^{-1}M_{DM}/M_\odot]}$")
    
    # fig.savefig("SMHM.png",bbox_inches='tight',pad_inches=0.1)
    # plt.close(fig)   
    
    
def main():
    plot_smf()


if __name__=="__main__":
    main()
    
