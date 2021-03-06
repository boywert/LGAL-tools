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
    filter['Len'] = True
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
        rangen = (0,1.)
        bins = 20
        step = (rangen[1]-rangen[0])/bins
        count[index] = numpy.zeros(bins,dtype=numpy.int64)
        total = long(0)
        firstgal = numpy.where(gal[index]["Type"] == 0)[0]
        for ii in range(len(firstgal)-1):
            for j in range(firstgal[ii+1]-firstgal[ii]+1):
                this_gal = firstgal[ii]+j
                distance = numpy.sqrt((gal[index][this_gal]['Pos'][0] - gal[index][firstgal[ii]]['Pos'][0])**2.+(gal[index][this_gal]['Pos'][1] - gal[index][firstgal[ii]]['Pos'][1])**2.+(gal[index][this_gal]['Pos'][2] - gal[index][firstgal[ii]]['Pos'][2])**2.)/(1.+float(z))
                if ((distance < gal[index][firstgal[ii]]['Rvir']) & (gal[index][this_gal]["Type"] > 0) & (gal[index][this_gal]["Len"] >= 20) ):
                    slot = int(distance/gal[index][firstgal[ii]]['Rvir']/step)
                    count[index][slot]+=1

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(model_names)):
        index = model_names[i]
        r = numpy.arange(0,1,1./bins)+1./bins/2.
        ax.plot(r,count[index].astype(numpy.float64)/numpy.sum(count[index])/step,color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    ax.set_xlabel(r"$r/R_{200c}$")
    ax.set_ylabel(r"$\mathrm{pdf}/(\Delta r/R_{200c})$")
    leg = ax.legend(loc='upper left', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    fig.savefig("r_dist.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)   

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(model_names)):
        index = model_names[i]
        r = numpy.arange(0,1,1./bins)+1./bins/2.
        ax.plot(r,count[index].astype(numpy.float64)/step,color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='upper left', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$r/R_{200c}$")
    ax.set_ylabel(r"$N/(\Delta r/R_{200c})$")
    
    fig.savefig("r_dist_total.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)   
    
    
def main():
    plot_smf()


if __name__=="__main__":
    main()
    
