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
    filter['StellarMass'] = True
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
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile[i],lastfile[i],filter[i],dt[i],1)

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(model_names)):
        index = model_names[i]
	print index, len(gal[index])
        pgal = numpy.where(gal[index]["StellarMass"]>0.)[0]
	print len(pgal),len(gal[index]["StellarMass"])
        ax.scatter(numpy.log10(1e10*gal[index]["Mvir"][pgal]/hubble_h),numpy.log10(1e10*gal[index]["StellarMass"][pgal]/hubble_h),s=3,color=model_plot_colors[i],label=model_labels[i])
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    x = numpy.arange(xmin,xmax,0.1)
    A = numpy.log10(1.0e2/4.0e10)/numpy.log10(2.0e6/2.0e12)
    B = numpy.log10(1.0e2)-A*numpy.log10(2.0e6)
    y = A*x+B
    ax.plot(x,y,'k--',linewidth=1,label="Behroozi+ (2013)")
    A = numpy.log10(1.0e2/4.0e10)/numpy.log10(5.0e7/2.0e12)
    B = numpy.log10(1.0e2)-A*numpy.log10(5.0e7)
    y = A*x+B
    ax.plot(x,y,'k.',linewidth=1,label="Garrison-Kimmel+ (2014)")
    ax.set_ylabel(r"$\mathrm{\log_{10}[M_*/M_\odot]}$")
    ax.set_xlabel(r"$\mathrm{\log_{10}[M_{DM}/M_\odot]}$")
    leg = ax.legend(loc='upper left', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    fig.savefig("SMHM.png",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)    
    
def main():
    plot_smf()


if __name__=="__main__":
    main()
