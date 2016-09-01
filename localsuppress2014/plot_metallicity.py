from mass_fn import *
import matplotlib
matplotlib.use('pdf') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct
import add_observations
sys.path.append("../python/")
from globalconf import *
import read_lgal_advance as read_lgal
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
    filter['Type'] = True
    filter['Metallicity'] = True
    filter['HaloM_Crit200'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)


dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)



pylab.rc('text', usetex=True)

zlist = open(zlistfile,"r").readlines()

def metallicity_plot(z):
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

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        logf = -2.5*numpy.log10(gal[index]["Sfr"])
        a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
        x = a[1][0:len(a[1])-1]+0.25-offset
        y = a[0]/47.**3/0.5
        ax.plot(x,y,model_plot_patterns[i],label=model_labels[i])
    
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"numbers $\mathrm{Mpc^{-3} Mag^-1}$")
    ax.set_yscale("log")
    fig.savefig("uv_l_z"+z+".pdf",bbox_inches='tight',pad_inches=0)



def main():
    metallicity_plot("6.00")
    
if __name__=="__main__":
    main()
