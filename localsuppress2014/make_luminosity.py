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
import uv_luminosity
sys.path.append("../python/")
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
    filter['Sfr'] = True
    filter['HaloM_Crit200'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)


h0 = 0.7
gadgetmass = 1.e10
model_names = ["okamoto","noreionization","patchy_I"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/n306/sams/43000.00/inputs/LGalaxyStruct.py"]

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/n306/sams/43000.00/"]
model_plot_patterns = ['r--','g--','b--']


pylab.rc('text', usetex=True)
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
#z = zlist[int(sys.argv[1])].strip()



z = "6.98"
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


import uv_luminosity
fig = plt.figure()
ax = fig.add_subplot(111)
uv_luminosity.add_obs_uv_z7("../../codes/47Mpc/observed_UVL/",ax)

for i in range(len(model_names)):
    index = model_names[i]
    if not index in gal:
        (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

    logf = -2.5*numpy.log10(gal[index]["Sfr"])
    a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
    x = a[1][0:len(a[1])-1]+0.25-18
    y = a[0]/47.**3/0.5
    ax.plot(x,y,label=model_labels[i],model_plot_patterns[i])
    
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
ax.set_xlabel(r"M1600 - 5log(h)")
ax.set_ylabel(r"numbers \mathrm{Mpc^{-3}} Mag$^-1$")
ax.set_yscale("log")
fig.savefig("uv_l_z7.pdf",bbox_inches='tight',pad_inches=0)

