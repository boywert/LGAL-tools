from mass_fn import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct

sys.path.append("../python/")
import read_lgal_advance as read_lgal

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


pylab.rc('text', usetex=True)
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
#z = zlist[int(sys.argv[1])].strip()
z = "6.09"
file_prefix = "SA_z"+z
firstfile = 0
lastfile = 127
config = {}

h0 = 0.7
gadgetmass = 1.e10
model_names = ["okamoto","noreionization","patchy_I"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/L-Galaxy_with_RT/sams/43000.00/inputs/LGalaxyStruct.py"]

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/L-Galaxy_with_RT/sams/43000.00/"]
model_plot_patterns = ['r--','g--','b--']

try:
    gal
except NameError:
    gal = {}
    nTrees = {}
    nGals = {}
    nTreeGals = {}
    star = {}
    totsfr = {}
    sfr = {}

for i in range(len(model_names)):
    index = model_names[i]
    if not index in gal:
        (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)




