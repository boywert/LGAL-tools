from mass_fn import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct

sys.path.append("../python/")
import read_lgal_advance as read_lgal

def loadfilter(structfile):
    sys.path.insert(0,"../tmp/")
    os.system("cp "+structfile+" ../tmp/LGalaxyStruct.py")
    os.system("rm -f ../tmp/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    filter['Type'] = True
    filter['Sfr'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)



#z = sys.argv[1]
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
z = zlist[sys.argv[1]]

file_prefix = "SA_z"+z
firstfile = 0
lastfile = 127
config = {}
model_names = ["okamoto","noreionization","patchy_I"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_001/sams/5500.00/LGalaxyStruct.py"]

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_001/sams/5500.00/"]
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
    #hotgas = {}
    #coldgas = {}
    #Blackhole = {}
    sfr = {}

for i in range(len(model_names)):
    index = model_names[i]
    if not index in gal:
        (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)


pylab.rc('text', usetex=True)


gal_type0 = {}
gal_type1 = {}
gal_type2 = {}

for i in range(len(model_names)):
    index = model_names[i]
    gal_type0[index] = gal[index][numpy.where(gal[index]["Type"] == 0)[0]]
    gal_type1[index] = gal[index][numpy.where(gal[index]["Type"] == 1)[0]]
    gal_type2[index] = gal[index][numpy.where(gal[index]["Type"] == 2)[0]]


sfr_type0 = {}
sfr_type1 = {}
sfr_type2 = {}
for i in range(len(model_names)):
    index = model_names[i]
    sfr_type0[index] = numpy.sum(gal_type0[index]["Sfr"],dtype = numpy.float64)
    sfr_type1[index] = numpy.sum(gal_type1[index]["Sfr"],dtype = numpy.float64)
    sfr_type2[index] = numpy.sum(gal_type2[index]["Sfr"],dtype = numpy.float64)

folder = "sfr/"
os.system("mkdir -p "+folder)
f = open(folder+"/"+z+".dat","w+")
for i in range(len(model_names)):
    index = model_names[i]
    print >> f,sfr_type0[index], sfr_type1[index],sfr_type2[index]
    f.close()
