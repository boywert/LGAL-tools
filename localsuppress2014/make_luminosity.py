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


pylab.rc('text', usetex=True)


exit()

gal_lomass ={}
gal_himass = {}
gal_type0 = {}
gal_type1 = {}
gal_type2 = {}


for i in range(len(model_names)):
    index = model_names[i]
    gal_type0[index]  = gal[index][numpy.where(gal[index]["Type"] == 0)[0]]
    gal_type1[index]  = gal[index][numpy.where(gal[index]["Type"] == 1)[0]]
    gal_type2[index]  = gal[index][numpy.where(gal[index]["Type"] == 2)[0]]
    gal_lomass[index] = gal[index][numpy.where(gal[index]["HaloM_Crit200"] < 0.1*h0)[0]]
    gal_himass[index] = gal[index][numpy.where(gal[index]["HaloM_Crit200"] > 0.1*h0)[0]]
sfr_type0 = {}
sfr_type1 = {}
sfr_type2 = {}
sfr_himass = {}
sfr_lomass = {}
for i in range(len(model_names)):
    index = model_names[i]
    sfr_type0[index] = numpy.sum(gal_type0[index]["Sfr"],dtype = numpy.float64)
    sfr_type1[index] = numpy.sum(gal_type1[index]["Sfr"],dtype = numpy.float64)
    sfr_type2[index] = numpy.sum(gal_type2[index]["Sfr"],dtype = numpy.float64)
    sfr_himass[index] = numpy.sum(gal_himass[index]["Sfr"],dtype = numpy.float64)
    sfr_lomass[index] = numpy.sum(gal_lomass[index]["Sfr"],dtype = numpy.float64)
    
folder = "sfr/"
os.system("mkdir -p "+folder)
f = open(folder+"/"+z+".dat","w+")
print >> f,"#type0 type1 type2 Lo-Mass Hi-Mass"
for i in range(len(model_names)):
    index = model_names[i]
    print >> f, sfr_type0[index], sfr_type1[index],sfr_type2[index],sfr_lomass[index],sfr_himass[index]
f.close()


