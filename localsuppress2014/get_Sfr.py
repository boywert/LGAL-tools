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
sys.path.append("../python/")
import read_lgal_advance as read_lgal

rank = sys.argv[1]
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
z = zlist[int(sys.argv[1])].strip()

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

file_prefix = "SA_z"+z

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


gal_lomass ={}
gal_himass = {}
gal_type0 = {}
gal_type1 = {}
gal_type2 = {}

sfr_hist = {}
fig = pylab.figure()
fig2 = pylab.figure()
ax = fig.add_subplot(111)
ax2 = fig2.add_subplot(111)
for i in range(len(model_names)):
    index = model_names[i]
    gal_type0[index]  = gal[index][numpy.where(gal[index]["Type"] == 0)[0]]
    gal_type1[index]  = gal[index][numpy.where(gal[index]["Type"] == 1)[0]]
    gal_type2[index]  = gal[index][numpy.where(gal[index]["Type"] == 2)[0]]
    gal_lomass[index] = gal[index][numpy.where(gal[index]["HaloM_Crit200"] < 0.1/h0)[0]]
    gal_himass[index] = gal[index][numpy.where(gal[index]["HaloM_Crit200"] > 0.1/h0)[0]]
    (sfr_bin_x,sfr_bin_y) = sfr_massbin_fn(gal[index],mass_min=1e8,mass_max=1.e12,nbins=20)
    (massfn_x,massfn_y) = M200c_mass_fn_gal(gal[index],mass_min=1e8,mass_max=1.e12,nbins=20)
    ax.plot(sfr_bin_x,numpy.cumsum(sfr_bin_y),label=model_labels[i])
    ax2.plot(massfn_x,massfn_y,label=model_labels[i])

ax.set_xlabel(r"$\mathrm{\log_{10}(M_{200c}/M_\odot)}$")
ax2.set_xlabel(r"$\mathrm{\log_{10}(M_{200c}/M_\odot)}$")
ax.set_ylabel(r"$\mathrm{\phi(Mpc^{-3} dex^{-1})}$")
ax2.set_ylabel(r"$\mathrm{\phi(Mpc^{-3} dex^{-1})}$")
ax2.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
leg2 = ax2.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg2.get_frame().set_linewidth(0)
fig.savefig("sum_sfr_bin_"+str(rank)+".pdf")
fig2.savefig("massfn_"+str(rank)+".pdf")

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


