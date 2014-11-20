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
    for index in filter.keys():
        filter[index] = True

    filter['Type'] = True
    filter['HaloIndex'] = True
    filter['Sfr'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    filter['ICM'] = True
    filter['HotGas'] = True
    filter['ColdGas'] = True
    filter['sfh_numbins'] = True
    filter['sfh_BulgeMass'] = True
    filter['sfh_DiskMass'] = True
    filter['GasDiskRadius'] = True
    filter['CoolingRate'] = True
    filter['EjectedMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

firstfile = 1
lastfile = 1
# read galtree
filelist = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zf = open(filelist,"r")
zlist = zf.readlines()
(f1,t1) = loadfilter("/mnt/lustre/scratch/cs390/codes/47Mpc/L-Galaxies_development/galtree/LGalaxyStruct.py")
(nGals_gt,galtree) = read_lgal.read_lgaltree_advance("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","SA_",firstfile,lastfile,f1,t1,0)
i = 0
for z in zlist:
    z = z.strip()
    print z
    gal1 = galtree[numpy.where(galtree["SnapNum"] == i)[0]]
    (f2,t2) = loadfilter("/mnt/lustre/scratch/cs390/codes/47Mpc/L-Galaxies_development/gal/LGalaxyStruct.py")
    (nTrees_g,nGals_g,nTreeGals_g,gal2) = read_lgal.readsnap_lgal_advance("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","SA_z"+z,firstfile,lastfile,f2,t2,0)
    print gal2["Sfr"][numpy.where(gal2["Sfr"] > 100.)[0]]
    # (f3,t3) = loadfilter("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py")
    # (nTrees_g,nGals_g,nTreeGals_g,gal3) = read_lgal.readsnap_lgal_advance("/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","SA_z"+z,firstfile,lastfile,f3,t3,0)
    
    # (f4,t4) = loadfilter("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py")
    # (nTrees_g,nGals_g,nTreeGals_g,gal4) = read_lgal.readsnap_lgal_advance("/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","SA_z"+z,firstfile,lastfile,f4,t4,0)
    i += 1
    print numpy.sum(gal1["Sfr"],dtype=numpy.float64),numpy.sum(gal2["Sfr"],dtype=numpy.float64)#,numpy.sum(gal3["Sfr"],dtype=numpy.float64),numpy.sum(gal4["Sfr"],dtype=numpy.float64)
    print len(gal1),len(gal2)#,len(gal3),len(gal4)
    
