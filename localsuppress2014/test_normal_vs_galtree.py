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

firstfile = 10
lastfile = 10
# read galtree

(f1,t1) = loadfilter("/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/LGalaxyStruct.py")
(nGals_gt,galtree) = read_lgal.read_lgaltree_advance("/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/","SA_",firstfile,lastfile,f1,t1,0)

gal1 = galtree[numpy.where(galtree["SnapNum"] == 75)[0]]
(f2,t2) = loadfilter("/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/sams/43000.00/LGalaxyStruct.py")
(nTrees_g,nGals_g,nTreeGals_g,gal2) = read_lgal.readsnap_lgal_advance("/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/sams/43000.00/","SA_z6.00",firstfile,lastfile,f2,t2,0)




