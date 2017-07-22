from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.25, Tcmb0=2.725)
import tables
from mass_fn import *
from globalconf import *
from math import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct
import add_observations
sys.path.append("../python/")
import read_lgal_advance as read_lgal
import timeit
from numpy.ctypeslib import ndpointer
from ctypes import CDLL, POINTER, c_int, c_float, c_double
#import test as mymodule
mymodule = CDLL('./test.so')
_twodimp = ndpointer(dtype=c_float,ndim=2)
_onedimp = ndpointer(dtype=c_float,ndim=1)
arg2 = ndpointer(ndim=2)
arg3 = ndpointer(shape=(10,10))
mymodule.make_sphere.argtypes = [c_int, c_float, _twodimp, _twodimp, _twodimp, _twodimp]
import healpy
from timeit import default_timer as timer
rank = "0"
os.system("mkdir -p ../tmp/"+rank)
db_struct = numpy.dtype([
('PosX'                      , numpy.float32),
('PosY'                      , numpy.float32),
('PosZ'                      , numpy.float32),
('PosR'                      , numpy.float32),
('PosTheta'                  , numpy.float32),
('PosPhi'                    , numpy.float32),
('VelX'                      , numpy.float32),
('VelY'                      , numpy.float32),
('VelZ'                      , numpy.float32),
('VelR'                      , numpy.float32),
('VelTheta'                  , numpy.float32),
('VelPhi'                    , numpy.float32),
('StellarMass'               , numpy.float32),
('ColdGas'                   , numpy.float32), 
('Healpix'                   , numpy.int32),
('Frequency'                 , numpy.float32),
('LuminosityDistance'        , numpy.float32),
('NeutralH'                  , numpy.float32),
('Intensity'                 , numpy.float32)])

def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False    
    filter['Pos'] = True
    filter['Vel'] = True
    filter['StellarMass'] = True
    filter['ColdGas'] = True
    filter['Mvir'] = True
    filter['FileUniqueGalID'] = True
    #filter['FileUniqueGalCentralID'] = True
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

NSIDE = 2048
f21cm  = 1420.4057517667 #MHz
def readgal(z,i_model,i_file):
        i = i_model
        index = model_names[i]
        if index[:4] == "lgal":
            zz = "%10.2f"%(z)
        elif index[:4] == "sage":
            zz = "%10.3f"%(z)
        file_prefix = "model_z"+zz.strip()
        (nTrees,nGals,nTreeGals,gal) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,i_file,i_file,filter[i],dt[i],1)
        pos = numpy.ascontiguousarray(gal['Pos'])
        vel = numpy.ascontiguousarray(gal['Vel'])
        pos_sphere = numpy.empty((nGals*8,3),dtype=numpy.float32)
        vel_R = numpy.empty((nGals*8,3),dtype=numpy.float32)
        mymodule.make_sphere(c_int(nGals),c_float(500.0),pos,vel,pos_sphere,vel_R)
        return nGals,gal,pos_sphere,vel_R
def nu_from_a(a): #MHz
    return a*f21cm
def a_from_nu(f):
    return f/f21cm
def nu_from_z(z):
    return f21cm/(1.+z)
def z_from_nu(f):
    return f21cm/f - 1.0
def t_from_a(a):
    return cosmo.age(1./a - 1.0)
def t_from_z(z):
    return cosmo.age(z)
def a_from_z(z):
    return 1./(z+1.)
def z_from_a(a):
    return 1./a - 1.0

alist_file =  "/lustre/HI_FAST/SAM_code/LGAL/input/zlists/zlist_MR.txt"


def read_lightcone(dataset,dataname,file):
    gal = numpy.load('model_%s_%d.npy'%(dataname,file))
    return gal
def dtype_struct_to_descr(dtype):
    desc = {}
    for element in dtype.descr:
        if element[1] == "<f4":
            desc[element[0]] = tables.Float32Col()
        if element[1] == "<i4":
            desc[element[0]] = tables.Int32Col()
    return desc
def main():
    # import sqlite3
    # print "Reading SQLite3 table"
    # start = timer()
    # conn = sqlite3.connect('/share/data2/VIMALA/Lightcone/example.db')
    # c = conn.cursor()
    # c.execute("SELECT * FROM lightcone WHERE PosPhi > %f" %(numpy.pi/2))

    # a = c.fetchone()
    # while (a != None):
    #     a = c.fetchone()
        
    # conn.close()
    # end = timer()
    # print "sqlite uses ",end-start

    db_desc = dtype_struct_to_descr(db_struct)
    print "PyTables HDF5 file"
    start = timer()
    h5f = tables.open_file('/share/data2/VIMALA/Lightcone/example.hdf5', 'r')
    tbl = h5f.root.table_name
    for row in tbl.where('(PosPhi > %f)' % (numpy.pi/2)):
        print row
    h5f.close()
    end = timer()
    print "/share/data2 uses ",(end-start)
    
    print "PyTables HDF5 file"
    start = timer()
    h5f = tables.open_file('/lustre/HI_FAST/VIMALA/Lightcone/example.hdf5', 'r')

    tbl = h5f.root.table_name
    for row in tbl.where('(PosPhi > %f)' % (numpy.pi/2)):
        print row[:]
    h5f.close()
    end = timer()
    print "/lustre uses ",(end-start)
if __name__ == "__main__":
    main()
