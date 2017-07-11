from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.25, Tcmb0=2.725)
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
mymodule.make_sphere.argtypes = [c_int, c_float, _twodimp, _twodimp, _twodimp, _onedimp]
import healpy
from timeit import default_timer as timer
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
    filter['Pos'] = True
    filter['Vel'] = True
    filter['FileUniqueGalID'] = True
    filter['FileUniqueGalCentralID'] = True
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
def readgal(z):
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    r = {}
    theta = {}
    phi = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if index[:4] == "lgal":
            zz = "%10.2f"%(z)
        elif index[:4] == "sage":
            zz = "%10.3f"%(z)
        file_prefix = "model_z"+zz.strip()
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,5,5,filter[i],dt[i],1)
            pos = numpy.ascontiguousarray(gal[index]['Pos'])
            vel = numpy.ascontiguousarray(gal[index]['Vel'])
            pos_sphere = numpy.empty((nGals[index]*8,3),dtype=numpy.float32)
            vel_R = numpy.empty(nGals[index]*8,dtype=numpy.float32)
            mymodule.make_sphere(c_int(nGals[index]),c_float(500.0),pos,vel,pos_sphere,vel_R)
        return gal[index]
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

def main():
    first_z = 0.0
    last_z = z_from_nu(1230.0)
    print "a", a_from_z(first_z), a_from_z(last_z)
    print "f", nu_from_z(first_z), nu_from_z(last_z)
    print "t", t_from_z(first_z), t_from_z(last_z)
    alist = numpy.loadtxt(alist_file)
    alist = alist[(alist >= a_from_z(last_z)) & (alist <= a_from_z(first_z))]
    gal = []
    for a in alist:
        z = "%10.3f" % (z_from_a(a))
        gal.append(readgal(float(z)))

    f_step = 0.5 #MHz
    fc_list = numpy.arange(nu_from_z(first_z),nu_from_z(last_z)-f_step,-1*f_step)
    Rc_list = numpy.empty(len(fc_list),dtype = numpy.float32)
    for i in range(len(Rc_list)):
        Rc_list[i] = cosmo.comoving_distance(z_from_nu(fc_list[i])).value*0.73
    print Rc_list

    # #track gals backward
    # for igal in gal[len(gal)-1]:
    #     id = igal['FileUniqueGalID']
    #     isnap = len(gal)-2
    #     while (id > -1) & (isnap > -1):
    #         listgal = gal[isnap][gal[isnap]['FileUniqueGalID'] == id]
    #         if len(listgal) == 0:
    #             id = -1
    #         isnap -= 1
if __name__ == "__main__":
    main()
