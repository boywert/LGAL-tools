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
from ctypes import CDLL, POINTER, c_int, c_float, c_double
mymodule = CDLL('./test.so')
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



def plot_coldgas(z):
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
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,0,5,filter[i],dt[i],1)

    for i in range(len(model_names)):
        index = model_names[i]
        pos =  gal[index]['Pos'].astype(numpy.float32)
        c = numpy.empty((nGals[index],3),dtype=numpy.float32)
        start = timer()
        i = 0
        while i<1000:
            c[:,0] = numpy.sqrt(pos[:,0]*pos[:,0]+pos[:,1]*pos[:,1]+pos[:,2]*pos[:,2])
            c[:,1] = numpy.acos(pos[:,2]/c[:,0])
            c[:,2] = numpy.atan(pos[:,1]/pos[:,0])
            i += 1
        end = timer()
        print "Python:", (end - start)/1000.0
        start = timer()
        mymodule.cart2sphere1(c_int(nGals[index]),pos.ctypes.data_as(POINTER(c_float)),c.ctypes.data_as(POINTER(c_float)))
        end = timer()
        print "fortran:", (end - start)/1000.0
        start = timer()
        mymodule.cart2sphere2(c_int(nGals[index]),pos.ctypes.data_as(POINTER(c_float)),c.ctypes.data_as(POINTER(c_float)))
        end = timer()
        print "mkl vector fortran:", (end - start)/1000.0
        start = timer()
        mymodule.cart2sphere3(c_int(nGals[index]),pos.ctypes.data_as(POINTER(c_float)),c.ctypes.data_as(POINTER(c_float)))
        end = timer()
        print "mkl vector fortran (tweaked):", (end - start)/1000.0
        
def main():
    plot_coldgas(0.0)


if __name__=="__main__":
    main()
