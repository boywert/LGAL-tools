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
#import test as mymodule
mymodule = CDLL('./test.so')
mymodule.make_sphere.argtypes = [c_int, c_float, POINTER(c_float),POINTER(c_float)]
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
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,0,511,filter[i],dt[i],1)

        R = numpy.empty(nGals[index]*8,dtype=c_float)
        pix = numpy.empty(nGals[index]*8,dtype=numpy.int64)
        pixmap = numpy.zeros(healpy.nside2npix(NSIDE),dtype=numpy.float64)
        pos =  gal[index]['Pos']
        print pos
        index_out = 0
        N = nGals[index]
        print N
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    index_in = 0
                    
                    while (index_in < N):
                        #print pos_tmp[index][0:3] , pos[index_in][0:3]
                        pos_tmp = pos[index_in,0:3]-500.*numpy.array([i,j,k])
                        R[index_out] = numpy.sqrt(pos_tmp[0]*pos_tmp[0]+pos_tmp[1]*pos_tmp[1]+pos_tmp[2]*pos_tmp[2])
                        pix[index_out] = healpy.pixelfunc.vec2pix(NSIDE,pos_tmp[0],pos_tmp[1],pos_tmp[2])
                        if ((R[index_out] >35) & (R[index_out] < 70)):
                            pixmap[pix[index_out]] += 1.0
                        index_in += 1
                        index_out += 1
        healpy.write_map("my_map_full.fits", pixmap/numpy.sum(pixmap))
        
        
        #mymodule.make_sphere(c_int(nGals[index]),c_float(500.0),pos.ctypes.data_as(POINTER(c_float)),pos_sphere.ctypes.data_as(POINTER(c_float)))
        # pos_sphere[:,0] = numpy.sqrt(pos_tmp[:,0]*pos_tmp[:,0]+pos_tmp[:,1]*pos_tmp[:,1]+pos_tmp[:,2]*pos_tmp[:,2]) 
        # pos_sphere[:,1] = numpy.arccos(pos_tmp[:,2]/pos_sphere[:,0])
        # pos_sphere[:,2] = numpy.arctan(pos_tmp[:,1]/pos_tmp[:,0])
        #print pos_sphere

        #healpix
        
        
        

def main():
    plot_coldgas(0.0)


if __name__=="__main__":
    main()
