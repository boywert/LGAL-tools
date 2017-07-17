from astropy.cosmology import FlatLambdaCDM
import sqlite3
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
mymodule.make_sphere.argtypes = [c_int, c_float, _twodimp, _twodimp, _twodimp, _twodimp]
import healpy
import sqlite3
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


def main():
    
    first_z = 0.0
    last_z = z_from_nu(1220.0)
    f_step = 0.5 #MHz
    print "a", a_from_z(first_z), a_from_z(last_z)
    print "f", nu_from_z(first_z), nu_from_z(last_z)
    print "t", t_from_z(first_z), t_from_z(last_z)
    print "d", cosmo.comoving_distance(first_z)*.73,cosmo.comoving_distance(last_z)*.73
    #construct table for lookup f-d
    f_array = numpy.arange(nu_from_z(first_z),nu_from_z(last_z)-f_step,-0.1)
    d_array = numpy.empty(len(f_array),dtype=numpy.float32)
    d_array[:] = cosmo.comoving_distance(z_from_nu(f_array[:])).value*0.73
    print f_array
    print d_array

    
    alist = numpy.loadtxt(alist_file)
    alist = alist[(alist >= a_from_z(last_z)) & (alist <= a_from_z(first_z))]
    alist.sort()
    alist = alist[::-1]
    print alist,a_from_z(last_z),a_from_z(first_z)

    fc_list = numpy.arange(nu_from_z(first_z),nu_from_z(last_z)-f_step,-1*f_step)
    fb_list = numpy.empty(len(fc_list)-1,dtype = numpy.float32)
    for i in range(len(fc_list)-1):
        fb_list[i] = 0.5*(fc_list[i]+fc_list[i+1])
    ngals = []
    gals = []
    start_r = 0.0
    totalNgals = []
    for i in range(len(alist)):
        z = "%10.3f" % (z_from_a(alist[i]))
        if i < len(alist)-1:
            alist_distance = cosmo.comoving_distance(z_from_a(alist[i+1])).value*0.73
        else:
            alist_distance = cosmo.comoving_distance(last_z).value*0.73
        ngal_i,gal_i,pos_i,vR_i = readgal(float(z),0,5)
        #ngals.append(ngal_i)
        #pos.append(pos_i)
        #vR.append(vR_i)
        fullgal = numpy.empty(ngal_i*8,dtype = gal_i.dtype)
        if "FileUniqueGalID" in gal_i.dtype.names:
            for j in range(8):
                fullgal[ngal_i*j:ngal_i*(j+1)] = gal_i
                fullgal[ngal_i*j:ngal_i*(j+1)]['FileUniqueGalID'] += ngal_i*j

        else:
            for j in range(8):
                fullgal[ngal_i*j:ngal_i*(j+1)] = gal_i      

        del(gal_i)
        gallist = numpy.where((pos_i[:,0] >= start_r) & (pos_i[:,0] <= alist_distance))[0]
        print "z = ",z,"a=",alist[i],"r = ",start_r,"-",alist_distance
        #store data
        ogal = numpy.empty(len(gallist),dtype=db_struct)
        ogal['PosX'] = fullgal['Pos'][gallist,0]
        ogal['PosY'] = fullgal['Pos'][gallist,1]
        ogal['PosZ'] = fullgal['Pos'][gallist,2]
        ogal['VelX'] = fullgal['Vel'][gallist,0]
        ogal['VelY'] = fullgal['Vel'][gallist,1]
        ogal['VelZ'] = fullgal['Vel'][gallist,2]
        ogal['StellarMass'] = fullgal['StellarMass'][gallist]
        ogal['ColdGas'] = fullgal['ColdGas'][gallist]
        coldtostellar =  ogal['ColdGas']/ogal['StellarMass']
        ogal['PosR'] = pos_i[gallist,0]
        ogal['PosTheta'] = pos_i[gallist,1]
        ogal['PosPhi'] = pos_i[gallist,2]
        ogal['VelR'] = vR_i[gallist,0]
        ogal['VelTheta'] = vR_i[gallist,1]
        ogal['VelPhi'] = vR_i[gallist,2]
        ogal['Frequency'] = numpy.interp(ogal['PosR'],d_array,f_array)
        ogal['LuminosityDistance'] = ogal['PosR']*(z_from_nu(ogal['Frequency'][:])+1)
        ogal['NeutralH'] = ogal['ColdGas']*0.41/(numpy.power(coldtostellar,-0.52)+numpy.power(coldtostellar,0.56))/0.73
        ogal['Intensity'] = ogal['NeutralH']/49.8*numpy.power(ogal['LuminosityDistance'],-2)
        gals.append(ogal)
        start_r = alist_distance
        totalNgals.append(len(gallist))
    totalNgals = numpy.array(totalNgals)
    db_gal = numpy.empty(numpy.sum(totalNgals),dtype=db_struct)
    first_gal = 0
    for i in range(len(gal)):
        db_gal[first_gal:first_gal+totalNgals[i]] = gals[i]
        first_gal += totalNgals[i]
    print db_gal
    return


    # i = len(alist)-1
    # z = "%10.3f" % (z_from_a(alist[i]))
    
    # ngal_i,gal_i,pos_i,vR_i = readgal(float(z))
    # gallist = numpy.where((pos_i[:,0] >= start_r) & (pos_i[:,0] <= alist_distance))[0]
    # print "z = ",z,"a=",alist[i],"r = ",start_r,"-",alist_distance
    # start_r = alist_distance
    
                        
    # return


    # Rb_list = numpy.empty(len(fb_list),dtype = numpy.float32)
    # Rb_list[:] = cosmo.comoving_distance(z_from_nu(fb_list[:])).value*0.73

    # for i in range(len(Rb_list)-1):
    #     r_check = len(alist_distance)-1
    #     toggle = 0
    #     while (toggle==0) & (r_check >= 0):
    #         if alist_distance[r_check] > Rb_list[i]:
    #             toggle = 1
    #             break
    #         r_check -= 1
    #     gallist = numpy.where((pos[r_check][:,0] >= Rb_list[i]) & (pos[r_check][:,0] <= Rb_list[i+1]))[0]
    #     print len(gallist)

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
