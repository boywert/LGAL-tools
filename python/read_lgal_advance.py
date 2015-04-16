#import LGalaxyStruct
import numpy
import os
import sys

struct_lgalinput = numpy.dtype([
('Descendant',numpy.int32,1),
('FirstProgenitor',numpy.int32,1),
('NextProgenitor',numpy.int32,1),
('FirstHaloInFOFgroup',numpy.int32,1),
('NextHaloInFOFgroup',numpy.int32,1),
('Len',numpy.int32,1),
('M_Mean200',numpy.float32,1),
('M_Crit200',numpy.float32,1),
('M_TopHat',numpy.float32,1),
('Pos',numpy.float32,3),
('Vel',numpy.float32,3),
('VelDisp',numpy.float32,1),
('Vmax',numpy.float32,1),
('Spin',numpy.float32,3),
('MostBoundID',numpy.int64,1),
('SnapNum',numpy.int32,1),
('FileNr',numpy.int32,1),
('SubhaloIndex',numpy.int32,1),
('SubHalfMass',numpy.int32,1)
])

tree_properties_used = {}
for el in struct_lgalinput.names:
    tree_properties_used[el] = False

def read_lgal_input_tree(folder,file_prefix,firstfile,lastfile,filter_arr,verbose):
    dt = struct_lgalinput
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    filter_tuple = []
    for prop in dt.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,dt[prop]))
    filter_dtype = numpy.dtype(filter_tuple)
    output_Galaxy = numpy.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
        nTrees += this_nTrees
        this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
        nHalos += this_nHalos
        if(verbose):
            print "File ", ifile," nHalos = ",this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
        this_addedHalos = numpy.fromfile(f,dt,this_nHalos)
        addedHalos = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in dt.names:
            if(filter_arr[prop] is True):
                addedHalos[prop] = this_addedHalos[prop]
        output_Halos = numpy.append(output_Halos,addedHalos)
        f.close()
    return (nTrees,nHalos,nTreeHalos,output_Halos)

def read_lgaltree_advance(folder,file_prefix,firstfile,lastfile,filter_arr,dt,verbose):
    nHalos = 0
    filter_tuple = []
    for prop in dt.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,dt[prop]))
    filter_dtype = numpy.dtype(filter_tuple)
    output_Galaxy = numpy.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"galtree_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        one = dummy[0]
        dummy = numpy.fromfile(f,numpy.int32,1)
        structsize = dummy[0]
        if(structsize != dt.itemsize):
            print "size mismatch:",structsize,dt.itemsize
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        f.seek(structsize, os.SEEK_SET) 
        print "File ", ifile," nGals = ",this_nHalos
        this_addedGalaxy = numpy.fromfile(f,dt,this_nHalos)
        addedGalaxy = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in dt.names:
            if(filter_arr[prop] is True):
                addedGalaxy[prop] = this_addedGalaxy[prop]
        output_Galaxy = numpy.append(output_Galaxy,addedGalaxy)
       
      
        f.close()

    return (nHalos,output_Galaxy)


# This function return (nTrees,nHalos,nTreeHalos,Galaxy)
# The input are (folder,file_prefix,firstfile,lastfile [,filter_arr])
def readsnap_lgal_advance(folder,file_prefix,firstfile,lastfile,filter_arr,dt,verbose):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    filter_tuple = []
    for prop in dt.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,dt[prop]))
    filter_dtype = numpy.dtype(filter_tuple)
    output_Galaxy = numpy.array([],dtype=filter_dtype)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nTrees =  dummy[0]
        nTrees += this_nTrees
        dummy = numpy.fromfile(f,numpy.int32,1)
        this_nHalos = dummy[0]
        nHalos += this_nHalos
        if(verbose):
            print "File ", ifile," nGals = ",this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
        this_addedGalaxy = numpy.fromfile(f,dt,this_nHalos)
        addedGalaxy = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in dt.names:
            if(filter_arr[prop] is True):
                addedGalaxy[prop] = this_addedGalaxy[prop]
        output_Galaxy = numpy.append(output_Galaxy,addedGalaxy)
       
      
        f.close()
    return (nTrees,nHalos,nTreeHalos,output_Galaxy)

