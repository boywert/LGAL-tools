#import LGalaxyStruct
import numpy
import os
import sys
import time
import cPickle as pickle
import os.path
import hashlib

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

struct_lgaldbidsinput = numpy.dtype([
    ('HaloID',numpy.int64,1),
    ('FileTreeNr',numpy.int64,1),
    ('FirstProgenitor',numpy.int64,1),
    ('LastProgenitor',numpy.int64,1),
    ('NextProgenitor',numpy.int64,1),
    ('Descendant',numpy.int64,1),
    ('FirstHaloInFOFgroup',numpy.int64,1),
    ('NextHaloInFOFgroup',numpy.int64,1),
    ('Redshift',numpy.float64,1),
    ('PeanoKey',numpy.int32,1),
    ('dummy',numpy.int32,1)
])

struct_lgaldbidsinput_MRII = numpy.dtype([
    ('HaloID',numpy.int64,1),
    ('FileTreeNr',numpy.int64,1),
    ('FirstProgenitor',numpy.int64,1),
    ('LastProgenitor',numpy.int64,1),
    ('NextProgenitor',numpy.int64,1),
    ('Descendant',numpy.int64,1),
    ('FirstHaloInFOFgroup',numpy.int64,1),
    ('NextHaloInFOFgroup',numpy.int64,1),
    ('MainLeafID',numpy.int64,1),
    ('Redshift',numpy.float64,1),
    ('PeanoKey',numpy.int32,1),
    ('dummy',numpy.int32,1)
])

tree_properties_used = {}
for el in struct_lgalinput.names:
    tree_properties_used[el] = False

tree_dbid_properties_used = {}
for el in struct_lgaldbidsinput.names:
    tree_dbid_properties_used[el] = False

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
    output_Halos = numpy.array([],dtype=filter_dtype)
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

def read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose):
    nTrees = 0
    nHalos = 0
    tree_index = numpy.zeros(lastfile-firstfile+1,dtype=numpy.int32)
    halo_index = numpy.zeros(lastfile-firstfile+1,dtype=numpy.int32)
    i = 0
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/trees_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
        f = open(filename,"rb")
        tree_index[i] = numpy.fromfile(f,numpy.int32,1)[0]
        halo_index[i] = numpy.fromfile(f,numpy.int32,1)[0]
        f.close()
        i+=1
    tree_findex = numpy.cumsum(tree_index)-tree_index
    halo_findex = numpy.cumsum(halo_index)-halo_index
    
    nTreeHalos = numpy.zeros(numpy.sum(tree_index),dtype=numpy.int32)
    output_Halos = numpy.zeros(numpy.sum(halo_index),dtype=struct_lgalinput)
    output_HaloIDs = numpy.zeros(numpy.sum(halo_index),dtype=struct_lgaldbidsinput)
    i=0
    for ifile in range(firstfile,lastfile+1):
        start = time.time()
        filename = folder+'/trees_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
        nTrees += this_nTrees
        this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
        nHalos += this_nHalos
        if(verbose):
            print "File ", ifile,filename,"nTrees = ",this_nTrees," nHalos = ",this_nHalos
        #addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)

        nTreeHalos[tree_findex[i]:tree_findex[i]+tree_index[i]] = numpy.fromfile(f,numpy.int32,this_nTrees)
        #this_addedHalos = numpy.fromfile(f,struct_lgalinput,this_nHalos)
        output_Halos[halo_findex[i]:halo_findex[i]+halo_index[i]] = numpy.fromfile(f,struct_lgalinput,this_nHalos)
        f.close()
        filename = folder+'/tree_dbids_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
        f = open(filename,"rb")
        #this_addedHalos = numpy.fromfile(f,struct_lgaldbidsinput,this_nHalos)
        output_HaloIDs[halo_findex[i]:halo_findex[i]+halo_index[i]] = numpy.fromfile(f,struct_lgaldbidsinput,this_nHalos)
        f.close()
        end = time.time()
        if(verbose == 2):
            print end-start,"s"
        i += 1
    return (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs)

def read_lgal_input_fulltrees_withids_advance(folder,lastsnap,firstfile,lastfile,trees_filter,tree_dbids_filter,verbose):
    nTrees = 0
    nHalos = 0
    nTreeHalos = numpy.array([],dtype=numpy.int32)
    output_Halos = numpy.array([],dtype=struct_lgalinput)
    output_HaloIDs = numpy.array([],dtype=struct_lgaldbidsinput)
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/trees_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
        nTrees += this_nTrees
        this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
        nHalos += this_nHalos
        if(verbose):
            print "File ", ifile," nHalos = ",this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        
        nTreeHalos = numpy.append(nTreeHalos,addednTreeHalos)
        this_addedHalos = numpy.fromfile(f,struct_lgalinput,this_nHalos)
        output_Halos = numpy.append(output_Halos,this_addedHalos)
        f.close()
        filename = folder+'/tree_dbids_'+"%03d"%(lastsnap)+'.'+"%d"%(ifile)
        f = open(filename,"rb")
        this_addedHalos = numpy.fromfile(f,struct_lgaldbidsinput,this_nHalos)
        output_HaloIDs = numpy.append(output_HaloIDs,this_addedHalos)
        f.close()
    return (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs)

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



def get_filter_array_to_string(filter_arr,dt):
    string = ""
    for prop in dt.names:
        if(filter_arr[prop] is True):
            string = string+"1"
        else:
            string = string+"0"
    return string

# This function return (nTrees,nHalos,nTreeHalos,Galaxy)
def readsnap_lgal_advance(folder,file_prefix,firstfile,lastfile,filter_arr,dt,verbose=1,cache_on=False):
    startx = time.time()
    cache_filename = hashlib.sha1(folder+"_"+file_prefix).hexdigest()+"_"+str(firstfile)+"_"+str(lastfile)+"_"+ str(int(get_filter_array_to_string(filter_arr,dt),2)) +".pickle"
    if cache_on is True:
        if os.path.isfile(cache_filename):
            print "Read data from Pickled file: "+cache_filename
            (nTrees,nHalos,nTreeHalos,output_Galaxy) = pickle.load( open( cache_filename, "rb" ) )
            endx = time.time()
            if(verbose > 0):
                print "Read ",folder,"file",firstfile,"-",lastfile,":",endx-startx,"s"
                return (nTrees,nHalos,nTreeHalos,output_Galaxy)
    nTrees = 0
    nHalos = 0
    filter_tuple = []
    for prop in dt.names:
        if(filter_arr[prop] is True):
            filter_tuple.append((prop,dt[prop]))
    filter_dtype = numpy.dtype(filter_tuple)

    tree_index = numpy.zeros(lastfile-firstfile+1,dtype=numpy.int32)
    halo_index = numpy.zeros(lastfile-firstfile+1,dtype=numpy.int32)
    i = 0
    for ifile in range(firstfile,lastfile+1):
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        tree_index[i] = numpy.fromfile(f,numpy.int32,1)[0]
        halo_index[i] = numpy.fromfile(f,numpy.int32,1)[0]
        f.close()
        i+=1
    tree_findex = numpy.cumsum(tree_index)-tree_index
    halo_findex = numpy.cumsum(halo_index)-halo_index
    
    nTreeHalos = numpy.zeros(numpy.sum(tree_index),dtype=numpy.int32)
    output_Galaxy = numpy.zeros(numpy.sum(halo_index),dtype=filter_dtype)
    i = 0
    for ifile in range(firstfile,lastfile+1):
        start = time.time()
        filename = folder+'/'+file_prefix+"_"+"%d"%(ifile)
        f = open(filename,"rb")
        this_nTrees = numpy.fromfile(f,numpy.int32,1)[0]
        nTrees += this_nTrees
        this_nHalos = numpy.fromfile(f,numpy.int32,1)[0]
        nHalos += this_nHalos
        if(verbose==2):
            print "File ", ifile," nGals = ",this_nHalos
        addednTreeHalos = numpy.fromfile(f,numpy.int32,this_nTrees)
        nTreeHalos[tree_findex[i]:tree_findex[i]+tree_index[i]] = addednTreeHalos
        this_addedGalaxy = numpy.fromfile(f,dt,this_nHalos)
        addedGalaxy = numpy.zeros(this_nHalos,dtype=filter_dtype)
        for prop in dt.names:
            if(filter_arr[prop] is True):
                addedGalaxy[prop] = this_addedGalaxy[prop]
        output_Galaxy[halo_findex[i]:halo_findex[i]+halo_index[i]] = addedGalaxy
        f.close()
        i += 1
        end = time.time()
        if(verbose == 2):
            print end-start,"s"
    endx = time.time()
    pickle.dump( (nTrees,nHalos,nTreeHalos,output_Galaxy), open( cache_filename, "wb" ))
    if(verbose > 0):
        print "Read ",folder,"file",firstfile,"-",lastfile,":",endx-startx,"s"
    return (nTrees,nHalos,nTreeHalos,output_Galaxy)



