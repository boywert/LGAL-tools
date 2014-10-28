#import LGalaxyStruct
import numpy
import os
import sys
# This function return (nTrees,nHalos,nTreeHalos,Galaxy)
# The input are (folder,file_prefix,firstfile,lastfile [,filter_arr])
def readsnap_lgal_advance(folder,file_prefix,firstfile,lastfile,filter_arr,dt):
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

