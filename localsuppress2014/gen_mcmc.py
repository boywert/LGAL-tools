import matplotlib
matplotlib.use('Agg')
import numpy
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("../python/")
from  read_lgal_advance import *
import mass_fn

hubble_h = 0.7

zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
lastsnap = 75
selected_file = 127
def main():
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_%03d." % (lastsnap)
    flist = [selected_file]
    for i in flist:
        firstfile = i
        lastfile = i
        (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=False)
        rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
        #halos = output_Halos[index]
        #haloids = output_Halos[haloindex]
        a = numpy.histogram(numpy.log10(output_Halos["M_Crit200"]*1e10/hubble_h),bins=50,range=(8.,10.))
        print a 
     
    return 0

if __name__=="__main__":
    main()
    
