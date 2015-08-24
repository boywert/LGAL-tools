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
        firsthalo = numpy.cumsum(nTreeHalos)-nTreeHalos
        haloindex = numpy.where(output_Halos['SnapNum'] == lastsnap)
        halos = output_Halos[haloindex]
        haloids = output_Halos[haloindex]
        firsthalo2 = numpy.where((output_HaloIDs["FirstHaloInFOFgroup"] == output_HaloIDs["HaloID"]) & (output_Halos['SnapNum'] == lastsnap))[0]
        print len(firsthalo),len(firsthalo2)
        j = 0
        while firsthalo[j] == firsthalo2[j]:
            print j, firsthalo[j],firsthalo2[j]
            j += 1
    return 0

if __name__=="__main__":
    main()
    
