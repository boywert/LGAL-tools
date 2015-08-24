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
selected_file = 100
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
        print firsthalo
    return 0

if __name__=="__main__":
    main()
    
