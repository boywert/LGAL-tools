import numpy
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("../python/")
from  read_lgal_advance import *
import random
hubble_h = 0.7
gadget_m_conv = 1.e10
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
lastsnap = 75
selected_file = 127
min_m = 8.
max_m = 10.
nbins = 50
delta_logm = int(round((max_m-min_m)/nbins))
def main():
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_%03d." % (lastsnap)
    flist = [selected_file]
    for i in flist:
        firstfile = i
        lastfile = i
        (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=False)
        rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
        for j in range(nbins):
            min = min_m+j*delta_logm
            max = min+delta_logm
            print numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv)
            #choose_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) >=min) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) <=max))
            #print choose_list

    return 0

if __name__=="__main__":
    main()
    
