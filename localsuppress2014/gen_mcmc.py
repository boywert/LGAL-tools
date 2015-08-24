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
min_m = 8.
max_m = 10.
nbins = 20
delta_logm = (max_m-min_m)/nbins
def main():
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_%03d." % (lastsnap)

    firstfile = 127
    lastfile = 127
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=False)
    rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
    for j in range(nbins):
        lbound = min_m+j*delta_logm
        rbound = lbound+delta_logm
        print lbound,rbound
        choose_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) >=lbound))[0]
        choose_list = random.sample(choose_list,min(len(choose_list),20))
        print rootindex[choose_list]

    return 0

if __name__=="__main__":
    main()
    
