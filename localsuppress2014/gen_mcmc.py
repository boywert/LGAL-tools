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
nbins = 30
delta_logm = (max_m-min_m)/nbins
sample_bin = 10
def main():
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_%03d." % (lastsnap)
    firstfile = 40
    lastfile = 40
    tot_ntrees = 0
    tot_nhalos = 0
    tot_ntreehalos = numpy.array([],dtype=numpy.int32)
    tot_output_halos = numpy.array([],dtype=struct_lgalinput)
    tot_output_haloids_mcmc  = numpy.array([],dtype=numpy.int64)
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=True)
    rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
    for j in range(nbins):
        lbound = min_m+j*delta_logm
        rbound = lbound+delta_logm
        print lbound,rbound
        choose_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) >=lbound))[0]
        choose_list = random.sample(choose_list,min(len(choose_list),sample_bin))
        for h in choose_list:
            tot_ntrees += 1
            tot_nhalos += nTreeHalos[h]
            tot_ntreehalos = numpy.append(tot_ntreehalos,nTreeHalos[h])
            tot_output_halos = numpy.append(tot_output_halos,output_Halos[rootindex[h]:rootindex[h+1]])
            tot_output_haloids_mcmc = numpy.append(tot_output_haloids_mcmc,output_HaloIDs[rootindex[h]]["HaloID"])

    print tot_ntrees
    print tot_nhalos
    print tot_ntreehalos
    return 0

if __name__=="__main__":
    main()
    
