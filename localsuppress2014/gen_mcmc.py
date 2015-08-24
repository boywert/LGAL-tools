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
sample_bin = 1
nFiles = 128
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
    tot_weight = numpy.array([],dtype=numpy.float64)
    tot_nbins = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    tot_count = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)

    for ifile in range(100,nFiles):
        firstfile = ifile
        lastfile = ifile
        (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=True)
        rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
        print "Making total table ..."
        # compute weight table
        percent = 0
        for i in range(lastsnap+1):
            for j in range(nbins):
                lbound = min_m+j*delta_logm
                rbound = lbound+delta_logm
                t_list = numpy.where((numpy.log10(output_Halos['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(output_Halos['M_Crit200']*gadget_m_conv/hubble_h) >=lbound) & (output_Halos['SnapNum'] == i) & (output_HaloIDs["HaloID"] == output_HaloIDs["FirstHaloInFOFgroup"]))[0]
                tot_nbins[j,i] += len(t_list)
                percent += 100./(lastsnap+1)/nbins
                print "%d \% completed " % (percent)
        # sample data
        print "Sampling data ..."
        for j in range(nbins):
            lbound = min_m+j*delta_logm
            rbound = lbound+delta_logm
            print lbound,rbound
            r_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) >=lbound))[0]
            choose_list = random.sample(r_list,min(len(r_list),sample_bin))
            for h in choose_list:
                tot_ntrees += 1
                tot_nhalos += nTreeHalos[h]
                tot_ntreehalos = numpy.append(tot_ntreehalos,nTreeHalos[h])
                tot_output_halos = numpy.append(tot_output_halos,output_Halos[rootindex[h]:rootindex[h]+nTreeHalos[h]])
                tot_output_haloids_mcmc = numpy.append(tot_output_haloids_mcmc,output_HaloIDs[rootindex[h]:rootindex[h]+nTreeHalos[h]]["HaloID"])

    
    print tot_ntrees
    print tot_nhalos
    print tot_ntreehalos
    print tot_output_halos
    print tot_output_haloids_mcmc
    print tot_weight
    return 0

if __name__=="__main__":
    main()
    
