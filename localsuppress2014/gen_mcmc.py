import numpy
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("../python/")
from  read_lgal_advance import *
import random
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

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
    tot_ntrees = 0
    tot_nhalos = 0
    tot_ntreehalos = numpy.array([],dtype=numpy.int32)
    tot_output_halos = numpy.array([],dtype=struct_lgalinput)
    tot_output_haloids_mcmc  = numpy.array([],dtype=struct_lgaldbidsinput)
    tot_nbins = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    tot_count = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    all_list = numpy.array(range(120,nFiles))
    filelist = numpy.array_split(all_list,size)[rank]
    for ifile in filelist:
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
                tot_output_haloids_mcmc = numpy.append(tot_output_haloids_mcmc,output_HaloIDs[rootindex[h]:rootindex[h]+nTreeHalos[h]])
                
    for i in range(lastsnap+1):
        for j in range(nbins):
            lbound = min_m+j*delta_logm
            rbound = lbound+delta_logm
            c_list = numpy.where((numpy.log10(tot_output_halos['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(tot_output_halos['M_Crit200']*gadget_m_conv/hubble_h) >=lbound) & (tot_output_halos['SnapNum'] == i) & (tot_output_haloids_mcmc["HaloID"] == tot_output_haloids_mcmc["FirstHaloInFOFgroup"]))[0]
            tot_count[j,i] += len(c_list)    

    

    comm.Barrier()
    comm.Reduce(None, [tot_ntrees], op=MPI.SUM, root=0)
    comm.Reduce(None, [tot_nhalos], op=MPI.SUM, root=0)
    tot_ntreehalos = comm.Gather(tot_ntreehalos, root=0)
    tot_output_halos = comm.Gather(tot_output_halos, root=0)
    tot_output_haloids_mcmc = comm.Gather(tot_output_haloids_mcmc, root=0)
    comm.Reduce(None, tot_count, op=MPI.SUM, root=0)
    comm.Reduce(None, tot_nbins, op=MPI.SUM, root=0)
    if rank==0:
        print tot_ntrees 
        print tot_nhalos 
        print tot_ntreehalos
        print tot_output_halos
        print tot_output_haloids_mcmc
        print tot_nbins
        print tot_count
    return 0

if __name__=="__main__":
    main()
    
