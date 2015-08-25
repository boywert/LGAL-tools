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
lastsnap = 75
min_m = 8.
max_m = 10.
nbins = 20
delta_logm = (max_m-min_m)/nbins
sample_bin = 1
nFiles = 128
def main(argv):
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    #file_prefix = "trees_%03d." % (lastsnap)
    tot_ntrees = numpy.zeros(1)
    tot_nhalos = numpy.zeros(1)
    tot_ntreehalos = numpy.array([],dtype=numpy.int32)
    tot_output_halos = numpy.array([],dtype=struct_lgalinput)
    tot_output_haloids_mcmc  = numpy.array([],dtype=struct_lgaldbidsinput)
    tot_nbins = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    tot_count = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    f_tot_nbins = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    f_tot_count = numpy.zeros((nbins,lastsnap+1),dtype=numpy.int64)
    all_list = range(124,nFiles)
    random.shuffle(all_list)
    all_list = numpy.array(all_list)
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
            r_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) <=rbound) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv/hubble_h) >=lbound))[0]
            choose_list = random.sample(r_list,min(len(r_list),sample_bin))
            for h in choose_list:
                tot_ntrees[0] += 1
                tot_nhalos[0] += nTreeHalos[h]
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
    a = numpy.zeros(1,dtype=numpy.int32)
    b = numpy.zeros(1,dtype=numpy.int32)
    comm.Reduce(tot_ntrees, a, op=MPI.SUM, root=0)
    comm.Reduce(tot_nhalos, b, op=MPI.SUM, root=0)
   
    tot_ntreehalos = comm.gather(tot_ntreehalos, root=0)
    tot_output_halos = comm.gather(tot_output_halos, root=0)
    tot_output_haloids_mcmc = comm.gather(tot_output_haloids_mcmc, root=0)
    comm.Reduce(tot_count, f_tot_count, op=MPI.SUM, root=0)
    comm.Reduce(tot_nbins, f_tot_nbins, op=MPI.SUM, root=0)
    if rank==0:
        print a 
        print b
        f_tot_ntreehalos = numpy.concatenate(tot_ntreehalos)
        f_tot_output_halos = numpy.concatenate(tot_output_halos)
        f_tot_output_haloids_mcmc = numpy.concatenate(tot_output_haloids_mcmc)
        print f_tot_ntreehalos
        print f_tot_output_halos
        print f_tot_output_haloids_mcmc
        print f_tot_nbins
        print f_tot_count
        fp_halo = open("trees","wb")
        fp_haloids = open("tree_dbids","wb")
        a.tofile(fp_halo)
        b.tofile(fp_halo)
        f_tot_ntreehalos.tofile(fp_halo)
        f_tot_output_halos.tofile(fp_halo)
        fp_halo.close()
        f_tot_output_haloids_mcmc.tofile(fp_haloids)
        fp_haloids.close()
    return 0

if __name__=="__main__":
    main(sys.argv)
    
