import matplotlib 
import sys
from numpy import *
matplotlib.use("Agg")
import pylab
from matplotlib import gridspec
import matplotlib.pyplot as plt
import time
sys.path.append("../python/")
from  read_lgal_advance import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
matplotlib.rc('text', usetex=True)
matplotlib.rc('lines', linewidth=2)
matplotlib.rcParams['ytick.major.size'] = 8
matplotlib.rcParams['xtick.major.size'] = 8
def find_nearest(array,value,x):
    idx = numpy.searchsorted(array, value, side="left")
    v = (x[idx+1]-x[idx])*(value-array[idx])/(array[idx+1]-array[idx]) + x[idx]
    return v
def main(argv):
    gadget_m_conv = 1.e10
    hubble_h = 0.7
    a_list = numpy.loadtxt(a_list_file)
    num_bin = 100
    bin_size = 2.0/num_bin
    hist_x = numpy.linspace(-1.0, 1.0, num=num_bin, endpoint=False)
    for i in range(num_bin-1):
        hist_x[i] += bin_size*0.5
    hist_y = numpy.zeros(num_bin,dtype=numpy.int64)
    numfiles = lastfile - firstfile + 1
    if rank == 0:
        pool = numpy.arange(numfiles,dtype=numpy.int32)
        pool_list = numpy.array_split(pool,size)
        t_hist_y = numpy.zeros_like(hist_y)
    else:
        t_hist_y = None
        pool_list = None
    pool_list = comm.bcast(pool_list, root=0)
    for ifile in range(numfiles):
        if ifile in pool_list[rank]:
            (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,ifile,ifile,verbose=2)
            rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
            lastleafindex = numpy.cumsum(nTreeHalos)
            for i in range(len(nTreeHalos)):
                halos = output_Halos[rootindex[i]:lastleafindex[i]]
                haloIDs = output_HaloIDs[rootindex[i]:lastleafindex[i]]
                for this_h in numpy.where((halos['SnapNum'] == lastsnap))[0]:
                    c_mass = halos[this_h]['M_Crit200']*gadget_m_conv
                    if (haloIDs[this_h]["HaloID"] == haloIDs[this_h]["FirstHaloInFOFgroup"]) & (c_mass == 0.0):
                        c_mass = halos[this_h]["Len"]*mass_part*gadget_m_conv
                    c_prog = halos[this_h]['FirstProgenitor']
                    c_a = a_list[halos[this_h]["SnapNum"]]
                    while c_prog > -1:
                        p_mass = halos[c_prog]['M_Crit200']*gadget_m_conv
                        if (haloIDs[c_prog]["HaloID"] == haloIDs[c_prog]["FirstHaloInFOFgroup"]) & (p_mass == 0.0):
                            p_mass = halos[c_prog]["Len"]*mass_part*gadget_m_conv
                        p_prog = halos[c_prog]['FirstProgenitor']
                        #print "%g, %g"%(p_mass,c_mass)
                        p_a = a_list[halos[c_prog]["SnapNum"]]
                        if p_prog > -1:
                            pp_mass = halos[p_prog]['M_Crit200']*gadget_m_conv
                            if (haloIDs[p_prog]["HaloID"] == haloIDs[p_prog]["FirstHaloInFOFgroup"]) & (pp_mass == 0.0):
                                pp_mass = halos[p_prog]["Len"]*mass_part*gadget_m_conv
                            #print "%g, %g, %g"%(pp_mass,p_mass,c_mass)
                            pp_a = a_list[halos[p_prog]["SnapNum"]]
                            if(c_mass > 1.e8/hubble_h) & (p_mass > 1.e8/hubble_h) & (pp_mass > 1.e8/hubble_h):
                                alpha_a = numpy.log10(c_mass/p_mass)/numpy.log10(c_a/p_a)
                                alpha_b = numpy.log10(c_mass/pp_mass)/numpy.log10(c_a/pp_a)
                                epsilon = (numpy.arctan(alpha_b)-numpy.arctan(alpha_a))/numpy.pi
                                #print epsilon
                                slot = (epsilon + 1.0)/bin_size
                                hist_y[slot] += 1
                        c_prog = p_prog
                        c_mass = p_mass
                        c_time = 0
    comm.Barrier()
    sys.stdout.flush()
    comm.Reduce(
        [hist_y, MPI.LONG],
        [t_hist_y, MPI.LONG],
        op = MPI.SUM,
        root = 0
    )
    comm.Barrier()
    if rank == 0:
        c = 1.0*numpy.cumsum(t_hist_y,dtype=numpy.int64)/numpy.sum(t_hist_y,dtype=numpy.int64)
        a05 = find_nearest(c,0.05,hist_x)
        a50 = find_nearest(c,0.50,hist_x)
        a95 = find_nearest(c,0.95,hist_x)

        fig = plt.figure(figsize=(8, 6)) 
        gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        plt.subplots_adjust(hspace = 0)
        
        ax1.plot(hist_x,t_hist_y)
        ax0.plot(hist_x,c)
        ax0.grid()
        ax1.grid()
        ax1.axvline(x = a05,color='k',ls='dashed')
        ax1.axvline(x = a50,color='k',ls='dashed')
        ax1.axvline(x = a95,color='k',ls='dashed')
        
        ax1.set_xlabel(r"$\epsilon$")
        ax0.set_ylabel(r"$P( \le \epsilon)$")
        ax1.set_ylabel(r"$dN/d\epsilon(\epsilon)$")

        ax1.set_yscale("log")
        fig.savefig("fluc.pdf",bbox_inches='tight',pad_inches=0.1)
        plt.close(fig)
    comm.Barrier()
    return 0

if __name__ == "__main__":
    from folder import *
    main(sys.argv)
  
