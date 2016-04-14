import matplotlib as plt
import sys
from numpy import *
plt.use("Agg")
import pylab 
import time
sys.path.append("../python/")
from  read_lgal_advance import *
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
plt.rc('text', usetex=True)
plt.rc('lines', linewidth=2)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8

def main(argv):
    gadget_m_conv = 1.e10
    hubble_h = 0.7
    a_list = numpy.loadtxt(a_list_file)
    num_bin = 100
    bin_size = 2.0/num_bin
    hist_x = numpy.linspace(-1.0, 1.0, num=num_bin, endpoint=False)
    for i in range(num_bin-1):
        hist_x[i] = 0.5*(hist_x[i]+hist_x[i+1])
    hist_x[num_bin-1] = (hist_x[num_bin-1]+1.0)*0.5
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
        print t_hist_y
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        ax.plot(hist_x,t_hist_y)
        #leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
        #leg.get_frame().set_linewidth(0)
        ax.set_xlabel(r"$\pi^{-1}\epsilon$")
        ax.set_ylabel(r"$P(\epsilon)$")
        ax.set_yscale("log")
        fig.savefig("fluc.pdf",bbox_inches='tight',pad_inches=0.1)
        pylab.close(fig)
    comm.Barrier()
    return 0

if __name__ == "__main__":
    from folder import *
    main(sys.argv)
  
