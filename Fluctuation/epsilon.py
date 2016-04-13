import matplotlib as plt
import sys
from numpy import *
plt.use("Agg")
from pylab import *
import time
sys.path.append("../python/")
from  read_lgal_advance import *

def main(argv):
    gadget_m_conv = 1.e10
    hubble_h = 0.7
    for ifile in range(1):
        (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,ifile,ifile,verbose=True)
        rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
        lastleafindex = numpy.cumsum(nTreeHalos)
        for i in range(len(nTreeHalos)):
            halos = output_Halos[rootindex[i]:lastleafindex[i]]
            haloIDs = output_HaloIDs[rootindex[i]:lastleafindex[i]]
            for this_h in numpy.where((halos['SnapNum'] == lastsnap))[0]:
                c_mass = halos[this_h]['M_Crit200']
                c_prog = halos[this_h]['FirstProgenitor']
                c_time = 0
                while c_prog > -1:
                    p_mass = halos[c_prog]['M_Crit200']
                    #alpha_1 = numpy.arctan(numpy.log10(c_mass/p_mass)/numpy.log10(c_time/p_time))
                    p_prog = halos[c_prog]['FirstProgenitor']
                    print c_prog,p_prog
                    if p_prog > -1:
                        pp_mass = halos[p_prog]['M_Crit200']
                        pp_time = 0
                    c_prog = p_prog
                    c_mass = p_mass
                    c_time = 0
    return 0

if __name__ == "__main__":
    from folder import *
    main(sys.argv)
  
