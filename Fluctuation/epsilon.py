import matplotlib as plt
import sys
from numpy import *
plt.use("Agg")

def main(argv):
    gadget_m_conv = 1.e10
    hubble_h = 0.7
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=True)
    rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
    lastleafindex = numpy.cumsum(nTreeHalos)
    for i in range(len(nTreeHalos)):
        halos = output_Halos[rootindex[i]:lastleafindex[i]]
        haloIDs = output_HaloIDs[rootindex[i]:lastleafindex[i]]
        for this_h in numpy.where((halos['SnapNum'] == lastsnap)):
            c_mass = halos[this_h]['M_Crit200']
            c_prog = halos[this_h]['FirstProgenitor']
            if c_prog > -1:
                p_mass = halos[c_prog]['M_Crit200']
                p_prog = halos[c_prog]['FirstProgenitor']
                if c_prog > -1:
                    pp_mass = halos[p_prog]['M_Crit200']

    return 0

if __name__ == "__main__":
    from folder import *
    main(sys.argv)
  
