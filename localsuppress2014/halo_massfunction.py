import matplotlib
matplotlib.use('Agg')
import numpy
import os
import matplotlib.pyplot as plt
sys.path.append("../python/")
import read_lgal_advance as read_lgal
import mass_fn

hubble_h = 0.7

matplotlib.rc('text', usetex=True)
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
def main():
    filter = read_lgal.tree_properties_used
    filter['M_Crit200'] = True
    filter['SnapNum'] = True
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_075."
    firstfile = 0
    firstfile = 127
    (nTrees,nHalos,nTreeHalos,output_Halos) = read_lgal_input_tree(folder,file_prefix,firstfile,lastfile,filter,verbose=False)
    for i in range(len(zlist)):
        halos = output_Halos[numpy.where(output_Halos['SnapNum'] == i)]
        (massftn_x,massftn_y) = M200c_mass_fn(halos,mass_min=1e8,mass_max=1.e13,nbins=50)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(massfn_x,massfn_y)
        ax.set_xlabel(r"$M_{200c}(h^{-1}M_\odot)$")
        ax.set_ylabel(r"numbers $\mathrm{Mpc^{-3} dex^-1}$")
        ax.set_yscale("log")
        fig.savefig(zlist[i]+".pdf",bbox_inches='tight',pad_inches=0)
    return 0

if __name__=="__main__":
    main()
    
