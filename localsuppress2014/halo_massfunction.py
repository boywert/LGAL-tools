import matplotlib
matplotlib.use('Agg')
import numpy
import os
import sys
import matplotlib.pyplot as plt
sys.path.append("../python/")
import read_lgal_advance as read_lgal
import mass_fn

hubble_h = 0.7

matplotlib.rc('text', usetex=True)
matplotlib.rc('savefig', dpi=300)
def main(folder,file_prefix,zlistfile):
    # zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
    zlist = open(zlistfile,"r").readlines()
    filter = read_lgal.tree_properties_used
    filter['M_Crit200'] = True
    filter['SnapNum'] = True
    # folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    # file_prefix = "trees_075."
    firstfile = 0
    lastfile = 127
    os.system("mkdir -p pdf")
    (nTrees,nHalos,nTreeHalos,output_Halos) = read_lgal.read_lgal_input_tree(folder,file_prefix,firstfile,lastfile,filter,verbose=True)
    for i in range(len(zlist)):
        print 'halo mass function: z = '+zlist[i]
        halos = output_Halos[numpy.where(output_Halos['SnapNum'] == i)]
        if(len(halos) > 0):
            (massfn_x,massfn_y) = mass_fn.M200c_mass_fn(halos,mass_min=1e8,mass_max=1.e13,nbins=50)
            f = open(zlist[i].strip()+"_massfn.dat","w")
            for j in range(len(massfn_x)):
                print >> f,massfn_x[j],massfn_y[j]
            f.close()
            # if(numpy.sum(massfn_y) > 0.00001):
            #     fig = plt.figure()
            #     ax = fig.add_subplot(111)
            #     ax.plot(massfn_x,massfn_y)
            #     ax.set_xlabel(r"$\mathrm{M_{200c}[}h^{-1}\mathrm{M_\odot]}$")
            #     ax.set_ylabel(r"$\mathrm{\Phi[}h^3\mathrm{Mpc^{-3} dex^-1]}$")
            #     ax.set_yscale("log")            
            #     fig.savefig("pdf/"+zlist[i].strip()+"_massfn.pdf",bbox_inches='tight',pad_inches=0)
    return 0

if __name__=="__main__":
    folder = sys.argv[1]
    prefix = sys.argv[2]
    zlist = sys.argv[3]
    main(folder,prefix,zlist)
    
