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
npart = 1728
matplotlib.rc('text', usetex=True)
matplotlib.rc('savefig', dpi=300)
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
def main():
    filter = read_lgal.tree_properties_used
    filter['M_Crit200'] = True
    filter['Len'] = True
    filter['SnapNum'] = True
    filter['FirstHaloInFOFgroup'] = True
    folder = "/mnt/lustre/scratch/cs390/47Mpc/treedata/"
    file_prefix = "trees_075."
    firstfile = 0
    lastfile = 127
    f = open("clumping.dat","w")
    (nTrees,nHalos,nTreeHalos,output_Halos) = read_lgal.read_lgal_input_tree(folder,file_prefix,firstfile,lastfile,filter,verbose=True)
   
    countroot = 0
    for i in range(nTrees):
        for j in range(nTreeHalos[i]):
            output_Halos[countroot+j]["M_Crit200"] = output_Halos[countroot+output_Halos[countroot+j]["FirstHaloInFOFgroup"]]["M_Crit200"]
        countroot += nTreeHalos[i]

    for i in range(len(zlist)):
        print 'halo mass function: z = '+zlist[i]
        halos = output_Halos[numpy.where(output_Halos['SnapNum'] == i)]
        if(len(halos) > 0):
            high = numpy.sum(halos[numpy.where(halos['M_Crit200'] >= 0.1/hubble_h)]['Len'],dtype=numpy.int64)
            high = float(high)/npart**3
            low = numpy.sum(halos[numpy.where((halos['M_Crit200'] < 0.1/hubble_h) & (halos['M_Crit200'] >= 0.01/hubble_h))]['Len'],dtype=numpy.int64)
            low = float(low)/npart**3
            mid = numpy.sum(halos[numpy.where(halos['M_Crit200'] < 0.01/hubble_h)]['Len'],dtype=numpy.int64)
            mid = float(mid)/npart**3
            print zlist[i].strip(),low,mid,high
            print >> f,zlist[i].strip(),low,mid,high
    f.close()
    return 0

if __name__=="__main__":
    main()
    
