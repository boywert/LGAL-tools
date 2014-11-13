from mass_fn import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct

sys.path.append("../python/")
import read_lgal_advance as read_lgal

def loadfilter(structfile):
    sys.path.insert(0,"../tmp/")
    os.system("cp "+structfile+" ../tmp/LGalaxyStruct.py")
    os.system("rm -f ../tmp/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for index in filter.keys():
        filter[index] = True

    filter['Type'] = True
    filter['HaloIndex'] = True
    filter['Sfr'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    filter['ICM'] = True
    filter['HotGas'] = True
    filter['ColdGas'] = True
    filter['sfh_numbins'] = True
    filter['sfh_BulgeMass'] = True
    filter['sfh_DiskMass'] = True
    filter['GasDiskRadius'] = True
    filter['CoolingRate'] = True
    filter['EjectedMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)



#z = sys.argv[1]
z = '6.00'

file_prefix = "SA_z"+z
file_prefix = "SA_"
firstfile = 0
lastfile = 0
config = {}
#model_names = ["okamoto","noreionization","patchy_I"]
model_names = ["okamoto","patchy_I"]
#struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/LGalaxyStruct.py"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/lgaltree/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/LGalaxyStruct.py"]
dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

#model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
model_labels = ["Okamoto et al. (2008)","Patchy Reionization (Gradual)"]
#model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/lgaltree/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/lgaltree","/mnt/lustre/scratch/cs390/47Mpc/couple/model_002/fullgaltree/43000.00/"]
model_plot_patterns = ['r--',"g--"]



try:
    gal
except NameError:
    gal = {}
    nTrees = {}
    nGals = {}
    nTreeGals = {}
    star = {}
    totsfr = {}
    #hotgas = {}
    #coldgas = {}
    #Blackhole = {}
    sfr = {}


for i in range(len(model_names)):
    index = model_names[i]
    if not index in gal:
        (nGals[index],gal[index]) = read_lgal.read_lgaltree_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        # totsfr[index] = numpy.sum(gal[index]["Sfr"],dtype=numpy.float64)
        # star[index] = stellar_mass_fn(gal[index],1.,1.e10,50)
        # sfr[index] = sfr_fn(gal[index])

gal_type0 = {}
for i in range(len(model_names)):
    index = model_names[i]
    gal_type0[index] = gal[index][numpy.where((gal[index]["Type"] == 0) & (gal[index]["SnapNum"] == 75))[0]]

cmp_sfr = {}

for i in range(len(model_names)):
    index = model_names[i]
    cmp_sfr[index] = []
    haloid = 0

    id = numpy.where(gal[index]["HaloID"] == haloid)[0][0]
    nextid = id
    while nextid > -1:
        #print id
        #print index,gal[index][id]["Sfr"],gal[index][id]["ColdGas"],gal[index][id]["HotGas"],gal[index][id]["EjectedMass"] 
        cmp_sfr[index].append([gal[index][id]["SnapNum"],gal[index][id]["Sfr"],gal[index][id]["ColdGas"],gal[index][id]["HotGas"],gal[index][id]["BulgeMass"]+gal[index][id]["DiskMass"]])
        nextgalid = gal[index][id]["FirstProgGal"]
        nextid = numpy.where(gal[index]["GalID"] == nextgalid)[0]
        if len(nextid) > 0:
            id = nextid[0]
        else:
            id = -1
    
print cmp_sfr
for i in range(len(cmp_sfr['okamoto'])):
    print i,cmp_sfr['okamoto'][i],cmp_sfr['patchy_I'][i]

fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
   index = model_names[i]
   ax.plot(cmp_sfr[index][:][0],cmp_sfr[index][:][1],model_plot_patterns[i],label=model_labels[i])
ax.set_yscale("log")
#ax.set_xscale("log")
#leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
#leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(SFR)$")
ax.set_ylabel(r"$\log(SFR)$")
fig.suptitle("SFR file "+str(firstfile)+"-"+str(lastfile))

#fig.show()
fig.savefig('cmp_sfr_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')
#pylab.close(fig)

# for i in range(len(model_names)):
#     index = model_names[i]
#     print nTrees[index],nGals[index],totsfr[index]

# pylab.rc('text', usetex=True)

# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(star[index][0],star[index][1],model_plot_patterns[i],label=model_labels[i])
# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)

# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Stellar Mass Function z = "+z+" file "+str(firstfile)+"-"+str(lastfile))


# fig.savefig('reion_star_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')
# pylab.close(fig)

# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(sfr[index][0],sfr[index][1],model_plot_patterns[i],label=model_labels[i])
# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)

# ax.set_xlabel(r"$\log(Sfr/M_\odot yrs)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Stellar formation rate z = "+z+" file "+str(firstfile)+"-"+str(lastfile))

# fig.savefig('reion_sfr_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')
# pylab.close(fig)



# gal_type0 = {}
# gal_type1 = {}
# gal_type2 = {}

# for i in range(len(model_names)):
#     index = model_names[i]
#     gal_type0[index] = gal[index][numpy.where(gal[index]["Type"] == 0)[0]]
#     gal_type1[index] = gal[index][numpy.where(gal[index]["Type"] == 1)[0]]
#     gal_type2[index] = gal[index][numpy.where(gal[index]["Type"] == 2)[0]]


# # star_type0 = {}
# for i in range(len(model_names)):
#     index = model_names[i]
#     star_type0[index] = stellar_mass_fn(gal_type0[index],1.,1.e10,50)

# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(star_type0[index][0],star_type0[index][1],model_plot_patterns[i],label=model_labels[i])
# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Stellar Mass Function z = "+z+" file "+str(firstfile)+"-"+str(lastfile))
# fig.savefig('star_type0_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')
# pylab.close(fig)

    
# # histogram tree by tree since the number of trees must be identcal
# try:
#     sfr_tree
# except NameError:
#     sfr_tree = {}

# for i in range(len(model_names)):
#     index = model_names[i]
#     if not index in sfr_tree:
#         cumsumntrees = numpy.cumsum(nTreeGals[index])
#         sfr_tree[index] = numpy.zeros(nTrees[index],dtype=numpy.float64)
#         for j in range(nTrees[index]):
#             sfr_tree[index][j] = numpy.sum(gal[index]["Sfr"][cumsumntrees[j]:cumsumntrees[j]+nTreeGals[index][j]],dtype=numpy.float64)
#         a = numpy.where(sfr_tree[index] > 0.)[0]
#         print len(a)
#         print sfr_tree[index][a]

# for i in range(len(model_names)): 
#     model_i = model_names[i]
#     for j in range(i+1,len(model_names)):
#         model_j = model_names[j]
#         minval = min(min(sfr_tree[model_i]),min(sfr_tree[model_j]))
#         maxval = min(max(sfr_tree[model_i]),max(sfr_tree[model_j]))
#         x = numpy.linspace(0.0,maxval,num=5)
#         fig = pylab.figure()
#         ax = fig.add_subplot(111)
#         ax.scatter(sfr_tree[model_i],sfr_tree[model_j], s=1.0)
#         ax.plot(x,x,"k-")
#         #ax.set_yscale("log")
#         #ax.set_xscale("log")
#         fig.suptitle("SFR-SFR tree by tree, "+model_labels[i]+" vs "+model_labels[j])
#         fig.savefig('sfr_vs_sfr_'+model_names[i]+"_vs_"+model_names[j]+".pdf",bbox_inches='tight')
#         pylab.close(fig)
