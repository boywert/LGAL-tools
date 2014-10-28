from mass_fn import *

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
    filter['Sfr'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)



z = sys.argv[1]
file_prefix = "SA_z"+z
firstfile = 0
lastfile = 127


config = {}
model_names = ["okamoto","noreionization","patchy_I"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_001/sams/5500.00/LGalaxyStruct.py"]
dt = []
filter = []
for i in range(len(struct_file)):
    filter.append(loadfilter(struct_file[i])[0])
    dt.append(loadfilter(struct_file[i])[1])
    print len(filter[i])

model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/okamoto/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_001/sams/5500.00/"]
model_plot_patterns = ['r--','g--','b--']

try:
    gal
except NameError:
    gal = {}
    nTrees = {}
    nGals = {}
    nTreeGals = {}
    star = {}
    #hotgas = {}
    #coldgas = {}
    #Blackhole = {}
    sfr = {}

for i in range(len(model_names)):
    index = model_names[i]
    if not index in gal:
        (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i])
        star[index] = stellar_mass_fn(gal[index],1.,1.e10,50)
        sfr[index] = sfr_fn(gal[index])



pylab.rc('text', usetex=True)

fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(star[index][0],star[index][1],model_plot_patterns[i],label=model_labels[i])
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(M/M_\odot h)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Stellar Mass Function z = "+z+" file "+str(firstfile)+"-"+str(lastfile))


pylab.savefig('reion_star_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')


fig = pylab.figure()
ax = fig.add_subplot(111)

for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(sfr[index][0],sfr[index][1],model_plot_patterns[i],label=model_labels[i])
ax.set_yscale("log")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)

ax.set_xlabel(r"$\log(Sfr/M_\odot yrs)$")
ax.set_ylabel(r"$N$")
fig.suptitle("Stellar Mass Function z = "+z+" file "+str(firstfile)+"-"+str(lastfile))

pylab.savefig('reion_sfr_'+str(firstfile)+'-'+str(lastfile)+'_'+file_prefix+'.pdf',bbox_inches='tight')




# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(hotgas[index][0],hotgas[index][1],model_plot_patterns[i],label=model_labels[i])

# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)

# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Hot Gass Mass Function z = 6 file "+str(firstfile)+"-"+str(lastfile))
# pylab.savefig('reion_hotgas_'+str(firstfile)+'-'+str(lastfile)+"_"+file_prefix+'.pdf',bbox_inches='tight')


# fig = pylab.figure()
# ax = fig.add_subplot(111)

# for i in range(len(model_names)):
#     index = model_names[i]
#     ax.plot(coldgas[index][0],coldgas[index][1],model_plot_patterns[i],label=model_labels[i])

# ax.set_yscale("log")
# leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
# leg.get_frame().set_linewidth(0)
# ax.set_xlabel(r"$\log(M/M_\odot h)$")
# ax.set_ylabel(r"$N$")
# fig.suptitle("Cold Gass Mass Function z = 6 file "+str(firstfile)+"-"+str(lastfile))

# pylab.savefig('reion_coldgas_'+str(firstfile)+'-'+str(lastfile)+'.pdf',bbox_inches='tight')


#fig = pylab.figure()
#ax = fig.add_subplot(111)

#for i in range(len(model_names)):
#    index = model_names[i]
#    ax.plot(Blackhole[index][0],Blackhole[index][1],model_plot_patterns[i],label=model_labels[i])

#leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
#leg.get_frame().set_linewidth(0)
#ax.set_xlabel(r"$\log(M/M_\odot h)$")
#ax.set_ylabel(r"$N$")
#fig.suptitle("Blackhole Mass Function z = 6 file "+str(ff)+"-"+str(lf))

#pylab.savefig('reion_bh_'+str(ff)+'-'+str(lf)+'.pdf',bbox_inches='tight')


