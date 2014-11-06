import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import os
import matplotlib.pyplot as plt
import numpy
#z = sys.argv[1]

zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
folder = "sfr/"
model_names = ["okamoto","noreionization","patchy_I"]
model_labels = ["Okamoto et al. (2008)","No Reionization","Patchy Reionization (Gradual)"]
color = ["r","g","b"]
pattern = ["-","--","-."]
sfr_t0 = {}
sfr_t1 = {}
sfr_t2 = {}
for i in range(len(model_names)):
    index = model_names[i]
    sfr_t0[index] = []
    sfr_t1[index] = []
    sfr_t2[index] = []
j = 0
for z in zlist:
    z = z.strip()
    zlist[j] = float(z)
    file = folder+"/"+z+".dat"
    data = numpy.loadtxt(file)
    for i in range(len(model_names)):
        index = model_names[i]
        sfr_t0[index].append(data[i][0])
        sfr_t1[index].append(data[i][1])
        sfr_t2[index].append(data[i][2])
    j += 1


fig = pylab.figure()
ax = fig.add_subplot(111)
for i in range(len(model_names)):
    index = model_names[i]
    ax.plot(zlist,sfr_t0[index],color[i]+pattern[0],label=model_labels[i]+",Type 0")
    #ax.plot(zlist,sfr_t1[index],color[i]+pattern[1],label=model_labels[i]+",Type 1")
    #ax.plot(zlist,sfr_t2[index],color[i]+pattern[2],label=model_labels[i]+",Type 2")
leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
#ax.set_yscale("log")
fig.suptitle("SFR")
fig.savefig("SFR_history.pdf",bbox_inches='tight')
pylab.close(fig)
