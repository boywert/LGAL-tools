import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import os
import matplotlib.pyplot as plt
import numpy
#z = sys.argv[1]
h0=0.7
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
folder = "sfr/"
model_names = ["okamoto","noreionization","patchy_I"]
model_labels = ["AS","NS","LS"]
color = ["#4477AA","#DDCC77","#CC6677"]
pattern = ["-.","--","-."]
sfr_t0 = {}
sfr_t1 = {}
sfr_t2 = {}
sfr_lm = {}
sfr_hm = {}
for i in range(len(model_names)):
    index = model_names[i]
    sfr_t0[index] = []
    sfr_t1[index] = []
    sfr_t2[index] = []
    sfr_lm[index] = []
    sfr_hm[index] = []
j = 0
for z in zlist[0:len(zlist)-1]:
    z = z.strip()
    zlist[j] = float(z)
    file = folder+"/"+z+".dat"
    data = numpy.loadtxt(file)
    for i in range(len(model_names)):
        index = model_names[i]
        sfr_t0[index].append(data[i][0])
        sfr_t1[index].append(data[i][1])
        sfr_t2[index].append(data[i][2])
        sfr_lm[index].append(data[i][3])
        sfr_hm[index].append(data[i][4])
    j += 1

pylab.rc('text', usetex=True)
fig = pylab.figure()
ax = fig.add_subplot(111)
for i in range(len(model_names)):
    index = model_names[i]
    x = zlist[0:len(zlist)-1]
    y = numpy.array(sfr_lm[index])/(47./h0)**3
    #y = numpy.array(sfr_t0[index])+numpy.array(sfr_t1[index])+numpy.array(sfr_t2[index])
    ax.plot(x,numpy.log10(y),color=color[i],linestyle=pattern[0],label="LMACH - "+model_labels[i], linewidth=1.5)
    y = numpy.array(sfr_hm[index])/(47./h0)**3
    #y = numpy.array(sfr_t0[index])+numpy.array(sfr_t1[index])+numpy.array(sfr_t2[index])
    ax.plot(x,numpy.log10(y),color=color[i],linestyle=pattern[1],label="HMACH - "+model_labels[i], linewidth=1.5)
  

leg = ax.legend(loc='best', handlelength = 9,ncol=1, fancybox=True, prop={'size':10})
leg.get_frame().set_linewidth(0)
ax.set_ylabel(r"log$_{\mathrm{10}}$(SFRD)[$\mathrm{M_\odot yr^{-1} Mpc^{-3}}$]")
ax.set_xlabel(r"z")
fig.savefig("SFR_history.pdf",bbox_inches='tight')
pylab.close(fig)
