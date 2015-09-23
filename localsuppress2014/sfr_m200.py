from mass_fn import *
from globalconf import *
import matplotlib
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import os
import matplotlib.pyplot as plt
os.system("cp dummy_dtype.py LGalaxyStruct.py")
import LGalaxyStruct
import add_observations
sys.path.append("../python/")
import read_lgal_advance as read_lgal
import timeit
rank = "0"
SEC_PER_YEAR = 3600*24*365.25
Msun2kg = 1.989e30
h_mass = 1.6737237e-27 #kg
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['NPhotReion'] = True
    filter['HaloM_Crit200'] = True
    filter['Sfr'] = True
    filter['Type'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

#filter model
filter_tmp = []
dt_tmp = []
model_names_tmp = []
struct_file_tmp = []
model_labels_tmp = []
model_paths_tmp = []
for i in range(len(use_model)):
    if use_model[i]:
        filter_tmp.append(filter[i])
        dt_tmp.append(dt[i])
        model_names_tmp.append(model_names[i])
        struct_file_tmp.append(struct_file[i])
        model_labels_tmp.append(model_labels[i])
        model_paths_tmp.append(model_paths[i])

filter = filter_tmp
dt = dt_tmp
model_names = model_names_tmp
struct_file = struct_file_tmp
model_labels = model_labels_tmp
model_paths = model_paths_tmp       


pylab.rc('text', usetex=True)

zlist = open(zlistfile,"r").readlines()


def plot_uv_z8():
    z = "6.00"
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sum_logphoton = {}
    sum_SFR = {}
    sum_SFR_sq = {}
    N = {}
    mean_SFR = {}
    mean_SFR_sq = {}
    mean_logphoton = {}
    m200c = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        rangen = (7.5,11.5)
        bins = 40
        gal[index] = gal[index][numpy.where((gal[index]["HaloM_Crit200"]>0.))]
        nummax= numpy.nanmax(gal[index]["NPhotReion"])
        gal[index]["NPhotReion"] = numpy.clip(gal[index]["NPhotReion"]+numpy.log10(SEC_PER_YEAR),0.0,nummax)
        print numpy.log10(gal[index]["Sfr"].astype(numpy.float64)*Msun2kg/h_mass)
        sum_logphoton[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["NPhotReion"]-numpy.log10(gal[index]["Sfr"].astype(numpy.float64)*Msun2kg/h_mass))
        ssfr = gal[index]["Sfr"]/(gal[index]["HaloM_Crit200"]*1.e10/hubble_h)
        ssfr = numpy.nan_to_num(ssfr)
        sum_SFR[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=ssfr)
        sum_SFR_sq[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["Sfr"]**2)
        N[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins)
        mean_SFR[index] = sum_SFR[index][0]/N[index][0]
        mean_logphoton[index] = sum_logphoton[index][0] /N[index][0]
        print sum_logphoton[index]
        mean_SFR_sq[index]= sum_SFR_sq[index][0]/N[index][0]
        m200c[index] = []
        for i in range(len(sum_SFR[index][0])):
            m200c[index].append(0.5*(sum_SFR[index][1][i]+sum_SFR[index][1][i+1]))
        del(gal[index])
        del(nTreeGals[index])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(m200c[index],mean_SFR[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$M_{200c}[M_\odot]$")
    ax.set_ylabel(r"$\mathrm{SFR [M_\odot/year]}$")
    ax.set_yscale("log")
    fig.savefig("SFRvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(m200c[index],mean_logphoton[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$M_{200c}[M_\odot]$")
    ax.set_ylabel(r"$\mathrm{NPHOT}$")
    fig.savefig("NPHOTvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
def main():
    #plot_uv_z6()
    #plot_uv_z7()
    plot_uv_z8()

if __name__=="__main__":
    main()
