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
import random
rank = "0"
SEC_PER_YEAR = 3600*24*365.25
Msun2kg = 1.989e30
h_mass = 1.6737237e-27 #kg
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    ranki = str(random.randint(0,1000000))
    sys.path.insert(0,"../tmp/"+ranki)
    os.system("cp "+structfile+" ../tmp/"+ranki+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+ranki+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['NPhotReion'] = True
    filter['Mvir'] = True
    filter['HaloM_Crit200'] = True
    filter['HotGas'] = True
    filter['ColdGas'] = True
    filter['EjectedMass'] = True
    filter['StellarMass'] = True
    filter['ICM'] = True
    filter['BlackHoleGas'] = True
    filter['BlackHoleMass'] = True
    filter['Sfr'] = True
    filter['Type'] = True
    #filter['CumulativeSFR'] = True
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


def plot_z(z):
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sum_logphoton = {}
    sum_sSFR = {}
    sum_SFR = {}
    sum_hotgas = {}
    sum_coldgas= {}
    sum_SFR_sq = {}
    sum_stellarmass  = {}
    sum_stellarratio  = {}
    sum_ejectedmass = {}
    sum_ejectedratio = {}
    sum_baryons = {}
    sum_gammaratio = {}
    N = {}
    mean_SFR = {}
    mean_SFR_sq = {}
    sum_ASFR = {}
    mean_logphoton = {}
    m200c = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        rangen = (7.5,11.5)
        bins = 40
        gal[index] = gal[index][numpy.where((gal[index]["HaloM_Crit200"] >0.))]
        
        # nummax= numpy.nanmax(gal[index]["NPhotReion"])
        # gal[index]["NPhotReion"] = numpy.clip(gal[index]["NPhotReion"]+numpy.log10(SEC_PER_YEAR),0.0,nummax)
        # total_sfr =  numpy.log10(gal[index]["Sfr"].astype(numpy.float64)*Msun2kg/h_mass)
        # nummax2= numpy.nanmax(total_sfr)
        # total_sfr = numpy.clip(total_sfr,0.0,nummax2)
        # avg = numpy.sum(gal[index]["NPhotReion"] - total_sfr,dtype=numpy.float64)/len(total_sfr)
        # print index,"avg = ",10.**avg
        sum_baryons[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=(gal[index]["StellarMass"]+gal[index]["EjectedMass"]+gal[index]["ColdGas"]+gal[index]['HotGas']+gal[index]["ICM"]+gal[index]["BlackHoleMass"])/gal[index]["HaloM_Crit200"]/0.166)
        #sum_logphoton[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins,weights=numpy.float64(1)*10.**gal[index]["NPhotReion"].astype(numpy.float64) -1.)
        sum_logphoton[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["Sfr"].astype(numpy.float64)*Msun2kg/h_mass*11.5e6*2800)
        sum_hotgas[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["HotGas"].astype(numpy.float64)*1.e10 )
        sum_coldgas[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["ColdGas"].astype(numpy.float64)*1.e10)
        sum_stellarmass[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["StellarMass"].astype(numpy.float64)*1.e10)
        sum_stellarratio[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["StellarMass"].astype(numpy.float64)/gal[index]["HaloM_Crit200"])
        sum_sSFR[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["Sfr"]/(gal[index]["StellarMass"].astype(numpy.float64)/hubble_h*1.e10))
        sum_ejectedmass[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["EjectedMass"].astype(numpy.float64)*1.e10)
        sum_ejectedratio[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["EjectedMass"].astype(numpy.float64)/gal[index]["Mvir"])
        ssfr = gal[index]["Sfr"]/(gal[index]["HaloM_Crit200"]*1.e10/hubble_h)
        ssfr = numpy.nan_to_num(ssfr)
        sum_SFR[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=gal[index]["Sfr"]*11.6e6)
        #sum_ASFR[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins,weights=gal[index]["CumulativeSFR"].astype(numpy.float64)*1.e10/hubble_h*Msun2kg/h_mass*500.)
        #sum_ASFR[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins,weights=gal[index]["CumulativeSFR"].astype(numpy.float64)*500./gal[index]['Mvir'])
        sum_SFR_sq[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins,weights=gal[index]["Sfr"]**2)
        N[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins)
        mean_SFR[index] = sum_SFR[index][0]/N[index][0]
        mean_logphoton[index] = sum_logphoton[index][0]/N[index][0]
        mean_SFR_sq[index]= sum_SFR_sq[index][0]/N[index][0]
        m200c[index] = []
        for i in range(len(sum_SFR[index][0])):
            m200c[index].append(0.5*(sum_SFR[index][1][i]+sum_SFR[index][1][i+1]))
        del(gal[index])
        del(nTreeGals[index])
        #print sum_logphoton[index][0]
        #print sum_logphoton[index][0].dtype
        print sum_logphoton[index][0]*SEC_PER_YEAR/(sum_SFR[index][0].astype(numpy.float64)*Msun2kg/h_mass)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_baryons[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$BaryonMass/(M_{200c}\Omega_b\Omega_m^{-1})$")
    # #ax.set_yscale("log")
    # fig.savefig("baryonsratiovsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_logphoton[index][0]*SEC_PER_YEAR*11.6e6/(sum_SFR[index][0].astype(numpy.float64)*Msun2kg/h_mass),model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{sSFR[yr^{-1}]}$")
    # ax.set_yscale("log")
    # fig.savefig("gamma_peratomvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_ASFR[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{ASFR[]}$")
    # ax.set_yscale("log")
    # fig.savefig("ASFRvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_stellarmass[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{Stellar Mass[h^{-1}M_\odot]}$")
    # ax.set_yscale("log")
    # fig.savefig("StarvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_ejectedratio[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"Ejected Mass / M$_{200c}$")
    # ax.set_yscale("log")
    # fig.savefig("EjectRatiovsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(m200c[index],sum_stellarratio[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    ax.set_ylabel(r"Stellar Mass / M$_{200c}$")
    ax.set_yscale("log")
    ax.set_ylim([1.e-6,1.e-2])
    fig.savefig("StarRatiovsM_z"+"%08.2f"%(float(z))+".png",bbox_inches='tight')
    plt.close(fig)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_ejectedmass[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{Ejected Mass[h^{-1}M_\odot]}$")
    # ax.set_yscale("log")
    # fig.savefig("EjectedvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)


    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_coldgas[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{ColdGas[M_\odot/h}$")
    # ax.set_yscale("log")
    # fig.savefig("ColdgasvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],sum_hotgas[index][0]/N[index][0],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{HotGas[M_\odot/h}$")
    # ax.set_yscale("log")
    # fig.savefig("HotgasvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],mean_SFR[index],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{SFR [M_\odot/year]}$")
    # ax.set_yscale("log")
    # fig.savefig("SFRvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(m200c[index],mean_logphoton[index],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$M_{200c}[h^{-1}M_\odot]$")
    # ax.set_ylabel(r"$\mathrm{NPHOT}$")
    # ax.set_yscale("log")
    # fig.savefig("NPHOTvsM_z"+str(z)+".pdf",bbox_inches='tight',pad_inches=0)
    
def main():
    zlist = open(zlistfile).readlines()
    zi = zlist[long(sys.argv[1])].strip()
    plot_z(zi)
    #plot_z("7.96")

if __name__=="__main__":
    main()
