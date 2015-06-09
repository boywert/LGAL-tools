from mass_fn import *
import matplotlib
matplotlib.use('pdf') 
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
rank = "0"
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    filter['Type'] = True
    filter['Mag'] = True
    filter['Sfr'] = True
    filter['HaloM_Crit200'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    filter['HotGas'] = True
    filter['ColdGas'] = True
    filter['EjectedMass'] = True
    filter['MetalsDiskMass'] = True
    filter['MetalsBulgeMass'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)


h0 = 0.7
gadgetmass = 1.e10
model_names = ["noreionization","patchy8","patchy9","patchyG"]
struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_7/sams/14000.00/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_8/sams/14000.00/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5/sams/14000.00/LGalaxyStruct.py"]

use_model = [True,True,True,True]
dt = []
filter = []
for i in range(len(struct_file)):
    (f,t) = loadfilter(struct_file[i])
    filter.append(f)
    dt.append(t)

model_labels = ["No Reionization","Patchy Reionization (cutoff at 8)","Patchy Reionization (cutoff at 9)","Patchy Reionization (Gradual)"]
model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_7/sams/14000.00/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_8/sams/14000.00/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5/sams/14000.00/"]
model_plot_patterns = ['r--','g--','b--','y--']


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
zlistfile = "/mnt/lustre/scratch/cs390/47Mpc/snap_z.txt"
zlist = open(zlistfile,"r").readlines()
offset = 18.0


def plot_uv_z8():
    z = "7.96"
    file_prefix = "SA_z"+z
    firstfile = 0
    lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    uvlf_x = {}
    uvlf_y = {}
    luvlf_x = {}
    luvlf_y = {}
    sfr_x = {}
    sfr_y = {}
    smf_x = {}
    smf_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        logf = -2.5*numpy.log10(gal[index]["Sfr"])
        a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
        uvlf_x[index] = a[1][0:len(a[1])-1]+0.25-offset
        uvlf_y[index] = a[0]/47.**3/0.5
        (sfr_x[index],sfr_y[index]) =  sfr_density_fn(gal[index],mass_min=10**-0.5,mass_max=10.**3,nbins=10)
        (smf_x[index],smf_y[index]) =  stellar_mass_fn(gal[index],mass_min=10**7,mass_max=10.**12,nbins=50)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-23.,max=-17,nbins=12)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        

    # metals
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(metalicity_x[index],metalicity_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"Z")
    ax.set_ylabel(r"$\mathrm{\Phi (Mpc^{-3} Mag^-1)}$")
    ax.set_yscale("log")
    fig.savefig("metal_z8.pdf",bbox_inches='tight',pad_inches=0)
    
    # UVLF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_uv_z8("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(uvlf_x[index],uvlf_y[index],model_plot_patterns[i],label=model_labels[i])
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi (Mpc^{-3} Mag^-1)}$")
    ax.set_yscale("log")
    fig.savefig("uv_l_z8.pdf",bbox_inches='tight',pad_inches=0)

    # # SFR
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # add_observations.add_obs_sfr_z8("observations/SFR/",ax)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(sfr_x[index],sfr_y[index],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    # ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    # ax.set_yscale("log")
    # fig.savefig("sfr_z8.pdf",bbox_inches='tight',pad_inches=0)

    # # SMF
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # add_observations.add_obs_smf_z8("observations/SMF/",ax)
    # for i in range(len(model_names)):
    #     index = model_names[i]
    #     ax.plot(smf_x[index],smf_y[index],model_plot_patterns[i],label=model_labels[i])
    # leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    # leg.get_frame().set_linewidth(0)
    # ax.set_xlabel(r"$\mathrm{\log_{10} (M/M_\odot)}$")
    # ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    # ax.set_yscale("log")
    # fig.savefig("uv_l_z8.pdf",bbox_inches='tight',pad_inches=0)


def plot_uv_z7():
    z = "6.98"
    file_prefix = "SA_z"+z
    firstfile = 0
    lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    uvlf_x = {}
    uvlf_y = {}
    luvlf_x = {}
    luvlf_y = {}
    sfr_x = {}
    sfr_y = {}
    smf_x = {}
    smf_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        logf = -2.5*numpy.log10(gal[index]["Sfr"])
        a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
        uvlf_x[index] = a[1][0:len(a[1])-1]+0.25-offset
        uvlf_y[index] = a[0]/47.**3/0.5
        (sfr_x[index],sfr_y[index]) =  sfr_density_fn(gal[index],mass_min=10**-0.5,mass_max=10.**3,nbins=10)
        (smf_x[index],smf_y[index]) =  stellar_mass_fn(gal[index],mass_min=10**7,mass_max=10.**12,nbins=50)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-23.,max=-17,nbins=12)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        

    # metals
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(metalicity_x[index],metalicity_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"Z")
    ax.set_ylabel(r"$\mathrm{\Phi (Mpc^{-3} Mag^-1)}$")
    ax.set_yscale("log")
    fig.savefig("metal_z7.pdf",bbox_inches='tight',pad_inches=0)
    
    # UVLF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_uv_z7("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(uvlf_x[index],uvlf_y[index],model_plot_patterns[i],label=model_labels[i])
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi (Mpc^{-3} Mag^-1)}$")
    ax.set_yscale("log")
    fig.savefig("uv_l_z7.pdf",bbox_inches='tight',pad_inches=0)

    # SFR
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_sfr_z7("observations/SFR/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("sfr_z7.pdf",bbox_inches='tight',pad_inches=0)

    # SMF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_smf_z7("observations/SMF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(smf_x[index],smf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} (M/M_\odot)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("smf_z7.pdf",bbox_inches='tight',pad_inches=0)
    
def plot_uv_z6():
    z = "6.06"
    file_prefix = "SA_z"+z
    firstfile = 0
    lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    uvlf_x = {}
    uvlf_y = {}
    luvlf_x = {}
    luvlf_y = {}
    sfr_x = {}
    sfr_y = {}
    smf_x = {}
    smf_y = {}
    hmf_x = {}
    hmf_y = {}
    coldmf_x = {}
    coldmf_y = {}
    hotmf_x = {}
    hotmf_y = {}
    ejectedmf_x = {}
    ejectedmf_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        logf = -2.5*numpy.log10(gal[index]["Sfr"])
        a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
        uvlf_x[index] = a[1][0:len(a[1])-1]+0.25-offset-5*numpy.log10(0.7)
        uvlf_y[index] = a[0]/47.**3/0.5
        (sfr_x[index],sfr_y[index]) =  sfr_density_fn(gal[index],mass_min=10.**-0.5,mass_max=1.e3,nbins=10)
        (smf_x[index],smf_y[index]) =  stellar_mass_fn(gal[index],mass_min=1.e7,mass_max=1.e12,nbins=50)
        (hmf_x[index],hmf_y[index]) =  M200c_mass_fn_gal(gal[index],mass_min=1.e7,mass_max=1.e12,nbins=50)
        (coldmf_x[index],coldmf_y[index]) =  coldgas_mass_fn(gal[index],mass_min=1.e7,mass_max=1.e12,nbins=50)
        (hotmf_x[index],hotmf_y[index]) =  hotgas_mass_fn(gal[index],mass_min=1.e7,mass_max=1.e12,nbins=50)
        (ejectedmf_x[index],ejectedmf_y[index]) =  ejected_mass_fn(gal[index],mass_min=1.e7,mass_max=1.e12,nbins=50)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-23.,max=-17,nbins=12)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        

    # metals
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(metalicity_x[index],metalicity_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"Z")
    ax.set_ylabel(r"$\mathrm{\Phi Mpc^{-3} Mag^-1}$")
    ax.set_yscale("log")
    fig.savefig("metal_z6.pdf",bbox_inches='tight',pad_inches=0)
    
    # UVLF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_uv_z6("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        #ax.plot(uvlf_x[index],uvlf_y[index],model_plot_patterns[i],label=model_labels[i])
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"numbers $\mathrm{Mpc^{-3} Mag^-1}$")
    ax.set_yscale("log")
    fig.savefig("uv_l_z6.pdf",bbox_inches='tight',pad_inches=0)

    # SFR
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_sfr_z6("observations/SFR/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("sfr_z6.pdf",bbox_inches='tight',pad_inches=0)

    # SMF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_smf_z6("observations/SMF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(smf_x[index],smf_y[index],model_plot_patterns[i],label=model_labels[i])
        #ax.plot(hmf_x[index],hmf_y[index],model_plot_patterns[i],label=model_labels[i])
        #ax.plot(hotmf_x[index],hotmf_y[index],"b-",label=model_labels[i])
        #ax.plot(coldmf_x[index],coldmf_y[index],"k-",label=model_labels[i])
        #ax.plot(ejectedmf_x[index],ejectedmf_y[index],"r-",label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} (M/M_\odot)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("smf_z6.pdf",bbox_inches='tight',pad_inches=0)


def plot_uv_z12():
    z = "12.05"
    file_prefix = "SA_z"+z
    firstfile = 0
    lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    uvlf_x = {}
    uvlf_y = {}
    luvlf_x = {}
    luvlf_y = {}
    sfr_x = {}
    sfr_y = {}

    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)

        logf = -2.5*numpy.log10(gal[index]["Sfr"])
        a = numpy.histogram(logf,bins=9,range=(-3.0,1.5))
        uvlf_x[index] = a[1][0:len(a[1])-1]+0.25-offset-5*numpy.log10(0.7)
        uvlf_y[index] = a[0]/47.**3/0.5
        (sfr_x[index],sfr_y[index]) =  sfr_density_fn(gal[index],mass_min=10.**-0.5,mass_max=1.e3,nbins=10)

    # SFR
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(sfr_x[index],sfr_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"$\mathrm{\log_{10} SFR(M_\odot/year)}$")
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    ax.set_yscale("log")
    fig.savefig("sfr_z12.pdf",bbox_inches='tight',pad_inches=0)

    
def main():
    #plot_uv_z6()
    #plot_uv_z7()
    #plot_uv_z8()
    plot_uv_z12()
if __name__=="__main__":
    main()
