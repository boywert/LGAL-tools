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
os.system("mkdir -p ../tmp/"+rank)
def loadfilter(structfile):
    sys.path.insert(0,"../tmp/"+rank)
    os.system("cp "+structfile+" ../tmp/"+rank+"/LGalaxyStruct.py")
    os.system("rm -f ../tmp/"+rank+"/LGalaxyStruct.pyc")
    reload(LGalaxyStruct)
    filter = LGalaxyStruct.properties_used
    for fi in filter:
        fi = False
    
    filter['Mag'] = True
    filter['Sfr'] = True
    filter['DiskMass'] = True
    filter['BulgeMass'] = True
    filter['EjectedMass'] = True
    filter['MetalsDiskMass'] = True
    filter['MetalsBulgeMass'] = True
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
    z = "7.96"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sfr_x = {}
    sfr_y = {}
    luvlf_x = {}
    luvlf_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (sfr_x[index],sfr_y[index]) = sfr_density_fn(gal[index],mass_min=10**-4.5,mass_max=10.**4,nbins=40)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-13,nbins=24)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        
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
        ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} dex^{-1}}]$")
        ax.set_yscale("log")
        fig.savefig("sfr_z8.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)
    
    # metals
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(metalicity_x[index],metalicity_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"Z")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^-1]}$")
    ax.set_yscale("log")
    fig.savefig("metal_z8.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)
    
    # UVLF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_uv_z8("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi [ Mpc^{-3} Mag^{-1}]$")
    ax.set_yscale("log")
    ax.set_xlim([-20.5,-16.])
    fig.savefig("uv_l_z8.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)


def plot_uv_z7():
    z = "6.98"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sfr_x = {}
    sfr_y = {}
    luvlf_x = {}
    luvlf_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (sfr_x[index],sfr_y[index]) = sfr_density_fn(gal[index],mass_min=10**-4.5,mass_max=10.**4,nbins=40)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-15,nbins=20)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        
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
    plt.close(fig)
                                                                                
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
    plt.close(fig)
    
    # UVLF
    fig = plt.figure()
    ax = fig.add_subplot(111)
    add_observations.add_obs_uv_z7("observations/UVLF/",ax)
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    ax.set_xlim([-21,-16])
    fig.savefig("uv_l_z7.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)


    
def plot_uv_z6():
    z = "6.00"
    file_prefix = "SA_z"+z
    #firstfile = 0
    #lastfile = 127
    config = {}

    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}

    luvlf_x = {}
    luvlf_y = {}
    sfr_x = {}
    sfr_y = {}
    metalicity_x = {}
    metalicity_y = {}
    for i in range(len(model_names)):
        index = model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],0)
        (sfr_x[index],sfr_y[index]) = sfr_density_fn(gal[index],mass_min=10**-4.5,mass_max=10.**4,nbins=40)
        (luvlf_x[index],luvlf_y[index]) = uv_luminosity_fn(gal[index],min=-25.,max=-15,nbins=20)
        (metalicity_x[index],metalicity_y[index]) = metallicity_fn(gal[index],mass_min=1.e-5,mass_max=1.,nbins=20)
        
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
    plt.close(fig)
                                                                                
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
        ax.plot(luvlf_x[index],luvlf_y[index],model_plot_patterns[i],label=model_labels[i])
    leg = ax.legend(loc='best', handlelength = 7,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel(r"M1600 - 5log(h)")
    ax.set_ylabel(r"$\mathrm{\Phi[Mpc^{-3} Mag^{-1}]}$")
    ax.set_yscale("log")
    ax.set_xlim([-22.,-15.5])
    fig.savefig("uv_l_z6.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)


    
def main():
    #plot_uv_z6()
    #plot_uv_z7()
    plot_uv_z8()

if __name__=="__main__":
    main()
