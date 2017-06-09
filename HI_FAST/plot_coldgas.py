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
    filter['ColdGas'] = True
    filter['Pos'] = True
    filter['Vel'] = True
    filter['Spin'] = True
    filter['Rvir'] = True
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
pylab.rc('lines', linewidth=2)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8
#zlist = open(zlistfile,"r").readlines()



def plot_coldgas(z,ax,sam):
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

    smf_x = {}
    smf_y = {}

    for i in range(len(model_names)):
        index = model_names[i]
        if index[:4] == "lgal":
            z = "%10.2f"%(z)
        elif index[:4] == "sage":
            z = "%10.3f"%(z)
        file_prefix = "model_z"+z.strip()
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
            (smf_x[index],smf_y[index]) = coldgas_mass_fn(gal[index],mass_min=1.e4,mass_max=1e11,nbins=40)
            
    for i in range(len(model_names)):
        index = model_names[i]
        ax.plot(smf_x[index],smf_y[index],color=model_plot_colors[i],linestyle=model_plot_patterns[i],label=model_labels[i])
    ax.set_xlabel(r"$\mathrm{\log_{10}[m_{\rm cold}/M_\odot]}$")
    ax.set_yscale("log")
    ax.set_ylim([1.e-2,1e1])
    ax.set_xlim([6,8])
    ax.text(0.1, 0.1, 'z = %d'%(int(float(z)+0.5)),
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=15)
    leg = ax.legend(loc='upper right', handlelength = 10,ncol=1, fancybox=True, prop={'size':10})
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_alpha(0)
    ax.set_ylabel(r"$\mathrm{\Phi(Mpc^{-3} dex^{-1}})$")
    labels = ["",r"$6.5$",r"$7.0$",r"$7.5$",r"$8.0$",r"$8.5$",r"$9.5$",r"10.0"]
    ax.xaxis.set_ticklabels(labels)
    ax.yaxis.set_ticklabels([])
    ax.set_axisbelow(True)
    ax.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax.yaxis.grid(True,linestyle='-', color='#C0C0C0')

    
def main():
    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(1,1,1)
    #ax3 = fig.add_subplot(3,1,3)
    plot_coldgas(0.0,ax1,"l")
    fig.savefig("CG00.pdf",bbox_inches='tight',pad_inches=0.1)
    plt.close(fig)

if __name__=="__main__":
    main()
