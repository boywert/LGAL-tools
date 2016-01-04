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
ranki = str(random.randint(0,1000000))
os.system("mkdir -p ../tmp/"+ranki)
pylab.rc('text', usetex=True)
pylab.rc('lines', linewidth=2)
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.major.size'] = 8

def loadfilter(structfile):
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
    filter['CumulativeSFR'] = True
    dt = LGalaxyStruct.struct_dtype
    return (filter,dt)

def setfilter(models):
    dt = []
    filter = []
    for i in range(len(models.struct_file)):
        (f,t) = loadfilter(models.struct_file[i])
        filter.append(f)
        dt.append(t)
    return dt,filter


# filter_tmp = []
# dt_tmp = []
# model_names_tmp = []
# struct_file_tmp = []
# model_labels_tmp = []
# model_paths_tmp = []
# for i in range(len(use_model)):
#     if use_model[i]:
#         filter_tmp.append(filter[i])
#         dt_tmp.append(dt[i])
#         model_names_tmp.append(model_names[i])
#         struct_file_tmp.append(struct_file[i])
#         model_labels_tmp.append(model_labels[i])
#         model_paths_tmp.append(model_paths[i])

# filter = filter_tmp
# dt = dt_tmp
# model_names = model_names_tmp
# struct_file = struct_file_tmp
# model_labels = models.model_labels_tmp
# models.model_paths = models.model_paths_tmp
def plot_z(z,models,ax,pos):    
    dt,filter = setfilter(models)
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sum_baryons_sq = {}
    sum_baryons = {}
    N = {}
    m200c = {}
    for i in range(len(models.model_names)):
        index = models.model_names[i]
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(models.model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        rangen = (7.5,11.5)
        bins = 40
        gal[index] = gal[index][numpy.where((gal[index]["HaloM_Crit200"] >0.))]
        
        sum_baryons[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=((numpy.float64(1)*gal[index]["StellarMass"]+gal[index]["EjectedMass"]+gal[index]["ColdGas"]+gal[index]['HotGas']+gal[index]["ICM"]+gal[index]["BlackHoleMass"]+gal[index]["BlackHoleGas"])/gal[index]["HaloM_Crit200"]/0.165))
        sum_baryons_sq[index] = numpy.histogram(numpy.log10(gal[index]["HaloM_Crit200"]*1.e10),range=rangen,bins=bins,weights=((numpy.float64(1)*gal[index]["StellarMass"]+gal[index]["EjectedMass"]+gal[index]["ColdGas"]+gal[index]['HotGas']+gal[index]["ICM"]+gal[index]["BlackHoleMass"]+gal[index]["BlackHoleGas"])/gal[index]["HaloM_Crit200"]/0.165)**2)
        N[index] = numpy.histogram(numpy.log10(gal[index]["Mvir"]*1.e10),range=rangen,bins=bins)
        m200c[index] = []
        for i in range(len(sum_baryons[index][0])):
            m200c[index].append(0.5*(sum_baryons[index][1][i]+sum_baryons[index][1][i+1]))
        del(gal[index])
        del(nTreeGals[index])

    for i in range(len(models.model_names)):
        index = models.model_names[i]
        ax.plot(m200c[index],sum_baryons[index][0]/N[index][0],color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i],label=models.model_labels[i])
        mean = sum_baryons[index][0]/N[index][0]
        sd =  numpy.sqrt(numpy.fabs(sum_baryons_sq[index][0]/N[index][0] - mean**2))
        #print mean,sd
        ax.fill_between(m200c[index], sum_baryons[index][0]/N[index][0] - sd, sum_baryons[index][0]/N[index][0] + sd, alpha=0.25, edgecolor='#CC4F1B', facecolor=models.model_plot_colors[i],linewidth=0)
    if pos == "r":
        leg = ax.legend(loc=4, handlelength = 10,ncol=1, fancybox=True, prop={'size':14})
        leg.get_frame().set_linewidth(0)
        ax.yaxis.set_ticklabels([])
        labels = [item.get_text() for item in ax.xaxis.get_ticklabels()]
        print labels
        labels[0] = ""
        print labels
        ax.xaxis.set_ticklabels(labels)
    ax.set_ylim([0,1.4])
    ax.set_xlabel(r"$M_{200c}[h^{-1}\mathrm{M_\odot}]$")
    if pos == "l":
        ax.set_ylabel(r"$f_{\mathrm{Baryon}}/f_b$")
    #ax.set_yscale("log")


    
def main():
    zlist = open(zlistfile).readlines()
    zi = zlist[long(sys.argv[1])].strip()
    fig = plt.figure(figsize=(16, 6))
    plt.subplots_adjust(wspace = 0)
    import model1 as model1
    ax1 = fig.add_subplot(121)
    plot_z(zi,model1,ax1,"l")
    import model2 as model2
    ax2 = fig.add_subplot(122)
    fig.canvas.draw()
    plot_z(zi,model2,ax2,"r")
    fig.savefig("Baryons.pdf",bbox_inches='tight',pad_inches=0)
    plt.close(fig)

if __name__=="__main__":
    main()
