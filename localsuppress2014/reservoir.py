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
    filter['Mvir'] = True
    filter['Rvir'] = True
    filter['HaloM_Crit200'] = True
    filter['HotGas'] = True
    filter['ColdGas'] = True
    filter['EjectedMass'] = True
    filter['StellarMass'] = True
    filter['ICM'] = True
    filter['Pos'] = True
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
    round_z = "%4.2f" % (float(int(float(z)+0.5)))
    #fmod_file = "fmods/"+round_z+".txt"
    #fmod_model = numpy.loadtxt(fmod_file)
    #fmod_model[:,0] = numpy.log10(fmod_model[:,0]/hubble_h)
    dt,filter = setfilter(models)
    file_prefix = "SA_z"+z
    try:
        gal
    except NameError:
        gal = {}
        nTrees = {}
        nGals = {}
        nTreeGals = {}
    sum_hotgas = {}
    sum_coldgas = {}
    sum_ejectedmass = {}
    N = {}
    m200c = {}
    for i in range(len(models.model_names)):
        index = models.model_names[i]
        print index
        if not index in gal:
            (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(models.model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)
        rangen = (7.5,11.5)
        bins = 40
        total_baryon = numpy.float64(1)*(gal[index]["StellarMass"]+gal[index]["EjectedMass"]+gal[index]["ColdGas"]+gal[index]['HotGas']+gal[index]["ICM"]+gal[index]["BlackHoleMass"]+gal[index]["BlackHoleGas"])
        firstgal = numpy.where(gal[index]["Type"] == 0)[0]
        cenmass = numpy.zeros(len(firstgal),dtype=numpy.float64)
        hotgas = numpy.zeros(len(firstgal),dtype=numpy.float64)
        coldgas = numpy.zeros(len(firstgal),dtype=numpy.float64)
        ejectedmass = numpy.zeros(len(firstgal),dtype=numpy.float64)
        for ii in range(len(firstgal)-1):
            for j in range(len(total_baryon[firstgal[ii]:firstgal[ii+1]])):
                #print total_baryon[firstgal[i]:firstgal[i+1]]
                this_gal = firstgal[ii]+j
                distance = numpy.sqrt((gal[index][this_gal]['Pos'][0] - gal[index][firstgal[ii]]['Pos'][0])**2.+(gal[index][this_gal]['Pos'][1] - gal[index][firstgal[ii]]['Pos'][1])**2.+(gal[index][this_gal]['Pos'][2] - gal[index][firstgal[ii]]['Pos'][2])**2.)/(1.+float(z))
                if ( distance < gal[index][firstgal[ii]]['Rvir']):
                    cenmass[ii] += total_baryon[this_gal]
                    hotgas[ii] += numpy.float64(1)*gal[index][this_gal]['HotGas']
                    coldgas[ii] += numpy.float64(1)*gal[index][this_gal]['ColdGas']
                    ejectedmass[ii] += numpy.float64(1)*gal[index][this_gal]['EjectedMass']
        cenhalomass = gal[index]["Mvir"][firstgal]
        cond = numpy.where(~numpy.isnan(cenmass) & ~numpy.isnan(cenhalomass) & (cenhalomass > 0.) & (gal[index]["Type"][firstgal] == 0))[0]
        cenmass = cenmass[cond]
        cenhalomass = cenhalomass[cond]
        hotgas = hotgas[cond]
        coldgas = coldgas[cond]
        ejectedmass = ejectedmass[cond]
        print index,len(cenmass),len(cenhalomass),len(firstgal)
	#print "min",numpy.amin(cenmass/cenhalomass/0.165), "max", numpy.amax(cenmass/cenhalomass/0.165)
        sum_hotgas[index] = numpy.histogram(numpy.log10(cenhalomass*1.e10/hubble_h),range=rangen,bins=bins,weights=(hotgas/cenhalomass/0.165))
        sum_coldgas[index] = numpy.histogram(numpy.log10(cenhalomass*1.e10/hubble_h),range=rangen,bins=bins,weights=(coldgas/cenhalomass/0.165))
        sum_ejectedmass[index] = numpy.histogram(numpy.log10(cenhalomass*1.e10/hubble_h),range=rangen,bins=bins,weights=(ejectedmass/cenhalomass/0.165))
        N[index] = numpy.histogram(numpy.log10(cenhalomass*1.e10/hubble_h),range=rangen,bins=bins)
        m200c[index] = []
        for ii in range(len(N[index][0])):
            m200c[index].append(0.5*(N[index][1][ii]+N[index][1][ii+1]))
        del(gal[index])
        del(nTreeGals[index])
        m200c[index] = numpy.array(m200c[index])

    for i in range(len(models.model_names)):
        index = models.model_names[i]
        mean_hotgas = sum_hotgas[index][0]/N[index][0]
        mean_coldgas = sum_coldgas[index][0]/N[index][0]
        mean_ejectedmass = sum_ejectedmass[index][0]/N[index][0]
        cond = ~numpy.isnan(1./N[index][0])
        mean_hotgas = mean_hotgas[cond]
        mean_coldgas = mean_coldgas[cond]
        mean_ejectedmass = mean_ejectedmass[cond]
        m200c[index] = m200c[index][cond]
        ax.plot(m200c[index],mean_hotgas,color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i],label="Hot "+models.model_labels[i])
        ax.plot(m200c[index],mean_coldgas,color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i],label="Cold "+models.model_labels[i])
        ax.plot(m200c[index],mean_ejectedmass,color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i],"Ejected "+label=models.model_labels[i])
    if pos == "r":
        # leg = ax.legend(loc=4, handlelength = 10,ncol=1, fancybox=True, prop={'size':12})
        # leg.get_frame().set_linewidth(0)
        # leg.get_frame().set_alpha(0)
        ax.yaxis.set_ticklabels([])
        labels = ["",r"$8.5$",r"$9.0$",r"$9.5$",r"$10.0$",r"$10.5$",r"$11.0$",r"$11.0$",r"$11.5$"]
        #ax.xaxis.set_ticklabels(labels)
    ax.set_ylim([0,1.0])
    #ax.set_yscale('log')
    ax.set_xlim([8,11])
    ax.set_xlabel(r"$\log_{10}(M_{200c}/\mathrm{M_\odot})$")
    if pos == "l":
        ax.set_ylabel(r"$f_{\mathrm{baryon}}/f_b$")
    if pos == "l":
        ax.text(0.5, 0.25, 'stripping 0',
                verticalalignment='bottom', horizontalalignment='center',
                transform=ax.transAxes, fontsize=14)
        ax.text(0.5, 0.9, 'z = '+str(int(float(z)+0.5)),
                verticalalignment='bottom', horizontalalignment='center',
                transform=ax.transAxes, fontsize=14)
    else:
        ax.text(0.5, 0.25, 'stripping 1',
                verticalalignment='bottom', horizontalalignment='center',
                transform=ax.transAxes, fontsize=14)
    ax.set_axisbelow(True)
    ax.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax.yaxis.grid(True,linestyle='-', color='#C0C0C0')

        
def main():
    import globalconf as model1
    zi = "6.00"
    fig = plt.figure(figsize=(16, 6))
    plt.subplots_adjust(wspace = 0)
    ax1 = fig.add_subplot(121)
    plot_z(zi,model1,ax1,"l")
    ax2 = fig.add_subplot(122)
    fig.canvas.draw()
    zi = "9.03"
    plot_z(zi,model1,ax2,"r")
    fig.savefig("reservoir69.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
