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
try:
   import cPickle as pickle
except:
   import pickle
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
    filter['BulgeMass'] = True
    filter['DiskMass'] = True
    filter['HaloM_Crit200'] = True
    filter['Rvir'] = True
    filter['Pos'] = True
    filter['Mvir'] = True
    filter['CumulativeSFR'] = True
    # filter['StellarMass'] = True
    # filter['ICM'] = True
    # filter['BlackHoleGas'] = True
    # filter['BlackHoleMass'] = True
    # filter['Sfr'] = True
    filter['Type'] = True
    filter['StellarMass'] = True
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
def plot_z(z,models,ax,pos,label=0,bottom=0,top=0):
    round_z = "%4.2f" % (float(int(float(z)+0.5)))
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
        cachefile = index+"_"+round_z+"_gal.pickle"
        rangen = (6.0,13)
        bins = 50
        if os.path.isfile(cachefile) == False: 
            if not index in gal:
                (nTrees[index],nGals[index],nTreeGals[index],gal[index]) = read_lgal.readsnap_lgal_advance(models.model_paths[i],file_prefix,firstfile,lastfile,filter[i],dt[i],1)

            firstgal = numpy.where(gal[index]["Type"] == 0)[0]
            star = numpy.zeros(len(firstgal),dtype=numpy.float64)
            stargross =  numpy.zeros(len(firstgal),dtype=numpy.float64)
            for ii in range(len(firstgal)-1):
                for j in range(firstgal[ii+1]-firstgal[ii]):
                    #print total_baryon[firstgal[i]:firstgal[i+1]]
                    this_gal = firstgal[ii]+j
                    distance = numpy.sqrt((gal[index][this_gal]['Pos'][0] - gal[index][firstgal[ii]]['Pos'][0])**2.+(gal[index][this_gal]['Pos'][1] - gal[index][firstgal[ii]]['Pos'][1])**2.+(gal[index][this_gal]['Pos'][2] - gal[index][firstgal[ii]]['Pos'][2])**2.)/(1.+float(z))
                    if ( distance < gal[index][firstgal[ii]]['Rvir']):
                        star[ii] += gal[index][this_gal]['StellarMass']
                        stargross[ii] += gal[index][this_gal]['CumulativeSFR']
            gal[index][firstgal]['StellarMass'] = star
            gal[index][firstgal]['CumulativeSFR'] = stargross
            fp = open(cachefile,'wb')
            pickle.dump(gal[index],fp)
            fp.close()
        else:
            fp = open(cachefile,'rb')
            gal[index] = pickle.load(fp)
            fp.close()
       	gal[index] = gal[index][numpy.where((gal[index]["Type"]==0)&((gal[index]["StellarMass"]) >0.))]
        #gal[index] = gal[index][gal[index]["Type"]==0]
	mass = gal[index]['HaloM_Crit200']# (gal[index]["BulgeMass"]+gal[index]["DiskMass"])
        sum_baryons[index] = numpy.histogram(numpy.log10(mass*1.e10/hubble_h),range=rangen,bins=bins,weights=numpy.log10(numpy.float64(1)*(gal[index]["StellarMass"])*1.e10/hubble_h))
        sum_baryons_sq[index] = numpy.histogram(numpy.log10(mass*1.e10/hubble_h),range=rangen,bins=bins,weights=numpy.log10(numpy.float64(1)*(gal[index]["StellarMass"])*1.e10/hubble_h)**2)
        N[index] = numpy.histogram(numpy.log10(mass*1.e10/hubble_h),range=rangen,bins=bins)
        m200c[index] = []
        for i in range(len(sum_baryons[index][0])):
            m200c[index].append(0.5*(sum_baryons[index][1][i]+sum_baryons[index][1][i+1]))
        del(gal[index])
        del(nTreeGals[index])
        m200c[index] = numpy.array(m200c[index])
  
    for i in range(len(models.model_names)):
        index = models.model_names[i]
        mean = sum_baryons[index][0]/N[index][0]
        sd =  numpy.sqrt(numpy.fabs(sum_baryons_sq[index][0]/N[index][0] - mean**2))
        cond = ~numpy.isnan(mean)
        mean = mean[cond]
        sd = sd[cond]
        m200c[index] = m200c[index][cond]
    
        for j in range(len(m200c[index])):
            print m200c[index][j],mean[j]
        if pos == "l":
            ax.plot(m200c[index],mean,color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i])
        elif pos == "r":
            print models.model_labels[i]
            ax.plot(m200c[index],mean,color=models.model_plot_colors[i],linestyle=models.model_plot_patterns[i],label=models.model_labels[i])
        ax.fill_between(m200c[index], mean - sd, mean + sd, alpha=0.25, edgecolor='#CC4F1B', facecolor=models.model_plot_colors[i],linewidth=0)
    xplot = numpy.arange(0,20.)
    ref = (-9.25+float(z)/24)+xplot*1.7-0.8
    ax.plot(xplot,ref,'k--', label = r'$m_{\mathrm{*}} \propto M_{\mathrm{200c}}^{1.7}$')
  
    ax.set_ylim([3,10])
    ax.set_xlim([8.0,11.])
    ax.set_xlabel(r"$\log_{10}(M_{\mathrm{200c}}/\mathrm{M_\odot})$")
    
    if pos == "r":
        labels = ["",r"$8.5$",r"$9.0$",r"$9.5$",r"$10.0$",r"$10.5$",r"$11.0$",r"$11.5$",r"$12.0$"]
        ax.xaxis.set_ticklabels(labels)
    if bottom != 1:
        ax.xaxis.set_ticklabels([])
    if label == 1:
        leg = ax.legend(loc=4, handlelength = 10,ncol=1, fancybox=True, prop={'size':12})
        leg.get_frame().set_linewidth(0)
        leg.get_frame().set_alpha(0)
    if pos == "l":
        ax.set_ylabel(r"$\log_{10}(m_{\mathrm{*}}/\mathrm{M_\odot})$")
    if top == 1:
        ax2 = ax.twiny()
        ax2.set_xlim([8.0,11.])
        ax2.set_axisbelow(True)
        ax2.xaxis.grid(True,linestyle='-', color='#C0C0C0')
        ax2.yaxis.grid(True,linestyle='-', color='#C0C0C0')
        if pos =="r":
            labels = ["",r"$8.5$",r"$9.0$",r"$9.5$",r"$10.0$",r"$10.5$",r"$11.0$",r"$11.5$",r"$12.0$"]
            ax2.xaxis.set_ticklabels(labels)
    if top == 0:
        labels = [r"$3$",r"$4$",r"$5$",r"$6$",r"$7$",r"$8$",r"$9$",""]
        ax.yaxis.set_ticklabels(labels)
    if pos == "r":
        ax.yaxis.set_ticklabels([])
#     if pos == "l":
#         ax.text(0.9, 0.95, 'stripping 0',
#                 verticalalignment='bottom', horizontalalignment='right',
#                 transform=ax.transAxes, fontsize=14)
    ax.text(0.5, 0.9, 'z = '+str(int(float(z)+0.5)),
            verticalalignment='bottom', horizontalalignment='center',
            transform=ax.transAxes, fontsize=14)
    # else:
#         ax.text(0.9, 0.95, 'stripping 1',
#                 verticalalignment='bottom', horizontalalignment='right',
#                 transform=ax.transAxes, fontsize=14)

    ax.set_axisbelow(True)
    ax.xaxis.grid(True,linestyle='-', color='#C0C0C0')
    ax.yaxis.grid(True,linestyle='-', color='#C0C0C0')


    
def main():
    zlist = open(zlistfile).readlines()
    #zi = zlist[long(sys.argv[1])].strip()
    fig = plt.figure(figsize=(16, 18))
    plt.subplots_adjust(wspace = 0,hspace = 0)
    import model2 as model1
    ax1 = fig.add_subplot(321)
    zi = zlist[75].strip()
    plot_z(zi,model1,ax1,"l",top=1)
    ax2 = fig.add_subplot(322)
    fig.canvas.draw()
    zi = zlist[49].strip()
    plot_z(zi,model1,ax2,"r",label=1,top=1)
    ax2 = fig.add_subplot(323)
    fig.canvas.draw()
    zi = zlist[34].strip()
    plot_z(zi,model1,ax2,"l")
    ax2 = fig.add_subplot(324)
    fig.canvas.draw()
    zi = zlist[24].strip()
    plot_z(zi,model1,ax2,"r")
    ax2 = fig.add_subplot(325)
    fig.canvas.draw()
    zi = zlist[18].strip()
    plot_z(zi,model1,ax2,"l",bottom=1)
    ax2 = fig.add_subplot(326)
    fig.canvas.draw()
    zi = zlist[13].strip()
    plot_z(zi,model1,ax2,"r",bottom=1)
    fig.savefig("mstar_m200c2.pdf",bbox_inches='tight',pad_inches=0.05)
    plt.close(fig)

if __name__=="__main__":
    main()
