import numpy
import pylab
from pylab import *
from matplotlib import gridspec
from matplotlib.colors import LogNorm
class xfrac:
    grid = 0
    data = 0

def read_xfrac(filename,doubleflag):
    f = open(filename,"rb")
    output = xfrac()
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    output.grid = numpy.fromfile(f,numpy.int32,3)
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    if(doubleflag == 1):
        output.data = numpy.fromfile(f,numpy.float64,output.grid[0]**3).reshape(( output.grid[0], output.grid[1], output. grid[2]))
    else:
        output.data = numpy.fromfile(f,numpy.float32,output.grid[0]**3).reshape(( output.grid[0], output.grid[1], output. grid[2]))
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    return output

def get_space(data3d,x,y,z):
    minx = x[0]
    maxx = x[1]
    miny = y[0]
    maxy = y[1]
    minz = z[0]
    maxz = z[1]
    a = data3d[minx:maxx,miny:maxy,minz:maxz]
    b = numpy.zeros(shape=((maxx-minx),(maxy-miny)), dtype=numpy.float64)
    for i in range(len(a)):
        for j in range(len(a[i])):
            b[i,j] = numpy.sum(a[i,j], dtype=numpy.float64)/(maxz-minz)
    return b

def get_plot(filename,doubleflag,x,y,z):
    xfrac = read_xfrac(filename,doubleflag)
    data_plot = get_space(xfrac.data,x,y,z) 
    return data_plot

#global x,y,z



def plot_reionized(nrow,ncol,filelist,labellist,doubleflaglist,frac):
    fig = pylab.figure(figsize=(4*ncol, 4*nrow+1))
    plt.subplots_adjust(wspace = 0.01)
    plt.subplots_adjust(hspace = 0.10)
    #fig.suptitle(r"$x_{\mathrm{HII}} = %3.1f$" % (frac))
    gs_width_ratios = []
    gs_height_ratios = []
    for i in range(nrow):
        gs_height_ratios.append(1.)
    for i in range(ncol):
        gs_width_ratios.append(1.)
    gs = gridspec.GridSpec(nrow, ncol, width_ratios=gs_width_ratios, height_ratios = gs_height_ratios) 
    ax = []
    im = []
    ifile = 0
    for i in range(nrow):
        ii = i
        print "Plotiing row %d" %(ii+1)
        for j in range(ncol):
            print "\tPlotiing column %d, file %d of %d" %(j+1,ii*ncol+j+1,len(filelist))
            if(ii*ncol+j < len(filelist)):
                ax.append(pylab.subplot(gs[ii,j]))
                filename = filelist[ifile] #+"/xfrac3d_"+redshift+".bin"
                data_plot = get_plot(filename,doubleflaglist[ifile],x,y,z)
                im.append(ax[ifile].imshow(data_plot+1, cmap=plt.get_cmap("Blues"), norm=LogNorm(vmin=1.0,  vmax=2.0), extent=[x[0], x[1], y[0], y[1]]))
                ax[ifile].axis("on")
		ax[ifile].set_xlabel(labellist[ifile],fontsize=10)
            	im[ifile].set_interpolation('bilinear')
                ax[ifile].yaxis.set_ticklabels([])
                ax[ifile].xaxis.set_ticklabels([])
                if j == 0:
                    ax[ifile].set_ylabel(r"47 Mpc/h",fontsize=10)
                ifile += 1
    outfile = "%3.1f_pic.pdf" % (frac)
    fig.savefig(outfile, bbox_inches='tight')
    plt.close(fig)

x = (0,306)
y = (0,306)
z = (90,120)

nrow = 3
ncol = 2

# 30%
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_8.515.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_8.283.bin"]
labellist = ["No suppression, stripping 0 (z = 8.76)",
            "No suppression, stripping 1 (z = 8.76)",
            "Homogeneous, stripping 0 (z = 8.76)",
            "Homogeneous, stripping 1 (z = 8.76)",
            "Patchy suppression, stripping 0 (z = 8.51)",
            "Patchy suppression, stripping 1 (z = 8.28)"]

doubleflaglist =[0,0,0,0,0,0]
plot_reionized(nrow,ncol,filelist,labellist,doubleflaglist,0.3)
# 70%
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_7.859.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_7.760.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_7.859.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_7.760.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_7.664.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_7.305.bin"]
labellist = ["No suppression, stripping 0 (z = 7.86)",
            "No suppression, stripping 1 (z = 7.76)",
            "Homogeneous, stripping 0 (z = 7.86)",
            "Homogeneous, stripping 1 (z = 7.76)",
            "Patchy suppression, stripping 0 (z = 7.66)",
            "Patchy suppression, stripping 1 (z = 7.30)"]
doubleflaglist =[0,0,0,0,0,0]
plot_reionized(nrow,ncol,filelist,labellist,doubleflaglist,0.7)


# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.938")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.457")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.026")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"8.515")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"7.960")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"7.480")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"6.981")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"6.483")

