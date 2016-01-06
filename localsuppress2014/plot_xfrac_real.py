import numpy
import pylab
from pylab import *
from matplotlib import gridspec

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



def plot_reionized(nrow,ncol,filelist,doubleflaglist):
    fig = pylab.figure(figsize=(8*nrow, 8*ncol))
    plt.subplots_adjust(wspace = 0)
    plt.subplots_adjust(hspace = 0)
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
        for j in range(ncol):
            if(i*nrow+j < len(filelist)):
                ax.append(pylab.subplot(gs[i,j]))
                filename = filelist[ifile] #+"/xfrac3d_"+redshift+".bin"
                data_plot = get_plot(filename,doubleflaglist[ifile],x,y,z)
                im.append(ax[ifile].imshow(data_plot, cmap=cm.RdBu, vmin=0.0, vmax=1.0, extent=[x[0], x[1], y[0], y[1]]))
                ax[ifile].axis("on")
		#ax[ifile].set_xlabel(filelist[ifile])
                #ifile += 1
            	im[ifile].set_interpolation('bilinear')
		ifile += 1
    fig.savefig("0.3_pic.pdf", bbox_inches='tight')
    

x = (0,306)
y = (0,306)
z = (200,230)

nrow = 3
ncol = 2

# 30%
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_8.762.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_8.515.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_8.283.bin"]
doubleflaglist =[0,0,0,0]
plot_reionized(nrow,ncol,filelist,doubleflaglist)
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.938")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.457")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"9.026")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"8.515")
# plot_reionized(nrow,ncol,filelist,doubleflaglist,"7.960")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"7.480")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"6.981")
# #plot_reionized(nrow,ncol,filelist,doubleflaglist,"6.483")

