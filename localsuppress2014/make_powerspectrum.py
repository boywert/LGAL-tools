import numpy
import pylab
from pylab import *
from matplotlib import gridspec
import power_spectrum

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
        output.data = numpy.fromfile(f,numpy.float64,output.grid[0]**3).reshape(( output.grid[0], output.grid[1], output.grid[2]))
    else:
        output.data = numpy.fromfile(f,numpy.float32,output.grid[0]**3).reshape(( output.grid[0], output.grid[1], output.grid[2]))
    padd = numpy.fromfile(f,numpy.int32,1)[0]
    return output



def do_powerspectrum(filename,doubleflag):
    xfrac = read_xfrac(filename,doubleflag).data
    ps,bins = power_spectrum.power_spectrum_1d(xfrac, [47.,47.,47.], kbins=153)
    return ps

def plot_powerspectrum(filelist,doubleflaglist,redshift):    
    fig = pylab.figure()
    ax = fig.add_subplot(111)
    for ifile in range(len(filelist)):
        filename = filelist[ifile]+"/xfrac3d_"+redshift+".bin"
        ps = do_powerspectrum(filename,doubleflaglist[ifile])
        ax.plot(ps)
    ax.set_yscale('log')
    fig.savefig(redshift+"_ps.pdf", bbox_inches='tight')
    

filelist = ["/mnt/lustre/scratch/cs390/codes/ionz_codes/nosupwithnohist/43000.00/","/mnt/lustre/scratch/cs390/codes/ionz_codes/okamotowithnohist/43000.00/","/mnt/lustre/scratch/cs390/47Mpc/couple/n306/xfrac/43000.00/"]
doubleflaglist =[0,0,0]
plot_powerspectrum(filelist,doubleflaglist,"9.938")
# plot_powerspectrum(filelist,doubleflaglist,"9.026")
# plot_powerspectrum(filelist,doubleflaglist,"7.960")
# plot_powerspectrum(filelist,doubleflaglist,"6.981")

