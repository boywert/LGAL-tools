import matplotlib as plt
import numpy

hubble_h = 0.7
def add_obs_uv_z8(observe_folder,ax):
    
    bouwens2011_file = observe_folder+"/bouwens2011_z8.txt"
    bouwens2011 = numpy.loadtxt(bouwens2011_file)
    bouwens2011_x = bouwens2011[:,0]-5.*numpy.log10(hubble_h)
    bouwens2011_y = (10.**bouwens2011[:,1])/hubble_h**3.

    bouwens2011_errorup = (10.**(bouwens2011[:,1] + bouwens2011[:,5]) - 10.**bouwens2011[:,1])/hubble_h**3.
    bouwens2011_errordown = (10.**bouwens2011[:,1] - 10.**(bouwens2011[:,1] + bouwens2011[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2011_x,bouwens2011_y,yerr=bouwens2011_y/10., fmt='o',label="Bouwens et al. (2011)") 

    mclure2010_file = observe_folder+"/mclure2010_z8.txt"
    mclure2010= numpy.loadtxt(mclure2010_file)
    mclure2010_x = mclure2010[:,0]-5.*numpy.log10(hubble_h)
    mclure2010_y = (10.**mclure2010[:,1])/hubble_h**3.

    mclure2010_errorup = (10.**(mclure2010[:,1] + mclure2010[:,5]) - 10.**mclure2010[:,1])/hubble_h**3.
    mclure2010_errordown = (10.**mclure2010[:,1] - 10.**(mclure2010[:,1] + mclure2010[:,4]))/hubble_h**3.
    ax.errorbar(mclure2010_x,mclure2010_y,yerr=mclure2010_y/10., fmt='o',label="Mclure et al. (2011)") 
  
    bouwens2010_file = observe_folder+"/bouwens2010_z8.txt"
    bouwens2010 = numpy.loadtxt(bouwens2011_file)
    bouwens2010_x = bouwens2010[:,0]-5.*numpy.log10(hubble_h)
    bouwens2010_y = (10.**bouwens2010[:,1])/hubble_h**3.

    bouwens2010_errorup = (10.**(bouwens2010[:,1] + bouwens2010[:,5]) - 10.**bouwens2010[:,1])/hubble_h**3.
    bouwens2010_errordown = (10.**bouwens2010[:,1] - 10.**(bouwens2010[:,1] + bouwens2010[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2010_x,bouwens2010_y,yerr=bouwens2010_y/10., fmt='o',label="Bouwens et al. (2010)")

    

    return ax

def add_obs_uv_z6(observe_folder,ax):
    data_file = observe_folder+"bouwens2007_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = (10.**data[:,1])/hubble_h**3.

    data_errorup = (10.**(data[:,1] + data[:,3]) - 10.**data[:,1])/hubble_h**3.
    data_errordown = (10.**data[:,1] - 10.**(data[:,1] + data[:,2]))/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_y/10., fmt='o',label="Bouwens et al. (2007)") 
    return ax


