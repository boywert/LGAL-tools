import matplotlib as plt
import numpy

hubble_h = 0.7

def add_obs_smf_z8(observe_folder,ax):
    data_file = observe_folder+"/song2016_z8.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = 10.**data[:,1]
    data_errorup = 10.**(data[:,2]+data[:,1])-10.**data[:,1]
    data_errordown = 10.**data[:,1] - 10.**(data[:,1]-data[:,3])
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Song et al. (2016)")
    return ax

def add_obs_uv_z8(observe_folder,ax):
    bouwens2011_file = observe_folder+"/bouwens2011_z8.txt"
    bouwens2011 = numpy.loadtxt(bouwens2011_file)
    bouwens2011_x = bouwens2011[:,0]-5.*numpy.log10(hubble_h)
    bouwens2011_y = (10.**bouwens2011[:,1])/hubble_h**3.

    bouwens2011_errorup = (10.**(bouwens2011[:,1] + bouwens2011[:,5]) - 10.**bouwens2011[:,1])/hubble_h**3.
    bouwens2011_errordown = (10.**bouwens2011[:,1] - 10.**(bouwens2011[:,1] + bouwens2011[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2011_x,bouwens2011_y,yerr=[ bouwens2011_errordown, bouwens2011_errorup], fmt='o',label="Bouwens et al. (2011)") 

    data_file = observe_folder+"/bouwens2014_z8.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = data[:,1]/hubble_h**3.
    data_error = data[:,2]/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Bouwens et al. (2015)")
    
    mclure2010_file = observe_folder+"/mclure2010_z8.txt"
    mclure2010= numpy.loadtxt(mclure2010_file)
    mclure2010_x = mclure2010[:,0]-5.*numpy.log10(hubble_h)
    mclure2010_y = (10.**mclure2010[:,1])/hubble_h**3.

    mclure2010_errorup = (10.**(mclure2010[:,1] + mclure2010[:,3]) - 10.**mclure2010[:,1])/hubble_h**3.
    mclure2010_errordown = (10.**mclure2010[:,1] - 10.**(mclure2010[:,1] + mclure2010[:,2]))/hubble_h**3.
    ax.errorbar(mclure2010_x,mclure2010_y,yerr=[mclure2010_errordown,mclure2010_errorup], fmt='o',label="Mclure et al. (2011)") 
  
    bouwens2010_file = observe_folder+"/bouwens2010_z8.txt"
    bouwens2010 = numpy.loadtxt(bouwens2010_file)
    bouwens2010_x = bouwens2010[:,0]-5.*numpy.log10(hubble_h)
    bouwens2010_y = (10.**bouwens2010[:,1])/hubble_h**3.

    bouwens2010_errorup = (10.**(bouwens2010[:,1] + bouwens2010[:,3]) - 10.**bouwens2010[:,1])/hubble_h**3.
    bouwens2010_errordown = (10.**bouwens2010[:,1] - 10.**(bouwens2010[:,1] + bouwens2010[:,2]))/hubble_h**3.
    ax.errorbar(bouwens2010_x,bouwens2010_y,yerr=[bouwens2010_errordown,bouwens2010_errorup], fmt='o',label="Bouwens et al. (2010)")
    return ax

##################################################################################################################
# z = 7

def add_obs_uv_z7(observe_folder,ax):
    hubble_h=0.7

    bouwens2011_file = observe_folder+"/bouwens2011_z7.txt"
    bouwens2011 = numpy.loadtxt(bouwens2011_file)
    bouwens2011_x = bouwens2011[:,0]-5.*numpy.log10(hubble_h)
    bouwens2011_y = (10.**bouwens2011[:,1])/hubble_h**3.

    bouwens2011_errorup = (10.**(bouwens2011[:,1] + bouwens2011[:,5]) - 10.**bouwens2011[:,1])/hubble_h**3.
    bouwens2011_errordown = (10.**bouwens2011[:,1] - 10.**(bouwens2011[:,1] + bouwens2011[:,4]))/hubble_h**3.
    ax.errorbar(bouwens2011_x,bouwens2011_y,yerr=[ bouwens2011_errordown, bouwens2011_errorup], fmt='o',label="Bouwens et al. (2011)") 
    
    data_file = observe_folder+"/bouwens2014_z7.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = data[:,1]/hubble_h**3.
    data_error = data[:,2]/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Bouwens et al. (2015)")
    
    mclure2010_file = observe_folder+"/mclure2010_z8.txt"
    mclure2010= numpy.loadtxt(mclure2010_file)
    mclure2010_x = mclure2010[:,0]-5.*numpy.log10(hubble_h)
    mclure2010_y = (10.**mclure2010[:,1])/hubble_h**3.

    mclure2010_errorup = (10.**(mclure2010[:,1] + mclure2010[:,3]) - 10.**mclure2010[:,1])/hubble_h**3.
    mclure2010_errordown = (10.**mclure2010[:,1] - 10.**(mclure2010[:,1] + mclure2010[:,2]))/hubble_h**3.
    ax.errorbar(mclure2010_x,mclure2010_y,yerr=[mclure2010_errordown,mclure2010_errorup], fmt='o',label="Mclure et al. (2011)") 
  
    oesch2010_file = observe_folder+"/oesch2010_z7.txt"
    oesch2010 = numpy.loadtxt(oesch2010_file)
    oesch2010_x = oesch2010[:,0]-5.*numpy.log10(hubble_h)
    oesch2010_y = (10.**oesch2010[:,1])/hubble_h**3.

    oesch2010_errorup = (10.**(oesch2010[:,1] + oesch2010[:,3]) - 10.**oesch2010[:,1])/hubble_h**3.
    oesch2010_errordown = (10.**oesch2010[:,1] - 10.**(oesch2010[:,1] + oesch2010[:,2]))/hubble_h**3.
    ax.errorbar(oesch2010_x,oesch2010_y,yerr=[ oesch2010_errordown, oesch2010_errorup], fmt='o',label="Oesch et al. (2010)")
    return ax


def add_obs_sfr_z7(observe_folder,ax):
    data_file = observe_folder+"/smit2012_z7.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]+numpy.log10(0.55)
    data_y = data[:,1]
    data_error = data[:,2]
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Smit et al. (2012)")
    
    data_file = observe_folder+"/Duncan14_SFRF_SED_z7.cat"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = data[:,1]

    data_errorup = data[:,3]
    data_errordown = data[:,2]
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Duncan et al. (2014)")
    
def add_obs_smf_z7(observe_folder,ax):
    data_file = observe_folder+"/gonzalez2011_z7.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = 10.**data[:,1]
    data_errorup = 10.**(data[:,2]+data[:,1])-10.**data[:,1]
    data_errordown = 10.**data[:,1] - 10.**(data[:,1]+data[:,3])
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Gonzalez et al. (2011)")
    
    data_file = observe_folder+"/Duncan14_MF_z7.cat"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = data[:,1]
    data_errorup = data[:,3]
    data_errordown = data[:,2]
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Duncan et al. (2014)")

    data_file = observe_folder+"/song2016_z7.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = 10.**data[:,1]
    data_errorup = 10.**(data[:,2]+data[:,1])-10.**data[:,1]
    data_errordown = 10.**data[:,1] - 10.**(data[:,1]-data[:,3])

    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Song et al. (2016)")
    

#####################################################################################################################
# z = 6  

def add_obs_uv_z6(observe_folder,ax):
    data_file = observe_folder+"/bouwens2007_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = (10.**data[:,1])/hubble_h**3.
    data_errorup = (10.**(data[:,1] + data[:,3]) - 10.**data[:,1])/hubble_h**3.
    data_errordown = (10.**data[:,1] - 10.**(data[:,1] + data[:,2]))/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Bouwens et al. (2007)") 

    data_file = observe_folder+"/bouwens2014_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = data[:,1]/hubble_h**3.
    data_error = data[:,2]/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Bouwens et al. (2015)")

    data_file = observe_folder+"/Duncan14_LF_z6.txt"
    data = numpy.loadtxt(data_file)
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = (data[:,1])/hubble_h**3.
    data_errorup = data[:,3]/hubble_h**3.
    data_errordown = data[:,2]/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Duncan et al. (2014)")

    data_file = observe_folder+"/bowler2014_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]-5.*numpy.log10(hubble_h)
    data_y = data[:,1]/hubble_h**3.
    data_error = data[:,2]/hubble_h**3.
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Bowler et al. (2014)")

    
    return ax

def add_obs_sfr_z6(observe_folder,ax):
    data_file = observe_folder+"/smit2012_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]+numpy.log10(0.55)
    data_y = data[:,1]
    data_error = data[:,2]
    ax.errorbar(data_x,data_y,yerr=data_error, fmt='o',label="Smit et al. (2012)")
    
    data_file = observe_folder+"/Duncan14_SFRF_SED_z6.cat"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = data[:,1]

    data_errorup = data[:,3]
    data_errordown = data[:,2]
    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Duncan et al. (2014)")

    
def add_obs_smf_z6(observe_folder,ax):
    # data_file = observe_folder+"/gonzalez2011_z6.txt"
    # data = numpy.loadtxt(data_file)
    # data_x = data[:,0]
    # data_y = 10.**data[:,1]
    # data_errorup = 10.**(data[:,2]+data[:,1])-10.**data[:,1]
    # data_errordown = 10.**data[:,1] - 10.**(data[:,1]+data[:,3])
    # ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Gonzalez et al. (2011)")
    
    # data_file = observe_folder+"/Duncan14_MF_z6.cat"
    # data = numpy.loadtxt(data_file)
    # data_x = data[:,0]
    # data_y = data[:,1]
    # data_errorup = data[:,3]
    # data_errordown = data[:,2]
    # ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Duncan et al. (2014)")
    
    data_file = observe_folder+"/song2016_z6.txt"
    data = numpy.loadtxt(data_file)
    data_x = data[:,0]
    data_y = 10.**data[:,1]
    data_errorup = 10.**(data[:,2]+data[:,1])-10.**data[:,1]
    data_errordown = 10.**data[:,1] - 10.**(data[:,1]-data[:,3])

    ax.errorbar(data_x,data_y,yerr=[data_errordown,data_errorup], fmt='o',label="Song et al. (2016)")
    


