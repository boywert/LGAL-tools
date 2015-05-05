import numpy

gadget2msun=1.e10
boxsize = 47.0
hubble_h = 0.7


def uv_luminosity_fn(gal,min=-23.,max=-17.,nbins=12):
    massf = gal["Mag"][:,5] - numpy.log10(hubble_h)
    stellarmass = numpy.histogram(massf,nbins,(min,max))
    massftn_y = stellarmass[0]/(boxsize)**3/((max-min)/nbins)
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)


def metallicity_fn(gal,mass_min=1.e-5,mass_max=1.,nbins=20):
    massf = numpy.log10((gadget2msun*gal['MetalsDiskMass']+gadget2msun*gal['MetalsBulgeMass'])/(gadget2msun*gal['DiskMass']+gadget2msun*gal['BulgeMass']))
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]/(boxsize/hubble_h)**3/(numpy.log10(mass_max/mass_min)/nbins)
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)


def ssfr_density_fn(gal,mass_min=0.1,mass_max=1.,nbins=20):
    massf = numpy.log10(gal['Sfr']/(gal['DiskMass']+gal['BulgeMass']))
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]/(boxsize/hubble_h)**3/(numpy.log10(mass_max/mass_min)/nbins)
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def sfr_density_fn(gal,mass_min=0.1,mass_max=1000.,nbins=20):
    massf = numpy.log10(gal['Sfr'])
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]/(boxsize/hubble_h)**3/(numpy.log10(mass_max/mass_min)/nbins)
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def M200c_mass_fn(halos,mass_min=1e8,mass_max=1.e15,nbins=20):
    massf = numpy.log10(gadget2msun*halos['M_Crit200'])
    mass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = mass[0]/(boxsize)**3/(numpy.log10(mass_max/mass_min)/nbins)
    massftn_x = []
    for i in range(len(mass[0])):
        massftn_x.append((mass[1][i]+mass[1][i+1])/2.)
    return (massftn_x,massftn_y)


def stellar_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = numpy.log10(gadget2msun*(gal['DiskMass']+gal['BulgeMass'])/hubble_h*10.)
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]/(boxsize/hubble_h)**3/(numpy.log10(mass_max/mass_min)/nbins)
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def hotgas_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = gadget2msun*gal['HotGas']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def coldgas_mass_fn(gal,mass_min=1.,mass_max=1.e20,nbins=20):
    massf = gadget2msun*gal['ColdGas']
    stellarmass = numpy.histogram(numpy.log10(massf),nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    return (massftn_x,massftn_y)

def bh_mass_fn(gal,mass_min=1.0,mass_max=1.e8,nbins=50):
    massf = numpy.log10(gal['BlackHoleMass'])
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    for i in range(len(massftn_x)):
        massftn_x[i] = pow(10.,massftn_x[i])
    return (massftn_x,massftn_y)


def sfr_fn(gal,mass_min=0.000001,mass_max=1.,nbins=50):
    massf = numpy.log10(gal['Sfr'])
    stellarmass = numpy.histogram(massf,nbins,(numpy.log10(mass_min),numpy.log10(mass_max)))
    massftn_y = stellarmass[0]
    massftn_x = []
    for i in range(len(stellarmass[0])):
        massftn_x.append((stellarmass[1][i]+stellarmass[1][i+1])/2.)
    for i in range(len(massftn_x)):
        massftn_x[i] = pow(10.,massftn_x[i])
    return (massftn_x,massftn_y)
