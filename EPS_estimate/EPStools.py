import matplotlib as plt
plt.use("Agg")
from numpy import *
from pylab import *
import time
import sys
sys.path.append("../python/")
from  read_lgal_advance import *
omegam = 0.3
Mpc2m = 3.08567758e22
m2Mpc = 1./Mpc2m
Msun2kg = 1.98855e30
kg2Msun = 1./Msun2kg
m2km = 0.001
G = 6.674e-11   # SI
H0 = 100.0        # km/s / (Mpc/h)
Msun2Gadget = 1e-10
Gadget2Msun = 1./Msun2Gadget
G = G * (m2km**2.) * (m2Mpc) / (kg2Msun * Msun2Gadget) # (Mpc/h) (km/s)^2 / (1e10Msun/h)
rho_crit_0 = 3.* H0**2 / (8.*pi*G)  # (1e10 Msun/h)/(Mpc/h)^3
rho_m_0 = omegam*rho_crit_0
overdensity = 200.0
psdoc = "PS/billennium_powspec.txt"
ps = loadtxt(psdoc)
ps = 10.**ps
lastsnap = 79

def find_closest(m_list,m):
    l = 0
    u = len(m_list)-1
    c = (u+l)/2
    while m_list[l] < m_list[c]:
        if m > m_list[c]:
            l = c
        elif m < m_list[c]:
            u = c
        else:
            return c
        c = (u+l)/2
    return c
def cal_radius(mass):
    mass = Msun2Gadget * 10.**mass
    return ((mass*3/(4*pi*overdensity*rho_crit_0))**(1./3))
def cal_variance(k_min,k_max):
    var = 0.
    for i in range(len(ps[:,0])-1):
        k = ps[i,0]
        k_s = ps[i+1,0]
        if (k >= k_min) & (k <= k_max):
            if (k_min >= k) & (k_min <= k_s):
                y_min = ps[i,1] + (ps[i+1,1] - ps[i,1])/(ps[i+1,0] - ps[i,0])*(k_min - ps[i,0])
                var += 0.5*(ps[i+1,1]+y_min)*log10(ps[i+1,0]/k_min)
            if (k_max >= k) & (k_max <= k_s):
                y_max = ps[i,1] + (ps[i+1,1] - ps[i,1])/(ps[i+1,0] - ps[i,0])*(k_max - ps[i,0])
                var += 0.5*(ps[i+1,1]+y_max)*log10(k_max/ps[i,0])
            else:
                var += 0.5*(ps[i,1]+ps[i+1,1])*log10(ps[i+1,0]/ps[i,0])
    return var

def map_Correa2015(logmass,boxsize,z):
    radius = cal_radius(logmass)
    z_f = -0.0064*logmass**2 + 0.0237*logmass +1.8837
    q = 4.137*z_f**(-0.9476)
    radius_q = cal_radius(logmass-log10(q))
    f = []
    D = -1.5/e
    ss = []
    sq = []
    for i in range(len(logmass)):
        S = cal_variance(2*pi/boxsize,2*pi/radius[i])
        ss.append(S)
        Sq = cal_variance(2*pi/boxsize,2*pi/radius_q[i])
        sq.append(Sq)
        f.append(1./sqrt(Sq - S))
    f = array(f)
    alpha = f*(1.686*sqrt(2./pi)*D + 1.)
    beta = -1.*f
    M0 = 10.**logmass
    Mz = M0*(1+z)**alpha * exp(beta*z)
    return (log10(M0),log10(Mz))
    
def Correa2015(M0,boxsize,z):
    logmass = log10(M0)
    radius = cal_radius(logmass)
    z_f = -0.0064*logmass**2 + 0.0237*logmass +1.8837
    q = 4.137*z_f**(-0.9476)
    radius_q = cal_radius(logmass-log10(q))
    D = -1.5/e
    S = cal_variance(2*pi/boxsize,2*pi/radius)
    Sq = cal_variance(2*pi/boxsize,2*pi/radius_q)
    f = 1./sqrt(Sq-S)
    alpha = f*(1.686*sqrt(2./pi)*D + 1.)
    beta = -1.*f
    Mz = M0*(1+z)**alpha * exp(beta*z)
    return Mz

def mz_Correa2015(m0,z0,zlist,boxsize):
    min_lgmass = 5.0
    max_lgmass = 16.0
    step = 0.01
    logmass = arange(min_lgmass,max_lgmass,step)
    m = map_Correa2015(logmass,boxsize,z0)
    M0 = 10.**m[0][find_closest(m[1],m0)]
    mz = zeros(len(zlist))
    for i in range(len(zlist)):
        mz[i] = Correa2015(M0,boxsize,zlist[i])
    return mz
def main(argv):
    z = 6.0
    boxsize = 112.0 #Mpc/h
    folder = "/scratch/01937/cs390/test_4_1440_112/trees/treedata/"
    #folder = "/Volumes/Backup/Work/test_4_1440_112/treedata/"
    snapfile = "/scratch/01937/cs390/test_4_1440_112/trees/snap_z3.txt"
    #snapfile = "/Volumes/Backup/Work/test_4_1440_112/snap_z3.txt"
    z_list_lgal = loadtxt(snapfile)
    firstfile  = 0
    lastfile = 7
    gadget_m_conv = 1.e10
    hubble_h = 0.7
    m6 = arange(9.5,12.5,0.5)
    zlist = arange(6.0,18.0,0.1)
    rc('text', usetex=True)
    fig = figure()
    ax = fig.add_subplot(111)
    for t_m6 in m6:
        mz = mz_Correa2015(t_m6,z,zlist,boxsize)
        ax.plot(zlist,log10(mz))
    (nTrees,nHalos,nTreeHalos,output_Halos,output_HaloIDs) = read_lgal_input_fulltrees_withids(folder,lastsnap,firstfile,lastfile,verbose=True)
    rootindex = numpy.cumsum(nTreeHalos)-nTreeHalos
    for t_m6 in m6:
        mass = zeros(len(z_list_lgal),dtype=float64)
        count = zeros(len(z_list_lgal),dtype=int64)
        mask = ones(len(z_list_lgal),dtype=int32)
        l_m = t_m6-0.1
        r_m = t_m6+0.1
        r_list = numpy.where((numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv) <=r_m) & (numpy.log10(output_Halos[rootindex]['M_Crit200']*gadget_m_conv) >=l_m))[0]
        for i in r_list:
            root = rootindex[i]
            M0 = output_Halos[root]["M_Mean200"]*Gadget2Msun
            nexthaloid = output_Halos[root]['FirstProgenitor']
            mass[output_Halos[root]["SnapNum"]] += M0
            count[output_Halos[root]["SnapNum"]] += 1
            while nexthaloid > -1:
                nexthalo = output_Halos[root+nexthaloid]
                nextprogid = nexthalo['NextProgenitor']
                if nextprogid == -1:
                    mass[nexthalo["SnapNum"]] += nexthalo["M_Mean200"]*Gadget2Msun
                    count[nexthalo["SnapNum"]] += 1
                    nexthaloid = nexthalo['FirstProgenitor']
                else:
                    nexthaloid = -1
        
        mask = count > count[len(z_list_lgal)-1]/5
        mass = mass*mask
        ax.plot(z_list_lgal,log10(mass/count))
    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$\log(hM/M_\odot)$")
    savefig("test.pdf")
    return 0
if __name__ == "__main__":
    main(sys.argv)
