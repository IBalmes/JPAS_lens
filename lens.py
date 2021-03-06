################################################################################
import matplotlib.pyplot as plt
import numpy as np
import math as math
from scipy.interpolate import griddata

"""Note:
A quasar absolute magnitude is typically between -29 and -20.
Here is a table translating that to apparent magnitude, function of redshift:
redshift // min app mag (-29) // max app mag (-20)
0 // -29 // -20
0.1 // 8.8 // 17.8
0.3 // 10.8 // 19.8
0.5 // 11.4 // 20.4 
1 // 12 // 21
2 // 12.1 // 21.1
3 // 11.9 // 20.9
"""

def quasarluminosity(M,z):
    """Compute dPhi/dM for QSOs at redshift z.

    The magnitude here is the absolute i-band magnitude.
    Arguments:
    M -- ABSOLUTE magnitude at which to do the computation.
    z -- redshift of the quasars.
    Outputs:
    diff -- dPhi/dM in absolute i-band magnitude in (h/Mpc)**3.
    """
    # parameters of the fonction
    Phi_star = 5.34e-6 # units (h/Mpc)^3
    zeta = 2.98
    ksi = 4.05
    z_star = 1.60
    beta = -1.45
    if z > 3:
        alpha = -2.58
    else:
        alpha = -3.31

    h = 0.72 # Hubble's constant value taken in Oguri & Marshall

    # value of the break absolute magnitude at this redshift
    f = np.exp(zeta*z)*(1+np.exp(ksi*z_star))/(np.sqrt(np.exp(ksi*z)) \
                                               +np.sqrt(np.exp(ksi*z_star)))**2
    M_star = -20.90+5*np.log10(h)-2.5*np.log10(f)

    diff = Phi_star/(10**(0.4*(alpha+1)*(M-M_star)) \
                     +10**(0.4*(beta+1)*(M-M_star)))

    # the result is in (h/Mpc)**3.

    return diff


def velocitydispersion(v):
    """Compute the velocity function of early-type galaxies.

    Arguments:
    v -- velocity at which to compute the number of galaxies in km/s.
    Outputs:
    diff -- number density of lenses at v in 1/(Mpc/h)^3/(km/s).
    """
    # parameters
    Phi_star = 8e-3 # units (h/Mpc)^3
    v_star = 161. # units km/s
    alpha = 2.32
    beta = 2.67

    x = v/v_star
    norm = Phi_star*beta/math.gamma(alpha/beta) # units (h/Mpc)^3 

    diff = norm*(x**alpha)*np.exp(-(x**beta))/v # units 1/(Mpc/h)^3/(km/s)

    # There seems to be a factor 3 difference between this and 
    # the result in Choi et al 2007. WHY???
    # It also appears to drop much too quickly.
    # -> Choi et al 2007 plot the log function, so diff*v*np.log(10).

    return diff


def lensingcrosssection(v,zl,zs,M):
    """Find the biased lensing cross-section for v, zl, zs.

    Interpolate from a pre-computed table.
    Arguments:
    v -- velocity dispersion of the lens.
    zl -- redshift of the lens.
    zs -- redshift of the source.
    M -- ABSOLUTE magnitude of the source.
    Outputs:
    sigma -- biased lensing cross-section, precomputed with gravlens, in square
    radians
    """
    result = np.load('crosssection.npz')

    zstab = result['x'][0]
    zltab = result['x'][1]
    vtab = result['x'][2]
    Mtab = result['x'][3]
    sigtable = result['x'][4]

    sigma = griddata((zstab,zltab,vtab,Mtab), sigtable, (zs,zl,v,M),\
                     method='nearest')

    return sigma


def lensingprobability(zs,Mapp):
    """Compute the lensing probability for a QSO at redshift zs, magnitude M.

    Arguments:
    zs -- redshift of the source quasar.
    Mapp -- APPARENT magnitude of the lensed quasar.
    Outputs:
    p -- probability that the quasar is lensed.
    """

    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    h = 0.72 # *100 km/s/Mpc = Hubble's constant
    c = 299792.458 # km/s

    # integration over velocity dispersion
    nbin_v = 30
    minv = 10**1.6
    maxv = 10**2.6
    vel = np.linspace(minv,maxv,nbin_v)

    # this array will contain the function to integrate over v
    integrand = []

    # for the integration over redshift zl
    # quantities independant from v
    nbin = 30
    # in order to avoid nans (dls = 0), zl must never be equal to zs
    zl = np.arange(nbin)*np.float(zs)/nbin

    hl = 100*np.sqrt(Omega_m*(1+zl)**3+Omega_L*(1+zl)**(-3*(1+w))) 
    # hl in km/s/(Mpc/h)
    volume = (1+zl)**2*c/hl # in Mpc/h

    ds = distance(zs)
    dls = []
    for z in zl:
        dls.append(distance(zs,z))
    dls = np.array(dls)

    for v in vel:
        vdisp = velocitydispersion(v) # km/s/(Mpc/h)**3

        # biased lensing cross-section
        # can lensingcrosssection take an array for zl???
        # It seems that yes
        sigma_l = []
        for i in range(len(zl)):
            dist = distance(zs)
            # converting magnitude to absolute magnitude
            Mabs = Mapp-5*(np.log10(dist*1e6)-1)
            # a is in square radian
            a = lensingcrosssection(v,zl[i],zs,Mabs)
            # sigma_l is in Mpc/h
            sigma_l.append(a*dist**2)
        sigma_l = np.array(sigma_l)

        # result of the integration over the redshift
        integrand_z = volume*vdisp*sigma_l 
        integral = np.trapz(integrand_z,zl)

        integrand.append(integral)

    integral_vel = np.trapz(integrand,vel)

    return integral_vel


def distance(z1,z2=0):
    """Compute the angular diameter distance between redshifts z1 and z2 in Mpc.

    If z2 has no explicit value, computes the distance between 0 and z1.
    Arguments:
    z1 -- final redshift. Mandatory argument. 
    z2 -- initial redshift. Optional argument, default value 0.
    Outputs:
    d -- angular diameter distance from z1 to z2 in Mpc/h.
    """
    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    h = 0.72 # *100 km/s/Mpc = Hubble's constant
    c = 299792.458 # km/s

    astart = 1./(1+z1)
    aend = 1./(1+z2)

    # fonction to integrate
    nbin = 1000
    a = astart+np.arange(nbin)*(aend-astart)/(nbin-1)
    rgral = np.sqrt(Omega_m*a+Omega_L*a**(3*(1+w)+4))

    d = np.trapz(1./rgral,x=a)*astart*c #in Mpc/h

    return d


def lensingbyredshift(zs,Mmax):
    """Compute the expected number of lenses in a slice of redshift.

    Arguments:
    zs -- source redshift.
    Mmax -- maximum APPARENT magnitude of the survey. Typically of order 23.
    Outputs:
    nlens -- dN/dzs, number of expected lenses by interval of source redshift.
    """
    nbin = 10#200
    Mmin = -30. 
    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    c = 299792.458 # km/s
    h = 0.72 

    # the density of quasars there should be below 1e-9 per
    # magnitude per Mpc^3.
    # M is the APPARENT magnitude of the lensed quasar.
    Mapp = np.arange(nbin)*(Mmax-Mmin)/(nbin-1)+Mmin
    # the magnitude has to be converted to ABSOLUTE
    dist = distance(zs) # in Mpc/h
    Mabs = Mapp-5*(np.log10(dist*1e6)-1)

    # luminosity function
    luminosity = quasarluminosity(Mabs,zs) # in (Mpc/h)**-3
    
    # volume factor
    area = 2.5 # survey area in steradian. 8000 to 9000 square degrees
    hs = 100*np.sqrt(Omega_m*(1+zs)**3+Omega_L*(1+zs)**(-3*(1+w)))
    volume = area*dist**2*(1+zs)**2*c/hs # in Mpc/h
    
    p = []
    for mag in Mapp:
        # lensingprobability takes the APPARENT magnitude as input.
        p.append(lensingprobability(zs,mag))
    
    integrand = luminosity*volume*p

    nlens = np.trapz(integrand,Mapp)
    
    return nlens


def numberoflenses(zmax,Mmax):
    """Compute the total number of lenses expected in the JPAS survey.

    Arguments:
    zmax -- maximum redshift at which objects will be detected. ~6
    Mmax -- maximum APPARENT magnitude at which objects will be detected.
    Outputs:
    N -- number of lenses expected to be observed by JPAS.
    """
    nbin = 10#20
    zmin = 0.1 #zs must not be 0
    zs = zmin+np.arange(nbin)*np.float(zmax-zmin)/(nbin-1)
    
    integrand = []
    k = 0
    for z in zs:
        integrand.append(lensingbyredshift(z,Mmax))
        k = k+1
        print k,'/',nbin

    plt.plot(zs,integrand)
    plt.show()

    N = np.trapz(integrand,zs)

    return N
