################################################################################
import matplotlib.pyplot as plt
import numpy as np
import math as math
from scipy.interpolate import griddata

def quasarluminosity(M,z):
    """Compute dPhi/dM for QSOs at redshift z.

    The magnitude here is the absolute i-band magnitude.
    Arguments:
    M -- absolute magnitude at which to do the computation.
    z -- redshift of the quasars.
    Outputs:
    diff -- dPhi/dM in absolute i-band magnitude.
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

    # the result is in h**3/mag/Mpc**3. Is it necessary to multiply it by h**3 
    # before using it???

    return diff

def velocitydispersion(v):
    """Compute the velocity function of early-type galaxies.

    Arguments:
    v -- velocity at which to compute the number of galaxies.
    Outputs:
    diff -- number density of lenses at v.
    """
    # parameters
    Phi_star = 8e-3 # units (h/Mpc)^3
    v_star = 161. # units km/s
    alpha = 2.32
    beta = 2.67

    x = v/v_star
    norm = Phi_star*beta/math.gamma(alpha/beta)

    diff = norm*(x**alpha)*np.exp(-(x**beta))/v

    # There seems to be a factor 3 difference between this and 
    # the result in Choi et al 2007. WHY???
    # It also appears to drop much too quickly.

    return diff

def lensingcrosssection(v,zl,zs,M):
    """Find the biased lensing cross-section for v, zl, zs.

    Interpolate from a pre-computed table.
    Arguments:
    v -- velocity dispersion of the lens.
    zl -- redshift of the lens.
    zs -- redshift of the source.
    M -- magnitude of the source.
    Outputs:
    sigma -- biased lensing cross-section, precomputed with gravlens.
    """
    file = 'crosssection.dat'
    f = open(file,'r')
    zltab = []
    zstab = []
    Mtab = []
    vtab = []
    sigtable = []
    for line in f:
        zstab.append(float(line.split()[0]))
        zltab.append(float(line.split()[1]))
        vtab.append(float(line.split()[2]))
        Mtab.append(float(line.split()[3]))
        sigtable.append(float(line.split()[4]))

    sigma = griddata((zstab,zltab,vtab,Mtab), sigtable, (zs,zl,v,M),\
                     method='nearest')

    return sigma

def lensingprobability(zs,M):
    """Compute the lensing probability for a QSO at redshift zs, magnitude M.

    Arguments:
    zs -- redshift of the source quasar.
    M -- magnitude of the lensed quasar.
    Outputs:
    p -- probability that the quasar is lensed.
    """

    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    h = 0.72 # *100 km/s.Mpc = Hubble's constant
    c = 299792.458 # km/s

    # conversion from arcseconds to radian
    arcsectorad = 4.85e-6

    # integration over Theta
    nbin_theta = 100
    theta_min = 1. # PSF of JPAS is about 1 arcsecond
    theta_max = 4. # we don't expect new lenses 
                  # with separation above 4 arcseconds
    theta = theta_min+np.arange(nbin_theta)*(theta_max-theta_min)/(nbin_theta-1)
    # theta in radians
    theta = arcsectorad*theta
    #dtheta = theta[1]-theta[0]

    # this array will contain the function to integrate over theta
    integrand = []

    # for the integration over redshift zl
    # quantities independant from theta
    nbin = 100
    # in order to avoid nans (dls = 0), zl must never be equal to zs
    zl = np.arange(nbin)*np.float(zs)/nbin

    hl = h*100*np.sqrt(Omega_m*(1+zl)**3+Omega_L*(1+zl)**(-3*(1+w)))
    volume = (1+zl)**2*c/hl

    ds = distance(zs)
    dls = []
    for z in zl:
        dls.append(distance(zs,z))
    dls = np.array(dls)

    for t in theta:
        # The relation between the velocity dispersion
        # and the angular separation comes (I assume) from the 
        # Einstein radius relation
        # Theta_E = 4pi(v/c)**2*D_ls/D_s
        # Therefore dv/dTheta is given by the following relation
        diffvtheta = c/2*np.sqrt(ds/(4*np.pi*dls*t))

        # velocity corresponding to the given theta
        # same units as c, km/s
        v = c*np.sqrt(ds*t/(4*np.pi*dls))
        # dv corresponding to dtheta
        #dv = diffvtheta*dtheta
        vdisp = velocitydispersion(v)#,dv)

        # biased lensing cross-section
        # can lensingcrosssection take an array for zl???
        # It seems that yes
        sigma_l = []
        for z in zl:
            sigma_l.append(lensingcrosssection(v,z,zs,M))
        sigma_l = np.array(sigma_l)

        # result of the integration over the redshift
        integrand_z = volume*vdisp*diffvtheta*sigma_l 
        integral = np.trapz(integrand_z,zl)

#        print integral
#        plt.plot(zl,integrand_z)
#        plt.show()

        integrand.append(integral)

    integral_theta = np.trapz(integrand,theta)

#    print integral_theta
#    plt.plot(theta,integrand)
#    plt.show()

    return integral_theta

def distance(z1,z2=0):
    """Compute the angular diameter distance between redshifts z1 and z2.

    If z2 has no explicit value, computes the distance between 0 and z1.
    Arguments:
    z1 -- final redshift. Mandatory argument. 
    z2 -- initial redshift. Optional argument, default value 0.
    Outputs:
    D -- angular diameter distance from z1 to z2.
    """
    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    h = 0.72 # *100 km/s.Mpc = Hubble's constant
    c = 299792.458 # km/s

    astart = 1./(1+z1)
    aend = 1./(1+z2)

    # fonction to integrate
    nbin = 1000
    a = astart+np.arange(nbin)*(aend-astart)/(nbin-1)
    rgral = np.sqrt(Omega_m*a+Omega_L*a**(3*(1+w)+4))

    d = np.trapz(1./rgral,x=a)*astart*c/(100*h) #in Mpc

    return d

def lensingbyredshift(zs,Mmax):
    """Compute the expected number of lenses in a slice of redshift.

    Arguments:
    zs -- source redshift.
    Mmax -- maximum magnitude of the survey.
    Outputs:
    nlens -- dN/dzs, number of expected lenses by interval of source redshift.
    """
    nbin = 200
    Mmin = -30. 
    # cosmological parameters
    Omega_m = 0.3
    Omega_L = 1-Omega_m
    w = -1
    c = 299792.458 # km/s
    h = 0.72 

    # the density of quasars there should be below 1e-9 per
    # magnitude per Mpc^3.
    # M is the magnitude of the lensed quasar.
    M = np.arange(nbin)*(Mmax-Mmin)/(nbin-1)+Mmin
    
    # luminosity function
    luminosity = quasarluminosity(M,zs)
    
    # volume factor
    area = 2.5 # survey area in steradian. 8000 to 9000 square degrees
    dist = distance(zs)
    hs = h*100*np.sqrt(Omega_m*(1+zs)**3+Omega_L*(1+zs)**(-3*(1+w)))
    volume = area*dist**2*(1+zs)**2*c/hs
    
    p = []
    for mag in M:
        p.append(lensingprobability(zs,mag))
    
    integrand = luminosity*volume*p

    nlens = np.trapz(integrand,M)
    
    return nlens

def numberoflenses(zmax,Mmax):
    """Compute the total number of lenses expected in the JPAS survey.

    Arguments:
    zmax -- maximum redshift at which objects will be detected. ~6
    Mmax -- maximum magnitude at which objects will be detected.
    Outputs:
    N -- number of lenses expected to be observed by JPAS.
    """
    nbin = 20
    zmin = 0.1 #zs must not be 0
    zs = zmin+np.arange(nbin)*np.float(zmax-zmin)/(nbin-1)
    
    integrand = []
    k = 0
    for z in zs:
        integrand.append(lensingbyredshift(z,Mmax))
        k = k+1
        print k,'/',nbin
    N = np.trapz(integrand,zs)

    return N
