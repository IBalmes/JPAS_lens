################################################################################
import matplotlib.pyplot as plt
import numpy as np
import math as math

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
    if z < 3:
        alpha = -3.31
    else:
        alpha = -2.58

    h = 0.72 # Hubble's constant value taken in Oguri & Marshall

    # value of the break absolute magnitude at this redshift
    f = np.exp(zeta*z)*(1+np.exp(ksi*z_star))/(np.sqrt(np.exp(ksi*z)) \
                                               +np.sqrt(np.exp(ksi*z_star)))**2
    M_star = -20.90+5*np.log10(h)-2.5*np.log10(f)

    diff = np.log10(Phi_star/(10**(0.4*(alpha+1)*(M-M_star)) \
                     +10**(0.4*(beta+1)*(M-M_star))))
    
    return diff

def velocitydispersion(v,dv):
    """Compute the velocity function of early-type galaxies.

    Arguments:
    v -- velocity at which to compute the number of galaxies.
    dv -- velocity interval.
    Outputs:
    diff -- dn/dv.
    """
    # parameters
    Phi_star = 8e-3 # units (h/Mpc)^3
    v_star = 161 # units km/s
    alpha = 2.32
    beta = 2.67

    diff = Phi_star*(v/v_star)**alpha*np.exp(-(v/v_star)**beta)* \
        beta/math.gamma(alpha/beta)*dv/v

    return diff

def lensingcrosssection():
    """ """

def lensingprobability(zs,M):
    """Compute the lensing probability for a QSO at redshift zs, magnitude M.

    Arguments:
    zs -- redshift of the source quasar.
    M -- magnitude of the quasar.
    Outputs:
    p -- probability that the quasar is lensed.
    """
    
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

    astart = 1/(1+z1)
    aend = 1/(1+z2)

    # fonction to integrate
    
    rgral = np.sqrt(Omega_m*a+Omega_L*a**(3*(1+w)+4))

def lensingbyredshift(zs,Mmax):
    """Compute the expected number of lenses in a slice of redshift.

    Arguments:
    zs -- source redshift.
    Outputs:
    nlens -- dN/dzs, number of expected lenses by interval of source redshift.
    """
    nbin = 200
    Mmin = -30 
    # the density of quasars there should be below 1e-9 per
    # magnitude per Mpc^3.
    M = np.arange(nbin)*(Mmax-Mmin)/(nbin-1)+Mmin
    
    # luminosity function
    luminosity = quasarluminosity(M,zs)

    # volume factor
    area = 2.5 # survey area in steradian. 8000 to 9000 square degrees
    dist = distance(zs)
    volume = area*dist**2*c*(1+zs)**3.....
    
    p = lensingprobability(zs,M)

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
    nbin = 200
    zs = np.arange(nbin)*zmax/(nbin-1)
    
    integrand = lensingbyredshift(zs,Mmax)
    N = np.trapz(integrand,zs)

    return N
