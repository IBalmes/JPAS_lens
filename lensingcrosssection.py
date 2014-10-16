################################################################################
import numpy as np
import matplotlib.pyplot as plt
import draw_glafic as gl
import lens as l
import os

def main():
    """Fills the table of biased lensing cross-section.

    Varies magnitude and redshift of the source, velocity dispersion and
    redshift of the lens to fill a complete table of biased lensing 
    cross-section, to be used by lens.py.
    Arguments:
    Outputs:
    """
    nbin = 10#000

    maxmag = -20
    minmag = -29
    mag = [-25]#np.linspace(minmag,maxmag,nbin)

    maxz = 5
    minz = 0.1
    z = np.linspace(minz,maxz,nbin)

    maxv = 10**2.6
    minv = 10**1.6
    vel = np.linspace(minv,maxv,nbin)

    # writing the script for galfic once and for all
    gl.write_script()

    c = 299792.458 # km/s
    arcsectorad = 4.85e-6
    thetaE = []
    sig = []

    file = 'crosssection.dat'
    f = open(file,'w')
    for zl in z:
        zrange = z[np.where(z > zl)]
        for zs in zrange:
            for v in vel:
                gl.write_initfile(v,zl,zs)
                os.system("./script_gl")
                thetaE.append(4*np.pi*(v/c)**2*l.distance(zs,zl)/l.distance(zs))
                for M in mag:
                    sigma = gl.analyse_output(M,zs)
                    sig.append(sigma)
                    line = str(zs)+' '+str(zl)+' '+str(v)+' '+str(M)+' '\
                           +str(sigma)+'\n'
                    f.write(line)

    f.close()
                    
    plt.plot(thetaE,sig,'.')
    plt.show()
