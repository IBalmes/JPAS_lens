################################################################################
import numpy as np
import matplotlib.pyplot as plt
import draw_glafic as gl
import lens as l
import os
#import multiprocessing

def main():
    """Fills the table of biased lensing cross-section.
    
    Varies magnitude and redshift of the source, velocity dispersion and
    redshift of the lens to fill a complete table of biased lensing 
    cross-section, to be used by lens.py.
    Arguments:
    Outputs:
    """
    nbin = 20
    
    maxmag = -20
    minmag = -29
    mag = np.linspace(minmag,maxmag,nbin)    
    minz = 0.1
    maxz = 5
    z = np.linspace(minz,maxz,nbin)
    
    minv = 10**1.6
    maxv = 10**2.6
    vel = np.linspace(minv,maxv,nbin)
    
    # writing the script for galfic once and for all
    gl.write_script()
    
    file = 'crosssection.dat'
    f = open(file,'w')
    for zl in z:
        zrange = z[np.where(z > zl)]
        for zs in zrange:
            for v in vel:
                gl.write_initfile(v,zl,zs)
                os.system('./script_gl > /dev/null 2>&1')
                for M in mag:
                    sigma = gl.analyse_output(M,zs)
                    line = str(zs)+' '+str(zl)+' '+str(v)+' '+str(M)+' '\
                           +str(sigma)+'\n'
                    f.write(line)

    f.close()
