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
    nbin = 6#0
    nbinM = 10#0
    
    maxmag = -20
    minmag = -29
    mag = np.linspace(minmag,maxmag,nbinM)    
    minz = 0.1
    maxz = 5
    z = np.linspace(minz,maxz,nbin)
    
    minv = 10**1.6
    maxv = 10**2.6
    vel = np.linspace(minv,maxv,nbin)
    
    # writing the script for galfic once and for all
    gl.write_script()

    result = [[] for i in range(5)]
    
    for zl in z:
        zrange = z[np.where(z > zl)]
        for zs in zrange:
            for v in vel:
                gl.write_initfile(v,zl,zs)
                os.system('./script_gl > /dev/null 2>&1')
                sigma = gl.analyse_output(mag,zs,zl,v)
                for i in range(nbinM):
                    result[0].append(zs)
                    result[1].append(zl)
                    result[2].append(v)
                    result[3].append(mag[i])
                    result[4].append(sigma[i])

    np.savez('crosssection.npz',x=result)

def test():
    """Test the smoothness of the function in all directions"""

    file = 'crosssection.dat'
    f = open(file,'r')
    lines = f.readlines()
    nline = len(lines)
    points = np.zeros(shape=(nline,4))
    sigtable = np.zeros(nline)
    for i in range(nline):
        points[i,0] = float(lines[i].split()[0])
        points[i,1] = float(lines[i].split()[1])
        points[i,2] = float(lines[i].split()[2])
        points[i,3] = float(lines[i].split()[3])
        sigtable[i] = float(lines[i].split()[4])

    nbin = 60
    npts = nline/nbin

    # checking lensing cross section against magnitude
    '''
    for i in range(npts):
        plt.plot(points[i*nbin:(i+1)*nbin,3],sigtable[i*nbin:(i+1)*nbin])
        plt.show()
    '''
    npts = npts/nbin

    # checking lensing cross section against velocity dispersion
    '''
    for i in range(nline):
        mask, = np.where((points[:,1]==points[i,1])&(points[:,0]==points[i,0])\
                        &(points[:,3]==points[i,3]))
        vel = points[mask,2]
        sigma = sigtable[mask]
        plt.plot(vel,sigma)
        plt.show()
    '''

    # checking lensing cross section against lens redshift
    #'''
    for i in range(3000,nline):
        mask, = np.where((points[:,1]==points[i,1])&(points[:,2]==points[i,2])\
                        &(points[:,3]==points[i,3]))
        print mask
        zl = points[mask,0]
        sigma = sigtable[mask]
        plt.plot(zl,sigma)
        plt.show()
    #'''

    # checking lensing cross section against source redshift
    for i in reversed(range(nline)):
        mask, = np.where((points[:,0]==points[i,0])&(points[:,2]==points[i,2])\
                        &(points[:,3]==points[i,3]))
        print mask
        zs = points[mask,1]
        sigma = sigtable[mask]
        plt.plot(zs,sigma)
        plt.show()
