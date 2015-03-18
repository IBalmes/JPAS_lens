################################################################################
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
import draw_glafic as gl
import lens as l
import os

def write_initfile(nlens,thetaE,zl,zs,ell=True):
    """Write the initial files for gravlens for a given Einstein radius.

    Arguments:
    nlens -- number of lenses
    thetaE -- Einstein radius in arcseconds
    zl -- lens redshift
    zs -- source redshift
    ell -- if True, lenses are elliptical. If False, spherical. Default True.
    Outputs:
    tonnes of files
    """
    c = 299792.458 # km/s
    arcsectorad = 4.85e-6
    # converting thetaE to radian
    thetaE = thetaE*arcsectorad
    v = c*np.sqrt(thetaE/(4*np.pi)*l.distance(zs)/l.distance(zs,zl))
    e,gamma,theta_g,x,y = gl.draw_lensparam(nlens,thetaE)
    if not ell:
        e = e*0

    for i in range(nlens):
        file = 'init_test/lens_'+str(i)+'.in'
        f = open(file,'w')
        f.write('#parameters\n')
        f.write('omega 0.26\n')
        f.write('lambda 0.74\n')
        f.write('weos -1\n')
        f.write('hubble 0.72\n')
        line = 'zl '+str(zl)+'\n'
        f.write(line)
        f.write('\n')
        f.write('#defining the lens model\n')
        f.write('startup 2 0 1\n')
        line = 'lens sie '+str(v)+' 0 0 '+str(e[i])+' 0 0 0\n'
        f.write(line)
        line = 'lens pert '+str(zs)+' 0 0 0 0 0 0\n'
        f.write(line)
        line = 'point '+str(zs)+' '+str(x[i])+' '+str(y[i])+'\n'
        f.write(line)
        f.write('end_startup\n')
        f.write('\n')
        f.write('start_command\n')
        f.write('findimg\n')
        f.write('\n')
        f.write('quit')
        f.close()

def write_script(nlens):
    """Write the script to run gravlens with the initial files.

    Arguments:
    nlens -- number of lenses to run.
    """
    file = 'script_gl'
    f = open(file,'w')
    for i in range(nlens):
        line = './glafic init_test/lens_'+str(i)+'.in > init_test/lens_'\
               +str(i)+'.out\n'
        f.write(line)

def run(nlens,thetaE,zl,zs,ell=True):
    """Creates the init files and runs the glafic script.

    Arguments:
    nlens -- number of lenses to run.
    thetaE -- Einstein radius in arcseconds
    zl -- lens redshift
    zs -- source redshift
    ell -- if True, lenses are elliptical. If False, spherical. Default True.
    """
    write_initfile(nlens,thetaE,zl,zs,ell)
    os.system('./script_gl > /dev/null 2>&1')
    
def analyze(nlens):
    """Analyze the outputs of glafic. Compute the separation angle.

    Arguments:
    nlens -- number of lenses to run.
    Outputs:
    sep -- array of all separation angle between images.
    """
    # array containing the magnification, and number of multiply imaged lenses
    mu = []
    mutmp = []
    nmult = 0
    for i in range(start,start+nlens):
        file = 'init_gl/lens_'+str(i)+'.out'
        f = open(file,'r')
        line = f.readline()
        image_number = float(line.split()[2])
        if image_number == 2:
            # double lens case
            image1 = f.readline()
            image2 = f.readline()
            mutmp.append(abs(float(image1.split()[8])))
            mutmp.append(abs(float(image2.split()[8])))
            mutmp.sort()
            #if mutmp[0]/mutmp[1]>0.1: 
            # in the end, no condition is necessary there
            nmult = nmult+1
            mu.append(mutmp[0]) # faintest image
            mutmp = []
        if image_number == 3:
            # naked cusp case
            image1 = f.readline()
            image2 = f.readline()
            image3 = f.readline()
            mutmp.append(abs(float(image1.split()[8])))
            mutmp.append(abs(float(image2.split()[8])))
            mutmp.append(abs(float(image3.split()[8])))
            mutmp.sort()
            nmult = nmult+1
            mu.append(mutmp[0]) # third brightest image
            mutmp = []
        if image_number == 4:
            # quad lens case
            image1 = f.readline()
            image2 = f.readline()
            image3 = f.readline()
            image4 = f.readline()
            mutmp.append(abs(float(image1.split()[8])))
            mutmp.append(abs(float(image2.split()[8])))
            mutmp.append(abs(float(image3.split()[8])))
            mutmp.append(abs(float(image4.split()[8])))
            mutmp.sort()
            nmult = nmult+1
            mu.append(mutmp[1]) # third brightest image
            mutmp = []
