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
    #c = 299792.458 # km/s
    #arcsectorad = 4.85e-6
    # converting thetaE to radian
    #thetaE = thetaE*arcsectorad
    #v = c*np.sqrt(thetaE/(4*np.pi)*l.distance(zs)/l.distance(zs,zl))
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
        f.write('startup 1 0 1\n')
        #line = 'lens powpot '+str(zs)+' 0 0 '+str(e[i])+' 0 '+str(thetaE)+\
        #       ' 2\n'
        line = 'lens nfwpot 1e14 0 0 '+str(e[i])+' 0 5 0\n'
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
    file = 'script_test'
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
    os.system('./script_test > /dev/null 2>&1')
    
def analyze(nlens):
    """Analyze the outputs of glafic. Compute the separation angle.

    Arguments:
    nlens -- number of lenses to run.
    Outputs:
    sep -- array of all separation angle between images.
    """
    # array containing the separation angle
    sep = []
    nmult = 0
    for i in range(nlens):
        file = 'init_test/lens_'+str(i)+'.out'
        f = open(file,'r')
        line = f.readline()
        image_number = float(line.split()[2])
        if image_number == 2:
            # double lens case
            image1 = f.readline()
            image2 = f.readline()
            angle12 = computeangle(image1,image2)
            sep.append(angle12)
            nmult = nmult+1
        if image_number == 3:
            # naked cusp case
            image1 = f.readline()
            image2 = f.readline()
            image3 = f.readline()
            angle12 = computeangle(image1,image2)
            sep.append(angle12)
            angle13 = computeangle(image1,image3)
            sep.append(angle13)
            angle23 = computeangle(image2,image3)
            sep.append(angle23)
            nmult = nmult+1
        if image_number == 4:
            # quad lens case
            image1 = f.readline()
            image2 = f.readline()
            image3 = f.readline()
            image4 = f.readline()
            angle12 = computeangle(image1,image2)
            sep.append(angle12)
            angle13 = computeangle(image1,image3)
            sep.append(angle13)
            angle14 = computeangle(image1,image4)
            sep.append(angle14)
            angle23 = computeangle(image2,image3)
            sep.append(angle23)
            angle24 = computeangle(image2,image4)
            sep.append(angle24)
            angle34 = computeangle(image3,image4)
            sep.append(angle34)
            nmult = nmult+1

    print nmult

    sep = np.array(sep)

    return sep


def computeangle(image1,image2):
    """Compute the separation angle between two images.

    Arguments:
    image1,image2 -- strings, outputs of glafic.
    Outputs:
    angle -- separation angle between the two images.
    """
    x1 = float(image1.split()[2])
    y1 = float(image1.split()[5])
    x2 = float(image2.split()[2])
    y2 = float(image2.split()[5])
    angle = np.sqrt((x1-x2)**2+(y1-y2)**2)

    return angle

def main(nlens,thetaE,zl,zs):
    """Run all glafic test and graph the histogram of separation angles.

    Arguments:
    nlens -- number of lenses to use.
    thetaE -- Einstein radius in arcseconds
    zl -- lens redshift
    zs -- source redshift
    Outputs:
    """
    write_script(nlens)
    # running with ellipticity in the lens
    run(nlens,thetaE,zl,zs)
    angle_ell = analyze(nlens)
    np.savez('angle_ell.npz',x=angle_ell)

    # running without ellipticity
    run(nlens,thetaE,zl,zs,ell=False)
    angle_noell = analyze(nlens)
    np.savez('angle_noell.npz',x=angle_noell)

    angles = [angle_ell,angle_noell]

    plt.hist(angles,50,histtype='bar')
    plt.show()
