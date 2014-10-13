################################################################################
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand

def draw_lensparam(n):
    """Draw ellipticity, shear, shear angle and source positions.

    Draw n sets of lens parameters from the given distributions (see Oguri &
    Marshall 2010)
    Arguments:
    n -- number of sets to draw
    Outputs:
    e -- ellipticity
    gamma -- shear strength
    theta_g -- shear angle
    x,y -- position of the source in arcseconds
    """
    e = rand.normal(0.3,0.16,n)
    for i in range(n):
        while ((e[i]<0.0) or (e[i]>0.9)):
            e[i] = rand.normal(0.3,0.16)

    gamma = rand.lognormal(0.05,0.2,n)

    theta_g = rand.uniform(0,180,n)

    x = rand.uniform(-1,1,n)
    y = rand.uniform(-1,1,n)

    return e,gamma,theta_g,x,y

def write_initfile(theta):
    """Write the initial files for gravlens for a given Einstein radius.

    Arguments:
    theta -- Einstein radius
    Outputs:
    tonnes of files
    """
    nlens = 500
    e,gamma,theta_g,x,y = draw_lensparam(nlens)

    for i in range(nlens):
        file = 'init_gl/lens_'+str(i)+'.in'
        f = open(file,'w')
        f.write('#defining the lens model\n')
        f.write('startup 1 1\n')
        line = 'alpha '+str(theta)+' 0 0 '+str(e[i])+' 0 '+ \
            str(gamma[i])+' '+str(theta_g[i])+' 0 0 1\n'
        f.write(line)
        f.write('0 0 0 0 0 0 0 0 0 0\n')
        f.write('\n')
        lineb = 'findimg '+str(x[i])+' '+str(y[i])+'\n'
        f.write(lineb)
        f.close()

def write_script():
    """ Write the script to run gravlens with the previous initial files."""
    nlens = 500

    file = 'init_gl/script_gl'
    f = open(file,'w')
    for i in range(nlens):
        line = '/Users/irenebalmes/Tout/gravlens/gl1.99o/gravlens lens_'+ \
            str(i)+'.in > lens_'+str(i)+'.out\n'
        f.write(line)
