################################################################################
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rand
import lens as l

def draw_lensparam(n,thetaE):
    """Draw ellipticity, shear, shear angle and source positions.

    Draw n sets of lens parameters from the given distributions (see Oguri &
    Marshall 2010)
    Arguments:
    n -- number of sets to draw
    thetaE -- typical Einstein radius of the lens
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

    # Note that lognormal takes as parameters the mean and sigma of the 
    # _underlying_ normal distribution.
    # This implies the following relation between mu and sigma (to feed to
    # rand.lognormal) and M and S (from the Oguri & Marshall article)
    M = 0.05
    S = 0.2
    sigma = np.sqrt(np.log(S/M**2+1))
    mu = np.log(M)-sigma**2/2
    gamma = rand.lognormal(mu,sigma,n)
    for i in range(n):
        while (gamma[i]>1.0):
            gamma[i] = rand.lognormal(mu,sigma)

    theta_g = rand.uniform(0,180,n)

    x = rand.uniform(-thetaE,thetaE,n)
    y = rand.uniform(-thetaE,thetaE,n)

    return e,gamma,theta_g,x,y

def write_initfile(v,zl,zs):
    """Write the initial files for gravlens for a given Einstein radius.

    Arguments:
    v -- velocity dispersion of the lens, in km/s
    zl -- lens redshift
    zs -- source redshift
    Outputs:
    tonnes of files
    """
    nlens = 10#000#0
    c = 299792.458 # km/s
    arcsectorad = 4.85e-6
    thetaE = 4*np.pi*(v/c)**2*l.distance(zs,zl)/l.distance(zs)
    # the above formula gives thetaE in radian
    # converting to arcseconds
    thetaE = thetaE/arcsectorad
    e,gamma,theta_g,x,y = draw_lensparam(nlens,thetaE)
    print thetaE

    for i in range(nlens):
        file = 'init_gl/lens_'+str(i)+'.in'
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
        line = 'lens pert '+str(zs)+' 0 0 '+str(gamma[i])+' '+ \
               str(theta_g[i])+' 0 0\n'
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

def write_script():
    """Write the script to run gravlens with the initial files."""
    nlens = 10#000#0

    file = 'script_gl'
    f = open(file,'w')
    for i in range(nlens):
        line = './glafic init_gl/lens_'+str(i)+'.in > init_gl/lens_'\
               +str(i)+'.out\n'
        f.write(line)

def analyse_output(M,zs):
    """Using the files output by glafic, compute the lensing cross-section.

    Arguments:
    M -- magnitude of the source.
    zs -- redshift of the source.
    nlens -- number of lenses to take into account
    Outputs:
    sigma -- biased lensing cross-section computed from glafic output files.
    """
    nlens = 10#000#0

    # array containing the magnification, and number of multiply imaged lenses
    mu = []
    mutmp = []
    nmult = 0
    for i in range(nlens):
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
            if mutmp[0]/mutmp[1]>0.1:
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

    mu = np.array(mu)
    mag = M+2.5*np.log10(mu)

    # absolute i-band magnitude luminosity function
    lumfunction = l.quasarluminosity(mag,zs)
    lumfunction_nomag = l.quasarluminosity(M,zs)

    sigma = 0
    sigtmp = []
    var = []
    for i in range(nmult):
        lumratio = lumfunction[i]/lumfunction_nomag
        var.append(lumratio/mu[i])
        sigma = sigma+lumratio/mu[i]
        sigtmp.append(sigma/(i+1))
    if nmult == 0:
        sigma = 0
    else:
        sigma = sigma/nmult

    return sigma
