import sys
import numpy as np
from scipy.linalg.lapack import dgges
from numpy.polynomial.legendre import leggauss

# Read in command line args from bash
if len(sys.argv) != 7:
    sys.exit('Error, needs nPart cubeMult yp NXNYMAX nq root')

#SLy4
T0 = -2488.91
T1 = 486.82
T2 = -546.39
T3 = 13777.0
X0 = 0.834
X1 = -0.344
X2 = -1.000
X3 = 1.354
ALPHA = 1./6.

hbar2m = 20.72125 #hbar^2/2m (197.3^2 MeV^2 fm^2)/2(931.5 MeV)
chargeSquared = 1.44 #natural units of proton/electron charge (197.3 MeV fm)/137 fine stucture constant
COUNTLIMIT = 20 #arbitratry hard stop if convergence not reached
PI = np.pi

#Command line arguments from shell script
nMult = int(sys.argv[1])
cubeMult = int(sys.argv[2])
yp = int(sys.argv[3]) #int value for labeling
NXNYMAX = int(sys.argv[4]) #int value for scanning degeneracies in x and y direction
nq = int(sys.argv[5]) #number of oscillations
root = sys.argv[6].strip() #scratch path to be used for file output, excluding final E/A which is stdout

nParticles = nMult*cubeMult**3
nPoints = 600 #or points along z axis
rho0_tot = 0.160 #Need to adjust, if saturation density, for yp
tau0_tot = (3./5.)*(3*PI**2)**(2./3.)*rho0_tot**(5./3.)
boxLength = (nParticles/rho0_tot)**(1./3.)
dz = boxLength/nPoints
I = np.identity(nPoints)
zVec = np.linspace(0., boxLength, nPoints, endpoint=False) # up to but not including L, which is just 0


EF = hbar2m*((3*PI*PI*rho0_tot)**(2./3.))
vq_index = 0.125
vq = vq_index*EF #v_q from Sam's paper is actually 2v_q, his code also labels files as 2v_q
phiTol = 1.e-6 #convergence

#neutron & protons
ypCoef = float(yp/100.0)
protPart = int(ypCoef*nParticles)

#For guagleg integration over fixed r
x, w = leggauss(nPoints)
x = 0.5*(x + 1)*boxLength
w = 0.5*w*boxLength

#external potential, if nq is zero, return 0 array
def extPot():
    if nq == 0:
        return np.zeros(nPoints)
    else:
        q = (2.*PI*nq)/boxLength
        return 2.*vq*np.cos(q*zVec)

#Dummy function required in dgges
def dselect(alphar, alphai, beta):
    return 0

extPotVec = extPot() #intitilize external potential from parameters module, will not change during program
