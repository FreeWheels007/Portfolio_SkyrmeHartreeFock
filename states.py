import numpy as np
from math import floor
from parameters import *

# State object to hold individual eigenstates calculated by DGGES
# with corrosponding eigenvalue, eigenvector, nx ny component and resulting degeneracy
# vector is normalized with creation of object
class State:
	def __init__(self, eigVal, eigVec, nx, ny, degeneracy, eigValIndex):
		self.eigVal = eigVal
		self.nx = nx
		self.ny = ny
		self.deg = degeneracy
		self.eigValIndex = eigValIndex
		
		norm = np.sqrt(self.gaulegint(eigVec**2))
		
		#normalize eigenvector if normalization const is not 0
		if not np.isclose(norm, 0.):
			self.eigVec = eigVec/norm
		else:
			self.eigVec = eigVec
		
		self.compVec = 0 + 0j

	
	
	def gaulegint(self, vec):
		fsum = 0.
		for i in range(len(vec)):
			interlow = floor(x[i]/dz)
			interhigh = interlow + 1
			interweight = x[i]/dz - floor(x[i]/dz)
			
			if (interlow == (nPoints-1)):
				fsum += w[i]*(vec[nPoints-1]*(1. - interweight) + vec[0]*interweight)
			else:
				fsum += w[i]*(vec[interlow]*(1. - interweight) + vec[interhigh]*interweight)
			
			
		return fsum
