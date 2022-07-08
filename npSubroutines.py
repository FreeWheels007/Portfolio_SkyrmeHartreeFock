from parameters import *
from states import State

def firder(arr):
	arrder = np.zeros(nPoints)
	
	for i in range(nPoints):
		arrder[i] = (arr[(i+3)%nPoints] - 9*arr[(i+2)%nPoints] + 45*arr[(i+1)%nPoints] - 45*arr[(i-1)%nPoints] + 9*arr[(i-2)%nPoints] - arr[(i-3)%nPoints])/(60.*dz)
	
	return arrder

def secder(arr):
	arr2der = np.zeros(nPoints)
	
	for i in range(nPoints):
		arr2der[i] = (arr[(i+3)%nPoints] - 13.5*arr[(i+2)%nPoints] + 135*arr[(i+1)%nPoints] - 245*arr[i] + 135*arr[(i-1)%nPoints] - 13.5*arr[(i-2)%nPoints] + arr[(i-3)%nPoints])/(90.*dz*dz)
	
	return arr2der
	
def firderC(arr):
	arrder = np.zeros(nPoints) + 1j*np.zeros(nPoints)
	
	for i in range(nPoints):
		arrder[i] = (arr[(i+3)%nPoints] - 9*arr[(i+2)%nPoints] + 45*arr[(i+1)%nPoints] - 45*arr[(i-1)%nPoints] + 9*arr[(i-2)%nPoints] - arr[(i-3)%nPoints])/(60.*dz)
	
	return arrder

	
def absvec(x):
	return np.sqrt(x*x)

#density function
def rhoFunc(eigPairList, nParticles_q):
	rhoVec = np.zeros(nPoints)
	
	#sum over all particle
	pCount = 0
	for i in range(len(eigPairList)):
		#temp variables, for ease of reading final sum loop
		#phiVec = eigPairList[i].eigVec
		phiVec = eigPairList[i].compVec
		deg = eigPairList[i].deg
		
		if ((pCount + deg) >= nParticles_q): #ensure last degeneracy does not go over particle limit
			rhoVec += (nParticles_q - pCount)*(abs(phiVec)**2)
			break
			
		rhoVec += deg*(abs(phiVec)**2)
		pCount += deg
	
	#print(rhoVec)
	
	return rhoVec/(boxLength*boxLength)

#kinetic density calculator
def tauFunc(eigPairList, nParticles_q):
	tauVec = np.zeros(nPoints)
	
	#sum over all particles
	pCount = 0
	for i in range(len(eigPairList)):
		#temp variables, for ease of reading final sum loop
		#phiVec = eigPairList[i].eigVec
		phiVec = eigPairList[i].compVec
		dphiVec = firderC(phiVec)
		nx, ny = eigPairList[i].nx, eigPairList[i].ny
		deg = eigPairList[i].deg
		
		if ((pCount + deg) >= nParticles_q): #ensure last degeneracy does not go over particle limit
			tauVec += (nParticles_q - pCount)*((boxLength*boxLength*abs(dphiVec)**2)/(4.*PI*PI) + (nx*nx +ny*ny)*abs(phiVec)**2)
			break
		
		tauVec += deg*((boxLength*boxLength*abs(dphiVec)**2)/(4.*PI*PI) + (nx*nx +ny*ny)*abs(phiVec)**2)
		pCount += deg
	
	return (4.*PI*PI/boxLength**4)*tauVec

#set complex vectors from degenerate eigenvectors
def createComplexVec(eigPairList):
	listLen = len(eigPairList)
	i = 0
	
	while i < listLen:
		if (i == listLen-1) or not(np.isclose(eigPairList[i].eigVal, eigPairList[i+1].eigVal) and eigPairList[i].nx == eigPairList[i+1].nx and eigPairList[i].ny == eigPairList[i+1].ny):
			eigPairList[i].compVec = eigPairList[i].eigVec + 1j*np.zeros(nPoints)
			i = i + 1			
		else:
			eigPairList[i].compVec = (eigPairList[i].eigVec + 1j*eigPairList[i+1].eigVec)/np.sqrt(2.)
			eigPairList[i+1].compVec = (eigPairList[i].eigVec - 1j*eigPairList[i+1].eigVec)/np.sqrt(2.)
			i = i + 2


#effective mass
def effmass(rho, rho_q):
	return hbar2m + 0.125*(T1*(2.+X1)+T2*(2.+X2))*rho - 0.125*(T1*(1.+2*X1)-T2*(1.+2*X2))*rho_q

#skyrme potential
def skyrmePot(rho, tau, rho_q, tau_q, nx, ny):
	term1 = 0.5*T0*((2.+X0)*rho - (1.+2*X0)*rho_q)
	
	rho_qnp = rho - rho_q #opposite q
	
	term2 = (1./24.)*T3*((2.+X3)*(2.+ALPHA)*rho**(ALPHA+1.) - (2*X3+1.)*(2*rho_q*rho**ALPHA + ALPHA*(rho_q**2 + rho_qnp**2)*rho**(ALPHA-1.)))
	
	term3 = 0.125*((T1*(2.+X1) + T2*(2.+X2))*tau + (T2*(1.+2*X2) - T1*(1.+2*X1))*tau_q)
	
	term4 = 0.1875*((T2*(2.+X2) - 3*T1*(2.+X1))*secder(rho) + (3*T1*(1.+2*X1) + T2*(1.+2*X2))*secder(rho_q))
	
	return term1 + term2 + term3 + term4


#A Function
def AFunc(rho, rho_q):
	return -effmass(rho, rho_q)

#B Function
def BFunc(rho, rho_q):
	return -firder(effmass(rho, rho_q))

#C Function
def CFunc(rho, tau, rho_q, tau_q, nx, ny, vPotz, q):
	return skyrmePot(rho, tau, rho_q, tau_q, nx, ny) + extPotVec + q*vPotz + effmass(rho, rho_q)*(nx*nx + ny*ny)*((4.*PI*PI)/(boxLength*boxLength))

	
#Abar for matrix
def Abar(rho, rho_q):
	return AFunc(rho, rho_q)/(90.*dz*dz)

#Bbar for matrix
def Bbar(rho, rho_q):
	return BFunc(rho, rho_q)/(60.*dz)


#M matrix creator and populator
def mMatrixConstructor(rho, tau, rho_q, tau_q, nx, ny, vPotz, q):
	AVec = Abar(rho, rho_q)
	BVec = Bbar(rho, rho_q)
	CVec = CFunc(rho, tau, rho_q, tau_q, nx, ny, vPotz, q)
	
	mMatrix = np.zeros((nPoints, nPoints))
	
	for j in range(nPoints):
		if (j == 0):
			#first line
			mMatrix[j, 0] = CVec[j] - 245*AVec[j]
			mMatrix[j, 1] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, 2] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, 3] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-3] = AVec[j] - BVec[j]
			mMatrix[j, nPoints-2] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, nPoints-1] = 135*AVec[j] - 45*BVec[j]
		elif (j == 1):
			#second
			mMatrix[j, 0] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, 1] = CVec[j] - 245*AVec[j]
			mMatrix[j, 2] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, 3] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, 4] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-2] = AVec[j] - BVec[j]
			mMatrix[j, nPoints-1] = -(13.5*AVec[j] - 9*BVec[j])
		elif (j == 2):
			#third
			mMatrix[j, 0] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, 1] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, 2] = CVec[j] - 245*AVec[j]
			mMatrix[j, 3] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, 4] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, 5] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-1] = AVec[j] - BVec[j]
		elif (j == nPoints-3):
			#third last
			mMatrix[j, 0] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-6] = AVec[j] - BVec[j]
			mMatrix[j, nPoints-5] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, nPoints-4] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, nPoints-3] = CVec[j] - 245*AVec[j]
			mMatrix[j, nPoints-2] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, nPoints-1] = -(13.5*AVec[j] + 9*BVec[j])
		elif (j == nPoints-2):
			#second last
			mMatrix[j, 0] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, 1] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-5] = AVec[j] - BVec[j]
			mMatrix[j, nPoints-4] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, nPoints-3] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, nPoints-2] = CVec[j] - 245*AVec[j]
			mMatrix[j, nPoints-1] = 135*AVec[j] + 45*BVec[j]
		elif (j == nPoints-1):
			#last
			mMatrix[j, 0] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, 1] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, 2] = AVec[j] + BVec[j]
			mMatrix[j, nPoints-4] = AVec[j] - BVec[j]
			mMatrix[j, nPoints-3] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, nPoints-2] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, nPoints-1] = CVec[j] - 245*AVec[j]
		else:
			#rest in between
			mMatrix[j, j-3] = AVec[j] - BVec[j]
			mMatrix[j, j-2] = -(13.5*AVec[j] - 9*BVec[j])
			mMatrix[j, j-1] = 135*AVec[j] - 45*BVec[j]
			mMatrix[j, j] = CVec[j] - 245*AVec[j]
			mMatrix[j, j+1] = 135*AVec[j] + 45*BVec[j]
			mMatrix[j, j+2] = -(13.5*AVec[j] + 9*BVec[j])
			mMatrix[j, j+3] = AVec[j] + BVec[j]
	
	return mMatrix

def degenFactor(nx, ny):
	
	#MY Degeneracies
	if (nx == 0) and (ny == 0):
		deg = 1
	elif (nx == 0) or (ny == 0):
		deg = 4
	elif (nx == ny):
		deg = 4
	else:
		deg = 8
	
	return deg*2

def hamiltonian(rho_tot, tau_tot, rho_n, tau_n, rho_p, tau_p, vPotz):
	term1 = hbar2m*tau_tot
	term2 = (1./4.)*T0*((2.+X0)*rho_tot**2 - (2*X0+1.)*(rho_n**2 + rho_p**2))
	term3 = (1./24.)*T3*((2.+X3)*rho_tot**2 - (2*X3+1.)*(rho_n**2 + rho_p**2))*rho_tot**ALPHA
	term4 = (1./8.)*((T1*(2.+X1) + T2*(2.+X2))*rho_tot*tau_tot + (T2*(2*X2+1) - T1*(2*X1+1))*(rho_n*tau_n + rho_p*tau_p))
	term5 = (1./32.)*((3*T1*(2.+X1) - T2*(2.+X2))*(firder(rho_tot))**2 - (3*T1*(2*X1+1.) + T2*(2*X2+1.))*((firder(rho_n))**2 + (firder(rho_p))**2))
	term6 = extPotVec*rho_tot
	term7 = vPotz*rho_p
	
	return term1 + term2 + term3 + term4 + term5 + term6 + term7 
