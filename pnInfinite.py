from parameters import *
from states import State
from npSubroutines import *
import sys
from fcoul2py import coulomb_potential

if __name__ == '__main__':
        isConverged = False
        count = 0
        rho_tot = np.full((nPoints), rho0_tot)
        tau_tot = np.full((nPoints), tau0_tot)
        rho_n = np.full((nPoints), (1. - ypCoef)*rho0_tot)
        tau_n = np.full((nPoints), (1. - ypCoef)*tau0_tot)
        rho_p = np.full((nPoints), ypCoef*rho0_tot)
        tau_p = np.full((nPoints), ypCoef*tau0_tot)
        prevRho_tot = np.zeros(nPoints)
        prevTau_tot = np.zeros(nPoints)
        prevRho_n = np.zeros(nPoints)
        prevTau_n = np.zeros(nPoints)
        prevRho_p = np.zeros(nPoints)
        prevTau_p = np.zeros(nPoints)
        # with rho_p starting as flat, assume only exchange term for initial condition
        vPotz = -0.75*chargeSquared*((3./PI)**(1./3.))*(rho_p**(1./3.))
        
        
        energy = hamiltonian(rho_tot, tau_tot, rho_n, tau_n, rho_p, tau_p, vPotz)
        prevTotalEnergy = (boxLength*boxLength)*State.gaulegint(State, energy)
        
        while not isConverged:
                nx, ny = 0, 0
                eigPairList = []
                nxnyFinished = False
                
                # Neutron Matrix and Eigenfunction
                
                while not nxnyFinished:
                        mMatrix = mMatrixConstructor(rho_tot, tau_tot, rho_n, tau_n, nx, ny, vPotz, 0)
                        
                        outDgges = dgges(dselect, mMatrix, I, jobvsl=0, jobvsr=1, sort_t=0, ldvsr=nPoints, lwork=8*nPoints+16)
                        a, b, sdim, alphar, alphai, beta, vsl, vsr, work, info = outDgges
                        #NOTE, vsr is phi
                        
                        if (info != 0):
                                sys.exit('Error dgges: ', info)
                        
                        #create finished eigenvalues from alpha real and beta and add them to list
                        eigValues = alphar/beta
                        for i in range(len(eigValues)):
                                eigPairList.append(State(np.around(eigValues[i],7), vsr.T[i], nx, ny, degenFactor(nx, ny), i))
                        
                        
                        #my deg
                        if (nx == NXNYMAX) and (ny == NXNYMAX):
                                nxnyFinished = True
                        elif (nx == ny):
                                nx = 0
                                ny += 1
                        else:
                                nx += 1
                        
                
                # sort eigPairs, from least to greatest, by eigenvalue (eigenenergy), then by nx ny 
                eigPairList.sort(key=lambda eigPair: (eigPair.eigVal, eigPair.nx, eigPair.ny))
                
                
                # make complex vectors from degenerate eigenvectors
                createComplexVec(eigPairList)
                
                # ensure the explored nx ny space produced enough eigenstates
                if (len(eigPairList) < (nParticles-protPart)):
                        sys.exit('Error, nx/ny searchspace too small for NEUTRONS')
                
                # with updated eigenvectors (from phi) and degeneracies, update rho and tau
                rho_n = rhoFunc(eigPairList, nParticles-protPart)
                tau_n = tauFunc(eigPairList, nParticles-protPart)
                
                
                # END of Neutron Calc clear list
                eigPairList = []
                nxnyFinished = False
                nx, ny = 0, 0
                
                # Proton Matrix and Eigenfunction
                
                while not nxnyFinished:
                        mMatrix = mMatrixConstructor(rho_tot, tau_tot, rho_p, tau_p, nx, ny, vPotz, 1)
                        
                        outDgges = dgges(dselect, mMatrix, I, jobvsl=0, jobvsr=1, sort_t=0, ldvsr=nPoints, lwork=8*nPoints+16)
                        a, b, sdim, alphar, alphai, beta, vsl, vsr, work, info = outDgges
                        #NOTE, vsr is phi
                        
                        if (info != 0):
                                sys.exit('Error dgges: ', info)
                        
                        # create finished eigenvalues from alpha real and beta and add them to list
                        eigValues = alphar/beta
                        for i in range(len(eigValues)):
                                eigPairList.append(State(np.around(eigValues[i],7), vsr.T[i], nx, ny, degenFactor(nx, ny), i))
                        
                        
                        # my deg
                        if (nx == NXNYMAX) and (ny == NXNYMAX):
                                nxnyFinished = True
                        elif (nx == ny):
                                nx = 0
                                ny += 1
                        else:
                                nx += 1
                        
                
                # code sort eigPairs
                eigPairList.sort(key=lambda eigPair: (eigPair.eigVal, eigPair.nx, eigPair.ny))
                
                # make complex vectors from degenerate eigenvectors
                createComplexVec(eigPairList)
                
                      
                if (len(eigPairList) < protPart):
                        sys.exit('Error, nx/ny searchspace too small for PROTONS')
                        
                # with updated eigenvectors (from phi) and degeneracies, update rho and tau
                rho_p = rhoFunc(eigPairList, protPart)
                tau_p = tauFunc(eigPairList, protPart)

                # calculate new coulomb potential
                vPotz = coulomb_potential(zvec=zVec, rho=rho_p, d2rho=secder(rho_p), chargefactor=chargeSquared, boxlength=boxLength, dz=dz, npoints=nPoints)
                
                # Combine debsities
                rho_tot = rho_n + rho_p
                tau_tot = tau_n + tau_p
                
                # END of proton calc
                
                # Calculate system energy from updated hamiltonian
                energy = hamiltonian(rho_tot, tau_tot, rho_n, tau_n, rho_p, tau_p, vPotz)
                totalEnergy = (boxLength*boxLength)*State.gaulegint(State, energy)
                
                #exit loop if convergence is reached
                if abs(totalEnergy - prevTotalEnergy) < phiTol*prevTotalEnergy:
                        isConverged = True
                        #print('reached convergence')
                #or exit if minimum is passed
                elif (count > 0) and (totalEnergy > prevTotalEnergy):
                        isConverged = True
                        totalEnergy = prevTotalEnergy
                        rho_tot = np.copy(prevRho_tot)
                        tau_tot = np.copy(prevTau_tot)
                        rho_n = np.copy(prevRho_n)
                        tau_n = np.copy(prevTau_n)
                        rho_p = np.copy(prevRho_p)
                        tau_p = np.copy(prevTau_p)
                        energy = np.copy(prevEnergy)
                        #print('reached minima')
                elif (count >= COUNTLIMIT):
                        print('no converge')
                        break
                else:
                        prevRho_tot = np.copy(rho_tot)
                        prevTau_tot = np.copy(tau_tot)
                        prevRho_n = np.copy(rho_n)
                        prevTau_n = np.copy(tau_n)
                        prevRho_p = np.copy(rho_p)
                        prevTau_p = np.copy(tau_p)
                        prevEnergy = np.copy(energy)
                        prevTotalEnergy = totalEnergy
                        count += 1
        
        # Output file naming code
        rhoFileStr = '{}/rhoAtZvq{:03.2f}_n{}_acc{}_nxny{}_part{}.dat'.format(root, vq_index, nq, nPoints, NXNYMAX, nParticles)
        rhonFileStr = '{}/rhonAtZvq{:03.2f}_n{}_acc{}_nxny{}_part{}.dat'.format(root, vq_index, nq, nPoints, NXNYMAX, nParticles)
        rhopFileStr = '{}/rhopAtZvq{:03.2f}_n{}_acc{}_nxny{}_part{}.dat'.format(root, vq_index, nq, nPoints, NXNYMAX, nParticles)
        energyFileStr = '{}/energyAtZvq{:03.2f}_n{}_acc{}_nxny{}_part{}.dat'.format(root, vq_index, nq, nPoints, NXNYMAX, nParticles)
        
        
        rhoFile, rhonFile, rhopFile, energyFile = open(rhoFileStr, 'w'), open(rhonFileStr, 'w'), open(rhopFileStr, 'w'), open(energyFileStr, 'w')
        try:
            np.savetxt(rhoFile, np.column_stack((zVec, rho_tot)), delimiter=' ', fmt='%s')
            np.savetxt(rhonFile, np.column_stack((zVec, rho_n)), delimiter=' ', fmt='%s')
            np.savetxt(rhopFile, np.column_stack((zVec, rho_p)), delimiter=' ', fmt='%s')
            np.savetxt(energyFile, np.column_stack((zVec, energy)), delimiter=' ', fmt='%s')
        finally:
            rhoFile.close()
            rhonFile.close()
            rhopFile.close()
            energyFile.close()

        # Let slurm define output db location
        print('{},{},{},{},{},{}\n'.format(nParticles, nMult, cubeMult, nPoints, totalEnergy/nParticles, isConverged))

