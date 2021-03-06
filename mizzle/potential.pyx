"""Cython module with objective function:

.. math:: 
    D(e^{-2\sigma(r-r_0)}-2e^{-\sigma(r-r_0)})+
    \\frac{s(-2De^{-2\\sigma r_0}(e^{\\sigma r_0}-1))}{3}e^{-3\sigma(r-r_0)}+\\
    \sum_{i,j}\\frac{A}{r_{ij}^4}+\sum_{i,k}\\frac{B}{r_{ik}^4}

Notes
-----
Contains the objective function which is minimized during optimization

"""
import cython
import numpy as np
cimport numpy as np
from radish import Topologizer
from libc.math cimport exp 

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t
@cython.cdivision(True)
@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)

def potential(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours,
                np.ndarray[np.float64_t, ndim=1] boxVectors,
                np.ndarray[np.float64_t, ndim=1] bondlengths,
                np.ndarray[ITYPE_t, ndim=1] numSpecies):

    """Potential function

    Parameters
    ----------
    solvateCoords : numpy ndarray
        Independant variables of the potential function
    centers : ndarray
        Coordinates for centers which binds solvate
    topol : Topologizer instance
        Topologizer instance of system
    centerNeighbours : ndarray(float)
        coordinates of neighbours to solvate
    centerNumNeighbours : array(int)
        number of neighbours to each solvate
    boxVectors : numpy array
        Contains the box vectors

    Returns
    -------
    sumPot
        Value of the potential at given solvate coordinates

    """
    cdef int i = 0
    cdef int u = 0
    cdef int r = 0
    cdef int k = 0
    cdef int d = 0
    cdef double sumPot = 0

    #Does not drop duplicates(which is a good thing)
    cdef np.ndarray[np.float32_t, ndim=2] centersXYZ =\
                                 topol.trj.xyz[0][centers]*10
    cdef np.ndarray[dtype=double, ndim=1] tempNeighbour = np.empty([0], dtype=np.float64)
    cdef double distance = 0
    cdef int solvateLen = len(solvateCoords)/3
    cdef float eqDist = 2.2
    cdef float A = 20
    cdef float B = 5
    cdef float D = 40
    cdef float a = 1
    cdef float s = 0.01
    cdef float p = 4

    for r in range(solvateLen):
        # Harmonic potential
        distance = ((solvateCoords[r*3] - centersXYZ[i][0])**2 +\
                    (solvateCoords[r*3+1] - centersXYZ[i][1])**2 +\
                    (solvateCoords[r*3+2] - centersXYZ[i][2])**2)**(1./2)
        #sumPot += D*(exp(-2*a*(distance-eqDist))-2*exp(-a*(distance-eqDist)))+\
        #        s*(-2/3*D*exp(-2*a*eqDist)*(exp(a*eqDist)-1))*exp(-3*a*(distance-eqDist))
        #if(r > numSpecies[0]):
        #    sumPot += D*(exp(-2*a*(distance-bondlengths[1]))-2*exp(-a*(distance-bondlengths[1])))+\
        #            s*(-2/3*D*exp(-2*a*bondlengths[1])*(exp(a*bondlengths[1])-1))*exp(-3*a*(distance-bondlengths[1]))
        #else:
        #    sumPot += D*(exp(-2*a*(distance-bondlengths[0]))-2*exp(-a*(distance-bondlengths[0])))+\
        #            s*(-2/3*D*exp(-2*a*bondlengths[0])*(exp(a*bondlengths[0])-1))*exp(-3*a*(distance-bondlengths[0]))
        sumPot += D*(exp(-2*a*(distance-bondlengths[r]))-2*exp(-a*(distance-bondlengths[r])))+\
                s*(-2/3*D*exp(-2*a*bondlengths[r])*(exp(a*bondlengths[r])-1))*exp(-3*a*(distance-bondlengths[r]))
        #sumPot += 5*(distance - 2.2)**2

        # Solvate-neighbours pair-potential
        for d in range(centerNumNeighbours[r]):
            tempNeighbour = centerNeighbours[k]

            distance = ((solvateCoords[r*3] - tempNeighbour[0])**2 +\
                        (solvateCoords[r*3+1] - tempNeighbour[1])**2 +\
                        (solvateCoords[r*3+2] -\
                         tempNeighbour[2])**2)**(1./2)

            #sumPot += 1/(distance**p)
            sumPot += A/(distance**p)
            k += 1

        # Solvate - solvate pair-potential
        for u in range(solvateLen):
            if(i != u): #Don't include solvate[i] - solvate[i]
                tempNeighbour[0] = solvateCoords[u*3]
                tempNeighbour[1] = solvateCoords[u*3+1]
                tempNeighbour[2] = solvateCoords[u*3+2]
                if(boxVectors[0] > 0):
                    if(tempNeighbour[0] - solvateCoords[r*3] > boxVectors[0]/2):
                        tempNeighbour[0] = tempNeighbour[0] - boxVectors[0]
                    elif(tempNeighbour[0] - solvateCoords[r*3] < -boxVectors[0]/2):
                        tempNeighbour[0] = tempNeighbour[0] + boxVectors[0]
                    if(tempNeighbour[1] - solvateCoords[r*3 + 1] > boxVectors[1]/2):
                        tempNeighbour[1] = tempNeighbour[1] - boxVectors[1]
                    elif(tempNeighbour[1] - solvateCoords[r*3 + 1] < -boxVectors[1]/2):
                        tempNeighbour[1] = tempNeighbour[1] + boxVectors[1]
                    if(tempNeighbour[2] - solvateCoords[r*3 + 2] > boxVectors[2]/2):
                        tempNeighbour[2] = tempNeighbour[2] - boxVectors[2]
                    elif(tempNeighbour[2] - solvateCoords[r*3 + 2] < -boxVectors[2]/2):
                        tempNeighbour[2] = tempNeighbour[2] + boxVectors[2]

                distance = ((solvateCoords[r*3] - tempNeighbour[0])**2 +\
                            (solvateCoords[r*3+1] - tempNeighbour[1])**2 +\
                            (solvateCoords[r*3+2] -\
                            tempNeighbour[2])**2)**(1./2)

                #if(distance != 0.0):
                #sumPot += 1/(distance**p)
                sumPot += B/(distance**p)
        i += 1
    return sumPot

@cython.cdivision(True)
@cython.boundscheck(False)
def potential_jac(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours,
                np.ndarray[np.float64_t, ndim=1] boxVectors,
                np.ndarray[np.float64_t, ndim=1] bondlengths,
                np.ndarray[ITYPE_t, ndim=1] numSpecies):

    """Jacobian of potential

    Parameters
    ----------
    Same as potential.potential_c
        
    Returns
    -------
    jac
        Jacobian of potential.potential_c

    """
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef int d = 0
    cdef np.ndarray[np.float32_t, ndim=2] centersXYZ
    cdef double denom = 0
    cdef int solvateLen = len(solvateCoords)/3
    cdef np.ndarray[np.float64_t, ndim=1] jac = np.empty(solvateLen*3, dtype=np.float64)
    cdef np.ndarray[dtype=double, ndim=1] tempNeighbour = np.empty([0], dtype=np.float64)
    centersXYZ = topol.trj.xyz[0][centers]*10

    cdef float eqDist = 2.2
    cdef float A = 20
    cdef float B = 5
    cdef float D = 40
    cdef float a = 1
    cdef float s = 0.01
    cdef float p = 4

    for i in range(solvateLen):
        # Harmonic potential
        denom = ((solvateCoords[3*i] - centersXYZ[i][0])**2 +\
                 (solvateCoords[3*i+1] - centersXYZ[i][1])**2 +\
                 (solvateCoords[3*i+2] - centersXYZ[i][2])**2)**(1./2)

        distance = ((solvateCoords[3*i] - centersXYZ[i][0])**2 +\
                    (solvateCoords[3*i+1] - centersXYZ[i][1])**2 +\
                    (solvateCoords[3*i+2] - centersXYZ[i][2])**2)**(1./2)
        #if(i > numSpecies[0]):
        #    jac[3*i] = (2*a*(exp(bondlengths[1]*a)-1)*D*s*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-3*a*(distance - bondlengths[1])-2*bondlengths[1]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-a*(distance - bondlengths[1])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-2*a*(distance - bondlengths[1])))/distance)

        #    jac[3*i+1] = (2*a*(exp(bondlengths[1]*a)-1)*D*s*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-3*a*(distance - bondlengths[1])-2*bondlengths[1]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-a*(distance - bondlengths[1])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-2*a*(distance - bondlengths[1])))/distance)

        #    jac[3*i+2] = (2*a*(exp(bondlengths[1]*a)-1)*D*s*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-3*a*(distance - bondlengths[1])-2*bondlengths[1]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-a*(distance - bondlengths[1])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-2*a*(distance - bondlengths[1])))/distance)
        #else:
        #    jac[3*i] = (2*a*(exp(bondlengths[0]*a)-1)*D*s*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-3*a*(distance - bondlengths[0])-2*bondlengths[0]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-a*(distance - bondlengths[0])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-2*a*(distance - bondlengths[0])))/distance)

        #    jac[3*i+1] = (2*a*(exp(bondlengths[0]*a)-1)*D*s*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-3*a*(distance - bondlengths[0])-2*bondlengths[0]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-a*(distance - bondlengths[0])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-2*a*(distance - bondlengths[0])))/distance)

        #    jac[3*i+2] = (2*a*(exp(bondlengths[0]*a)-1)*D*s*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-3*a*(distance - bondlengths[0])-2*bondlengths[0]*a))/distance +\
        #    D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-a*(distance - bondlengths[0])))/distance) -\
        #    D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-2*a*(distance - bondlengths[0])))/distance)  


        jac[3*i] = (2*a*(exp(bondlengths[i]*a)-1)*D*s*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-3*a*(distance - bondlengths[i])-2*bondlengths[i]*a))/distance +\
        D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-a*(distance - bondlengths[i])))/distance) -\
        D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-2*a*(distance - bondlengths[i])))/distance)

        jac[3*i+1] = (2*a*(exp(bondlengths[i]*a)-1)*D*s*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-3*a*(distance - bondlengths[i])-2*bondlengths[i]*a))/distance +\
        D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-a*(distance - bondlengths[i])))/distance) -\
        D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-2*a*(distance - bondlengths[i])))/distance)

        jac[3*i+2] = (2*a*(exp(bondlengths[i]*a)-1)*D*s*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-3*a*(distance - bondlengths[i])-2*bondlengths[i]*a))/distance +\
        D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-a*(distance - bondlengths[i])))/distance) -\
        D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-2*a*(distance - bondlengths[i])))/distance)  



        #jac[3*i] = (2*a*(exp(eqDist*a)-1)*D*s*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-3*a*(distance - eqDist)-2*eqDist*a))/distance +\
        #D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-a*(distance - eqDist)))/distance) -\
        #D*((2*a*(solvateCoords[3*i] - centersXYZ[i][0])*exp(-2*a*(distance - eqDist)))/distance)

        #jac[3*i+1] = (2*a*(exp(eqDist*a)-1)*D*s*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-3*a*(distance - eqDist)-2*eqDist*a))/distance +\
        #D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-a*(distance - eqDist)))/distance) -\
        #D*((2*a*(solvateCoords[3*i+1] - centersXYZ[i][1])*exp(-2*a*(distance - eqDist)))/distance)

        #jac[3*i+2] = (2*a*(exp(eqDist*a)-1)*D*s*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-3*a*(distance - eqDist)-2*eqDist*a))/distance +\
        #D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-a*(distance - eqDist)))/distance) -\
        #D*((2*a*(solvateCoords[3*i+2] - centersXYZ[i][2])*exp(-2*a*(distance - eqDist)))/distance)

        #jac[3*i] = 10*(solvateCoords[3*i]-centersXYZ[i][0])*(denom-2.2)/denom
        #jac[3*i+1] = 10*(solvateCoords[3*i+1]-centersXYZ[i][1])*(denom-2.2)/denom
        #jac[3*i+2] = 10*(solvateCoords[3*i+2]-centersXYZ[i][2])*(denom-2.2)/denom       
        
        # Solvate - neighbours pair-potential
        for d in range(centerNumNeighbours[i]):
            tempNeighbour = centerNeighbours[k]

            denom = ((solvateCoords[3*i] - tempNeighbour[0])**2 +\
                     (solvateCoords[3*i+1] - tempNeighbour[1])**2 +\
                     (solvateCoords[3*i+2] - tempNeighbour[2])**2)**(p-1)

            jac[3*i] += -p*A*(solvateCoords[3*i] - tempNeighbour[0])/denom
            jac[3*i+1] += -p*A*(solvateCoords[3*i+1] - tempNeighbour[1])/denom
            jac[3*i+2] += -p*A*(solvateCoords[3*i+2] - tempNeighbour[2])/denom

            k += 1

        #Solvate - solvate pair-potential
        for j in range(solvateLen):
            if(i != j): #Don't include solvate[i] - solvate[i]
                tempNeighbour[0] = solvateCoords[j*3]
                tempNeighbour[1] = solvateCoords[j*3+1]
                tempNeighbour[2] = solvateCoords[j*3+2]
                if(boxVectors[0] > 0):
                    if(tempNeighbour[0] - solvateCoords[i*3] > boxVectors[0]/2):
                        tempNeighbour[0] = tempNeighbour[0] - boxVectors[0]
                    elif(tempNeighbour[0] - solvateCoords[i*3] < -boxVectors[0]/2):
                        tempNeighbour[0] = tempNeighbour[0] + boxVectors[0]
                    if(tempNeighbour[1] - solvateCoords[i*3 + 1] > boxVectors[1]/2):
                        tempNeighbour[1] = tempNeighbour[1] - boxVectors[1]
                    elif(tempNeighbour[1] - solvateCoords[i*3 + 1] < -boxVectors[1]/2):
                        tempNeighbour[1] = tempNeighbour[1] + boxVectors[1]
                    if(tempNeighbour[2] - solvateCoords[i*3 + 2] > boxVectors[2]/2):
                        tempNeighbour[2] = tempNeighbour[2] - boxVectors[2]
                    elif(tempNeighbour[2] - solvateCoords[i*3 + 2] < -boxVectors[2]/2):
                        tempNeighbour[2] = tempNeighbour[2] + boxVectors[2]

                denom = ((solvateCoords[3*i] - tempNeighbour[0])**2 +\
                         (solvateCoords[3*i+1] - tempNeighbour[1])**2 +\
                         (solvateCoords[3*i+2] - tempNeighbour[2])**2)**(p-1)

                jac[3*i] += (-p)*B*(solvateCoords[3*i] - tempNeighbour[0])/denom
                jac[3*i+1] += (-p)*B*(solvateCoords[3*i+1] - tempNeighbour[1])/denom
                jac[3*i+2] += (-p)*B*(solvateCoords[3*i+2] - tempNeighbour[2])/denom

                denom = ((tempNeighbour[0] - solvateCoords[3*i])**2 +\
                         (tempNeighbour[1] - solvateCoords[3*i+1])**2 +\
                         (tempNeighbour[2] - solvateCoords[3*i+2])**2)**(p-1)

                jac[3*i] += p*B*(tempNeighbour[0] - solvateCoords[3*i])/denom
                jac[3*i+1] += p*B*(tempNeighbour[1] - solvateCoords[3*i+1])/denom
                jac[3*i+2] += p*B*(tempNeighbour[2] - solvateCoords[3*i+2])/denom
    return jac      