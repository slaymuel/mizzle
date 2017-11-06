"""Example NumPy style docstrings.

This Cython module contains the objective function and the Jacobian


Notes
-----
The objective function is minimized during optimization

"""

#cython: cdivision=True
#cython: boundscheck=False, wraparound=False, nonecheck=False
import cython
import numpy as np
cimport numpy as np
from radish import Topologizer

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t

#@cython.boundscheck(False)

def potential_c(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours,
                np.ndarray[np.float64_t, ndim=1] boxVectors):

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

#unsigned int
    cdef int i = 0
    cdef int u = 0
    cdef int r = 0
    cdef int k = 0
    cdef int d = 0
    cdef double sumPot = 0

    #Does not drop duplicates(which is a good thing)
    cdef np.ndarray[np.float32_t, ndim=2] centersXYZ =\
                                 topol.trj.xyz[0][centers]*10
    cdef tempNeighbour = np.array([0, 0, 0], dtype=float)
    cdef double distance = 0
    cdef int solvateLen = len(solvateCoords)/3

    for r in range(solvateLen):
        # Harmonic potential
        distance = ((solvateCoords[r*3] - centersXYZ[i][0])**2 +\
                    (solvateCoords[r*3+1] - centersXYZ[i][1])**2 +\
                    (solvateCoords[r*3+2] - centersXYZ[i][2])**2)**(1./2)
        sumPot += 5*(distance - 2.2)**2

        # Solvate-neighbours pair-potential
        for d in range(centerNumNeighbours[r]):
            tempNeighbour = centerNeighbours[k]

            distance = ((solvateCoords[r*3] - tempNeighbour[0])**2 +\
                        (solvateCoords[r*3+1] - tempNeighbour[1])**2 +\
                        (solvateCoords[r*3+2] -\
                         tempNeighbour[2])**2)**(1./2)

            sumPot += 1/(distance**4)
            k += 1

        # Solvate - solvate pair-potential
        for u in range(solvateLen):
            if(i != u): #Don't include solvate[i] - solvate[i]
                tempNeighbour[0] = solvateCoords[u*3]
                tempNeighbour[1] = solvateCoords[u*3+1]
                tempNeighbour[2] = solvateCoords[u*3+2]

                if(tempNeighbour[0] - solvateCoords[r*3] > boxVectors[0]/2):
                    tempNeighbour[0] = tempNeighbour[0] - boxVectors[0]
                elif(tempNeighbour[0] - solvateCoords[r*3] <= -boxVectors[0]/2):
                    tempNeighbour[0] = tempNeighbour[0] + boxVectors[0]
                if(tempNeighbour[1] - solvateCoords[r*3 + 1] > boxVectors[1]/2):
                    tempNeighbour[1] = tempNeighbour[1] - boxVectors[1]
                elif(tempNeighbour[1] - solvateCoords[r*3 + 1] <= -boxVectors[1]/2):
                    tempNeighbour[1] = tempNeighbour[1] + boxVectors[1]
                if(tempNeighbour[2] - solvateCoords[r*3 + 2] > boxVectors[2]/2):
                    tempNeighbour[2] = tempNeighbour[2] - boxVectors[2]
                elif(tempNeighbour[2] - solvateCoords[r*3 + 2] <= -boxVectors[2]/2):
                    tempNeighbour[2] = tempNeighbour[2] + boxVectors[2]

                distance = ((solvateCoords[r*3] - tempNeighbour[0])**2 +\
                            (solvateCoords[r*3+1] - tempNeighbour[1])**2 +\
                            (solvateCoords[r*3+2] -\
                            tempNeighbour[2])**2)**(1./2)

                if(distance != 0.0):
                    sumPot += 1/(distance**4)
        i += 1
    return sumPot


def potential_c_jac(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours,
                np.ndarray[np.float64_t, ndim=1] boxVectors):

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
    cdef tempNeighbour = np.array([0, 0, 0], dtype=float)
    centersXYZ = topol.trj.xyz[0][centers]*10

    for i in range(solvateLen):
        # Harmonic potential
        denom = ((solvateCoords[3*i] - centersXYZ[i][0])**2 +\
                 (solvateCoords[3*i+1] - centersXYZ[i][1])**2 +\
                 (solvateCoords[3*i+2] - centersXYZ[i][2])**2)**(1./2)

        jac[3*i] = 10*(solvateCoords[3*i]-centersXYZ[i][0])*(denom-2.2)/denom
        jac[3*i+1] = 10*(solvateCoords[3*i+1]-centersXYZ[i][1])*(denom-2.2)/denom
        jac[3*i+2] = 10*(solvateCoords[3*i+2]-centersXYZ[i][2])*(denom-2.2)/denom       
        
        # Solvate - neighbours pair-potential
        for d in range(centerNumNeighbours[i]):
            tempNeighbour = centerNeighbours[k]

            denom = ((solvateCoords[3*i] - tempNeighbour[0])**2 +\
                     (solvateCoords[3*i+1] - tempNeighbour[1])**2 +\
                     (solvateCoords[3*i+2] - tempNeighbour[2])**2)**3

            jac[3*i] += -4*(solvateCoords[3*i] - tempNeighbour[0])/denom
            jac[3*i+1] += -4*(solvateCoords[3*i+1] - tempNeighbour[1])/denom
            jac[3*i+2] += -4*(solvateCoords[3*i+2] - tempNeighbour[2])/denom

            k += 1

        #Solvate - solvate pair-potential
        for j in range(solvateLen):
            if(i != j): #Don't include solvate[i] - solvate[i]
                tempNeighbour[0] = solvateCoords[j*3]
                tempNeighbour[1] = solvateCoords[j*3+1]
                tempNeighbour[2] = solvateCoords[j*3+2]

                if(tempNeighbour[0] - solvateCoords[i*3] > boxVectors[0]/2):
                    tempNeighbour[0] = tempNeighbour[0] - boxVectors[0]
                elif(tempNeighbour[0] - solvateCoords[i*3] <= -boxVectors[0]/2):
                    tempNeighbour[0] = tempNeighbour[0] + boxVectors[0]
                if(tempNeighbour[1] - solvateCoords[i*3 + 1] > boxVectors[1]/2):
                    tempNeighbour[1] = tempNeighbour[1] - boxVectors[1]
                elif(tempNeighbour[1] - solvateCoords[i*3 + 1] <= -boxVectors[1]/2):
                    tempNeighbour[1] = tempNeighbour[1] + boxVectors[1]
                if(tempNeighbour[2] - solvateCoords[i*3 + 2] > boxVectors[2]/2):
                    tempNeighbour[2] = tempNeighbour[2] - boxVectors[2]
                elif(tempNeighbour[2] - solvateCoords[i*3 + 2] <= -boxVectors[2]/2):
                    tempNeighbour[2] = tempNeighbour[2] + boxVectors[2]

                denom = ((solvateCoords[3*i] - tempNeighbour[0])**2 +\
                         (solvateCoords[3*i+1] - tempNeighbour[1])**2 +\
                         (solvateCoords[3*i+2] - tempNeighbour[2])**2)**3

                jac[3*i] += (-4)*(solvateCoords[3*i] - tempNeighbour[0])/denom
                jac[3*i+1] += (-4)*(solvateCoords[3*i+1] - tempNeighbour[1])/denom
                jac[3*i+2] += (-4)*(solvateCoords[3*i+2] - tempNeighbour[2])/denom

                denom = ((tempNeighbour[0] - solvateCoords[3*i])**2 +\
                         (tempNeighbour[1] - solvateCoords[3*i+1])**2 +\
                         (tempNeighbour[2] - solvateCoords[3*i+2])**2)**3

                jac[3*i] += 4*(tempNeighbour[0] - solvateCoords[3*i])/denom
                jac[3*i+1] += 4*(tempNeighbour[1] - solvateCoords[3*i+1])/denom
                jac[3*i+2] += 4*(tempNeighbour[2] - solvateCoords[3*i+2])/denom
    return jac      