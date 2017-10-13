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
                np.ndarray[np.float32_t, ndim=1] atoms,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours):

    """Potential between solvate molecules and their
    neighbours. Minimization of this potential yields suitable 
    coordinates for solvate molecules.

    Parameters
    ----------
    solvate : 3N*array(float)
        Independant variables of the potential function
    centers : ndarray(float)
        Coordinates for centers which binds solvate

    Returns
    -------
    sumPot
        Value of the potential
    """
#unsigned int
    cdef int i = 0
    cdef int u = 0
    cdef int r = 0
    cdef int k = 0
    cdef int d = 0
    cdef double sumPot = 0
    cdef np.ndarray[np.float32_t, ndim=2] centersXYZ = topol.trj.xyz[0][centers]*10   #Does not drop duplicates(which is a good thing)
    cdef double distance = 0
    cdef int solvateLen = len(solvateCoords)/3

    for r in range(solvateLen):
        distance = ((solvateCoords[r*3] - centersXYZ[i][0])**2 +\
                    (solvateCoords[r*3+1] - centersXYZ[i][1])**2 +\
                    (solvateCoords[r*3+2] - centersXYZ[i][2])**2)**(1./2)
        sumPot += 5*(distance - 2.2)**2

        #for d in range(len(atoms)/3):
        #    distance = ((solvateCoords[r*3] - atoms[d*3])**2 +\
        #               (solvateCoords[r*3+1] - atoms[d*3+1])**2 +\
        #               (solvateCoords[r*3+2] - atoms[d*3+2])**2)**(1./2) 
        #    sumPot += 1/(distance**4)

        for d in range(centerNumNeighbours[r]):
            distance = ((solvateCoords[r*3] - centerNeighbours[k][0])**2 +\
                        (solvateCoords[r*3+1] - centerNeighbours[k][1])**2 +\
                        (solvateCoords[r*3+2] - centerNeighbours[k][2])**2)**(1./2) 
            sumPot += 1/(distance**4)
            k += 1

        for u in range(solvateLen):
            distance = ((solvateCoords[r*3] - solvateCoords[u*3])**2 +\
                        (solvateCoords[r*3+1] - solvateCoords[u*3+1])**2 +\
                        (solvateCoords[r*3+2] - solvateCoords[u*3+2])**2)**(1./2)
    
            if(distance != 0.0 and distance < 5):
                sumPot += 1/(distance**4)
        i += 1
    return sumPot


def potential_c_jac(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol,
                np.ndarray[np.float32_t, ndim=1] atoms,
                np.ndarray[np.float64_t, ndim=2] centerNeighbours,
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours):

    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef int d = 0
    cdef np.ndarray[np.float32_t, ndim=2] centersXYZ
    cdef double denom = 0
    cdef int solvateLen = len(solvateCoords)/3
    cdef np.ndarray[np.float64_t, ndim=1] jac = np.empty(solvateLen*3, dtype=np.float64)

    centersXYZ = topol.trj.xyz[0][centers]*10

    for i in range(solvateLen):
        denom = ((solvateCoords[3*i] - centersXYZ[i][0])**2 +\
                 (solvateCoords[3*i+1] - centersXYZ[i][1])**2 +\
                 (solvateCoords[3*i+2] - centersXYZ[i][2])**2)**(1./2)

        jac[3*i] += 2*(solvateCoords[3*i]-centersXYZ[i][0])*(denom-2.2)/denom
        jac[3*i+1] += 2*(solvateCoords[3*i+1]-centersXYZ[i][1])*(denom-2.2)/denom
        jac[3*i+2] += 2*(solvateCoords[3*i+2]-centersXYZ[i][2])*(denom-2.2)/denom       
        

        for d in range(centerNumNeighbours[i]):
            denom = ((solvateCoords[3*i] - centerNeighbours[k][0])**2 +\
                     (solvateCoords[3*i+1] - centerNeighbours[k][1])**2 +\
                     (solvateCoords[3*i+2] - centerNeighbours[k][2])**2)**3

            jac[3*i] += -4*(solvateCoords[3*i] - centerNeighbours[k][0])/denom
            jac[3*i+1] += -4*(solvateCoords[3*i+1] - centerNeighbours[k][1])/denom
            jac[3*i+2] += -4*(solvateCoords[3*i+2] - centerNeighbours[k][2])/denom

            k += 1

        for j in range(solvateLen):
            if(i != j):
                denom = ((solvateCoords[3*i] - solvateCoords[3*j])**2 +\
                         (solvateCoords[3*i+1] - solvateCoords[3*j+1])**2 +\
                         (solvateCoords[3*i+2] - solvateCoords[3*j+2])**2)**3

                jac[3*i] += (-4)*(solvateCoords[3*i] - solvateCoords[3*j])/denom
                jac[3*i+1] += (-4)*(solvateCoords[3*i+1] - solvateCoords[3*j+1])/denom
                jac[3*i+2] += (-4)*(solvateCoords[3*i+2] - solvateCoords[3*j+2])/denom

                denom = ((solvateCoords[3*j] - solvateCoords[3*i])**2 +\
                         (solvateCoords[3*j+1] - solvateCoords[3*i+1])**2 +\
                         (solvateCoords[3*j+2] - solvateCoords[3*i+2])**2)**3

                jac[3*i] += 4*(solvateCoords[3*j] - solvateCoords[3*i])/denom
                jac[3*i+1] += 4*(solvateCoords[3*j+1] - solvateCoords[3*i+1])/denom
                jac[3*i+2] += 4*(solvateCoords[3*j+2] - solvateCoords[3*i+2])/denom
    return jac      