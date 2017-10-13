import numpy as np
cimport numpy as np
from radish import Topologizer


FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t

def potential_c(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol, 
                np.ndarray[np.float64_t, ndim=1] centerNeighboursCoords,
                np.ndarray[ITYPE_t, ndim=1] centerCoordination,
                np.ndarray[np.float32_t, ndim=1] atoms):

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

        cdef int i = 0
        cdef int u = 0
        cdef int r = 0
        cdef int k = 0
        cdef int d = 0
        cdef float sumPot = 0
        cdef np.ndarray[np.float32_t, ndim=2] centersXYZ
        cdef float distance = 0
        cdef int solvateLen = len(solvateCoords)/3

        centersXYZ = topol.trj.xyz[0][centers]*10   #Does not drop duplicates(which is a good thing)

        for r in range(solvateLen):
            sumPot += 5*(((solvateCoords[r*3] - centersXYZ[i][0])**2 +\
                      (solvateCoords[r*3+1] - centersXYZ[i][1])**2 +\
                      (solvateCoords[r*3+2] - centersXYZ[i][2])**2)**(1./2) -\
                      2.2)**2

            for d in range(len(atoms)/3):
                distance = ((solvateCoords[r*3] - atoms[d*3])**2 +\
                           (solvateCoords[r*3+1] - atoms[d*3+1])**2 +\
                           (solvateCoords[r*3+2] - atoms[d*3+2])**2)**(1./2) 
                sumPot += 1/(distance**4)

            #for d in range(centerCoordination[r]):
            #    distance = ((solvateCoords[r*3] - centerNeighboursCoords[k*3])**2 +\
            #               (solvateCoords[r*3+1] - centerNeighboursCoords[k*3+1])**2 +\
            #               (solvateCoords[r*3+2] - centerNeighboursCoords[k*3+2])**2)**(1./2) 
            #    sumPot += 1/(distance**4)

                k += 1

            for u in range(solvateLen):
                distance = ((solvateCoords[r*3] - solvateCoords[u*3])**2 +\
                            (solvateCoords[r*3+1] - solvateCoords[u*3+1])**2 +\
                            (solvateCoords[r*3+2] - solvateCoords[u*3+2])**2)**(1./2)
        
                if(distance != 0.0 and distance < 5):
                    sumPot += 2/(distance**4)
            i += 1
        return sumPot