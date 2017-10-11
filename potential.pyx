import numpy as np
cimport numpy as np
from radish import Topologizer
cimport cython

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t

def potential_c(np.ndarray[np.float64_t, ndim=1] solvateCoords, 
                np.ndarray[ITYPE_t, ndim=1] centers, 
                object topol, 
                np.ndarray[np.float32_t, ndim=1] centerNeighboursXYZ, 
                ITYPE_t numCenters, 
                np.ndarray[ITYPE_t, ndim=1] centerCoordination):

        cdef int i = 0
        cdef int u = 0
        cdef int r = 0
        cdef float sumPot = 0
        cdef np.ndarray[np.float32_t, ndim=2] centersXYZ
        cdef float distance = 0
        cdef int solvateLen = len(solvateCoords)/3

        centersXYZ = topol.trj.xyz[0][centers]*10

        for r in range(solvateLen):
            sumPot += (((solvateCoords[r*3] - centersXYZ[i][0])**2 +\
                      (solvateCoords[r*3+1] - centersXYZ[i][1])**2 +\
                      (solvateCoords[r*3+2] - centersXYZ[i][2])**2)**(1./2) -\
                      2.2)**2

            for d in range(centerCoordination[r]):
                distance = ((solvateCoords[r*3] - centerNeighboursXYZ[(r+d)*3])**2 +\
                           (solvateCoords[r*3+1] - centerNeighboursXYZ[(r+d)*3+1])**2 +\
                           (solvateCoords[r*3+2] - centerNeighboursXYZ[(r+d)*3+2])**2)**(1./2) 
                sumPot += 1/(distance**4)

            for u in range(solvateLen):
                distance = ((solvateCoords[r*3] - solvateCoords[u*3])**2 +\
                            (solvateCoords[r*3+1] - solvateCoords[u*3+1])**2 +\
                            (solvateCoords[r*3+2] - solvateCoords[u*3+2])**2)**(1./2)
                
                if(distance != 0.0):
                    pass
                    #sumPot += 2/(distance**4)
            i += 1
        print sumPot
        return sumPot