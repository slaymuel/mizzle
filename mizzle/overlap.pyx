"""Calculates minimum distance between atoms

"""

import cython
import numpy as np
cimport numpy as np
from radish import Topologizer

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t
@cython.cdivision(True)

def shortest_distance(np.ndarray[np.float64_t, ndim=2] solvCoords, np.ndarray elements, np.ndarray[np.float32_t, ndim=2] structCoords):
    cdef int i = 0
    cdef int j = 0
    cdef int k = 0
    cdef double minODist = 9999.999
    cdef double minHDist = 9999.999
    cdef double minStructDist = 9999.999
    cdef double distance = 0

    while i < len(solvCoords):
        j = i + 1
        while j < len(solvCoords):
            if(elements[i] == 'O' and elements[j] == 'O'):
                distance = ((solvCoords[i][0] - solvCoords[j][0])**2 + \
                            (solvCoords[i][1] - solvCoords[j][1])**2 + \
                            (solvCoords[i][2] - solvCoords[j][2])**2)**(1.0/2)
                if(distance < minODist):
                    minODist = distance

            elif(elements[i] == 'H' and elements[j] == 'H'):
                distance = ((solvCoords[i][0] - solvCoords[j][0])**2 + \
                            (solvCoords[i][1] - solvCoords[j][1])**2 + \
                            (solvCoords[i][2] - solvCoords[j][2])**2)**(1.0/2)
                if(distance < minHDist):
                    minHDist = distance

            j += 1
        i += 1
    i = 0
    while i < len(solvCoords):
        while k < len(structCoords):
            distance = ((solvCoords[i][0] - structCoords[k][0])**2 + \
                (solvCoords[i][1] - structCoords[k][1])**2 + \
                (solvCoords[i][2] - structCoords[k][2])**2)**(1.0/2)
            if(distance < minStructDist):
                minStructDist = distance
            k += 1
        i += 1
    return minODist, minHDist, minStructDist