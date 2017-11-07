import cython
import numpy as np
cimport numpy as np
from radish import Topologizer

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t
@cython.cdivision(True)

def shortest_distance(np.ndarray[np.float64_t, ndim=2] coords, np.ndarray elements):
    cdef int i = 0
    cdef int j = 0
    cdef double minODist = 9999.999
    cdef double minHDist = 9999.999
    cdef double distance = 0

    while i < len(coords):
        j = i + 1
        while j < len(coords):
            if(elements[i] == 'O' and elements[j] == 'O'):
                distance = ((coords[i][0] - coords[j][0])**2 + \
                            (coords[i][1] - coords[j][1])**2 + \
                            (coords[i][2] - coords[j][2])**2)**(1.0/2)
                if(distance < minODist):
                    minODist = distance

            elif(elements[i] == 'H' and elements[j] == 'H'):
                distance = ((coords[i][0] - coords[j][0])**2 + \
                            (coords[i][1] - coords[j][1])**2 + \
                            (coords[i][2] - coords[j][2])**2)**(1.0/2)
                if(distance < minHDist):
                    minHDist = distance

            j += 1
        i += 1
    return minODist, minHDist