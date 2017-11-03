#cython: cdivision=True
#cython: boundscheck=False, wraparound=False, nonecheck=False
import cython
import numpy as np
cimport numpy as np

FTYPE = np.float32
ITYPE = np.int_
ctypedef np.float_t FTYPE_t
ctypedef np.int_t ITYPE_t

def get(np.ndarray[np.float32_t, ndim=2] centers, np.ndarray[np.float64_t, ndim=2] neighbours, np.ndarray[np.float64_t, ndim=1] boxVectors):
    cdef int i=0
    cdef int j=0

    for i in range(len(centers)):
        for j in range(len(neighbours)):
            if(neighbours[j][0] - centers[i][0] > boxVectors[0]/2):
                neighbours[j][0] = neighbours[j][0] - boxVectors[0]
            elif(neighbours[j][0] - centers[i][0] <= -boxVectors[0]/2):
                neighbours[j][0] = neighbours[j][0] + boxVectors[0]
            if(neighbours[j][1] - centers[i][1] > boxVectors[1]/2):
                neighbours[j][1] = neighbours[j][1] - boxVectors[1]
            elif(neighbours[j][1] - centers[i][1] <= -boxVectors[1]/2):
                neighbours[j][1] = neighbours[j][1] + boxVectors[1]
            if(neighbours[j][2] - centers[i][2] > boxVectors[2]/2):
                neighbours[j][2] = neighbours[j][2] - boxVectors[2]
            elif(neighbours[j][2] - centers[i][2] <= -boxVectors[2]/2):
                neighbours[j][2] = neighbours[j][2] + boxVectors[2]
    return neighbours