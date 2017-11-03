"""Example NumPy style docstrings.

This Cython module contains the potential function and the Jacobian

Example
-------
Examples can be given using either the ``Example`` or ``Examples``
sections. Sections support any reStructuredText formatting, including
literal blocks::

    $ python example_numpy.py


Section breaks are created with two blank lines. Section breaks are also
implicitly created anytime a new section starts. Section bodies *may* be
indented:

Notes
-----
    This is an example of an indented section. It's like any other section,
    but the body is indented to help it stand out from surrounding text.

If a section is indented, then a section break is created by
resuming unindented text.

Attributes
----------
module_level_variable1 : int
    Module level variables may be documented in either the ``Attributes``
    section of the module docstring, or in an inline docstring immediately
    following the variable.

    Either form is acceptable, but the two should not be mixed. Choose
    one convention to document module level variables and be consistent
    with it.


.. _NumPy Documentation HOWTO:
   https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

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
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours):

    """Potential function

    Parameters
    ----------
    solvate : 3N*array(float)
        Independant variables of the potential function
    centers : ndarray(float)
        Coordinates for centers which binds solvate
    topol : Topologizer instance
        Topologizer instance of system
    centerNeighbours : ndarray(float)
        coordinates of neighbours to solvate
    centerNumNeighbours : array(int)
        number of neighbours to each solvate

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
                np.ndarray[ITYPE_t, ndim=1] centerNumNeighbours):

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