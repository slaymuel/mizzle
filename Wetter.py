'''
Questions:
What is TER at the end of pdb file, only for proteins/aminoacid chains?
Ideas for 4-coordinated centers
'''


'''
To-do list: 
. Find close atoms (recursive function and bondgraph?, use networkX?) (maybe not needed)
. Particle swarm optimization, checkout pyswarm
. Provide Jacobian to minimize to improve speed
. Not enough with only neighbours for minimization (maybe it is, seems to work when changed potential function)
. 4-coordinated splits 1 and adsorbs one
. Topologizer.extract throws error if no atoms found
. Boundary conditions, affects distance between atoms
. Tkinter GUI?

'''


"""Wetter module

This module demonstrates documentation as specified by the `NumPy
Documentation HOWTO`_. Docstrings may extend over multiple lines. Sections
are created with a section header followed by an underline of equal length.

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
# pymatgen
# ase
import numpy as np
import pandas as pd
from pdbExplorer import append_atoms
import AtomCreator
import random
#import quaternion as quat
from pyquaternion import Quaternion
from radish import Topologizer
from tqdm import tqdm
from IPython import embed
from scipy.optimize import minimize
import sys
from IPython import embed
from timeit import default_timer as timer


class Wetter:

    def __init__(self, verbose, topol, theta=104.5, waterFrac=1, 
        hydroxylFrac=0, MOBondlength=2.2, HOHBondlength=1, OHBondlength=1, 
        centers=['Ti'], Nmax=6):

        #Input parameters
        self.theta = theta
        self.waterFrac = waterFrac
        self.hydroxylFrac = hydroxylFrac
        self.verbose = verbose
        self.MOBondlength = MOBondlength
        self.HOHBondlength = HOHBondlength
        self.OHBondlength = OHBondlength
        self.theta = theta/360*2*np.pi
        self.centers = centers
        self.Nmax = Nmax
        #Prepare input file for processing

        #Use radish (RAD-algorithm) to compute the coordination of each atom
        #self.topol = Topologizer.from_coords(file)
        self.topol = topol
        #self.topol.topologize()

        #Set Nmax which is the maximum coordination (coordination in bulk)
        #self.Nmax = self.topol.bondgraph['i'].value_counts().max()
        #Format float precision(will remove from here later maybe....)
        self.float_format = lambda x: "%.3f" % x

    def get_center_neighbours(self, coordination):
        """Construct neighbourgraph where columns are metal centers and rows 
        their coordination shell

        Parameters
        ----------
        coordination : int
            The coordination for which to construct the neighbourgraph and
            the indices for all such metal centers

        Returns
        -------
        centerIndices
            Indices for all metal centers if there are any, otherwise empty
        neighbourgrap
            neighbourgraph if any atoms found, empty otherwise
        """
        try:
            centerIndices = self.topol.extract(self.centers[0], environment =
                {'O': self.Nmax - coordination}).index.get_level_values(1)

            neighbourgraph = self.topol.bondgraph.loc[self.topol.\
                bondgraph['j'].isin(centerIndices)].\
                sort_values(['j'])[['i', 'j']].pivot(columns= 'j').\
                apply(lambda x: pd.Series(x.dropna().values)).\
                apply(np.int64)['i']

            return (centerIndices, neighbourgraph)

        except IndexError:  #If no atoms found, extract() throws IndexError
            return [], []

    def overlap(self, point, points):
        i = 0
        for atom in points:
            distance = np.linalg.norm(point - atom)
            if(distance < 2):
                i += 1
        if(i > 0):
            print("Atom overlaps with " + str(i) + " other atoms.")
            return True
        else:
            return False

    def find_close_atoms(self, point, atomIndices):
        i = 0
        closeAtoms = np.empty([0, 3], dtype = float)
        
        for atom in atomIndices:
            distance = np.linalg.norm(self.topol.trj.xyz[0][point] -\
                                      self.topol.trj.xyz[0][atom])
            if(distance < 5):
                np.vstack((closeAtoms, atom))
        return closeAtoms

    def potential(self, solvate, centers):
        """Repulsion potential between solvate molecules and their
        neighbours. Minimization of this potential yields suitable 
        coordinates for solvate molecules.

        Parameters
        ----------
        solvate : 3N*array(float)
            Independant variables in the potential function
        centers : ndarray(float)
            Coordinates for centers which binds solvate

        Returns
        -------
        sumPot
            Value of the potential
        """
        A = 1
        cutoff = 4  #Not needed since only neighbours are included
        sumPot = 0

        solvate = np.reshape(solvate, (-1, 3))

        i = 0
        for r in solvate:
            centerNeighbours = np.array(self.topol.bondgraph.loc[self.topol.bondgraph['j'] == centers[i]]['i'])
            centerNeighbours = self.topol.trj.xyz[0][centerNeighbours]*10
            center = self.topol.trj.xyz[0][centers[i]]*10

            sumPot += (np.linalg.norm(r - center) - 2.2)**2  #Harmonic potential, places solvate on correct distance from M
            sumPot += np.sum(1/((np.sum((centerNeighbours - r)**2,axis=-1)**(1./2))**4))    #M-O pairpotential
            #for d in centerNeighbours:
            #    sumPot += 1/(np.linalg.norm(r - d)**4)
            sumPot += np.sum(1/((np.sum((solvate[np.arange(len(solvate)) != i] - r)**2, axis=-1)**(1./2))**4))  #O-O pairpotential
            i += 1

        # i = 0
        # j = 0
        # while i < len(solvate):
        #     j = i + 1
        #     while j < len(solvate):
        #         sumPot += 2/(np.linalg.norm(solvate[i] - solvate[j])**4)
        #         j += 1
        #     i += 1

        return sumPot

    def minimization_callback(self, x):
        print("One iteration done!")
        #print("sumPot: " + str(self.potential(x, self.topol.trj.xyz[0])))
        pass

    def optimize(self, coords, centers, addedSolvate = [], cons = []):
        """Set up for minimization scheme by the SLSQP algorithm

        Parameters
        ----------
        coords : 3*N array(float)
            The coordination for which to construct the neighbourgraph and
            the indices for all such metal centers
        centers : array(int)
            Array containing all the centers which binds solvate

        Raises
        ------
        ValueError
            If one or more coordinates are missing centers

        Returns
        -------
        coords
            Optimized oxygen coordinates
        vectors
            Optimized M-O vectors
        """
        vectors = np.empty([0, 3], dtype = float)

        centerNeighbours = np.array(self.topol.bondgraph.loc[self.topol.bondgraph['j'].isin(centers)]['i'])
        centerNeighbours = self.topol.trj.xyz[0][centerNeighbours]*10

        M = self.topol.trj.xyz[0][centers]*10

        cons = ({'type': 'eq',
                        'fun' : lambda x: np.sum(np.absolute(np.apply_along_axis(np.linalg.norm, 1, np.reshape(x, (-1, 3)) - M) - self.MOBondlength))})
        if(len(M) != len(coords)):
            raise ValueError("Optimization failed: some solvate molecules were not assigned to a center.")
        start = timer()
        self.potential(coords.flatten(), centers)
        end = timer()
        print("Potential takes: " + str(end-start) + " seconds to calculate")

        #print("Before: " + str(coords))
        
        res = minimize(self.potential, coords,
                        #constraints = cons,
                        args = centers,#centerNeighbours,#np.vstack((addedSolvate, centerNeighbours)),
                        callback = self.minimization_callback,
                        method = 'BFGS',
                        options={'disp': True, 'gtol': 1e-03, 'eps': 1.4901161193847656e-08, 'return_all': False, 'maxiter': None, 'norm': np.inf})
                        #options={'disp': True, 'iprint': 1, 'eps': 1.4901161193847656e-08, 'maxiter': 100, 'ftol': 1e-06})
        print(res)
        coords = np.reshape(res.x, (-1, 3))

        i = 0
        for coord in coords:
            vectors = np.vstack((vectors, (coord-M[i])/np.linalg.norm(coord-M[i])))
            i += 1

        return coords, vectors


    #For Nmax - 2 coordinated centers calculate pair-vectors
    def calculate_pair_vectors(self):
        """Calculates the sum of vectors from oxygen to 4-coordinated metal
        atoms. These vectors are used as first guesses of the positions for
        the hydration molecules.

        Parameters
        ----------
        addedVectors : numpy array
            Description bla bla.

        Returns
        -------
        bool
            True always
        """
        centerIndices, neighbourgraph = self.get_center_neighbours(2)

        if(len(centerIndices) > 0):
            centers = np.empty([0], dtype=int)
            print("Found " + str(len(centerIndices)) +\
                " centers with Nmax - 2 coordination")

            vectors = np.empty([0, 3], dtype = float)
            coordinates = np.empty([0, 3], dtype = float)

            randIndices = random.sample(range(0, len(centerIndices)),
                                        int((self.waterFrac +\
                                        self.hydroxylFrac)*\
                                        float(len(centerIndices))))
            indices = centerIndices[randIndices]

            for center in indices:
                points = np.empty([0, 3], dtype=float)
                tempVectors = np.empty([0, 3], dtype=float)
                pairVectors = np.empty([0, 3], dtype=float)
                sumVec = np.array([0, 0, 0], dtype=float)
                points = np.vstack((points, self.topol.trj.xyz[0][center]*10))

                for neighbour in neighbourgraph[center]:
                    points = np.vstack((points, 
                                        self.topol.trj.xyz[0][neighbour]*10))
                    tempVec = self.topol.trj.xyz[0][center] -\
                        self.topol.trj.xyz[0][neighbour]
                    tempVec = tempVec/np.linalg.norm(tempVec)

                    tempVectors = np.vstack((tempVectors, tempVec))
                    sumVec += tempVec

                #Normalize the sum of vectors
                sumVec = sumVec/np.linalg.norm(sumVec)

                for vector in tempVectors:
                    
                    tempVector = sumVec - vector
                    pairVectors = np.vstack((pairVectors, tempVector))

                #Sort vectors with increasing norm
                pairVectors = sorted(pairVectors, 
                                    key = lambda x: np.sqrt(x[0]**2 + x[1]**2\
                                                            + x[2]**2))

                pairVec1 = pairVectors[0]/np.linalg.norm(pairVectors[0])
                pairVec2 = pairVectors[1]/np.linalg.norm(pairVectors[1])

                #If there are only 2 neighbours, rotate Pi/2 around 
                #directional vector to maximize distance to neighbours

                if(len(tempVectors) == 2):
                    axis = sumVec
                    angle = np.pi/2
                    pairVec1 = Quaternion(axis = axis,angle = angle).\
                        rotate(pairVec1)
                    pairVec2 = Quaternion(axis = axis,angle = angle).\
                        rotate(pairVec2)


                crossProd = np.cross(pairVec1, pairVec2)
                dotProdVec1 = np.dot(sumVec, pairVec1)
                dotProdVec2 = np.dot(sumVec, pairVec2)
                if(dotProdVec1 < 0):
                    pairVec1 = Quaternion(axis = crossProd,angle = -np.pi/7).\
                        rotate(pairVec1)
                    pairVec2 = Quaternion(axis = crossProd,angle = np.pi/7).\
                        rotate(pairVec2)
                else:
                    pairVec2 = Quaternion(axis = crossProd,angle = -np.pi/7).\
                        rotate(pairVec2)
                    pairVec1 = Quaternion(axis = crossProd,angle = np.pi/7).\
                        rotate(pairVec1)


                coord1 = self.topol.trj.xyz[0][center]*10 +\
                         pairVec1*self.MOBondlength
                coord2 = self.topol.trj.xyz[0][center]*10 +\
                         pairVec2*self.MOBondlength

                #Minimize potential function
                M = self.topol.trj.xyz[0][center]*10
                cons = ({'type': 'eq',
                     'fun' : lambda x: np.array(np.linalg.norm([x[0], x[1], x[2]] - M) - self.MOBondlength)},
                     {'type': 'eq',
                     'fun' : lambda x: np.array(np.linalg.norm([x[3], x[4], x[5]] - M) - self.MOBondlength)})
                #[coord1, coord2], [pairVec1, pairVec2] = self.optimize(np.vstack((coord1, coord2)), [center, center], coordinates, cons)
                # #embed()
                # res = minimize(self.potential_function, np.hstack((coord1, coord2)), constraints = cons, args=(np.vstack((points, coordinates))))

                # if(np.array_equal(coord1, res.x[:3]) and np.array_equal(coord2, res.x[3:])):
                #     print("minimization did nothing...")
                
                # coord1 = res.x[:3]
                # coord2 = res.x[3:]

                # pairVec1 = (coord1-M)/np.linalg.norm(coord1-M)
                # pairVec2 = (coord2-M)/np.linalg.norm(coord2-M)

                #Save relevant vectors
                vectors = np.vstack((vectors, pairVec1))
                vectors = np.vstack((vectors, pairVec2))
                centers = np.append(centers, center)
                centers = np.append(centers, center)
                coordinates = np.vstack((coordinates, coord1))
                coordinates = np.vstack((coordinates, coord2))

            return vectors, coordinates, centers

        else:
            return (np.empty([0, 3], dtype=float),
                    np.empty([0, 3], dtype=float),
                    np.empty([0], dtype=int))

    #Calculate M-O vectors which determine the positions of hydration molecules
    def calculate_vectors(self):

        vectors = np.empty([0, 3], dtype=float)
        coords = np.empty([0, 3], dtype=float)
        vec = np.array([0,0,0])

        #Get indices for metal centers with coordination Nmax - 1
        centerIndices, neighbourgraph = self.get_center_neighbours(1)

        if(len(centerIndices) > 0):
            centers = np.empty([0], dtype=int)

            print("Found " + str(len(centerIndices)) + " centers with Nmax - 1 coordination")

            #Calculate only for user specified fraction
            randIndices = random.sample(range(0, len(centerIndices)), int((self.waterFrac + self.hydroxylFrac) * float(len(centerIndices))))
            indices = centerIndices[randIndices]

            #Calculate M-O vectors
            for center in tqdm(indices, ascii=True, desc='Calculating directional vectors'):
                vec = [0, 0, 0]
                for neighbour in neighbourgraph[center]:

                    #M-O vector
                    tempVec = self.topol.trj.xyz[0][center] -\
                              self.topol.trj.xyz[0][neighbour]

                    #normalize
                    tempVec = tempVec/np.linalg.norm(tempVec)

                    #Sum vectors
                    #vec = [vec[0] + tempVec[0], vec[1] + tempVec[1], vec[2] + tempVec[2]]
                    vec = vec + tempVec


                if(np.linalg.norm(vec) < 0.1):
                    print("\n\nWarning: directional vector almost 0, will " +\
                          "not hydrate at atom with index: " +\
                           str(center + 1))
                    print("Probably due to symmetric center..........  \n")

                else:
                    centers = np.append(centers, center)
                    vec = vec/np.linalg.norm(vec)
                    vectors = np.vstack((vectors, vec))

                    #Set vector length to get coordinate
                    vec = vec*self.MOBondlength

                    #Coordinates where oxygen should be placed: vec + coordinate of metal atom multiplied with bondlength
                    coord = self.topol.trj.xyz[0][center]*10 + vec

                    #Save coords and correct units
                    coords = np.vstack((coords, coord))

            return vectors, coords, centers
        else:
            return (np.empty([0, 3], dtype=float),
                    np.empty([0, 3], dtype=float),
                    np.empty([0], dtype=int))

    def wet(self):
        #A sum of all vectors pointing from each neighbour of each center to the center itself
        #get coordinates where oxygen should be placed

        vectors, coords, centers = self.calculate_vectors()
        pairVectors, pairCoords, pairCenters = self.calculate_pair_vectors()
        #vectors, coords, centers = self.calculate_pair_vectors()

        vectors = np.concatenate((vectors, pairVectors), axis=0)
        coords = np.concatenate((coords, pairCoords), axis=0)
        centers = np.concatenate((centers, pairCenters), axis=0)

        start = timer()
        coords, vectors = self.optimize(coords, centers)
        end = timer()

        print("Optimization took: " + str(end-start))

        #Check for overlap

        np.set_printoptions(suppress=True) #print whole number and not 1.234 e+4 etc, will remove later....
        
        randIndices = random.sample(range(0, len(coords)), 
                                    int(self.hydroxylFrac *\
                                        float(len(coords))))

        #Add hydroxide to user specified fraction
        hydCoords = coords[randIndices]
        hydVectors = vectors[randIndices]
        hydCoords, hydElements = AtomCreator.add_hydroxyl(hydCoords, 
                                                          hydVectors, 
                                                          self.theta)

        print("Adding: " + str(len(hydCoords)) + " atoms ("+\
              str(len(hydCoords)/3) +" hydroxyl molecules)")

        #Create mask for the selection of water molecules
        mask = np.ones(len(coords), np.bool)
        mask[randIndices] = 0

        #Create water molecules at the inverse selection of where hydroxide was placed
        watCoords = coords[mask]
        watVectors = vectors[mask]
        watCoords, watElements = AtomCreator.add_water(watCoords, 
                                                       watVectors, 
                                                       self.theta)

        print("Adding: " + str(len(watCoords)) + " atoms ("+\
              str(len(watCoords)/3) +" water molecules)")

        #Concatenate water and hydroxide coordinates
        coords = np.concatenate((watCoords, hydCoords))
        elements = np.concatenate((watElements, hydElements))

        #append_atoms(file = self.file, coords = atoms, elements = elements)
        return coords, elements
        #self.add_water(pairCoords, pairVectors)
        #embed()