'''
Questions:
What is TER at the end of pdb file, only for proteins/aminoacid chains?
Ideas for 4-coordinated centers
'''


'''
To-do list: 
. Rewrite config file for better readability
. Boundary conditions, affects distance between atoms
. Tkinter GUI?

'''

'''
Notes:
Repulsion from center is not included
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

    $ python main.py  config.wet TiO_110.pdb


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
import pyximport; pyximport.install()
import potential
import mdtraj as md

class Wetter:

    def __init__(self, verbose, topol, theta=104.5, highWaterFrac=0.5, 
        highHydroxylFrac=0.5, lowFrac = 0.5, MOBondlength=2.2, HOHBondlength=1, OHBondlength=1, 
        center='Ti', Nmax=6):

        #Input parameters
        self.theta = theta

        self.verbose = verbose
        self.MOBondlength = MOBondlength
        self.HOHBondlength = HOHBondlength
        self.OHBondlength = OHBondlength
        self.theta = theta/360*2*np.pi
        self.center = center
        self.Nmax = Nmax
        self.lowFrac = float(lowFrac)
        self.highWaterFrac = float(highWaterFrac)
        self.highHydroxylFrac = float(highHydroxylFrac)

        #Use radish (RAD-algorithm) to compute the coordination of each atom
        self.topol = topol
        #self.topol.topologize()

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
            centerIndices = self.topol.extract(self.center, environment =
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


    def optimize(self, coords, centers):
        """Set up for minimization scheme for the L-BFGS-B algorithm

        Parameters
        ----------
        coords : 3*N array(float)
            The coordination for which to construct the neighbourgraph and
            the indices for all such metal centers
        centers : array(int)
            Array containing all the centers which has bound each solvate

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
        cutoff = 6.0
        M = self.topol.trj.xyz[0][centers]*10

        if(len(M) != len(coords)):
            raise ValueError("Optimization failed: some solvate molecules were not assigned to a center.")        

        # Drops duplicates
        centerCoordination = np.array(self.topol.bondgraph.\
                             loc[self.topol.bondgraph['j'].\
                             isin(centers)]['j'].value_counts())
        missing = len(centers) - len(centerCoordination)

        # Add removed duplicates
        if(missing > 0):
            centerCoordination = np.append(centerCoordination,\
                                           centerCoordination[\
                                           len(centerCoordination)-missing:])

        centerNeighbours = np.empty([0, 3], dtype = float)
        centerNumNeighbours = np.empty([0], dtype = int)

        # Get neighbours to each center
        for center in centers:
            neighbours = md.compute_neighbors(self.topol.trj, cutoff/10.0,\
                                              np.array([center]))
            centerNumNeighbours = np.hstack((centerNumNeighbours,\
                                             len(neighbours[0])))
            centerNeighbours = np.vstack((centerNeighbours,\
                                    self.topol.trj.xyz[0][neighbours[0]]*10))


        start = timer()
        potential.potential_c(coords.flatten(), centers, self.topol,\
                              centerNeighbours, centerNumNeighbours)
        end = timer()
        print("Potential takes: " + str(end-start) + " seconds to calculate")

        start = timer()
        potential.potential_c_jac(coords.flatten(), centers, self.topol,\
                                  centerNeighbours, centerNumNeighbours)
        end = timer()
        print("Potential Jacobian takes: " + str(end-start) +\
              " seconds to calculate")

        print("Minimizing potential......")
        # Run minimization
        res = minimize(potential.potential_c, coords,
                        args = (centers, self.topol, centerNeighbours,\
                                centerNumNeighbours),
                        jac = potential.potential_c_jac,
                        method = 'L-BFGS-B',
                        options={'disp': False, 'gtol': 1e-05, 'iprint': 0,\
                                 'eps': 1.4901161193847656e-05,\
                                 'maxiter': 1000})
        print res
        if(res.success):
            print ("\nFound adequate local minima!\n")
        else:
            print("\nDid not end up in local minima but probably the hydration went well! :)\n")

        # Since minimization returns flat array we need to reshape
        coords = np.reshape(res.x, (-1, 3))

        # Recalculate directional vectors
        i = 0
        for coord in coords:
            vectors = np.vstack((vectors, (coord-M[i])/\
                                 np.linalg.norm(coord-M[i])))
            i += 1

        return coords, vectors


    def calculate_pair_vectors(self):
        """Calculate coordinates and directional vectors

        Notes
        -----
        Calculates the sum of vectors from oxygen to 4-coordinated metal
        atoms. These vectors are used as first guesses of the positions for
        the hydration molecules.

        Returns
        -------
        vectors : ndarray(float)
            Directional vectors for each solvate
        coord : ndarray(float)
            Coordinates for oxygen atom in each solvate
        centers : ndarray(int)
            List of each center that each solvate belongs to
        """

        centerIndices, neighbourgraph = self.get_center_neighbours(2)

        if(len(centerIndices) > 0):
            centers = np.empty([0], dtype=int)
            print("Found " + str(len(centerIndices)) +\
                " centers with Nmax - 2 coordination")

            vectors = np.empty([0, 3], dtype = float)
            coordinates = np.empty([0, 3], dtype = float)

            randIndices = random.sample(range(0, len(centerIndices)),
                                        int(self.lowFrac*\
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

                # If there are only 2 neighbours, rotate Pi/2 around
                # directional vector to maximize distance to neighbours
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


    def calculate_vectors(self):
        """ See Wetter.calculate_pair_vectors
        """
        vectors = np.empty([0, 3], dtype=float)
        coords = np.empty([0, 3], dtype=float)
        vec = np.array([0,0,0])

        #Get indices for metal centers with coordination Nmax - 1
        centerIndices, neighbourgraph = self.get_center_neighbours(1)

        if(len(centerIndices) > 0):
            centers = np.empty([0], dtype=int)

            print("Found " + str(len(centerIndices)) +\
                  " centers with Nmax - 1 coordination")

            #Calculate only for user specified fraction
            randIndices = random.sample(range(0, len(centerIndices)),\
                                        int((self.highWaterFrac +\
                                            self.highHydroxylFrac)*\
                                            float(len(centerIndices))))
            indices = centerIndices[randIndices]

            #Calculate M-O vectors
            for center in tqdm(indices, ascii=True,\
                               desc='Calculating directional vectors'):
                vec = [0, 0, 0]
                for neighbour in neighbourgraph[center]:

                    #M-O vector
                    tempVec = self.topol.trj.xyz[0][center] -\
                              self.topol.trj.xyz[0][neighbour]
                    tempVec = tempVec/np.linalg.norm(tempVec)
                    vec = vec + tempVec

                #If norm of directional vector too small
                if(np.linalg.norm(vec) < 0.1):
                    print("\n\nWarning: directional vector almost 0, will " +\
                          "not hydrate at atom with index: " +\
                           str(center + 1))
                    print("Probably due to symmetric center..........  \n")

                else:
                    centers = np.append(centers, center)
                    vec = vec/np.linalg.norm(vec)
                    vectors = np.vstack((vectors, vec))
                    vec = vec*self.MOBondlength

                    #Save coordinates and correct units
                    coord = self.topol.trj.xyz[0][center]*10 + vec
                    coords = np.vstack((coords, coord))

            return vectors, coords, centers
        else:
            return (np.empty([0, 3], dtype=float),
                    np.empty([0, 3], dtype=float),
                    np.empty([0], dtype=int))

    def wet(self):
        """ Main algorithm

        Returns
        -------
        coords : ndarray(float)
            coordinates for each solvate atom
        element : array(char)
            List of elements
        """

        vectors, coords, centers = self.calculate_vectors()
        pairVectors, pairCoords, pairCenters = self.calculate_pair_vectors()
        numHighCoord = len(centers)
        numLowCoord = len(pairCenters)
        #vectors, coords, centers = self.calculate_pair_vectors()

        vectors = np.concatenate((vectors, pairVectors), axis=0)
        coords = np.concatenate((coords, pairCoords), axis=0)
        centers = np.concatenate((centers, pairCenters), axis=0)

        start = timer()
        coords, vectors = self.optimize(coords, centers)
        end = timer()

        print("Optimization took: " + str(end-start) + " seconds")

        np.set_printoptions(suppress=True)

        hydCoords = coords[:int(self.highHydroxylFrac*numHighCoord)]
        hydVectors = vectors[:int(self.highHydroxylFrac*numHighCoord)]
        watCoords = coords[int(self.highHydroxylFrac*numHighCoord):numHighCoord]
        watVectors = vectors[int(self.highHydroxylFrac*numHighCoord):numHighCoord]

        
        hydCoords = np.concatenate((hydCoords, coords[numHighCoord::2]), axis=0)
        hydVectors = np.concatenate((hydVectors, vectors[numHighCoord::2]), axis=0)
        watCoords = np.concatenate((watCoords, coords[numHighCoord + 1::2]), axis=0)
        watVectors = np.concatenate((watVectors, vectors[numHighCoord + 1::2]), axis=0)

        #randIndices = random.sample(range(0, len(coords)), 
        #                            int(self.hydroxylFrac *\
        #                                float(len(coords))))

        #Add hydroxide to user specified fraction
        #hydCoords = coords[randIndices]
        #hydVectors = vectors[randIndices]
        hydCoords, hydElements = AtomCreator.add_hydroxyl(hydCoords, 
                                                          hydVectors, 
                                                          self.theta)

        print("Adding: " + str(len(hydCoords)) + " atoms ("+\
              str(len(hydCoords)/2) +" hydroxyl molecules)")

        #Create mask for the selection of water molecules
        #mask = np.ones(len(coords), np.bool)
        #mask[randIndices] = 0

        #watCoords = coords[mask]
        #watVectors = vectors[mask]
        watCoords, watElements = AtomCreator.add_water(watCoords, 
                                                       watVectors, 
                                                       self.theta)

        print("Adding: " + str(len(watCoords)) + " atoms ("+\
              str(len(watCoords)/3) +" water molecules)")

        #Concatenate water and hydroxide coordinates and elements
        coords = np.concatenate((watCoords, hydCoords))
        elements = np.concatenate((watElements, hydElements))

        return coords, elements