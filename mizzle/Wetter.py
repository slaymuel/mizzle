"""Hydrates metal oxide surfaces

This is the main module that contains the methods for hydration of
metal oxide surfaces.

Example
-------
Using Wetter module in existing project::

    import mizzle.Wetter as mw
    wet = mw(verbose, radish_instance)
    wet.add_solvate({'Nmax': Nmax, 'element': element,\
'coordination': coordination, 'OH': hydroxylFrac, 'OH2': waterFrac, 'O':0.05})
    wet.optimize()
    wet.wet()
    wet.append_atoms('outputfile.pdb', 'residuename')

Attributes
----------
    hydVectors : ndarray
        Array that holds directional vectors for hydroxyl molecules.
    hydCoords
        Array that holds oxygen coordinates for hydroxyl molecules.
    hydCenters
        Array that holds the corresponding metal atom for each 
        hydroxyl molecule.
    watVectors : ndarray
        see hydVectors
    watCoords
        see hydCoords
    watCenters
        see hydCenters
    oxyVectors : ndarray
        see hydVectors
    oxyCoords
        hydCoords
    oxyCenters
        see hydCenters

Parameters
----------
    theta : double, optional
        HOH bond angle. Default is 104.5.
    MOBondlength : double, optional
        M-O bond length. Default is 2.2
    HOHBondlength : double, optional
        H-O bond length in water. Default is 1
    OHBondlength : double, optional
        O-H bond length in hydroxyl, Default is 1
"""

from __future__ import print_function
import numpy as np
import pandas as pd
import random
from radish import Topologizer
from tqdm import tqdm
from scipy.optimize import minimize
import sys
from timeit import default_timer as timer
import mdtraj as md
from IPython import embed
from scipy.optimize.optimize import _approx_fprime_helper
# mizzle imports
from mizzle import potential,ghosts
from mizzle.pdbExplorer import append_atoms,remove_low_coordinated
import mizzle.AtomCreator as mac
from mizzle.overlap import shortest_distance
from mizzle.puts import puts

if sys.version_info[0] == 2:
    # workaround for Sphinx autodoc bug
    import __builtin__
    def print(*args, **kwargs):
        __builtin__.print(*args, **kwargs)

class Wetter:

    def __init__(self, silent, topol,theta=104.5, MOBondlength=2.2,\
                 HOHBondlength=1, OHBondlength=1):

        self.hydVectors = np.empty([0, 3], dtype=float)
        self.hydCoords = np.empty([0, 3], dtype=float)
        self.hydCenters = np.empty([0], dtype=float)

        self.watVectors = np.empty([0, 3], dtype=float)
        self.watCoords = np.empty([0, 3], dtype=float)
        self.watCenters = np.empty([0], dtype=float)

        self.oxyVectors = np.empty([0, 3], dtype=float)
        self.oxyCoords = np.empty([0, 3], dtype=float)
        self.oxyCenters = np.empty([0], dtype=float)

        #Input parameters
        self.theta = theta
        self.MOBondlength = MOBondlength
        self.HOHBondlength = HOHBondlength
        self.OHBondlength = OHBondlength
        self.theta = theta/360*2*np.pi
        self.lowFrac = 1
        self.highWaterFrac = 1
        self.highHydroxylFrac = 1
        self.silent = silent
        #Use radish (RAD-algorithm) to compute the coordination of each atom
        self.topol = Topologizer.from_coords(topol)

        if self.topol.trj.unitcell_lengths is None:
            raise RuntimeError("No box vectors in PDB file")
        else:
            self.boxVectors = 10*self.topol.trj.unitcell_lengths.reshape(-1).astype(float)

        #Format float precision(will remove from here later maybe....)
        self.float_format = lambda x: "%.3f" % x

        #Set up verbose print function
        self.verboseprint = print if not silent else lambda *a, **k: None
        self.verboseputs = puts if not silent else lambda *a, **k: None
        self.i = 0
        self.potTime = None
        self.jacPotTime = None
        self.start = None
        self.end = None

        #temp
        self.centers = None
        self.centerNeighbours = None
        self.centerNumNeighbours = None

    def show_progress(self, x):
        #jac = potential.potential_c_jac(x, self.centers, self.topol,\
        #                          self.centerNeighbours, self.centerNumNeighbours, self.boxVectors)
        #print(np.linalg.norm(jac))
        self.end = timer()
        self.i += 1
        if(self.i % 10 == 0 or self.i == 1):
            self.verboseprint('\r', end='')
            sys.stdout.flush()
            #self.verboseprint("Maximum time left: " + str(self.end - self.start)*(500-self.i) + 's', end='')
            self.verboseprint("Maximum time left: " + str(int((500 - self.i*2)*self.jacPotTime*2 + (500 - self.i*2)*self.potTime*2)) + 's', end='')
        #_approx_fprime_helper(x, potential.potential_c, 1.4901161193847656e-8, args=(self.centers, self.topol, self.centerNeighbours, self.centerNumNeighbours, self.boxVectors))
        #print(potential.potential_c_jac(x))
        self.start = timer()

    def optimizer(self, coords, centers):
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
            raise ValueError("Internal error: some solvate molecules were\
                              not assigned to a center.")        

        # Drops duplicates
        centerCoordination = np.array(self.topol.bondgraph.\
                             loc[self.topol.bondgraph['j'].\
                             isin(centers)]['j'].value_counts())
        missing = len(centers) - len(centerCoordination)

        # Readd duplicates
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
            #centerNeighbours = np.vstack((centerNeighbours,\
            #                        self.topol.trj.xyz[0][neighbours[0]]*10))

            centerArray = np.repeat(center, len(neighbours[0]))
            dispArray = np.vstack((neighbours, centerArray))
            dispArray = np.transpose(dispArray)
            dispVectors = md.compute_displacements(self.topol.trj, dispArray, periodic = True)[0]
            centerNeighbours = np.vstack((centerNeighbours, -dispVectors*10 + self.topol.trj.xyz[0][center]*10))


        #temp
        self.centers = centers
        self.centerNeighbours = centerNeighbours
        self.centerNumNeighbours = centerNumNeighbours
        #Time objective function
        start = timer()
        potential.potential_c(coords.flatten(), centers, self.topol,\
                              centerNeighbours, centerNumNeighbours, self.boxVectors)
        end = timer()
        self.potTime = end-start
        self.verboseprint("Potential takes: " + str(end-start) +\
                          " seconds to calculate")

        start = timer()
        potential.potential_c_jac(coords.flatten(), centers, self.topol,\
                                  centerNeighbours, centerNumNeighbours, self.boxVectors)
        end = timer()
        self.jacPotTime = end-start
        self.verboseprint("Potential Jacobian takes: " + str(end-start) +\
              " seconds to calculate\n")

        self.verboseprint("Minimizing potential with " +\
                          str(len(coords.flatten())) + " free variables...")

        # Run minimization
        self.start = timer()
        res = minimize(potential.potential_c, coords,
                        args = (centers, self.topol, centerNeighbours,\
                                centerNumNeighbours, self.boxVectors),
                        jac = potential.potential_c_jac,
                        method = 'SLSQP',
                        callback = self.show_progress,
                        #options={'disp': False, 'gtol': 1e-5, 'iprint': 0,\
                        #         'eps': 1.4901161193847656e-5,\
                        #         'maxiter': 5000})
                        options={'disp': True, 'ftol':1e-10, 'iprint': 1,\
                                 'eps': 1.4901161193847656e-8, 'maxiter':500})
        print(res)
        if(res.success):
            self.verboseprint("\n")
            if(not self.silent):
                self.verboseputs.success("Successfully minimized potential!\n")
        else:
            print("\nOptimization failed...\n")

        coords = np.reshape(res.x, (-1, 3))# Since minimization returns flat
                                           # array we need to reshape

        # Recalculate directional vectors
        i = 0
        for coord in coords:
            vectors = np.vstack((vectors, (coord-M[i])/\
                                 np.linalg.norm(coord-M[i])))
            i += 1

        return coords, vectors


    def get_center_neighbours(self, coordination, center):
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

        #try:
        centerIndices = self.topol.extract(center, environment =
                                           {'O': coordination})\
                                         .filter("i").squeeze()

        neighbourgraph = self.topol.bondgraph.loc[self.topol.\
                bondgraph['j'].isin(centerIndices)].\
                sort_values(['j'])[['i', 'j']].pivot(columns= 'j').\
                apply(lambda x: pd.Series(x.dropna().values)).\
                apply(np.int64)['i']

        return (centerIndices, neighbourgraph)

        #except IndexError:  #If no atoms found, extract() throws IndexError
        #    return [], []


    def calculate_pair_vectors(self, coordination, O_frac,\
                               OH_frac, OH2_frac, center):
        """Calculate coordinates for solvate on for 4-coordinated atoms

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

        centerIndices, neighbourgraph = self.get_center_neighbours(\
                                            coordination, center)
        
        if(len(centerIndices) > 0):
            centers = np.empty([0], dtype=int)
            self.verboseprint(("Found {} centers that are Nmax - 2 = {}-fold "
                               "coordinated\n".format(len(centerIndices),coordination)))


            vectors = np.empty([0, 3], dtype = float)
            coordinates = np.empty([0, 3], dtype = float)

            randIndices = random.sample(range(0, len(centerIndices)),
                                        int((OH_frac+OH2_frac)*\
                                        float(len(centerIndices))))
            indices = centerIndices[randIndices]

            for center in indices:
                points = np.empty([0, 3], dtype=float)
                tempVectors = np.empty([0, 3], dtype=float)
                pairVectors = np.empty([0, 3], dtype=float)
                sumVec = np.array([0, 0, 0], dtype=float)
                points = np.vstack((points, self.topol.trj.xyz[0][center]*10))

                neighbours = np.array(neighbourgraph[center])
                centerArray = np.repeat(center, len(neighbours))
                dispArray = np.vstack((neighbours, centerArray))
                dispArray = np.transpose(dispArray)

                dispVectors = md.compute_displacements(self.topol.trj, dispArray, periodic = True)[0]

                for vector in dispVectors:
                    if(np.linalg.norm(vector) < 0.1):
                        puts.warning("Directional vector less than 0.1,"+\
                              " will not hydrate atom with index: " +\
                              str(center) + ". Please verify structure"+\
                              " before using it.\n")
                        continue
                    vector = vector/np.linalg.norm(vector)
                    tempVectors = np.vstack((tempVectors, vector))
                    sumVec += vector # Sum M-O vectors

                sumVec = sumVec/np.linalg.norm(sumVec)

                for vector in tempVectors:
                    tempVector = sumVec - vector
                    pairVectors = np.vstack((pairVectors, tempVector))

                #Sort vectors with increasing norm
                pairVectors = sorted(pairVectors, 
                                    key = lambda x: np.sqrt(x[0]**2 + x[1]**2\
                                                            + x[2]**2))

                # Choose the two vectors with smallest norm
                pairVec1 = pairVectors[0]/np.linalg.norm(pairVectors[0])
                pairVec2 = pairVectors[1]/np.linalg.norm(pairVectors[1])

                # If there are only 2 neighbours, rotate Pi/2 around
                # directional vector to maximize distance to neighbours
                if(len(tempVectors) == 2):
                    axis = sumVec
                    angle = np.pi/2
                    pairVec1 = np.dot(mac.rotate_around_vec(axis, angle), pairVec1)
                    pairVec2 = np.dot(mac.rotate_around_vec(axis, angle), pairVec2)

                # Rotate the vectors towards each other
                # (away from center of bulk)
                crossProd = np.cross(pairVec1, pairVec2)
                dotProdVec1 = np.dot(sumVec, pairVec1)
                dotProdVec2 = np.dot(sumVec, pairVec2)
                if(dotProdVec1 < 0):
                    pairVec1 = np.dot(mac.rotate_around_vec(crossProd, np.pi/7), pairVec1)
                    pairVec2 = np.dot(mac.rotate_around_vec(crossProd, -np.pi/7), pairVec2)
                    #pairVec1 = Quaternion(axis = crossProd,angle = -np.pi/7).\
                    #    rotate(pairVec1)
                    #pairVec2 = Quaternion(axis = crossProd,angle = np.pi/7).\
                    #    rotate(pairVec2)
                else:
                    angle = -np.pi/7
                    pairVec1 = np.dot(mac.rotate_around_vec(crossProd, np.pi/7), pairVec1)
                    pairVec2 = np.dot(mac.rotate_around_vec(crossProd, -np.pi/7), pairVec2)
                    #pairVec2 = Quaternion(axis = crossProd,angle = -np.pi/7).\
                    #    rotate(pairVec2)
                    #pairVec1 = Quaternion(axis = crossProd,angle = np.pi/7).\
                    #    rotate(pairVec1)

                # Calculate coordinates
                coord1 = self.topol.trj.xyz[0][center]*10 +\
                         pairVec1*self.MOBondlength
                coord2 = self.topol.trj.xyz[0][center]*10 +\
                         pairVec2*self.MOBondlength

                # Save relevant coordinates
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


    def calculate_vectors(self, coordination, O_frac,\
                          OH_frac, OH2_frac, center):
        """ Calculate coordinates for solvate on for 5-coordinated atoms

            Notes
            -----
            See calculate_pair_vectors
        """

        vectors = np.empty([0, 3], dtype=float)
        coords = np.empty([0, 3], dtype=float)
        
        vec = np.array([0,0,0])
        tempNeighbour = np.array([0,0,0])
        #Get indices for metal centers with coordination Nmax - 1
        centerIndices, neighbourgraph = self.get_center_neighbours(\
                                            coordination, center)

        if(len(centerIndices) > 0): # If there are any Nmax-1 centers
            centers = np.empty([0], dtype=int)

            self.verboseprint(("Found {} centers that are Nmax - 1 = {}-fold "
                               "coordinated.\n".format(len(centerIndices), coordination)))

            #Calculate only for user specified fraction
            randIndices = random.sample(range(0, len(centerIndices)),\
                                        int((OH2_frac + OH_frac)*\
                                            float(len(centerIndices))))
            indices = centerIndices[randIndices]

            #Calculate M-O vectors
            for center in indices:
                vec = [0, 0, 0]

                neighbours = np.array(neighbourgraph[center])
                centerArray = np.repeat(center, len(neighbours))
                dispArray = np.vstack((neighbours, centerArray))
                dispArray = np.transpose(dispArray)

                dispVectors = md.compute_displacements(self.topol.trj, dispArray, periodic = True)[0]

                #for neighbour in neighbourgraph[center]:
                for vector in dispVectors:
                     vector = vector/np.linalg.norm(vector)
                     vec = vec + vector

                #If norm of directional vector too small
                if(np.linalg.norm(vec) < 0.1):
                    puts.warning("Directional vector almost 0, will " +\
                          "not hydrate at atom with index: " +\
                           str(center + 1))
                    print("Probably due to symmetric center...\n")

                else:
                    centers = np.append(centers, center)
                    vec = vec/np.linalg.norm(vec)
                    vectors = np.vstack((vectors, vec))
                    vec = vec*self.MOBondlength

                    #Save coordinates and correct units
                    coord = self.topol.trj.xyz[0][center]*10 + vec
                    coords = np.vstack((coords, coord))

            return vectors, coords, centers

        else: # Else return empty arrays
            return (np.empty([0, 3], dtype=float),
                    np.empty([0, 3], dtype=float),
                    np.empty([0], dtype=int))


    def solvate(self,params):
        """Calculates coordinates of solvate for specified metal center type

        Parameters
        ----------
        params : dict with keys
        element : string
            Element symbol
        coordination : int
            coordination number of element
        OH_frac : float
            fraction of element atoms that should be hydrated with
            hydroxyl.
        OH2_frac : float
            fraction of element atoms that should be hydrated with
            water.
        Nmax : int
            Maximum coordination number

        Returns
        -------
        vectors : ndarray(float)
            Directional vectors for each solvate
        coord : ndarray(float)
            Coordinates for oxygen atom in each solvate
        centers : ndarray(int)
            List of each center that each solvate belongs to
        """

        element = params['element']
        coordination = params['coordination']
        OH_frac = params['OH']
        OH2_frac = params['OH2']
        O_frac = params['O']
        Nmax = params['Nmax']
 
        if(coordination == Nmax - 1):
            vectors, coords, centers = self.calculate_vectors(\
                                           coordination,\
                                           O_frac, OH_frac, OH2_frac,
                                           element)
            self.hydCoords = np.vstack((self.hydCoords,\
                                        coords[:int(OH_frac*len(vectors))]))
            self.hydVectors = np.vstack((self.hydVectors,\
                                         vectors[:int(OH_frac*len(vectors))]))
            self.hydCenters = np.hstack((self.hydCenters,\
                                         centers[:int(OH_frac*len(vectors))]))
            self.watCoords = np.vstack((self.watCoords,\
                                        coords[int(OH_frac*len(vectors)):]))
            self.watVectors = np.vstack((self.watVectors,\
                                         vectors[int(OH_frac*len(vectors)):]))
            self.watCenters = np.hstack((self.watCenters,\
                                         centers[int(OH_frac*len(vectors)):]))

        elif(coordination == Nmax - 2):
            vectors, coords, centers = self.calculate_pair_vectors(\
                                           coordination,\
                                           O_frac, OH_frac, OH2_frac,
                                           element)
            if(np.isnan(coords).any()):
                raise ValueError("Some coordinates are NaN, aborting....")

            randIndices = random.sample(range(0, len(coords)),\
                                        int((OH_frac)*\
                                            float(len(coords))))

            self.hydCoords = np.vstack((self.hydCoords, coords[randIndices]))
            self.hydVectors = np.vstack((self.hydVectors, vectors[randIndices]))
            self.hydCenters = np.hstack((self.hydCenters, centers[randIndices]))
            mask = np.ones(len(coords), np.bool)
            mask[randIndices] = 0
            self.watCoords = np.vstack((self.watCoords, coords[mask]))
            self.watVectors = np.vstack((self.watVectors, vectors[mask]))
            self.watCenters = np.hstack((self.watCenters, centers[mask]))

            # Take every other element as hydroxyl and the others as water since
            # each 4-coordinated metal site should have one water and one 
            # hydroxyl
            # self.hydCoords = np.vstack((self.hydCoords, coords[::2]))
            # self.hydVectors = np.vstack((self.hydVectors, vectors[::2]))
            # self.hydCenters = np.hstack((self.hydCenters, centers[::2]))
            # self.watCoords = np.vstack((self.watCoords, coords[1::2]))
            # self.watVectors = np.vstack((self.watVectors, vectors[1::2]))
            # self.watCenters = np.hstack((self.watCenters, centers[1::2]))

        else:
            raise ValueError('Can only hydrate Nmax - 1 and Nmax - 2 centers.\
                             You tried to hydrate ' + str(Nmax) + ' - ' +\
                             str(coordination-Nmax) + ' centers. To solve this\
                                                       issue edit config.wet.')


    def optimize(self):
        """Maximize distance between solvate and it's neighbours

        Notes
        -------
        Collects all coordinates and calls optimizer()
        """

        vectors = np.empty([0, 3], dtype=float)
        coords = np.empty([0, 3], dtype=float)
        centers= np.empty([0], dtype=int)

        #Stack coordinates and centers and call optimize()
        coords = np.vstack((coords, self.hydCoords))
        coords = np.vstack((coords, self.watCoords))

        centers = np.hstack((centers, self.hydCenters.astype(int)))
        centers = np.hstack((centers, self.watCenters.astype(int)))

        coords, vectors = self.optimizer(coords, centers)

        #Set attributes to optimized values
        self.hydCoords = coords[:len(self.hydCoords)]
        self.hydVectors = vectors[:len(self.hydCoords)]
        self.watCoords = coords[len(self.hydCoords):]
        self.watVectors = vectors[len(self.hydCoords):]


    def wet(self):
        """Creates atoms at the calculated coordinates

        Returns
        -------
        coords : ndarray
            2D array containing coordinates for all atoms in solvate
        elements : array
            Array containing the element symbols
        """

        hydCoords, hydElements = mac.add_hydroxyl(self.hydCoords, 
                                                  self.hydVectors, 
                                                  self.theta)
        
        watCoords, watElements = mac.add_water(self.watCoords, 
                                               self.watVectors, 
                                               self.theta)

        coords = np.concatenate((watCoords, hydCoords))
        elements = np.concatenate((watElements, hydElements))
        self.coords = coords
        self.elements = elements
        self.verboseprint("Generating output...")

        minODist, minHDist, minStructDist = shortest_distance(coords, elements, self.topol.trj.xyz[0]*10)
        self.verboseprint("Shortest O-O distance in solvent: " + str(minODist) + " Angstrom.")
        self.verboseprint("Shortest H-H distance in solvent: " + str(minHDist) + " Angstrom.")
        self.verboseprint("Shortest distance between solvent and structure: " + str(minStructDist) + " Angstrom.\n")

    def remove_low_coordinated(self, Nmax, element):
        self.topol = remove_low_coordinated(self.topol, Nmax, element, self.silent)

    def save(self, fileWet, resname = 'SOL'):
        self.topol.trj.save(fileWet, force_overwrite = True)
        append_atoms(file = fileWet, coords = self.coords, elements = self.elements, resname = resname)