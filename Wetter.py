'''
To-do list: 
. Add from_dataframe to radish
. Check for overlap
. 4-coordinated splits 1 and adsorbs one
. np.trigonometric radians
. Topologizer.extract throws error if no atoms found
. Restructure (Make the code more effective etc)
. Boundary conditions, affects distance between atoms
. Add support different environments
. Add support for different metals
. Readme
. Timer
. Tkinter GUI?

Notes:
Add water to water+hydroxyl fraction then remove one hydrogen from totalAdded - water fraction

'''

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
import sys

class Wetter:

    def __init__(self, file, verbose, theta = 104.5, waterFrac=0.5, hydroxylFrac = 0.5, MOBondlength = 2.2, MEnvironment = 5, HOHBondlength = 1, OHBondlength = 1, centers = ['Ti']):
        #Input parameters
        self.theta = theta
        self.waterFrac = waterFrac
        self.hydroxylFrac = hydroxylFrac
        self.file = file
        self.verbose = verbose
        self.MOBondlength = MOBondlength
        self.HOHBondlength = HOHBondlength
        self.OHBondlength = OHBondlength
        self.theta = theta/360*2*np.pi
        self.centers = centers

        #Prepare input file for processing

        #Use radish (RAD-algorithm) to compute the coordination of each atom
        self.topol = Topologizer.from_coords(file)
        self.topol.topologize()
        #Set Nmax which is the maximum coordination (çoordination in bulk)
        self.Nmax = self.topol.bondgraph['i'].value_counts().max()

        #Format float precision(will remove from here later maybe....)
        self.float_format = lambda x: "%.3f" % x


    #Construct neighbourgraph where columns are metal centers and rows their coordination shell
    #Sort_values unecessary but nice for readability, will remove later....
    def get_center_neighbours(self, coordination):
        try:
            centerIndices = self.topol.extract(self.centers[0], environment = {'O': self.Nmax - coordination}).index.get_level_values(1)
            neighbourgraph = self.topol.bondgraph.loc[self.topol.bondgraph['j'].isin(centerIndices)].sort_values(['j'])[['i', 'j']].pivot(columns= 'j').apply(lambda x: pd.Series(x.dropna().values)).apply(np.int64)['i']
            return (centerIndices, neighbourgraph)

        except IndexError:
            return [], []


    #For Nmax - 2 coordinated centers calculate pair-vectors
    def calculate_pair_vectors(self):
        centerIndices, neighbourgraph = self.get_center_neighbours(2)

        vectors = np.empty([0, 3], dtype=float)
        tempVectors = np.empty([0, 3], dtype=float)
        pairVectors = np.empty([0, 3], dtype=float)
        tempVector = np.empty([3], dtype=float)
        coordinates = np.empty([0, 3], dtype=float)
        coord = np.empty([3], dtype=float)

        randIndices = random.sample(range(0, len(centerIndices)), int((self.waterFrac + self.hydroxylFrac) * float(len(centerIndices))))
        indices = centerIndices[randIndices]

        for center in indices:

            sumVec = np.array([0, 0, 0], dtype=float)

            for neighbour in neighbourgraph[center]:
                vec_x = self.topol.trj.xyz[0][center][0] - self.topol.trj.xyz[0][neighbour][0]
                vec_y = self.topol.trj.xyz[0][center][1] - self.topol.trj.xyz[0][neighbour][1]
                vec_z = self.topol.trj.xyz[0][center][2] - self.topol.trj.xyz[0][neighbour][2]

                mag = np.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
                tempVectors = np.vstack((tempVectors, [vec_x/mag, vec_y/mag, vec_z/mag]))
                sumVec += [vec_x/mag, vec_y/mag, vec_z/mag]

            #Normalize the sum of vectors
            norm = np.sqrt(sumVec[0]**2 + sumVec[1]**2 + sumVec[2]**2)
            sumVec = sumVec/norm

            for vector in tempVectors:
                
                tempVector = sumVec - vector
                pairVectors = np.vstack((pairVectors, tempVector))

            #Sort vectors with increasing norm
            pairVectors = sorted(pairVectors, key=lambda x: np.sqrt(x[0]**2 + x[1]**2 + x[2]**2))

            #Normalize the two vectors with largest norm that are to be saved
            norm = np.sqrt(pairVectors[2][0]**2+pairVectors[2][1]**2+pairVectors[2][2]**2)
            pairVectors[2] = pairVectors[2]/norm
            norm = np.sqrt(pairVectors[3][0]**2+pairVectors[3][1]**2+pairVectors[3][2]**2)
            pairVectors[3] = pairVectors[3]/norm

            #Save relevant vectors
            vectors = np.vstack((vectors, pairVectors[2]))
            vectors = np.vstack((vectors, pairVectors[3]))

            #Calculate and save coordinates of where oxygen should be placed
            coord = [self.topol.trj.xyz[0][center][0]*10 + pairVectors[2][0]*self.MOBondlength, 
                    self.topol.trj.xyz[0][center][1]*10 + pairVectors[2][1]*self.MOBondlength,
                    self.topol.trj.xyz[0][center][2]*10 + pairVectors[2][2]*self.MOBondlength]
            coordinates = np.vstack((coordinates, coord))

            coord = [self.topol.trj.xyz[0][center][0]*10 + pairVectors[3][0]*self.MOBondlength, 
                    self.topol.trj.xyz[0][center][1]*10 + pairVectors[3][1]*self.MOBondlength,
                    self.topol.trj.xyz[0][center][2]*10 + pairVectors[3][2]*self.MOBondlength]
            coordinates = np.vstack((coordinates, coord))

            tempVectors = np.empty([0, 3], dtype=float)
            pairVectors = pairVectors = np.empty([0, 3], dtype=float)

        return vectors, coordinates

    #Calculate M-O vectors which determine the positions of hydration molecules
    def calculate_vectors(self):

        vectors = np.empty([0, 3], dtype=float)
        coords = np.empty([0, 3], dtype=float)
        vec = np.array([0,0,0])

        #Get indices for metal centers with coordination Nmax - 1
        centerIndices, neighbourgraph = self.get_center_neighbours(1)

        #Calculate only for user specified fraction
        randIndices = random.sample(range(0, len(centerIndices)), int((self.waterFrac + self.hydroxylFrac) * float(len(centerIndices))))
        indices = centerIndices[randIndices]

        #Calculate M-O vectors
        for center in tqdm(indices, ascii=True, desc='Calculating directional vectors'):
            vec = [0, 0, 0]
            for neighbour in neighbourgraph[center]:

                #M-O Vector components
                vec_x = self.topol.trj.xyz[0][center][0] - self.topol.trj.xyz[0][neighbour][0]
                vec_y = self.topol.trj.xyz[0][center][1] - self.topol.trj.xyz[0][neighbour][1]
                vec_z = self.topol.trj.xyz[0][center][2] - self.topol.trj.xyz[0][neighbour][2]

                #M-O vector
                tempVec = np.array([vec_x, vec_y, vec_z])

                #normalize
                mag = np.sqrt(tempVec[0]**2 + tempVec[1]**2 + tempVec[2]**2)
                tempVec = tempVec/mag

                #Sum vectors
                vec = [vec[0] + tempVec[0], vec[1] + tempVec[1], vec[2] + tempVec[2]]

            #normalize
            mag = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
            vec = [vec[0]/mag, vec[1]/mag, vec[2]/mag]
            vectors = np.vstack((vectors, vec))

            #Set vector length to get coordinate
            vec = [vec[0]*self.MOBondlength, vec[1]*self.MOBondlength, vec[2]*self.MOBondlength]

            #Coordinates where oxygen should be placed: vec + coordinate of metal atom multiplied with bondlength
            coord = np.array([self.topol.trj.xyz[0][center][0]*10 + vec[0],
                                self.topol.trj.xyz[0][center][1]*10 + vec[1],
                                self.topol.trj.xyz[0][center][2]*10 + vec[2]])

            #Save coords and correct units
            coords = np.vstack((coords, coord))

        return vectors, coords


    def wet(self):
        #A sum of all vectors pointing from each neighbour of each center to the center itself
        #get coordinates where oxygen should be placed

        vectors, coords = self.calculate_vectors()

        #vectors, coords = self.calculate_pair_vectors()
        
        #vectors = np.concatenate((vectors, pairVectors), axis=0)
        #coords = np.concatenate((coords, pairCoords), axis=0)
        #Check for overlap

        #print whole number and not 1.234 e+4 etc, will remove later....
        np.set_printoptions(suppress=True)
        
        randIndices = random.sample(range(0, len(coords)), int(self.hydroxylFrac * float(len(coords))))

        #Add hydroxide to user specified fraction
        hydCoords = coords[randIndices]
        hydVectors = vectors[randIndices]
        #hydCoords, hydElements = self.add_hydroxyl(hydCoords, hydVectors)
        hydCoords, hydElements = AtomCreator.add_hydroxyl(hydCoords, hydVectors, self.theta)
        print("Adding: " + str(len(hydCoords)) + " atoms (hydroxyl)")

        #Create mask for the selection of water molecules
        mask = np.ones(len(coords), np.bool)
        mask[randIndices] = 0

        #Create water molecules at the inverse selection of where hydroxide was placed
        watCoords = coords[mask]
        watVectors = vectors[mask]
        #watCoords, watElements = self.add_water(watCoords, watVectors)
        watCoords, watElements = AtomCreator.add_water(watCoords, watVectors, self.theta)
        print("Adding: " + str(len(watCoords)) + " atoms (water)")

        #Concatenate water and hydroxide coordinates
        atoms = np.concatenate((watCoords, hydCoords))
        elements = np.concatenate((watElements, hydElements))

        append_atoms(file = self.file, coords = atoms, elements = elements)

        #self.add_water(pairCoords, pairVectors)
        #embed()