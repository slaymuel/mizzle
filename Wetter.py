'''
To-do list: 
. Remove one coordinated oxygens?
. Add water molecules (calculate coordinates, angles etc)
. Put 2 molecules on Nmax-2 coordinated centers
. Restructure (Make the code more effective etc)
. Boundary conditions, affects distance between atoms
. Add support different environments
. Add support for different metals
. Tkinter GUI?
. remove low (N < Nmax-1) coordinated atoms

Notes:
Add water to water+hydroxyl fraction then remove one hydrogen from totalAdded - water fraction

'''

# pymatgen
# ase
import numpy as np
import pandas as pd
from pdbExplorer import append_atoms
import random
#import quaternion as quat
from pyquaternion import Quaternion
from radish import Topologizer
from tqdm import tqdm

class Wetter:

	def __init__(self, file, verbose, theta = 104.5, waterFrac=0.5, hydroxylFrac = 0.1, MOBondlength = 2, MEnvironment = 5, HOHBondlength = 1, OHBondlength = 1):
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

		#Prepare input file for processing
		self.topol = Topologizer.from_coords(file)
		self.topol.topologize()
		#Set Nmax
		self.Nmax = self.topol.bondgraph['i'].value_counts().max()
		#Get centers for environment
		self.centerIndices = self.topol.extract('Ti', environment = {'O': self.Nmax - 1}).index.get_level_values(1)
		#Construct neighbourgraph where columns are metal centers and rows their coordination shell
		#Sort_values unecessary but nice for readability, will remove later....
		self.neighbourgraph = self.topol.bondgraph.loc[self.topol.bondgraph['j'].isin(self.centerIndices)].sort_values(['j'])[['i', 'j']].pivot(columns= 'j').apply(lambda x: pd.Series(x.dropna().values)).apply(np.int64)['i']

		#Format float precision(will remove from here later maybe....)
		self.float_format = lambda x: "%.3f" % x
		

	#NOT USED
	def find_under_coordinated(self, atom):
		for neighbour in self.topol.bondgraph.loc[self.topol.bondgraph['j'] == atom, 'i']:
			#print(neighbour)
			if(self.topol.top.loc[neighbour]['element'] == 'Ti'):
				if(self.topol.bondgraph['i'].value_counts().loc[neighbour] - 1 <= self.Nmax - 3):
					print(self.topol.top.loc[neighbour]['element'])
					print("found undercoordinated")
					return self.findUnderCoordinated(neighbour).append(neighbour)

			elif(self.topol.top.loc[neighbour]['element'] == 'O'):
				if(self.topol.bondgraph['i'].value_counts().loc[neighbour] < 2):
					pass
					#Remove this oxygen?
			else:
				return []

		#removeFromPDB(atom)

	#Rotates around x-axis
	def xRotate(self, vector, angle):
		rotMatrix = [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]
		dotProd = np.dot(rotMatrix, vector)
		return (dotProd)

	def random_rotate(self, vector, angle):
		rotMatrix = np.array([[np.cos(angle)+vector[0]**2*(1-np.cos(angle)), vector[0]*vector[1]*(1 - np.cos(angle))-vector[2]*np.sin(angle), vector[0]*vector[2]*(1-np.cos(angle))+vector[1]*np.sin(angle)],
					[vector[1]*vector[0]*(1 - np.cos(angle))+vector[2]*np.sin(angle), np.cos(angle)+vector[1]**2*(1-np.cos(angle)), vector[1]*vector[2]*(1-np.cos(angle))-vector[0]*np.sin(angle)],
					[vector[2]*vector[0]*(1 - np.cos(angle))-vector[1]*np.sin(angle), vector[2]*vector[1]*(1 - np.cos(angle))+vector[0]*np.sin(angle), np.cos(angle)+vector[2]**2*(1-np.cos(angle))]])
		return rotMatrix

	#Calculate rotation matrix by rotating z-axis (0, 0, 1) to align with MO-vector
	def align(self, vec1, vec2):
		I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		crossProd = np.cross(vec1, vec2)	#sine
		dotProd = np.dot(vec1, vec2)	#cosine
		#print("crossprod: " + str(crossProd))
		#print("dotprod: " + str(dotProd))

		#Edge cases: If vectors are parallell since dot product of normalized vectors is the cosine of the angle.
		if(dotProd < 1.01 and dotProd > 0.99):
			#print("dotprod is 1")
			rotMatrix = I

		elif(dotProd < -0.99 and dotProd > -1.01):
			#print("dotprod is -1")
			#Find orthonormal vector to both (cross product) and rotate pi
			#mag = np.sqrt(crossProd[0]**2 + crossProd[1]**2 + crossProd[2]**2)
			#ortVec = crossProd/mag

			#Will always be the same? Are there edge cases.....?
			#rotMatrix = -np.array([[-1+2*ortVec[0]**2, 2*ortVec[0]*ortVec[1], 2*ortVec[0]*ortVec[2]],
			#					[2*ortVec[0]*ortVec[1], -1+2*ortVec[1]**2, 2*ortVec[1]*ortVec[2]],
			#					[2*ortVec[0]*ortVec[2], 2*ortVec[1]*ortVec[2], -1+2*ortVec[2]**2]])
			rotMatrix = -I 		#np.array([[1, 0, 0],[0, -1, 0],[0, 0, -1]])

		#Need to take into account when vectors are orthagonal i.e when dot product is 0?
		else:
			#skew-symmetric: transpose equals negative. Used to represent cross product as matrix multiplication
			skewSym	= np.array([[0, -crossProd[2], crossProd[1]],[crossProd[2], 0, -crossProd[0]],[-crossProd[1], crossProd[0], 0]])
			prod = np.matmul(skewSym, skewSym)*(1/(1+dotProd))
			rotMatrix = np.add(np.add(I, skewSym), prod)

		return rotMatrix

	def add_atoms(self, coords, vectors):
		atoms = np.empty([0, 3], dtype=float)
		elements = []
		i=0

		#Loop over all coordinates where oxygen can be placed
		while i < len(coords):
			O = np.array([0, 0, 0])
			H1 = np.array([np.sin(self.theta/2), 0, np.cos(self.theta/2)])
			H2 = np.array([-np.sin(self.theta/2), 0, np.cos(self.theta/2)])

			#No need to rotate O since it lies on the x-axis
			angle = np.arccos(np.cos(115)/np.cos(104.5/2))
			H1 = self.xRotate(H1, angle)
			H2 = self.xRotate(H2, angle)

			#Align z axis to the directional vector
			rotMatrix = self.align([0, 0, 1], vectors[i])
			O = np.dot(rotMatrix, O)
			H1 = np.dot(rotMatrix, H1)
			H2 = np.dot(rotMatrix, H2)



			#Random rotation along directional vector
			#randRotMatrix = self.randomRot(vectors[i], random.uniform(0.1, 2*np.pi))
			#randRotMatrix = self.randomRot(vectors[i], 1.2)

			#Rotate randomly around the directional vector
			axis = vectors[i]
			theta = random.uniform(0.1, 2*np.pi)
			#O = Quaternion(axis=axis,angle=theta).rotate(O)
			H1 = Quaternion(axis=axis,angle=theta).rotate(H1)
			H2 = Quaternion(axis=axis,angle=theta).rotate(H2)

			#O = randRotMatrix.dot(O)
			#H1 = randRotMatrix.dot(H1)
			#H1 = randRotMatrix.dot(H2)

			#Translate to correct coordinates
			transVector = [coords[i][0] - O[0], coords[i][1] - O[1], coords[i][2] - O[2]]
			O = np.array([O[0] + transVector[0], O[1] + transVector[1], O[2] + transVector[2]])
			H1 = np.array([H1[0] + transVector[0], H1[1] + transVector[1], H1[2] + transVector[2]])
			H2 = np.array([H2[0] + transVector[0], H2[1] + transVector[1], H2[2] + transVector[2]])


			i += 1

			#Save atoms to be added to pdb file
			atoms = np.vstack((atoms, O))
			elements.extend('O')
			atoms = np.vstack((atoms, H1))
			elements.extend('H')
			atoms = np.vstack((atoms, H2))
			elements.extend('H')
		#Append atoms to .pdb file
		append_atoms(file = self.file, coords = atoms, elements = elements)
		#self.appendAtoms(coords = [O, H1, H2], elements=['O', 'H', 'H'])


	#Calculate M-O vectors
	def calculate_vectors(self):
		if(self.verbose):
			print("Ind:")
			print self.ind
			print("neighbourgraph:")
			print self.neighbourgraph

		vectors = np.empty([0, 3], dtype=float)
		#pairVectors = vectors = np.empty([0, 3], dtype=float)
		coords = np.empty([0, 3], dtype=float)
		vec = np.array([0,0,0])
		#pairVec = np.array([0,0,0])

		#Calculate only for user specified fraction
		randIndices = random.sample(range(0, len(self.centerIndices)), int(self.waterFrac * float(len(self.centerIndices))))
		indices = self.centerIndices[randIndices]
		#indices = self.centerIndices

		#Calculate M-O vectors
		print(indices)
		for center in tqdm(indices, ascii=True, desc='Calculating directional vectors'):
			#Only calculate for the fraction that will actually get hydrated (change this later to actually choose that fraction)
			vec = [0, 0, 0]
			for neighbour in self.neighbourgraph[center]:
				# for otherNeighbour in self.neighbourgraph[center]:
				# 	vec_x = self.topol.trj.xyz[0][center][0] - self.topol.trj.xyz[0][neighbour][0]
				# 	vec_y = self.topol.trj.xyz[0][center][1] - self.topol.trj.xyz[0][neighbour][1]
				# 	vec_z = self.topol.trj.xyz[0][center][2] - self.topol.trj.xyz[0][neighbour][2]
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

			#Set vector length
			vec = [vec[0]*self.MOBondlength, vec[1]*self.MOBondlength, vec[2]*self.MOBondlength]

			#Coordinates where oxygen should be placed: vec + coordinate of metal atom multiplied with bondlength
			coord = np.array([self.topol.trj.xyz[0][center][0]*10 + vec[0],
								self.topol.trj.xyz[0][center][1]*10 + vec[1],
								self.topol.trj.xyz[0][center][2]*10 + vec[2]])

			#Save coords and correct units
			coords = np.vstack((coords, coord))

		return vectors, coords


	def wet(self):
		#get coordinates where oxygen should be placed
		#self.topol = pdbExplorer.removeLowerCoordinated(self.topol, self.file)

		vectors, coords = self.calculate_vectors()
		i = 0

		#Check for overlap
		#print whole number and not 1.234 e+4 etc, will remove later....
		np.set_printoptions(suppress=True)

		self.add_atoms(coords, vectors)

		if(self.verbose):
			for vector in vectors:
				print("Vector: " + str(i))
				print vector
				i = i + 1
			i=0
			
			for coord in coords:
				print (coord)
				i = i + 1