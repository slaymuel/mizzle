'''
To-do list: 
1. correctly structure lines to be appended to pdb-file (index of atom...)
2. Add water molecules (calculate coordinates, angles etc)
3. Restructure (Make the code more effective etc)
4. Add support different environments
5. Add support for different metals

. remove low (N < Nmax-1) coordinated atoms

'''
import numpy as np
import pandas as pd
from radish import Topologizer
from tqdm import tqdm

class Wetter:

	def __init__(self, file, verbose, theta = 104.5, frac=0.4, MOBondlength = 2, MEnvironment = 5):
		#Input parameters
		self.theta = theta
		self.frac = frac
		self.file = file
		self.verbose = verbose
		self.MOBondlength = MOBondlength

		#Prepare input file for processing
		self.topol = Topologizer.from_coords(file)
		self.topol.topologize()
		self.ind = self.topol.extract('Ti', environment = {'O': MEnvironment}).index.get_level_values(1)

		#Construct neighbourgraph where columns are metal centers and rows their coordination shell
		#Sort_values unecessary but nice for readability, will remove later....
		self.neighbourgraph = self.topol.bondgraph.loc[self.topol.bondgraph['j'].isin(self.ind)].sort_values(['j'])[['i', 'j']].pivot(columns= 'j').apply(lambda x: pd.Series(x.dropna().values)).apply(np.int64)['i']

		#Format float precision
		self.float_format = lambda x: "%.3f" % x


		#Not used.......
	def createWater(self, coord, r):
		O = [coord[0], coord[1], coord[2]]
		H1 = [r*np.sin(self.theta/2), 0, r*np.cos(self.theta/2)]
		H2 = [r*np.sin(-self.theta/2), 0, r*np.cos(-self.theta/2)]

		#Append atoms to .pdb file
	def appendAtoms(self, atomList = [], element='O', coords=[[0,0,0],[1,1,1]]):
		f = open(self.file, "r")
		content = f.readlines()
		f.close()

		#get indices for all entries in list 'content' where the substring 'ATOM' occurs
		indices = [index for index, line in enumerate(content) if 'ATOM' in line]

		#get number of atoms
		nrAtoms = self.topol.trj.top.n_atoms

		#Prepare list of atoms to be appended to pdb file
		for coord in coords:
			atomList.append("ATOM   "+str(nrAtoms+1)+"  "+element+"   TiO A   1      "+self.float_format(coord[0])+"   "+self.float_format(coord[1])+"   "+self.float_format(coord[2])+"   1.00  0.00          "+element)
			nrAtoms += 1
	   	#"ATOM      1 Ti   TiO A   1       0.000   0.000  36.728  1.00  0.00          Ti"

	   	if(content[indices[-1] + 1][:3] == 'TER'):
	   		print 'found'

	   	new_content = content[:indices[-1] + 1]
	   	
	   	#Append new atoms
	   	new_content.extend(atomList)
	   	new_content.extend(content[indices[-1] + 1:])

	   	#Print to file

	   	#for line in new_content:
	   	#	print line

	#Calculate M-O vectors
	def calculateVectors(self):
		if(self.verbose):
			print("Ind:")
			print self.ind
			print("neighbourgraph:")
			print self.neighbourgraph

		#Initialize array to hold vectors
		vectors = np.empty([0, 4], dtype=float)	#change initializer....
		coords = np.array([[-1, -1, -1]])

		#Calculate M-O vectors
		for center in tqdm(self.ind, ascii=True, desc='Progress'):
			#Only calculate for the fraction that will actually get hydrated
			if(np.random.rand(1) < self.frac):
				for neighbour in self.neighbourgraph[center]:
					#M-O Vector components
					vec_x = self.topol.trj.xyz[0][center][0] - self.topol.trj.xyz[0][neighbour][0]
					vec_y = self.topol.trj.xyz[0][center][1] - self.topol.trj.xyz[0][neighbour][1]
					vec_z = self.topol.trj.xyz[0][center][2] - self.topol.trj.xyz[0][neighbour][2]

					#normalize. also save center to keep track of which atom the vector belongs to
					mag = np.sqrt(vec_x**2 + vec_y**2 + vec_z**2)
					vec = np.array([center, vec_x/mag, vec_y/mag, vec_z/mag])

					#Coordinates where oxygen should be placed: vec + coordinate of metal atom
					coord = np.array([self.topol.trj.xyz[0][center][0] + vec[1],
										self.topol.trj.xyz[0][center][1] + vec[2],
										self.topol.trj.xyz[0][center][2] + vec[3]])
					
					#Save coords and correct units
					coords = np.vstack((coords, coord*10))
					vectors = np.vstack((vectors, vec))

		return vectors, coords


	def wet(self):
		#get coordinates where oxygen should be placed
		vectors, coords = self.calculateVectors()
		i = 0

		#print whole number and not 1.234 e+4 etc, will remove later....
		np.set_printoptions(suppress=True)
		self.appendAtoms(coords = coords)

		if(self.verbose):
			for vector in vectors:
				print("Vector: " + str(i))
				print vector
				i = i + 1
			i=0
			
			#"ATOM      1 Ti   TiO A   1       0.000   0.000  36.728  1.00  0.00          Ti"
			#"ATOM      1 Ti   TiO A   1       float_format(coord[0])   float_format(coord[1])  float_format(coord[2])  1.00  0.00          Ti"
			#self.appendAtoms(atomList)
			for coord in coords:
				#print("coord: " + str(i))
				print (coord)
				i = i + 1


