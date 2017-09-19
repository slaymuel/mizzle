'''
To-do list: 
1. correctly structure lines to be appended to pdb-file (index of atom...)
2. Add water molecules (calculate coordinates, angles etc)
3. Restructure (Make the code more effective etc)
4. Add support different environments
5. Add support for different metals

. remove low (N < Nmax-1) coordinated atoms

'''

# pymatgen
# ase
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
	def appendAtoms(self, atomList = [], element='H', coords=[[0,0,0],[1,1,1]]):
		f = open(self.file, "r")
		content = f.readlines()
		f.close()

		#get indices for all entries in list 'content' where the substring 'ATOM' occurs
		indices = [index for index, line in enumerate(content) if 'ATOM' in line]

		#get number of atoms
		nrAtoms = self.topol.trj.top.n_atoms

		#Variables for the .pdb format
		residueName = 'TiO'
		chainIdentifier = 'A'
		residueSequenceNumber = '1'
		occupancy = '1.00'
		temperatureFactor = '0.00'

		#Prepare list of atoms to be appended to pdb file
		for coord in coords:
			
			#nrAtoms += 1
			i = 4
			nrAtoms += 1
			tempString = "ATOM"

			#Format tempString to .pdb format
			while(i <= 80):
				if(i + len(str(nrAtoms)) == 11):
					tempString += str(nrAtoms)
					i += len(str(nrAtoms))

				elif(i + len(element) == 14):
					tempString += element
					i += len(element)

				elif(i + len(residueName) == 20):
					tempString += residueName
					i += len(residueName)

				elif(i +len(chainIdentifier) == 22):
					tempString += chainIdentifier
					i += 1

				elif(i + len(residueSequenceNumber) == 26):
					tempString += residueSequenceNumber
					i += len(residueSequenceNumber)

				elif(i + len(str(self.float_format(coord[0]))) == 38):
					tempString += self.float_format(coord[0])
					i += len(str(self.float_format(coord[0])))

				elif(i + len(str(self.float_format(coord[1]))) == 46):
					tempString += self.float_format(coord[1])
					i += len(str(self.float_format(coord[1])))

				elif(i + len(str(self.float_format(coord[2]))) == 54):
					tempString += self.float_format(coord[2])
					i += len(str(self.float_format(coord[2])))

				elif(i + len(occupancy) == 60):
					tempString += occupancy
					i += len(occupancy)

				elif(i + len(temperatureFactor) == 66):
					tempString += temperatureFactor
					i += len(temperatureFactor)

				elif(i == 76):
					tempString += element
					i += len(element)

				tempString += " "
				i += 1

			#Append formatted tempString
			atomList.append(tempString)

##################################    .pdb format     ##########################################
	   	#"ATOM      1 Ti   TiO A   1       0.000   0.000  36.728  1.00  0.00          Ti"

	   	'''
			 1 -  6        Record name   "ATOM  "
			 7 - 11        Integer       serial       Atom  serial number.
			13 - 16        Atom          name         Atom name.
			17             Character     altLoc       Alternate location indicator.
			18 - 20        Residue name  resName      Residue name.
			22             Character     chainID      Chain identifier.
			23 - 26        Integer       resSeq       Residue sequence number.
			27             AChar         iCode        Code for insertion of residues.
			31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
			39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
			47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			55 - 60        Real(6.2)     occupancy    Occupancy.
			61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			77 - 78        LString(2)    element      Element symbol, right-justified.
			79 - 80        LString(2)    charge       Charge  on the atom.
	   	'''
#################################################################################################

		#Don't know what to do with this yet.......
	   	if(content[indices[-1] + 1][:3] == 'TER'):
	   		print 'found'

	   	new_content = content[:indices[-1] + 1]
	   	#Append new atoms
	   	new_content.extend(atomList)
	   	new_content.extend(content[indices[-1] + 1:])

	   	#Print to file
		file = open('test.pdb', 'w')
		for line in new_content:
			file.write("%s\n" % line.rstrip())	#also remove newline characters
		file.close()

	#Calculate M-O vectors
	def calculateVectors(self):
		if(self.verbose):
			print("Ind:")
			print self.ind
			print("neighbourgraph:")
			print self.neighbourgraph

		vectors = np.empty([0, 3], dtype=float)
		coords = np.empty([0, 3], dtype=float)
		vec = np.array([0,0,0])

		#Calculate M-O vectors
		for center in tqdm(self.ind, ascii=True, desc='Progress'):
			#Only calculate for the fraction that will actually get hydrated (change this later to actually choose that fraction)
			vec = [0, 0, 0]
			for neighbour in self.neighbourgraph[center]:
				#M-O Vector components
				vec_x = self.topol.trj.xyz[0][center][0] - self.topol.trj.xyz[0][neighbour][0]
				vec_y = self.topol.trj.xyz[0][center][1] - self.topol.trj.xyz[0][neighbour][1]
				vec_z = self.topol.trj.xyz[0][center][2] - self.topol.trj.xyz[0][neighbour][2]

				#M-O vector
				tempVec = np.array([vec_x, vec_y, vec_z])

				#Sum vectors
				vec = [vec[0] + tempVec[0], vec[1] + tempVec[1], vec[2] + tempVec[2]]

			#normalize
			mag = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
			vec = [vec[0]/mag, vec[1]/mag, vec[2]/mag]

			#Coordinates where oxygen should be placed: vec + coordinate of metal atom multiplied with bondlength
			coord = np.array([self.topol.trj.xyz[0][center][0] + vec[0],
								self.topol.trj.xyz[0][center][1] + vec[1],
								self.topol.trj.xyz[0][center][2] + vec[2]])
			
			#print(coord)

			#Save coords and correct units
			coords = np.vstack((coords, coord*10))
			vectors = np.vstack((vectors, vec))

		print("Added " + str(len(coords)) + " atoms")
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


