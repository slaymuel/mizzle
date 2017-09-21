'''
To-do list: 
. Add water molecules (calculate coordinates, angles etc)
. Edge cases, directional vector = (0, 0, 1)
. Restructure (Make the code more effective etc)
. Add support different environments
. Add support for different metals

. remove low (N < Nmax-1) coordinated atoms

'''

# pymatgen
# ase
import numpy as np
import pandas as pd
from radish import Topologizer
from tqdm import tqdm

class Wetter:

	def __init__(self, file, verbose, theta = 104.5, frac=0.4, MOBondlength = 2, MEnvironment = 5, HOHBondlength = 1, OHBondlength = 1):
		#Input parameters
		self.theta = theta
		self.frac = frac
		self.file = file
		self.verbose = verbose
		self.MOBondlength = MOBondlength
		self.HOHBondlength = HOHBondlength
		self.OHBondlength = OHBondlength

		#Prepare input file for processing
		self.topol = Topologizer.from_coords(file)
		self.topol.topologize()
		self.ind = self.topol.extract('Ti', environment = {'O': MEnvironment}).index.get_level_values(1)

		#Construct neighbourgraph where columns are metal centers and rows their coordination shell
		#Sort_values unecessary but nice for readability, will remove later....
		self.neighbourgraph = self.topol.bondgraph.loc[self.topol.bondgraph['j'].isin(self.ind)].sort_values(['j'])[['i', 'j']].pivot(columns= 'j').apply(lambda x: pd.Series(x.dropna().values)).apply(np.int64)['i']

		#Format float precision
		self.float_format = lambda x: "%.3f" % x

	#Rotates around x-axis
	def xRotate(self, vector, angle):
		rotMatrix = [[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]
		dotProd = np.dot(rotMatrix, vector)
		return (dotProd)

	#Calculate rotation matrix by rotating z-axis (0, 0, 1) to align with MO-vector
	def align(self, vec1, vec2):
		I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		crossProd = np.cross(vec1, vec2)	#sine
		dotProd = np.dot(vec1, vec2)	#cosine
		print("crossprod: " + str(crossProd))
		print("dotprod: " + str(dotProd))

		#Edge cases: If vectors are parallell since dot product of normalized vectors is the cosine of the angle.
		if(dotProd < 1.01 and dotProd > 0.99):
			print("dotprod is 1")
			rotMatrix = -I

		elif(dotProd < -0.99 and dotProd > -1.01):
			print("dotprod is -1")
			#Find orthonormal vector to both (cross product) and rotate pi
			mag = np.sqrt(crossProd[0]**2 + crossProd[1]**2 + crossProd[2]**2)
			ortVec = crossProd/mag

			rotMatrix = -np.array([[-1+2*ortVec[0]**2, 2*ortVec[0]*ortVec[1], 2*ortVec[0]*ortVec[2]],
								[2*ortVec[0]*ortVec[1], -1+2*ortVec[1]**2, 2*ortVec[1]*ortVec[2]],
								[2*ortVec[0]*ortVec[2], 2*ortVec[1]*ortVec[2], -1+2*ortVec[2]**2]])

		else:
			#skew-symmetric: transpose equals negative. Used to represent cross product as matrix multiplication
			skewSym	= np.array([[0, -crossProd[2], crossProd[1]],[crossProd[2], 0, -crossProd[0]],[-crossProd[1], crossProd[0], 0]])
			prod = np.matmul(skewSym, skewSym)*(1/(1+dotProd))
			rotMatrix = np.add(np.add(I, skewSym), prod)

		return rotMatrix

	def createWater(self, coords, vectors):
		O = [0, 0, 0]
		H1 = [np.sin(self.theta/2), 0, np.cos(self.theta/2)]
		H2 = [-np.sin(self.theta/2), 0, np.cos(self.theta/2)]

		H1 = self.xRotate(H1, np.pi/5)
		H2 = self.xRotate(H2, np.pi/5)
		#O = self.xRotate(O, np.pi/5)
		i = 49
		rotMatrix = self.align([0, 0, 1], [vectors[i][0], vectors[i][1], vectors[i][2]])
		print("Rotmatrix: \n" + str(rotMatrix))
		O = np.dot(rotMatrix, O)
		H1 = np.dot(rotMatrix, H1)
		H2 = np.dot(rotMatrix, H2)

		#Translate to correct coordinates
		transVector = [coords[i][0] - O[0], coords[i][1] - O[1], coords[i][2] - O[2]]
		O = [O[0] + transVector[0], O[1] + transVector[1], O[2] + transVector[2]]
		H1 = [H1[0] + transVector[0], H1[1] + transVector[1], H1[2] + transVector[2]]
		H2 = [H2[0] + transVector[0], H2[1] + transVector[1], H2[2] + transVector[2]]

		self.appendAtoms(coords = [O, H1, H2])

		#Append atoms to .pdb file
	def appendAtoms(self, atomList = [], element='H', coords=[[0,0,0],[1,1,1]], elements = []):
		f = open(self.file, "r")
		content = f.readlines()
		f.close()

		#get indices for all entries in list 'content' where the substring 'ATOM' occurs
		indices = [index for index, line in enumerate(content) if 'ATOM' in line]

		#get number of atoms
		nrAtoms = len(indices)#self.topol.trj.top.n_atoms

		#Variables for the .pdb format
		residueName = 'TiO'
		chainIdentifier = 'A'
		residueSequenceNumber = '1'
		occupancy = '1.00'
		temperatureFactor = '0.00'

		#Prepare list of atoms to be appended to pdb file
		for coord in coords:
			
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

	   	#Get old content of file until last atom entry
	   	new_content = content[:indices[-1] + 1]
	   	#Append new atoms
	   	new_content.extend(atomList)
	   	#append lines after the final atom entry
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

		print("Added " + str(len(coords)) + " atoms")
		return vectors, coords


	def wet(self):
		#get coordinates where oxygen should be placed
		vectors, coords = self.calculateVectors()
		i = 0

		#print whole number and not 1.234 e+4 etc, will remove later....
		np.set_printoptions(suppress=True)

		#self.appendAtoms(coords = coords)
		self.createWater(coords, vectors)
		if(self.verbose):
			for vector in vectors:
				print("Vector: " + str(i))
				print vector
				i = i + 1
			i=0
			
			for coord in coords:
				#print("coord: " + str(i))
				print (coord)
				i = i + 1


