##################################    .pdb format     ##########################################
#
#		ATOM      1 Ti   TiO A   1       0.000   0.000  36.728  1.00  0.00          Ti"
#
#
#			 1 -  6        Record name   "ATOM  "
#			 7 - 11        Integer       serial       Atom  serial number.
#			13 - 16        Atom          name         Atom name.
#			17             Character     altLoc       Alternate location indicator.
#			18 - 20        Residue name  resName      Residue name.
#			22             Character     chainID      Chain identifier.
#			23 - 26        Integer       resSeq       Residue sequence number.
#			27             AChar         iCode        Code for insertion of residues.
#			31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
#			39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
#			47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
#			55 - 60        Real(6.2)     occupancy    Occupancy.
#			61 - 66        Real(6.2)     tempFactor   Temperature  factor.
#			77 - 78        LString(2)    element      Element symbol, right-justified.
#			79 - 80        LString(2)    charge       Charge  on the atom.
#
#################################################################################################

def append_atoms(file, atomList = [], element='H', coords=[[0,0,0],[1,1,1]], elements = []):
	float_format = lambda x: "%.3f" % x
	f = open(file, "r")
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

			elif(i + len(str(float_format(coord[0]))) == 38):
				tempString += float_format(coord[0])
				i += len(str(float_format(coord[0])))

			elif(i + len(str(float_format(coord[1]))) == 46):
				tempString += float_format(coord[1])
				i += len(str(float_format(coord[1])))

			elif(i + len(str(float_format(coord[2]))) == 54):
				tempString += float_format(coord[2])
				i += len(str(float_format(coord[2])))

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
	
	print("Added " + str(len(coords)) + " atoms to pdb file.")