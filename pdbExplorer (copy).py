#Handles pdb file

##################################    .pdb format     ############################################
#                                                                                                #
#       ATOM      1 Ti   TiO A   1       0.000   0.000  36.728  1.00  0.00          Ti"          #
#                                                                                                #
#                                                                                                #
#            1 -  6        Record name   "ATOM  "                                                #
#            7 - 11        Integer       serial       Atom  serial number.                       #
#           13 - 16        Atom          name         Atom name.                                 #
#           17             Character     altLoc       Alternate location indicator.              #
#           18 - 20        Residue name  resName      Residue name.                              #
#           22             Character     chainID      Chain identifier.                          #
#           23 - 26        Integer       resSeq       Residue sequence number.                   #
#           27             AChar         iCode        Code for insertion of residues.            #
#           31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms. #
#           39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms. #
#           47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms. #
#           55 - 60        Real(6.2)     occupancy    Occupancy.                                 #
#           61 - 66        Real(6.2)     tempFactor   Temperature  factor.                       #
#           77 - 78        LString(2)    element      Element symbol, right-justified.           #
#           79 - 80        LString(2)    charge       Charge  on the atom.                       #
#                                                                                                #
##################################################################################################
"""Example NumPy style docstrings.

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
from radish import Topologizer
import re
import numpy as np
from shutil import copyfile
import mdtraj as md
import pandas as pd
from timeit import default_timer as timer

#Appends atoms to pdb file
def append_atoms(file, coords=[], elements = []):
    atomList = []
    fileName = file
    float_format = lambda x: "%.3f" % x
    f = open(file, "r")
    content = f.readlines()
    f.close()

    #get indices for all entries in list 'content' where the substring 'ATOM' occurs
    indices = [index for index, line in enumerate(content) if 'ATOM' in line]

    #get number of atoms
    nrAtoms = len(indices)#self.topol.trj.top.n_atoms

    #Variables for the .pdb format
    residueName = 'SOL'
    chainIdentifier = 'A'
    residueSequenceNumber = '1'
    occupancy = '1.00'
    temperatureFactor = '0.00'

    j = 0
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

            elif(i + len(elements[j]) == 14):
                tempString += elements[j]
                i += len(elements[j])

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
                tempString += elements[j]
                i += len(elements[j])

            tempString += " "
            i += 1
        j += 1

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
    file = open(file, 'w')
    for line in new_content:
        file.write("%s\n" % line.rstrip())  #also remove newline characters
    file.close()

    print("Added " + str(len(coords)) + " atoms to " + fileName)

#Remove atoms with coordination less than Nmax - 3 from pdb file
def remove_lower_coordinated(topol, Nmax):

    #topol = Topologizer.from_coords(file)
    topol.topologize()

    while(1):
        indices = []

        #Create inputs for new Topologizer instance
        xyz = topol.trj.xyz[0]
        topology = topol.trj.topology.to_dataframe()

        try:
            i = 3
            while(i <= Nmax):
                centerIndices = topol.extract('Ti', environment = {'O': Nmax - i}).index.get_level_values(1)
                indices.extend(centerIndices)
                i += 1
        except IndexError:
            pass

        try:
#######################################     FIX     ################################################################
            #pass
            oxygenIndices = topol.extract('O', environment = {'Ti': 1}).index.get_level_values(1)
            indices.extend(oxygenIndices)
####################################################################################################################
        except IndexError:
            pass

        if(len(indices) < 1):
            print("All low coordinated atoms removed")
            return topol

        #Remove atoms
        topology = topology[0].drop(topology[0].index[indices])

        #Rewrite indices
        topology.reset_index(drop = True, inplace = True)
        #Rewrite serial column
        topology['serial'] = topology.index + 1

        #Remove atoms
        indices = sorted(indices, reverse = True)   #Sort in decreasing order to allow loop to find correct index
        for index in indices:
            xyz = np.vstack((xyz[:index], xyz[index + 1:]))

        print("Removed " + str(len(indices)) + " low coordinated atoms")
        topology = md.Topology.from_dataframe(topology, bonds = None)
        trajectory = md.Trajectory(xyz, topology)
        topol = Topologizer.from_mdtraj(trajectory)

        topol.topologize()

'''
    content = []
    f = open(file, 'r')
    content = f.readlines()
    f.close()
    print("Length of file: " + str(len(content)))

    #Remove low coordinated center atoms
    while(1):
        indices = []

        try:
            i = 3
            while(i <= Nmax):
                centerIndices = topol.extract('Ti', environment = {'O': Nmax - i}).index.get_level_values(1)
                indices.extend(centerIndices)
                i += 1
        except IndexError:
            #print("All low coordinated Ti removed")
            pass

        try:
#######################################     FIX     ################################################################
            pass
            #oxygenIndices = topol.extract('O', environment = {'Ti': 1}).index.get_level_values(1)
            #indices.extend(oxygenIndices)
####################################################################################################################
        except IndexError:
            print("All low coordinated oxygen removed")
            pass

        if(len(indices) < 1):
            return file

        #Create regular expressions to find undercoordinated atoms by their index
        regexes = []
        atomRegex = re.compile("ATOM")

        for index in indices:
            regexes.append(re.compile("ATOM\s+" + str(index+1) + " "))

        print("Atoms to be removed: " + str(len(regexes)))
        newContent = []

        i = 0
        j = 0
        for line in content:
            if(atomRegex.match(line)):
                j += 1
            if any(regex.match(line) for regex in regexes):
                i += 1
            else:
                if(atomRegex.match(line)):
                    k = len(str(j-i))
                    tempStr = ""
                    while(k < len(str(j))):
                        tempStr += " "
                        k += 1
                    tempStr += str(j-i)
                    newLine = line.replace(str(j), tempStr, 1)
                else:
                    newLine = line
                newContent.append(newLine)
            #rewrite indices

        print("Removed " + str(i))

        if(i != len(indices)):
            print("Error!")

        #file = fileWet
        #file = 'test.pdb'

        with open(file, 'w') as f:
            for line in newContent:
                f.write("%s\n" % line.rstrip()) #also remove newline characters

        print("Length of file after removal: " + str(len(newContent)))


        content = newContent


        topol = Topologizer.from_coords(file)
        topol.topologize()
'''
    #return file