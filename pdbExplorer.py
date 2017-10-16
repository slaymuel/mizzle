"""pdbExplorer

This module handles formatting of output pdb-file and preparation of
input pdb-file.

Notes
-----
The append_atoms method appends new atoms to an existing pdb file and
the remove_lower_coordianted method removes low coordinated atoms from a
Topology instance.
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

        #Remove atoms from topology
        topology = topology[0].drop(topology[0].index[indices])

        #Rewrite indices and serial
        topology.reset_index(drop = True, inplace = True)
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