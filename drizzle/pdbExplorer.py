"""Formatting of output pdb-file and preparation of input pdb-file.

"""

from __future__ import print_function
from radish import Topologizer
import numpy as np
import mdtraj as md


def append_atoms(file, coords=[], elements = []):
    """Append atoms to pdb file

    Parameters
    ----------
    file : string
        pdb file to which atoms will be appended
    coords : array
        1D Array of coordinates of all atoms
    elements : array
        Array containing element symbols for all atoms
    verbose : boolean
        sets verbosity
    """

    atomList = []
    fileName = file
    float_format = lambda x: "%.3f" % x
    f = open(file, "r")
    content = f.readlines()
    f.close()
    fileExt = fileName.rsplit('.', 1)[1]

    if(fileExt == "pdb"):
        # Get indices for all entries in list 'content' where the substring
        # 'ATOM' occurs
        indices = [ind for ind, line in enumerate(content) if 'ATOM' in line]

        # Get number of atoms
        nrAtoms = len(indices)

        # Variables for the .pdb format
        residueName = 'SOL'
        chainIdentifier = 'A'
        residueSequenceNumber = '2'
        occupancy = '1.00'
        temperatureFactor = '0.00'

        j = 0
        # Prepare list of atoms to be appended to pdb file
        for coord in coords:
            
            i = 4
            nrAtoms += 1
            tempString = "ATOM"

            # Format tempString to .pdb format
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

            # Append formatted tempString
            atomList.append(tempString)

        # Termination sequence
        if(content[indices[-1] + 1][:3] == 'TER'):
            pass

        # Get old content of file until last atom entry
        new_content = content[:indices[-1] + 1]
        # Append new atoms
        new_content.extend(atomList)
        # Append lines after the final atom entry
        new_content.extend(content[indices[-1] + 1:])

        # Print to file
        file = open(file, 'w')
        for line in new_content:
            file.write("%s\n" % line.rstrip())
        file.close()

        print("Added " + str(len(coords)) + " atoms to " + fileName)

    # mdtraj doesnt like xyz so no use
    elif(fileExt == "xyz"):
        i = 0
        file = open(file, 'w')
        file.write(content)
        for element in elements:
            file.write("%s\n" % element + "   " + float_format(coord[0]) +\
                       "   " + float_format(coord[1]) + "   " +\
                       float_format(coord[2]))
            i += 1
        file.close()


def remove_low_coordinated(topol, Nmax, element, verbose):
    """Removes low coordinated (<Nmax-3) atoms

    Parameters
    ----------
    topol : radish.Topologizer instance
        Topologizer instance of metal oxide crystal
    Nmax : int
        Maximum coordination for element
    element : string
        Element symbol for metal atom
    verbose : boolean
        sets verbosity

    Returns
    -------
    topol
        radish.Topologizer instance with low coordinated atoms removed
    """

    verboseprint = print if verbose else lambda *a, **k: None
    topol.topologize()
    verboseprint("\nRemoving low coordinated atoms....")

    while(1):
        indices = []

        #Create inputs for new Topologizer instance
        xyz = topol.trj.xyz[0]
        topology = topol.trj.topology.to_dataframe()

        # Get low coordinated atoms using radish
        i = 3   # Coordination less than or equal to Nmax-3
        while(i <= Nmax):
            try:
                centerIndices = topol.extract(element,\
                                              environment = {'O': Nmax - i}).\
                                                      index.get_level_values(1)
                indices.extend(centerIndices)
            # If no indices found Topologizer throws IndexError
            except IndexError:
                pass

            i += 1

        i = 0
        while i < 2:
            try:
                oxygenIndices = topol.extract('O',\
                                              environment = {element: i}).\
                                                   index.get_level_values(1)
                indices.extend(oxygenIndices)

            except IndexError:
                pass
            i += 1

        # Get uncoordinated atoms
        topologyIndices = np.array([atom.index for atom in topol.trj.topology.atoms])
        bondgraphIndices = np.array(topol.bondgraph['i'].unique())
        mask = np.ones(len(topologyIndices), np.bool)
        mask[bondgraphIndices] = 0
        indices.extend(topologyIndices[mask])

        if(len(indices) < 1):
            verboseprint("All low coordinated atoms removed\n")
            return topol

        #Remove atoms from topology
        topology = topology[0].drop(topology[0].index[indices])

        #Rewrite indices and serial
        topology.reset_index(drop = True, inplace = True)
        topology['serial'] = topology.index + 1

        #Remove atoms
        indices = sorted(indices, reverse = True)   #Sort in decreasing order
                                                    #to allow loop to find 
                                                    #correct index
        for index in indices:
            xyz = np.vstack((xyz[:index], xyz[index + 1:]))

        verboseprint("Removed " + str(len(indices)) + " low coordinated atoms")
        topology = md.Topology.from_dataframe(topology, bonds = None)
        trajectory = md.Trajectory(xyz, topology)
        topol = Topologizer.from_mdtraj(trajectory)

        topol.topologize()