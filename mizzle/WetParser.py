"""Parses the configuration file and extract the maximum coordination for the metals

"""

def parse_config(file):
    """Parses the configuration in *file*

    """

    f = open(file, "r")
    content = f.readlines()
    f.close()

    keyWords = ["atom", "molecule", "resname", "end"]
    subKeyWords = ["defect", "surface", "HOH bond length",
                   "OH bond length", "MOH bond length",
                   "bond angle", "water", "hydroxyl", "fraction"]
    atoms = []
    resname = None
    solver = []
    i = 0

    while(i < len(content)):
        if("atom" in content[i]):
            keyWordIndex = content[i].index("atom")
            colonIndex = content[i].index(":")
            element = content[i][keyWordIndex + 5:colonIndex]
            tempDict = {"element": element}
            coordination = content[i][colonIndex+2:].rstrip()
            tempDict["coordination"] = coordination

            i += 1

            while(not any(substring in content[i] for substring in keyWords)):
                for subKeyWord in subKeyWords:
                    if(subKeyWord in content[i]):
                        colonIndex = content[i].index(":")
                        tempDict[subKeyWord] = content[i][colonIndex + 2:].rstrip()
                i += 1
            atoms.append(tempDict)

            continue

        elif("resname" in content[i]):
            colonIndex = content[i].index(":")
            resname = content[i][colonIndex+1:].strip()
            if(len(resname) > 3):
                raise LookupError(("Error in config file: resname is \
                                    restricted to 3 characters"))

        elif(content[i].strip() == "end"):
            break

        i += 1
    for atom in atoms:
        assert float(atom.get('water', 0)) +\
                     float(atom.get('hydroxyl', 0)) == 1,\
                     "Water and hydroxyl fractions in config file does\
                      not sum to 1!"

    return atoms, resname

def parse_data(element):
    """Finds metal properties needed by the wetting algorithm metals.data

    """

    foundMax = False

    import os
    import mizzle.WetParser as wp

    atomList = ["Ti", "Fe", "Si"]
    metalData = {}
    maxNf = os.path.join(os.path.dirname(wp.__file__), "metals.data")
    i = 0
    foundDMOH = False
    foundDMOH2 = False
    foundMax = False

    f = open(maxNf)
    content = f.readlines()
    f.close()

    while(i < len(content)):
        if(element in content[i]):
            i += 1
            while(not any(atom in content[i] for atom in atomList)):

                if("Nmax:" in content[i]):
                    colonIndex = content[i].index(":")
                    metalData['Nmax'] = int(content[i][colonIndex + 1:])
                    foundMax = True
                
                elif("d_MOH:" in content[i]):
                    colonIndex = content[i].index(":")
                    metalData['d_MOH'] = float(content[i][colonIndex + 1:])
                    foundDMOH = True
                elif("d_MOH2:" in content[i]):
                    colonIndex = content[i].index(":")
                    metalData['d_MOH2'] = float(content[i][colonIndex + 1:])
                    foundDMOH2 = True

                if(all([foundDMOH, foundDMOH2, foundMax])):
                    return metalData
                i += 1
        i += 1

    if(not (foundMax and foundDMOH and foundDMOH2)):
        raise ValueError(("Data is missing for "
                          "element {} in {}.".format(element, maxNf)))
