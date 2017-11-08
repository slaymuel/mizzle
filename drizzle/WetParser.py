"""Parses Wetter config file and finds Nmax

"""

def parse(file):
    """Parses config file
    
    """

    f = open(file, 'r')
    content = f.readlines()
    f.close()

    keyWords = ['atom', 'molecule', 'resname', 'end']
    subKeyWords = ['defect', 'surface', 'HOH bond length', 'OH bond length', 'MOH bond length', 'bond angle', 'water', 'hydroxyl', 'fraction']
    atoms = []
    resname = None
    i = 0

    while(i < len(content)):
        if('atom' in content[i]):
            keyWordIndex = content[i].index('atom')
            colonIndex = content[i].index(':')
            element = content[i][keyWordIndex + 5:colonIndex]
            tempDict = {'element': element}
            coordination = content[i][colonIndex+2:].rstrip()
            tempDict['coordination'] = coordination

            i += 1

            while(not any(substring in content[i] for substring in keyWords)):
                for subKeyWord in subKeyWords:
                    if(subKeyWord in content[i]):
                        colonIndex = content[i].index(':')
                        tempDict[subKeyWord] = content[i][colonIndex + 2:].rstrip()
                i += 1
            atoms.append(tempDict)

            continue

        elif('resname' in content[i]):
            colonIndex = content[i].index(':')
            resname = content[i][colonIndex+1:].strip()
            if(len(resname) > 3):
                raise LookupError("Error in config file: resname can be no longer than 3 characters")

        elif(content[i].strip() == 'end'):
            break

        i += 1
    #assert()
    return atoms, resname

def get_max_coordination(element):
    """Finds bulk coordination in MaxCoordinations.data

    """

    foundMax = False

    import os
    import WetParser

    maxNf = os.path.join(os.path.dirname(WetParser.__file__), "MaxCoordinations.data")
    f = open(maxNf)
    
    content = f.readlines()
    f.close()
    for line in content:
        if(element in line):
            colonIndex = line.index(':')
            maxCoordination = int(line[colonIndex + 1:])
            foundMax = True
            return maxCoordination
    if(not foundMax):
        raise ValueError(("Could not find maximum coordination number for "
                          " {} element in {}.".format(element, maxNf)))
