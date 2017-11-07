"""Parses Wetter config file and finds Nmax

"""

def parse(file):
    """Parses config file
    
    """

    f = open(file, 'r')
    content = f.readlines()
    f.close()

    keyWords = ['atom', 'molecule', 'end']
    subKeyWords = ['defect', 'surface', 'HOH bond length', 'OH bond length', 'MOH bond length', 'bond angle', 'water', 'hydroxyl', 'fraction']
    atoms = []

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

        elif('end' in content[i]):
            break

        i += 1
    return atoms

def get_max_coordination(element):
    """Finds bulk coordination in MaxCoordinations.data

    """

    foundMax = False
    f = open("drizzle/MaxCoordinations.data")
    content = f.readlines()
    f.close()
    for line in content:
        if(element in line):
            colonIndex = line.index(':')
            maxCoordination = int(line[colonIndex + 1:])
            foundMax = True
            return maxCoordination
    if(not foundMax):
        raise ValueError("Could not find maximum coordination number for \"" + element + "\" in MaxCoordinations.lib")