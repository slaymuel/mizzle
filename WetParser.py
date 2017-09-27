def parse(file):
	f = open(file, 'r')
	content = f.readlines()
	f.close()

	keyWords = ['atom', 'molecule']
	subKeyWords = ['coordination', 'HOH bond length', 'OH bond length', 'MOH bond length', 'bond angle', 'water', 'hydroxyl']
	atoms = []
	molecules = []

	i = 0
	while(i < len(content)):
		if('atom' in content[i]):
			keyWordIndex = content[i].index('atom')
			colonIndex = content[i].index(':')
			element = content[i][keyWordIndex + 5:colonIndex]
			print("Center: ")
			print(element)
			tempDict = {'element': element}

			i += 1

			while(not any(substring in content[i] for substring in keyWords)):
				for subKeyWord in subKeyWords:
					if(subKeyWord in content[i]):
						colonIndex = content[i].index(':')
						tempDict[subKeyWord] = content[i][colonIndex + 2:].rstrip()
				i += 1
			atoms.append(tempDict)
			print tempDict


		elif('molecule' in content[i]):
			keyWordIndex = content[i].index('molecule')
			colonIndex = content[i].index(':')
			molecule = content[i][keyWordIndex + 9:colonIndex]
			print("molecule: ")
			print(molecule)

			i += 1

			while(not any(substring in content[i] for substring in keyWords)):
				for subKeyWord in subKeyWords:
					if(subKeyWord in content[i]):
						colonIndex = content[i].index(':')
						tempDict[subKeyWord] = content[i][colonIndex + 2:].rstrip()
				i += 1
			molecules.append(tempDict)
			print tempDict
		tempDict = {}
		i += 1
		
