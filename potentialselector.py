def selectmolecule():


	selectmolecule = ''

	print('\nWelcome to Energy Minimiser v1.0 \n')

	while selectmolecule.upper() != 'A' and selectmolecule.upper() != 'B':
		selectmolecule = input('\nChoose potential (A for Morse, B for Lennard Jones): ')

	if selectmolecule.upper() == 'A':
		import morse

	if selectmolecule.upper() == 'B':
		import lennardjones



selectmolecule()
