from numpy import *
import math
from os import remove, environ, getcwd
from ase import Atoms
from ase.calculators.gulp import GULP, Conditions
from ase import io
from pathlib import Path


def GULPRelaxation(inputStr):
	##for gulp to work. gulp folder must be in the same folder as the script
	environ["GULP_LIB"] = "\"" + getcwd() + "/Gulp_40/Libraries/\""
	environ["ASE_GULP_COMMAND"] = "\"" + getcwd() + "/Gulp_40/Exe/gulp.exe\" < PREFIX.gin > PREFIX.got"
	
	##takes string inputStr and turns it into stepNb, shouldRestartFromScratch listOfAtoms
	if len(inputStr)%4 != 2 :
		print("Args nb " + str(len(inputStr)))
		print(inputStr)
		print("Not enough args. Args order : stepNb, shouldRestartFromScratch (0 or 1), \n atom list formatted as Name Xpos Ypos Zpos for each atom (example CO2 would be C 0 0 0 O 0 1 1 O 0 1 -1)") 
		return "", False
	
	atomNames = []
	atomXYZ = []  
	
	stepNb = int(inputStr[0])
	shouldRestartFromScratch = int(inputStr[1])
	
	#if restart from scratch, destroy previously known trajectory
	if (shouldRestartFromScratch == 1) :
		my_file = Path('tmp.pckl')
		if my_file.is_file() :
			remove('tmp.pckl')		
	
	#create Atoms object to optimize
	for i in range(0, math.floor(len(inputStr) / 4)) :
		atomNames.append(inputStr[4 * i + 2])
		atomXYZ.append((float(inputStr[4 * i + 3]), float(inputStr[4 * i + 4]), float(inputStr[4 * i + 5])))
		
	atoms = Atoms(atomNames, positions = atomXYZ)
	c = Conditions(atoms)
	
	
	#attach calculator
	calc = GULP(keywords='conp',  library="./Gulp_40/Libraries/reaxff.lib")
	atoms.calc = calc
	calc.set(keywords='conp opti')
	
	#optimize
	opt = calc.get_optimizer(atoms)
	opt.run(fmax=0.05)
		
	#io.write('../MoleculeLibrary/tmp.xyz', atoms, format='xyz')
	
	pos = atoms.get_positions()
	outStr = ""
	
	for i in range(0, len(atomNames)):
		outStr += atomNames[i] + " "
		outStr += str(pos[i][0]) + " "
		outStr += str(pos[i][1]) + " "
		outStr += str(pos[i][2]) + " "
	
	
	return outStr, True