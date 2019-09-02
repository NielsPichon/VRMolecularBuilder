import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.offline import plot
import time
from operator import attrgetter
#import GPUtil
from copy import deepcopy
#from scipy.sparse import coo_matrix, bmat
GPU = False
try:
	import cupy as cp
	from cupy import cuda
	cp.get_default_memory_pool().free_all_blocks()
	cp.get_default_pinned_memory_pool().free_all_blocks()
	GPU = True
except:
	GPU = False
#GPU = False
vpp = 2.7
extra = 1.2
cap = False

#from collections import deque
values = [0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67, 1.55, 1.39, 1.26, 1.16, 1.11, 1.03, 0.99, 0.96, 1.96, 1.71, 1.48, 1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.1, 1.12, 1.18, 1.24, 1.21, 1.21, 1.16, 1.14, 1.17, 2.1, 1.85, 1.63, 1.54, 1.47, 1.38, 1.28, 1.25, 1.25, 1.2, 1.28, 1.36, 1.42, 1.4, 1.4, 1.36, 1.33, 1.31, 2.32, 1.96, 1.8, 1.63, 1.76, 1.74, 1.73, 1.72, 1.68, 1.69, 1.68, 1.67, 1.66, 1.65, 1.64, 1.7, 1.62, 1.52, 1.46, 1.37, 1.31, 1.29, 1.22, 1.23, 1.24, 1.33, 1.44, 1.44, 1.51, 1.45, 1.47, 1.42, 2.23, 2.01, 1.86, 1.75, 1.69, 1.7, 1.71, 1.72, 1.66, 1.66, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57]
keys = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
CovRadDict = dict(zip(keys, values))

#orbitalArrays = [coo_matrix([[2.7]]) for n in range(len(keys))]
#OrbitalDict = dict(zip(keys, orbitalArrays))
''' The data-structure is mainly atom objects that contain all the important information.
    The molecule then contain lists of these atom objects, the full molecule and the leads seperately.
	In some functions i use get___List to get a list of a certain property of all the atoms and use that
	because im lazy. Then i update the atom objects in the end, like in "sortMolecule".'''

class Atom:
	def __init__(self, element, position, neighbours, lead, leadDist = 0, neighbourAtoms = [], discovered = False, index = int):
		self.element = element #string, example "C" for carbon
		self.covalentRadius = CovRadDict[self.element]
		self.position = position #3 floats, example [1.0, 0.0, 0.0] in angstrom
		self.neighbourIndices = neighbours #indices of neighbours, varying nr. of int
		self.neighbourAtoms = neighbourAtoms #This is filled for a given molecule
		self.index = index
		self.lead = lead
		self.leadDist = leadDist
		self.discovered = discovered
		#self.orbitals = OrbitalDict[self.element] # Orbitals intended for Slater-Koster implementation
			
	def copy(self):
		# All but neighbourAtoms is deepcopied. DO NOT deepcopy neighbourAtoms, the structure is so deep it will reach recursion depth error.
		return Atom(deepcopy(self.element), deepcopy(self.position), deepcopy(self.neighbourIndices), deepcopy(self.lead), deepcopy(self.leadDist), self.neighbourAtoms, deepcopy(self.discovered), deepcopy(self.index))
	
	
class Molecule:
	def __init__(self):
		self.molecule = []
		self.nrOfAtoms = 0
		self.leftLead = []
		self.rightLead = []
		
	def addAtom(self, atom):
		""" Adds atom to molecule, the indices of atoms """
		self.molecule.append(atom)
		if atom.lead == -1:
			self.leftLead.append(atom)
		if atom.lead == 1:
			self.rightLead.append(atom)
		atom.index = self.nrOfAtoms
		self.nrOfAtoms += 1
	
	def getLeadList(self):
		return [atom.lead for atom in self.molecule]
	
	def getLeadDistList(self, default = None):
		if default == None:
			default = self.molecule
		return [atom.leadDist for atom in default]
		
	def getPositionList(self, default = None):
		if default == None:
			default = self.molecule
		return [atom.position for atom in default]

	def updateNeighbourAtoms(self):
		""" Converts the neighbour indices to the neighbour's Atom-objects.
			This is not a general neighbour updater routine, but a one time function to launch
			after all atoms are imported."""
		for atom in self.molecule:
			atom.neighbourAtoms = [self.molecule[i] for i in atom.neighbourIndices]
	
	def makeHamiltonian(self, V, cap):
		""" Makes hamiltonian by looping through Atom-objects """
		hamiltonian = np.zeros([self.nrOfAtoms, self.nrOfAtoms], dtype=np.complex128)
		
		for atom in self.molecule:			
			if not (atom.lead == 0) and cap:
				hamiltonian[atom.index, atom.index] = complex(0,self.CAP(atom))
			
			for neigh in atom.neighbourAtoms:
				if neigh not in self.molecule:
					continue
				
				#r = distanceToObject(atom,neigh)
				
				if neigh != atom:
					
					# Create list of second nearest neighbour elements
					secondNeighElems = [secNN.element for secNN in neigh.neighbourAtoms]
					
					# "fluor model" is changing hopping to carbon to 0 if second nearest neighbour is not carbon
					if neigh.element != "C" or atom.element != "C":
						hamiltonian[atom.index, neigh.index] = 0
					elif any(secNN != "C" for secNN in secondNeighElems):
						hamiltonian[atom.index, neigh.index] = 0
						hamiltonian[atom.index, atom.index] = -1
					else:
						hamiltonian[atom.index, neigh.index] = V#*r

	
		return hamiltonian
	
	def sortMolecule(self, leftCopyAxis, rightCopyAxis):
		#First sort for the leads
		self.molecule = sorted(self.molecule, key=attrgetter('lead'))
		#Then sort the leads after the closest atom to an arbitrary atom in the other lead
		
		# Here the leads are sorted after their projected lead Distance
		self.leftLead = sorted(self.leftLead, key=attrgetter('leadDist'), reverse = True)
		self.rightLead = sorted(self.rightLead, key=attrgetter('leadDist'))
		
		#Replace the respective leads in the molecule.
		self.molecule[0:len(self.leftLead)] = self.leftLead #leftLead
		self.molecule[-len(self.rightLead):] = self.rightLead #rightLead


		#Set each atom's index as the index in the list
		for idx, atom in enumerate(self.molecule):
			atom.index = idx
			
	
	def periodicityDirection(self, lead, plane):
		# This function calculated the direction in which leads will be replicated
		if self.leftLead == [] and self.rightLead == []:
			print('No leads were found')
			return 
		if lead == 1:
			lead = self.rightLead
		elif lead == -1:
			lead = self.leftLead
		else:
			print('Chosen lead is not an option')
			return
		
		# LeadEdges are calculated from the perspective of all atoms of given lead type
		leadEdgeAtoms, deviceEdgeAtoms = calcLeadEdge(lead)
		
		#If only one atom in the lead is connected to the device, the direction line
		#will be fitted using positions of all atoms in lead + the edge atom in device
		if len(leadEdgeAtoms) == 1:
			# Uses average of device edge atoms and whole lead to fit line
			coords = np.array(self.getPositionList(deviceEdgeAtoms))
			coords = np.append(coords, self.getPositionList(lead),axis = 0)
			meanCoord = coords.mean(axis=0)
			lineVector = np.linalg.svd(coords-meanCoord)[-1][0]
			
			
			# Fix the direction of the normal vector, using vector from random lead edge to random device edge
			deviceEdge = deviceEdgeAtoms[0].position
			if np.dot(lineVector, coords[1]-deviceEdge) < 0:
				lineVector = -lineVector
			
			lineVector = np.append(lineVector, np.dot(lineVector,coords[1]))
			return lineVector, leadEdgeAtoms
		else:
			coords = np.array(self.getPositionList(leadEdgeAtoms))
			coords = projToPlane(plane, coords)
			
			# Here we take take two random coordinates from the lead edge
			v1 = coords[-1:]-coords[0]
			v2 = plane[:3]
			planeNormVector = np.cross(v1, v2)[0]
			planeNormVector = planeNormVector/np.linalg.norm(planeNormVector)
			
			# Fix the direction of the normal vector, using vector from random lead edge to random device edge
			deviceEdge = deviceEdgeAtoms[0].position
			if np.dot(planeNormVector, coords[0]-deviceEdge) < 0:
				planeNormVector = -planeNormVector
				
			planeNormVector = np.append(planeNormVector, np.dot(planeNormVector,coords[0]))
			return planeNormVector, leadEdgeAtoms
		
	
	def leadCopy(self, lead, copyAxis, princPlane ,repeats = 1):
		if lead == 1:
			lead = self.rightLead
			print('Copying right lead...')
		elif lead == -1:
			lead = self.leftLead
			print('Copying left lead...')
		else:
			print('Chosen lead is not an option')
			return
		
		# Calculate the projected distance for all atoms in the lead
		print('Calculating projected distance')
		for atom in lead:
			inPPlane = projToPlane(princPlane, [atom.position])[0]
			atom.leadDist = distanceToObject(inPPlane, copyAxis)
		
		print('Finding uniques')
		# Append the unique projected positions to projDistList in order
		projDistList = []
		atomDistList = [atom.leadDist for atom in lead]
		threshold = 0.01
		for item in sorted(atomDistList):
			if len(projDistList) == 0 or item > projDistList[-1] + threshold:
				projDistList.append(item)
		
		print('Finding average segment length')
		# Finding length of each consecutive segment and modified average
		if len(projDistList) > 1:
			projSegList = [projDistList[i + 1]-projDistList[i] for i in range(len(projDistList) - 1)]
			averSegLen = (sum(projSegList)+2*lead[0].covalentRadius)/(len(projSegList)+1)
			
		else:
			###### Warning if leads of multiple material is implemented #######
			#-0.001 will move the point to the very edge of the valid range
			averSegLen = 2*lead[0].covalentRadius-0.001
			
		# The projected tile distance should always be smaller than the bond length
		tileDist = max(projDistList) - min(projDistList) + averSegLen
		if tileDist > 2*lead[0].covalentRadius*extra:
			tileDist = 2*lead[0].covalentRadius*extra
	
		# Lists to save atom and neightbour indices, inside the unit cell
		insideAtoms = np.array([], dtype = int)
		insideNeigh = np.array([], dtype = int)
		
		print('Creating lead repeats')
		# Searches the lead for internal neighbours
		for idxA, atom in enumerate(lead):
			for nn in atom.neighbourAtoms:
				for idxN, insideAtom in enumerate(lead):
					if nn == insideAtom:
						insideAtoms= np.append(insideAtoms, idxA)
						insideNeigh = np.append(insideNeigh, idxN)
		
		# Manual copy function is used to avoid the recursive property of deepcopy. 
		# because the neighbourAtoms indirectly reference each other, which make a very deep structure
		newLeadMaster = lead
		
		#Create manual copy and apply neighbour indices to append new neighbours
		for n in range(repeats-1):
			newLeadSlave = [atom.copy() for atom in newLeadMaster]
			if n == 0:
				# newtonRaphson will search for neighbour unit cells and optimise tiling distance
				insideAtoms_, outsideNeigh, newTileDist, newLeadMaster = self.newtonRaphson(newLeadSlave, lead, tileDist, copyAxis, insideAtoms, insideNeigh)
			else:
				for idx, atom in enumerate(newLeadSlave):
					atom.position += newTileDist*copyAxis[:3]
					atom.leadDist += newTileDist
					atom.neighbourAtoms = []
					
					#Appends neighbours inside the newLeadSlave
					for newNeigh in insideNeigh[insideAtoms == idx]:
						atom.neighbourAtoms.append(newLeadSlave[newNeigh])
						
					# First iteration calculates the neighbours outside the unitcell
					# Subsequent iterations apply the same outside neighbours
					for newNeigh in outsideNeigh[insideAtoms_ == idx]:
						atom.neighbourAtoms.append(newLeadMaster[newNeigh])
						newLeadMaster[newNeigh].neighbourAtoms.append(atom)
			
			# Appends the atom to the lead afterwards so they aren't included in the 'outside' search
			for atom in newLeadSlave:
				self.addAtom(atom)
			
			newLeadMaster[:] = newLeadSlave
		
	def CAP(self, atom):
		hbar = 1 
		m = 0.511e6 # eV
		c = 2.62
		if atom.lead == 1:
			projList = self.getLeadDistList(self.rightLead)
		elif atom.lead == -1:
			projList = self.getLeadDistList(self.leftLead)
		z1 = min(projList)
		z2 = max(projList)+0.1
		dz_CAP = z2 - z1
		z = abs(atom.leadDist)
		fz = (4/(c**2)) * ((dz_CAP/(z2-2*z1+z))**2 + (dz_CAP/(z2-z))**2 - 2 )
		Wz = ((hbar**2)/(2*m)) * (2*np.pi/(dz_CAP/2000))**2 * fz
		return Wz
	
	
	def newtonRaphson(self, newLeadSlave, lead, tileDist, copyAxis, insideAtoms, insideNeigh):
		# DYI Newton's method
		def calcNeighbours(addDist, newLeadSlave):
			nonlocal neighbourDists, theoryDists, insideAtoms_, outsideNeigh, lead, tileDist, copyAxis, insideAtoms, insideNeigh
			for idx, atom in enumerate(newLeadSlave):
				atom.position += addDist*copyAxis[:3]
				atom.leadDist += addDist
				atom.neighbourAtoms = []
				
				#Appends neighbours inside the newLeadSlave
				for newNeigh in insideNeigh[insideAtoms == idx]:
					atom.neighbourAtoms.append(newLeadSlave[newNeigh])
					
				# First iteration calculates the neighbours outside the unitcell
				for idxN, outsideAtom in enumerate(lead):
					#neighbourCheck also appends the atoms to each others neighbour list
					check, r = nearestNeighbourCheck(atom, outsideAtom)
					if check:
						neighbourDists = np.append(neighbourDists, r)
						theoryDists = np.append(theoryDists,atom.covalentRadius+outsideAtom.covalentRadius)
						insideAtoms_ = np.append(insideAtoms_, idx)
						outsideNeigh = np.append(outsideNeigh, idxN)
			return newLeadSlave
		
		# Prepare empty matrices
		neighbourDists = np.array([], dtype = np.float64)
		theoryDists = np.array([], dtype = np.float64)
		insideAtoms_ = np.array([], dtype = int)
		outsideNeigh = np.array([], dtype = int)
		addDist = tileDist*extra 
		
		# Calculate unit cell neighbours once
		newLeadSlave = calcNeighbours(addDist, newLeadSlave)
		totalDist = np.array([addDist], dtype = np.float64)
		
		# Conduct Newton raphson 
		if neighbourDists.size > 0:
			func = 1/len(neighbourDists)*sum(neighbourDists-theoryDists)
			while abs(func) > 0.1:
				neighbourDists = np.array([], dtype = np.float64)
				theoryDists = np.array([], dtype = np.float64)
				insideAtoms_ = np.array([], dtype = int)
				outsideNeigh = np.array([], dtype = int)
				
				newLeadSlave = calcNeighbours(-func/2, newLeadSlave)
				totalDist = np.append(totalDist, -func/2)
		
				func = 1/len(neighbourDists)*sum(neighbourDists-theoryDists)
		else:
			print('Something went wrong; no unit cell neighbours found in newton Rapson')
		
		return insideAtoms_, outsideNeigh, sum(totalDist), newLeadSlave
				
def calcLeadEdge(lead):
	leadEdgeAtoms, deviceEdgeAtoms = [], []
	for atom in lead:
		for neighbour in atom.neighbourAtoms:
			if neighbour not in lead:
				if neighbour not in deviceEdgeAtoms:
					deviceEdgeAtoms.append(neighbour)
				if atom not in leadEdgeAtoms:
					leadEdgeAtoms.append(atom)

	return leadEdgeAtoms, deviceEdgeAtoms        
	
def distanceToObject(atom, obj):
	# Calculate distance between an atoms object and either another atom, line or plane
	if type(obj)==Atom:
		# This is atom-case
		try:
			return np.linalg.norm(np.array(atom.position)-np.array(obj.position))
		except:
			print("Couldn't calculate distance between two atoms")
	elif type(obj)==np.ndarray and len(obj) == 3:
		# This is line-case
		if type(atom) == Atom:
			atom = atom.position
		return abs(np.dot(atom, obj))
	elif type(obj)==np.ndarray and len(obj) == 4:
		# This is plane-case
		if type(atom) == Atom:
			atom = atom.position
		return abs((np.dot(atom, obj[:3])-obj[-1]))
	
def nearestNeighbourCheck(i, k):
	result = False
	dist = distanceToObject(i, k)
	maxDist = ((i.covalentRadius + k.covalentRadius)*extra)
	# Calculate whether two atoms objects are neighest neighbours depending on their covalent radius
	if dist < maxDist and not (i == k) and i.lead == k.lead:
		if k not in i.neighbourAtoms:
			i.neighbourAtoms.append(k)
			result = True
		if i not in k.neighbourAtoms:
			k.neighbourAtoms.append(i)
			result = True
		return result, dist

	else:
		return False, dist
	
def projToPlane(plane, oldCoords):
	# Takes a list of coordinates even if only one coordinate (eg. [[x,y,z]])
	planeNorm = plane[:3]
	newCoords = np.array([])
	if len(oldCoords) == 1:
		newCoords = np.hstack((newCoords,oldCoords[0]-(np.dot(planeNorm,oldCoords[0])-plane[3])*planeNorm))
	else:
		for i in oldCoords:
			newCoords = np.hstack((newCoords,i-(np.dot(planeNorm,i)-plane[3])*planeNorm))
	
	newCoords = newCoords.reshape((int(len(newCoords)/3),3))
	newCoords = np.round(newCoords,4)
	return newCoords
	
'''TODO: This function is not needed for implementation, only for sanity checks.'''
def plotMeanPlane(coords, planeModel):
	fig = go.Figure()
	x,y,z = np.transpose(coords)
	
	#plotly library plotting
	fig.add_trace(go.Scatter3d(x=x,y=y,z=z,
	mode='markers',
	marker=dict(
		size=7,
		color=z,                # set color to an array/list of desired values
		colorscale='Viridis',   # choose a colorscale
		opacity=0.8
		)
	))

	if len(planeModel) == 3:
		a,b,c = planeModel
		d = 0
	else:
		a, b, c, d = planeModel
	
	# calculate plane equation inside specified area
	xx, yy = np.mgrid[min(x):max(x), min(y):max(y)]
	z1 = (d - a * xx - b * yy) / c
	
	fig.add_trace(go.Surface(x=xx, y=yy, z=z1, showscale=False , opacity=0.4, hoverinfo='skip'))
	fig.update_layout(margin=dict(l=0, r=0, b=0, t=0), scene_aspectmode='data')
	
	plot(fig)

def PlanefitPCA(molecule):
	# PCA as found in scikit-learn
	xyz = molecule.getPositionList()
	
	# move coordinates around the the mean value
	xyz -= np.mean(xyz, axis = 0)
	
	# compute Singular Value Decomposition of the coordinates, V will contain
	# the three principle vectors, ordered so the third is the plane normal vector
	U, S, V = np.linalg.svd(xyz, full_matrices=False)
	U, V = svd_flip(U, V)
	
	# calculate distance to the plane for the plane equation
	d=np.dot(V[2], molecule.getPositionList()[0])
	return np.append(V[2],d)

def svd_flip(u, v, u_based_decision=True):
	"""Sign correction to ensure deterministic output from SVD.

	Adjusts the columns of u and the rows of v such that the loadings in the
	columns in u that are largest in absolute value are always positive.
	
	As found in Scikit-learn
	"""
	if u_based_decision:
		# columns of u, rows of v
		max_abs_cols = np.argmax(np.abs(u), axis=0)
		signs = np.sign(u[max_abs_cols, range(u.shape[1])])
		u *= signs
		v *= signs[:, np.newaxis]
	else:
		# rows of v, columns of u
		max_abs_rows = np.argmax(np.abs(v), axis=1)
		signs = np.sign(v[range(v.shape[0]), max_abs_rows])
		u *= signs
		v *= signs[:, np.newaxis]
	return u, v
    
def inputDecoder(inputStr):
    #Takes the input string and construct output lists used in the transmission calculation
	eMin, eMax, eSteps = inputStr[0], inputStr[1], inputStr[2]
	inputStr = inputStr[3:]
	elementList = []
	positionList = []
	leadList = []
	neighbourList = []
	for curIndex, text in enumerate(inputStr):
		# The search is seperated by 'alphabet string', after which subsequent numbers is identified
		# It is safe to add more numbers before the first letter.
		if text.isalpha():
			elementList.append(text)
			# fixed length inputs
			try:
				positionList.append([float(inputStr[curIndex+1]), float(inputStr[curIndex+2]),float(inputStr[curIndex+3])])
				leadList.append(int(inputStr[curIndex+4]))
			except Exception as e:
				return None, None, None, None, None, None, None, e
			
			# Neighbours is variable length input, reads until an alphabet string is reached
			i = 5
			neighbourRow = []
			try:
				inputStr[curIndex+i]
			except IndexError:
				neighbourList.append(neighbourRow)
				break
            
			while (not inputStr[curIndex+i].isalpha()):
				neighbourRow.append(int(inputStr[curIndex+i]))
				if curIndex+i == len(inputStr)-1:
					break
				i += 1
			neighbourList.append(neighbourRow)
	return eMin, eMax, eSteps, elementList, positionList, leadList, neighbourList, False

def clusterSearch(startMolecule):
	# Make list of all atoms in the startMolecule that are not leads
	G = [atom for atom in startMolecule.molecule if atom.lead == 0]
	clusterList = []
	
	# Search all undiscovered atoms, .discovered value is handled in function 'buildCluster'
	for atom in G:
		if (not atom.discovered):
			newCluster = buildCluster(atom)
			if newCluster != None:
				clusterList.append(newCluster)
			
	return clusterList

def buildCluster(startAtom):
	# Modified depth-first search algorithm, where input is a non-lead atom
	containLlead = False
	containRlead = False
	discoveredLeads = []
	
	cluster = Molecule()
	cluster.addAtom(startAtom)
	startAtom.discovered = True
	Q = [startAtom]
	while len(Q) > 0:
		atom = Q.pop()
		for nAtom in atom.neighbourAtoms:
		
			# This ensures that the search go into the leads, but cannot leave back to non-leads
			if not atom.lead == 0 and not nAtom.lead == atom.lead:
				continue
			
			# This is the depth-first search part
			if nAtom.discovered == False:
				nAtom.discovered = True
				cluster.addAtom(nAtom)
				Q.append(nAtom)
				
			# Here we test for validity of the given cluster 
			if nAtom.lead == -1:
				containLlead = True
				discoveredLeads.append(nAtom)
			elif nAtom.lead == 1:
				containRlead = True
				discoveredLeads.append(nAtom)
	
	# Reset the discovered status of leads so the next cluster will also search them
	for leadAtom in discoveredLeads:
		leadAtom.discovered = False
	
	# If the cluster contains both leads it is valid and is returned.
	if (containLlead and containRlead):
		return cluster
	else:
		return None

def clusterMerge(validClusterList):
	# Merges all the valid clusters assuming they share the same lead
	try:
		combinedMolecule = validClusterList[0]
	except:
		print('No valid clusters were present in cluster merge')
		return []
		
	# Get all atoms from the leads of the original molecule
	if len(validClusterList) > 1:
		for cluster in validClusterList[1:]:
			for atom in cluster.molecule:
				if atom not in combinedMolecule.molecule:
					combinedMolecule.addAtom(atom)
					
		# All neighbours will be updated according to the cluster
		for atom in combinedMolecule.molecule:
			atom.neighbourAtoms[:] = [neigh for neigh in atom.neighbourAtoms if neigh in combinedMolecule.molecule]
	else:
		for atom in combinedMolecule.molecule:
			atom.neighbourAtoms[:] = [neigh for neigh in atom.neighbourAtoms if neigh in combinedMolecule.molecule]

	return combinedMolecule

def greensFunctionInvert(M):
	# Old implementation to test numba @jit
	result = np.empty(np.shape(M))
	for idx, matrix in enumerate(M[:,:]):
		result[idx,:,:] = np.linalg.inv(matrix)
	return result

###############################################################################
# The following functions are used to compute the inverse of a batch of       #
# function on the GPU using CuPy.                                             #
# The code is modified from https://github.com/cupy/cupy/issues/1647          #

def _as_batch_mat(x):
	return x.reshape(len(x), x.shape[1], -1)


def _mat_ptrs(a):
	if len(a) == 1:
		return cp.full((1,), a.data.ptr, dtype=np.uintp)
	else:
		stride = a.strides[0]
		ptr = a.data.ptr
		out = cp.arange(ptr, ptr + stride * len(a), stride, dtype=np.uintp)
		return out


def _get_ld(a):
	strides = a.strides[-2:]
	trans = np.argmin(strides)
	return trans, int(max(a.shape[trans - 2], max(strides) // a.itemsize))


def inv_gpu(b):
	# We do a batched LU decomposition on the GPU to compute the inverse
	# Change the shape of the array to be size=1 minibatch if necessary
	# Also copy the matrix as the elments will be modified in-place
	a = _as_batch_mat(b).copy()
	
	n = a.shape[1]
	n_matrices = len(a)
	# Pivot array
	p = cp.empty((n, n_matrices), dtype=np.int32)
	# Output array
	c = cp.empty_like(a)
	# These arrays hold information on the execution success
	# or if the matrix was singular
	info = cp.empty(n_matrices, dtype=np.int32)
	ap = _mat_ptrs(a)
	cpv = _mat_ptrs(c)
	_, lda = _get_ld(a)
	_, ldc = _get_ld(c)
	handle = cuda.Device().cublas_handle
	if b.dtype == np.float32:
		cuda.cublas.sgetrfBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, info.data.ptr, n_matrices)
		cuda.cublas.sgetriBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, cpv.data.ptr, ldc,
			info.data.ptr, n_matrices)
	elif b.dtype == np.float64:
		cuda.cublas.dgetrfBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, info.data.ptr, n_matrices)
		cuda.cublas.dgetriBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, cpv.data.ptr, ldc,
			info.data.ptr, n_matrices)
	elif b.dtype == np.complex64:
		cuda.cublas.cgetrfBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, info.data.ptr, n_matrices)
		cuda.cublas.cgetriBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, cpv.data.ptr, ldc,
			info.data.ptr, n_matrices)
	elif b.dtype == np.complex128:
		cuda.cublas.zgetrfBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, info.data.ptr, n_matrices)
		cuda.cublas.zgetriBatched(
			handle, n, ap.data.ptr, lda, p.data.ptr, cpv.data.ptr, ldc,
			info.data.ptr, n_matrices)
	else:
		print('Matrix had invalid dtype. Valid dtypes: float32, float64, complex64, complex128')
	return c, info

###############################################################################
# Put in TCPphysics to precompile c++ functions and reserve GPU memory
# H = np.eye(x, dtype='complex64')
# Where x is the expected limit for GPU accelerated
# TransportCalculation(H,H,H,[0])

def transportCalculation_GPU(hamiltonian, rightGamma, leftGamma, energy, seR, seL, eta = 0):
	# Large _gpu objects are stored in gpu memory and is not accessible(not even in variable explorer)
	# until they are transfered back to the cpu via .asnumpy()
	
	M = np.stack([(e+eta*1j)*np.eye(len(hamiltonian))-hamiltonian for e in energy])
	
	# Change arrays into cupy arrays which put them in GPU memory.
	M_gpu = cp.array(M).astype(cp.complex64)
	rG_gpu = cp.array(rightGamma).astype(cp.complex64)
	lG_gpu = cp.array(leftGamma).astype(cp.complex64)
	
	if len(rightGamma.shape) == 2:
		G_gpu = inv_gpu(M_gpu)[0]
		result = cp.einsum('ij,ojk,kl,oil->o', rG_gpu, G_gpu, lG_gpu, cp.conj(G_gpu))
	elif len(rightGamma.shape) == 3:
		seR_gpu = cp.array(seR).astype(cp.complex64)
		seL_gpu = cp.array(seL).astype(cp.complex64)
		G_gpu = inv_gpu(M_gpu-seR_gpu-seL_gpu)[0]
		result = cp.einsum('oij,ojk,okl,oil->o', rG_gpu, G_gpu, lG_gpu, cp.conj(G_gpu))
	# Example of alternative to einsum
	# result = cp.trace(cp.matmul(cp.matmul(cp.matmul(rG_gpu, G_gpu),lG_gpu),cp.conj(G_gpu)),axis1 = 1, axis2 = 2)
	return cp.asnumpy(result)

def transportCalculation(hamiltonian, rightGamma, leftGamma, energy, seR, seL, eta=0.001):
	# Change to single precision for improved speed
	hamiltonian = hamiltonian.astype('complex64')
	rightGamma = rightGamma.astype('complex64')
	leftGamma = leftGamma.astype('complex64')
	
	# stack matrices into a 3D ndarray
	M = np.stack([(e+eta*1j)*np.eye(len(hamiltonian))-hamiltonian for e in energy])
	
	
	if len(rightGamma.shape) == 2:
		# linalg.inv will parallelize the inversion of all matrices in the M stack.
		G = np.linalg.inv(M)
		
		# A guide to convoluted einsum notation http://ajcr.net/Basic-guide-to-einsum/
		# Notice, dot products are done like: 'ij,jk,kl,li', but li is transposed to il.
		# The start and ending index are both 'i' which causes the trace to be calculated.
		# The index 'o' is the 'energy' index, and is the index we preseve by ->o 
		# 'optimal' algorithm: small overhead for small matrices, essential for speed for medium and large matrices
		result = np.einsum('ij,ojk,kl,oil->o', rightGamma, G, leftGamma, G.conj(), optimize = 'optimal')
	elif len(rightGamma.shape) == 3:
		G = np.linalg.inv(M-seR-seL)
		result = np.einsum('oij,ojk,okl,oil->o', rightGamma, G, leftGamma, G.conj(), optimize = 'optimal')
	return result

def capHamiltonian(molecule, en):
	# Create the full hamiltonian
	print("Creating the CAP Hamiltonian")
	hamiltonian = molecule.makeHamiltonian(vpp, True)
	
	# Slice the left and right self-energy from the full hamiltonian
	print("Slicing matrices")
	leftSE = np.zeros_like(hamiltonian)
	rightSE = np.zeros_like(hamiltonian)
	leftIdcs = slice(0,len(molecule.leftLead),1)
	rightIdcs = slice(-len(molecule.rightLead),None,1)
	leftSE[leftIdcs, leftIdcs] = hamiltonian[leftIdcs, leftIdcs]
	rightSE[rightIdcs, rightIdcs] = hamiltonian[rightIdcs, rightIdcs]
	
    # Calculate Gamma matrices
	leftGamma = -2*leftSE.imag
	rightGamma = -2*rightSE.imag
	return hamiltonian, rightGamma, leftGamma, None, None

def cap2Hamiltonian(molecule, en):
	# Attempted implementation of https://aip.scitation.org/doi/pdf/10.1063/1.4898729
	# Notice; tiles should be much larger than normal.
	
	# Create the full hamiltonian
	print("Creating the CAP Hamiltonian")
	ham = molecule.makeHamiltonian(vpp, True)
	global hamiltonian
	hamiltonian = ham
	# Creating slices
	leftFullLen = len(molecule.leftLead)
	rigthFullLen = len(molecule.rightLead)
	n = 15
	Nl = n
# =============================================================================
# 	Nr = n
# =============================================================================
	
	leftSlice  = slice(leftFullLen-Nl,leftFullLen,1)
	deviceSlice = slice(leftFullLen,-rigthFullLen,1)
# =============================================================================
# 	rightSlice = slice(-rigthFullLen,rigthFullLen+Nr,1)
# =============================================================================
	
	hL = ham[leftSlice, leftSlice]
	vDL= ham[deviceSlice, leftSlice]
	hD = ham[deviceSlice, deviceSlice]
# =============================================================================
# 	vDR= ham[rightSlice, deviceSlice]
# 	hR = ham[rightSlice, rightSlice]
# =============================================================================
	eg = time.time()
	ekL,eigenVectorL = np.linalg.eig(hL)
	#ekStarL,eigenWectorL = np.linalg.eig(np.conj(hL))
	print('eig took {}'.format(time.time()-eg))
# =============================================================================
# 	ekR,eigenVectorR = np.linalg.eig(hR)
# 	ekStarR,eigenWectorR = np.linalg.eig(np.conj(hR))
# =============================================================================
	
	gL = np.empty((len(en), len(eigenVectorL), len(eigenVectorL)), dtype=np.complex128)
	for k, ek in enumerate(ekL):
		gL = np.add(gL,np.stack([np.outer(eigenVectorL[:,k], np.conj(eigenVectorL[:,k]))/(E-ek) for E in en]))
	
# =============================================================================
# 	gR = np.empty((len(en), len(eigenVectorR), len(eigenWectorR)), dtype=np.complex128)
# 	for k, ek in enumerate(ekR):
# 		gR = np.add(gR, np.stack([np.outer(eigenVectorR[:,k], np.conj(eigenWectorR[:,k]))/(E-ek) for E in en]))
# =============================================================================
	
	seL = (vDL) @ gL @ dagger(vDL)
# =============================================================================
# 	seR = dagger(vDR) @ gR @ (vDR)
# =============================================================================
	
	# Calculate Gamma functions
	leftGamma = -2*seL.imag
# =============================================================================
# 	rightGamma = -2*seR.imag
# =============================================================================
	return hD, leftGamma, leftGamma, seL, seL
	

def recurveHamiltonian(molecule, tiles, en):
	print("Creating the recurve Hamiltonian")
	# Creating hamiltonian with no CAP
	ham = molecule.makeHamiltonian(vpp, False)
	
	# Creating slices
	leftFullLen = len(molecule.leftLead)
	rigthFullLen = len(molecule.rightLead)
	leftLen = int(leftFullLen/tiles)
	rightLen = int(rigthFullLen/tiles)

	leftSlice1  = slice(leftFullLen-2*leftLen,leftFullLen-leftLen,1)
	leftSlice2 = slice(leftFullLen-leftLen,leftFullLen,1)
	deviceSlice = slice(leftFullLen,-rigthFullLen,1)
	rightSlice2 = slice(-rigthFullLen, -rigthFullLen+rightLen,1)
	rightSlice1 = slice(-rigthFullLen+rightLen,-rigthFullLen+2*rightLen,1)
	
	#Picking out the device and lead arrays closest to the device region
	hL = ham[leftSlice1, leftSlice1]
	vL = ham[leftSlice1, leftSlice2]
	vDL= ham[deviceSlice, leftSlice2]
	hD = ham[deviceSlice, deviceSlice]
	vDR= ham[rightSlice2, deviceSlice]
	hR = ham[rightSlice1, rightSlice1]
	vR = ham[rightSlice2, rightSlice1]
	
	# Calculate Green's functions
	gL = GreenFunctionIteration(hL, vL, en)
	gR = GreenFunctionIteration(hR, dagger(vR), en)
	
	# Calculate self energies, robust to varying matrix sizes
	seL = (vDL) @ gL @ dagger(vDL)
	seR = dagger(vDR) @ gR @ (vDR)
	
	# Calculate Gamma functions
	leftGamma = -2*seL.imag
	rightGamma = -2*seR.imag
	return hD, rightGamma, leftGamma, seR, seL
	
def GreenFunctionIteration(h, v, energy, eta = 0.01):
	# input: h, v er 2D arrays, energy is 1D array
	# output: G is 3D array with energy as the first axis
	z = np.stack([(e+eta*1j)*np.eye(len(h)) for e in energy])
	bj = v
	aj = dagger(v)
	ej = h
	esj = h
	while np.max(np.abs(aj))>1e-6:
		gj = np.linalg.inv(np.subtract(z,ej))
		ag = aj @ gj
		bg = bj @ gj
		ej = ej + ag @ bj + bg @ aj
		esj = esj + ag @ bj
		aj = ag @ aj
		bj = bg @ bj

	G = np.linalg.inv(np.subtract(z,esj))
	return G

def dagger(Matrix):
	# Hermitian conjugate of matrix, and of 2 last dimensions for higher dim arrays
	if len(Matrix.shape) > 2:
		return np.transpose(np.conj(Matrix), axes = np.append([n for n in range(len(Matrix.shape)-2)], [len(Matrix.shape), len(Matrix.shape)-1]))
	else:
		return np.transpose(np.conj(Matrix))

class Calculation:
	def __init__(self, inputStr):
		self.hamiltonian = None
		self.molecule = None
		self.leftGamma = None
		self.rightGamma = None
		self.inputStr = inputStr
		self.transmission = np.array([])
		self.energy = np.array([])
		self.eMin, self.eMax, self.eSteps = 0,0,0
		self.noMolecule = False
		self.tiles = int

	def setupMatrices(self, tiles = 16, showStructure = False):
		print("Decoding input data")
		self.eMin, self.eMax, self.eSteps, elementList, positionList, leadList, neighbourList, error = inputDecoder(self.inputStr)
		self.tiles = tiles
		if error:
			print('Exited inputDecoder with error: \n{}\nInput should be "Str float float float int (int ..)"'.format(error))
			return True

	    # Define empty molecule object
		allAtoms = Molecule()
	    
	    # Define atom objects from input string and add to the molecule
		print("Creating atoms and molecule")
		try:
			for i in range(len(elementList)):
				atom = Atom(elementList[i], positionList[i], neighbourList[i], leadList[i])
				allAtoms.addAtom(atom)
		except Exception as e:
			print('Creating the Atom objects failed with message: \n{}'.format(e))
			return True
		
	    # The data structure of neighbour index is put into the corresponding neighbour atom object
		print("Registering atom neighbours")
		try:
			allAtoms.updateNeighbourAtoms()
		except Exception as e:
			print('updateNeighbourAtoms failed with message: \n{}'.format(e))
			return True
		
		# Do Depth-first search of device part of molecule to check for clusters seperated by leads
		print("Clustering molecule islands")
		try:
			clusterList = clusterSearch(allAtoms)
		except Exception as e:
			print('ClusterMerge failed with message: \n{}'.format(e))
			return True
		
		# In case leads are connected to multiple valid devices, merges the molecules, otherwise only update neighbours
		print("Merging valid clusters")
		try:
			self.molecule = clusterMerge(clusterList)
		except Exception as e:
			print('ClusterMerge failed with message: \n{}'.format(e))
			return True
		
		# NoMolecule exception assumes leads directly connected and returns 1 for all energies
		if not self.molecule:
			self.noMolecule = True
			return False
		
	    # Calculate mean Plane model of the molecule.
		print("Fitting the PCA plane")
		try:
			planeModel = PlanefitPCA(self.molecule)
		except Exception as e:
			print('The principle component analysis failed with message: \n{}'.format(e))
			return True

	    # Calculate the direction in which the leads will be tiled
		print("Calculating the direction of the leads")
		try:
			rightDir, rLeadEdge = self.molecule.periodicityDirection(1, planeModel)
			leftDir, lLeadEdge = self.molecule.periodicityDirection(-1, planeModel)
		except Exception as e:
			print('The periodicity direction could not be calculated with message: \n{}'.format(e))
			return True
		
	    # Tile the leads in the given direction
		print("Tiling the leads")		
		try:
			self.molecule.leadCopy(-1, leftDir, planeModel, self.tiles)
			self.molecule.leadCopy(1, rightDir, planeModel, self.tiles)
		except Exception as e:
			print('The lead could not be copied with message: \n{}'.format(e))
			return True
		
	    # Run sort to properly sort all the new atoms added
		try:
			self.molecule.sortMolecule(leftDir, rightDir)
		except Exception as e:
			print('The molecule could not be sorted with message: \n{}'.format(e))
			return True
		
		#print(self.molecule.getPositionList())
		if showStructure:
			plotMeanPlane(self.molecule.getPositionList(), planeModel)
		return False

	def calcTransport(self, data=None):
		# we use self.* value if none is specified.
		if data==None:
			eMin = self.eMin
			eMax = self.eMax
			eSteps = self.eSteps
		else:
			eMin = data[0]
			eMax = data[1]
			eSteps = data[2]
		
		enList = np.linspace(float(eMin),float(eMax), int(eSteps))
		if len(enList) == 0:
			print("No energy points were given in calcTransport")
			return True
		
		#if no valid molecule is present, it is only made of leads due to logic in ue4.
		if self.noMolecule:
			print('Lead to lead transmission is enforced')
			self.energy = enList
			self.transmission = np.ones(len(enList))
			return False
		
		#Logic to ensure refinement doesn't calculate more points than necessary
		self.energy = [en for en in enList if round(en,3) not in np.round(self.energy,3)]
		
		start = time.time()
		global cap
		if cap:
			try:
				self.hamiltonian, self.rightGamma, self.leftGamma, seRP, seLP = capHamiltonian(self.molecule, self.energy)
			except Exception as e:
				print('Creating CAP Hamiltonian failed with message:\n{}'.format(e))
				return True
		else:
			try:
				self.hamiltonian, self.rightGamma, self.leftGamma, seRP, seLP = recurveHamiltonian(self.molecule, self.tiles, self.energy)
			except Exception as e:
				print('Creating recurve Hamiltonian failed with message:\n{}'.format(e))
				return True
		
		# Calculate the Transmission using self.hamiltonian, rightGamma and leftGamma
		print('Calculating transport') 
		h1, h2 = self.hamiltonian.shape
		if int(eSteps) * h1 * h2 > 3.8e7:
			global GPU
			GPU = False

		if GPU:
			self.transmission = np.empty_like(self.energy, dtype=np.complex128)
			parts = int(int(eSteps)/25) # Parts ensure 25 energies per batch
			mul = int(len(self.transmission)/parts)
			try:
				# Calculation using GPU is cut into smaller batches for efficiency and memory concerns
				for i in range(parts):
					if len(self.rightGamma.shape) == 3:
						rightGamma = self.rightGamma[i*mul:i*mul+mul,:,:]
						leftGamma = self.leftGamma[i*mul:i*mul+mul,:,:]
						seR = seRP[i*mul:i*mul+mul,:,:]
						seL = seLP[i*mul:i*mul+mul,:,:]
					else:
						rightGamma = self.rightGamma
						leftGamma = self.leftGamma
						seR = seRP
						seL = seLP
						
					self.transmission[i*mul:i*mul+mul] = transportCalculation_GPU(self.hamiltonian, rightGamma, leftGamma, self.energy[i*mul:i*mul+mul], seR, seL, eta = 0)	
			except Exception as e:
				print('GPU calculations failed... with message:\n{}\nTrying CPU calculation'.format(e))
				try:
					self.transmission = transportCalculation(self.hamiltonian, self.rightGamma, self.leftGamma, self.energy, seR, seL, eta = 0)
				except Exception as e:
					print('CPU calculations also failed, with message:\n{}'.format(e))
					return True
				else:
					print('CPU calculations succeeded.')
			else:
				print('GPU calculations succeeded.')
		else:
			try:
				self.transmission = transportCalculation(self.hamiltonian,self.rightGamma, self.leftGamma, self.energy, seRP, seLP, eta = 0)
				print('CPU calculations succeeded.')
			except Exception as e:
				print('CPU calculations failed with message:\n{}'.format(e))
				return True
		print('Calculation took: {} seconds'.format(time.time()-start))
		
		
		# Check for numerical zero imaginary parts of the transmission
		if all(abs(self.transmission.imag) < 10**(-5)):
			self.transmission = np.round(self.transmission.real,3)
		else:
			print('Something went wrong the transmission is imaginary')
			print(abs(self.transmission.imag))
			return True
			
		# Sort the energy and transmission with respect to self.energy
		try:
			self.energy, self.transmission = zip(*sorted(zip(self.energy, self.transmission)))
			self.energy = np.array(self.energy)
			self.transmission = np.array(self.transmission)
		except Exception as e:
			print('Sorting of transmission failed with message: {}'.format(e))
			return True

		return False

	

if __name__ == "__main__":
	def Plot(x,y,fig = None):
		# Expects a list of x- and y-lists.
		if fig == None:
			fig = plt.figure(figsize=(5, 5), dpi=216,facecolor='w', edgecolor='w')
		style = ['b-', 'r--', 'g-.']
		idx = 0
		for en, tm in zip(x,y):
			plt.plot(en, tm, style[idx], linewidth=1.5)
			idx += 1
		plt.legend(['Recursive','CAP','DFT'], loc = 'lower right')
		plt.xlabel('Energy (eV)', fontsize=12)
		plt.ylabel('Transmission', fontsize=12)
		plt.xticks(fontsize=10)
		plt.yticks(fontsize=10)
		plt.xlim(min(x[0]), max(x[0]))
		plt.ylim(min(y[0])-0.05*max(y[0]), max(y[0])*1.05)
		plt.gca().set_aspect('auto', adjustable='box')
		plt.grid(which='both')
		plt.minorticks_on()
		return fig
	
	def SaveImage(name, fig, x= None, y = None):
		from pathlib import Path
		#outpath = Path.cwd().parent / 'Images'
		outpath = Path.home() / 'Pictures' / 'DeViNa'
		if not outpath.is_dir():
			print('Making new directory..')
			outpath.mkdir()
	
		imgCount = len([path for path in list(outpath.glob('*')) if path.suffix == '.png']) + 1
		
		try:
			string = str(imgCount) + name
			if type(x) == np.ndarray and type(y) == np.ndarray:
				np.savez(outpath / string, x=x, y=y)
			fig.savefig(outpath / string)
			plt.close(fig)
		except Exception as e:
			print('Filename is invalid with message:\n{}'.format(e))
			
			
	# Lead to lead test
	testInputStr0 = '''-7.0 7.0 50 C 1.73 1 0 1 1 C 0.866 1.5 0 -1 0'''.split()
	testInputStr01 = '''-7.0 7.0 50 C 1.73 1 0 1 1 C 0.866 1.5 0 -1 0 C 0.866 2.5 0 0'''.split()
	
	# With and Without Fluor Test
	testInputStr1 = '''-7.0 7.0 50 C -6.62 -12.20 68.81 1 1 C -5.12 -12.32 68.96 0 0 2 3 C -3.67 -12.04 69.26 -1 1 F -5.12 -13.62 68.96 0 1'''.split()
	testInputStr2 = '''-7.0 7.0 50 C -6.62 -12.20 68.81 1 1 C -5.12 -12.32 68.96 0 0 2 C -3.67 -12.04 69.26 -1 1'''.split()
	
	# Nanoribbon tests
	testInputStr3 = '''-7.0 7.0 50 C -1.24 -10.41 59.60 1 1 9 C -2.64 -10.94 59.65 1 3 0 C 1.56 -13.69 59.69 -1 6 8 4 C -2.89 -12.44 59.76 1 1 7 C 2.76 -14.71 60.06 -1 2 C -0.68 -13.09 59.67 0 7 9 6 C 0.19 -14.11 59.88 0 5 2 C -2.03 -13.62 59.89 1 3 5 C 2.02 -12.37 59.49 -1 2 10 11 C -0.24 -11.64 59.63 0 0 5 10 C 1.14 -11.27 59.63 0 9 8 C 3.47 -11.98 59.37 -1 8'''.split()
	testInputStr4 = '''-7.0 7.0 50 C -1.36 -9.15 68.08 -1 1 C -2.56 -9.89 68.75 -1 2 4 0 C -2.53 -10.98 69.92 -1 7 1 3 C -1.11 -11.45 70.66 -1 2 C -4.12 -9.26 68.35 0 5 1 C -5.07 -10.18 69.20 0 4 6 8 C -4.79 -11.07 70.22 0 5 7 11 C -3.50 -11.74 70.69 0 2 6 C -6.54 -9.63 68.61 1 5 9 C -7.61 -10.15 69.63 1 10 8 C -7.37 -11.34 70.66 1 9 11 C -5.84 -11.91 71.14 1 6 10'''.split()

	# Larger nanoribbons (#9 is largest structure)
	testInputStr5 = '''-5.0 5.0 50 C -11.73 -1.55 57.89 1 1 19 C -10.94 -2.88 57.80 1 2 9 0 C -9.29 -2.74 57.68 1 1 3 23 C -8.58 -1.60 57.72 1 2 C -13.80 -5.48 58.20 0 5 6 15 C -13.10 -4.19 58.04 0 4 9 20 C -15.17 -5.55 58.53 0 4 14 18 C -9.41 -5.41 57.52 0 8 13 23 C -10.91 -5.38 57.71 0 7 9 10 C -11.65 -4.18 57.84 0 1 5 8 C -11.66 -6.64 57.78 -1 8 11 15 C -10.90 -7.74 57.63 -1 10 12 C -9.45 -7.92 57.41 -1 11 13 C -8.58 -6.63 57.34 -1 7 12 C -15.91 -4.39 58.67 0 6 21 C -13.04 -6.67 58.03 -1 4 10 16 C -13.57 -8.12 58.12 -1 15 17 C -15.00 -7.97 58.51 -1 16 18 C -15.97 -6.99 58.78 -1 6 17 C -13.17 -1.56 58.00 1 20 0 C -13.80 -2.99 58.14 1 5 19 21 C -15.23 -2.85 58.45 1 14 20 22 C -15.85 -1.73 58.57 1 21 C -8.66 -4.05 57.46 0 2 7'''.split()
	testInputStr6 = '''-5.0 5.0 50 C -11.73 -1.55 57.89 1 1 19 C -10.94 -2.88 57.80 1 2 9 0 C -9.29 -2.74 57.68 1 1 3 23 C -8.58 -1.60 57.72 1 2 C -13.80 -5.48 58.20 0 5 6 15 C -13.10 -4.19 58.04 0 4 9 20 C -15.17 -5.55 58.53 0 4 14 18 C -9.41 -5.41 57.52 0 8 13 23 C -10.91 -5.38 57.71 0 7 9 10 C -11.65 -4.18 57.84 0 1 5 8 C -11.66 -6.64 57.78 -1 8 11 15 C -10.90 -7.74 57.63 -1 10 12 C -9.45 -7.92 57.41 -1 11 13 C -8.58 -6.63 57.34 -1 7 12 C -15.91 -4.39 58.67 0 6 21 C -13.04 -6.67 58.03 -1 4 10 16 C -13.57 -8.12 58.12 -1 15 17 C -15.00 -7.97 58.51 -1 16 18 C -15.97 -6.99 58.78 -1 6 17 C -13.17 -1.56 58.00 1 20 0 C -13.80 -2.99 58.14 1 5 19 21 C -15.23 -2.85 58.45 1 14 20 22 C -15.85 -1.73 58.57 1 21 C -8.66 -4.05 57.46 0 2 7'''.split()
	testInputStr7 = '''-5.0 5.0 50 C -14.15 -1.45 61.73 1 1 20 C -13.95 -2.51 61.05 1 0 2 10 C -13.51 -2.51 59.46 1 1 3 5 C -13.24 -1.44 58.69 1 2 C -14.52 -4.82 64.11 0 6 7 16 C -13.51 -3.92 59.23 0 2 8 C -14.44 -3.61 63.29 0 4 10 21 C -14.66 -4.76 65.51 0 4 15 19 C -13.69 -5.17 59.77 0 5 9 14 C -14.02 -5.00 61.23 0 8 10 11 C -14.15 -3.74 61.85 0 1 6 9 C -14.19 -6.18 62.07 -1 9 12 16 C -14.07 -7.35 61.42 -1 11 13 C -13.80 -7.67 60.00 -1 12 14 C -13.58 -6.47 59.04 -1 8 13 C -14.78 -3.54 66.15 0 7 22 C -14.42 -6.08 63.46 -1 4 11 17 C -14.51 -7.47 64.12 -1 16 18 C -14.62 -7.19 65.58 -1 17 19 C -14.69 -6.12 66.49 -1 7 18 C -14.51 -0.98 63.08 1 0 21 C -14.58 -2.35 63.86 1 6 20 22 C -14.77 -2.08 65.30 1 15 21 23 C -14.87 -0.90 65.81 1 22'''.split()
	testInputStr8 = '''-5.0 5.0 50 C -9.28 -6.12 60.19 0 5 9 11 C -7.28 -4.80 60.07 1 6 9 C -9.55 -3.62 60.07 1 6 5 C -11.57 -5.08 60.19 1 4 5 13 C -12.23 -3.84 60.15 1 7 3 C -10.16 -4.97 60.15 1 2 3 0 C -8.03 -3.51 60.03 1 2 1 C -13.79 -3.94 60.19 1 4 8 C -14.40 -5.28 60.28 1 7 12 C -7.93 -6.01 60.15 0 0 1 10 C -7.11 -7.20 60.19 0 9 17 C -10.00 -7.35 60.27 0 0 14 16 C -13.48 -6.46 60.31 0 8 13 15 C -12.09 -6.30 60.27 0 3 12 14 C -11.40 -7.49 60.31 0 11 13 21 C -14.22 -7.67 60.39 0 12 20 C -9.14 -8.48 60.30 0 11 17 19 C -7.75 -8.43 60.26 0 10 16 18 C -6.94 -9.62 60.30 0 17 25 C -9.78 -9.73 60.38 0 16 22 24 C -13.32 -8.84 60.43 0 15 21 23 C -11.90 -8.73 60.39 0 14 20 22 C -11.19 -9.93 60.43 0 19 21 30 C -14.00 -10.07 60.51 0 20 28 C -8.98 -10.89 60.42 0 19 25 27 C -7.57 -10.88 60.38 0 18 24 26 C -6.76 -12.06 60.42 0 25 33 C -9.58 -12.15 60.50 0 24 29 32 C -13.14 -11.27 60.54 0 23 30 31 C -11.00 -12.38 60.54 0 27 30 37 C -11.72 -11.18 60.50 0 22 28 29 C -13.80 -12.51 60.62 0 28 36 C -8.81 -13.32 60.53 0 27 33 35 C -7.39 -13.33 60.50 0 26 32 34 C -6.58 -14.51 60.54 0 33 41 C -9.38 -14.58 60.61 0 32 38 40 C -12.96 -13.71 60.66 0 31 37 39 C -11.53 -13.62 60.62 0 29 36 38 C -10.81 -14.83 60.66 0 35 37 45 C -13.62 -14.95 60.74 0 36 44 C -8.63 -15.76 60.65 0 35 41 43 C -7.21 -15.77 60.62 0 34 40 42 C -6.40 -16.95 60.65 0 41 49 C -9.20 -17.03 60.73 0 40 46 48 C -12.78 -16.16 60.78 0 39 45 47 C -11.35 -16.07 60.74 0 38 44 46 C -10.63 -17.28 60.78 0 43 45 53 C -13.43 -17.40 60.86 0 44 52 C -8.45 -18.21 60.77 -1 43 49 51 C -7.03 -18.22 60.73 -1 42 48 50 C -6.21 -19.40 60.77 -1 49 C -9.01 -19.47 60.85 -1 48 54 C -12.59 -18.61 60.90 -1 47 53 55 C -11.17 -18.52 60.85 -1 46 52 54 C -10.44 -19.73 60.90 -1 51 53 C -13.25 -19.85 60.97 -1 52'''.split()
	testInputStr9 = '''-7.0 7.0 50 C 5.37 -4.06 59.04 -1 1 3 C 5.31 -5.58 58.97 -1 0 2 C 3.93 -6.17 58.93 -1 1 13 C 4.11 -3.29 59.07 -1 0 6 12 C 5.56 0.18 59.24 -1 5 7 C 5.48 -1.37 59.17 -1 4 6 C 4.16 -1.89 59.13 -1 3 5 17 C 4.30 0.94 59.26 -1 4 10 16 C 5.75 4.43 59.44 -1 9 11 C 5.67 2.87 59.37 -1 8 10 C 4.35 2.35 59.33 -1 7 9 21 C 4.48 5.19 59.46 -1 8 20 C 2.86 -4.02 59.02 0 3 13 15 C 2.81 -5.39 58.96 0 2 12 14 C 1.53 -6.05 58.92 0 13 25 C 1.72 -3.17 59.05 0 12 18 24 C 3.04 0.20 59.22 0 7 17 19 C 3.01 -1.19 59.15 0 6 16 18 C 1.76 -1.76 59.12 0 15 17 29 C 1.92 1.06 59.25 0 16 22 28 C 3.21 4.42 59.41 0 11 21 23 C 3.20 3.02 59.35 0 10 20 22 C 1.94 2.47 59.31 0 19 21 33 C 2.09 5.30 59.45 0 20 32 C 0.49 -3.88 59.01 0 15 25 27 C 0.38 -5.27 58.94 0 14 24 26 C -0.88 -5.93 58.90 0 25 37 C -0.66 -3.11 59.03 0 24 30 36 C 0.69 0.34 59.20 0 19 29 31 C 0.57 -1.08 59.14 0 18 28 30 C -0.68 -1.67 59.10 0 27 29 41 C -0.45 1.12 59.23 0 28 34 40 C 0.82 4.54 59.40 0 23 33 35 C 0.76 3.12 59.33 0 22 32 34 C -0.50 2.56 59.30 0 31 33 45 C -0.31 5.36 59.43 0 32 44 C -1.91 -3.75 58.99 0 27 37 39 C -2.06 -5.15 58.93 0 26 36 38 C -3.33 -5.82 58.88 0 37 49 C -3.09 -3.02 59.01 0 36 42 48 C -1.71 0.47 59.19 0 31 41 43 C -1.86 -0.96 59.12 0 30 40 42 C -3.13 -1.57 59.08 0 39 41 53 C -2.88 1.20 59.21 0 40 46 52 C -1.60 4.66 59.38 0 35 45 47 C -1.68 3.23 59.32 0 34 44 46 C -2.96 2.66 59.28 0 43 45 57 C -2.76 5.46 59.41 0 44 56 C -4.34 -3.64 58.97 1 39 49 51 C -4.51 -5.04 58.91 1 38 48 50 C -5.78 -5.71 58.87 1 49 C -5.53 -2.91 59.00 1 48 54 C -4.14 0.58 59.17 1 43 53 55 C -4.32 -0.85 59.10 1 42 52 54 C -5.58 -1.46 59.06 1 51 53 C -5.31 1.30 59.19 1 52 58 C -4.05 4.77 59.36 1 47 57 59 C -4.13 3.34 59.30 1 46 56 58 C -5.41 2.76 59.26 1 55 57 C -5.21 5.57 59.39 1 56'''.split()
	
	# Cluster isolation test (remove excess garbage)
	testInputStr10 = '''-7.0 7.0 50 C -11.66 1.61 61.35 0 10 1 3 C -11.06 1.09 60.10 0 0 4 3 C -9.28 -0.11 58.16 0 4 C -12.36 1.09 59.52 0 0 1 C -9.89 0.43 59.43 0 1 2 C -6.66 1.57 61.58 1 11 6 C -5.18 0.89 60.42 0 5 7 C -5.35 0.43 59.00 0 6 8 C -6.38 0.51 57.91 0 7 9 C -7.55 1.19 57.29 0 8 C -10.76 1.68 62.52 0 12 0 C -8.06 2.09 61.87 0 5 12 C -9.32 1.86 62.66 -1 11 10'''.split()
    
	# Carbon ring, cluster merging test
	testInputStr11 = '''-7.0 7.0 50 C 1.73 1 0 0 1 2 C 0.866 1.5 0 -1 0 3 C 2.598 1.5 0 1 0 4 C 0.866 2.5 0 -1 1 5 C 2.598 2.5 0 1 2 5 C 1.73 3 0 0 3 4'''.split()

	# index order error (these should both be 1 transmission)
	testInputStr12 = '''-5.0 5.0 50 C -14.20 -4.13 63.94 -1 2 C -14.00 -1.15 63.43 1 2 C -14.04 -2.65 63.63 0 1 0'''.split()
	testInputStr13 = '''-5.0 5.0 50 C -14.20 -4.13 63.94 -1 2 3 C -14.00 -1.15 63.43 1 2 C -14.04 -2.65 63.63 0 1 0 C -13.90 -5.56 64.22 0 0'''.split()
	
	#2 atom leads, 1 atom device:
	#90 Deg of one lead:
	testInputStr14 = ''' -5.0 5.0 50 C -4.96 0.41 62.73 0 3 2 C -2.52 -1.30 62.41 -1 3 C -6.26 0.94 63.26 1 0 4 C -3.72 -0.45 62.77 -1 1 0 C -7.60 1.52 63.66 1 2'''.split()

	#larger than 90 Deg of one lead:
	testInputStr15 = ''' -5.0 5.0 50 C -4.96 0.41 62.73 0 3 2 C -2.52 -1.30 62.41 -1 3 C -6.26 0.94 63.26 1 0 4 C -3.72 -0.45 62.77 -1 1 0 C -6.95 0.91 61.75 1 2'''.split()

	testInputStr16 = ''' -5.0 5.0 50 C -4.96 0.41 62.73 0 3 2 C -2.52 -1.30 62.41 -1 3 C -6.26 0.94 63.26 1 0 4 C -3.72 -0.45 62.77 -1 1 0 C -7.59 1.20 62.22 1 2'''.split()

	#0 transmission (full length):
	testInputStr17 = '''-5.0 5.0 50 C -4.96 0.41 62.73 0 3 2 C -2.52 -1.30 62.41 -1 3 C -6.26 0.94 63.26 1 0 4 C -3.72 -0.45 62.77 -1 1 0 C -7.75 1.16 61.93 1 2'''.split()
    
	#full length:
	testInputStr18 = '''-5.0 5.0 50 C -12.66 4.33 64.03 0 2 1 C -14.53 3.57 64.20 -1 0 C -10.72 4.89 63.83 0 0 3 C -8.73 5.23 63.53 1 2'''.split()
	testInputStr19 = '''-5.0 5.0 50 C -12.66 4.33 64.03 0 2 1 C -14.53 3.57 64.20 -1 0 C -10.72 4.89 63.83 1 0 3 C -8.73 5.23 63.53 1 2'''.split()
	
	# Carbonhagen tests:
	testInputStr20 = '''-2.6 2.6 100 C 4.40 3.34 61.51 -1 1 3 C 4.53 1.82 61.50 -1 0 2 C 3.24 1.07 61.51 -1 1 9 C 3.05 3.94 61.53 -1 0 6 8 C 4.06 7.58 61.54 -1 5 7 C 4.17 6.02 61.53 -1 4 6 C 2.93 5.35 61.54 -1 3 5 13 C 2.71 8.18 61.56 -1 4 12 C 1.90 3.06 61.54 0 3 9 11 C 2.03 1.70 61.53 0 2 8 10 C 0.84 0.88 61.53 0 9 17 C 0.67 3.76 61.55 0 8 14 16 C 1.54 7.26 61.57 0 7 13 15 C 1.71 5.87 61.56 0 6 12 14 C 0.52 5.16 61.56 0 11 13 21 C 0.32 7.99 61.58 0 12 20 C -0.45 2.90 61.56 0 11 17 19 C -0.39 1.51 61.55 0 10 16 18 C -1.57 0.70 61.56 0 17 25 C -1.70 3.53 61.58 0 16 22 24 C -0.83 7.08 61.59 0 15 21 23 C -0.72 5.66 61.58 0 14 20 22 C -1.91 4.95 61.59 0 19 21 29 C -2.07 7.75 61.61 0 20 28 C -2.86 2.73 61.59 0 19 25 27 C -2.84 1.32 61.58 0 18 24 26 C -4.01 0.51 61.58 0 25 33 C -4.12 3.32 61.60 0 24 30 32 C -3.26 6.89 61.62 0 23 29 31 C -3.16 5.47 61.61 0 22 28 30 C -4.36 4.74 61.61 0 27 29 37 C -4.51 7.54 61.63 0 28 36 C -5.29 2.54 61.61 0 27 33 35 C -5.28 1.13 61.60 0 26 32 34 C -6.46 0.31 61.61 0 33 41 C -6.56 3.11 61.63 0 32 38 40 C -5.71 6.70 61.64 0 31 37 39 C -5.61 5.27 61.63 0 30 36 38 C -6.81 4.54 61.64 0 35 37 45 C -6.96 7.34 61.66 0 36 44 C -7.73 2.35 61.64 1 35 41 43 C -7.73 0.94 61.63 1 34 40 42 C -8.91 0.11 61.63 1 41 C -9.00 2.91 61.65 1 40 46 C -8.16 6.50 61.67 1 39 45 47 C -8.06 5.07 61.66 1 38 44 46 C -9.26 4.34 61.66 1 43 45 C -9.41 7.15 61.68 1 44'''.split()
	testInputStr21 = '''-2.6 2.6 100 C 4.40 3.34 61.51 -1 1 3 C 4.53 1.82 61.50 -1 0 2 C 3.24 1.07 61.51 -1 1 9 C 3.05 3.94 61.53 -1 0 6 8 C 4.06 7.58 61.54 -1 5 7 C 4.17 6.02 61.53 -1 4 6 C 2.93 5.35 61.54 -1 3 5 13 C 2.71 8.18 61.56 -1 4 12 C 1.90 3.06 61.54 0 3 9 11 C 2.03 1.70 61.53 0 2 8 10 C 0.84 0.88 61.53 0 9 17 C 0.67 3.76 61.55 0 8 14 16 C 1.54 7.26 61.57 0 7 13 15 C 1.71 5.87 61.56 0 6 12 14 C 0.52 5.16 61.56 0 11 13 21 C 0.32 7.99 61.58 0 12 20 C -0.35 2.88 61.43 0 11 17 19 C -0.39 1.51 61.55 0 10 18 16 C -1.57 0.70 61.56 0 17 25 C -1.74 3.58 62.16 0 48 22 16 24 C -0.83 7.08 61.59 0 15 21 23 C -0.72 5.66 61.58 0 14 20 22 C -1.96 4.99 61.45 0 19 21 29 C -2.07 7.75 61.61 0 20 28 C -2.82 2.66 61.44 0 19 25 27 C -2.84 1.32 61.58 0 18 26 24 C -4.01 0.51 61.58 0 25 33 C -4.12 3.32 61.60 0 30 32 24 C -3.26 6.89 61.62 0 23 29 31 C -3.16 5.47 61.61 0 28 30 22 C -4.36 4.74 61.61 0 27 29 37 C -4.51 7.54 61.63 0 28 36 C -5.29 2.54 61.61 0 27 33 35 C -5.28 1.13 61.60 0 26 32 34 C -6.46 0.31 61.61 0 33 41 C -6.56 3.11 61.63 0 32 38 40 C -5.71 6.70 61.64 0 31 37 39 C -5.61 5.27 61.63 0 30 36 38 C -6.81 4.54 61.64 0 35 37 45 C -6.96 7.34 61.66 0 36 44 C -7.73 2.35 61.64 1 35 41 43 C -7.73 0.94 61.63 1 34 40 42 C -8.91 0.11 61.63 1 41 C -9.00 2.91 61.65 1 40 46 C -8.16 6.50 61.67 1 39 45 47 C -8.06 5.07 61.66 1 38 44 46 C -9.26 4.34 61.66 1 43 45 C -9.41 7.15 61.68 1 44 F -1.76 3.47 63.52 0 19'''.split()
	# DTF data 1 fluor
	energy = '''[-1.98000, -1.94000, -1.90000, -1.86000, -1.82000, -1.78000, -1.74000, -1.70000, -1.66000, -1.62000, -1.58000, -1.54000, -1.50000, -1.46000, -1.42000, -1.38000, -1.34000, -1.30000, -1.26000, -1.22000, -1.18000, -1.14000, -1.10000, -1.06000, -1.02000, -0.98000, -0.94000, -0.90000, -0.86000, -0.82000, -0.78000, -0.74000, -0.70000, -0.66000, -0.62000, -0.58000, -0.54000, -0.50000, -0.46000, -0.42000, -0.38000, -0.34000, -0.30000, -0.26000, -0.22000, -0.18000, -0.14000, -0.10000, -0.06000, -0.02000,  0.02000, 0.06000, 0.10000, 0.14000, 0.18000, 0.22000, 0.26000, 0.30000, 0.34000, 0.38000, 0.42000, 0.46000, 0.50000, 0.54000, 0.58000, 0.62000, 0.66000, 0.70000, 0.74000, 0.78000, 0.82000, 0.86000, 0.90000, 0.94000, 0.98000, 1.02000, 1.06000, 1.10000, 1.14000, 1.18000, 1.22000, 1.26000, 1.30000, 1.34000, 1.38000, 1.42000, 1.46000, 1.50000, 1.54000, 1.58000, 1.62000, 1.66000, 1.70000, 1.74000, 1.78000, 1.82000, 1.86000, 1.90000, 1.94000, 1.98000]'''
	trans = '''[0.17940512E+01 0.17481945E+01 0.16365986E+01 0.83681250E+00 0.96348977E+00 0.98676282E+00 0.99797022E+00 0.99963760E+00 0.99421054E+00 0.98336381E+00 0.96827465E+00 0.94978911E+00 0.92852050E+00 0.90491265E+00 0.87928247E+00 0.85185122E+00 0.82276690E+00 0.79212213E+00 0.75996774E+00 0.72632456E+00 0.69119316E+00 0.65456313E+00 0.61642230E+00 0.57676625E+00 0.53560913E+00 0.49299595E+00 0.44901738E+00 0.40382674E+00 0.35766026E+00 0.31086046E+00 0.26390171E+00 0.21741791E+00 0.17222886E+00 0.12936281E+00 0.90069525E-01 0.55817436E-01 0.28267233E-01 0.92151817E-02 0.50218950E-03 0.38912978E-02 0.20925997E-01 0.52794170E-01 0.10023560E+00 0.16354537E+00 0.24275179E+00 0.33813709E+00 0.45168680E+00 0.59290171E+00 0.89097857E+00 0.33747137E+01 0.32676508E+01 0.49114135E+00 0.93080294E+00 0.94610500E+00 0.94939321E+00 0.95019549E+00 0.95027733E+00 0.95015347E+00 0.95000440E+00 0.94989729E+00 0.94985503E+00 0.94988245E+00 0.94997656E+00 0.95013160E+00 0.95034069E+00 0.95059723E+00 0.95089513E+00 0.95122880E+00 0.95159352E+00 0.95198518E+00 0.95240045E+00 0.95283639E+00 0.95329070E+00 0.95376152E+00 0.95424759E+00 0.95474786E+00 0.95526189E+00 0.95578957E+00 0.95633119E+00 0.95688754E+00 0.95745993E+00 0.95805007E+00 0.95866060E+00 0.95929462E+00 0.95995653E+00 0.96065164E+00 0.96138716E+00 0.96217227E+00 0.96301913E+00 0.96394420E+00 0.96496999E+00 0.96612865E+00 0.96746790E+00 0.96906334E+00 0.97104633E+00 0.97368282E+00 0.97771794E+00 0.17218781E+01 0.19364297E+01 0.19709136E+01]'''
	trans= trans.replace(' ', ',')
	en1 = eval(energy)
	tm1 = eval(trans)
	# DTF data pristine ribbon
	energy2 = '''[-1.98000 -1.94000 -1.90000 -1.86000 -1.82000 -1.78000 -1.74000 -1.70000 -1.66000 -1.62000 -1.58000 -1.54000 -1.50000 -1.46000 -1.42000 -1.38000 -1.34000 -1.30000 -1.26000 -1.22000 -1.18000 -1.14000 -1.10000 -1.06000 -1.02000 -0.98000 -0.94000 -0.90000 -0.86000 -0.82000 -0.78000 -0.74000 -0.70000 -0.66000 -0.62000 -0.58000 -0.54000 -0.50000 -0.46000 -0.42000 -0.38000 -0.34000 -0.30000 -0.26000 -0.22000 -0.18000 -0.14000 -0.10000 -0.06000 -0.02000 0.02000 0.06000 0.10000 0.14000 0.18000 0.22000 0.26000 0.30000 0.34000 0.38000 0.42000 0.46000 0.50000 0.54000 0.58000 0.62000 0.66000 0.70000 0.74000 0.78000 0.82000 0.86000 0.90000 0.94000 0.98000 1.02000 1.06000 1.10000 1.14000 1.18000 1.22000 1.26000 1.30000 1.34000 1.38000 1.42000 1.46000 1.50000 1.54000 1.58000 1.62000 1.66000 1.70000 1.74000 1.78000 1.82000 1.86000 1.90000 1.94000 1.98000]'''
	energy2= energy2.replace(' ', ',')
	en2 = eval(energy2)
	transmission2 ='''[0.29999998E+01 0.29999993E+01 0.29999967E+01 0.10150543E+01 0.10000017E+01 0.10000004E+01 0.10000001E+01 0.10000001E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000002E+01 0.10000010E+01 0.10015843E+01 0.39999831E+01 0.39996142E+01 0.10000188E+01 0.10000008E+01 0.10000002E+01 0.10000001E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000000E+01 0.10000001E+01 0.10000001E+01 0.10000005E+01 0.10000032E+01 0.29999511E+01 0.29999988E+01 0.29999998E+01]'''
	transmission2= transmission2.replace(' ', ',')
	tm2 = eval(transmission2)
	
	# Multi input
	inputList = [testInputStr1, testInputStr2, testInputStr3, testInputStr4,
			  testInputStr5, testInputStr6, testInputStr7, testInputStr8,
			  testInputStr9, testInputStr10, testInputStr11, testInputStr12,
			  testInputStr13, testInputStr14, testInputStr15, testInputStr16,
			  testInputStr17, testInputStr18, testInputStr19]
	
	#Single input
	inputList = [testInputStr21]
	
	cap = False
	energies = np.zeros(int(inputList[0][2]))
	transmission = np.zeros(int(inputList[0][2]))
	atimes = []
	aconvergence = []
	aatoms = []
	ainputAtoms = []
	for idx, testInputStr in enumerate(inputList):
		for i in range(1):
			start = time.time()
			calc = Calculation(testInputStr)
			#error = calc.setupMatrices(tiles = (i+1))
			error = calc.setupMatrices(tiles = 25, showStructure = False)
			if not error:
				
				error = calc.calcTransport()
				
				if error:
					print('An error has occured in the transport calculation')
					break
			else:
				print('An error has occured in matrix setup')
				break
			
			print('Full program execution time: {} seconds'.format(time.time()-start))
			transmission = np.vstack((transmission, np.array(calc.transmission)))
			energies = np.vstack((energies, np.array(calc.energy)))
			atimes.append(time.time()-start)
			fig = Plot([calc.energy], [calc.transmission])
			SaveImage('Transport', fig, calc.energy, calc.transmission)
			#plt.close(fig)
# =============================================================================
# 			if max(abs(transmission[i,:]-transmission[i+1,:])) < 0.1 and i > 0:
# 				print('Convergence was reached after {} iterations'.format(i))
# 				aconvergence.append(i)
# 				aatoms.append(len(calc.hamiltonian))
# 				ainputAtoms.append(sum([i.isalpha() for i in testInputStr]))
# 				d = [max(abs(transmission[i,:]-transmission[i+1,:])) for i in range(len(transmission)-1)]
# 				break
# =============================================================================
		if not error:	
			if i > 300:
				fig2 = go.Figure()
				fig2.add_trace(go.Surface(x=energies[1:], z=transmission[1:], showscale=False , opacity=0.8, hoverinfo='skip'))
				fig2.update_layout(margin=dict(l=0, r=0, b=0, t=0), scene_aspectmode='data')
				plot(fig2)
			else:
				pass
				plt.figure(idx)
				plt.plot(calc.energy,calc.transmission)
				
			energies = np.zeros(50)
			transmission = np.zeros(50)	
	
	if not error:
		print('The mean time took {}, with min-max range of {} - {}'.format(np.mean(atimes), min(atimes),max(atimes)))
	
# =============================================================================
# 	cap = True
# 	energies = np.zeros(int(inputList[0][2]))
# 	transmission = np.zeros(int(inputList[0][2]))
# 	atimes = []
# 	aconvergence = []
# 	aatoms = []
# 	ainputAtoms = []
# 	for idx, testInputStr in enumerate(inputList):
# 		for i in range(30):
# 			start = time.time()
# 			calc = Calculation(testInputStr)
# 			error = calc.setupMatrices(tiles = (i+1), showStructure = False)
# 			#error = calc.setupMatrices(tiles = 13)
# 			if not error:    
# 				error = calc.calcTransport()
# 				if error:
# 					print('An error has occured in the transport calculation')
# 					break
# 			else:
# 				print('An error has occured in matrix setup')
# 				break
# 			
# 			print('Full program execution time: {} seconds'.format(time.time()-start))
# 			transmission = np.vstack((transmission, np.array(calc.transmission)))
# 			energies = np.vstack((energies, np.array(calc.energy)))
# 			atimes.append(time.time()-start)
# 			#fig = Plot(calc.energy, calc.transmission)
# 			#SaveImage('Transport', fig)
# 			if max(abs(transmission[i,:]-transmission[i+1,:])) < 0.01 and i > 0:
# 				print('Convergence was reached after {} iterations'.format(i))
# 				aconvergence.append(i)
# 				aatoms.append(len(calc.hamiltonian))
# 				ainputAtoms.append(sum([i.isalpha() for i in testInputStr]))
# 				d = [max(abs(transmission[i,:]-transmission[i+1,:])) for i in range(len(transmission)-1)]
# 				break
# 		if not error:	
# 			if i > 0:
# 				fig2 = go.Figure()
# 				fig2.add_trace(go.Surface(x=energies[1:], z=transmission[1:], showscale=False , opacity=0.8, hoverinfo='skip'))
# 				fig2.update_layout(scene = dict(xaxis_title='Energy [eV]', yaxis_title='Lead Repetitions', zaxis_title='Transmission'), margin=dict(l=0, r=0, b=0, t=0), scene_aspectmode='data')
# 				plot(fig2)
# 			else:
# 				plt.figure(idx)
# 				plt.plot(calc.energy,calc.transmission)
# 				
# 			energies = np.zeros(50)
# 			transmission = np.zeros(50)	
# 	
# 	if not error:
# 		print('The mean time took {}, with min-max range of {} - {}'.format(np.mean(atimes), min(atimes),max(atimes)))
# =============================================================================
	
	
	
	
	#fig = go.Figure(go.Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1]))
	#fig.show()
	
	#plotly.io.write_html
	#GPU = GPUtil.getGPUs()[0]
	#cp.get_default_memory_pool().free_all_blocks()
	#cp.get_default_pinned_memory_pool().free_all_blocks()
	#mempool.free_all_blocks()
	#pinned_mempool.free_all_blocks()
	# =============================================================================
# 		print('Total memory: {}'.format(GPU.memoryTotal))
# 		print('Free memory: {}'.format(GPU.memoryFree))
# 		print('Used memory: {}'.format(GPU.memoryUsed))
# 		print('Memory percentage: {}'.format((100-(GPU.memoryTotal-GPU.memoryUsed)/GPU.memoryTotal*100)))
# 		print('Cuda memory size: {}'.format(cp.cuda.Memory.size))
# =============================================================================