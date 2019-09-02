import socket
import ASEscriptOnServer
import Transport
import threading
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

s=None
TCP_PORT=3939
TCP_IP = 'localhost'

#first word tells what the message is about
ECHO = 'ECHO' #debug and conformity check
RUN = 'RUN' #check if relaxation is running	
ERROR = 'ERROR' #signal an error
RELAX = 'RELAX' #start relaxation
STOP = 'STOP' #stop all computation 
TRANSPORT = 'TRANSPORT' #start transport calculation
REFINE_TRANSPORT = 'REFINE_TRANSPORT' #refine transport calculation by adding the specified extra points
TRANSPORT_COMPLETE = 'TRANSPORT_COMPLETE' #send completed transport calculation results
RELAX_COMPLETE = 'RELAX_COMPLETE' #send completed relaxation calculation results
IMAGE = 'IMAGE'
IMAGE_COMPLETE = 'IMAGE_COMPLETE'


#store last sent message to check conformity
lastSentMessage = bytearray()
lastSentMsgType = ECHO

#flag for stopping threads
exitFlag = 0

#################################################################################################
def openConnection(TCP_PORT=3939, TCP_IP = '127.0.0.1'):
	global s
	
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.connect((TCP_IP, TCP_PORT))
    

def sendMessage(msg, msgType = 'ECHO'):
	bMessage = formatMessage(msg, msgType)
	s.sendall(bMessage)

def processData(data):
	#because UE4 TCP implementation uses variable length messages with 3 termination bytes all equal to 255, we check where these bytes are and strip everything afterwards
	i = 1	
	while ( (i < len(data) - 4) & ((data[i] != 0xff) | (data[i + 1] != 0xff) | (data[i + 2] != 0xff)| (data[i + 3] != 0xff))):
		i += 1
	if ( data[i] == 0xff ):
		data = data[0 : i-1]

	#decode data
	strData = data.decode("utf-8")

	#remove potential spaces at the end of the string
	while (strData.endswith(" ")):
		strData = strData[:-1]
	
	#extract msgType
	words = strData.split()
	msgType = words[0]
	
	#remove msgType from data
	strData = strData[len(msgType) + 1 : len(strData) ]
	
	return strData, msgType
	
def receiveMessage():
	#We use a fixed size buffer on python side with too large of a buffer (8kB) and then cut at the termination character
	BUFFER_SIZE = 8192
	data = s.recv(BUFFER_SIZE)	
	
	#transform the byte array into usable data, as well as extract message type flag
	strData, msgType = processData(data) 
	
	print ("client received data:", strData, "with type", msgType)
	return strData, msgType

def formatMessage(message, msgType = 'ECHO'):
	
	while (message.endswith(" ")):
		message = message[:-1]
	
	bMessage = bytearray((msgType + " ").encode('utf-8'))	
	
	#convert the message to bytes and append it to the message type
	bMessage.extend(bytearray(message.encode("utf-8")))
		
	#add the termination bytes
	bMessage.append(0xff)
	bMessage.append(0xff)
	bMessage.append(0xff)
	bMessage.append(0xff)

	return bMessage

def closeConnection():
    s.close()


def checkConformity(data, lastSentMessage, lastSentMsgType):
	if (len(lastSentMessage) > 0):
		if (data == lastSentMessage):
			sendMessage("All Good", ECHO)
		else:
			sendMessage("Not conform", ERROR)
			sendMessage(lastSentMessage, lastSentMsgType)
	else: 
		sendMessage("No previously sent message", ERROR)
	
###############################################################################################
class relaxationThread (threading.Thread):
	def __init__(self, threadID, name, data):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name
		self.data = data
	def run(self):
		print ("Starting " + self.name)
		output, ok = ASEscriptOnServer.GULPRelaxation(data.split())
		if (ok):
			print("Relaxation thread %s complete" %self.name)
			sendMessage(output, RELAX_COMPLETE)
		else:
			print("something went wrong in relaxation")
			sendMessage("something went wrong in relaxation", ERROR)
		

##############################################################################################
class transportThread (threading.Thread):
	def __init__(self, threadID, name, data):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name
		self.data = data
        
	def run(self):
		print ("Starting " + self.name)
		if not refine:
			#Setting up a new calculation object, data is saved in self.calc
			self.calc = Transport.Calculation(self.data.split())
			print("Setting up matrices")
			
			# Modifying structure, repeating leads and making matrices
			error = self.calc.setupMatrices()
			if not error:
				print("Starting calculation")
				# Calculating the electronic transport
				error = self.calc.calcTransport()
		else:
			print("Starting calculation")
			#Refinement reuses already defined Calculation object with new input data
			try:
				error = self.calc.calcTransport(self.data.split())
			except:
				print('Transport refinement failed.')
		
		if (not error):
			print("Transport thread %s complete" %self.name)
			msg = formatTransport(self.calc.energy, self.calc.transmission)
			sendMessage(msg, TRANSPORT_COMPLETE)
			
			fig = Plot([self.calc.energy], [self.calc.transmission])
			filepath, err = SaveImage('Transport', fig, self.calc.energy, self.calc.transmission)
			if (not err):
				sendMessage(filepath.absolute().as_posix(), IMAGE_COMPLETE)
		else:
			print("something went wrong in transport calculations")
			sendMessage("something went wrong in transport calculations", ERROR)
			

def formatTransport(energy, transmission):
	msg = ""
	for idx, e in enumerate(energy):
		msg = msg + str(e) + " " + str(transmission[idx]) + " "
		
	msg = msg[0 : len(msg) - 1]
	return msg

def Plot(x,y,fig = None):
	# Expects a list of x- and y-lists.
	if fig == None:
		fig = plt.figure(figsize=(5, 5), dpi=216,facecolor='w', edgecolor='w')
	for en, tm in zip(x,y):
		plt.plot(en, tm, linewidth=1.5)
	#plt.legend(['CAP','Recursive'], loc = 'lower right')
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
	try:
		from pathlib import Path
		outpath = Path.home() / 'Pictures' / 'DeViNa'
		if not outpath.is_dir():
			print('Making new directory..')
			outpath.mkdir()

		imgCount = len([path for path in list(outpath.glob('*')) if path.suffix == '.png']) + 1
	except:
		return None, True
	
	try:
		string = name + "_" + str(imgCount)
		if type(x) == np.ndarray and type(y) == np.ndarray:
			np.savez(outpath / string, x=x, y=y)
		fig.savefig(outpath / string)
		plt.close(fig)
		string = string + ".png"
		return (outpath / string), False
	except Exception as e:
		print('Filename is invalid with message:\n{}'.format(e))
		return None, True
		
	
##############################################################################################
	
class imageThread(threading.Thread):
	def __init__(self, threadID, name, data):
		threading.Thread.__init__(self)
		self.threadID = threadID
		self.name = name
		self.data = data
		#self.openImages = []
		
	def run(self):
		# Something to open images
		#outPath = Path.cwd().parent / 'Images'
		try:
			outPath = Path.home() / 'Pictures' / 'DeViNa' #Default picture folder on windows, not sure about linux
			# sendMessage(outPath.absolute().as_posix(), IMAGE_COMPLETE) # Send only folder
			
			#imgPaths = [path for path in list(outpath.glob('*')) if path.suffix == '.png' and path not in self.openImages]
			imgPaths = [path for path in list(outPath.glob('*')) if path.suffix == '.png']
			for path in imgPaths:
				sendMessage(path.absolute().as_posix(), IMAGE_COMPLETE) # Send all image directories
				#self.openImages.append(path)
		except:
			sendMessage("Error sending file path", ERROR)


#open connection
openedConnection = False
while (not openedConnection):
	try:
		print("opening connection")
		openConnection(TCP_PORT, TCP_IP)
		openedConnection = True
		break
	except:
		print(sys.exc_info()[0])


threadID = 3
#initialize relaxation thread but don't run
relaxThread = relaxationThread(0, "relaxationThread", "")	
transThread = transportThread(1, "transportThread", "")
imageThread = imageThread(2, "imageThread", "")
	
print("Running.")

#receive message and deal with its content
try:
	while True:
		#receive incoming message
		data, msgType = receiveMessage()
		
		#echo Message back
		sendMessage(data, ECHO)	
		
		#deal with data
		#if echo msg was received and previous message was of COMPLETE type, check if conform to data sent
		if (msgType == ECHO):
			if (lastSentMsgType == COMPLETE):
				checkConformity(data, lastSentMessage)
		#check if thread is still running
		elif (msgType == RUN):
			if (relaxThread.isAlive()):
				sendMessage("Running", RUN)
			else:
				sendMessage("Not running", ERROR)	
		#create relaxation thread and start running it
		elif (msgType == RELAX):
			if (not relaxThread.isAlive()):
				relaxThread = relaxationThread(threadID, "relaxationThread", data)
				relaxThread.start()
				threadID += 1
			else:
				sendMessage("Relaxation is already running", ERROR)
        #create transport thread and start running it
		elif (msgType == TRANSPORT):
			refine = 0
			if (transThread.isAlive()):
				transThread.stop()				
			transThread.data = data
			transThread.run()
		#create refine transport
		elif (msgType == REFINE_TRANSPORT):
			refine = 1
			if (transThread.isAlive()):
				transThread.stop()
			transThread.data = data
			transThread.run()
		elif (msgType == IMAGE):
			if (imageThread.isAlive()):
				imageThread.stop()				
			imageThread.data = data
			imageThread.run()		
		#kill thread if running
		elif (msgType == STOP):
			if (relaxThread.isAlive()):
				relaxThread.exit()
				sendMessage("Computation stopped", ECHO)
			else:
				sendMessage("Was already stopped", ECHO)
		
		#Otherwise, just echo that msg type doesn't exist
		else:
			sendMessage("Unknown msg Type", ERROR)			
		
except KeyboardInterrupt:
	if (relaxThread.isAlive()):
		relaxThread.exit()
	try:
		sys.exit(0)
	except SystemExit:
		os._exit(0)

	
closeConnection()	