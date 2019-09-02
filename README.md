# VRMolecularBuilder
This project is a VR software for teaching, research and outreach in molecular physics and 2D materials.

At the atomic scale, correlation between form and function is key. Electron transport is no exception to the rule. 
Any subtle change in lattice structure, defect or external contamination can change radically the transmission of a given low-dimensional system. 
However, making accurate theoretical predictions of the transport properties of a single configuration can be very computation heavy, 
which prohibits rapid testing. This ultimately makes the development of an intuition for the phenomena involved at the atomic-scale an extremely slow proces.

The goal of this project is to leverage the interactivity and intuitivity of VR, to allow the user to experiment and obtain a grasp of how physical properties of 2D materials and nano-scaled structures 
evolve with geometry and various interaction. Ultimately, we would like to offer a comprehensive package of tools to learn and research electronic transport, photonics, and molecular mechanics. 

At the current stage we have implemented relaxation and a simplified version of electronic transport.

The project was developped using the Unreal Engine 4 game engine. The engine handles the grpahics and VR implementation, as well atoms binding.
All the physics simulation (relaxation and transport) are implemented in Python 3. It runs seperately on a TCP/IP client (also coded in python), the server being the executable made in UE4. 
This allows for a high modularity, as any researcher can easily plug in there python code this way.

#Requirements
Hardware:
- a HTC Vive (this is currently the only supported VR equipment.)
- a computer capable of running VR, typically a stationnary computer with i7, and a NVDIA GTX 1060 graphics card.
- A CUDA enable GPU if you want to run transport in GPU

Software:
- Windows 10
- UE4.20.3 (Using another version of the engine might not work as some functions may get deprecated).
- python3. We recommend the Anaconda distribution as it comes with a handful of usefull libraries, in use in this project.
- The following Add-ons for UE4.20.3 : 
		- Object Deliverer (for TCP/IP communication)
- The following libraries for python:
		- numpy
		- ASE (https://wiki.fysik.dtu.dk/ase/index.html)
		- matplotlib
		- plotly
		- cupy (if you want to run transport on GPU. Might not work on all platforms)
- GULP for windows if you want to be able to use relaxation. This can be downloaded for free for academics at https://nanochemistry.curtin.edu.au/gulp/request.cfm?rel=download
Once downloaded, add the folder to ".\VRmolecularBuilderV2\Content\ExternalScripts\" and edit the "environ" commands in "ASEscriptOnServer.py" to match your version of Gulp.
-Git LFS : A lot of the files in intermediate and binaries are fairly large. We use Git LFS to work around the issue.


#Usage
- To simply use as is, you can go download the latest build here: 
- If you simply want to contribute to the physics part, you can also download the latest build, using the link above. Then you can edit the python files in ".\VRmolecularBuilderV2\Content\ExternalScripts\"
- To contribute to the whole experience, simply clone the repository. (WARNING: Some things might appear very slow or laggy in the UE4 editor window, such as live refresh of graphs).

We are working on a more detailed documentation. Bear with us :)

If you have any question, want to contribute, or for commercial applications, please contact us at 
niels.pichon@outlook.com



		
