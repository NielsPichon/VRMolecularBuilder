// Fill out your copyright notice in the Description page of Project Settings.

#include "Public/importExportLib.h"
#include <fstream> //for file data stream
#include <iostream> //for input output data stream
//#include <Windows.h> //for directory check and creation //conflicts with ue4 libraries and prevents compilation
#include "Windows/MinWindows.h" //ue4 replacement for windows.h "because they didn't like the original lol"
#include <sstream>  //for std::istringstream
#include <iterator> //for std::istream_iterator
#include <vector>   //for std::vector 
#include "PlatformFilemanager.h" //for getting file list in folder
#include "LocalTimestampDirectoryVisitor.h"//for getting file list in folder
#include "Runtime/Core/Public/Misc/MonitoredProcess.h" //interactable proc
#include "FileManagerGeneric.h" //for checking last nmodified date


//ue4 file manager system
#include "paths.h"
#include "FileManager.h"


FMonitoredProcess relaxationInteractableProc(*FString(""), *FString(""), true);
FDateTime ASEstartTime;

bool UimportExportLib::exportXYZ(const TArray<FVector>& positions, const TArray<FString>& names, FString moleculeName)
{
	//Exports the position and name of the molecules'atoms into .xyz format. Note that no check is made on the prior 
	//existence of the export file. A warning to the user should be made prior to this.

	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();

	//Open 
	std::string strPath = std::string(TCHAR_TO_UTF8(*moleculeName)); //transform UE4 strings into real std strings
	std::ofstream fileStream(WorkDirectory + "/MoleculeLibrary/" + strPath + ".xyz"); // open data stream for writing the xyz file

	if (!fileStream)
	{
		return false; //If cannot open the specified file for some reason, return false and the exception will be handled in the blueprints
	}


	try
	{
		std::string strNbOfAtoms = std::to_string(positions.Num());

		fileStream << strNbOfAtoms << "\n";
		fileStream << "\\\\" << strPath << "\n";

		for (int i = 0; i < names.Num(); i++)
		{
			FString tmpName = names[i];
			std::string strName = std::string(TCHAR_TO_UTF8(*tmpName));
			std::string strX = std::to_string(positions[i].X);
			std::string strY = std::to_string(positions[i].Y);
			std::string strZ = std::to_string(positions[i].Z);
			fileStream << strName << " " << strX << " " << strY << " " << strZ << "\n";
		}
	}
	catch(...){
		fileStream.close();
		return false;
	}

	fileStream.close();

	return true;
}

bool UimportExportLib::importXYZ(TArray<FVector>& positions, TArray<FString>& names, FString moleculeName, TArray<FVector>& unitVectors, bool isUnitCell)
{
	//
	//takes a path to a xyz and reads it. The result is stored in the parsed referenced array inputs
	//

	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();

	std::string strPath = std::string(TCHAR_TO_UTF8(*moleculeName)); //transform UE4 strings into real std strings
	std::ifstream fileStream(WorkDirectory + "/MoleculeLibrary/" + strPath + ".xyz"); // open data stream for reading the xyz file
	if (!fileStream)
	{
		return false; //If cannot open the specified file for some reason, return false and the exception will be handled in the blueprints
	}


	try {
		//read first line to get total atom number
		std::string line;
		std::getline(fileStream, line);
		int atomNb = std::stoi(line);

		if (isUnitCell)
		{
			std::getline(fileStream, line);

			//split the line into tokens
			std::istringstream ss(line);
			std::istream_iterator<std::string> begin(ss), end;

			//putting all the tokens in the vector
			std::vector<std::string> splitLine(begin, end);


			//Extract coordiantes of all 3 unit vectors
			FVector unitVector;
			unitVector.X = std::stof(splitLine[1]);
			unitVector.Y = std::stof(splitLine[2]);
			unitVector.Z = std::stof(splitLine[3]);
			unitVectors.Add(unitVector);

			unitVector.X = std::stof(splitLine[5]);
			unitVector.Y = std::stof(splitLine[6]);
			unitVector.Z = std::stof(splitLine[7]);
			unitVectors.Add(unitVector);

			unitVector.X = std::stof(splitLine[9]);
			unitVector.Y = std::stof(splitLine[10]);
			unitVector.Z = std::stof(splitLine[11]);
			unitVectors.Add(unitVector);					   
		}
		else
		{
			//skip comment line which is basically useless comment
			std::getline(fileStream, line);
		}

		//for each atom, read the line, split it into usefull information and then store it
		for (int i = 0; i < atomNb; i++)
		{
			std::getline(fileStream, line);


			//split the line into tokens
			std::istringstream ss(line);
			std::istream_iterator<std::string> begin(ss), end;

			//putting all the tokens in the vector
			std::vector<std::string> splitLine(begin, end);

			//first element is atom symbol/name
			FString fstringName(splitLine[0].c_str());
			names.Add(fstringName);

			//the 3 next elements are the x y z components in angstrom
			FVector atomPosition;
			atomPosition.X = std::stof(splitLine[1]);
			atomPosition.Y = std::stof(splitLine[2]);
			atomPosition.Z = std::stof(splitLine[3]);

			////Logging instructions
			//UE_LOG(LogTemp, Warning, TEXT("%f"), atomPosition.X);
			//UE_LOG(LogTemp, Warning, TEXT("%f"), atomPosition.Y);
			//UE_LOG(LogTemp, Warning, TEXT("%f"), atomPosition.Z);

			positions.Add(atomPosition);
		}
	}
	catch (...) {
		return false;
	}


	return true;
}

bool UimportExportLib::moleculeFileExistsTest(FString moleculeName)
{

	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();

	//Create path from moleculeName
	std::string strPath = std::string(TCHAR_TO_UTF8(*moleculeName));
	strPath = WorkDirectory + "/MoleculeLibrary/" + strPath + ".xyz";

	//Check if the MoleculeLibrary directory exists first. If it doesn't create it and return false (no directory means no file ofc)
	std::string dirName = WorkDirectory + "/MoleculeLibrary/";
	DWORD ftyp = GetFileAttributesA(dirName.c_str());
	if (ftyp == INVALID_FILE_ATTRIBUTES)
	{
		CreateDirectoryA(dirName.c_str(), NULL);
		return false;
	}

	//test file prior existence
	return UimportExportLib::fileExistsTest(FString(strPath.c_str()));
}

bool UimportExportLib::fileExistsTest(FString path)
{
	struct stat buffer;
	return (stat(TCHAR_TO_UTF8(*path), &buffer) == 0);;
}

TArray<FString> UimportExportLib::GetFileList()
{
	const FString onlyFilesStartingWith = "";
	const FString onlyFilesWithExtension = "xyz";
	const bool fullPath = false;

	//Get molecule library directory
	std::string libPath = GetGameDirectoryPath() + "/MoleculeLibrary/";
	const FString directory = UTF8_TO_TCHAR(libPath.c_str());
	
	TArray<FString> directoriesToSkip;
	IPlatformFile &PlatformFile = FPlatformFileManager::Get().GetPlatformFile();
	FLocalTimestampDirectoryVisitor Visitor(PlatformFile, directoriesToSkip, directoriesToSkip, false);
	PlatformFile.IterateDirectory(*directory, Visitor);
	TArray<FString> files;

	for (TMap<FString, FDateTime>::TIterator TimestampIt(Visitor.FileTimes); TimestampIt; ++TimestampIt)
	{
		const FString filePath = TimestampIt.Key();
		const FString fileName = FPaths::GetCleanFilename(filePath);
		bool shouldAddFile = true;

		// Check if filename starts with required characters
		if (!onlyFilesStartingWith.IsEmpty())
		{
			const FString left = fileName.Left(onlyFilesStartingWith.Len());

			if (!(fileName.Left(onlyFilesStartingWith.Len()).Equals(onlyFilesStartingWith)))
				shouldAddFile = false;
		}

		// Check if file extension is required characters
		if (!onlyFilesWithExtension.IsEmpty())
		{
			if (!(FPaths::GetExtension(fileName, false).Equals(onlyFilesWithExtension, ESearchCase::IgnoreCase)))
				shouldAddFile = false;
		}

		// Add full path to results (without extension)
		if (shouldAddFile)
		{
			files.Add(fullPath ? filePath : FPaths::GetBaseFilename(filePath));
		}
	}
	
	//ignore tmp.xyz files which are data used for relaxation
	files.Remove(FString("tmp"));

	return files;
}

bool UimportExportLib::exportToNapuraXR(const TArray<FVector>& positions,
	const TArray<FString>& fullNameElements,
	const TArray<FVector>& bonds,
	const TArray<int>& owningMolecule,
	const TArray<int>& MM3Type,
	const float minTemperature,
	const float maxTemperature,
	const float eqTemperature,
	FString simName, 
	bool useOpenMMForceField)
{
	//Get working directory
	//std::string WorkDirectory = UimportExportLib::GetGameDirectoryPath();

	//Get path to simulation folder
	std::string* paths = UimportExportLib::GetPaths();
	if (!paths)
	{
		return false;
	}

	//Open 
	std::string strPath = std::string(TCHAR_TO_UTF8(*simName)); //transform UE4 strings into real std strings
	//std::ofstream fileStream(WorkDirectory + "/MoleculeLibrary/" + strPath + ".xml"); // open data stream for writing the xyz file
	std::ofstream fileStream(paths[0] + strPath + ".xml"); // open data stream for writing the xyz file

	if (!fileStream)
	{
		return false; //If cannot open the specified file for some reason, return false and the exception will be handled in the blueprints
	}

	//Ugly file writer. Could be improved with XML serialisation probably
	try {
		/* XML Header */
		fileStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << "\n";
		fileStream << "<Simulation Name = \"" << strPath << "\">" << "\n";

		/* System Properties */
		fileStream << "\t" << "<!-- Define the System Properties analogously to SystemProperties.cs-->" << "\n";
		fileStream << "\t" << "<SystemProperties>" << "\n";
		fileStream << "\t" << "\t" << "<SimulationBoundary SimulationBox=\"5, 5, 5\" MinimumBoxVolume=\"8\" />" << "\n";
		fileStream << "\t" << "\t" << "<Thermostat Type=\"BerendsenThermostat\" EquilibriumTemperature=\""
			<< eqTemperature << "\" MaximumTemperature=\"" << maxTemperature << "\" BerendsenCoupling=\"0.003\" />" << "\n";
		fileStream << "\t" << "\t" << "<Integrator Type=\"VelocityVerlet\" TimeStep=\"0.0005\" />" << "\n";
		fileStream << "\t" << "</SystemProperties>" << "\n";

		/* Topology */
		fileStream << "\t" << "<Topology>" << "\n";
		fileStream << "\t" << "\t" << "<Templates>" << "\n";

		int i = 0;
		int j = 0;
		int moleculeIdxStart = 0;
		int moleculeIdxEnd = 0;
		int moleculeNumber = owningMolecule[0];

		int templateCount = 1;

		fileStream << "\t" << "\t" << "\t" << "<Residue Name=\"" << moleculeNumber << "\">" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "<Atoms>" << "\n";

		while (i < fullNameElements.Num())
		{
			if (moleculeNumber != owningMolecule[i])
			{
				moleculeIdxEnd = i - 1;
				fileStream << "\t" << "\t" << "\t" << "\t" << "</Atoms>" << "\n";
				//Add bonds
				fileStream << "\t" << "\t" << "\t" << "\t" << "<Bonds>" << "\n";
				while (j < bonds.Num() && bonds[j].Z == moleculeNumber)
				{
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Bond A = \""
						<< std::to_string(int(bonds[j].X)) << "\" B = \"" << std::to_string(int(bonds[j].Y)) << "\" />" << "\n";
					j++;
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "</Bonds>" << "\n";

				//Add force field
				fileStream << "\t" << "\t" << "\t" << "\t" << "<ForceFields>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<InteractiveGaussianForceField GradientScaleFactor=\"500\" />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3ForceField>" << "\n";
				//fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<DataFile File=\"^ / Assets / Data / mm3.xml\" />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMappings>" << "\n";
				for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
				{
					std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMapping AtomPath = \"" << j << "\" Type = \"" + typeNb + "\" />" << "\n";
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3AtomMappings>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<AdditionalTerms />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3ForceField>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesForceField CutOffDistance=\"0.5\">" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMappings>" << "\n";
				for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
				{
					std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMapping AtomPath = \"" << j << "\" Sigma=\"0\" Epsilon=\"0\" MM3Type=\"" + typeNb + "\" />" << "\n";
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesAtomMappings>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesForceField>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "</ForceFields>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "</Residue>" << "\n";

				//new residue
				moleculeNumber = owningMolecule[i];
				fileStream << "\t" << "\t" << "\t" << "<Residue Name=\"" << moleculeNumber << "\">" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "<Atoms>" << "\n";

				moleculeIdxStart = i;
				templateCount++;
			}

			std::string element = std::string(TCHAR_TO_UTF8(*(fullNameElements[i])));

			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Atom Element = \""
				<< element << "\" Position = \""
				<< std::to_string(positions[i].X) << ", " << std::to_string(positions[i].Y) << ", " << std::to_string(positions[i].Z) << "\" />" << "\n";

			i++;
		}


		moleculeIdxEnd = i - 1;
		fileStream << "\t" << "\t" << "\t" << "\t" << "</Atoms>" << "\n";

		//add bonds
		fileStream << "\t" << "\t" << "\t" << "\t" << "<Bonds>" << "\n";
		while (j < bonds.Num() && bonds[j].Z == moleculeNumber)
		{
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Bond A = \""
				<< std::to_string(int(bonds[j].X)) << "\" B = \"" << std::to_string(int(bonds[j].Y)) << "\" />" << "\n";
			j++;
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "</Bonds>" << "\n";

		//add force field
		fileStream << "\t" << "\t" << "\t" << "\t" << "<ForceFields>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<InteractiveGaussianForceField GradientScaleFactor=\"500\" />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3ForceField>" << "\n";
		//fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<DataFile File=\"^ / Assets / Data / mm3.xml\" />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMappings>" << "\n";
		for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
		{
			std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMapping AtomPath = \"" << j << "\" Type = \"" + typeNb + "\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3AtomMappings>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<AdditionalTerms />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3ForceField>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesForceField CutOffDistance=\"0.5\">" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMappings>" << "\n";
		for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
		{
			std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMapping AtomPath = \"" << j << "\" Sigma=\"0\" Epsilon=\"0\" MM3Type=\"" + typeNb + "\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesAtomMappings>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesForceField>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "</ForceFields>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "</Residue>" << "\n";

		fileStream << "\t" << "\t" << "</Templates>" << "\n";

		/* Spawners */
		fileStream << "\t" << "\t" << "<Spawners>" << "\n";
		for (int j = 0; j < templateCount; j++)
		{
			fileStream << "\t" << "\t" << "\t" << "<Spawner Name=\"" << j << "\" Template=\"" << j << "\" Count=\"1\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "</Spawners>" << "\n";

		fileStream << "\t" << "</Topology>" << "\n";
		fileStream << "</Simulation>" << "\n";

	}
	catch (...)
	{
		fileStream.close();
		return false;
	}

	fileStream.close();

	////copy the file to the simulation directory from Napura
	//std::string* paths = UimportExportLib::GetPaths();

	//FString srcStr((WorkDirectory + "/MoleculeLibrary/" + strPath + ".xml").c_str());
	//const TCHAR * src = *srcStr;

	//FString tgtStr((paths[0] + strPath + ".xml").c_str());
	//const TCHAR * tgt = *srcStr;

	return true;
}

void UimportExportLib::exportToPDB(FString moleculeName, bool makeForceField, bool runSimulation, int stepNb)
{
	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();

	//Open 
	std::string fileName = std::string(TCHAR_TO_UTF8(*moleculeName)); //transform UE4 strings into real std strings
	std::string pdbPath = WorkDirectory + "/MoleculeLibrary/" + fileName + ".pdb"; // path to the corresponding pdb file to create

	std::string param = fileName + ".xyz" + " -O " + fileName + ".pdb";
	FPlatformProcess::CreateProc(*FString("obabel"), *FString(param.c_str()), true, false, false, 0, 0, *FString((WorkDirectory + "/MoleculeLibrary/").c_str()), nullptr);

	if (makeForceField)
		makeOpenMMForceField(pdbPath, runSimulation, stepNb);

	return;
}

void UimportExportLib::ASErelaxation(TArray<FString> atomNames, 
	TArray<FVector> atomPositions,
	bool startFromScratch,
	int stepNb)
{
	std::string WorkDirectory = GetGameDirectoryPath();
	
	//format data to parse to ASE
	FString Finput = UimportExportLib::formatDataForASE(atomNames, atomPositions, startFromScratch, stepNb);
	std::string input = std::string(TCHAR_TO_UTF8(*Finput));


	//debug
	std::ofstream fileStream(WorkDirectory + "relaxationDebug.txt"); // open data stream for writing the xyz file
	fileStream << input;
	fileStream.close();
	//

	//UE_LOG(LogTemp, Warning, TEXT("Running ASE Script"));

	//register start time to compare to tmp.xyz last modification time
	ASEstartTime = FDateTime::UtcNow();

	FString timeToString = ASEstartTime.ToString();
	//UE_LOG(LogTemp, Warning, TEXT("Proc creation time %s"), *timeToString);
	FPlatformProcess::CreateProc(*FString("python"), *FString(("ASEscript.py " + input).c_str()), true, false, false, 0, 0, *FString((WorkDirectory + "/ExternalScripts/").c_str()), nullptr);
	/*relaxationInteractableProc = FMonitoredProcess(*FString("\"C:/ProgramData/Anaconda3/python.exe\""), *FString(input.c_str()), false, true);
	if (!relaxationInteractableProc.Init())
	{
		UE_LOG(LogTemp, Warning, TEXT("ASE initialisation failed"));
	}
	else
	{
		relaxationInteractableProc.Run();

		if (relaxationInteractableProc.IsRunning())
		{
			UE_LOG(LogTemp, Warning, TEXT("ASE succesfully launched"));
		}
		else
		{
			int returnCode = relaxationInteractableProc.GetReturnCode();
			FString code = FString((std::to_string(returnCode)).c_str());
			UE_LOG(LogTemp, Warning, TEXT("Something went wrong in ASE. Return code : %s"), *code);
		}
	}*/
}

FString UimportExportLib::formatDataForASE(TArray<FString> atomNames, TArray<FVector> atomPositions, bool startFromScratch, int stepNb)
{

	std::string input = "";
	input += std::to_string(stepNb);
	if (startFromScratch)
		input += " 1 ";
	else
		input += " 0 ";

	for (int i = 0; i < atomNames.Num(); i++)
	{
		input += std::string(TCHAR_TO_UTF8(*(atomNames[i]))) + " ";
		input += std::to_string(atomPositions[i].X) + " ";
		input += std::to_string(atomPositions[i].Y) + " ";
		input += std::to_string(atomPositions[i].Z) + " ";
	}

	return FString(input.c_str());
}

bool UimportExportLib::checkRelaxationFinished()
{

	if (UimportExportLib::fileExistsTest("tmp"))
	{
		//Get working directory
		std::string WorkDirectory = GetGameDirectoryPath();

		FString fileName = FString(WorkDirectory.c_str()) + "/MoleculeLibrary/tmp.xyz"; 

		FFileManagerGeneric fm;
		FDateTime mod = fm.GetTimeStamp(*fileName);
	/*	FString timeToString = mod.ToString();
		UE_LOG(LogTemp, Warning, TEXT("tmp.xyz last modification %s"), *timeToString);*/

		if (mod > ASEstartTime)
		{
			return true;
		}
		else
		{
			//UE_LOG(LogTemp, Warning, TEXT("tmp hasn't yet be overwritten"));
			return false;
		}
	}
	else 
	{
		//UE_LOG(LogTemp, Warning, TEXT("tmp File does not exist"));
		return false;
	}

	//return !(relaxationInteractableProc.IsRunning());
}

void UimportExportLib::getRelaxationData(TArray<FVector>& atomPositions)
{
	//relaxationInteractableProc.Exit();

	//get data
	TArray<FString> names;
	TArray<FVector> unitVectors;
	UimportExportLib::importXYZ(atomPositions, names, FString("tmp"), unitVectors);

	//Center barycenter of atoms
	FVector center = FVector(0, 0, 0);
	for (int i = 0; i < atomPositions.Num(); i++)
	{
		center += atomPositions[i];
	}

	center /= (float)(atomPositions.Num());

	for (int i = 0; i < atomPositions.Num(); i++)
	{
		atomPositions[i] -= center;
	}

	return;
}

void UimportExportLib::forceKillRelaxation()
{
	system("taskkill /IM python.exe /F");
	return;
}

bool UimportExportLib::startRemoteClient()
{
	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();


	try {
		//start server
		FPlatformProcess::CreateProc(*FString("python"), *FString("TCPPhysics.py"), true, false, false, 0, 0, *FString((WorkDirectory + "/ExternalScripts/").c_str()), nullptr);
	}
	catch (...) {
		return false;
	}

	return true;
}

void UimportExportLib::makeOpenMMForceField(std::string fileName, bool runSimulation, int stepNb)
{
	if (!runSimulation)
		stepNb = 0;

	UE_LOG(LogTemp, Warning, TEXT("Running Python Script"));

	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();
	std::string inputPath = "\"" + WorkDirectory + "/MoleculeLibrary/" + fileName + ".pdb\"";
	std::string outputPath = "\"" + WorkDirectory + "/MoleculeLibrary/" + fileName + "_output.pdb\"";
	std::string params = WorkDirectory + "/ExternalScripts/NewOpenMMsim.py " + std::to_string(stepNb) + " " + inputPath + " " + outputPath;

	//FPlatformProcess::CreateProc(*FString("C:\ProgramData\Anaconda3\python.exe"), *FString(params.c_str()), true, false, false, 0, 0, nullptr, nullptr);
}

std::string UimportExportLib::GetGameDirectoryPath()
{
	FString RelativePath = FPaths::GameContentDir();

	FString FullPath = IFileManager::Get().ConvertToAbsolutePathForExternalAppForRead(*RelativePath);
	//Log
	//UE_LOG(LogTemp, Warning, TEXT("%s"), *FullPath);

	std::string strPath(TCHAR_TO_UTF8(*FullPath));

	return strPath;
}

std::string * UimportExportLib::GetPaths()
{
	//Get working directory
	std::string WorkDirectory = GetGameDirectoryPath();

	std::string* paths = new std::string[2];


	std::ifstream fileStream(WorkDirectory + "/Paths.txt"); // open data stream for reading the xyz file
	if (!fileStream)
	{
		return false; //If cannot open the specified file for some reason, return false and the exception will be handled in the blueprints
	}


	try {
		// read first line to get total atom number
		std::string line;

		//skip one line
		std::getline(fileStream, line);

		//store path to simulation folder
		std::getline(fileStream, line);
		paths[0] = line;

		//skip 2 lines
		std::getline(fileStream, line);
		std::getline(fileStream, line);

		//store path to NapuraXR
		std::getline(fileStream, line);
		paths[1] = line;

		fileStream.close();
	}
	catch (...) {
		return new std::string[2];
	}

	return paths;
}

//I don't really know why but it seems that this method is linked to something in the engine somehow and I can't delete it even if it is not used.
bool UimportExportLib::exportToNapuraSimulationXML(const TArray<FVector>& positions,
	const TArray<FString>& fullNameElements,
	const TArray<FVector>& bonds,
	const TArray<int>& owningMolecule,
	const TArray<int>& MM3Type,
	const float minTemperature,
	const float maxTemperature,
	const float eqTemperature,
	FString simName)
{
	//Get working directory
	std::string WorkDirectory = UimportExportLib::GetGameDirectoryPath();

	//Open 
	std::string strPath = std::string(TCHAR_TO_UTF8(*simName)); //transform UE4 strings into real std strings
	std::ofstream fileStream(WorkDirectory + "/MoleculeLibrary/" + strPath + ".xml"); // open data stream for writing the xyz file

	if (!fileStream)
	{
		return false; //If cannot open the specified file for some reason, return false and the exception will be handled in the blueprints
	}


	try {
		//Ugly file writer. Could be improved with XML serialisation probably

		/* XML Header */
		fileStream << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << "\n";
		fileStream << "<Simulation Name = \"" << strPath << "\">" << "\n";

		/* System Properties */
		fileStream << "\t" << "<!-- Define the System Properties analogously to SystemProperties.cs-->" << "\n";
		fileStream << "\t" << "<SystemProperties>" << "\n";
		fileStream << "\t" << "\t" << "<SimulationBoundary SimulationBox=\"2, 2, 2\" MinimumBoxVolume=\"8\" />" << "\n";
		fileStream << "\t" << "\t" << "<Thermostat Type=\"BerendsenThermostat\" EquilibriumTemperature=\""
			<< eqTemperature << "\" MaximumTemperature=\"" << maxTemperature << "\" BerendsenCoupling=\"0.003\" />" << "\n";
		fileStream << "\t" << "\t" << "<Integrator Type=\"VelocityVerlet\" TimeStep=\"0.0005\" />" << "\n";
		fileStream << "\t" << "</SystemProperties>" << "\n";

		/* Topology */
		fileStream << "\t" << "<Topology>" << "\n";
		fileStream << "\t" << "\t" << "<Templates>" << "\n";

		int i = 0;
		int j = 0;
		int moleculeIdxStart = 0;
		int moleculeIdxEnd = 0;
		int moleculeNumber = owningMolecule[0];

		int templateCount = 1;

		fileStream << "\t" << "\t" << "\t" << "<Residue Name=\"" << moleculeNumber << "\">" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "<Atoms>" << "\n";

		while (i < fullNameElements.Num())
		{
			if (moleculeNumber != owningMolecule[i])
			{
				moleculeIdxEnd = i - 1;
				fileStream << "\t" << "\t" << "\t" << "\t" << "</Atoms>" << "\n";
				//Add bonds
				fileStream << "\t" << "\t" << "\t" << "\t" << "<Bonds>" << "\n";
				while (j < bonds.Num() && bonds[j].Z == moleculeNumber)
				{
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Bond A = \""
						<< std::to_string(int(bonds[j].X)) << "\" B = \"" << std::to_string(int(bonds[j].Y)) << "\" />" << "\n";
					j++;
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "</Bonds>" << "\n";

				//Add force field
				fileStream << "\t" << "\t" << "\t" << "\t" << "<ForceFields>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<InteractiveGaussianForceField GradientScaleFactor=\"500\" />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3ForceField>" << "\n";
				//fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<DataFile File=\"^ / Assets / Data / mm3.xml\" />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMappings>" << "\n";
				for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
				{
					std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMapping AtomPath = \"" << j << "\" Type=\"" + typeNb + "\" />" << "\n";
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3AtomMappings>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<AdditionalTerms />" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3ForceField>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesForceField CutOffDistance=\"0.5\">" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMappings>" << "\n";
				for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
				{
					std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
					fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMapping AtomPath = \"" << j << "\" Sigma=\"0\" Epsilon=\"0\" MM3Type=\"" + typeNb + "\" />" << "\n";
				}
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesAtomMappings>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesForceField>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "</ForceFields>" << "\n";
				fileStream << "\t" << "\t" << "\t" << "</Residue>" << "\n";

				//new residue
				moleculeNumber = owningMolecule[i];
				fileStream << "\t" << "\t" << "\t" << "<Residue Name=\"" << moleculeNumber << "\">" << "\n";
				fileStream << "\t" << "\t" << "\t" << "\t" << "<Atoms>" << "\n";

				moleculeIdxStart = i;
				templateCount++;
			}

			std::string element = std::string(TCHAR_TO_UTF8(*(fullNameElements[i])));

			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Atom Element = \""
				<< element << "\" Position = \""
				<< std::to_string(positions[i].X) << ", " << std::to_string(positions[i].Y) << ", " << std::to_string(positions[i].Z) << "\" />" << "\n";

			i++;
		}


		moleculeIdxEnd = i - 1;
		fileStream << "\t" << "\t" << "\t" << "\t" << "</Atoms>" << "\n";

		//add bonds
		fileStream << "\t" << "\t" << "\t" << "\t" << "<Bonds>" << "\n";
		while (j < bonds.Num() && bonds[j].Z == moleculeNumber)
		{
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<Bond A = \""
				<< std::to_string(int(bonds[j].X)) << "\" B = \"" << std::to_string(int(bonds[j].Y)) << "\" />" << "\n";
			j++;
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "</Bonds>" << "\n";

		//add force field
		fileStream << "\t" << "\t" << "\t" << "\t" << "<ForceFields>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<InteractiveGaussianForceField GradientScaleFactor=\"500\" />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3ForceField>" << "\n";
		//fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<DataFile File=\"^ / Assets / Data / mm3.xml\" />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMappings>" << "\n";
		for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
		{
			std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<MM3AtomMapping AtomPath = \"" << j << "\" Type=\"" + typeNb + "\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3AtomMappings>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<AdditionalTerms />" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</MM3ForceField>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesForceField CutOffDistance=\"0.5\">" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMappings>" << "\n";
		for (int j = 0; j <= moleculeIdxEnd - moleculeIdxStart; j++)
		{
			std::string typeNb = std::to_string(MM3Type[j + moleculeIdxStart]);
			fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "<LennardJonesAtomMapping AtomPath = \"" << j << "\" Sigma=\"0\" Epsilon=\"0\" MM3Type=\"" + typeNb + "\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesAtomMappings>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "\t" << "</LennardJonesForceField>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "\t" << "</ForceFields>" << "\n";
		fileStream << "\t" << "\t" << "\t" << "</Residue>" << "\n";

		fileStream << "\t" << "\t" << "</Templates>" << "\n";

		/* Spawners */
		fileStream << "\t" << "\t" << "<Spawners>" << "\n";
		for (int j = 0; j < templateCount; j++)
		{
			fileStream << "\t" << "\t" << "\t" << "<Spawner Name=\"" << j << "\" Template=\"" << j << "\" Count=\"1\" />" << "\n";
		}
		fileStream << "\t" << "\t" << "</Spawners>" << "\n";

		fileStream << "\t" << "</Topology>" << "\n";
		fileStream << "</Simulation>" << "\n";
	}
	catch (...) {
		fileStream.close();
		return false;
	}

	fileStream.close();
	return true;
}




