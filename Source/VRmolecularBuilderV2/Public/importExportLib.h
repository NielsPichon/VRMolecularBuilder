// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include <string>

//always include this last
#include "importExportLib.generated.h"

/**
 * 
 */
UCLASS()
class VRMOLECULARBUILDERV2_API UimportExportLib : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
	
		UFUNCTION(BlueprintCallable)
		static bool exportXYZ(const TArray<FVector>& positions, const TArray<FString>& names, FString path);

		UFUNCTION(BlueprintCallable)
		static bool importXYZ(TArray<FVector>& positions, TArray<FString>& names, FString moleculeName, TArray<FVector>& unitVectors, bool isUnitCell = false);

		UFUNCTION(BlueprintCallable)
		static bool exportToNapuraSimulationXML(const TArray<FVector>& positions,
			const TArray<FString>& fullNameElements,
			const TArray<FVector>& bonds,
			const TArray<int>& owningMolecule,
			const TArray<int>& MM3Type,
			const float minTemperature,
			const float maxTemperature,
			const float eqTemperature,
			FString moleculeName);

		UFUNCTION(BlueprintCallable)
		static bool moleculeFileExistsTest(FString moleculeName);


		UFUNCTION(BlueprintCallable)
		static TArray<FString> GetFileList();

		UFUNCTION(BlueprintCallable)
		static bool exportToNapuraXR(const TArray<FVector>& positions,
			const TArray<FString>& fullNameElements,
			const TArray<FVector>& bonds,
			const TArray<int>& owningMolecule,
			const TArray<int>& MM3Type,
			const float minTemperature,
			const float maxTemperature,
			const float eqTemperature,
			FString simName,
			bool useOpenMMForceField = false);

		UFUNCTION(BlueprintCallable)
		static void exportToPDB(FString MoleculeName, 
			bool makeForceField,
			bool runSimulation = false, 
			int stepNb = 0);

		UFUNCTION(BlueprintCallable)
		static void ASErelaxation(TArray<FString> atomNames,
			TArray<FVector> atomPositions, 
			bool startFromScratch = true, 
			int stepNb = 10000);

		UFUNCTION(BlueprintCallable)
		static FString formatDataForASE(TArray<FString> atomNames,
			TArray<FVector> atomPositions,
			bool startFromScratch = true,
			int stepNb = 10000);

		UFUNCTION(BlueprintCallable)
		static bool checkRelaxationFinished();

		UFUNCTION(BlueprintCallable)
		static void getRelaxationData(TArray<FVector>& atomPositions);

		UFUNCTION(BlueprintCallable)
		static void forceKillRelaxation();

		UFUNCTION(BlueprintCallable)
		static bool startRemoteClient();

		static void makeOpenMMForceField(std::string path, bool runSimulation = false , int stepNb = 0);

	static std::string GetGameDirectoryPath();

public:
	static std::string* GetPaths();
	static bool fileExistsTest(FString path);
};
