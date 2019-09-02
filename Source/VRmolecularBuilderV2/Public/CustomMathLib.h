// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "CustomMathLib.generated.h"

/**
 * 
 */
UCLASS()
class VRMOLECULARBUILDERV2_API UCustomMathLib : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
	
	UFUNCTION(BlueprintCallable)
	static void interpolate(TArray<float>& Xpoints,
		TArray<float>& Ypoints,
		const TArray<float>& Xdata,
		const TArray<float>& Ydata,
		bool smooth);
	
	
};
