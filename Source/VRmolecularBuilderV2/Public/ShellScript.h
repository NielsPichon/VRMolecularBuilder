// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Kismet/BlueprintFunctionLibrary.h"
#include "ShellScript.generated.h"

/**
 * 
 */
UCLASS()
class VRMOLECULARBUILDERV2_API UShellScript : public UBlueprintFunctionLibrary
{
	GENERATED_BODY()
	
	UFUNCTION(BlueprintCallable)
	static bool goToNapura();
};
