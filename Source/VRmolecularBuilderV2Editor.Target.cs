// Fill out your copyright notice in the Description page of Project Settings.

using UnrealBuildTool;
using System.Collections.Generic;

public class VRmolecularBuilderV2EditorTarget : TargetRules
{
	public VRmolecularBuilderV2EditorTarget(TargetInfo Target) : base(Target)
	{
		Type = TargetType.Editor;

		ExtraModuleNames.AddRange( new string[] { "VRmolecularBuilderV2" } );
	}
}
