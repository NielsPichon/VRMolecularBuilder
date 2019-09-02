// Fill out your copyright notice in the Description page of Project Settings.

#include "ShellScript.h"
#include "Public/importExportLib.h"

bool UShellScript::goToNapura()
{
	std::string* paths = UimportExportLib::GetPaths();

	if (!UimportExportLib::fileExistsTest(FString(paths[1].c_str())))
		return false;

	try {
		system(paths[1].c_str());
	}
	catch (...) {
		return false;
	}

	return true;
}

