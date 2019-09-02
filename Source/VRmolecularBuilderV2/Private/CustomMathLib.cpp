// Fill out your copyright notice in the Description page of Project Settings.

#include "CustomMathLib.h"

void UCustomMathLib::interpolate(TArray<float>& Xpoints, TArray<float>& Ypoints, const TArray<float>& Xdata, const TArray<float>& Ydata, bool smooth)
{
	int idx = 0;

	for (int i = 0; i < Xpoints.Num(); i++) 
	{
		if (Xpoints[i] < Xdata[0])
			Ypoints[i] = Ydata[0];
		else
		{
			while (Xpoints[i] < Xdata[idx])
			{
				idx++;
				if (idx >= Xdata.Num() - 1)
					break;
			}
			if (idx >= Xdata.Num() - 1)
				Ypoints[i] = Ydata.Last();
			else
			{
				if (smooth) {
					float dxn = 0;
					float dxp = 0;
					if (i > 0)
						dxn = (Xpoints[i - 1] + Xpoints[i]) * 0.5f;
					if (i < Xpoints.Num() - 1)
						dxp = (Xpoints[i + 1] + Xpoints[i]) * 0.5f;

					float buffer = Ydata[idx];

					int k = 0;
					while (Xdata[idx + k] > Xpoints[i] - dxn )
					{

					}


				}
				else {
					Ypoints[i] = (Ydata[idx] + Ydata[idx + 1]) * 0.5f;
				}
			}
			
		}

	}
}
