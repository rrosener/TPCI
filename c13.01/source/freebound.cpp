/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "freebound.h"

void freeBound::Reset()
{
	PopLTE = 0.;
	ColIoniz = 0.;
	gamnc = -DBL_MAX;
	RecomInducRate = -DBL_MAX;
	DepartCoef = -DBL_MAX;
	RateLevel2Cont = 0.;
	RateCont2Level = 0.;
	ConOpacRatio = 1.;
	RadRecCon = 0.;
	RadRecCoolCoef = 0.;
	RadRecomb[ipRecRad] = 0.;
	RadRecomb[ipRecNetEsc] = 1.;
	RadRecomb[ipRecEsc] = 1.;
	DielecRecomb = 0.;
	PhotoHeat = 0.;

	return;
};

