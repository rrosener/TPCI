/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DEUTERIUM_H_
#define DEUTERIUM_H_

class t_deuterium
{
public:
	t_deuterium()
	{
		lgElmtOn = false;
		gas_phase = 0.f;
		xIonDense[0] = 0.;
		xIonDense[1] = 0.;
		xMolecules = 0.f;
		fractionation = 0.f;
	}
	bool lgElmtOn;
	realnum gas_phase;
	double xIonDense[2];
	realnum xMolecules;
	realnum fractionation;
};

extern t_deuterium deut;

void ScaleDensitiesDeuterium( const realnum &factor );
void SetDeuteriumFractionation( const realnum &frac );
void SetGasPhaseDeuterium( const realnum &Hdensity );
void SetDeuteriumIonization( const double &xNeutral, const double &xIonized );
void InitDeuteriumIonization();

#endif /* DEUTERIUM_H_ */
