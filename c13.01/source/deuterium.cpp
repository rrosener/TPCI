/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "deuterium.h"
#include "dense.h"
#include "mole.h"

t_deuterium deut;

void ScaleDensitiesDeuterium( const realnum &factor )
{
	deut.gas_phase *= factor;
	deut.xMolecules *= factor;
	deut.xIonDense[0] *= factor;
	deut.xIonDense[1] *= factor;
	return;
}

void InitDeuteriumIonization()
{
	deut.xIonDense[0] = deut.gas_phase;
	deut.xIonDense[1] = 0.;
	return;
}

void SetDeuteriumIonization( const double &xNeutral, const double &xIonized )
{
	if( !deut.lgElmtOn )
		return;

	total_molecule_deut( deut.xMolecules );

	realnum total = deut.gas_phase - deut.xMolecules;

	fixit(); // try to let D ionization be independent of H
#if 1
	realnum neut = total * xNeutral/( xNeutral + xIonized );
	realnum ionz = total * xIonized/( xNeutral + xIonized );
#else
	realnum src = mole.findrk("D+,e-=>D,PHOTON") * dense.eden;
	realnum snk = mole.findrk("D,PHOTON=>D+,e-");
	realnum ion_ratio = xIonized/xNeutral;
	if( src > SMALLFLOAT )
		ion_ratio = snk/src;
	realnum neut = total / ( ion_ratio + 1.f );
	realnum ionz = total * ion_ratio / ( ion_ratio + 1.f );
#endif
	if( total > 1e-4 * deut.gas_phase )
	{
		deut.xIonDense[0] = neut;
		deut.xIonDense[1] = ionz;
	}

	//if( iteration >= 3 )
	//	fprintf( ioQQQ, "DEBUGGG SetDeuteriumIonization %li\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\t%.12e\n",
	//	nzone, deut.gas_phase, deut.xMolecules, total, old0, deut.xIonDense[0], old1, deut.xIonDense[1] );
	
	return;
}

void SetDeuteriumFractionation( const realnum &frac )
{
	if( !deut.lgElmtOn )
		return;
	deut.fractionation = frac;
	return;
}

void SetGasPhaseDeuterium( const realnum &Hdensity )
{
	if( !deut.lgElmtOn )
		return;
	deut.gas_phase = Hdensity * deut.fractionation;
}

