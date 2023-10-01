/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "collision.h"

#include "dense.h"
#include "h2.h"
#include "hmi.h"

ColliderList::ColliderList()
{
	DEBUG_ENTRY(" ColliderList::ColliderList()");

	list.resize( ipNCOLLIDER );
	list[ipELECTRON].charge = -1;
	list[ipELECTRON].mass_amu = ELECTRON_MASS/ATOMIC_MASS_UNIT;
	
	list[ipPROTON].charge = 1;
	list[ipPROTON].mass_amu = dense.AtomicWeight[0];
	
	list[ipHE_PLUS].charge = 1;
	list[ipHE_PLUS].mass_amu = dense.AtomicWeight[1];
	
	list[ipALPHA].charge = 2;
	list[ipALPHA].mass_amu = dense.AtomicWeight[1];
	
	list[ipATOM_H].charge = 0;
	list[ipATOM_H].mass_amu = dense.AtomicWeight[0];
	
	list[ipATOM_HE].charge = 0;
	list[ipATOM_HE].mass_amu = dense.AtomicWeight[1];
	
	list[ipH2_ORTHO].charge = 0;
	list[ipH2_ORTHO].mass_amu = 2.f;
	
	list[ipH2_PARA].charge = 0;
	list[ipH2_PARA].mass_amu = 2.f;
	
	list[ipH2].charge = 0;
	list[ipH2].mass_amu = 2.f;
}

void ColliderList::init()
{
	DEBUG_ENTRY(" ColliderList::init()");
	colliders.list[ipELECTRON].density = &(dense.EdenHCorr);
	colliders.list[ipPROTON].density = &(dense.xIonDense[ipHYDROGEN][1]);
	colliders.list[ipHE_PLUS].density = &(dense.xIonDense[ipHELIUM][1]);
	colliders.list[ipALPHA].density = &(dense.xIonDense[ipHELIUM][2]);
	colliders.list[ipATOM_H].density = &(dense.xIonDense[ipHYDROGEN][0]);
	colliders.list[ipATOM_HE].density = &(dense.xIonDense[ipHELIUM][0]);
	colliders.list[ipH2_ORTHO].density = &(h2.ortho_density);
	colliders.list[ipH2_PARA].density = &(h2.para_density);
	colliders.list[ipH2].density = &(hmi.H2_total);
}
ColliderList colliders;

/*CollisionJunk set all elements of transition struc to dangerous values */
void CollisionJunk( const CollisionProxy & t )
{

	DEBUG_ENTRY( "CollisionJunk()" );

	/* Coll->cooling and Coll->heating due to collisional excitation */
	t.cool() = -FLT_MAX;
	t.heat() = -FLT_MAX;

	/* collision strengths for transition */
	t.col_str() = -FLT_MAX;

	for( long i=0; i<ipNCOLLIDER; i++ )
		t.rate_coef_ul_set()[i] = 0.f;

	t.rate_lu_nontherm_set() = 0.f;

	return;
}

/*CollisionZero zeros out the structure */
void CollisionZero( const CollisionProxy & t )
{

	DEBUG_ENTRY( "CollisionZero()" );

	/* Coll->cooling and Coll->heating due to collisional excitation */
	t.cool() = 0.;
	t.heat() = 0.;

	return;
}

