/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_hydro put H-iso seq into line intensity stack */
#include "cddefines.h"
#include "radius.h"
#include "thermal.h"
#include "dense.h"
#include "lines_service.h"
#include "grainvar.h"
#include "lines.h"

void lines_grains(void)
{
	double 
	  dhtot,
	  hold;
	long i;

	DEBUG_ENTRY( "lines_grains()" );

	if( !gv.lgGrainPhysicsOn )
	{
		return;
	}

	i = StuffComment( "grains" );
	linadd( 0., (realnum)i , "####", 'i',
		"the grain output");

	/* find total grain heating */
	dhtot = 0.;
	for( size_t nd=0; nd < gv.bin.size(); nd++ )
	{
		/* add heating due to all grain species that are included */
		dhtot += gv.bin[nd]->GasHeatPhotoEl;
	}

	/* total heating due to dust integrated over model */
	gv.TotalDustHeat += (realnum)(dhtot*radius.dVeffAper);
	/* largest fraction of local heating due to grains photo */
	gv.dphmax = MAX2((realnum)(dhtot/thermal.htot),gv.dphmax);
	/* largest local cooling of gas by collisions with grains */
	gv.dclmax = MAX2(gv.dclmax,(realnum)(gv.GasCoolColl/thermal.htot));

	/* largest relative number of electrons donated by grains */
	hold = SDIV(dense.EdenTrue);
	gv.GrnElecDonateMax = 
		(realnum)MAX2( gv.GrnElecDonateMax , gv.TotalEden/hold );

	/* largest relative number of electrons on grain surface */
	gv.GrnElecHoldMax = 
		(realnum)MAX2( gv.GrnElecHoldMax , -gv.TotalEden/hold );

	linadd(dhtot,0,"GrGH",'h',
		" gas heating by grain photoionization");

	linadd(thermal.heating[0][25],0,"GrTH",'h',
		" gas heating by thermionic emissions of grains ");

	linadd(MAX2(0.,gv.GasCoolColl),0,"GrGC",'c',
		"gas cooling by collisions with grains ");	

	linadd(MAX2(0.,-gv.GasCoolColl),0,"GrGC",'c',
		" gas heating by collisions with grains ");	

	linadd(gv.GrainHeatSum,0,"GraT",'i',
		" total grain heating by all sources, lines, collisions, incident continuum ");

	linadd(gv.GrainHeatInc,0,"GraI",'i',
		" grain heating by incident continuum ");

	linadd(gv.GrainHeatLya,1216,"GraL",'i',
		" grain heating due to destruction of Ly alpha  ");

	linadd(gv.GrainHeatCollSum,0,"GraC",'i',
		" grain heating due to collisions with gas ");

	linadd(gv.GrainHeatDif,0,"GraD",'i',
		" grain heating due to diffuse fields, may also have grain emission ");
	return;

}
