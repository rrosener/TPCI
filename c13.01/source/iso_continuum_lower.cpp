/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_continuum_lower - limit max prin. quan. no. due to continuum lowering processes	*/
#include "cddefines.h"
#include "phycon.h"
#include "dense.h"
#include "conv.h"
#include "iso.h"
#include "hydrogenic.h"
#include "trace.h"

void iso_continuum_lower( long ipISO, long nelem )
{
	double a;
	long int np, nd, ns, nc;
	long eff_charge;
	
	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* size of rate matrices will be defined according to the n calculated here	*/

	ASSERT( dense.xNucleiTotal < MAX_DENSITY );
	ASSERT( nelem < LIMELM );
	/* this may change at a future date. */
	ASSERT( ipISO <= 1 );

	eff_charge = nelem + 1 - ipISO;

	/* Particle packing - the density here should be density of all nuclei in the plasma */
	/* This one is just nuclear charge, which is independent of iso, and always nelem+1. */
	a = sqrt( 1.8887E8 * (nelem+1.) / pow((double)dense.xNucleiTotal, 0.333) );
	ASSERT( a > 0. );
	np = (long)( MIN2((double)INT16_MAX,floor(a)) );

	/* Debye shielding - the density here is electron density	*/
	/* This one depends on effective charge. */
	a = 2.6E7 * eff_charge * eff_charge * pow( phycon.te/dense.eden, 0.25);
	ASSERT( a > 0. );
	nd = (long)( MIN2((double)INT16_MAX,floor(a)) );

	/* Stark broadening - this should be the density of singly charged ions, 
	 * both positive and negative.  The sum of protons, electrons, and HeII should be
	 * good enough.	*/
	/* This one depends on effective charge. */
	a = 3171. * pow( (double)eff_charge, 0.8 ) * pow( dense.eden + (double)dense.xIonDense[ipHYDROGEN][1]
		+ (double)dense.xIonDense[ipHELIUM][1], -0.1333);
	ASSERT( a > 0. );
	ns = (long)( MIN2((double)INT16_MAX,floor(a)) );

	nc = MIN3(np, nd, ns);
	/* Don't allow continuum to be lowered below n=3. */
	nc = MAX2( nc, 3 );

	if( nc <= sp->n_HighestResolved_max)
	{
		sp->lgLevelsLowered = true;
		sp->lgLevelsEverLowered = true;
		sp->lgMustReeval = true;
		sp->n_HighestResolved_local = nc;
		sp->nCollapsed_local = 0;
		sp->numLevels_local = iso_get_total_num_levels( ipISO, nc, 0 );
	}
	/* Here is the case where the critical n lies among the one or more collapsed levels */
	/* we just get rid of any that are too high. */
	else if( nc <= sp->n_HighestResolved_max + sp->nCollapsed_max )
	{
		sp->lgLevelsLowered = true;
		sp->lgLevelsEverLowered = true;
		sp->lgMustReeval = true;
		sp->n_HighestResolved_local = sp->n_HighestResolved_max;
		sp->nCollapsed_local = nc - sp->n_HighestResolved_local;
		sp->numLevels_local = 
			iso_get_total_num_levels( ipISO, sp->n_HighestResolved_max, sp->nCollapsed_local );
	}
	/* This is usually where control will flow, because in most conditions the continuum will not be lowered.
	* Nothing changes in this case. */
	else
	{
		sp->numLevels_local = sp->numLevels_max;
		sp->nCollapsed_local = sp->nCollapsed_max;
		sp->n_HighestResolved_local = sp->n_HighestResolved_max;
		
		/* if levels were lowered on last pass but are not now, must reeval */
		if( sp->lgLevelsLowered )
		{
			sp->lgMustReeval = true;
		}
		else
		{
			sp->lgMustReeval = false;
		}

		sp->lgLevelsLowered = false;
	}
	
	if( !conv.nTotalIoniz )
		sp->lgMustReeval = true;

	/* None of these can be greater than that which was originally malloc'd. */
	ASSERT( sp->numLevels_local <= sp->numLevels_max );
	ASSERT( sp->nCollapsed_local <= sp->nCollapsed_max );
	ASSERT( sp->n_HighestResolved_local <= sp->n_HighestResolved_max );

	/* Lyman lines can not be greater than original malloc or critical pqn. */
	iso_ctrl.nLyman[ipISO] = MIN2( nc, iso_ctrl.nLyman_malloc[ipISO]);

	// zero out cooling and heating terms involving unused levels
	for( long ipHi=sp->numLevels_local; ipHi < sp->numLevels_max; ++ipHi )
	{
		for( long ipLo=0; ipLo < ipHi; ++ipLo )
			CollisionZero( sp->trans(ipHi,ipLo).Coll() );
	}

	if( trace.lgTrace && (trace.lgHBug||trace.lgHeBug)  )
	{
		fprintf( ioQQQ,"     iso_continuum_lower: ipISO %li nelem %li nc %li (np:%li,nd:%li,ns:%li) numLevels %li nCollapsed %li n_HighestResolved %li \n",
			ipISO, 
			nelem,
			nc,
			np, nd, ns, 
			sp->numLevels_local,
			sp->nCollapsed_local,
			sp->n_HighestResolved_local
			);
	}

	return;
}
