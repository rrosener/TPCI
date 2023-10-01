/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*dense_tabden interpolate on table of points for density with dlaw table command, by K Volk */

/* CHANGES: (M. Salz 14.05.2013)
 *  - add the if clause for interpolation on the linear
 *    values, not the log10 values 
 *  - change table lookup from linear to bisection 
 *  - add if clause for nzone = 0 
 *  - changed interpolation on log depth to interpolation on
      on delta log: log(depth) - log(rinner)
      --> needed to import radius.h */
/* CHANGES: (M. Salz 17.05.2013)
 *  - superseeded by interpol_tabulated()
      --> NOT called anymore */

#include "cddefines.h"
#include "dense.h"
#include "radius.h"

double dense_tabden(double r0, 
  double depth)
{
	long int jlow, jmid, jhigh;
	double frac, 
	  tabden_v, 
	  x;

	DEBUG_ENTRY( "dense_tabden()" );
	/*interpolate on table of points for density with dlaw table command, by K Volk
	 *each line is log radius and H density per cc. */

	if( nzone == 0 ){ /* radius values are not set correctly yet */
		tabden_v = dense.tabval[0];
	}
	else
	{
		/*begin sanity check */
		if( r0 < 0. || depth <= 0. )
		{
			fprintf( ioQQQ, " dense_tabden called with insane depth, radius, =%10.2e%10.2e\n", 
			  depth, r0 );
		}
		/* end sanity check */
		
		/* interpolate linear or on logscale? */
		if ( dense.lgTabLinear )
		{
			/* interpolate on radius or depth? */
			if( dense.lgTabDepth )
			{
				/* depth key appeared = we want depth */
				x = depth;
			}
			else
			{
				/* use radius */
				x = r0;
			}
		}
		else
		{
			/* interpolate on radius or depth? */
			if( dense.lgTabDepth )
			{
				/* depth key appeared = we want depth */
				x = log10(depth/radius.rinner);
			}
			else
			{
				/* use radius */
				x = log10(r0);
			}
		}

		/* set to impossible value, will crash if not reset */
		tabden_v = -DBL_MAX;

		if( x < dense.tabrad[0] || x >= dense.tabrad[dense.nvals-1] )
		{
			fprintf( ioQQQ, " requested radius outside range of dense_tabden\n" );
			if ( dense.lgTabLinear )
			{
				fprintf( ioQQQ, " radius/depth was %10.2e min, max=%10.2e %10.2e\n",
				  x, dense.tabrad[0], dense.tabrad[dense.nvals-1] );
			}
			else
			{
				fprintf( ioQQQ, " log10(radius/depth) was %10.2f min, max=%10.2f %10.2f\n",
				  x, dense.tabrad[0], dense.tabrad[dense.nvals-1] );
			}
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			/* table lookup */
			jlow = 0;
			jhigh = dense.nvals - 1;
			while (jlow != (jhigh-1))
			{
				jmid = 0.5*(jlow+jhigh);
				if (x <= dense.tabrad[jmid])
				{
					jhigh = jmid;
				}
				else if (x > dense.tabrad[jmid])
				{
					jlow = jmid;
				}
			}
			/* interpolation */
			frac = (x - dense.tabrad[jlow])/(dense.tabrad[jhigh] - dense.tabrad[jlow]);
			tabden_v = dense.tabval[jlow] + frac*(dense.tabval[jhigh] - dense.tabval[jlow]);
		}
	}


	/* got it, now return value, not log of density */
	if ( dense.lgTabLinear )
	{
		tabden_v *= 1.;
	}
	else
	{
		tabden_v = pow(10.,tabden_v);
	}

	return( tabden_v );
}