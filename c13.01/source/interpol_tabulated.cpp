
/* CHANGES: (M. Salz 17.05.2013)
 *  - moved from dense_tabden() to here,
 *    because it is a multiparameter interpolation for
 *    density, temperature and velocity
 *  - add the if clause for interpolation on the linear
 *    values, not the log10 values 
 *  - change table lookup from linear to bisection 
 *  - add if clause for nzone = 0 
 *  - changed interpolation on log depth to interpolation on
 *    on delta log: log(depth) - log(rinner)
 *    --> must import radius.h */

#include "cddefines.h"
#include "radius.h"

double interpol_tabulated( double r0, double depth,
	bool lgDepth, bool lgLinear,
	realnum* tbrad, realnum* tbval,
	long int numvals )
{
	long int jlow, jmid, jhigh;
	double frac, 
	  interp_v, 
	  x;

	DEBUG_ENTRY( "interpol_tabulated()" );
	/*interpolate on table of points for density with dlaw table command, by K Volk
	 *each line is log radius and H density per cc. */

	if( nzone == 0 ){ /* radius values are not set correctly yet */
		interp_v = tbval[0];
	}
	else
	{
		/*begin sanity check */
		if( r0 < 0. || depth <= 0. )
		{
			fprintf( ioQQQ, " interpol_tabulated called with insane depth, radius, =%10.2e%10.2e\n", 
			  depth, r0 );
		}
		/* end sanity check */
		
		/* interpolate linear or on logscale? */
		if ( lgLinear )
		{
			/* interpolate on radius or depth? */
			if( lgDepth )
			{
				x = depth;
			}
			else
			{
				x = r0;
			}
		}
		else
		{
			/* interpolate on radius or depth? */
			if( lgDepth )
			{
				x = log10(depth/radius.rinner);
			}
			else
			{
				x = log10(r0);
			}
		}

		/* set to impossible value, will crash if not reset */
		interp_v = -DBL_MAX;

		if( x < tbrad[0] || x >= tbrad[numvals-1] )
		{
			fprintf( ioQQQ, " requested radius outside range of dense_tabden\n" );
			if ( lgLinear )
			{
				fprintf( ioQQQ, " radius/depth was %10.2e min, max=%10.2e %10.2e\n",
				  x, tbrad[0], tbrad[numvals-1] );
			}
			else
			{
				fprintf( ioQQQ, " log10(radius/depth) was %10.2f min, max=%10.2f %10.2f\n",
				  x, tbrad[0], tbrad[numvals-1] );
			}
			cdEXIT(EXIT_FAILURE);
		}
		else
		{
			/* table lookup */
			jlow = 0;
			jhigh = numvals - 1;
			while (jlow != (jhigh-1))
			{
				jmid = 0.5*(jlow+jhigh);
				if (x <= tbrad[jmid])
				{
					jhigh = jmid;
				}
				else if (x > tbrad[jmid])
				{
					jlow = jmid;
				}
			}
			/* interpolation */
			frac = (x - tbrad[jlow])/(tbrad[jhigh] - tbrad[jlow]);
			interp_v = tbval[jlow] + frac*(tbval[jhigh] - tbval[jlow]);
		}
	}


	/* got it, now return value, not the log */
	if ( lgLinear )
	{
		interp_v *= 1.;
	}
	else
	{
		interp_v = pow(10.,interp_v);
	}

	return( interp_v );
}