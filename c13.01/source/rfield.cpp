/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "rfield.h"

#include "opacity.h"

t_rfield rfield;

const realnum *t_rfield::getCoarseTransCoef()
{
	// average opacity transmission coefficient fine to coarse
	if( opac.lgScatON && trans_coef_total_stale)
	{
		/* sum over coarse continuum */
		for( long i=0; i < nflux-1; i++ )
		{
			// find transmission coefficient if lower and upper bounds 
			// of coarse continuum is within boundaries of fine continuum 
			// unity is default
			if( ipnt_coarse_2_fine[i] && ipnt_coarse_2_fine[i+1] )
			{
				// first branch is normal case, where fine continuum is finer than
				// coarse continuum.  But, when end temp is very high, fine continuum is
				// very coarse, so may be just one cell, and following will not pass
				if( ipnt_coarse_2_fine[i+1]>ipnt_coarse_2_fine[i] )
				{
					trans_coef_total[i] = 0.;
					for( long j=ipnt_coarse_2_fine[i]; j<ipnt_coarse_2_fine[i+1]; ++j )
						trans_coef_total[i] += sexp(fine_opt_depth[j]);
					trans_coef_total[i] /= (ipnt_coarse_2_fine[i+1]-ipnt_coarse_2_fine[i]);
				}
				else
				{
					// in case where fine is coarser than coarse, 
					// just use first cell
					trans_coef_total[i] = sexp(fine_opt_depth[ipnt_coarse_2_fine[i]]);
				}
			}
		}
		trans_coef_total_stale = false;
	}
	return trans_coef_total;
}
