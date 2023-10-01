/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolSave save coolants */
#include "cddefines.h"
#include "thermal.h"
#include "dynamics.h"
#include "radius.h"
#include "conv.h"
#include "phycon.h"
#include "save.h"
#include "grainvar.h"
#include "hmi.h"
#include "coolheavy.h"


/* this is limit to number of coolants to print out */
static const int IPRINT = 100;

/*CoolSave save coolants */
void CoolSave(FILE * io, char chJob[])
{
	long int i, 
	  ip, 
	  is;

	int nFail;

	double cset, 
		cool_total,
		heat_total;

	realnum
		*csav,
		*sgnsav;
	long int *index;

	DEBUG_ENTRY( "CoolSave()" );

	/* cannot do one-time init since thermal.ncltot can change */
	index = (long int *)CALLOC((size_t)thermal.ncltot,sizeof(long int));
	csav = (realnum *)CALLOC((size_t)thermal.ncltot,sizeof(realnum));
	sgnsav = (realnum *)CALLOC((size_t)thermal.ncltot,sizeof(realnum));

	cool_total = thermal.ctot;
	heat_total = thermal.htot;

	/* >>chng 06 mar 17, comment out following block and replace with this 
	 * removing dynamics heating & cooling and report only physical
	 * heating and cooling 
	 * NB the heating and cooling as punched no longer need be
	 * equal for a converged model */
	cool_total -= dynamics.Cool();
	heat_total -= dynamics.Heat();
#	if 0
	if(dynamics.Cool > dynamics.Heat()) 
	{
		cool_total -= dynamics.Heat();
		heat_total -= dynamics.Heat();
	} 
	else
	{
		cool_total -= dynamics.Cool;
		heat_total -= dynamics.Cool;
	}
#	endif

	/* cset will be weakest cooling to consider
	 * WeakHeatCool set with 'set weakheatcool' command
	 * default is 0.05 */
	cset = cool_total*save.WeakHeatCool;

	/* first find all strong lines, both + and - sign */
	ip = thermal.ncltot;

	for( i=0; i < ip; i++ )
	{
		csav[i] = (realnum)(MAX2(thermal.cooling[i],thermal.heatnt[i])/
			SDIV(cool_total));

		/* save sign to remember if heating or cooling line */
		if( thermal.heatnt[i] == 0. )
		{
			sgnsav[i] = 1.;
		}
		else
		{
			sgnsav[i] = -1.;
		}
	}

	/* order strongest to weakest */
	/* now sort by decreasing importance */
	/*spsort netlib routine to sort array returning sorted indices */
	spsort(
		  /* input array to be sorted */
		  csav, 
		  /* number of values in x */
		  ip, 
		  /* permutation output array */
		  index, 
		  /* flag saying what to do - 1 sorts into increasing order, not changing
		   * the original routine */
		  -1, 
		  /* error condition, should be 0 */
		  &nFail);

	/* warn if tcovergence failure occurred */
	if( !conv.lgConvTemp )
	{
		fprintf( io, "#>>>>  Temperature not converged.\n" );
	}
	else if( !conv.lgConvEden )
	{
		fprintf( io, "#>>>>  Electron density not converged.\n" );
	}
	else if( !conv.lgConvIoniz() )
	{
		fprintf( io, "#>>>>  Ionization not converged.\n" );
	}
	else if( !conv.lgConvPres )
	{
		fprintf( io, "#>>>>  Pressure not converged.\n" );
	}

	if( strcmp(chJob,"EACH") == 0 )
	{
		/* begin the print out with zone number, total heating and cooling */
		fprintf( io, "%.5e\t%.4e\t%.4e", 
				 radius.depth_mid_zone, 
				 phycon.te, 
				 cool_total );
		double debug_ctot = 0.;
		
		for( int i=0 ; i <= LIMELM ; i++)
		{
			fprintf( io, "\t%.4e", thermal.elementcool[i] );
			debug_ctot += thermal.elementcool[i];
		}
		
		/*fprintf( io, "\t%.4e", thermal.dust );
		fprintf( io, "\t%.4e", thermal.H2cX );
		fprintf( io, "\t%.4e", thermal.CT_C );
		fprintf( io, "\t%.4e", thermal.H_fb );
		fprintf( io, "\t%.4e", thermal.H2ln );
		fprintf( io, "\t%.4e", thermal.HDro );
		fprintf( io, "\t%.4e", thermal.H2p );
		fprintf( io, "\t%.4e", thermal.FF_c );
		fprintf( io, "\t%.4e", thermal.hvFB );
		fprintf( io, "\t%.4e", thermal.eeff );
		fprintf( io, "\t%.4e", thermal.adve );
		fprintf( io, "\t%.4e", thermal.Comp );
		fprintf( io, "\t%.4e", thermal.Extr );
		fprintf( io, "\t%.4e", thermal.Expn );
		fprintf( io, "\t%.4e", thermal.Cycl );
		fprintf( io, "\t%.4e", thermal.Hvin );
		fprintf( io, " \n" );*/

		fprintf( io, "\t%.4e", MAX2(0.,gv.GasCoolColl) );
		debug_ctot += MAX2(0.,gv.GasCoolColl);
		
		fprintf( io, "\t%.4e", MAX2(0.,-hmi.HeatH2Dexc_used) );
		debug_ctot += MAX2(0.,-hmi.HeatH2Dexc_used);
		
		fprintf( io, "\t%.4e", thermal.char_tran_cool );
		debug_ctot += thermal.char_tran_cool;
		
		fprintf( io, "\t%.4e", hmi.hmicol );
		debug_ctot += hmi.hmicol;
		
		fprintf( io, "\t%.4e", CoolHeavy.h2line );
		debug_ctot += CoolHeavy.h2line;
		
		fprintf( io, "\t%.4e", CoolHeavy.HD );
		debug_ctot += CoolHeavy.HD;
		
		fprintf( io, "\t%.4e", CoolHeavy.H2PlsCool );
		debug_ctot += CoolHeavy.H2PlsCool;
		
		fprintf( io, "\t%.4e", MAX2(0.,CoolHeavy.brems_cool_net) );
		debug_ctot += MAX2(0.,CoolHeavy.brems_cool_net);
		
		fprintf( io, "\t%.4e", CoolHeavy.heavfb );
		debug_ctot += CoolHeavy.heavfb;
		
		fprintf( io, "\t%.4e", CoolHeavy.eebrm );
		debug_ctot += CoolHeavy.eebrm;
		
		fprintf( io, "\t%.4e", CoolHeavy.tccool );
		debug_ctot +=  CoolHeavy.tccool;
		
		fprintf( io, "\t%.4e", CoolHeavy.cextxx );
		debug_ctot += CoolHeavy.cextxx;
		
		fprintf( io, "\t%.4e", CoolHeavy.expans );
		debug_ctot += CoolHeavy.expans;
		
		fprintf( io, "\t%.4e", CoolHeavy.cyntrn );
		debug_ctot += CoolHeavy.cyntrn;
		
		fprintf( io, "\t%.4e", CoolHeavy.colmet );
		debug_ctot += CoolHeavy.colmet;
		
		fprintf( io, "\t%.4e", thermal.dima );
		debug_ctot += thermal.dima;
		fprintf( io, " \n" );
		if( fabs( (debug_ctot - cool_total)/cool_total ) > 1e-10 )
		{
			fprintf( ioQQQ , "PROBLEM with the SAVE EACH COOLING output\n" );
			fprintf( ioQQQ , "PROBLEM One or more coolants have been lost, the sum of the reported cooling is %.4e\n", debug_ctot );
			fprintf( ioQQQ , "PROBLEM The total cooling is %.4ee\n", cool_total );
			fprintf( ioQQQ , "PROBLEM The difference is %.4e\n", cool_total - debug_ctot );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if( strcmp(chJob,"COOL") == 0 )
		{
		/*>>chng 06 jun 06, change start of save to give same info as heating 
		 * as per comment by Yumihiko Tsuzuki */
		/* begin the print out with zone number, total heating and cooling */
		fprintf( io, "%.5e\t%.4e\t%.4e\t%.4e", 
				 radius.depth_mid_zone, 
				 phycon.te, 
				 heat_total, 
				 cool_total );

		/* print only up to IPRINT, which is defined above */
		ip = MIN2( ip , IPRINT );

		/* now print the coolants 
		 * keep sign of coolant, for strong negative cooling 
		 * order is ion, wavelength, fraction of total */
		for( is=0; is < ip; is++ )
		{
			if(is > 4 && (thermal.cooling[index[is]] < cset && thermal.heatnt[index[is]] < cset))
				break;
			fprintf( io, "\t%s %.1f\t%.7f", 
					 thermal.chClntLab[index[is]], 
					 thermal.collam[index[is]], 
					 sign(csav[index[is]],sgnsav[index[is]]) );
		}
		fprintf( io, " \n" );
	}
	else
		TotalInsanity();

	/* finished, now free space */
	free(sgnsav);
	free(csav);
	free(index);

	return;
}

