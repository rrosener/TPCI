/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ipoint returns pointer to any energy within energy mesh */
/*ipFine()Cont returns array index within fine energy mesh */
/*ipContEnergy generate pointer to energy within continuum array
 * continuum energy in Rydbergs */
/*ipLineEnergy generate pointer to line energy within energy mesh
 * line energy in Rydbergs */
#include "cddefines.h"
#include "continuum.h"
#include "prt.h"
#include "rfield.h"
#include "ipoint.h"

/*ipoint returns pointer to any energy within energy mesh */
long ipoint(double energy_ryd)
{
	long int i, 
	  ipoint_v;

	DEBUG_ENTRY( "ipoint()" );

	// make sure mesh is set up
	ASSERT( continuum.nrange > 0 );

	if( energy_ryd < continuum.filbnd[0] || energy_ryd > continuum.filbnd[continuum.nrange] )
	{
		fprintf( ioQQQ, " ipoint:\n" );
		fprintf( ioQQQ, " The energy_ryd array is not defined at nu=%11.3e. The bounds are%11.3e%11.3e\n", 
		  energy_ryd, continuum.filbnd[0], continuum.filbnd[continuum.nrange] );
		fprintf( ioQQQ, " ipoint is aborting to get trace, to find how this happened\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	for( i=0; i < continuum.nrange; i++ )
	{
		if( energy_ryd >= continuum.filbnd[i] && energy_ryd <= continuum.filbnd[i+1] )
		{

			/* this is on the Fortran scale of array index counting, so is one above the
			 * c scale.  later the -1 will appear on any index references */
			ipoint_v = (long int)(log10(energy_ryd/continuum.filbnd[i])/continuum.fildel[i] + 
			  1.0 + continuum.ifill0[i]);

			ASSERT( ipoint_v >= 0 );
			/* recall on F scale */
			ipoint_v = MIN2( rfield.nupper , ipoint_v );

			if( ipoint_v < rfield.nflux-2 && ipoint_v>2 )
			{
				// possible will need to adjust index if cells have been fiddled with
				if( energy_ryd > rfield.anu[ipoint_v-1]+rfield.widflx[ipoint_v-1]/2. )
					++ipoint_v;
				if( energy_ryd < rfield.anu[ipoint_v-1]-rfield.widflx[ipoint_v-1]/2. )
					--ipoint_v;
				{
					enum {DEBUG_LOC=false};
					if( DEBUG_LOC )
					{
						if( energy_ryd > rfield.anu[ipoint_v]+rfield.widflx[ipoint_v]/2. )
						{
							fprintf(ioQQQ,"DISASTER ipoint hi bounds about to throw "\
								  "energy_ryd=%e hi bound=%e ipoint_v=%li nflux=%li\n",
								  energy_ryd , rfield.anu[ipoint_v]+rfield.widflx[ipoint_v]/2.,
								  ipoint_v,rfield.nflux);
						}
						if( energy_ryd < rfield.anu[ipoint_v-2]-rfield.widflx[ipoint_v-2]/2. )
						{
							fprintf(ioQQQ,"DISASTER ipoint low bounds about to throw "\
								  "energy_ryd=%e low bound=%e ipoint_v=%li nflux=%li\n",
								  energy_ryd , rfield.anu[ipoint_v-2]-rfield.widflx[ipoint_v-2]/2.,
								  ipoint_v,rfield.nflux);
						}
					}
				}

				ASSERT( energy_ryd <= rfield.anu[ipoint_v]+rfield.widflx[ipoint_v]/2. );
				ASSERT( energy_ryd >= rfield.anu[ipoint_v-2]-rfield.widflx[ipoint_v-2]/2. );
			}
			return ipoint_v;
		}
	}

	/* this exit not possible, here to shut up some compilers */
	fprintf( ioQQQ, " IPOINT logic error, energy=%.2e\n", 
	  energy_ryd );
	cdEXIT(EXIT_FAILURE);
}

/*ipContEnergy generate pointer to energy within continuum array */
long ipContEnergy(
  /* continuum energy in Rydbergs */
  double energy, 
  /* 4 char label for continuum, like those returned by chLineLbl */
  const char *chLabel)
{
	long int ipConSafe_v;

	DEBUG_ENTRY( "ipContEnergy()" );

	ipConSafe_v = ipoint(energy);

	/* write label in this cell if not anything there yet */
	if( strcmp(rfield.chContLabel[ipConSafe_v-1],"    ") == 0 )
	{
		strcpy( rfield.chContLabel[ipConSafe_v-1], chLabel );
	}

	/* this is a quick way to find all continua that occur within a given cell,
	 * recall the off-by-one aspect of C */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			/* recall the off-by-one aspect of C - the number below is
			 * on the physics scale because this returns the index
			 * on that scale, so the first is 1 for [0] */
			if( ipConSafe_v == 23 )
				fprintf(ioQQQ,"%s\n", chLabel );
		}
	}
	return ipConSafe_v;
}

/*ipLineEnergy generate pointer to line energy within energy mesh
 * line energy in Rydbergs -
 * return value is array index on the physics or Fortran scale, counting from 1 */
long ipLineEnergy(double energy, 
  /* 4 char label for  line, like those returned by chLineLbl */
  const char *chLabel , 
  /* will make sure energy is < this array index if greater than 0, does nothing if <= 0*/
  long ipIonEnergy )
{
	long int ipLine_ret;

	DEBUG_ENTRY( "ipLineEnergy()" );

	ipLine_ret = ipoint(energy);
	ASSERT( ipLine_ret );
	/* make sure pointer is below next higher continuum if positive */
	if( ipIonEnergy > 0 )
	{
		ipLine_ret = MIN2( ipLine_ret , ipIonEnergy-1 );
	}

	ASSERT( ipLine_ret > 0 );
	/* stuff in a label if none there,
	 * note that this is offset one below actual number to be on C scale of arrays */
	/* >>chng 06 nov 23, faster to use line_count index rather than checking 5 chars 
	 * first call will have zero so false */
	/*if( strcmp(rfield.chLineLabel[ipLine_ret-1],"    ")==0 )*/
	if( !rfield.line_count[ipLine_ret-1] )
	{
		strcpy( rfield.chLineLabel[ipLine_ret-1], chLabel );
	}
	/* this keeps track of the number of lines within this cell */
	++rfield.line_count[ipLine_ret-1];

	/* this is a quick way to find all lines that occur within a given cell,
	 * recall the off-by-one aspect of C */
	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			/* recall the off-by-one aspect of C - the numbers is one more 
			 * than the index on the C scale */
			if( ipLine_ret == 23 )
				fprintf(ioQQQ,"%s\n", chLabel );
		}
	}

	/* this implements the print continuum indices command */
	if( prt.lgPrtContIndices )
	{
		/* print header if first time */
		static bool lgFirst = true;
		if( lgFirst )
		{
			/* print header and set flag so do not do again */
			fprintf(ioQQQ , "\n\noutput from print continuum indices command follows.\n");
			fprintf(ioQQQ , "cont ind (F scale)\tenergy(ryd)\tlabel\n");
			lgFirst = false;
		}
		if( energy >= prt.lgPrtContIndices_lo_E && energy <= prt.lgPrtContIndices_hi_E )
		{
			/* use varying formats depending on order of magnitude of energy 
			 * NB - the ipLine_ret is the index on the physical or Fortran scale,
			 * and is 1 greater than the array element used in the code, which is
			 * on the C scale */
			if( energy < 1. )
			{
				fprintf(ioQQQ , "%li\t%.3e\t%s\n" , ipLine_ret , energy , chLabel);
			}
			else if( energy < 10. )
			{
				fprintf(ioQQQ , "%li\t%.3f\t%s\n" , ipLine_ret , energy , chLabel);
			}
			else if( energy < 100. )
			{
				fprintf(ioQQQ , "%li\t%.2f\t%s\n" , ipLine_ret , energy , chLabel);
			}
			else
			{
				fprintf(ioQQQ , "%li\t%.1f\t%s\n" , ipLine_ret , energy , chLabel);
			}
		}
	}

	if( prt.lgPrnLineCell )
	{
		/* print line cell option - number on command line is cell on Physics scale */
		if( prt.nPrnLineCell == ipLine_ret )
		{
			static bool lgMustPrintHeader = true;
			if( lgMustPrintHeader )
				fprintf(ioQQQ, "Lines within cell %li (physics scale) \nLabel\tEnergy(Ryd)\n",prt.nPrnLineCell );
			lgMustPrintHeader = false;
			fprintf(ioQQQ,"%s\t%.3e\n" , chLabel , energy );
		}
	}
	return ipLine_ret;
}

/*ipFine()Cont returns array index within fine energy mesh */
long ipFineCont(
	/* energy in Ryd */
	double energy_ryd )
{
	long int ipoint_v;

	DEBUG_ENTRY( "ipFine()Cont()" );

	if( energy_ryd < rfield.fine_ener_lo || energy_ryd > rfield.fine_ener_hi )
	{
		return -1;
	}

	/* this is on the Fortran scale of array index counting, so is one above the
	 * c scale.  later the -1 will appear on any index references
	 * 
	 * 07 Jun 22 testing done to confirm that energy grid is correct:  did 
	 * same sim with standard fine continuum resolution, and another with 200x
	 * finer resolution, and confirmed that these line up correctly.  */
	ipoint_v = (long int)(log10(energy_ryd*(1.-rfield.fine_resol/2.) /
		rfield.fine_ener_lo)/log10(1.+rfield.fine_resol));

	ASSERT( ipoint_v >= 0 && ipoint_v< rfield.nfine_malloc );
	return ipoint_v;
}
