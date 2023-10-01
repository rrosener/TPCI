
/* CHANGES: (M. Salz 17.05.2013)
 *  - moved from parse_dlaw to here, because it is used for reading
      tabulated density, temperature and velocity
    - introduced the linear keyword for linear interpolation
    - added sanity checks to prevent calculations if invalid
      tables ranges are provided by the user */

#include "cddefines.h"
#include "radius.h"
#include "parser.h"
#include "dense.h"

void ParseTabulated(Parser &p, bool* lgDepth, bool* lgLinear, realnum* tbrad, realnum* tbval,long int* numvals )
{
	bool lgEnd;

	DEBUG_ENTRY( "ParseTabulated()" );
	
	if( p.nMatch("DEPT") )
	{
		*lgDepth = true;
	}
	else
	{
		*lgDepth = false;
	}

	/* set true if linear interpolation to be used */
	if( p.nMatch("LINE") )
	{
		*lgLinear = true;
	}
	else
	{
		*lgLinear = false;
	}

	p.getline();
	tbrad[0] = (realnum)p.FFmtRead();
	tbval[0] = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		fprintf( ioQQQ, " No pairs entered - can\'t interpolate.\n Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* sanity check for first point in dlaw table */
	if( *lgDepth )
	{
		if( *lgLinear )
		{
			if( abs(tbrad[0]) > 0.000001 )
			{
				fprintf( ioQQQ, " First point on dlaw linear depth table should be zero.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			if( tbrad[0] >= -30 )
			{
				fprintf( ioQQQ, " First point on dlaw logarithmic depth table shoule be < -30.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	else
	{
		if( radius.lgRadiusKnown )
		{
			if( *lgLinear )
			{
				if( tbrad[0] >= radius.Radius )
				{
					fprintf( ioQQQ, " First point on dlaw radius table shoule be smaller than starting radius. Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
				else if( tbrad[0] < 0.01*radius.Radius )
				{
					fprintf( ioQQQ, " First point on dlaw radius table was much smaller than starting radius. Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else 
			{
				if( pow(10.,tbrad[0]) >= radius.Radius )
				{
					fprintf( ioQQQ, " First point on dlaw radius table shoule be smaller than starting radius. Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
				else if( pow(10.,tbrad[0]) < 0.01*radius.Radius )
				{
					fprintf( ioQQQ, " First point on dlaw radius table was much smaller than starting radius. Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
		}
	}

	*numvals = 2;
	lgEnd = false;

	/* read pairs of numbers until we find line starting with END */
	/* >>chng 04 jan 27, loop to LIMTABDLAW from LIMTABD, as per
	 * var definitions, caught by Will Henney */
	while( !lgEnd && *numvals < LIMTABDLAW )
	{
		p.getline();
		lgEnd = p.m_lgEOF;
		if( !lgEnd )
		{
			if( p.strcmp("END") == 0 )
				lgEnd = true;
		}

		if( !lgEnd )
		{
			tbrad[*numvals-1] = (realnum)p.FFmtRead();
			tbval[*numvals-1] = (realnum)p.FFmtRead();
			*numvals += 1;
		}
	}
	*numvals -= 1;

	for( long i=1; i < *numvals; i++ )
	{
		/* the radius values are assumed to be strictly increasing */
		if( tbrad[i] <= tbrad[i-1] )
		{
			fprintf( ioQQQ, " Radii must be in increasing order.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* sanity check for last point in dlaw table if thickness of cloud is defined */
	if( radius.StopThickness[0] < 1.e30 )
	{
		if( *lgDepth )
		{
			if( *lgLinear )
			{
				if( tbrad[*numvals - 1] <= radius.StopThickness[0] )
				{
					fprintf( ioQQQ, " Defined thickness is larger than dlaw table.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else
			{
				if( pow(10.,tbrad[*numvals - 1]) <= radius.StopThickness[0] )
				{
					fprintf( ioQQQ, " Defined thickness is larger than dlaw table.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
		}
		else
		{
			if( *lgLinear )
			{
				if( (tbrad[*numvals - 1] - radius.Radius) <= radius.StopThickness[0] )
				{
					fprintf( ioQQQ, " Defined thickness is larger than dlaw table.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else
			{
				if( (pow(10.,tbrad[*numvals - 1])- radius.Radius) <= radius.StopThickness[0] )
				{
					fprintf( ioQQQ, " Defined thickness is larger than dlaw table.  Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
		}
	}
	
	return;
}