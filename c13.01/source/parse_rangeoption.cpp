/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseRangeOption parse the range option on the luminosity command */
#include "cddefines.h"
#include "physconst.h"
#include "rfield.h"
#include "parser.h"

void ParseRangeOption(
	/* the input line image */
	Parser &p)
{
	bool lgP1Absent, lgP2Absent;
	double p1, 
	  p2;

	DEBUG_ENTRY( "ParseRangeOption()" );

	if( p.nMatch("TOTA") )
	{
		rfield.range[p.m_nqh][0] = rfield.emm;
		rfield.range[p.m_nqh][1] = rfield.egamry;
	}
	else if( p.nMatch("RANG") )
	{
		p.set_point(0);
		/* first call return the luminosity on the card, ignore it */
		p1 = p.FFmtRead();
		lgP1Absent = p.lgEOL();

		/* read lower limit and upper limit*/
		p1 = p.FFmtRead();
		lgP1Absent = p.lgEOL();
		p2 = p.FFmtRead();
		lgP2Absent = p.lgEOL();

		/* option to enter log if first energy is neg */
		if( p1 < 0. || p.nMatch(" LOG") )
		{
			p1 = pow(10.,p1);
			p2 = pow(10.,p2);
		}

		if( lgP1Absent )
			p1 = rfield.emm;
		if( lgP2Absent )
			p2 = rfield.egamry;

		/* make sure that energies are within array bounds */
		rfield.range[p.m_nqh][0] = MAX2((realnum)p1,rfield.emm);
		rfield.range[p.m_nqh][1] = MIN2((realnum)p2,rfield.egamry);
		if( rfield.range[p.m_nqh][0] >= rfield.range[p.m_nqh][1] )
		{
			fprintf( ioQQQ, " Range MUST be in increasing order - sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		/* range not specified, use default - total ionizing field
		 * >>chng 96 dec 18, from 1.001 to 1 Ryd for H mass nuc */
		rfield.range[p.m_nqh][0] = HIONPOT;
		rfield.range[p.m_nqh][1] = rfield.egamry;
	}
	return;
}
