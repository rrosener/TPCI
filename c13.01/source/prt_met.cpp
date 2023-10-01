/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*prtmet print all line optical depths at end of iteration */
#include "cddefines.h"
#include "taulines.h"
#include "h2.h"
#include "iso.h"
#include "lines_service.h"
#include "dense.h"
#include "prt.h"
#include "mole.h"
#include "trace.h"

/*prtmet print all line optical depths at end of iteration */
void prtmet(void)
{
	long int i,
		nelem , 
		ipHi , 
		ipLo , 
		ipISO;

	DEBUG_ENTRY( "prtmet()" );

	/* default is to not print optical depths, turn on with
	 * print optical depths on command */
	if( prt.lgPrtTau  || (trace.lgTrace && trace.lgOptcBug) )
	{
		fprintf( ioQQQ, "\n\n                                                 Mean Line Optical Depths\n");

		// initialize  -- TauLines[0] is just a dummy here 
		prme(true, TauLines[0]);

		/* iso sequences */
		for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
		{
			for( nelem=ipISO; nelem < LIMELM; nelem++ )
			{
				if( dense.lgElmtOn[nelem] )
				{
					/* print Lyman, Balmer, Paschen, etc sequence optical depths */	
					for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local-1; ipLo++ )
					{
						for( ipHi=ipLo+1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
						{
							prme(false,iso_sp[ipISO][nelem].trans(ipHi,ipLo));
						}
					}
				}
			}
		}

		/* print main lines optical depths */
		for( i=1; i <= nLevel1; i++ )
		{
			prme(false,TauLines[i]);
		}

		for( i=0; i < nWindLine; i++ )
		{
			if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
			{
				prme(false,TauLine2[i]);
			}
		}

		for( i=0; i < nUTA; i++ )
		{
			prme(false,UTALines[i]);
		}

		/* print H2 line optical depths */
		for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
			(*diatom)->H2_Prt_line_tau();

		for( i=0; i < nHFLines; i++ )
		{
			prme(false,HFLines[i]);
		}

		/* data base lines */
		for (int ipSpecies=0; ipSpecies < nSpecies; ++ipSpecies)
		{
			for( EmissionList::iterator em=dBaseTrans[ipSpecies].Emis().begin();
				  em != dBaseTrans[ipSpecies].Emis().end(); ++em)
			{
				prme(false,(*em).Tran());
			}
		}

		fprintf( ioQQQ, "\n");
	}
	return;
}

/* prme - print line optical depth */
void prme(
  const bool lgReset,
  const TransitionProxy &t)
{
	static long int n ;

	DEBUG_ENTRY( "prme()" );

	if( lgReset )
		n = 0;

	if( t.ipCont() <= 0 )
	{
		/* line is not transferred */
		return;
	}

	/* print optical depth if greater than lower limit, or significantly negative
	 * PrtTauFnt is threshold for printing it
	 * */
	if( t.Emis().TauIn()*SQRTPI > prt.PrtTauFnt || t.Emis().TauIn()*SQRTPI < -1e-5 )
	{
		fprintf( ioQQQ, "  %10.10s",chLineLbl(t));
		/*>> chng 12 jul 25, print mean optical depths, rather than line center */
		fprintf( ioQQQ, PrintEfmt("%9.2e", t.Emis().TauIn()*SQRTPI ));

		// throw CR after printing 6 numbers
		++n;
		if(n == 6)
		{
			n = 0;
			fprintf( ioQQQ, " \n");
		}
	}

	return;
}
