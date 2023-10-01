/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*PrtColumns print column densities of all elements */
#include "cddefines.h"
#include "cddrive.h"
#include "dense.h"
#include "elementnames.h"
#include "h2.h"
#include "save.h"
#include "taulines.h"
#include "molcol.h"
#include "prt.h"

void PrtColumns(
	/* this is stream used for io, is stdout when called by final,
	 * is save unit when save output generated */
	 FILE *ioMEAN ,
	 /* PRETTY for nice looking table in main output
	  * TABLE for table for save column density output, feeding into post processors
	  */
	 const char *chType ,
	 /* index of save in save array */
	 long int ioPun )
{
	long int i, 
	  nelem;

	double aa;

	DEBUG_ENTRY( "PrtColumns()" );

	if( strcmp( chType , "PRETTY")==0 )
	{
		/* print molecular element column densities */
		molcol("PRIN" , ioMEAN);

		fprintf( ioMEAN, "\n" );

		fprintf( ioMEAN, "\n         " );
		for( i=1; i <= 17; i++ )
		{
			fprintf( ioMEAN, "%7ld", i );
		}
		fprintf( ioMEAN, "\n\n" );

		/* ionization column densities  */
		for( nelem=0; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				bool lgDONE = false;

				fprintf( ioMEAN, " %10.10s", elementnames.chElementName[nelem]  );

				i = 1;
				while( !lgDONE )
				{
					if( cdColm(
						/* return value is zero if all ok, 1 if errors happened */
						/* 4-char + eol string that is first
						* 4 char of element name as spelled by cloudy */
						elementnames.chElementNameShort[nelem],

						/* integer stage of ionization, 1 for atom, 0 for CO or H2 */
						i,

						/* the theoretical column density derived by the code */
						&aa ) )
						TotalInsanity();

					if( aa == 0. )
					{
						aa = -30.;
					}
					else if( aa > 0. )
					{
						aa = log10(aa);
					}

					if( i == 18 )
					{
						fprintf( ioMEAN, "\n" );
					}
					fprintf( ioMEAN, "%7.3f", aa );

					/* increment counter and check if at upper limit */
					++i;
					/* MAX2 is to include H2 in H array */
					if( i > MAX2(3,nelem+2) )
						lgDONE = true;

					/* print title for this info if we are done with hydrogen */
					if( nelem==ipHYDROGEN && lgDONE )
						fprintf(ioMEAN," (H2)                Log10 Column density (cm^-2)");
				}

				fprintf( ioMEAN, "\n" );
			}
		}

		/* only print excited state column densities if level2 lines are included
		 * since they populated the upper level by UV pumping.   This process
		 * is not included if level2 lines are not considered, by introducing
		 * the "no level2" command */
		if( nWindLine>0 )
		{
			const int nEXCIT_COL = 12;
			// use strings instead of char[][] to avoid bogus errors with PGCC bounds checking
			string strExcit_Col[nEXCIT_COL]={
				"He1*","CII*","C11*","C12*","C13*","O11*","O12*","O13*","Si2*","C30*","C31*","C32*"};
			long int nprt = 0;
			/* print  excited level column densities */
			fprintf(ioMEAN," Exc state ");
			nprt = 12;
			for(i=0; i<nEXCIT_COL; ++i )
			{
				ASSERT( strExcit_Col[i].size() <= 4 );
				if( cdColm(
					/* return value is zero if all ok, 1 if errors happened */
					/* 4-char + eol string that is first 
					 * 4 char of element name as spelled by cloudy */
					strExcit_Col[i].c_str(),
					/* integer stage of ionization, 1 for atom, 0 for CO, OH, CII*, or H2 */
					0,
					/* the theoretical column density derived by the code */
					&aa ) )
						TotalInsanity();

				if( nprt > 120 )
				{
					fprintf(ioMEAN,"\n           ");
					nprt = 0;
				}
				fprintf(ioMEAN,"   %s%7.3f",
					strExcit_Col[i].c_str(),
					log10(SDIV(aa) ));
				nprt += 14;
			}
			fprintf(ioMEAN,"\n");
		}

		/* print column densities for H2 */
		h2.H2_Prt_column_density(ioMEAN);

		fprintf(ioMEAN,"\n");
	}
	else if( strcmp( chType , "TABLE")==0 )
	{
		ASSERT( ioPun>=0 );
		if( save.lgPunHeader[ioPun] )
		{
			fprintf(ioMEAN,"#Column densities %s",
				mole_global.list[0]->label.c_str() );
			for( i=1; i < mole_global.num_calc; i++ )
			{
				fprintf(ioMEAN,"\t%s",
					mole_global.list[i]->label.c_str() );
			}
			fprintf(ioMEAN,"\n");
			save.lgPunHeader[ioPun] = false;
		}
		fprintf(ioMEAN,"%.3e",
			mole.species[0].column	);
		for( i=1; i < mole_global.num_calc; i++ )
		{
			fprintf(ioMEAN,"\t%.3e",
				mole.species[i].column	);
		}
		fprintf(ioMEAN,"\n");
	}
	else
		TotalInsanity();
	return;
}
