/* This file is part of Cloudy and is copyright (C)1978-2010 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program calling cloudy to produce a table giving ionization vs temperature */
#include "cddefines.h"
#include "cddrive.h"

int main( void )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	try {
#		define NMAX 100
		/* NELEMI and NELEMF define the range of elements for which the cooling will be calculated */
#		define NELEMI 0
#		define NELEMF LIMELM
		double xCoolSave[NMAX][LIMELM], tesave[NMAX];
		double telog , teinc, te_last, te_first;
		double hden, element_abund ;
		long int nte , i;
		long int nelem;
		/* this is the list of element names used to make printout */
		char chElementName[LIMELM][12] =
			{ "Hydrogen\0" ,"Helium\0" ,"Lithium\0" ,"Beryllium\0" ,"Boron\0" ,
			  "Carbon\0" ,"Nitrogen\0" ,"Oxygen\0" ,"Fluorine\0" ,"Neon\0" ,
			  "Sodium\0" ,"Magnesium\0" ,"Aluminium\0" ,"Silicon\0" ,"Phosphorus\0" ,
			  "Sulphur\0" ,"Chlorine\0" ,"Argon\0" ,"Potassium\0" ,"Calcium\0" ,
			  "Scandium\0" ,"Titanium\0" ,"Vanadium\0" ,"Chromium\0" ,"Manganese\0" ,
			  "Iron\0" ,"Cobalt\0" ,"Nickel\0" ,"Copper\0" ,"Zinc\0"};

		FILE *ioRES ;
		char chLine[100];

		/* this will be limit to the number of command chLines we can still put in */
		long int nleft;

		/* calculation's results are saved here */
		if( (ioRES = fopen("collcool.txt","w")) == NULL )
		{
			printf(" could not open collcool.txt for writing.\n");
			cdEXIT(EXIT_FAILURE);
		}

		for( i = 0; i < NMAX; i++ )
		{
			for( nelem=0; nelem<LIMELM; ++nelem)
			{
				xCoolSave[i][nelem] = -1.;
			}
		}

		/* the first temperature */
		te_first = 4.;
		te_last = 9.;
		/* the increment in the temperature */
		teinc = 0.2;
		/* the log of the hydrogen density that will be used */
		hden = 0.;
		/* the log of the ratio of the abundance of the element in question to hydrogen */
		element_abund = 5.;

		int npoints = (te_last-te_first)/teinc + 1;
		

		for( nelem=NELEMI; nelem<NELEMF; ++nelem)
		{
			nte = 0;
			telog = te_first;

			while( telog <= te_last+0.01 && nte<NMAX)
			{
				/* initialize the code for this run */
				cdInit();
				cdTalk(false);
				/*cdNoExec( );*/
				printf("%s  Te = %g\n",chElementName[nelem],telog);
				
				sprintf(chLine,"coronal %3.1f ",telog);
				nleft = cdRead( chLine  );

				nleft = cdRead( "atom chianti hybrid levels max \"CloudyChiantiKurucz.ini\" "  );

				nleft = cdRead( "atom feii"  );

				nleft = cdRead( "no co molecules"  );

				nleft = cdRead( "atom lamda off"  );

				/* just do the first zone - only want ionization distribution */
				nleft = cdRead( "stop zone 1 "  );

				/* the hydrogen density entered as a log */
				sprintf(chLine,"hden %i ",int(hden));
				nleft = cdRead( chLine  );

				nleft = cdRead( "set dr 0"  );

				// cannot turn elements off and on during one run, so set them
				// all to small abundances, relative to H
				nleft = cdRead( "element helium abundance -9"  );
				nleft = cdRead( "element lithium abundance -9"  );
				nleft = cdRead( "element beryllium abundance -9"  );
				nleft = cdRead( "element boron abundance -9"  );
				nleft = cdRead( "element carbon abundance -9"  );
				nleft = cdRead( "element nitrogen abundance -9"  );
				nleft = cdRead( "element oxygen abundance -9"  );
				nleft = cdRead( "element fluorine abundance -9"  );
				nleft = cdRead( "element neon abundance -9"  );
				nleft = cdRead( "element Sodium abundance -9"  );
				nleft = cdRead( "element Magnesium abundance -9"  );
				nleft = cdRead( "element Aluminium abundance -9"  );
				nleft = cdRead( "element Silicon abundance -9"  );
				nleft = cdRead( "element Phosphorus abundance -9"  );
				nleft = cdRead( "element Sulphur abundance -9"  );
				nleft = cdRead( "element Chlorine abundance -9"  );
				nleft = cdRead( "element Argon abundance -9"  );
				nleft = cdRead( "element Potassium abundance -9"  );
				nleft = cdRead( "element Calcium abundance -9"  );
				nleft = cdRead( "element Scandium abundance -9"  );
				nleft = cdRead( "element Titanium abundance -9"  );
				nleft = cdRead( "element Vanadium abundance -9"  );
				nleft = cdRead( "element Chromium abundance -9"  );
				nleft = cdRead( "element Manganese abundance -9"  );
				nleft = cdRead( "element Iron abundance -9"  );
				nleft = cdRead( "element Cobalt abundance -9"  );
				nleft = cdRead( "element Nickel abundance -9"  );
				nleft = cdRead( "element Copper abundance -9"  );
				nleft = cdRead( "element Zinc abundance -9"  );

				sprintf(chLine,"element %s abundance %f ",chElementName[nelem],element_abund);
				nleft = cdRead( chLine  );

				nleft = cdRead( "set eden 0"  );

				/* actually call the code */
				if( cdDrive() )
					exit_status = ES_FAILURE;

				xCoolSave[nte][nelem] = cdCooling_last();
				if( nelem != 0 )
				{
					xCoolSave[nte][nelem] /= pow(10.,element_abund);
				}
				//printf("%6.2e\t\%6.2e\n",xCoolSave[nte][nelem],cdCooling_last());

				tesave[nte] = telog;
				telog += teinc;
				++nte;
			}
		}
			
		/* this generates large printout */
		fprintf(ioRES,"Te");
		for( nelem=NELEMI; nelem<NELEMF; ++nelem)
		{
			fprintf(ioRES,"\t%s",chElementName[nelem]);
		}
		fprintf(ioRES,"\n");
		
		for( i = 0; i < npoints; i++ )
		{
			fprintf(ioRES,"%6.2e",pow(10.,tesave[i]));
			for( nelem=NELEMI; nelem<NELEMF; ++nelem)
			{
				fprintf(ioRES,"\t%6.2e",xCoolSave[i][nelem]);
			}
			fprintf(ioRES,"\n");
		}
		fclose(ioRES);
		cdEXIT(exit_status);
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT || e.sig() == SIGQUIT )
			{
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
				exit_status = ES_USER_INTERRUPT;
			}
			else if( e.sig() == SIGTERM )
			{
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
				exit_status = ES_TERMINATION_REQUEST;
			}
			else if( e.sig() == SIGILL )
			{
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
				exit_status = ES_ILLEGAL_INSTRUCTION;
			}
			else if( e.sig() == SIGFPE )
			{
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
				exit_status = ES_FP_EXCEPTION;
			}
			else if( e.sig() == SIGSEGV )
			{
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
				exit_status = ES_SEGFAULT;
			}
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
			{
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
				exit_status = ES_BUS_ERROR;
			}
#			endif
			else
			{
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
				exit_status = ES_UNKNOWN_SIGNAL;
			}

		}
	}
#endif
	catch( cloudy_exit& e )
	{
		if( ioQQQ != NULL )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	cdPrepareExit(exit_status);

	return exit_status;
}
