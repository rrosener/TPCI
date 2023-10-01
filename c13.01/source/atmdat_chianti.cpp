/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines_service.h"
#include "physconst.h"
#include "taulines.h"
#include "trace.h"
#include "string.h"
#include "thirdparty.h"
#include "rfield.h"
#include "atmdat.h"

static const bool DEBUGSTATE = false;
// minimum energy difference (wavenumbers) we will accept
const double ENERGY_MIN_WN = 1e-10;

void atmdat_STOUT_readin( long intNS, char *chPrefix )
{
	DEBUG_ENTRY( "atmdat_STOUT_readin()" );

	long int nMolLevs;
	char chLine[FILENAME_PATH_LENGTH_2],
		chNRGFilename[FILENAME_PATH_LENGTH_2],
		chTPFilename[FILENAME_PATH_LENGTH_2],
		chCOLLFilename[FILENAME_PATH_LENGTH_2];

	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLTE = false;
	dBaseSpecies[intNS].lgLAMDA = false;

	strcpy( chNRGFilename , chPrefix );
	strcpy( chTPFilename , chPrefix );
	strcpy( chCOLLFilename , chPrefix );

	/*Open the energy levels file*/
	strcat( chNRGFilename , ".nrg");
	uncaps( chNRGFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chNRGFilename);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_STOUT_readin opening %s:",chNRGFilename);

	FILE *ioDATA;
	ioDATA = open_data( chNRGFilename, "r" );
	long int i = 0;
	bool lgEOL = false;
	long index = 0;
	double nrg = 0.0;
	double stwt = 0.0;

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chNRGFilename );
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	const long int iyr = 11, imo=10 , idy = 14;
	long iyrread, imoread , idyread;
	iyrread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chNRGFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	nMolLevs = 0;
	//Count number of levels
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* # - comment, *** ends data */
		if( chLine[0] != '#' && chLine[0] != '\n' && chLine[0] != '*' )
		{
			i = 1;
			long n = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			if( n < 0 )
				break;
			nMolLevs++;
		}
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not rewind %s.\n", chNRGFilename );
		cdEXIT(EXIT_FAILURE);
	}
	//Skip the magic numbers this time
	read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );

	long HighestIndexInFile = nMolLevs;

	nMolLevs = MIN3(atmdat.nStoutMaxLevels,nMolLevs, MAX_NUM_LEVELS );

	if( atmdat.lgStoutPrint == true)
	{
		fprintf( ioQQQ,"     Using STOUT model %s with %li levels of %li available.\n",
				dBaseSpecies[intNS].chLabel , nMolLevs , HighestIndexInFile );
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	/*Resize the States array*/
	dBaseStates[intNS].resize(nMolLevs);
	/*Zero out the maximum wavenumber for each species */
	dBaseSpecies[intNS].maxWN = 0.;

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = "Stout";

	//This is creating transitions that we don't have collision data for
	//Maybe use gbar or keep all of the Fe 2 even if it was assumed (1e-5)
	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}


	/* Check for an end of file sentinel */
	bool lgSentinelReached = false;

	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the energy level file.\n");
		cdEXIT(EXIT_FAILURE);
	}
	//Read the remaining lines of the energy level file
	do
	{
		i = 1;

		/* Stop on *** */
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}
		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			/* Skip blank lines */
			if( chLine[0] == '\n')
				continue;

			index = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) -1 ;

			if( index < 0 )
			{
				fprintf( ioQQQ, " PROBLEM An energy level index was less than 1 in file %s\n",chNRGFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			nrg = (double)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			stwt = (double)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

			if( lgEOL )
			{
				fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chNRGFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			if( index >= nMolLevs )
			{
				//Skip unused levels
				continue;
			}

			dBaseStates[intNS][index].energy().set(nrg,"cm^-1");
			dBaseStates[intNS][index].g() = stwt;

			/*information needed for label*/
			strcpy(dBaseStates[intNS][index].chLabel(), "    ");
			strncpy(dBaseStates[intNS][index].chLabel(),dBaseSpecies[intNS].chLabel, 4);
			dBaseStates[intNS][index].chLabel()[4] = '\0';
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chNRGFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column of that line.\n");
		cdEXIT(EXIT_FAILURE);
	}
	fclose(ioDATA);

	if( DEBUGSTATE )
	{
		/*Test whether energy levels are read in properly*/
		printf("DEBUG STOUT ENERGY LEVELS:\n");
		for( int i = 0; i< nMolLevs; i++ )
		{
			printf("DEBUG\t%i\t%11.4e\t%3.1f\n",i+1,dBaseStates[intNS][i].energy().WN(),dBaseStates[intNS][i].g());
		}
	}

	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	double fenergyWN = 0.;
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = (*tr).ipHi();
		int ipLo = (*tr).ipLo();
		ASSERT(ipHi > ipLo );
		fenergyWN = dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN();
		if( fenergyWN <= 0.)
		{
			printf("\nWARNING: The %s transition between level %i and %i has zero or negative energy.\n",
					dBaseStates[intNS][ipHi].chLabel(),ipLo+1,ipHi+1);
			printf("Check the Stout energy level data file (*.nrg) of this species for missing or duplicate energies.\n");
			//cdEXIT(EXIT_FAILURE);
		}
		(*tr).EnergyWN() = (realnum)MAX2(ENERGY_MIN_WN,fenergyWN);
		(*tr).WLAng() = (realnum)(1e+8/(*tr).EnergyWN()/RefIndex((*tr).EnergyWN()));
		dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,(*tr).EnergyWN());
	}

	/******************************************************
	 ************* Transition Probability File ************
	 ******************************************************/
	strcat( chTPFilename , ".tp");
	uncaps( chTPFilename );

	ioDATA = open_data( chTPFilename, "r" );

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chTPFilename );
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	iyrread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chTPFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	long numtrans = 0;
	//Count number of transitions
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* # is comment, *** is end of data */
		if( chLine[0] != '#' && chLine[0] != '\n' && chLine[0] != '*' )
		{
			i = 1;
			long n = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
			if( n < 0 )
				break;
			numtrans++;
		}
	}
	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not rewind %s.\n", chTPFilename );
		cdEXIT(EXIT_FAILURE);
	}
	//Skip the magic numbers this time
	read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );


	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the transition probability file.\n");
		cdEXIT(EXIT_FAILURE);
	}

	long ipLo = 0;
	long ipHi = 0;
	double Aul = 0.0;
	lgSentinelReached = false;
	//Read the remaining lines of the transition probability file
	do
	{
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}

		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			i = 1;

			/* skip null lines */
			if( chLine[0] == '\n')
				continue;

			//This means last data column has Aul.
			if( nMatch("A",chLine) )
			{
				/* reset read pointer */
				i = 1;

				ipLo= (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;
				ipHi = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;
				Aul = (double)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
				if( lgEOL )
				{
					fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chTPFilename);
					fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
					cdEXIT(EXIT_FAILURE);
				}

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these lines
					continue;
				}

				TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
				if( (*tr).hasEmis() )
				{
					fprintf(ioQQQ," PROBLEM duplicate transition found by atmdat_STOUT_readin in %s, "
							"wavelength=%f\n", chTPFilename,dBaseStates[intNS][ipHi].energy().WN()
							- dBaseStates[intNS][ipLo].energy().WN());
					cdEXIT(EXIT_FAILURE);
				}

				if( (*tr).EnergyWN() > ENERGY_MIN_WN )
				{
					(*tr).AddLine2Stack();
					(*tr).Emis().Aul() = Aul;
					(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());
				}
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chTPFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chTPFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column of that line.");
		cdEXIT(EXIT_FAILURE);
	}
	fclose(ioDATA);

	if( DEBUGSTATE )
	{
		/*Test whether stout A's are read in properly */
		printf("DEBUG STOUT A's:\n");
		for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
				 tr!= dBaseTrans[intNS].end(); ++tr)
		{
			printf("DEBUG.STOUT.A:\t%i\t%i\t%8.2e\n",(*tr).ipLo()+1,(*tr).ipHi()+1,(*tr).Emis().Aul());
		}
	}



	/******************************************************
	 ************* Collision Data File ********************
	 ******************************************************/

	strcat( chCOLLFilename , ".coll");
	uncaps( chCOLLFilename );

	ioDATA = open_data( chCOLLFilename, "r" );

	/* first line is a version number - now confirm that it is valid */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of %s.\n", chCOLLFilename );
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	iyrread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	imoread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	idyread = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

	if(( iyrread != iyr ) ||
	  (  imoread != imo ) ||
	  (  idyread != idy ) )
	{
		fprintf( ioQQQ,
			" PROBLEM atmdat_STOUT_readin: the version of %s is not the current version.\n", chCOLLFilename );
		fprintf( ioQQQ,
			" atmdat_STOUT_readin: I expected the magic numbers %li %li %li but found %li %li %li.\n",
			iyr, imo , idy ,iyrread, imoread , idyread  );
		cdEXIT(EXIT_FAILURE);
	}

	/****** Could add ability to count number of temperature changes, electron CS, and proton CS ****/


	//Read the first line of data
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_STOUT_readin could not read first line of the collision data file.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* Malloc space for collision strengths */
	StoutCollData[intNS] = (StoutColls***)MALLOC(nMolLevs *sizeof(StoutColls**));
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		StoutCollData[intNS][ipHi] = (StoutColls **)MALLOC((unsigned long)(nMolLevs)*sizeof(StoutColls *));
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			StoutCollData[intNS][ipHi][ipLo] =
				(StoutColls *)MALLOC((unsigned long)(ipNCOLLIDER)*sizeof(StoutColls ));

			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				StoutCollData[intNS][ipHi][ipLo][k].ntemps = -1;
				StoutCollData[intNS][ipHi][ipLo][k].temps = NULL;
				StoutCollData[intNS][ipHi][ipLo][k].collstrs = NULL;
				StoutCollData[intNS][ipHi][ipLo][k].lgIsRate = false;
			}
		}
	}

	ipLo = 0;
	ipHi = 0;
	int numpoints = 0;
	double *temps = NULL;
	long ipCollider = -1;
	lgSentinelReached = false;

	//Read the remaining lines of the collision data file
	do
	{
		/* Stop on *** */
		if( chLine[0] == '*' )
		{
			lgSentinelReached = true;
			break;
		}

		//Comments start with #, skip them
		if( chLine[0] != '#' )
		{
			i = 1;

			/* Skip blank lines */
			if( chLine[0] == '\n')
				continue;

			//This is a temperature line
			if( nMatch("TEMP",chLine) )
			{
				if( temps != NULL)
					free( temps );

				i = 1;
				numpoints = (int)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
				ASSERT( numpoints > 0 );

				temps = (double*)MALLOC(numpoints*sizeof(double));
				memset( temps, 0, (unsigned long)(numpoints)*sizeof(double) );
				for( int j = 0; j < numpoints; j++ )
				{
					temps[j] = (double)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
					ASSERT( temps[j] > 0 );
				}
			}
			else if( nMatch("CS",chLine) || nMatch("RATE",chLine) )
			{

				bool isRate = false;
				if( nMatch("RATE", chLine) )
					isRate = true;

				if( nMatch( "ELECTRON",chLine ) )
				{
					ipCollider = ipELECTRON;
				}
				else if( nMatch( "PROTON",chLine ) )
				{
					ipCollider = ipPROTON;
				}
				else
				{
					fprintf( ioQQQ,	" PROBLEM atmdat_STOUT_readin: Each line of the collision data"
							"file must specify the collider.\n");
					fprintf( ioQQQ,	" Possible colliders are: ELECTRON, PROTON\n");
					cdEXIT(EXIT_FAILURE);
				}

				if( temps == NULL )
				{
					fprintf( ioQQQ,	" PROBLEM atmdat_STOUT_readin: The collision "
							"file must specify temperatures before the collision strengths");
					cdEXIT(EXIT_FAILURE);
				}

				i = 1;
				ipLo= (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;
				ipHi = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL) - 1;

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these lines
					continue;
				}

				/* Set this as a collision strength not a collision rate */
				StoutCollData[intNS][ipHi][ipLo][ipCollider].lgIsRate = isRate;

				StoutCollData[intNS][ipHi][ipLo][ipCollider].ntemps = numpoints;
				ASSERT( numpoints > 0 );

				/*malloc the space here for the temps and collision strengths*/
				StoutCollData[intNS][ipHi][ipLo][ipCollider].temps =
					(double *)MALLOC((unsigned long)(numpoints)*sizeof(double));
				StoutCollData[intNS][ipHi][ipLo][ipCollider].collstrs =
					(double *)MALLOC((unsigned long)(numpoints)*sizeof(double));
				/* Loop over all but one CS value */
				for( int j = 0; j < numpoints; j++ )
				{
					StoutCollData[intNS][ipHi][ipLo][ipCollider].temps[j] = temps[j];
					ASSERT( StoutCollData[intNS][ipHi][ipLo][ipCollider].temps[j] > 0 );
					StoutCollData[intNS][ipHi][ipLo][ipCollider].collstrs[j] = (double)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
					ASSERT( StoutCollData[intNS][ipHi][ipLo][ipCollider].collstrs[j] > 0 );
				}
			}
			else
			{
				fprintf( ioQQQ, " PROBLEM File %s contains an invalid line.\n",chCOLLFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}

			if( lgEOL )
			{
				fprintf( ioQQQ, " PROBLEM End of line reached prematurely in file %s\n",chCOLLFilename);
				fprintf( ioQQQ, " The line being read is between the braces {%.*s}\n",int(strlen(chLine)-1),chLine);
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL );
	if( !lgSentinelReached )
	{
		fprintf( ioQQQ, " PROBLEM End of data sentinel was not reached in file %s\n",chCOLLFilename);
		fprintf( ioQQQ, " Check that stars (*****) appear after the last line of data and start in the first column.");
		cdEXIT(EXIT_FAILURE);
	}
	free(temps);
	fclose(ioDATA);

	return;
}

void atmdat_CHIANTI_readin( long intNS, char *chPrefix )
{
	DEBUG_ENTRY( "atmdat_CHIANTI_readin()" );

	int intsplinepts,intTranType,intxs;
	long int nMolLevs,nMolExpLevs,nElvlcLines;// number of experimental and total levels
	FILE *ioElecCollData=NULL, *ioProtCollData=NULL;
	realnum  fstatwt,fenergyWN,fWLAng,fenergy,feinsteina;
	double fScalingParam,fEnergyDiff,*xs,*spl,*spl2;

	char chLine[FILENAME_PATH_LENGTH_2] ,
		chEnFilename[FILENAME_PATH_LENGTH_2],
		chTraFilename[FILENAME_PATH_LENGTH_2],
		chEleColFilename[FILENAME_PATH_LENGTH_2],
		chProColFilename[FILENAME_PATH_LENGTH_2];

	realnum *dBaseStatesEnergy;
	long int *intNewIndex=NULL,*intOldIndex=NULL, *intExperIndex=NULL;
	int interror;
	bool lgProtonData=false,lgEneLevOrder;

	// this is the largest number of levels allowed by the chianti format, I3
	static const int MAX_NUM_LEVELS = 999;

	dBaseSpecies[intNS].lgMolecular = false;
	dBaseSpecies[intNS].lgLAMDA = false;
	dBaseSpecies[intNS].lgLTE = false;

	strcpy( chEnFilename , chPrefix );
	strcpy( chTraFilename , chPrefix );	
	strcpy( chEleColFilename , chPrefix );		
	strcpy( chProColFilename , chPrefix );			

	/*For the CHIANTI DATABASE*/
	/*Open the energy levels file*/
	strcat( chEnFilename , ".elvlc");
	uncaps( chEnFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chEnFilename);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEnFilename);

	fstream elvlcstream;
	open_data(elvlcstream, chEnFilename,mode_r);

	/*Open the transition probabilities file*/
	strcat( chTraFilename , ".wgfa");
	uncaps( chTraFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chTraFilename);

	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chTraFilename);

	fstream wgfastream;
	open_data(wgfastream, chTraFilename,mode_r);

	/*Open the electron collision data*/
	strcat( chEleColFilename , ".splups");
	uncaps( chEleColFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chEleColFilename);
	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chEleColFilename);

	ioElecCollData = open_data( chEleColFilename, "r" );

	/*Open the proton collision data*/
	strcat( chProColFilename , ".psplups");
	uncaps( chProColFilename );
	if(DEBUGSTATE)
		printf("The data file is %s \n",chProColFilename);
	/*Open the files*/
	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_CHIANTI_readin opening %s:",chProColFilename);

	/*We will set a flag here to indicate if the proton collision strengths are available */
	if( ( ioProtCollData = fopen( chProColFilename , "r" ) ) != NULL )
	{
		lgProtonData = true;
		//fclose( ioProtCollData );
		//ioProtCollData = NULL;
	}
	else
	{
		lgProtonData = false;
	}

	/*Loop finds how many theoretical and experimental levels are in the elvlc file */
	//eof_col is used get the first 4 charcters per line to find end of file
	const int eof_col = 5;
	//length (+1) of the nrg in the elvlc file
	const int lvl_nrg_col=16;
	//# of columns skipped from the left to get to nrg start
	const int lvl_skipto_nrg = 40;
	/* # of columns to skip from eof check to nrg start */
	const int lvl_eof_to_nrg = lvl_skipto_nrg - eof_col + 1;
	nElvlcLines = 0;
	nMolExpLevs = 1;
	lgEneLevOrder = true;
	if (elvlcstream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		char exptemp[lvl_nrg_col];
		double tempexpenergy = 0.;
		/*This loop counts the number of valid rows within the elvlc file
		  as well as the number of experimental energy levels.*/
		while(nj != -1)
		{
			elvlcstream.get(otemp,eof_col);
			nj = atoi(otemp);
			if( nj == -1)
				break;
			nElvlcLines++;

			elvlcstream.seekg(lvl_eof_to_nrg,ios::cur);
			elvlcstream.get(exptemp,lvl_nrg_col);
			tempexpenergy = (realnum) atof(exptemp);
			if( tempexpenergy != 0.)
				nMolExpLevs++;

			elvlcstream.ignore(INT_MAX,'\n');

		}
		elvlcstream.seekg(0,ios::beg);
	}

	if(DEBUGSTATE)
	{
		printf("DEBUG: The file %s contains %li experimental energy levels of the %li total levels. \n",chEnFilename,nMolExpLevs,nElvlcLines);
	}

	/* The total number of levels depends on the experimental Chianti switch */
	if( atmdat.lgChiantiExp )
	{
		if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
		{
			// Fe is special case with more levels
			nMolLevs = MIN3(nMolExpLevs, atmdat.nChiantiMaxLevelsFe, MAX_NUM_LEVELS );
		}
		else
		{
			nMolLevs = MIN3(nMolExpLevs, atmdat.nChiantiMaxLevels, MAX_NUM_LEVELS );
		}
	}
	else
	{
		if( tolower(dBaseSpecies[intNS].chLabel[0]) == 'f' && tolower(dBaseSpecies[intNS].chLabel[1]) == 'e')
		{
			// Fe is special case with more levels
			nMolLevs = MIN3(nElvlcLines, atmdat.nChiantiMaxLevelsFe,MAX_NUM_LEVELS );
		}
		else
		{
			nMolLevs = MIN3(nElvlcLines, atmdat.nChiantiMaxLevels,MAX_NUM_LEVELS );
		}
	}

	if( nMolLevs <= 0 )
	{
		fprintf( ioQQQ, "The number of energy levels is non-positive in datafile %s.\n", chEnFilename );
		fprintf( ioQQQ, "The file must be corrupted.\n" );
		fclose( ioProtCollData );
		cdEXIT( EXIT_FAILURE );
	}

	dBaseSpecies[intNS].numLevels_max = nMolLevs;
	dBaseSpecies[intNS].numLevels_local = dBaseSpecies[intNS].numLevels_max;

	if( atmdat.lgChiantiPrint == true)
	{
		if( atmdat.lgChiantiExp )
		{
			fprintf( ioQQQ,"     Using CHIANTI model %s with %li experimental energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel , nMolLevs , nMolExpLevs );
		}
		else
		{
			fprintf( ioQQQ,"     Using CHIANTI model %s with %li theoretical energy levels of %li available.\n",
				dBaseSpecies[intNS].chLabel , nMolLevs , nElvlcLines );
		}
	}

	/*Malloc the energy array to check if the energy levels are in order*/
	dBaseStatesEnergy = (realnum*)MALLOC((unsigned long)(nMolLevs)*sizeof(realnum));

	/*malloc the States array*/
	dBaseStates[intNS].resize(nMolLevs);

	/* allocate the Transition array*/
	ipdBaseTrans[intNS].reserve(nMolLevs);
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
		ipdBaseTrans[intNS].reserve(ipHi,ipHi);
	ipdBaseTrans[intNS].alloc();
	dBaseTrans[intNS].resize(ipdBaseTrans[intNS].size());
	dBaseTrans[intNS].states() = &dBaseStates[intNS];
	dBaseTrans[intNS].chLabel() = "Chianti";

	int ndBase = 0;
	for( long ipHi = 1; ipHi < nMolLevs; ipHi++)
	{
		for( long ipLo = 0; ipLo < ipHi; ipLo++)
		{
			ipdBaseTrans[intNS][ipHi][ipLo] = ndBase;
			dBaseTrans[intNS][ndBase].Junk();
			dBaseTrans[intNS][ndBase].setLo(ipLo);
			dBaseTrans[intNS][ndBase].setHi(ipHi);
			++ndBase;
		}
	}

	/*Reading in the energy and checking that they are in order*/
	// C++ io works in terms of cursor movement rather than position on line
	//# of columns to skip over the rydberg energy, we don't use it
	const int lvl_skip_ryd = 15;

	/*Keep track of which levels have experimental data and then create a vector
	which relates their indices to the default chianti energy indices.
	 */
	long ncounter = 0;

	if( atmdat.lgChiantiExp )
	{
		//Relate Chianti level indices to a set that only include experimental levels
		intExperIndex = (long int*)MALLOC((unsigned long)(nElvlcLines)* sizeof(long int));

		//Fill relational vector with bad values
		for(int i = 0; i < nElvlcLines; i++)
		{
			intExperIndex[i] = -1;
		}
	}

	//Read in nrg levels to see if they are in order
	for( long ipLev=0; ipLev<nElvlcLines; ipLev++ )
	{
		if(elvlcstream.is_open())
		{
			char thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
			elvlcstream.seekg(lvl_skipto_nrg,ios::cur);
			elvlcstream.get(thtemp,lvl_nrg_col);
			fenergy = (realnum) atof(thtemp);

			if( atmdat.lgChiantiExp )
			{
				/* Go through the entire level list selectively choosing only experimental level energies.
				 * Store them, not zeroes, in order using ncounter to count the index.
				 * Any row on the level list where there is no experimental energy, put a -1 in the relational vector.
				 * If it is a valid experimental energy level store the new ncounter index.
				 */

				//Once we collect enough experimental levels stop looking for more.
				if( ncounter == nMolLevs)
					break;
				ASSERT( ncounter < nMolLevs );

				if(fenergy != 0. || ipLev == 0 )
				{
					dBaseStatesEnergy[ncounter] = fenergy;
					intExperIndex[ipLev] = ncounter;
					ncounter++;
				}
				else
				{
					intExperIndex[ipLev] = -1;
				}

			}
			else
			{
				//Stop looking for levels when the array is full
				if( ipLev == nMolLevs )
					break;

				elvlcstream.seekg(lvl_skip_ryd,ios::cur);
				elvlcstream.get(obtemp,lvl_nrg_col);
				if(atof(obtemp) != 0.)
				{
					fenergy = (realnum) atof(obtemp);
				}
				//else
					//fprintf( ioQQQ," WARNING: Switched to theoretical energies for species %s, level %3li\n", dBaseSpecies[intNS].chLabel, ipLev );
				dBaseStatesEnergy[ipLev] = fenergy;
			}

			elvlcstream.ignore(INT_MAX,'\n');

		}
		else
		{
			fprintf( ioQQQ, " The data file %s is corrupted .\n",chEnFilename);
			fclose( ioProtCollData );
			cdEXIT(EXIT_FAILURE);
		}
	}

	for( long ipLev=1; ipLev<nMolLevs; ipLev++ )
	{
		/*To check if energy levels are in order*/
		/*If the energy levels are not in order a flag is set*/
		if (dBaseStatesEnergy[ipLev] < dBaseStatesEnergy[ipLev-1])
		{
			lgEneLevOrder = false;
		}
	}

	// malloc vector for new energy-sorted indices
	intNewIndex = (long int*)MALLOC((unsigned long)(nElvlcLines)* sizeof(long int));

	if(lgEneLevOrder == false)
	{
		/*If the energy levels are not in order and the flag is set(lgEneLevOrder=FALSE)
		then we sort the energy levels to find the relation matrix between the old and new indices*/
		intOldIndex = (long int*)MALLOC((unsigned long)(nMolLevs)* sizeof(long int));
		/*We now do a modified quick sort*/
		spsort(dBaseStatesEnergy,nMolLevs,intOldIndex,2,&interror);
		/*intNewIndex has the key to map*/
		for( long i=0; i<nMolLevs; i++ )
		{
			ASSERT( intOldIndex[i] < nMolLevs );
			intNewIndex[intOldIndex[i]] = i;	
		}

		if( nMolLevs < nElvlcLines )
		{
			/* If any chianti energy levels do not have experimental values
			 * the vector that relates experimental levels to the default will not be the same size
			 * as the stored experimental energy levels.
			 * Therefore this code is required to pad the rest of the sorted vector.
			 * Any index larger than nMolLevs set to -1 */
			for( long i=nMolLevs; i<nElvlcLines; i++)
			{
				intNewIndex[i] = -1;
			}
		}

		if(DEBUGSTATE)
		{
			for( long i=0; i<nMolLevs; i++ )
			{
				printf("The %ld value of the relation matrix is %ld \n",i,intNewIndex[i]);
			}
			for( long i=0; i<nMolLevs; i++ )
			{
				printf("The %ld value of the sorted energy array is %f \n",i,dBaseStatesEnergy[i]);
			}
		}
		free( intOldIndex );
	}
	else
	{
		/* if energies already in correct order, new index is same as original. */
		for( long i=0; i<nMolLevs; i++ )
		{
			intNewIndex[i] = i;	
		}

		if( nMolLevs < nElvlcLines )
		{
			/*Pad the sorted vector so that it has the correct number of elements
			 * Any index larger than nMolLevs set to -1.*/
			for( long i=nMolLevs; i<nElvlcLines; i++)
			{
				intNewIndex[i] = -1;
			}
		}
	}

	/* insist that there the intNewIndex values are unique */
	for( long i=0; i<nMolLevs; i++ )
	{
		for( long j=i+1; j<nMolLevs; j++ )
		{
			ASSERT( intNewIndex[i] != intNewIndex[j] );
		}
	}

	/*This will print out the relational matrix
	 * which includes the original index,the new experimental index,
	 * and the sorted experimental index.
	 * The numbers NOT be on the C scale, so 1 is the lowest level. */
	if( DEBUGSTATE && atmdat.lgChiantiExp )
	{
		printf("\n\n%s Energy Level Matrix\n",dBaseSpecies[intNS].chLabel);
		printf("(Original, Experimental, Sorted)\n");
		for( long ipLevOld=0; ipLevOld<nElvlcLines; ipLevOld++ )
		{
			if( intExperIndex[ipLevOld] != -1)
				printf("%li\t%li\t%li\n",ipLevOld+1,intExperIndex[ipLevOld]+1,intNewIndex[intExperIndex[ipLevOld]]+1);
		}
		printf("\n");
	}

	//After we read in the energies we rewind the energy levels file
	elvlcstream.seekg(0,ios::beg);

	//lvl_skipto_statwt is the # of columns to skip to statwt from left
	const int lvl_skipto_statwt = 37;
	//lvl_statwt_col is the length (+1) of statwt
	const int lvl_statwt_col = 4;
	//Read in stat weight and energy
	for( long ipLevOld=0; ipLevOld<nElvlcLines; ipLevOld++ )
	{
		long ipLevNew = 0;
		if( atmdat.lgChiantiExp )
		{
			/* We know that non-experimental levels are stored as -1
			 * in the observed/experimental vector.
			 * Use that to skip over those values. */

			if( intExperIndex[ipLevOld] == -1 )
			{
				elvlcstream.ignore(INT_MAX,'\n');
				continue;
			}
			else
			{
				/* Store values based on the sorted experimental indices */
				ipLevNew = intNewIndex[intExperIndex[ipLevOld]];
			}
		}
		else
		{
			/* With level trimming on it is possible that there can be rows that
			 * have to be skipped when using theoretical
			 * since the levels no longer exist */
			if( intNewIndex[ipLevOld] == -1 )
			{
				elvlcstream.ignore(INT_MAX,'\n');
				continue;
			}
			else
			{
				ipLevNew = intNewIndex[ipLevOld];
			}
		}

		char gtemp[lvl_statwt_col],thtemp[lvl_nrg_col],obtemp[lvl_nrg_col];
	
		/*information needed for label*/
		strcpy(dBaseStates[intNS][ipLevNew].chLabel(), "    ");
		strncpy(dBaseStates[intNS][ipLevNew].chLabel(),dBaseSpecies[intNS].chLabel, 4);
		dBaseStates[intNS][ipLevNew].chLabel()[4] = '\0';
		
		//Move cursor to and read statwt
		elvlcstream.seekg(lvl_skipto_statwt,ios::cur);
		elvlcstream.get(gtemp,lvl_statwt_col);
		fstatwt = (realnum)atof(gtemp);

		if(fstatwt <= 0.)
		{
			fprintf( ioQQQ, " WARNING: A positive non zero value is expected for the "
					 "statistical weight but was not found in %s"
					" level %li\n", chEnFilename,ipLevNew);
			fstatwt = 1.f;
			//cdEXIT(EXIT_FAILURE);
		}
		dBaseStates[intNS][ipLevNew].g() = fstatwt;

		//Read nrg
		elvlcstream.get(thtemp,lvl_nrg_col);

		fenergy = (realnum) atof(thtemp);

		/* If we are looking for theoretical energies
		 * move over a couple columns before reading in the energy*/
		if( !atmdat.lgChiantiExp )
		{
			elvlcstream.seekg(lvl_skip_ryd,ios::cur);
			elvlcstream.get(obtemp,lvl_nrg_col);
			if(atof(obtemp) != 0.)
			{
				fenergy = (realnum) atof(obtemp);
			}
		}		
		elvlcstream.ignore(INT_MAX,'\n');
		dBaseStates[intNS][ipLevNew].energy().set(fenergy,"cm^-1");
	}
	//Close the level stream
	elvlcstream.close();

	// highest energy transition in chianti
	dBaseSpecies[intNS].maxWN = 0.;
	/* fill in all transition energies, can later overwrite for specific radiative transitions */
	for(TransitionList::iterator tr=dBaseTrans[intNS].begin();
		 tr!= dBaseTrans[intNS].end(); ++tr)
	{
		int ipHi = (*tr).ipHi();
		int ipLo = (*tr).ipLo();
		fenergyWN = (realnum)MAX2( ENERGY_MIN_WN , dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN() );

		(*tr).EnergyWN() = fenergyWN;
		(*tr).WLAng() = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
		dBaseSpecies[intNS].maxWN = MAX2(dBaseSpecies[intNS].maxWN,fenergyWN);
	}

	if(DEBUGSTATE)
	{
		/* Print out the stored data for each level.
		 * Level index is NOT on C scale. */
		printf("\n%s Stored Energy Levels\n",dBaseSpecies[intNS].chLabel);
		printf("(Index,Energy in WN, Statistical Weight)\n");
		for( long ipLo = 0; ipLo < nMolLevs; ipLo++ )
		{
			printf("%li\t%e\t%.1f\n",ipLo+1,dBaseStates[intNS][ipLo].energy().WN(),dBaseStates[intNS][ipLo].g());
		}
	}

	/************************************************************************/
	/*** Read in the level numbers, Einstein A and transition wavelength  ***/
	/************************************************************************/

	//Count the number of rows first
	long wgfarows = -1;
	if (wgfastream.is_open())
	{
		int nj = 0;
		char otemp[eof_col];
		while(nj != -1)
		{
			wgfastream.get(otemp,eof_col);
			wgfastream.ignore(INT_MAX,'\n');
			nj = atoi(otemp);
			wgfarows++;
		}
		wgfastream.seekg(0,ios::beg);
	}
	else 
		fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);

	//line_index_col is the length(+1) of the level indexes in the WGFA file
	const int line_index_col = 6;
	//line_nrg_to_eina is the # of columns to skip from wavelength to eina in WGFA file
	const int line_nrg_to_eina =15;
	//line_eina_col is the length(+1) of einsteinA in WGFA
	const int line_eina_col = 16;
	char lvltemp[line_index_col];
	//Start reading WGFA file
	if (wgfastream.is_open())
	{
		for (long i = 0;i<wgfarows;i++)
		{
			wgfastream.get(lvltemp,line_index_col);
			long ipLo = atoi(lvltemp)-1;
			wgfastream.get(lvltemp,line_index_col);
			long ipHi = atoi(lvltemp)-1;

			if( atmdat.lgChiantiExp )
			{
				/* If either upper or lower index is -1 in the relational vector,
				 * skip that line in the wgfa file.
				 * Otherwise translate the level indices.*/
				if( intExperIndex[ipLo] == - 1 || intExperIndex[ipHi] == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipLo = intNewIndex[intExperIndex[ipLo]];
					ipHi = intNewIndex[intExperIndex[ipHi]];
				}
			}
			else
			{
				/* With level trimming on it is possible that there can be rows that
				 * have to be skipped when using theoretical
				* since the levels no longer exist */
				if( intNewIndex[ipLo] == - 1 || intNewIndex[ipHi] == -1 )
				{
					wgfastream.ignore(INT_MAX,'\n');
					continue;
				}
				else
				{
					ipLo = intNewIndex[ipLo];
					ipHi = intNewIndex[ipHi];
				}
			}

			if( ipLo >= nMolLevs || ipHi >= nMolLevs )
			{
				// skip these lines
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}
	
			if( ipHi == ipLo )
			{
				fprintf( ioQQQ," WARNING: Upper level = lower for a radiative transition in %s. Ignoring.\n", chTraFilename );
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

	
			ASSERT( ipHi != ipLo );
			ASSERT(ipHi >= 0);
			ASSERT(ipLo >= 0);

			// sometimes the CHIANTI datafiles list the highest index first as in the middle of these five lines in ne_10.wgfa:
			//    ...
			//    8   10       187.5299      0.000e+00      4.127e+05                 3d 2D1.5 -                  4s 2S0.5           E2
			//    9   10       187.6573      0.000e+00      6.197e+05                 3d 2D2.5 -                  4s 2S0.5           E2
			//   11   10   4842624.0000      1.499e-05      9.423e-06                 4p 2P0.5 -                  4s 2S0.5           E1
			//    1   11         9.7085      1.892e-02      6.695e+11                 1s 2S0.5 -                  4p 2P0.5           E1
			//    2   11        48.5157      6.787e-02      9.618e+10                 2s 2S0.5 -                  4p 2P0.5           E1
			//    ...
			// so, just set ipHi (ipLo) equal to the max (min) of the two indices.
			// NB NB NB it looks like this may depend upon whether one uses observed or theoretical energies.

			//Read in wavelengths
			char trantemp[lvl_nrg_col];
			wgfastream.get(trantemp,lvl_nrg_col);
			fWLAng = (realnum)atof(trantemp);

			/* \todo 2 CHIANTI labels the H 1 2-photon transition as z wavelength of zero.
			 * Should we just ignore all of the wavelengths in this file and use the
			 * difference of level energies instead. */

			if( atmdat.lgChiantiExp )
			{
				/* Sometimes the wgfa file lists the levels to where ipLo > ipHi.
				 * This seems to be related to the theoretical energies being out of order for those transitions.
				 * When this is true it seems that the associated wavelength is negative which means they are using theoretical energies,
				 * even though they have experimental energies.
				 * For now, print that the levels are switched and ignore the lines if the wavelength is negative. */

				if( ipHi < ipLo )
				{
					if( strcmp(dBaseSpecies[intNS].chLabel,"Fe 3") == 0)
					{
						long ipLoTemp = ipLo;
						long ipHiTemp = ipHi;
						ipHi = max( ipHiTemp, ipLoTemp );
						ipLo = min( ipHiTemp, ipLoTemp );
						if( atmdat.lgChiantiPrint )
						{
							fprintf( ioQQQ," WARNING: Swapped level indices for species %6s, indices %3li %3li with wavelength %e \n",
								dBaseSpecies[intNS].chLabel, ipLoTemp, ipHiTemp,fWLAng );
						}
					}
					else if( atmdat.lgChiantiPrint )
					{
						fprintf( ioQQQ," WARNING: Upper and Lower Indices are reversed.Species: %6s, Indices: %3li %3li Wavelength: %e \n",
							dBaseSpecies[intNS].chLabel, ipLo+1, ipHi+1,fWLAng );
					}
				}

				if( fWLAng < 0.)
				{
					// skip these lines
					wgfastream.ignore(INT_MAX,'\n');
					if( atmdat.lgChiantiPrint )
					{
						fprintf(ioQQQ,"WARNING: Skipping Transition %li to %li of %s.\n",ipLo+1,ipHi+1,dBaseSpecies[intNS].chLabel);
					}
					continue;
				}

			}
			else
			{
				/* If Cloudy is using theoretical energies, then just flip the levels. */
				if( ipHi < ipLo )
				{
					long ipLoTemp = ipLo;
					long ipHiTemp = ipHi;
					ipHi = max( ipHiTemp, ipLoTemp );
					ipLo = min( ipHiTemp, ipLoTemp );
					fprintf( ioQQQ," WARNING: Swapped level indices for species %6s, indices %3li %3li with wavelength %e \n",
							dBaseSpecies[intNS].chLabel, ipLoTemp, ipHiTemp,fWLAng );
				}
			}

			/* If the given wavelength is negative, then theoretical enegies are being used.
			 * Take the difference in stored theoretical energies.
			 * It should equal the absolute value of the wavelength in the wgfa file. */
			if( !atmdat.lgChiantiExp && fWLAng <= 0. )
			{
				//if( fWLAng < 0.)
					//fprintf( ioQQQ," WARNING: Negative wavelength for species %6s, indices %3li %3li \n", dBaseSpecies[intNS].chLabel, ipLo, ipHi);
				fWLAng = (realnum)(1e8/
					(dBaseStates[intNS][ipHi].energy().WN() - dBaseStates[intNS][ipLo].energy().WN()));
			}
			//Skip from end of Wavelength to Einstein A and read in
			wgfastream.seekg(line_nrg_to_eina,ios::cur);
			wgfastream.get(trantemp,line_eina_col);
			feinsteina = (realnum)atof(trantemp);
			if( feinsteina < 1e-20 )
			{
				static bool notPrintedYet = true;
				if( notPrintedYet && atmdat.lgChiantiPrint)
				{
					fprintf( ioQQQ," WARNING: Radiative rate(s) equal to zero in %s.\n", chTraFilename );
					notPrintedYet = false;
				}
				wgfastream.ignore(INT_MAX,'\n');
				continue;
			}

			fixit(); // may need to do something with these rates
			//Read in the rest of the line and look for auto
			wgfastream.getline(chLine,INT_MAX);
			TransitionList::iterator tr = dBaseTrans[intNS].begin()+ipdBaseTrans[intNS][ipHi][ipLo];
			if( nMatch("auto", chLine) )
			{
				if( (*tr).hasEmis() )
				{
					(*tr).Emis().AutoIonizFrac() =
						feinsteina/((*tr).Emis().Aul() + feinsteina);
					ASSERT( (*tr).Emis().AutoIonizFrac() >= 0. );
					ASSERT( (*tr).Emis().AutoIonizFrac() <= 1. );
				}
				continue;
			}

			if( (*tr).hasEmis() )
			{
				fprintf(ioQQQ," PROBLEM duplicate transition found by atmdat_chianti in %s, "
						"wavelength=%f\n", chTraFilename,fWLAng);
				fclose( ioProtCollData );
				cdEXIT(EXIT_FAILURE);
			}
			(*tr).AddLine2Stack();
			(*tr).Emis().Aul() = feinsteina;
			(*tr).WLAng() = fWLAng;
		 
			fenergyWN = (realnum)(1e+8/fWLAng);

			// TODO::Check the wavelength in the file with the difference in energy levels

			(*tr).EnergyWN() = fenergyWN;
			(*tr).WLAng() = (realnum)(1e+8/fenergyWN/RefIndex(fenergyWN));
			(*tr).Emis().gf() = (realnum)GetGF((*tr).Emis().Aul(), (*tr).EnergyWN(), (*(*tr).Hi()).g());

			if(DEBUGSTATE)
			{
				fprintf( ioQQQ, "The lower level is %ld \n",ipLo);
				fprintf( ioQQQ, "The upper level is %ld \n",ipHi);
				fprintf( ioQQQ, "The Einstein A is %f \n",(*tr).Emis().Aul());
				fprintf( ioQQQ, "The wavelength of the transition is %f \n",(*tr).WLAng() );
			}
		}
	}
	else fprintf( ioQQQ, " The data file %s is corrupted .\n",chTraFilename);
	wgfastream.close();

	/* Malloc space for splines */
	AtmolCollSplines[intNS] = (CollSplinesArray***)MALLOC(nMolLevs *sizeof(CollSplinesArray**));
	for( long ipHi=0; ipHi<nMolLevs; ipHi++ )
	{
		AtmolCollSplines[intNS][ipHi] = (CollSplinesArray **)MALLOC((unsigned long)(nMolLevs)*sizeof(CollSplinesArray *));
		for( long ipLo=0; ipLo<nMolLevs; ipLo++ )
		{
			AtmolCollSplines[intNS][ipHi][ipLo] = 
				(CollSplinesArray *)MALLOC((unsigned long)(ipNCOLLIDER)*sizeof(CollSplinesArray ));

			for( long k=0; k<ipNCOLLIDER; k++ )
			{
				/* initialize all spline variables */
				AtmolCollSplines[intNS][ipHi][ipLo][k].collspline = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].SplineSecDer = NULL;
				AtmolCollSplines[intNS][ipHi][ipLo][k].nSplinePts = -1; 
				AtmolCollSplines[intNS][ipHi][ipLo][k].intTranType = -1;
				AtmolCollSplines[intNS][ipHi][ipLo][k].EnergyDiff = BIGDOUBLE;
				AtmolCollSplines[intNS][ipHi][ipLo][k].ScalingParam = BIGDOUBLE;
			}
		}
	}

	/************************************/
	/*** Read in the collisional data ***/
	/************************************/

	// ipCollider 0 is electrons, 1 is protons
	for( long ipCollider=0; ipCollider<=1; ipCollider++ )
	{
		char chFilename[FILENAME_PATH_LENGTH_2];

		if( ipCollider == ipELECTRON )
		{
			strcpy( chFilename, chEleColFilename );
		}
		else if( ipCollider == ipPROTON )
		{
			if( !lgProtonData )
				break;
			fprintf( ioQQQ," WARNING: Chianti proton collision data have different format than electron data files!\n" );
			strcpy( chFilename, chProColFilename );
		}
		else
			TotalInsanity();

		/*Dummy string used for convenience*/
		strcpy(chLine,"A");

		fstream splupsstream;
		open_data(splupsstream, chFilename,mode_r);

		//cs_eof_col is the length(+1) of the first column used for finding the end of file
		const int cs_eof_col = 4;
		//cs_index_col is the length(+1) of the indexes in the CS file
		const int cs_index_col = 4;
		//cs_trantype_col is the length(+1) of the transition type in the CS file
		const int cs_trantype_col = 4;
		//cs_values_col is the length(+1) of the other values in the CS file
		//including: GF, nrg diff, scaling parameter, and spline points
		const int cs_values_col = 11;
		//Determine the number of rows in the CS file
		if (splupsstream.is_open())
		{
			int nj = 0;
			//splupslines is -1 since the loop runs 1 extra time
			long splupslines = -1;
			char otemp[cs_eof_col];
			while(nj != -1)
			{
				splupsstream.get(otemp,cs_eof_col);
				splupsstream.ignore(INT_MAX,'\n');
				nj = atoi(otemp);
				splupslines++;
			}
			splupsstream.seekg(0,ios::beg);
	
			for (int m = 0;m<splupslines;m++)
			{
				if( ipCollider == ipELECTRON )
				{
					splupsstream.seekg(6,ios::cur);
				}

				/* level indices */
				splupsstream.get(otemp,cs_index_col);
				long ipLo = atoi(otemp)-1;
				splupsstream.get(otemp,cs_index_col);
				long ipHi = atoi(otemp)-1;

				/* If either upper or lower index is -1 in the relational vector,
				* skip that line in the splups file.
				* Otherwise translate the level indices.*/
				if( atmdat.lgChiantiExp )
				{
					if( intExperIndex[ipLo] == - 1 || intExperIndex[ipHi] == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipLo = intNewIndex[intExperIndex[ipLo]];
						ipHi = intNewIndex[intExperIndex[ipHi]];
					}

				}
				else
				{
					/* With level trimming on it is possible that there can be rows that
					 * have to be skipped when using theoretical
					 * since the levels no longer exist */
					if( intNewIndex[ipLo] == - 1 || intNewIndex[ipHi] == -1 )
					{
						splupsstream.ignore(INT_MAX,'\n');
						continue;
					}
					else
					{
						ipLo = intNewIndex[ipLo];
						ipHi = intNewIndex[ipHi];
					}
				}

				if( ipLo >= nMolLevs || ipHi >= nMolLevs )
				{
					// skip these transitions
					splupsstream.ignore(INT_MAX,'\n');
					continue;
				}

				ASSERT( ipLo < nMolLevs );
				ASSERT( ipHi < nMolLevs );
				/*Transition Type*/
				splupsstream.get(otemp,cs_trantype_col);
				intTranType = atoi(otemp);
				char qtemp[cs_values_col];
				splupsstream.get(qtemp,cs_values_col);
				/*Energy Difference*/
				splupsstream.get(qtemp,cs_values_col);
				fEnergyDiff = atof(qtemp);
				/*Scaling Parameter*/
				splupsstream.get(qtemp,cs_values_col);
				fScalingParam = atof(qtemp);

				ASSERT( ipLo < nMolLevs );
				ASSERT( ipHi < nMolLevs );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline == NULL );
				ASSERT( AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer == NULL );

				const int CHIANTI_SPLINE_MAX=9, CHIANTI_SPLINE_MIN=5;
				STATIC_ASSERT(CHIANTI_SPLINE_MAX > CHIANTI_SPLINE_MIN);

				/*We malloc the space here*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer =
					(double *)MALLOC((unsigned long)(CHIANTI_SPLINE_MAX)*sizeof(double));

				/* always read at least CHIANTI_SPLINE_MIN */
				for( intsplinepts=0; intsplinepts<=CHIANTI_SPLINE_MAX; intsplinepts++ )
				{
					//Look at the next character to see if it is the end of line.
					char p = splupsstream.peek();
					if( p == '\n' )
					{
						/* there should be 5 or 9 spline points.  If we got EOL,
						 * insist that we were trying to read the 6th or 10th. */
						if( (intsplinepts != CHIANTI_SPLINE_MAX) && (intsplinepts != CHIANTI_SPLINE_MIN) )
						{
							static bool notPrintedYet = true;
							if( notPrintedYet )
							{
								fprintf( ioQQQ, " WARNING: Wrong number of spline points in %s.\n", chFilename);
								notPrintedYet = false;
							}
							for( long j=intsplinepts-1; j<CHIANTI_SPLINE_MAX; j++ )
								AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[j] = 0.;
						}
						else
							break;
					}
					else
					{
						if( intsplinepts >= CHIANTI_SPLINE_MAX )
						{
							fprintf( ioQQQ, " WARNING: More spline points than expected in %s, indices %3li %3li.  Ignoring extras.\n", chFilename, ipHi, ipLo );
							break;
						}
						ASSERT( intsplinepts < CHIANTI_SPLINE_MAX );
						double temp;
						//Store a single spline point then look for more
						splupsstream.get(qtemp,cs_values_col);
						temp = atof(qtemp);
						if(temp < 0)
						{
							temp = 0.;
						}
						AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intsplinepts] = temp;
					}

				}

				ASSERT( (intsplinepts == CHIANTI_SPLINE_MAX) || (intsplinepts == CHIANTI_SPLINE_MIN) );

				/*The zeroth element contains the number of spline points*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].nSplinePts = intsplinepts;
				/*Transition type*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].intTranType = intTranType;
				/*Energy difference between two levels in Rydbergs*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].EnergyDiff = fEnergyDiff;
				/*Scaling parameter C*/
				AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].ScalingParam = fScalingParam;

				/*Once the spline points have been filled,fill the second derivatives structure*/
				/*Creating spline points array*/
				xs = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));
				spl = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));
				spl2 = (double *)MALLOC((unsigned long)intsplinepts*sizeof(double));

				// should be able to just loop to intsplinepts -- ASSERT above protects
				if(intsplinepts == CHIANTI_SPLINE_MIN)
				{
					for(intxs=0;intxs<CHIANTI_SPLINE_MIN;intxs++)
					{
						xs[intxs] = 0.25*intxs;
						spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
					}
				}
				else if(intsplinepts == CHIANTI_SPLINE_MAX)
				{
					for(intxs=0;intxs<CHIANTI_SPLINE_MAX;intxs++)
					{
						xs[intxs] = 0.125*intxs;
						spl[intxs] = AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].collspline[intxs];
					}
				}
				else
					TotalInsanity();

				spline(xs, spl,intsplinepts,2e31,2e31,spl2);

				/*Filling the second derivative structure*/
				for( long i=0; i<intsplinepts; i++)
				{
					AtmolCollSplines[intNS][ipHi][ipLo][ipCollider].SplineSecDer[i] = spl2[i];
				}

				free( xs );
				free( spl );
				free( spl2 );
				splupsstream.ignore(INT_MAX,'\n');
			}
			splupsstream.close();
		}
	}

	// free some memory
	free( dBaseStatesEnergy );
	free( intNewIndex );
	if( atmdat.lgChiantiExp)
	{
		free( intExperIndex );
	}

	// close open file handles
	fclose( ioElecCollData );
	if( lgProtonData )
		fclose( ioProtCollData );

	return;
}
