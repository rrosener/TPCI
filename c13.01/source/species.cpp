/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "lines_service.h"
#include "taulines.h"
#include "trace.h"
#include "string.h"
#include "input.h"
#include "thirdparty.h"
#include "dense.h"
#include "atmdat.h"
#include "mole.h"
#include "elementnames.h"
#include "version.h"

/*File nemala.cpp was developed by Humeshkar B Nemala as a part of his thesis work during the Summer of 2007*/
/* Initially the code has been developed to read in energy levels,radiative and
 * collisional data from the CHIANTI and LEIDEN databases. The idea is to extend it to more databases.
 * In the case of the Leiden database there is a single .dat file which has the energy levels information,
 * radiative and collisional data, with the data corresponding to each collider coming one after the other.
 * In the case of CHIANTI, the energy levels data, radiative data and collision data are present in seperate files.
 * While LEIDEN gives collisional rate coefficients, CHIANTI gives collisional strengths.
 * In the case of CHIANTI only two colliders are used:electrons and protons. They appear as separate files.
 * The electron collision strengths files are always expected to be there. A flag is set and data processed 
 * if the file on proton collision strengths is available.*/

/* There is an initialization file called species.ini which tells Cloudy what kind of data is to be used */
/* Structures are created separately to hold the transition data,radiative and collisional data */
/* The collisional structures are different for different databases depending upon whether */
/* collisional strengths or collisional rate coefficients are used.Finally a superstructure is constructed to hold */
/* the total collisional rate obtained by considering all the colliders */
/* The colliders considered are electron,proton,Atomic Hydrogen,He,He+,He++,Ortho Molecular Hydrogen,Para Molecular Hydrogen and Molecular Hydrogen */
STATIC void states_popfill(void);
STATIC void states_nelemfill(void);
STATIC void database_prep(int);
STATIC void set_fractionation( species *sp );
STATIC void states_propprint(void);
/*SpeciesJunk set all elements of species struc to dangerous values */
STATIC void SpeciesJunk( species *sp );

#define DEBUGSTATE false
void database_readin( void )
{
	int i,intNoSp;

	FILE *ioMASTERLIST, *ioVERSION;

	char *chToken;

	char chLine[FILENAME_PATH_LENGTH_2],
		chDLine[FILENAME_PATH_LENGTH_2],
		chPath[FILENAME_PATH_LENGTH_2] = "";

	const int MAX_NUM_SPECIES = 1000;

	char chLabels[MAX_NUM_SPECIES][CHARS_SPECIES];
	char chLabelsOrig[MAX_NUM_SPECIES][CHARS_SPECIES];
	char chPaths[MAX_NUM_SPECIES][FILENAME_PATH_LENGTH_2];

	static int nCalled = 0;
	long nSpeciesLAMDA, nSpeciesSTOUT, nSpeciesCHIANTI;

	DEBUG_ENTRY( "database_readin()" );

	/* only do this once. */
	if( nCalled > 0 )
	{
		return;
	}

	/* this is first call, increment the nCalled counterso never do this again */
	++nCalled;

	// read masterlists, count number of species
	nSpecies = 0;

	////////////////////////////////////
	//  
	// Read LAMDA masterlist 
	//
	//////////////////////////////////

	/* count how many lines are in the file, ignoring all lines
	 * starting with '#':This would give the number of molecules */
	nSpeciesLAMDA = 0;

	if( atmdat.lgLamdaOn )
	{
		long numModelsNotUsed = 0;
		strcpy( chPath, "lamda" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );
		ioMASTERLIST = open_data( chPath, "r" );

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of LAMDA masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{	
			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{	
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \t\n");
				if( findspecies( chToken ) != null_mole  || 
					( chToken[1]=='-' && findspecies( chToken+2 ) != null_mole ) )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
					ASSERT( strlen(chToken) < CHARS_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

					// path is, for example, LAMDA/no.dat
					strcpy( chPaths[nSpecies], "lamda" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					chToken = strtok( NULL," \t\n" );
					strcat( chPaths[nSpecies], chToken );
					++nSpecies;
					++nSpeciesLAMDA;
				}
				else
					++numModelsNotUsed;
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		/* \todo 1 - save this and stuff as note since not really a "PROBLEM" but worth reporting */
		//if( !t_version::Inst().lgRelease && numModelsNotUsed > 0 )
		//	fprintf( ioQQQ, "\n PROBLEM - %li LAMDA models could not be found in chemistry network.\n\n\n", numModelsNotUsed );

		fclose(ioMASTERLIST);
	}

	//////////////////////////////////
	//  
	// Read CDMS/JPL masterlist
	//
	// These data files are in LAMDA format
	//
	//////////////////////////////////

	if( atmdat.lgCalpgmOn )
	{
		long numModelsNotUsed = 0;
		strcpy( chPath, "cdms+jpl" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );
		ioMASTERLIST = open_data( chPath, "r" );

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CDMS/JPL masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{	
			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{	
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \t\n");
				// hacks for alternative dialects...
				if( strcmp( chToken, "SH" ) == 0 )
					strcpy( chToken, "HS" );
				if( strcmp( chToken, "SH+" ) == 0 )
					strcpy( chToken, "HS+" );
				if( strcmp( chToken, "CCH" ) == 0 )
					strcpy( chToken, "C2H" );
				if( findspecies( chToken ) != null_mole  || 
					( chToken[1]=='-' && findspecies( chToken+2 ) != null_mole ) )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesLAMDA + 1 <= MAX_NUM_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					chLabels[nSpecies][CHARS_SPECIES-1] = '\0';

					strcpy( chPaths[nSpecies], "cdms+jpl" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					chToken = strtok( NULL," \t\n" );
					strcat( chPaths[nSpecies], chToken );
					++nSpecies;
					++nSpeciesLAMDA;
				}
				else
					++numModelsNotUsed;
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		// if( !t_version::Inst().lgRelease && numModelsNotUsed > 0 )
		//	fprintf( ioQQQ, "\nPROBLEM - %li CDMS/JPL models could not be found in chemistry network.\n\n",
		//		 numModelsNotUsed );

		fclose(ioMASTERLIST);
	}

	////////////////////////////////////
	//
	// Read STOUT masterlist and VERSION
	//
	///////////////////////////////////
	nSpeciesSTOUT = 0;
	if( atmdat.lgStoutOn )
	{
		// default location of Stout masterlist file
		strcpy( chPath, "stout" );
		strcat( chPath, input.chDelimiter );
		strcat( chPath, "masterlist" );
		strcat( chPath, input.chDelimiter );

		strcat( chPath, atmdat.chStoutFile );

		// first try local directory, then data/SED
		if( (ioMASTERLIST = fopen( atmdat.chStoutFile , "r" ) ) == NULL )
		{
			ioMASTERLIST = open_data( chPath, "r" );
		}

		// magic number
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of stout.ini.\n");
			cdEXIT(EXIT_FAILURE);
		}

		bool lgEOLST;
		long int ipST = 1;
		long int nYrRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOLST);
		long int nMonRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOLST);
		long int nDayRdST = (long)FFmtRead(chLine,&ipST,sizeof(chLine),&lgEOLST);

		static long int nYrST =11 , nMonST = 10, nDayST = 25;
		if( ( nYrRdST != nYrST ) || ( nMonRdST != nMonST ) || ( nDayRdST != nDayST ) )
		{
			fprintf( ioQQQ,
				" I expected to find the number %2.2li %2.2li %2.2li and got %2.2li %2.2li %2.2li instead.\n" ,
				nYrST , nMonST , nDayST , nYrRdST , nMonRdST , nDayRdST );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CHIANTI masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{
			strcpy(chDLine, chLine);
			chToken = strtok(chDLine," \n");
			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{
				ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
				ASSERT( nSpeciesSTOUT + 1 <= MAX_NUM_SPECIES );
				strcpy( chLabels[nSpecies], chToken );
				strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies] );

				char *chElement, chTokenTemp[7];
				strcpy( chTokenTemp, chToken );
				chElement = strtok(chTokenTemp," \n");
				chElement = strtok(chTokenTemp,"_");
				uncaps( chElement );

//					printf("STOUT:%s\n",chLabels[nSpecies]);

				// path is, for example, CHIANTI/ar/ar_10/ar_10
				// we will append extensions later
				strcpy( chPaths[nSpecies], "stout" );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chElement );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chLabels[nSpecies] );
				strcat( chPaths[nSpecies], input.chDelimiter );
				strcat( chPaths[nSpecies], chLabels[nSpecies] );

				ASSERT( isalpha(chToken[0]) );
				long cursor=0;
				chLabels[nSpecies][0] = chToken[0];
				if( isalpha(chToken[1]) )
				{
					chLabels[nSpecies][1] = chToken[1];
					cursor = 2;
				}
				else
				{
					chLabels[nSpecies][1] = ' ';
					cursor = 1;
				}

				ASSERT( chToken[cursor++]=='_' );
				ASSERT( isdigit(chToken[cursor]) );

				if( isdigit(chToken[cursor+1]) )
				{
					chLabels[nSpecies][2] = chToken[cursor++];
					chLabels[nSpecies][3] =	chToken[cursor++];
				}
				else
				{
					chLabels[nSpecies][2] = ' ';
					chLabels[nSpecies][3] =	chToken[cursor++];
				}
				chLabels[nSpecies][4] = '\0';
				ASSERT( chToken[cursor]=='\0' || chToken[cursor]=='d' );

				// now capitalize the first letter
				chLabels[nSpecies][0] = toupper( chLabels[nSpecies][0] );
				++nSpecies;
				++nSpeciesSTOUT;
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );
		fclose(ioMASTERLIST);
	}


	////////////////////////////////////
	//  
	// Read CHIANTI masterlist and VERSION
	//
	///////////////////////////////////

	nSpeciesCHIANTI = 0;

	if( atmdat.lgChiantiOn )
	{
		char chPathSave[FILENAME_PATH_LENGTH_2];
		strcpy( chPath, "chianti" );
		strcat( chPath, input.chDelimiter );
		//Preserve the path /chianti/ with chPathSave
		//Start reading in the chianti version number
		strcpy( chPathSave , chPath );
		strcat(chPath,"VERSION");
		ioVERSION = open_data(chPath,"r");
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioVERSION ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of the Chianti VERSION.\n");
			cdEXIT(EXIT_FAILURE);
		}
		fclose(ioVERSION);
		// chianti version - string since can contain letters
		strncpy(atmdat.chVersion,chLine,atmdat.iVersionLength);
		// remove newline we may captured
		long len = strlen(atmdat.chVersion);
		if( atmdat.chVersion[len-1] == '\n' )
			atmdat.chVersion[len-1] = '\0';
		// may have been null earlier in string, but make sure terminated at limit
		atmdat.chVersion[atmdat.iVersionLength-1] = '\0';
		//Restore the previous chPath
		strcpy(chPath,chPathSave);
		// Read in the masterlist
		strcat( chPath, "masterlist" );
		strcat( chPath, input.chDelimiter );
		// save copy
		strcpy( chPathSave , chPath );

		// our subset of Chianti
		strcat( chPath, atmdat.chCloudyChiantiFile );

		// first try local directory, then data/chianti
		if( (ioMASTERLIST = fopen( atmdat.chCloudyChiantiFile , "r" ) ) == NULL )
		{
			ioMASTERLIST = open_data( chPath, "r" );
		}

		// magic number
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CloudyChianti.ini.\n");
			cdEXIT(EXIT_FAILURE);
		}

		bool lgEOL;
		long int ip = 1;
		long int nYrRd = (long)FFmtRead(chLine,&ip,sizeof(chLine),&lgEOL);
		long int nMonRd = (long)FFmtRead(chLine,&ip,sizeof(chLine),&lgEOL);
		long int nDayRd = (long)FFmtRead(chLine,&ip,sizeof(chLine),&lgEOL);

		static long int nYr=11 , nMon = 10, nDay = 3;
		if( ( nYrRd != nYr ) || ( nMonRd != nMon ) || ( nDayRd != nDay ) )
		{
			fprintf( ioQQQ,
				" database_readin: the version of CloudyChianti.ini is not the current version.\n" );
			fprintf( ioQQQ,
				" database_readin obtain the current version from the Cloudy web site.\n" );
			fprintf( ioQQQ,
				" I expected to find the number %2.2li %2.2li %2.2li and got %2.2li %2.2li %2.2li instead.\n" ,
				nYr , nMon , nDay , nYrRd , nMonRd , nDayRd );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}

		if( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) == NULL )
		{
			fprintf( ioQQQ, " database_readin could not read first line of CHIANTI masterlist.\n");
			cdEXIT(EXIT_FAILURE);
		}

		do
		{

			if ((chLine[0]!='#') && (chLine[0]!='\n')&&(chLine[0]!='\t')&&(chLine[0]!='\r'))
			{
				strcpy(chDLine, chLine);
				chToken = strtok(chDLine," \n");

				fixit();  //insert logic here to exclude some ions (for example, iso sequences)
				// exclude for now the satellite lines (denoted by a "d" after the label
				//if( chToken[3]=='d' || chToken[4]=='d' || chToken[5]=='d' )
				if( chToken[3]!='d' && chToken[4]!='d' && chToken[5]!='d' )
				{
					ASSERT( nSpecies + 1 <= MAX_NUM_SPECIES );
					ASSERT( nSpeciesCHIANTI + 1 <= MAX_NUM_SPECIES );
					strcpy( chLabels[nSpecies], chToken );
					strcpy( chLabelsOrig[nSpecies], chLabels[nSpecies]);

					bool skipSpecies = false;

					//Check for duplicate species with Stout
					for( int j = nSpeciesLAMDA; j < nSpecies; j++)
					{
						if( strcmp( chLabelsOrig[j], chLabelsOrig[nSpecies] ) == 0)
						{
							printf("Skipping the Chianti version of %s, using Stout version\n",chLabels[nSpecies]);
							skipSpecies = true;
							break;
						}
					}
					if( skipSpecies)
						continue;

					char *chElement, chTokenTemp[7];
					strcpy( chTokenTemp, chToken );
					chElement = strtok(chTokenTemp," \n");
					chElement = strtok(chTokenTemp,"_");
					uncaps( chElement );

					// path is, for example, CHIANTI/ar/ar_10/ar_10
					// we will append extensions later
					strcpy( chPaths[nSpecies], "chianti" );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chElement );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chLabels[nSpecies] );
					strcat( chPaths[nSpecies], input.chDelimiter );
					strcat( chPaths[nSpecies], chLabels[nSpecies] );

					ASSERT( isalpha(chToken[0]) );
					long cursor=0;
					chLabels[nSpecies][0] = chToken[0];
					if( isalpha(chToken[1]) )
					{
						chLabels[nSpecies][1] = chToken[1];
						cursor = 2;
					}
					else
					{
						chLabels[nSpecies][1] = ' ';
						cursor = 1;
					}

					ASSERT( chToken[cursor++]=='_' );
					ASSERT( isdigit(chToken[cursor]) );

					if( isdigit(chToken[cursor+1]) )
					{
						chLabels[nSpecies][2] = chToken[cursor++];
						chLabels[nSpecies][3] =	chToken[cursor++];
					}
					else
					{
						chLabels[nSpecies][2] = ' ';
						chLabels[nSpecies][3] =	chToken[cursor++];
					}
					chLabels[nSpecies][4] = '\0';
					ASSERT( chToken[cursor]=='\0' || chToken[cursor]=='d' );

					// now capitalize the first letter
					chLabels[nSpecies][0] = toupper( chLabels[nSpecies][0] );
					++nSpecies;
					++nSpeciesCHIANTI;
				}
			}
		}
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioMASTERLIST ) != NULL );

		fclose(ioMASTERLIST);
	}

	/* no species found, nothing to do */
	if( nSpecies==0 )
		return;

	/*Initialization of the dBaseSpecies Structure*/
	dBaseSpecies = (species *)MALLOC( (unsigned long)nSpecies*sizeof(species));

	/*Initialization of the collisional rates array structure*/
	AtmolCollRateCoeff.reserve( nSpecies );
	AtmolCollSplines = (CollSplinesArray****)MALLOC((unsigned long)nSpecies *sizeof(CollSplinesArray***));
	StoutCollData = (StoutColls****)MALLOC((unsigned long)nSpecies *sizeof(StoutColls***));

	/*Mallocing here takes care of the number of colliders*/
	for( i=0; i<nSpecies; i++ )
	{
		AtmolCollRateCoeff.reserve( i, ipNCOLLIDER );
	}
	AtmolCollRateCoeff.alloc();

	// malloc state and transition arrays
	dBaseStates.resize(nSpecies);
	ipdBaseTrans.resize(nSpecies);

	for( i = 0; i < nSpecies; i++ )
	{
		dBaseTrans.push_back(TransitionList("dBaseTrans",&dBaseStates[i]));
		SpeciesJunk( &dBaseSpecies[i] );
		// label should be a minimum of 4 characters long
		size_t los = max(4,strlen(chLabels[i]));
		ASSERT( los >= 4 && los <= 7 );
		dBaseSpecies[i].chLabel = new char[los+1];
		strcpy(dBaseSpecies[i].chLabel,chLabels[i]);
		dBaseSpecies[i].lgActive = true;

		/* set type and isotopologue fractions */
		set_fractionation( &dBaseSpecies[i] );

		// set_fractionation trims off "p-","o-", etc.  Now have set label.  Check size.
		los = (int)strlen( dBaseSpecies[i].chLabel );
		ASSERT( los < CHARS_SPECIES );

		// pad label to at least four characters.
		if( dBaseSpecies[i].chLabel[2]=='\0' )
		{
			dBaseSpecies[i].chLabel[2]=' ';
			dBaseSpecies[i].chLabel[3]=' ';
			dBaseSpecies[i].chLabel[4]='\0';
		}
		else if( dBaseSpecies[i].chLabel[3]=='\0' )
		{
			dBaseSpecies[i].chLabel[3]=' ';
			dBaseSpecies[i].chLabel[4]='\0';
		}

		if( i<nSpeciesLAMDA )
		{
			// Read in LAMDA data files
			atmdat_LAMDA_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT )
		{
			atmdat_STOUT_readin( i, chPaths[i] );
		}
		else if( i < nSpeciesLAMDA + nSpeciesSTOUT + nSpeciesCHIANTI )
		{
			// Read in CHIANTI data files
			atmdat_CHIANTI_readin( i, chPaths[i] );
		}
		else
			TotalInsanity();
	}

	states_popfill();
	states_nelemfill();

	/*Setting nelem of the states to an arbitrary value*/
	/*Loop over species*/
	for( intNoSp=0; intNoSp<nSpecies; intNoSp++ )
	{
		database_prep(intNoSp);
		AllTransitions.push_back(dBaseTrans[intNoSp]);
	}

	/*To print the states*/
	if(DEBUGSTATE)
		states_propprint();
	return;
}

STATIC void set_fractionation( species *sp )
{
	DEBUG_ENTRY("set_fractionation()");

	char chToken[3];

	sp->fracIsotopologue = 1.f;
	//types include "p-", "o-", "e-", and "a-"
	strncpy( chToken, sp->chLabel, 2 );
	chToken[2] = '\0';
	if( strcmp( "p-", chToken )==0 )
		sp->fracType = 0.25f;
	else if( strcmp( "o-", chToken )==0 )
		sp->fracType = 0.75f;
	else if( strcmp( "e-", chToken )==0 )
		sp->fracType = 0.5f;
	else if( strcmp( "a-", chToken )==0 ) 
		sp->fracType = 0.5f;
	else
		sp->fracType = 1.0f;

	fixit(); // what fraction should e-type and a-type Methanol have?  Assume 50/50 for now.

	// Now scrape the type specifier off the label.
	if( sp->chLabel[1]=='-')
		memmove(sp->chLabel,sp->chLabel+2,strlen(sp->chLabel+2)+1);

	return;
}

/*This function zeros the population of all states */
STATIC void states_popfill( void)
{
	DEBUG_ENTRY( "states_popfill()" );

	for( long i=0; i<nSpecies; i++)
	{
		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++)
		{
			dBaseStates[i][j].Pop() = 0.;
		}
	}
	return;
}

/*This function fills the nelem and IonStg fields */
STATIC void states_nelemfill(void)
{
	DEBUG_ENTRY( "states_nelemfill()" );

	for( long i=0; i<nSpecies; i++ )
	{
		long nelem = 0, IonStg;
		char chLabelChemical[CHARS_SPECIES];

		if( dBaseSpecies[i].lgMolecular )
		{
			fixit(); 
			/* these should never be used if lgMolecular
			 *set to dangerous values instead of unity. */
			nelem = -1;
			IonStg = -1;
			strcpy( chLabelChemical, dBaseSpecies[i].chLabel );
		}
		else
		{
			char chToken[3];
			strncpy( chToken, dBaseSpecies[i].chLabel, 2 );
			chToken[2] = '\0';
			strcpy( chLabelChemical, chToken );
			if( chLabelChemical[1]==' ' )
				chLabelChemical[1] = '\0';
			for( long ipElement=0; ipElement<LIMELM; ipElement++ )
			{
				if( strcmp( elementnames.chElementSym[ipElement], chToken )==0 )
				{
					nelem = ipElement + 1;
					break;
				}
			}
			ASSERT( nelem > 0 && nelem <= LIMELM );
			strncpy( chToken, dBaseSpecies[i].chLabel + 2, 2 );
			IonStg = atoi(chToken);
			char chStage[5] = {'\0'};
			if( IonStg==2 )
				sprintf( chStage, "+" );
			else if( IonStg>1 )
				sprintf( chStage, "+%li", IonStg-1 );
			strcat( chLabelChemical, chStage );
			ASSERT( IonStg >= 1 && IonStg <= nelem+1 );
			//Prevent importing of iso-sequences from Chianti
			if( nelem - IonStg < NISO )
			{
				fprintf(ioQQQ, " PROBLEM: Cannot use Chianti model for %s%li\n",elementnames.chElementSym[nelem-1],IonStg);
				fprintf(ioQQQ, "  Iso-sequences are handled by our own model.\n");
				cdEXIT(EXIT_FAILURE);
			}
			dBaseSpecies[i].fmolweight = dense.AtomicWeight[nelem-1];
			// do not evaluate our cooling if we are using Chianti for this species

			if( dBaseTrans[i].chLabel() == "Chianti" )
			{
				dense.lgIonChiantiOn[nelem-1][IonStg-1] = true;
			}
			else if( dBaseTrans[i].chLabel() == "Stout" )
			{
				dense.lgIonStoutOn[nelem-1][IonStg-1] = true;
			}
			else
			{
				TotalInsanity();
			}

			if( atmdat.lgChiantiHybrid || atmdat.lgStoutHybrid )
			{
				// used in cool_dima to indicate whether to include line
				// with shorter wl than these databases
				dense.maxWN[nelem-1][IonStg-1] = dBaseSpecies[i].maxWN;
			}
			else
			{
				dense.maxWN[nelem-1][IonStg-1] = 0.;
			}
		}

		molecule *sp = findspecies(chLabelChemical);
		if( sp == null_mole )
		{
			dBaseSpecies[i].index = INT_MAX;
			if( nelem-1 >= ipHYDROGEN && dense.lgElmtOn[nelem-1] )
				fprintf(ioQQQ," PROBLEM: could not find species %li - %s\n",i,
					chLabelChemical );
		}
		else
		{
			dBaseSpecies[i].index = sp->index;
			mole.species[ sp->index ].levels = &dBaseStates[i];
			mole.species[ sp->index ].lines = &dBaseTrans[i];
		}

		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			dBaseStates[i][j].nelem() = nelem;
			dBaseStates[i][j].IonStg() = IonStg;
		}
	}
	return;
}

/*This function prints the various properties of states*/
STATIC void states_propprint(void)
{
	DEBUG_ENTRY( "states_propprint()" );

	for( long i=0; i<nSpecies; i++ )
	{
		printf("The species is %s \n",dBaseSpecies[i].chLabel);
		printf("The data output is in the following format \n");
		printf("Label Energy St.wt Pop Lifetime\n");

		for( long j=0; j<dBaseSpecies[i].numLevels_max; j++ )
		{
			printf("This is the %ld state \n",j);
			printf("%s %f %f %f %e \n",dBaseStates[i][j].chLabel(),
				dBaseStates[i][j].energy().WN(),
				dBaseStates[i][j].g(),
				dBaseStates[i][j].Pop(),
				dBaseStates[i][j].lifetime());
		}
	}
	return;
}


STATIC void database_prep(int intSpIndex)
{
	vector<realnum> fsumAs(dBaseSpecies[intSpIndex].numLevels_max,SMALLFLOAT);

	DEBUG_ENTRY( "database_prep()" );

	/*Get the lifetimes*/
	for( EmissionList::iterator em = dBaseTrans[intSpIndex].Emis().begin();
		  em != dBaseTrans[intSpIndex].Emis().end(); ++em) 
	{
		fsumAs[(*em).Tran().ipHi()] += (*em).Aul();
		(*em).iRedisFun() = ipPRD;
	}
	
	dBaseStates[intSpIndex][0].lifetime()= BIGFLOAT;
	for( int ipHi=1; ipHi < dBaseSpecies[intSpIndex].numLevels_max; ipHi++ )
	{
		dBaseStates[intSpIndex][ipHi].lifetime() = 1./fsumAs[ipHi];
	}
	return;
}

/*SpeciesJunk set all elements of species struc to dangerous values */
STATIC void SpeciesJunk( species *sp )
{
	sp->chLabel = NULL;
	set_NaN(sp->fmolweight);
	set_NaN(sp->fracIsotopologue );
	set_NaN(sp->fracType);
	sp->lgMolecular = false;
	sp->numLevels_local = -INT_MAX;
	sp->numLevels_max = -INT_MAX;

	return;
}
