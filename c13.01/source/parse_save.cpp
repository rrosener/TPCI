/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseSave parse the save command */
/*SaveFilesInit initialize save file pointers, called from cdInit */
/*CloseSaveFiles close all save files */
/*ChkUnits check for keyword UNITS on line, then scan wavelength or energy units if present */
#include "cddefines.h"
#include "cddrive.h"
#include "physconst.h"
#include "elementnames.h"
#include "input.h"
#include "geometry.h"
#include "prt.h"
#include "optimize.h"
#include "rfield.h"
#include "hcmap.h"
#include "atomfeii.h"
#include "h2.h"
#include "mole.h"
#include "hmi.h"
#include "version.h"
#include "grainvar.h"
#include "parse.h"
#include "grid.h"
#include "save.h"
#include "mpi_utilities.h"
#include "parser.h"

/* check for keyword UNITS on line, then scan wavelength or energy units if present */
STATIC void ChkUnits(Parser &p);

/* NB NB NB NB NB NB NB NB NB NB
 *
 * if any "special" save commands are added (commands that copy the file pointer
 * into a separate variable, e.g. like SAVE _DR_), be sure to add that file pointer
 * to SaveFilesInit and CloseSaveFiles !!!
 *
 * SAVE FILE POINTERS SHOULD NEVER BE ALTERED BY ROUTINES OUTSIDE THIS MODULE !!!
 *
 * hence initializations of save file pointers should never be included in zero() !!
 * the pointer might be lost without the file being closed
 * 
 * NB NB NB NB NB NB NB NB NB NB */

/* save file header headers - these are written into the string save.chHeader[save.nsave] when
 * the command is parsed
 * save.lgPunHeader[] determines whether header is saved 
 */


void ParseSave(Parser& p)
{
	char chLabel[INPUT_LINE_LENGTH] ,
		chFilename[INPUT_LINE_LENGTH] ,
		chSecondFilename[INPUT_LINE_LENGTH];
	bool lgSecondFilename;
	/* pointer to column across line image for free format number reader*/
	long int i,
	  nelem;

	char chTemp[MAX_HEADER_SIZE];

	DEBUG_ENTRY( "ParseSave()" );

	/* check that limit not exceeded */
	if( save.nsave >= LIMPUN )
	{
		fprintf( ioQQQ, 
			"The limit to the number of SAVE options is %ld.  Increase "
			"LIMPUN in save.h if more are needed.\nSorry.\n", 
		  LIMPUN );
		cdEXIT(EXIT_FAILURE);
	}

	/* initialize this flag, forced true for special cases below (e.g. for FITS files) */
	save.lgSaveToSeparateFiles[save.nsave] = p.nMatch("SEPA");

	/* LAST keyword is an option to save only on last iteration */
	save.lgPunLstIter[save.nsave] = p.nMatch("LAST");

	/* get file name for this save output.
	 * GetQuote does the following -
	 * first copy original version of file name into chLabel, 
	 * string does include null termination.
	 * set filename in OrgCard and second parameter to spaces so 
	 * that not picked up below as keyword
	 * last parameter says whether to abort if no quote found 	 */
	if( p.GetQuote( chLabel , true ) )
		/* this can't happen since routine would not return at all if no double quotes found */
		TotalInsanity();

	/* check that name is not same as opacity.opc, a special file */
	if( strcmp(chLabel , "opacity.opc") == 0 )
	{
		fprintf( ioQQQ, "ParseSave will not allow save file name %s, please choose another.\nSorry.\n",
			chLabel);
		cdEXIT(EXIT_FAILURE);
	}
	else if( chLabel[0]=='\0' )
	{
		fprintf( ioQQQ, "ParseSave found a null file name between double quotes, please check command line.\nSorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* now copy to chFilename, with optional grid prefix first */
	strcpy( chFilename , save.chGridPrefix.c_str() );
	/* this is optional prefix, normally a null string, set with set save prefix command */
	strcat( chFilename , save.chFilenamePrefix.c_str() );
	strcat( chFilename , chLabel );

	/* there may be a second file name, and we need to get it off the line
	 * before we parse options, last false parameter says not to abort if
	 * missing - this is not a problem at this stage */
	if( p.GetQuote( chSecondFilename , false ) )
		lgSecondFilename = false;
	else
		lgSecondFilename = true;

	/* CLOBBER clobber keyword is an option to overwrite rather than
	 * append to a given file */
	if( p.nMatch("CLOB") )
	{
		if( p.nMatch(" NO ") )
		{
			/* do not clobber files */
			save.lgNoClobber[save.nsave] = true;
		}
		else
		{
			/* clobber files */
			save.lgNoClobber[save.nsave] = false;
		}
	}

	/* put version number and title of model on output file, but only if
	 * this is requested with a "title" on the line*/
	/* >>chng 02 may 10, invert logic from before - default had been title */
	/* put version number and title of model on output file, but only if
	 * there is not a "no title" on the line*/
	if( !p.nMatch("NO TI") && p.nMatch("TITL"))
	{
		sprintf( save.chHeader[save.nsave], 
			"# %s %s\n", 
		  t_version::Inst().chVersion, input.chTitle );
	}

	/* usually results for each iteration are followed by a series of
	 * hash marks, ####, which fool excel.  This is an option to not print
	 * the line.  If the keyword NO HASH no hash appears the hash marks 
	 * will not occur */
	if( p.nMatch("NO HA") )
		save.lgHashEndIter[save.nsave] = false;

	/* save opacity must come first since elements appear as sub-keywords */
	if( p.nMatch("OPAC") )
	{
		/* check for keyword UNITS on line, then scan wavelength or energy units if present,
		 * units are copied into save.chConPunEnr */
		ChkUnits(p);

		strcpy( save.chSave[save.nsave], "OPAC" );

		/* "every" option to save this on every zone -
		 * not present then only last zone is saved */
		if( p.nMatch("EVER" ) )
		{
			/* save every zone */
			save.lgSaveEveryZone[save.nsave] = true;
			save.nSaveEveryZone[save.nsave] = 1;
		}
		else
		{
			/* only save last zone */
			save.lgSaveEveryZone[save.nsave] = false;
			save.nSaveEveryZone[save.nsave] = 1;
		}

		if( p.nMatch("TOTA") )
		{
			/* DoPunch will call save_opacity to parse the subcommands
			 * save total opacity */
			strcpy( save.chOpcTyp[save.nsave], "TOTL" );
			sprintf( save.chHeader[save.nsave], 
				"#nu/%s\tTot opac\tAbs opac\tScat opac\tAlbedo\telem\n",
				save.chConPunEnr[save.nsave]);
		}

		else if( p.nMatch("FIGU") )
		{
			/* do figure for hazy */
			strcpy( save.chOpcTyp[save.nsave], "FIGU" );
			sprintf( save.chHeader[save.nsave], 
				"#nu/%s\tH\tHe\ttot opac\n",
				save.chConPunEnr[save.nsave] );
		}

		else if( p.nMatch("FINE") )
		{
			/* save the fine opacity array */
			rfield.lgSaveOpacityFine = true;
			strcpy( save.chOpcTyp[save.nsave], "FINE" );
			/* check for keyword UNITS on line, then scan wavelength or energy units if present,
			 * units are copied into save.chConPunEnr */
			ChkUnits(p);

			sprintf( save.chHeader[save.nsave], 
				"#nu/%s\topac\n",
				save.chConPunEnr[save.nsave] );

			/* range option - important since so much data - usually want to
			 * only give portion of the continuum */
			if( p.nMatch("RANGE") ) 
			{
				/* get lower and upper range, eventually must be in Ryd */
				double Energy1 = p.FFmtRead();
				double Energy2 = p.FFmtRead();
				if( p.lgEOL() )
				{
					fprintf(ioQQQ,"There must be two numbers, the lower and upper energy range in Ryd.\nSorry.\n");
					cdEXIT(EXIT_FAILURE);
				}
				if( p.nMatch("UNIT" ) )
				{
					// apply units to range option
					const char *energyUnits = p.StandardEnergyUnit();
					Energy unitChange;
					unitChange.set(Energy1, energyUnits );
					Energy1 = unitChange.Ryd();
					unitChange.set(Energy2, energyUnits );
					Energy2 = unitChange.Ryd();
				}
				/* get lower and upper rang in Ryd */
				save.punarg[save.nsave][0] = (realnum)MIN2( Energy1 , Energy2 );
				save.punarg[save.nsave][1] = (realnum)MAX2( Energy1 , Energy2 );
				//fprintf(ioQQQ , "DEBUG units change fine %.3e %.3e\n" , save.punarg[save.nsave][0] ,
				//		save.punarg[save.nsave][1] );
				//cdEXIT(EXIT_FAILURE);
			}
			else
			{
				/* these mean full energy range */
				save.punarg[save.nsave][0] = 0.;
				save.punarg[save.nsave][1] = 0.;
			}
			/* optional last parameter - how many points to bring together */
			save.punarg[save.nsave][2] = (realnum)p.FFmtRead();

			/* default is to average together ten */
			if( p.lgEOL() )
				save.punarg[save.nsave][2] = 10;

			if( save.punarg[save.nsave][2] < 1 )
			{
				fprintf(ioQQQ,"The number of fine opacities to skip must be > 0 \nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("GRAI") )
		{
			/* save grain opacity command, give optical properties of gains in calculation */
			strcpy( save.chSave[save.nsave], "DUSO" );
			/* save grain opacity command in twice, here and above in opacity */
			sprintf( save.chHeader[save.nsave], 
				"#grain\tnu\tabs+scat*(1-g)\tabs\tscat*(1-g)\tscat\tscat*(1-g)/[abs+scat*(1-g)]\n" );
		}

		else if( p.nMatch("BREM") )
		{
			/* save bremsstrahlung opacity */
			strcpy( save.chOpcTyp[save.nsave], "BREM" );
			sprintf( save.chHeader[save.nsave], 
				"#nu\tbrem opac\n" );
		}

		else if( p.nMatch("SHEL") )
		{
			/* save shells, a form of the save opacity command for showing subshell crossections*/
			strcpy( save.chSave[save.nsave], "OPAC" );

			/* save subshell cross sections */
			strcpy( save.chOpcTyp[save.nsave], "SHEL" );

			/* this is element */
			save.punarg[save.nsave][0] = (realnum)p.FFmtRead();

			/* this is ion */
			save.punarg[save.nsave][1] = (realnum)p.FFmtRead();

			/* this is shell */
			save.punarg[save.nsave][2] = (realnum)p.FFmtRead();

			if( p.lgEOL() )
			{
				fprintf( ioQQQ, "There must be atom number, ion, shell\nSorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			sprintf( save.chHeader[save.nsave], 
				"#sub shell cross section\n" );
		}

		else if( p.nMatch("ELEM") )
		{
			/* save element opacity, produces n name.n files, one for each stage of 
			 * ionization.  the name is the 4-char version of the element's name, and
			 * n is the stage of ionization.  the file name on the card is ignored.
			 * The code stops in save_opacity after these files are produced. */

			/* this will be used as check that we did find an element on the command lines */
			/* nelem is -1 if an element was not found */
			if( (nelem = p.GetElem() ) < 0 )
			{
				fprintf( ioQQQ, "I did not find an element name on the opacity element command.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			/* copy string over */
			strcpy( save.chOpcTyp[save.nsave], elementnames.chElementNameShort[nelem] );
		}
		else
		{
			fprintf( ioQQQ, " I did not recognize a keyword on this save opacity command.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* save H2 has to come early since it has many suboptions */
	else if( p.nMatchErase(" H2 ") )
	{
		/* this is in mole_h2_io.c */
		h2.H2_ParseSave( p , save.chHeader[save.nsave] );
	}

	/* save HD has to come early since it has many suboptions */
	else if( p.nMatchErase(" HD ") )
	{
		/* this is in mole_h2_io.c */
		hd.H2_ParseSave( p , save.chHeader[save.nsave] );
	}

	/* save grain abundance will be handled later */
	else if( p.nMatch("ABUN") && !p.nMatch("GRAI") )
	{
		/* save abundances */
		strcpy( save.chSave[save.nsave], "ABUN" );
		sprintf( save.chHeader[save.nsave], 
			"#abund H" );
		for(nelem=ipHELIUM;nelem<LIMELM; ++nelem )
		{
			sprintf( chTemp, "\t%s",
				elementnames.chElementNameShort[nelem] );
			strcat( save.chHeader[save.nsave], chTemp );
		}
		strcat( save.chHeader[save.nsave], "\n");
	}

	else if( p.nMatch(" AGE") )
	{
		/* save ages */
		strcpy( save.chSave[save.nsave], "AGES" );
		sprintf( save.chHeader[save.nsave], 
			"#ages depth\tt(cool)\tt(H2 dest)\tt(CO dest)\tt(OH dest)\tt(H rec)\n" );
	}

	else if( p.nMatch(" AGN") )
	{
		/* save tables needed for AGN3 */
		strcpy( save.chSave[save.nsave], " AGN" );
		/* this is the AGN option, to produce a table for AGN */

		/* charge exchange rate coefficients */
		if( p.nMatch("CHAR") )
		{
			strcpy( save.chSave[save.nsave], "CHAG" );
			sprintf( save.chHeader[save.nsave], 
				"#charge exchange rate coefficnt\n" );
		}

		else if( p.nMatch("RECO") )
		{
			/* save recombination rates for AGN3 table */
			strcpy( save.chSave[save.nsave], "RECA" );
			sprintf( save.chHeader[save.nsave], 
				"#Recom rates for AGN3 table\n" );
		}

		else if( p.nMatch("OPAC") )
		{
			/* create table for appendix in AGN */
			strcpy( save.chOpcTyp[save.nsave], " AGN" );
			strcpy( save.chSave[save.nsave], "OPAC" );
		}

		else if( p.nMatch("HECS") )
		{
			/* create table for appendix in AGN */
			strcpy( save.chSaveArgs[save.nsave], "HECS" );
			sprintf( save.chHeader[save.nsave], 
				"#AGN3 he cs \n" );
		}

		else if( p.nMatch("HEMI") )
		{
			/* HEMIS - continuum emission needed for chap 4 of AGN3 */
			strcpy( save.chSaveArgs[save.nsave], "HEMI" );

			/* check for keyword UNITS on line, then scan wavelength or energy units if present,
			 * units are copied into save.chConPunEnr */
			ChkUnits(p);
		}
		else if( p.nMatch("RECC") )
		{
			/* recombination cooling, for AGN */
			strcpy( save.chSave[save.nsave], "HYDr" );
			sprintf( save.chHeader[save.nsave], 
				"#T\tbAS\tb1\tbB\n" );
		}
		else
		{
			fprintf( ioQQQ, " I did not recognize this option on the SAVE HYDROGEN command.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("AVER") )
	{
		/* save averages */
		strcpy( save.chSave[save.nsave], "AVER" );
		/* no need to print this standard line of explanation*/
		/*sprintf( save.chHeader[save.nsave], " asserts\n" );*/

		/* actually get the averages from the input stream, and malloc the 
		 * space in the arrays 
		 * save io unit not used in read */
		parse_save_average( p, save.nsave, save.chHeader[save.nsave] );
	}

	/* save charge transfer */
	else if( p.nMatch("CHAR") && p.nMatch("TRAN") )
	{
		/* NB in SaveDo only the first three characters are compared to find this option,
		 * search for "CHA" */
		/* save charge transfer */
		strcpy( save.chSave[save.nsave], "CHAR" );
		sprintf( save.chHeader[save.nsave], 
			"#charge exchange rate coefficient\n" );
	}

	// save chianti collision strengths in physical units
	else if( p.nMatch("CHIA"))
	{
		strcpy( save.chSave[save.nsave], "CHIA" );

	}

	else if( p.nMatch("CHEM") )
	{

		if( p.nMatch( "RATE" ) )
		{
			/* >>chng 06 May 30, NPA.  Save reaction rates for selected species */
			if( lgSecondFilename )
			{
				if( p.nMatch( "DEST" ) )
					strcpy( save.chSaveArgs[save.nsave], "DEST" );
				else if( p.nMatch( "CREA" ) )
					strcpy( save.chSaveArgs[save.nsave], "CREA" );
				else if( p.nMatch( "CATA" ) )	
					strcpy( save.chSaveArgs[save.nsave], "CATA" );
				else if( p.nMatch( "ALL" ) )
					strcpy( save.chSaveArgs[save.nsave], "ALL " );
				else
					strcpy( save.chSaveArgs[save.nsave], "DFLT" );
					
				strcpy( save.chSave[save.nsave], "CHRT" );
				save.optname[save.nsave] = chSecondFilename;
				// Haven't read chemistry database yet, so put off setting up header
				//sprintf( save.chHeader[save.nsave], "#");  
			} 

			else
			{
				fprintf(ioQQQ," A species label must appear within a second set of quotes (following the output filename).\n" );
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

		}
	}

	else if( p.nMatch("COMP") )
	{
		/* save Compton, details of the energy exchange problem */
		strcpy( save.chSave[save.nsave], "COMP" );
		sprintf( save.chHeader[save.nsave], 
			"#nu, comup, comdn\n" );
	}

	else if( p.nMatch("COOL") )
	{
		/* save cooling, actually done by routine cool_save */
		if( p.nMatch("EACH") )
		{
			strcpy( save.chSave[save.nsave], "EACH");
			sprintf( save.chHeader[save.nsave], 
					 "#depth(cm)\tTemp(K)\tCtot(erg/cm3/s)\t" );
			for( int i = 0 ; i < LIMELM ; i++ )
			{
				strcat(save.chHeader[save.nsave], elementnames.chElementSym[i] );
				strcat(save.chHeader[save.nsave], "\t" );
			}
			strcat(save.chHeader[save.nsave], "molecule\tdust\tH2cX\tCT C\tH-fb\tH2ln\tHDro\tH2+ \tFFcm\thvFB\teeff\tComp\tExtr\tExpn\tCycl\tHvin\tdima\n" );
		}
		else
		{
			strcpy( save.chSave[save.nsave], "COOL");
			/*>>chng 06 jun 06, revise to be same as save cooling */
			sprintf( save.chHeader[save.nsave],
					 "#depth cm\tTemp K\tHtot erg/cm3/s\tCtot erg/cm3/s\tcool fracs\n" );
		}
	}

	// punch the dominant rates for a given species
	else if( p.nMatch("DOMI") && p.nMatch("RATE"))
	{			
		if( !lgSecondFilename )
		{
			fprintf( ioQQQ,"This command requires two items in quotes (a filename and a species label).  Only one set of quotes was found.\nSorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* in this case the "second filename" is really the species label. */
		strncpy( save.chSpeciesDominantRates[save.nsave], chSecondFilename, CHARS_SPECIES );

		/* save dominant rates "species" */
		strcpy( save.chSave[save.nsave], "DOMI" );
		sprintf( save.chHeader[save.nsave], 
			"#depth cm\t%s col cm-2\tsrc s-1\tsnk s-1\n", 
			save.chSpeciesDominantRates[save.nsave] );
	}

	else if( p.nMatch("DYNA") )
	{
		/* save something dealing with dynamics 
		 * in SaveDo the DYN part of key is used to call DynaPunch,
		 * with the 4th char as its second argument.  DynaSave uses that
		 * 4th letter to decide the job */
		if( p.nMatch( "ADVE" ) )
		{
			/* save information relating to advection */
			strcpy( save.chSave[save.nsave], "DYNa");
			sprintf( save.chHeader[save.nsave], 
				"#advection depth\tHtot\tadCool\tadHeat\tdCoolHeatdT\t"
				"Source[hyd][hyd]\tRate\tEnthalph\tadSpecEnthal\n" );
		}
		else
		{
			fprintf( ioQQQ, " I did not recognize a sub keyword on this SAVE DYNAMICS command.\n" );
			fprintf( ioQQQ, " Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("ENTH") )
	{
		/* contributors to the total enthalpy */
		strcpy( save.chSave[save.nsave], "ENTH" );
		sprintf( save.chHeader[save.nsave], 
			"#depth\tTotal\tExcit\tIoniz\tBind\tKE\tther+PdV\tmag \n" );
	}

	else if( p.nMatch("EXEC") && p.nMatch("TIME") )
	{
		/* output the execution time per zone */
		strcpy( save.chSave[save.nsave], "XTIM" );
		sprintf( save.chHeader[save.nsave], 
			"#zone\tdTime\tElapsed t\n" );
	}

	else if( p.nMatch("FEII") || p.nMatch("FE II") )
	{
		/* something to do with FeII atom - several options
		 * FeII column densities */
		if( p.nMatch("COLU") && p.nMatch("DENS") )
		{
			/* save FeII column density */
			strcpy( save.chSave[save.nsave], "FENl" );

			/* file will give energy, statistical weight, and column density [cm-2] */
			sprintf( save.chHeader[save.nsave], 
				"#FeII: energy\tstat wght\tcol den\n" );
		}

		/* FeII continuum, only valid if large atom is set */
		else if( p.nMatch("CONT") )
		{
			// save FeII continuum, options are total (default), inward,
			// and outward
			if( p.nMatch("INWA") )
			{
				// inward continuum
				strcpy( save.chSave[save.nsave], "FEcI" );
				//sprintf( save.chHeader[save.nsave], 
				//	"#FeII inward: Wl(A)\tInt[erg cm-2 s-1]\n" );
			}
			else if( p.nMatch(" OUT") )
			{
				// outward continuum
				strcpy( save.chSave[save.nsave], "FEcO" );
				//sprintf( save.chHeader[save.nsave], 
				//	"#FeII outward: Wl(A)\tInt[erg cm-2 s-1]\n" );
			}
			else
			{
				// total continuum
				strcpy( save.chSave[save.nsave], "FEcT" );
				//sprintf( save.chHeader[save.nsave], 
				//	"#FeII total: Wl(A)\tInt[erg cm-2 s-1]\n" );
			}

			// default units of spectral energy are Ryd, this can change to
			// many other units
			ChkUnits(p);

			// by default give numbers in two columns, row keyword says to
			// write the numbers across as one long row 
			if( p.nMatch(" ROW") )
				save.punarg[save.nsave][0] = 1;
			else
				// the default, two columns
				save.punarg[save.nsave][0] = 2;
		}

		else if( p.nMatch("DEPA") )
		{
			/* save out departure coefficient for the large FeII atom */
			sprintf( save.chHeader[save.nsave], 
				"#FeII departure coefficient \n" );
			/* optional keyword all means do all levels, if not then just do subset */
			if( p.nMatch(" ALL") )
			{
				/* save all levels, calls routine FeIIPunDepart */
				strcpy( save.chSave[save.nsave], "FE2D" );
			}
			else
			{
				/* save a very few selected levels, calls routine FeIIPunDepart */
				strcpy( save.chSave[save.nsave], "FE2d" );
			}
		}

		else if( p.nMatch("LEVE") )
		{
			/* save level energies and statistical weights for large FeII atom */
			sprintf( save.chHeader[save.nsave], 
				"#FeII energy(wn)\tstat weight\n" );
			strcpy( save.chSave[save.nsave], "FE2l" );
		}

		else if( p.nMatch("LINE") )
		{
			/* save FeII lines command
			 * three optional parameters, threshold for faintest
			 * line to print, lower and upper energy bounds */

			/* save intensities from large FeII atom */
			strcpy( save.chSave[save.nsave], "FEli" );

			/* short keyword makes save half as big */
			if( p.nMatch("SHOR") )
			{
				FeII.lgShortFe2 = true;
			}
			else
			{
				FeII.lgShortFe2 = false;
			}

			/* first optional number changes the threshold of weakest line to print*/
			/* fe2thresh is intensity relative to normalization line,
			* normally Hbeta, and is set to zero in zero.c */
			FeII.fe2thresh = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				FeII.fe2thresh = 0.;
			}

			/* it is a log if negative */
			if( FeII.fe2thresh < 0. )
			{
				FeII.fe2thresh = (realnum)pow((realnum)10.f,FeII.fe2thresh);
			}

			/* check for energy range (Rydberg) of lines to be saved,
			 * this is to limit size of output file */
			FeII.fe2ener[0] = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				FeII.fe2ener[0] = 0.;
			}

			FeII.fe2ener[1] = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				FeII.fe2ener[1] = 1e8;
			}
			/* if either is negative then both are logs */
			if( FeII.fe2ener[0] < 0. || FeII.fe2ener[1] < 0. )
			{
				FeII.fe2ener[0] = (realnum)pow((realnum)10.f,FeII.fe2ener[0]);
				FeII.fe2ener[1] = (realnum)pow((realnum)10.f,FeII.fe2ener[1]);
			}

			/* entered in Ryd in above, convert to wavenumbers */
			FeII.fe2ener[0] /= (realnum)WAVNRYD;
			FeII.fe2ener[1] /= (realnum)WAVNRYD;

			/* these results are actually created by the FeIISaveLines routine
			 * that lives in the FeIILevelPops file */
			sprintf( save.chHeader[save.nsave], 
				"#FeII ipHi\tipLo\tWL(A)\tlog I\tI/Inorm\t\tTau\n" );
		}

		else if( p.nMatch("OPTI") && p.nMatch("DEPT") )
		{
			/* save optical depths for large FeII atom */
			sprintf( save.chHeader[save.nsave], 
				"#FeII hi\tlow\twl(A)\ttau\n" );
			strcpy( save.chSave[save.nsave], "FE2o" );
		}

		else if( p.nMatch("POPU") )
		{
			/* save out level populations for the large FeII atom */
			sprintf( save.chHeader[save.nsave], 
				"#FeII level populations [cm^-3]\n" );

			/* this is keyword RELATIVE that says to save relative to total Fe+, 
			 * default is actual density (cm-3) */
			if( p.nMatch("RELA") )
			{
				save.punarg[save.nsave][2] = 0.;
			}
			else
			{
				/* default is to save density (cm-3) */
				save.punarg[save.nsave][2] = 1.;
			}

			/* optional keyword all means do all levels, if not then just do subset */
			if( p.nMatch(" ALL") )
			{
				/* save all levels, calls routine FeIIPunPop */
				strcpy( save.chSave[save.nsave], "FE2P" );
				save.punarg[save.nsave][0] = 0.;
				save.punarg[save.nsave][1] = (realnum)NFE2LEVN;
			}

			/* optional keyword range means read lower and upper bound, do these */
			else if( p.nMatch("RANG") )
			{
				/* save range of levels, calls routine FeIIPunPop */
				strcpy( save.chSave[save.nsave], "FE2P" );
				save.punarg[save.nsave][0] = (realnum)p.FFmtRead();
				save.punarg[save.nsave][1] = (realnum)p.FFmtRead();
				if( p.lgEOL() || save.punarg[save.nsave][0] <0 ||
					save.punarg[save.nsave][0]>= save.punarg[save.nsave][1] )
				{
					fprintf( ioQQQ, "There must be two numbers on this save "
						"FeII populations range command.\n" );
					fprintf( ioQQQ, "These give the lower and upper levels "
						"for the range of FeII levels.\n" );
					fprintf( ioQQQ, "The first, %g, must be less than the second, %g.\n",
						save.punarg[save.nsave][0],
						save.punarg[save.nsave][1]);
					fprintf( ioQQQ, "Sorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}

			else
			{
				/* save a very few selected levels, calls routine FeIIPunPop */
				strcpy( save.chSave[save.nsave], "FE2p" );
			}
		}
		else
		{
			fprintf( ioQQQ, "There must be a second keyword on this SAVE FEII command.\n" );
			fprintf( ioQQQ, "The ones I know about are COLUmn, CONTinuum, "
				"DEPArture, LEVEls, LINE, OPTIcal DEPTh, and POPUlations.\n" );
			fprintf( ioQQQ, "Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* the save continuum command, with many options,
	 * the first 3 char of the chSave flag will always be "CON" 
	 * with the last indicating which one */
	else if( p.nMatch("CONT") && !p.nMatch("XSPE") )
	{
		/* this flag is checked in PrtComment to generate a caution
		 * if continuum is saved but iterations not performed */
		save.lgPunContinuum = true;

		/* check for keyword UNITS on line, then scan wavelength or energy units if present,
		 * units are copied into save.chConPunEnr */
		ChkUnits(p);

		if( p.nMatch("BINS") )
		{
			/* continuum binning */
			strcpy( save.chSave[save.nsave], "CONB" );

			sprintf( save.chHeader[save.nsave], 
				"#Continuum binning enrOrg/%s\tEnergy\twidth of cells\n",
				save.chConPunEnr[save.nsave] );
		}

		else if( p.nMatch("DIFF") )
		{
			/* diffuse continuum, the locally emitted lines and continuum */
			strcpy( save.chSave[save.nsave], "COND" );

			/* by default gives lines and continuum separately only for
			 * last zone.  The keyword ZONE says to give the total for every
			 * zone in one very low row */
			if( p.nMatch("ZONE") )
			{
				sprintf( save.chHeader[save.nsave], 
					"#energy/%s then emission per zone\n",
					save.chConPunEnr[save.nsave] );
				save.punarg[save.nsave][0] = 2.;

			}
			else
			{
				sprintf( save.chHeader[save.nsave], 
					"#energy/%s\tConEmitLocal\tDiffuseLineEmission\tTotal\n",
					save.chConPunEnr[save.nsave] );
				save.punarg[save.nsave][0] = 1.;
			}
		}

		else if( p.nMatch("EMIS") )
		{
			/* continuum volume emissivity and opacity as a function of radius */
			strcpy( save.chSave[save.nsave], "CONS" );

			double num = p.FFmtRead();
			if( p.lgEOL() )
				p.NoNumb( "continuum emissivity frequency" );
			save.emisfreq[save.nsave].set( num, save.chConPunEnr[save.nsave] );
			if( save.emisfreq[save.nsave].Ryd() < rfield.emm ||
			    save.emisfreq[save.nsave].Ryd() > rfield.egamry )
			{
				fprintf( ioQQQ, " The frequency is outside the Cloudy range\n Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}

			sprintf( save.chHeader[save.nsave], 
				 "#Radius\tdepth\tnujnu\tkappa_abs\tkappa_sct @ %e Ryd\n",
				 save.emisfreq[save.nsave].Ryd() );
		}

		else if( p.nMatch("EMIT") )
		{
			/* continuum emitted by cloud */
			strcpy( save.chSave[save.nsave], "CONE" );

			sprintf( save.chHeader[save.nsave], 
				"#Energy/%s\treflec\toutward\ttotal\tline\tcont\n",
				save.chConPunEnr[save.nsave] );
		}

		else if( p.nMatch("FINE" ) )
		{
			rfield.lgSaveOpacityFine = true;
			/* fine transmitted continuum cloud */
			strcpy( save.chSave[save.nsave], "CONf" );

			sprintf( save.chHeader[save.nsave], 
				"#Energy/%s\tTransmitted\n",
				save.chConPunEnr[save.nsave] );

			/* range option - important since so much data */
			if( p.nMatch("RANGE") ) 
			{
				/* get lower and upper range, eventually must be in Ryd */
				double Energy1 = p.FFmtRead();
				double Energy2 = p.FFmtRead();
				if( p.lgEOL() )
				{
					fprintf(ioQQQ,"There must be two numbers, the lower and upper energies in Ryd.\nSorry.\n");
					cdEXIT(EXIT_FAILURE);
				}
				if( p.nMatch("UNIT" ) )
				{
					// apply units to range option
					const char *energyUnits = p.StandardEnergyUnit();
					Energy unitChange;
					unitChange.set(Energy1, energyUnits );
					Energy1 = unitChange.Ryd();
					unitChange.set(Energy2, energyUnits );
					Energy2 = unitChange.Ryd();
				}
				/* get lower and upper rang in Ryd */
				save.punarg[save.nsave][0] = (realnum)MIN2( Energy1 , Energy2 );
				save.punarg[save.nsave][1] = (realnum)MAX2( Energy1 , Energy2 );
				//fprintf(ioQQQ , "DEBUG units change fine %.3e %.3e\n" , save.punarg[save.nsave][0] ,
				//		save.punarg[save.nsave][1] );
				//cdEXIT(EXIT_FAILURE);
			}
			else
			{
				/* these mean full energy range */
				save.punarg[save.nsave][0] = 0.;
				save.punarg[save.nsave][1] = 0.;
			}
			/* optional last parameter - how many points to bring together */
			save.punarg[save.nsave][2] = (realnum)p.FFmtRead();

			/* default is to bring together ten */
			if( p.lgEOL() )
				save.punarg[save.nsave][2] = 10;

			if( save.punarg[save.nsave][2] < 1 )
			{
				fprintf(ioQQQ,"The number of fine opacities to skip must be > 0 \nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}

		else if( p.nMatch("GRAI") )
		{
			/* save grain continuum in optically thin limit */
			strcpy( save.chSave[save.nsave], "CONG" );

			sprintf( save.chHeader[save.nsave], 
				"#energy\tgraphite\trest\ttotal\n" );
		}

		else if( p.nMatch("INCI") )
		{
			/* incident continuum */
			strcpy( save.chSave[save.nsave], "CONC" );

			sprintf( save.chHeader[save.nsave], 
				"#Incident Continuum, Enr\tnFn \n" );
		}

		else if( p.nMatch("INTE") )
		{
			/* continuum interactions */
			strcpy( save.chSave[save.nsave], "CONi" );

			sprintf( save.chHeader[save.nsave], 
				"#Continuum interactions, inc, otslin. otscon, ConInterOut, outlin \n" );
			/* this is option for lowest energy, if nothing then zero */
			save.punarg[save.nsave][0] = (realnum)p.FFmtRead();
		}

		else if( p.nMatch("IONI") )
		{
			/* save ionizing continuum*/
			strcpy( save.chSave[save.nsave], "CONI" );

			/* this is option for lowest energy, if nothing then zero */
			save.punarg[save.nsave][0] = (realnum)p.FFmtRead();

			/* this is option for smallest interaction to save, def is 1 percent */
			save.punarg[save.nsave][1] = (realnum)p.FFmtRead();
			if( p.lgEOL() )
				save.punarg[save.nsave][1] = 0.01f;

			/* "every" option to save this on every zone -
			 * not present then only last zone is saved */
			if( p.nMatch("EVER" ) )
			{
				/* save every zone */
				save.lgSaveEveryZone[save.nsave] = true;
				save.nSaveEveryZone[save.nsave] = 1;
			}
			else
			{
				/* only save last zone */
				save.lgSaveEveryZone[save.nsave] = false;
				save.nSaveEveryZone[save.nsave] = 1;
			}

			/* put the header at the top of the file */
			sprintf( save.chHeader[save.nsave], 
				"#cell(on C scale)\tnu\tflux\tflx*cs\tFinc\totsl\totsc\toutlin\toutcon\trate/tot\tintegral\tline\tcont\n" );
		}
#ifdef USE_NLTE7
		else if( p.nMatch("NLTE") )
		{
			/* continuum emitted by cloud */
			strcpy( save.chSave[save.nsave], "CONl" );

			sprintf( save.chHeader[save.nsave],
				"   spectrum1   NeXY6   XNUMX\n");
		}
#endif

		else if( p.nMatch("OUTW") )
		{
			/* outward only continuum */
			if( p.nMatch("LOCA") )
			{
				strcpy( save.chSave[save.nsave], "CONo" );
				sprintf( save.chHeader[save.nsave], 
					"#Local Out   ConInterOut+line SvOt*opc pass*opc\n" );
			}
			else
			{
				strcpy( save.chSave[save.nsave], "CONO" );
				sprintf( save.chHeader[save.nsave], 
					"#Out Con      OutIncid  OutConD   OutLinD   OutConS\n" );
			}
		}

		else if( p.nMatch("TRAN") )
		{
			/* transmitted continuum */
			strcpy( save.chSave[save.nsave], "CONT" );

			sprintf( save.chHeader[save.nsave], 
				"#ener\tTran Contin\ttrn coef\n" );
		}

		else if( p.nMatch(" TWO") )
		{
			/* total two photon continua rfield.TotDiff2Pht */
			strcpy( save.chSave[save.nsave], "CON2" );

			sprintf( save.chHeader[save.nsave], 
				"#energy\t n_nu\tnuF_nu \n" );
		}

		else if( p.nMatch(" RAW") )
		{
			/* "raw" continua */
			strcpy( save.chSave[save.nsave], "CORA" );

			sprintf( save.chHeader[save.nsave], 
				"#Raw Con anu\tflux\totslin\totscon\tConRefIncid\tConEmitReflec\tConInterOut\toutlin\tConEmitOut\tline\tcont\tnLines\n" );
		}

		else if( p.nMatch("REFL") )
		{
			/* reflected continuum */
			strcpy( save.chSave[save.nsave], "CONR" );

			sprintf( save.chHeader[save.nsave], 
				"#Reflected\tcont\tline\ttotal\talbedo\tConID\n" );
		}

		else
		{
			/* this is the usual save continuum command,
			 * ipType is index for continuum array to set either
			 * iteration or cumulative output */
			int ipType = 0;
			if( p.nMatch( "CUMU" ) )
				ipType = 1;
			save.punarg[save.nsave][0] = (realnum)ipType;
			strcpy( save.chSave[save.nsave], "CON " );
			char chHold[100];
			strcpy( chHold, "#Cont " );
			if( ipType > 0 )
				strcpy( chHold , "#Cumul " );
			sprintf( save.chHeader[save.nsave], 
				"%s nu\tincident\ttrans\tDiffOut\tnet trans\treflc\ttotal\treflin\toutlin\tlineID\tcont\tnLine\n" ,
				chHold );

			/* >>chng 06 apr 03, add "every" option to save this on every zone -
			 * if every is not present then only last zone is saved */
			if( p.nMatch("EVER" ) )
			{
				/* save every zone */
				save.lgSaveEveryZone[save.nsave] = true;
				/* option to say how many to skip */
				save.nSaveEveryZone[save.nsave] = (long)p.FFmtRead();
				if( p.lgEOL() )
					save.nSaveEveryZone[save.nsave] = 1;
			}
			else
			{
				/* only save last zone */
				save.lgSaveEveryZone[save.nsave] = false;
				save.nSaveEveryZone[save.nsave] = 1;
			}
		}
	}

	/* save information about convergence of this model 
	 * reason - why it did not converge an iteration
	 * error - zone by zone display of various convergence errors */
	else if( p.nMatch("CONV") )
	{
		if( p.nMatch("REAS") )
		{
			/* this does not count as a save option (really) */
			save.lgPunConv = true;
			/* this is done below */
			strcpy( save.chSave[save.nsave], "" );
			save.lgRealSave[save.nsave] = false;
		}
		else if( p.nMatch("ERRO") )
		{
			/* save zone by zone errors in pressure, electron density, and heating-cooling */
			/* convergence error */
			strcpy( save.chSave[save.nsave], "CNVE" );
			sprintf( save.chHeader[save.nsave], 
				"#depth\tnPres2Ioniz\tP(cur)\tP%%error\tNE(cor)\tNE(cur)\tNE%%error\tHeat\tCool\tHC%%error\n" );
		}
		else if( p.nMatch("BASE") )
		{
			/* save converged quantities in Converge base for each pass through
			 * solvers - not one pass per zone */
			strcpy( save.chSave[save.nsave], "CNVB" );
			strcpy( save.chSave[save.nsave], "" );
			save.lgRealSave[save.nsave] = false;
		}
		else
		{
			fprintf( ioQQQ, "There must be a second keyword on this command.\n" );
			fprintf( ioQQQ, "The ones I know about are REASON, ERROR, and BASE.\n" );
			fprintf( ioQQQ, "Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch(" DR ") )
	{
		/* first occurrence of save dr to follow choice in change of zoning */
		save.lgDROn = true;
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("ELEM") && !p.nMatch("GAMMA") && !p.nMatch("COOL") ) // do not trip on SAVE COOLING EACH ELEMENT
	{
		/* option to save ionization structure of some element
		 * will give each stage of ionization, vs depth */
		strcpy( save.chSave[save.nsave], "ELEM" );

		/* this returns element number on c scale */
		/* >>chng 04 nov 23, had converted to f scale, leave on c */
		nelem = p.GetElem();
		if( nelem < 0 || nelem >= LIMELM )
		{
			fprintf( ioQQQ, " I could not recognize a valid element name on this line.\n" );
			fprintf( ioQQQ, " Please check your input script. Bailing out...\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is the atomic number on the c scale */
		save.punarg[save.nsave][0] = (realnum)nelem;

		/* >>chng 04 nov 24, add DENSE option to print density rather than fraction */
		save.punarg[save.nsave][1] = 0;
		if( p.nMatch("DENS")  )
			save.punarg[save.nsave][1] = 1.;

		/* start printing header line - first will be the depth in cm */
		sprintf( save.chHeader[save.nsave], "#depth");

		/* next come the nelem+1 ion stages */
		for(i=0; i<=nelem+1;++i )
		{
			sprintf( chTemp, 
				"\t%2s%2li", elementnames.chElementSym[nelem],i+1);
			strcat( save.chHeader[save.nsave], chTemp );
		}

		/* finally some fine structure or molecular populations */
		/* >>chng 04 nov 23, add fs pops of C, O 
		 * >>chng 04 nov 25, add molecules */
		if( nelem==ipHYDROGEN )
		{
			sprintf( chTemp, "\tH2");
			strcat( save.chHeader[save.nsave], chTemp );
		}
		if( nelem==ipCARBON )
		{
			sprintf( chTemp, "\tC1\tC1*\tC1**\tC2\tC2*\tCO");
			strcat( save.chHeader[save.nsave], chTemp );
		}
		else if( nelem==ipOXYGEN )
		{
			sprintf( chTemp, "\tO1\tO1*\tO1**");
			strcat( save.chHeader[save.nsave], chTemp );
		}

		/* finally the new line */
		strcat( save.chHeader[save.nsave], "\n");
	}

	else if( p.nMatch("FITS") )
	{

#ifdef FLT_IS_DBL
		fprintf( ioQQQ, "Saving FITS files is not currently supported in double precision.\n" );
		fprintf( ioQQQ, "Please recompile without the FLT_IS_DBL option.\n" );
		fprintf( ioQQQ, "Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
#else
		/* say that this is a FITS file output */
		save.lgFITS[save.nsave] = true;
		/* concatenating files in a grid run would be illegal FITS */
		save.lgSaveToSeparateFiles[save.nsave] = true;
		save.lgPunLstIter[save.nsave] = true;
		save.FITStype[save.nsave] = NUM_OUTPUT_TYPES;

		strcpy( save.chSave[save.nsave], "FITS" );
#endif

	}

	else if( p.nMatch("FRED") )
	{
		/* save out some stuff for Fred's dynamics project */
		sprintf( save.chHeader[save.nsave], 
			"#Radius\tDepth\tVelocity(km/s)\tdvdr(cm/s)\thden\teden\tTemperature\tRadAccel line\tRadAccel con\t"
			"Force multiplier\ta(e thin)\t"
			"HI\tHII\tHeI\tHeII\tHeIII\tC2\tC3\tC4\tO1\t"
			"O2\tO3\tO4\tO5\tO6\tO7\tO8\t" 
			"HI\tHII\tHeI\tHeII\tHeIII\tC2\tC3\tC4\tO1\t"
			"O2\tO3\tO4\tO5\tO6\tO7\tO8\tMg2\tMg2\tOVI(1034) TauIn\tTauCon\n");

		strcpy( save.chSave[save.nsave], "FRED" );
	}

	else if( p.nMatch("GAMM") )
	{
		/* save all photoionization rates for all subshells */
		sprintf( save.chHeader[save.nsave], 
			"#Photoionization rates \n" );
		if( p.nMatch("ELEMENT") )
		{
			/* element keyword, find element name and stage of ionization, 
			 * will print photoionization rates for valence of that element */
			strcpy( save.chSave[save.nsave], "GAMe" );

			/* this returns element number on c scale */
			nelem = p.GetElem();
			/* this is the atomic number on the C scale */
			save.punarg[save.nsave][0] = (realnum)nelem;

			/* this will become the ionization stage on C scale */
			save.punarg[save.nsave][1] = (realnum)p.FFmtRead() - 1;
			if( p.lgEOL() )
				p.NoNumb("element ionization stage" );
			if( save.punarg[save.nsave][1]<0 || save.punarg[save.nsave][1]> nelem+1 )
			{
				fprintf(ioQQQ,"Bad ionization stage - please check Hazy.\nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			/* no element - so make table of all rates */
			strcpy( save.chSave[save.nsave], "GAMt" );
		}

	}
	else if( p.nMatch("GRAI") )
	{
		/* save grain ... options */
		if( p.nMatch("OPAC") )
		{
			/* check for keyword UNITS on line, then scan wavelength or energy units, 
			 * sets save.chConPunEnr*/
			ChkUnits(p);

			strcpy( save.chSave[save.nsave], "DUSO" );
			/* save grain opacity command in twice, here and above in opacity */
			sprintf( save.chHeader[save.nsave], 
				"#grain\tnu/%s\tabs+scat*(1-g)\tabs\tscat*(1-g)\tscat\tscat*(1-g)/[abs+scat*(1-g)]\n",
				save.chConPunEnr[save.nsave] );
		}
		else if( p.nMatch("ABUN") )
		{
			/* save grain abundance */
			strcpy( save.chSave[save.nsave], "DUSA" );
			sprintf( save.chHeader[save.nsave], 
				 "#grain\tdepth\tabundance (g/cm^3)\n" );
		}
		else if( p.nMatch("D/G ") )
		{
			/* save grain dust/gas mass ratio */
			strcpy( save.chSave[save.nsave], "DUSD" );
			sprintf( save.chHeader[save.nsave], 
				 "#grain\tdepth\tdust/gas mass ratio\n" );
		}
		else if( p.nMatch("PHYS") )
		{
			/* save grain physical conditions */
			strcpy( save.chSave[save.nsave], "DUSP" );
			sprintf( save.chHeader[save.nsave], 
				"#grain\tdepth\tpotential\n" );
		}
		else if( p.nMatch(" QS ") )
		{
			strcpy( save.chSave[save.nsave], "DUSQ" );
			sprintf( save.chHeader[save.nsave], 
				"#grain\tnu\tQ_abs\tQ_scat\n" );
		}
		else if( p.nMatch("TEMP") )
		{
			/* save temperatures of each grain species */
			strcpy( save.chSave[save.nsave], "DUST" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain temperature\n" );
		}
		else if( p.nMatch("DRIF") )
		{
			/* save drift velocity of each grain species */
			strcpy( save.chSave[save.nsave], "DUSV" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain drift velocity\n" );
		}
		else if( p.nMatch("EXTI") )
		{
			/* save grain extinction */
			strcpy( save.chSave[save.nsave], "DUSE" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#depth\tA_V(extended)\tA_V(point)\n" );
		}
		else if( p.nMatch("CHAR") )
		{
			/* save charge per grain (# elec/grain) for each grain species */
			strcpy( save.chSave[save.nsave], "DUSC" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain charge\n" );
		}
		else if( p.nMatch("HEAT") )
		{
			/* save heating due to each grain species */
			strcpy( save.chSave[save.nsave], "DUSH" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain heating\n" );
		}
		else if( p.nMatch("POTE") )
		{
			/* save floating potential of each grain species */
			strcpy( save.chSave[save.nsave], "DUSP" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain\tdepth\tpotential\n" );
		}
		else if( p.nMatch("H2RA") )
		{
			/* save grain H2rate - H2 formation rate for each type of grains */
			strcpy( save.chSave[save.nsave], "DUSR" );
			/* cannot save grain labels since they are not known yet */
			sprintf( save.chHeader[save.nsave], 
				"#grain H2 formation rates\n" );
		}
		else
		{
			fprintf( ioQQQ, "There must be a second key on this GRAIN command; The options I know about follow (required key in CAPS):\n");
			fprintf( ioQQQ, "OPACity, ABUNdance, D/G mass ratio, PHYSical conditions,  QS , TEMPerature, DRIFt velocity, EXTInction, CHARge, HEATing, POTEntial, H2RAtes\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("GAUN") )
	{
		strcpy( save.chSave[save.nsave], "GAUN" );
		sprintf( save.chHeader[save.nsave], 
			"#Gaunt factors.\n" );
	}
	else if( p.nMatch("GRID") )
	{
		strcpy( save.chSave[save.nsave], "GRID" );
		/* automatically generate no hash option */
		save.lgHashEndIter[save.nsave] = false;
	}
	else if( p.nMatch( "HIST" ) )
	{
		/* save pressure history of current zone */
		if( p.nMatch( "PRES") )
		{
			/* save pressure history - density - pressure for this zone */
			strcpy( save.chSave[save.nsave], "HISp" );
			sprintf( save.chHeader[save.nsave], 
				"#iter zon\tdensity\tpres cur\tpres error\n" );
		}
		/* save temperature history of current zone */
		else if( p.nMatch( "TEMP" ) )
		{
			/* save pressure history - density - pressure for this zone */
			strcpy( save.chSave[save.nsave], "HISt" );
			sprintf( save.chHeader[save.nsave], 
				"#iter zon\ttemperature\theating\tcooling\n" );
		}
	}

	else if( p.nMatch("HTWO") )
	{
		fprintf(ioQQQ," Sorry, this command has been replaced with the "
			"SAVE H2 CREATION and SAVE H2 DESTRUCTION commands.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* QHEAT has to come before HEAT... */
	else if( p.nMatch("QHEA") )
	{
		/* this is just a dummy clause, do the work below after parsing is over. 
		 * this is a no-nothing, picked up to stop optimizer */
		((void)0);
	}

	else if( p.nMatch("HEAT") )
	{
		/* save heating */
		strcpy( save.chSave[save.nsave], "HEAT" );
		/*>>chng 06 jun 06, revise to be same as save cooling */
		sprintf( save.chHeader[save.nsave], 
			"#depth cm\tTemp K\tHtot erg/cm3/s\tCtot erg/cm3/s\theat fracs\n" );
	}

	else if( p.nMatch("HELI") &&!( p.nMatch("IONI")))
	{
		/* save helium & helium-like iso sequence, but not save helium ionization rate
		 * save helium line wavelengths */
		if( p.nMatch("LINE") && p.nMatch("WAVE") )
		{
			strcpy( save.chSave[save.nsave], "HELW" );
			sprintf( save.chHeader[save.nsave], 
				"#wavelengths of lines from He-like ions\n" );
		}
		else
		{
			fprintf( ioQQQ, "save helium has options: LINE WAVElength.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
			/* no key */
		}
	}

	else if( p.nMatch("HUMM") )
	{
		strcpy( save.chSave[save.nsave], "HUMM" );
		sprintf( save.chHeader[save.nsave], 
			"#input to DHs routine.\n" );
	}

	else if( p.nMatch("HYDR") )
	{
		/* save hydrogen physical conditions */
		if( p.nMatch("COND") )
		{
			strcpy( save.chSave[save.nsave], "HYDc" );
			sprintf( save.chHeader[save.nsave], 
				"#depth\tTe\tHDEN\tEDEN\tHI/H\tHII/H\tH2/H\tH2+/H\tH3+/H\tH-/H\n" );
			/* save hydrogen ionization */
		}

		/* save information on 21 cm excitation processes - accept either keyword 21cm or 21 cm */
		else if( p.nMatch("21 CM") ||p.nMatch("21CM"))
		{
			/* save information about 21 cm line */
			strcpy( save.chSave[save.nsave], "21CM" );
			sprintf( save.chHeader[save.nsave], 
				"#depth\tT(spin)\tT(kin)\tT(Lya/21cm)\tnLo\tnHi\tOccLya\ttau(21cm)"
				"\ttau(Lya)\topac(21 cm)\tn/Ts\ttau(21)\tTex(Lya)\tN(H0)/Tspin"
				"\tSum_F0\tSum_F1\tSum_T21\n" );
		}

		else if( p.nMatch("IONI") )
		{
			/* save hydrogen ionization */
			strcpy( save.chSave[save.nsave], "HYDi" );
			sprintf( save.chHeader[save.nsave], 
				"#hion\tzn\tgam1\tcoll ion1\tRecTot\tHRecCaB\thii/hi\tSim hii/hi"
				"\time_Hrecom_long(esc)\tdec2grd\texc pht\texc col\trec eff\tsec ion\n" );
		}
		else if( p.nMatch("POPU") )
		{
			/* save hydrogen populations */
			strcpy( save.chSave[save.nsave], "HYDp" );
			sprintf( save.chHeader[save.nsave], 
				"#depth\tn(H0)\tn(H+)\tn(1s)\tn(2s)\tn(2p)\tetc\n" );
		}
		else if( p.nMatch("LINE") )
		{
			/* save hydrogen lines
			 * hydrogen line intensities and optical depths  */
			strcpy( save.chSave[save.nsave], "HYDl" );
			sprintf( save.chHeader[save.nsave], 
				"#nHi\tlHi\tnLo\tlLo\tE(ryd)\ttau\n" );
		}
		else if( p.nMatch(" LYA") )
		{
			/* save hydrogen Lya some details about Lyman alpha  */
			strcpy( save.chSave[save.nsave], "HYDL" );
			sprintf( save.chHeader[save.nsave], 
				"#depth\tTauIn\tTauTot\tn(2p)/n(1s)\tTexc\tTe\tTex/T\tPesc\tPdes\tpump\topacity\talbedo\n" );
		}
		else
		{
			fprintf( ioQQQ, "Save hydrogen has options: CONDitions, 21 CM, LINE, POPUlations, and IONIzation.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("IONI") )
	{
		if( p.nMatch("RATE") )
		{
			/* save ionization rates, search for the name of an element */
			if( (nelem = p.GetElem() ) < 0 )
			{
				fprintf( ioQQQ, "There must be an element name on the ionization rates command.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			save.punarg[save.nsave][0] = (realnum)nelem;
			strcpy( save.chSave[save.nsave], "IONR" );
			sprintf( save.chHeader[save.nsave], 
				"#%s depth\teden\tdynamics.Rate\tabund\tTotIonize\tTotRecom\tSource\t ... \n",
				elementnames.chElementSym[nelem]);
		}
		else
		{
			/* save table giving ionization means */
			strcpy( save.chSave[save.nsave], "IONI" );
			sprintf( save.chHeader[save.nsave], 
				"#Mean ionization distribution\n" );
		}
	}

	else if( p.nMatch(" IP ") )
	{
		strcpy( save.chSave[save.nsave], " IP " );
		sprintf( save.chHeader[save.nsave], 
			"#ionization potentials, valence shell\n" );
	}

	else if( p.nMatch("LEID") )
	{
		if( p.nMatch( "LINE" ) )
		{
			/* save Leiden lines
			 * final intensities of the Leiden PDR models */
			strcpy( save.chSave[save.nsave], "LEIL" );
			sprintf( save.chHeader[save.nsave], "#ion\twl\tInt\trel int\n");
		}
		else
		{
			/* save Leiden structure
			* structure of the Leiden PDR models */
			strcpy( save.chSave[save.nsave], "LEIS" );
			sprintf( save.chHeader[save.nsave], 
				/* 1-17 */
				"#Leid  depth\tA_V(extentd)\tA_V(point)\tTe\tH0\tH2\tCo\tC+\tOo\tCO\tO2\tCH\tOH\te\tHe+\tH+\tH3+\t"
				/* 18 - 30 */
				"N(H0)\tN(H2)\tN(Co)\tN(C+)\tN(Oo)\tN(CO)\tN(O2)\tN(CH)\tN(OH)\tN(e)\tN(He+)\tN(H+)\tN(H3+)\t"
				/* 31 - 32 */
				"H2(Sol)\tH2(FrmGrn)\tH2(photodiss)\t"
				/* 33 - 46*/
				"G0(DB96)\trate(CO)\trate(C)\theat\tcool\tGrnP\tGr-Gas-Cool\tGr-Gas-Heat\tCOds\tH2dH\tH2vH\tChaT\tCR H\tMgI\tSI\t"
				"Si\tFe\tNa\tAl\tC\tC610\tC370\tC157\tC63\tC146\n" );
		}
	}


	// save results for NLTE series of plasma comparison meetings, 
	// specifically NLTE7 2011 Dec
	else if( p.nMatch("NLTE") )
	{
# ifdef USE_NLTE7
		strcpy( save.chSave[save.nsave], "NLTE" );
# else
		fprintf(ioQQQ," PROBLEM You must enable the USE_NLTE7 macro at compile-time to use this command.\n");
		fprintf(ioQQQ," To do so, add EXTRA=\"-DUSE_NLTE7\" to the end of the make command.\n");
		fprintf(ioQQQ," An example for a quad core machine:\n make -j 4 EXTRA=\"-DUSE_NLTE7\" \n");
		fprintf(ioQQQ," in the sys_XXX folder that you want to use.\n\n\n");
		cdEXIT(EXIT_FAILURE);
# endif
	}

	/* FE2NRG, FE2TP, and FE2COLL write the internal Fe II data into Stout format */
	else if (p.nMatch("FE2NRG"))
	{
		strcpy( save.chSave[save.nsave], "LY1" );
	}

	else if (p.nMatch("FE2TP"))
	{
		strcpy( save.chSave[save.nsave], "LY2" );
	}

	else if (p.nMatch("FE2COLL"))
	{
		strcpy( save.chSave[save.nsave], "LY3" );
	}

	else if( (p.nMatch("LINE") && p.nMatch("LIST")) || p.nMatch("LINELIST") )
	{
		/* save line list "output file" "Line List file" */
		strcpy( save.chSave[save.nsave], "LLST" );

		/* 
		 * we parsed off the second file name at start of this routine
		 * check if file was found, use it if it was, else abort
		 */
		if( !lgSecondFilename )
		{
			fprintf(ioQQQ , "There must be a second file name between "
				"double quotes on the SAVE LINE LIST command.  This second"
				" file contains the input line list.  I did not find it.\nSorry.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* actually get the lines, and malloc the space in the arrays 
		 * cdGetLineList will look on path */
		if( save.ipPnunit[save.nsave] == NULL )
		{
			/* make sure we free any allocated space from a previous call */
			save.SaveLineListFree(save.nsave);

			save.nLineList[save.nsave] = cdGetLineList(chSecondFilename,
								   save.chLineListLabel[save.nsave],
								   save.wlLineList[save.nsave]);

			if( save.nLineList[save.nsave] < 0 )
			{
				fprintf(ioQQQ,"DISASTER could not open SAVE LINE LIST file %s \n",
					chSecondFilename );
				cdEXIT(EXIT_FAILURE);
			}
		}

		// check whether intrinsic or emergent line emissivity
		save.lgEmergent[save.nsave] = false;
		if( p.nMatch("EMER") )
			save.lgEmergent[save.nsave] = true;

		// check whether cumulative or specific line emission
		save.lgCumulative[save.nsave] = false;
		if( p.nMatch("CUMU") )
			save.lgCumulative[save.nsave] = true;

		/* ratio option, in which pairs of lines form ratios, first over
		 * second */
		if( p.nMatch("RATI") )
		{
			save.lgLineListRatio[save.nsave] = true;
			if( save.nLineList[save.nsave]%2 )
			{
				/* odd number of lines - cannot take ratio */
				fprintf(ioQQQ , "There must be an even number of lines to"
					" take ratios of lines.  There were %li, an odd number."
					"\nSorry.\n", save.nLineList[save.nsave]);
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			/* no ratio */
			save.lgLineListRatio[save.nsave] = false;
		}

		/* keyword absolute says to do absolute rather than relative intensities 
		 * relative intensities are the default */
		if( p.nMatch("ABSO") )
		{
			save.punarg[save.nsave][0] = 1;
		}
		else
		{
			save.punarg[save.nsave][0] = 0;
		}

		// check whether column or row (default)
		if( p.nMatch("COLUMN") )
		{
			save.punarg[save.nsave][1] = 1;
		}
		else
		{
			save.punarg[save.nsave][1] = 0;
		}

		/* give header line */
		sprintf( save.chHeader[save.nsave], "#lineslist" );
		// do header now if reporting rows of lines
		if( !save.punarg[save.nsave][1] )
		{
			for( long int j=0; j<save.nLineList[save.nsave]; ++j )
			{
				/* if taking ratio then put div sign between pairs */
				if( save.lgLineListRatio[save.nsave] && is_odd(j) )
					strcat( save.chHeader[save.nsave] , "/" );
				else
					strcat( save.chHeader[save.nsave] , "\t" );
				sprintf( chTemp, "%s ", save.chLineListLabel[save.nsave][j] );
				strcat( save.chHeader[save.nsave], chTemp );
				sprt_wl( chTemp, save.wlLineList[save.nsave][j] );
				strcat( save.chHeader[save.nsave], chTemp );
			}
		}
		strcat( save.chHeader[save.nsave], "\n" );
	}

	else if( p.nMatch("LINE") && !p.nMatch("XSPE") && !p.nMatch("NEAR"))
	{
		/* save line options -
		 * this is not save xspec lines and not linear option
		 * check for keyword UNITS on line, then scan wavelength or energy units, 
		 * sets save.chConPunEnr*/
		ChkUnits(p);

		/* save line emissivity, line intensity, line array,
		 * and line data */
		if( p.nMatch("STRU") )
		{
			fprintf(ioQQQ," The	SAVE LINES STRUCTURE command is now SAVE LINES "
				"EMISSIVITY.\n Sorry.\n\n");
			cdEXIT(EXIT_FAILURE);
		}

		else if( p.nMatch("PRES") )
		{
			/* save contributors to line pressure */
			strcpy( save.chSave[save.nsave], "PREL" );
			sprintf( save.chHeader[save.nsave], 
				"#P depth\tPtot\tPline/Ptot\tcontributors to line pressure\n" );
		}

		else if( p.nMatch("EMIS") )
		{
			/* this used to be the save lines structure command, is now
			 * the save lines emissivity command 
			 * give line emissivity vs depth */
			// check whether intrinsic or emergent line emissivity
			save.lgEmergent[save.nsave] = false;
			if( p.nMatch("EMER") )
				save.lgEmergent[save.nsave] = true;
			strcpy( save.chSave[save.nsave], "LINS" );
			sprintf( save.chHeader[save.nsave],
				"#");
			/* read in the list of lines to examine */
			parse_save_line(p,false, chTemp );
			strcat( save.chHeader[save.nsave], chTemp );
		}

		else if( p.nMatch(" RT " ) )
		{
			/* save line RT */
			strcpy( save.chSave[save.nsave], "LINR" );
			/* save some details needed for line radiative transfer 
			 * routine in save_line.cpp */
			Parse_Save_Line_RT(p);
		}

		else if( p.nMatch("CUMU") )
		{
			bool lgEOL;
			/* save lines cumulative 
			 * this will be integrated line intensity, function of depth */
			strcpy( save.chSave[save.nsave], "LINC" );
			// option for intrinsic (default) or emergent
			save.lgEmergent[save.nsave] = false;
			if( p.nMatch("EMER") )
				save.lgEmergent[save.nsave] = true;
			/* option for either relative intensity or abs luminosity */
			if( p.nMatch("RELA") )
			{
				lgEOL = true;
				sprintf( save.chHeader[save.nsave], "#" );
			}
			else
			{
				sprintf( save.chHeader[save.nsave], "#" );
				lgEOL = false;
			}
			/* read in the list of lines to examine */
			parse_save_line(p, lgEOL, chTemp );
			strcat( save.chHeader[save.nsave], chTemp );
		}

		else if( p.nMatch("DATA") )
		{
			/* save line data, done in SaveLineData */

			/* the default will be to make wavelengths like in the printout, called labels,
			 * if units appears then other units will be used instead */
			save.chConPunEnr[save.nsave] = "labl";

			/* check for keyword UNITS on line, then scan wavelength or energy units if present,
			 * units are copied into save.chConPunEnr */
			if( p.nMatch("UNIT") )
				ChkUnits(p);
			strcpy( save.chSave[save.nsave], "LIND" );
			sprintf( save.chHeader[save.nsave], 
				"#Emission line data.\n" );
		}

		else if( p.nMatch("ARRA") )
		{
			/* save line array -
			 * output energies and luminosities of predicted lines */
			strcpy( save.chSave[save.nsave], "LINA" );
			sprintf( save.chHeader[save.nsave], 
				"#enr\tID\tI(intrinsic)\tI(emergent)\ttype\n" );
		}

		else if( p.nMatch("LABE") )
		{
			/* save line labels */
			strcpy( save.chSave[save.nsave], "LINL" );
			sprintf( save.chHeader[save.nsave], 
				"#index\tlabel\twavelength\tcomment\n" );
			/* this controls whether we will print lots of redundant 
			 * info labels for transferred lines - if keyword LONG appears
			 * then do so, if does not appear then do not - this is default */
			if( p.nMatch("LONG") )
				save.punarg[save.nsave][0] = 1;
			else
				save.punarg[save.nsave][0] = 0;
		}

		else if( p.nMatch("OPTI") )
		{
			/* save line optical depths, done in SaveLineStuff */
			strcpy( save.chSave[save.nsave], "LINO" );

			/* the default will be to make wavelengths line in the printout, called labels,
			 * if units appears then other units will be used instead */
			save.chConPunEnr[save.nsave] = "labl";

			/* check for keyword UNITS on line, then scan wavelength or energy units if present,
			 * units are copied into save.chConPunEnr */
			if( p.nMatch("UNIT") )
				ChkUnits(p);

			sprintf( save.chHeader[save.nsave], 
				"#species\tenergy/%s\topt depth\tdamp\n",
				save.chConPunEnr[save.nsave] );

			/* this is optional limit to smallest optical depths */
			save.punarg[save.nsave][0] = (realnum)pow(10.,p.FFmtRead());
			/* this is default of 0.1 napier */
			if( p.lgEOL() )
			{
				save.punarg[save.nsave][0] = 0.1f;
			}
		}

		else if( p.nMatch("POPU") )
		{
			/* save line populations command - first give index and inforamtion
			 * for all lines, then populations for lines as a function of
			 * depth, using this index */
			strcpy( save.chSave[save.nsave], "LINP" );
			sprintf( save.chHeader[save.nsave], 
				"#population information\n" );
			/* this is optional limit to smallest population to save - always
			 * interpreted as a log */
			save.punarg[save.nsave][0] = (realnum)pow(10.,p.FFmtRead());

			/* this is default - all positive populations */
			if( p.lgEOL() )
				save.punarg[save.nsave][0] = 0.f;

			if( p.nMatch(" OFF") )
			{
				/* no lower limit - print all lines */
				save.punarg[save.nsave][0] = -1.f;
			}
		}

		else if( p.nMatch("INTE") )
		{
			/* this will be full set of line intensities */
			strcpy( save.chSave[save.nsave], "LINI" );
			sprintf( save.chHeader[save.nsave], 
				"#Emission line intrinsic intensities per unit inner area\n" );
			if( p.nMatch("COLU") )
				/* column is key to save single column */
				strcpy( save.chPunRltType, "column" );
			else
				/* array is key to save large array */
				strcpy( save.chPunRltType, "array " );

			save.punarg[save.nsave][0] = 0.;
			// ALL option - all lines, even zero intensities
			if( p.nMatch( " ALL" ) )
				save.punarg[save.nsave][0] = -1.;

			// check whether intrinsic or emergent line emissivity
			save.lgEmergent[save.nsave] = false;
			if( p.nMatch("EMER") )
				save.lgEmergent[save.nsave] = true;

			if( p.nMatch("EVER") )
			{
				save.LinEvery = (long int)p.FFmtRead();
				save.lgLinEvery = true;
				if( p.lgEOL() )
				{
					fprintf( ioQQQ, 
						"There must be a second number, the number of zones to print.\nSorry.\n" );
					cdEXIT(EXIT_FAILURE);
				}
			}
			else
			{
				save.LinEvery = geometry.nend[0];
				save.lgLinEvery = false;
			}
		}
		else
		{
			fprintf( ioQQQ, 
				"This option for SAVE LINE is something that I do not understand.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch(" MAP") )
	{
		strcpy( save.chSave[save.nsave], "MAP " );
		sprintf( save.chHeader[save.nsave], 
			"#te, heating, cooling.\n" );
		/* do cooling space map for specified zones
		 * if no number, or <0, do map and save out without doing first zone
		 * does map by calling punt(" map") 
		 */
		hcmap.MapZone = (long)p.FFmtRead();
		if( p.lgEOL() )
		{
			hcmap.MapZone = 1;
		}

		if( p.nMatch("RANG") )
		{
			bool lgLogOn;
			hcmap.RangeMap[0] = (realnum)p.FFmtRead();
			if( hcmap.RangeMap[0] <= 10. && !p.nMatch("LINE") )
			{
				hcmap.RangeMap[0] = (realnum)pow((realnum)10.f,hcmap.RangeMap[0]);
				lgLogOn = true;
			}
			else
			{
				lgLogOn = false;
			}

			hcmap.RangeMap[1] = (realnum)p.FFmtRead();
			if( lgLogOn )
				hcmap.RangeMap[1] = (realnum)pow((realnum)10.f,hcmap.RangeMap[1]);

			if( p.lgEOL() )
			{
				fprintf( ioQQQ, "There must be a zone number, followed by two temperatures, on this line.  Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}

	else if( p.nMatch("MOLE") )
	{
		/* molecules, especially for PDR calculations */
		strcpy( save.chSave[save.nsave], "MOLE" );
	}

	else if( p.nMatch("MONI") )
	{
		/* save monitors */
		strcpy( save.chSave[save.nsave], "MONI" );
	}

	else if( p.nMatch("OPTICAL") && p.nMatch("DEPTH") )
	{
		/* check for keyword UNITS on line, then scan wavelength or energy units if present,
		 * units are copied into save.chConPunEnr */
		ChkUnits(p);

		/* "every" option to save this on every zone -
		 * not present then only last zone is saved */
		if( p.nMatch("EVER" ) )
		{
			/* save every zone */
			save.lgSaveEveryZone[save.nsave] = true;
			save.nSaveEveryZone[save.nsave] = 1;
		}
		else
		{
			/* only save last zone */
			save.lgSaveEveryZone[save.nsave] = false;
			save.nSaveEveryZone[save.nsave] = 1;
		}
		
		if( p.nMatch("FINE") )
		{
			/* save fine continuum optical depths */
			rfield.lgSaveOpacityFine = true;
			strcpy( save.chSave[save.nsave], "OPTf" );
			sprintf( save.chHeader[save.nsave], "#energy/%s\tTau tot\topacity\n",
				save.chConPunEnr[save.nsave] );
			/* range option - important since so much data */
			if( p.nMatch("RANGE") ) 
			{
				/* get lower and upper range, eventually must be in Ryd */
				double Energy1 = p.FFmtRead();
				double Energy2 = p.FFmtRead();
				if( p.lgEOL() )
				{
					fprintf(ioQQQ,"There must be two numbers, the lower and upper energy range in Ryd.\nSorry.\n");
					cdEXIT(EXIT_FAILURE);
				}
				if( p.nMatch("UNIT" ) )
				{
					// apply units to range option
					const char *energyUnits = p.StandardEnergyUnit();
					Energy unitChange;
					unitChange.set(Energy1, energyUnits );
					Energy1 = unitChange.Ryd();
					unitChange.set(Energy2, energyUnits );
					Energy2 = unitChange.Ryd();
				}
				/* get lower and upper rang in Ryd */
				save.punarg[save.nsave][0] = (realnum)MIN2( Energy1 , Energy2 );
				save.punarg[save.nsave][1] = (realnum)MAX2( Energy1 , Energy2 );
				//fprintf(ioQQQ , "DEBUG units change fine %.3e %.3e\n" , save.punarg[save.nsave][0] ,
				//		save.punarg[save.nsave][1] );
				//cdEXIT(EXIT_FAILURE);
			}
			else
			{
				/* these mean full energy range */
				save.punarg[save.nsave][0] = 0.;
				save.punarg[save.nsave][1] = 0.;
			}
			/* optional last parameter - how many points to bring together */
			save.punarg[save.nsave][2] = (realnum)p.FFmtRead();
			/* default is to bring together ten */
			if( p.lgEOL() )
				save.punarg[save.nsave][2] = 10;
			if( save.punarg[save.nsave][2] < 1 )
			{
				fprintf(ioQQQ,"The number of fine opacities to skip must be > 0 \nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			/* save coarse continuum optical depths */
			strcpy( save.chSave[save.nsave], "OPTc" );
			sprintf( save.chHeader[save.nsave], 
				"#energy/%s\ttotal\tabsorp\tscat\n",
				save.chConPunEnr[save.nsave] );
		}

	}
	else if( p.nMatch(" OTS") )
	{
		strcpy( save.chSave[save.nsave], " OTS" );
		sprintf( save.chHeader[save.nsave], 
			"#otscon, lin, conOpac LinOpc\n" );
	}

	else if( p.nMatch("OVER") && p.nMatch(" OVE") )
	{
		/* save overview of model results */
		strcpy( save.chSave[save.nsave], "OVER" );
		sprintf( save.chHeader[save.nsave], 
			"#depth\tTe\tHtot\thden\teden\t2H_2/H\tHI\tHII\tHeI\tHeII\tHeIII\tCO/C\tC1\tC2\tC3\tC4\tO1\tO2\tO3\tO4\tO5\tO6\tH2O/O\tAV(point)\tAV(extend)\n" );
	}

	else if( p.nMatch(" PDR") )
	{
		strcpy( save.chSave[save.nsave], " PDR" );
		sprintf( save.chHeader[save.nsave], 
			"#depth\tH colden\tTe\tHI/HDEN\tH2/HDEN\tH2*/HDEN\tCI/C\tCO/C\tH2O/O\tG0\tAV(point)\tAV(extend)\tTauV(point)\n" );
	}

	else if( p.nMatch("PERF") )
	{
		/* output performance characteristics per zone */
		strcpy( save.chSave[save.nsave], "PERF" );
		sprintf( save.chHeader[save.nsave],
			"#zone\tdTime\tElapsed t\tnPres2Ioniz\n" );
	}

	else if( p.nMatch("PHYS") )
	{
		/* save physical conditions */
		strcpy( save.chSave[save.nsave], "PHYS" );
		sprintf( save.chHeader[save.nsave], 
			"#PhyC depth\tTe\tn(H)\tn(e)\tHtot\taccel\tfillfac\n" );
	}

	else if( p.nMatch("POIN") )
	{
		/* save out the pointers */
		save.lgPunPoint = true;
		/* this does not count as a save option (really) */
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("PRES") )
	{
		/* the save pressure command */
		strcpy( save.chSave[save.nsave], "PRES" );
		sprintf( save.chHeader[save.nsave], 
			"#P depth\tPerror%%\tPcurrent\tPIn+Pinteg\tPgas(r0)\tPgas\tPram"
			"\tPrad(line)\tPinteg\tV(wind km/s)\tcad(wind km/s)\tP(mag)\tV(turb km/s)"
			"\tP(turb)\tPgr_Int\tint thin elec\tconv?\n" );
	}

	else if( p.nMatch("RADI") )
	{
		/* the save radius command */
		sprintf( save.chHeader[save.nsave], "#NZONE\tradius\tdepth\tdr\n" );
		/* option to only save the outer radius */
		if( p.nMatch( "OUTE" ) )
		{
			/* only outer radius */
			strcpy( save.chSave[save.nsave], "RADO" );
		}
		else
		{
			/* all radii */
			strcpy( save.chSave[save.nsave], "RADI" );
		}
	}

	else if( p.nMatch("RECO") )
	{
		if( p.nMatch("COEF") )
		{
			/* recombination coefficients for everything */

			/* this is logical flag used in routine ion_recom to create the save output */
			save.lgioRecom = true;
			/* this does not count as a save option (really) */
			strcpy( save.chSave[save.nsave], "" );
			save.lgRealSave[save.nsave] = false;
		}

		else if( p.nMatch("EFFI") )
		{
			/* save recombination efficiency */
			strcpy( save.chSave[save.nsave], "RECE" );
			sprintf( save.chHeader[save.nsave], 
				"#Recom effic H, Heo, He+\n" );
		}

		else
		{
			fprintf( ioQQQ, "No option recognized on this save recombination command\n" );
			fprintf( ioQQQ, "Valid options are COEFFICIENTS, AGN, and EFFICIENCY\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* save results command, either as single column or wide array */
	else if( p.nMatch("RESU") )
	{
		strcpy( save.chSave[save.nsave], "RESU" );
		if( p.nMatch("COLU") )
		{
			/* column is key to save single column */
			strcpy( save.chPunRltType, "column" );
		}
		else
		{
			/* array is key to save large array */
			strcpy( save.chPunRltType, "array " );
		}

		/* do not change following, is used as flag in getlines */
		sprintf( save.chHeader[save.nsave], 
			"#results of calculation\n" );
	}

	else if( p.nMatch("SECO") )
	{
		/* save secondary ionization rate */
		strcpy( save.chSave[save.nsave], "SECO" );
		sprintf( save.chHeader[save.nsave], 
			"#depth\tIon(H^0)\tDiss(H_2)\tExcit(Lya)\n" );
	}

	else if( p.nMatch("SOUR") )
	{

		/* check for keyword UNITS on line, then scan wavelength or energy units if present,
		* units are copied into save.chConPunEnr */
		ChkUnits(p);

		if( p.nMatch("DEPT") )
		{
			/* print continuum source function as function of depth */
			strcpy( save.chSave[save.nsave], "SOUD" );
			sprintf( save.chHeader[save.nsave], 
				"#continuum source function vs depth\n" );
		}
		else if( p.nMatch("SPEC") )
		{
			/* print spectrum continuum source function at 1 depth */
			strcpy( save.chSave[save.nsave], "SOUS" );
			sprintf( save.chHeader[save.nsave], 
				"#continuum source function nu/%s\tConEmitLocal/widflx"
				"\tabs opac\tConSourceFcnLocal\tConSourceFcnLocal/plankf\tConSourceFcnLocal/flux\n",
				save.chConPunEnr[save.nsave] );
		}
		else
		{
			fprintf( ioQQQ, "A second keyword must appear on this line.\n" );
			fprintf( ioQQQ, "They are DEPTH and SPECTRUM.\n" );
			fprintf( ioQQQ, "Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}


	/* save spectrum the new form of the save continuum, will eventually replace the standard
	 * save continuum command */
	else if( p.nMatch("SPECTRUM") && !p.nMatch("XSPE") )
	{
		/* this flag is checked in PrtComment to generate a caution
		 * if continuum is saved but iterations not performed */
		save.lgPunContinuum = true;

		/* set flag for spectrum */
		strcpy( save.chSave[save.nsave], "CONN" );

		/* check for keyword UNITS on line, then scan wavelength or energy units if present,
		 * units are copied into save.chConPunEnr */
		ChkUnits(p);

		sprintf( save.chHeader[save.nsave], 
			"#Cont Enr/%s\tincid nFn\ttrans\tdiff\tlines \n",
				save.chConPunEnr[save.nsave] );
	}

	else if( p.nMatch("SPECIAL") )
	{
		/* save special, will call routine SaveSpecial */
		strcpy( save.chSave[save.nsave], "SPEC" );
		sprintf( save.chHeader[save.nsave], "#Special.\n" );
	}

	else if( p.nMatch("SPECIES") )
	{
		strcpy( save.chSave[save.nsave], "SPCS" );

		// option to save information about a particular species,
		// the "second filename" may really be the species label.  Rename here for clarity
		strcpy( chLabel, chSecondFilename );
		if( !lgSecondFilename )
		{
			strcpy( save.chSaveSpecies[save.nsave], "" );
		}
		else if( strlen(chLabel) >= CHARS_SPECIES )
		{
			fprintf( ioQQQ,"Species string is limited to %li characters.\nSorry.\n", (long)CHARS_SPECIES );
			cdEXIT(EXIT_FAILURE);
		}
		else
			strncpy( save.chSaveSpecies[save.nsave], chLabel, CHARS_SPECIES );

		if (p.nMatch( "COLUMN" ) )
		{
			/* column densities*/
			strcpy( save.chSaveArgs[save.nsave], "COLU" );
		}
		else if (p.nMatch( "ENERG" ) )
		{
			/* energy levels, default Rydbergs but option to change units */
			ChkUnits(p);
			strcpy( save.chSaveArgs[save.nsave], "ENER" );
		}
		else if( p.nMatch("LABELS") )
		{
			strcpy( save.chSaveArgs[save.nsave], "LABE" );
		}
		else if (p.nMatch( "LEVELS" ) )
		{
			/* the number of levels in this zone */
			strcpy( save.chSaveArgs[save.nsave], "LEVL" );
		}
		else if (p.nMatch( "POPUL" ) )
		{
			/* save species population fraction for 1 level*/
			strcpy( save.chSaveArgs[save.nsave], "POPU" );
		}
		else
		{
			fprintf( ioQQQ, "ParseSave cannot find a recognized keyword on this SAVE SPECIES command line.\n" );
			fprintf( ioQQQ, "I know about the keywords COLUMN DENSITIES, LABELS, LEVELS, and POPULATIONS.\nSorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	else if( p.nMatch("TEMP") )
	{
		/* save temperature command */
		strcpy( save.chSave[save.nsave], "TEMP" );
		sprintf( save.chHeader[save.nsave], 
			"#depth\tTe\tcC/dT\tdt/dr\td^2T/dr^2\n" );
	}

	else if( p.nMatch("TIME") && p.nMatch("DEPE") )
	{
		/* information about time dependent solutions */
		strcpy( save.chSave[save.nsave], "TIMD" );
		/* do not want to separate iterations with special character */
		save.lg_separate_iterations[save.nsave] = false;
		/* write header */
		sprintf( save.chHeader[save.nsave] , 
			"#elapsed time\ttime step \tscale cont\tn(H)\t<T>\t<H+/H rad>\t<H0/H rad>\t<H2/H rad>\t<He+/He rad>\t<CO/H>\t<redshift>\t<ne/nH>\n" );
	}

	else if( p.nMatch("TPRE") )
	{
		/* debug output from the temperature predictor in zonestart, 
		 * set with save tpred command */
		strcpy( save.chSave[save.nsave], "TPRE" );
		sprintf( save.chHeader[save.nsave], 
			"#zone  old temp,  guess Tnew, new temp    delta \n" );
	}

	else if( p.nMatch("WIND") )
	{
		strcpy( save.chSave[save.nsave], "WIND" );
		sprintf( save.chHeader[save.nsave], 
			"#radius\tdepth\tvel [cm/s]\tTot accel [cm s-2]\tLin accel [cm s-2]"
			"\tCon accel [cm s-2]\tforce multiplier\ta_gravity\n" );
		if( p.nMatch( "TERM" ) )
		{
			/* only save for last zone, the terminal velocity, for grids */
			save.punarg[save.nsave][0] = 0.;
		}
		else
		{
			/* one means save every zone */
			save.punarg[save.nsave][0] = 1.;
		}
	}

	else if( p.nMatch("XSPE") )
	{
		/* say that this is a FITS file output */
		save.lgFITS[save.nsave] = true;

		/* the save xspec commands */
		save.lgPunLstIter[save.nsave] = true;

		/* remember that a save xspec command was entered */
		grid.lgSaveXspec = true;

		/* range option - important since so much data */
		if( p.nMatch("RANGE") ) 
		{
			/* get lower and upper range, must be in keV */
			save.punarg[save.nsave][0] = (realnum)p.FFmtRead();
			save.punarg[save.nsave][1] = (realnum)p.FFmtRead();
			if( p.lgEOL() )
			{
				fprintf(ioQQQ,"There must be two numbers, the lower and upper energy range in keV.\nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}
			if( save.punarg[save.nsave][0] >=save.punarg[save.nsave][1] )
			{
				fprintf(ioQQQ,"The two energies for the range must be in increasing order.\nSorry.\n");
				cdEXIT(EXIT_FAILURE);
			}

			grid.LoEnergy_keV = save.punarg[save.nsave][0];
			grid.HiEnergy_keV = save.punarg[save.nsave][1];
		}
		else
		{
			/* these mean full energy range */
			save.punarg[save.nsave][0] = 0;
			save.punarg[save.nsave][1] = 0;
		}

		if( p.nMatch("ATAB") )
		{
			/* save xspec atable command */
			
			if( p.nMatch("TOTA") )	
			{
				/* total spectrum */
				strcpy( save.chSave[save.nsave], "XTOT" );
				grid.lgOutputTypeOn[0] = true;
				save.FITStype[save.nsave] = 0;
			}
			else if( p.nMatch("INCI") )
			{
				if( p.nMatch("ATTE") )
				{
					/* attenuated incident continuum */
					strcpy( save.chSave[save.nsave], "XATT" );
					grid.lgOutputTypeOn[2] = true;
					save.FITStype[save.nsave] = 2;
				}
				else if( p.nMatch("REFL") )
				{
					/* reflected incident continuum */
					strcpy( save.chSave[save.nsave], "XRFI" );
					grid.lgOutputTypeOn[3] = true;
					save.FITStype[save.nsave] = 3;
				}
				else
				{
					/* incident continuum */
					strcpy( save.chSave[save.nsave], "XINC" );
					grid.lgOutputTypeOn[1] = true;
					save.FITStype[save.nsave] = 1;
				}
			}
			else if( p.nMatch("DIFF") )
			{
				if( p.nMatch("REFL") )
				{
					/* reflected diffuse continuous emission */
					strcpy( save.chSave[save.nsave], "XDFR" );
					grid.lgOutputTypeOn[5] = true;
					save.FITStype[save.nsave] = 5;
				}
				else
				{
					/* diffuse continuous emission outward */
					strcpy( save.chSave[save.nsave], "XDFO" );
					grid.lgOutputTypeOn[4] = true;
					save.FITStype[save.nsave] = 4;
				}
			}
			else if( p.nMatch("LINE") )
			{
				if( p.nMatch("REFL") )
				{
					/* reflected lines */
					strcpy( save.chSave[save.nsave], "XLNR" );
					grid.lgOutputTypeOn[7] = true;
					save.FITStype[save.nsave] = 7;
				}
				else
				{
					/* outward lines */
					strcpy( save.chSave[save.nsave], "XLNO" );
					grid.lgOutputTypeOn[6] = true;
					save.FITStype[save.nsave] = 6;
				}
			}
			else if( p.nMatch("SPEC") )
			{
				if( p.nMatch("REFL") )
				{
					/* reflected spectrum */
					strcpy( save.chSave[save.nsave], "XREF" );
					grid.lgOutputTypeOn[9] = true;
					save.FITStype[save.nsave] = 9;
				}
				else
				{
					/* transmitted spectrum */
					strcpy( save.chSave[save.nsave], "XTRN" );
					grid.lgOutputTypeOn[8] = true;
					save.FITStype[save.nsave] = 8;
				}
			}
			else
			{
				/* transmitted spectrum */
				strcpy( save.chSave[save.nsave], "XTRN" );
				grid.lgOutputTypeOn[8] = true;
				save.FITStype[save.nsave] = 8;
			}
		}
		else if( p.nMatch("MTAB") )
		{
			/* save xspec mtable */
			strcpy( save.chSave[save.nsave], "XSPM" );
			grid.lgOutputTypeOn[10] = true;
			save.FITStype[save.nsave] = 10;
		}
		else
		{
			fprintf( ioQQQ, "Support only for xspec atable and xspec mtable.\n" );
			cdEXIT( EXIT_FAILURE );
		}
	}

	/* save column density has to come last so do not trigger specific column 
	 * densities, H2, FeII, etc.
	 * Need both keywords since column is also the keyword for one line per line */
	else if( p.nMatch("COLU") && p.nMatch("DENS") )
	{
		if( p.nMatch("SOME" ))
		{
			/* flag saying save some column densities */
			strcpy( save.chSave[save.nsave], "COLS" );
			parse_save_colden( p, save.chHeader[save.nsave] );
		}
		else
		{
			/* save column densities table */
			strcpy( save.chSave[save.nsave], "COLU" );
		}
	}
	else
	{
		fprintf( ioQQQ, 
			"ParseSave cannot find a recognized keyword on this SAVE command line.\nSorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* only open if file has not already been opened during a previous call */
	if( save.ipPnunit[save.nsave] == NULL )
	{
		string file_name;
		file_name += chFilename;
		string mode = "w";
		if( save.lgFITS[save.nsave] )
			mode += "b";

		/* open the file with the name and mode generated above */
		save.ipPnunit[save.nsave] = open_data( file_name.c_str(), mode.c_str(), AS_LOCAL_ONLY );
	
		/* option to set no buffering for this file.  The setbuf command may 
		 * ONLY be issued right after the open of the file.  Giving it after 
		 * i/o has been done may result in loss of the contents of the buffer, PvH */
		if( p.nMatch("NO BUFFER") )
			setbuf( save.ipPnunit[save.nsave] , NULL );
	}

	/***************************************************************/
	/*                                                             */
	/*  The following are special save options and must be done   */
	/*  after the parsing and file opening above.                  */
	/*                                                             */
	/*  NB:  these are ALSO parsed above.  Here we DO something.   */
	/*                                                             */
	/***************************************************************/

	if( p.nMatch("CONV") && p.nMatch("REAS") )
	{
		/* save reason model declared not converged
		 * not a true save command, since done elsewhere */
		save.ipPunConv = save.ipPnunit[save.nsave];
		save.lgPunConv_noclobber = save.lgNoClobber[save.nsave];
		save.lgPunConv = true;
		fprintf( save.ipPunConv, 
			"# reason for continued iterations\n" );
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("CONV") && p.nMatch("BASE") )
	{
		/* save some quantities we are converging */
		save.lgTraceConvergeBase = true;
		/* the second save occurrence - file has been opened,
		* copy handle, also pass on special no hash option */
		if( p.nMatch("NO HA") )
			save.lgTraceConvergeBaseHash = false;
		save.ipTraceConvergeBase = save.ipPnunit[save.nsave];
		/* set save last flag to whatever it was above */
		save.lgTraceConvergeBase_noclobber = save.lgNoClobber[save.nsave];
		static bool lgPrtHeader = true;
		if( lgPrtHeader )
			fprintf( save.ipTraceConvergeBase, 
			"#zone\theat\tcool\teden\n" );
		lgPrtHeader = false;
	}

	else if( p.nMatch(" DR ") )
	{
		static bool lgPrtHeader = true;
		/* the second save dr occurrence - file has been opened,
		 * copy handle to ipDRout, also pass on special no hash option */
		if( p.nMatch("NO HA") )
			save.lgDRHash = false;
		save.ipDRout = save.ipPnunit[save.nsave];
		/* set save last flag to whatever it was above */
		save.lgDRPLst = save.lgPunLstIter[save.nsave];
		save.lgDROn_noclobber = save.lgNoClobber[save.nsave];
		if( lgPrtHeader )
			fprintf( save.ipDRout, 
			"#zone\tdepth\tdr\tdr 2 go\treason \n" );
		lgPrtHeader = false;
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("QHEA") )
	{
		gv.QHSaveFile = save.ipPnunit[save.nsave];
		gv.lgQHPunLast = save.lgPunLstIter[save.nsave];
		save.lgQHSaveFile_noclobber = save.lgNoClobber[save.nsave];
		fprintf( gv.QHSaveFile, 
			"#Probability distributions from quantum heating routine.\n" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("POIN") )
	{
		/* save out the pointers */
		save.ipPoint = save.ipPnunit[save.nsave];
		save.lgPunPoint_noclobber = save.lgNoClobber[save.nsave];
		save.lgPunPoint = true;
		fprintf( save.ipPoint, 
			"#pointers. \n" );
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("RECO") && p.nMatch("COEF") )
	{
		/* recombination coefficients for everything
		 * save.lgioRecom set to false in routine zero, non-zero value
		 * is flag to save recombination coefficients. the output is actually
		 * produced by a series of routines, as they generate the recombination
		 * coefficients.  these include 
		 * diel supres, helium, hydrorecom, iibod, and makerecomb*/
		save.ioRecom = save.ipPnunit[save.nsave];
		save.lgioRecom_noclobber = save.lgNoClobber[save.nsave];
		/* this is logical flag used in routine ion_recom to create the save output */
		save.lgioRecom = true;
		fprintf( save.ioRecom, 
			"#recombination coefficients cm3 s-1 for current density and temperature\n" );
		strcpy( save.chSave[save.nsave], "" );
		save.lgRealSave[save.nsave] = false;
	}

	else if( p.nMatch("GRID") )
	{
		/* this enables saving GRID output outside the main SaveDo() loop */
		grid.pnunit = save.ipPnunit[save.nsave];
		save.lgSaveGrid_noclobber = save.lgNoClobber[save.nsave];
	}

	else if( p.nMatch(" MAP") )
	{
		/* say output goes to special save */
		ioMAP = save.ipPnunit[save.nsave];
	}

	/* check that string written into save.chHeader[save.nsave] can actually fit there
	 * we may have overrun this buffer, an internal error */
	/* check that there are less than nChar characters in the line */
	char *chEOL = strchr_s(save.chHeader[save.nsave] , '\0' );

	/* return null if input string longer than nChar, the longest we can read.
	* Print and return null but chLine still has as much of the line as 
	* could be placed in cdLine */
	if( (chEOL==NULL) || (chEOL - save.chHeader[save.nsave])>=MAX_HEADER_SIZE-1 )
	{
		fprintf( ioQQQ, "DISASTER save.chHeader[%li] has been overwritten "
			"with a line too long to be read.\n", save.nsave );
		cdEXIT(EXIT_FAILURE);
	}

	/* if lgPunHeader true and cdHeader has been set to a string then print header
	 * logic to prevent more than one header in grid calculation */
	if( save.lgPunHeader[save.nsave] && !nMatch(save.chHeader[save.nsave],save.chNONSENSE) )
	{
		fprintf( save.ipPnunit[save.nsave], "%s", save.chHeader[save.nsave] );
		save.lgPunHeader[save.nsave] = false;
	}

	/* increment total number of save commands, */
	++save.nsave;
	return;
}

/*SaveFilesInit initialize save file pointers, called from InitCoreload
 * called one time per core load 
 * NB KEEP THIS ROUTINE SYNCHED UP WITH THE NEXT ONE, CloseSaveFiles */
void SaveFilesInit()
{
	long int i;
	static bool lgFIRST = true;

	DEBUG_ENTRY( "SaveFilesInit()" );

	ASSERT( lgFIRST );
	lgFIRST = false;

	/* set lgNoClobber to not overwrite files, reset with clobber on save line
	 * if we are running a grid (grid command entered in cdRead) grid.lgGrid 
	 * true, is false if single sim.  For grid we want to not clobber files 
	 * by default, do clobber for optimizer since this was behavior before */
	bool lgNoClobberDefault = false;
	if( grid.lgGrid )
	{
		/* cdRead encountered grid command - do not want to clobber files */
		lgNoClobberDefault = true;
	}

	for( i=0; i < LIMPUN; i++ )
	{
		save.lgNoClobber[i] = lgNoClobberDefault;
	}
	save.lgPunConv_noclobber = lgNoClobberDefault;
	save.lgDROn_noclobber = lgNoClobberDefault;
	save.lgTraceConvergeBase_noclobber = lgNoClobberDefault;
	save.lgPunPoint_noclobber = lgNoClobberDefault;
	save.lgioRecom_noclobber = lgNoClobberDefault;
	save.lgQHSaveFile_noclobber = lgNoClobberDefault;
	save.lgSaveGrid_noclobber = lgNoClobberDefault;

	/* initialize chHeader strings with nonsense, compare later to see if we have any actual headers. */
	save.chNONSENSE = "ArNdY38dZ9us4N4e12SEcuQ";

	for( i=0; i < LIMPUN; i++ )
	{
		save.ipPnunit[i] = NULL;

		// is this a real save command?  set false with the dummy 
		// save commands like save dr
		save.lgRealSave[i] = true;

		// do we need to save header?
		save.lgPunHeader[i] = true;
		strcpy( save.chHeader[i], save.chNONSENSE );
	}

	save.lgTraceConvergeBase = false;

	save.ipDRout = NULL;
	save.lgDROn = false;

	save.ipTraceConvergeBase = NULL;
	save.lgTraceConvergeBase = false;

	save.ipPunConv = NULL;
	save.lgPunConv = false;

	save.ipPoint = NULL;
	save.lgPunPoint = false;

	gv.QHSaveFile = NULL;

	save.ioRecom = NULL;
	save.lgioRecom = false;

	grid.pnunit = NULL;

	ioMAP = NULL;

	return;
}

/*CloseSaveFiles close save files called from cdEXIT upon termination,
 * from cloudy before returning 
 * NB - KEEP THIS ROUTINE SYNCHED UP WITH THE PREVIOUS ONE, SaveFilesInit */
void CloseSaveFiles( bool lgFinal )
{
	long int i;

	DEBUG_ENTRY( "CloseSaveFiles()" );

	/* close all save units cloudy opened with save command,
	 * lgNoClobber is set false with CLOBBER option on save, says to
	 * overwrite the files */
	for( i=0; i < save.nsave; i++ )
	{
		/* if lgFinal is true, we close everything, no matter what. 
		 * this means ignoring "no clobber" options */
		if( save.ipPnunit[i] != NULL && ( !save.lgNoClobber[i] || lgFinal ) )
		{
			/* Test that any FITS files are the right size! */ 
			if( save.lgFITS[i] )
			{
				/* \todo 2 This overflows for file sizes larger (in bytes) than
				 * a long int can represent (about 2GB on most 2007 systems)  */
				fseek(save.ipPnunit[i], 0, SEEK_END);
				long file_size = ftell(save.ipPnunit[i]);
				if( file_size%2880 )
				{
					fprintf( ioQQQ, " PROBLEM  FITS file is wrong size!\n" );
				}
			}

			fclose( save.ipPnunit[i] );
			save.ipPnunit[i] = NULL;
		}
	}

	/* following file handles are aliased to ipPnunit which was already closed above */
	if( save.ipDRout != NULL && ( !save.lgDROn_noclobber || lgFinal ) )
	{
		save.ipDRout = NULL;
		save.lgDROn = false;
	}

	if( save.ipTraceConvergeBase != NULL && ( !save.lgTraceConvergeBase_noclobber || lgFinal ) )
	{
		save.ipTraceConvergeBase = NULL;
		save.lgTraceConvergeBase = false;
	}

	if( save.ipPunConv != NULL && ( !save.lgPunConv_noclobber || lgFinal ) )
	{
		save.ipPunConv = NULL;
		save.lgPunConv = false;
	}
	if( save.ipPoint != NULL && ( !save.lgPunPoint_noclobber || lgFinal ) )
	{
		save.ipPoint = NULL;
		save.lgPunPoint = false;
	}
	if( gv.QHSaveFile != NULL  && ( !save.lgQHSaveFile_noclobber || lgFinal ) )
	{
		gv.QHSaveFile = NULL;
	}
	if( save.ioRecom != NULL  && ( !save.lgioRecom_noclobber || lgFinal ) )
	{
		save.ioRecom = NULL;
		save.lgioRecom = false;
	}
	if( grid.pnunit != NULL && ( !save.lgSaveGrid_noclobber || lgFinal ) )
	{
		grid.pnunit = NULL;
	}
	ioMAP = NULL;

	return;
}

/*ChkUnits check for keyword UNITS on line, then scan wavelength or energy units if present,
 * units are copied into save.chConPunEnr - when doing output, the routine call
 * AnuUnit( energy ) will automatically return the energy in the right units,
 * when called to do save output */
STATIC void ChkUnits( Parser &p )
{

	DEBUG_ENTRY( "ChkUnits()" );

	/* option to set units for continuum energy in save output */
	if( p.nMatch("UNITS") )
	{
		// p.StandardEnergyUnit() will terminate if no unit was recognized
		save.chConPunEnr[save.nsave] = p.StandardEnergyUnit();
	}
	else
	{
		save.chConPunEnr[save.nsave] = StandardEnergyUnit(" RYD ");
	}
	return;
}
