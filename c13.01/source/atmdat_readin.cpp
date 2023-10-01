/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*atmdat_readin read in some data files, but only if this is very first call, 
 * called by Cloudy */
#include "cddefines.h"
#include "physconst.h"
#include "taulines.h"
#include "mewecoef.h"
#include "iterations.h"
#include "heavy.h"
#include "mole.h"
#include "input.h"
#include "h2.h"
#include "yield.h"
#include "trace.h"
#include "lines.h"
#include "lines_service.h"
#include "ionbal.h"
#include "struc.h"
#include "geometry.h"
#include "dense.h"
#include "dynamics.h"
#include "elementnames.h"
#include "hyperfine.h"
#include "atmdat.h"
#include "iso.h"
/* */
/* this was needed to get array to crash out of bounds if not set.
 * std limits on limits.h did not work with visual studio! */
const long INTBIG = 2000000000;

/* these are the individual pointers to the level 1 lines, they are set to
 * very large negative.  
 * NB NB NB!!
 * these occur two times in the code!!
 * They are declared in taulines.h, and defined here  */
long	ipT1656=INTBIG,    ipT9830=INTBIG,  
	ipT8727=INTBIG,    ipT1335=INTBIG,  ipT1909=INTBIG,
	ipT977=INTBIG,    ipT1550=INTBIG,  ipT1548=INTBIG, ipT386=INTBIG, 
	ipT310=INTBIG,    ipc31175=INTBIG, ipT291=INTBIG,  ipT280=INTBIG,
	ipT274=INTBIG,    ipT270=INTBIG,   ipT312=INTBIG,  ipT610=INTBIG, 
	ipT370=INTBIG,    ipT157=INTBIG,   ipT1085=INTBIG, 
	ipT990=INTBIG,    ipT1486=INTBIG,  ipT765=INTBIG,  ipT1243=INTBIG, 
	ipT1239=INTBIG,   ipT374g=INTBIG,  ipT374x=INTBIG, ipT1200=INTBIG,
	ipT2140=INTBIG,   ipT671=INTBIG,   ipT315=INTBIG,  ipT324=INTBIG, 
	ipT333=INTBIG,    ipT209=INTBIG,   ipT122=INTBIG,  ipT205=INTBIG,
	ipT57=INTBIG,     ipT6300=INTBIG,  ipT6363=INTBIG, ipT5577=INTBIG, 
	ipT834=INTBIG,    ipT1661=INTBIG,  ipT1666=INTBIG, ipT835=INTBIG,
	ipT789=INTBIG,    ipT630=INTBIG,   ipT1304=INTBIG, ipSi10_606=INTBIG,
	ipT1039=INTBIG,   ipT8446=INTBIG,  ipT4368=INTBIG, ipTOI13=INTBIG,
	ipTOI11=INTBIG,   ipTOI29=INTBIG,  ipTOI46=INTBIG, ipTO1025=INTBIG, 
	ipT304=INTBIG,    ipT1214=INTBIG,  ipT150=INTBIG,  ipT146=INTBIG,
	ipT63=INTBIG,     ipTO88=INTBIG,   ipT52=INTBIG,   ipT26=INTBIG, 
	ipT1032=INTBIG,   ipT1037=INTBIG,  ipT770=INTBIG, ipT780=INTBIG,
	ipxNe0676=INTBIG, ipT895=INTBIG,   ipT88=INTBIG, ipTNe13=INTBIG,
	ipTNe36=INTBIG,   ipTNe16=INTBIG,  ipTNe14=INTBIG, ipTNe24=INTBIG, 
	ipT5895=INTBIG,   ipfsNa373=INTBIG,ipfsNa490=INTBIG, ipfsNa421=INTBIG,
	ipxNa6143=INTBIG, ipxNa6862=INTBIG,ipxNa0746=INTBIG, ipMgI2853=INTBIG, 
	ipMgI2026=INTBIG, ipT2796=INTBIG,  ipT2804=INTBIG,
	ipT705=INTBIG,    ipT4561=INTBIG,  ipxMg51325=INTBIG, ipxMg52417=INTBIG, 
	ipxMg52855=INTBIG,ipxMg71190=INTBIG, ipxMg72261=INTBIG,
	ipxMg72569=INTBIG,ipxMg08303=INTBIG, ipTMg610=INTBIG, ipTMg625=INTBIG, 
	ipT58=INTBIG,     ipTMg4=INTBIG, ipTMg14=INTBIG, ipTMg6=INTBIG,
	ipfsMg790=INTBIG, ipfsMg755=INTBIG, ipAlI3957=INTBIG, ipAlI3090=INTBIG, 
	ipT1855=INTBIG,   ipT1863=INTBIG, ipT2670=INTBIG,
	ipAl529=INTBIG,   ipAl6366=INTBIG, ipAl6912=INTBIG, ipAl8575=INTBIG, 
	ipAl8370=INTBIG,  ipAl09204=INTBIG, ipT639=INTBIG,
	ipTAl550=INTBIG,  ipTAl568=INTBIG, ipTAl48=INTBIG, ipSii2518=INTBIG, 
	ipSii2215=INTBIG, ipT1808=INTBIG,
	ipT1207=INTBIG,   ipT1895=INTBIG, ipT1394=INTBIG, ipT1403=INTBIG, 
	ipT1527=INTBIG,   ipT1305=INTBIG, ipT1260=INTBIG, ipSi619=INTBIG,
	ipSi10143=INTBIG, ipTSi499=INTBIG, ipTSi521=INTBIG, ipTSi41=INTBIG, 
	ipTSi35=INTBIG,   ipTSi25=INTBIG, ipTSi65=INTBIG,
	ipTSi3=INTBIG,    ipTSi4=INTBIG, ipP0260=INTBIG, ipP0233=INTBIG, 
	ipP0318=INTBIG,   ipP713=INTBIG, ipT1256=INTBIG, ipT1194=INTBIG,
	ipTS1720=INTBIG,  ipT1198=INTBIG, ipT786=INTBIG,
	ipTS19=INTBIG, ipTS34=INTBIG,    ipTS11=INTBIG, ipfsCl214=INTBIG, ipfsCl233=INTBIG,
	ipCl04203=INTBIG, ipCl04117=INTBIG, ipCl973=INTBIG,
	ipTAr7=INTBIG,    ipTAr9=INTBIG, ipTAr22=INTBIG, ipTAr13=INTBIG,
	ipTAr8=INTBIG,    ipAr06453=INTBIG, ipKI7745=INTBIG, ipxK03462=INTBIG, ipxK04598=INTBIG,
	ipxK04154=INTBIG, ipxK07319=INTBIG, ipCaI4228=INTBIG, ipT3934=INTBIG,
	ipT3969=INTBIG,   ipT8498=INTBIG, ipT8542=INTBIG,
	ipT8662=INTBIG,   ipT7291=INTBIG, ipT7324=INTBIG, ipTCa3=INTBIG,
	ipSc05231=INTBIG, ipSc13264=INTBIG,
	ipFeI3884=INTBIG, ipFeI3729=INTBIG, ipFeI3457=INTBIG,
	ipFeI3021=INTBIG, ipFeI2966=INTBIG, ipTuv3=INTBIG,
	ipTr48=INTBIG,    ipTFe16=INTBIG, ipTFe26=INTBIG, ipTFe34=INTBIG, 
	ipTFe35=INTBIG,   ipTFe46=INTBIG, ipTFe56=INTBIG, ipT1122=INTBIG,
	ipT191=INTBIG,    ipCo11527=INTBIG;
long 	ipS4_1405=INTBIG,ipS4_1398=INTBIG,ipS4_1424=INTBIG,ipS4_1417=INTBIG,ipS4_1407=INTBIG,
	ipO4_1400=INTBIG,ipO4_1397=INTBIG,ipO4_1407=INTBIG,ipO4_1405=INTBIG,ipO4_1401=INTBIG,
	ipN3_1749=INTBIG,ipN3_1747=INTBIG,ipN3_1754=INTBIG,ipN3_1752=INTBIG,ipN3_1751=INTBIG,
	ipC2_2325=INTBIG,ipC2_2324=INTBIG,ipC2_2329=INTBIG,ipC2_2328=INTBIG,ipC2_2327=INTBIG,
	ipSi2_2334=INTBIG,ipSi2_2329=INTBIG,ipSi2_2350=INTBIG,ipSi2_2344=INTBIG,ipSi2_2336=INTBIG,
	ipS1_25m=INTBIG,ipS1_56m=INTBIG,ipCl1_11m=INTBIG,
	ipFe1_24m=INTBIG,ipFe1_35m=INTBIG,ipFe1_54m=INTBIG,ipFe1_111m=INTBIG,
	ipNi1_7m=INTBIG ,ipNi1_11m=INTBIG,ipSi1_130m=INTBIG,ipSi1_68m=INTBIG,
	ipNI_pumpDirect[NI_NDP]={INTBIG, INTBIG, INTBIG, INTBIG, INTBIG, INTBIG, INTBIG, INTBIG, INTBIG},
	ipNI_pumpIndirect=INTBIG;

/* above are the individual pointers to the level 1 lines, they are set to
 * very large negative.  
 * NB NB NB!!
 * these occur two times in the code!!
 * They are declared in taulines.h, and defined here  */

/* definition for whether level 2 lines are enabled, will be set to -1 
 * with no level2 command */
/*long nWindLine = NWINDDIM;*/
/*realnum TauLine2[NWINDDIM][NTA];*/
/*realnum **TauLine2;*/
#include "atomfeii.h"

// NB NB - IS_TOP should always be the last entry of this enum!
typedef enum { IS_NONE, IS_K_SHELL, IS_L1_SHELL, IS_L2_SHELL, IS_TOP } exc_type;

struct t_BadnellLevel
{
	string config;
	int irsl;
	int S;
	int L;
	int g; // 2*J+1
	realnum energy; // in cm^-1
	bool lgAutoIonizing;
	exc_type WhichShell;
	t_BadnellLevel() : irsl(0), S(0), L(0), g(0), energy(0.f), lgAutoIonizing(false), WhichShell(IS_NONE) {}
};

// read Hummer and Storey 98 He1 photoionization cross-section data
STATIC void read_SH98_He1_cross_sections(void);

STATIC void read_Helike_cross_sections(void);

// read autoionization data from Badnell data file
STATIC void ReadBadnellAIData(const string& fnam,      // filename containing the Badnell data
			      long nelem,              // nelem is on C scale
			      long ion,                // ion is on C scale
			      TransitionList& UTA,     // UTA lines will be pushed on this stack
			      bitset<IS_TOP> Skip);    // option to skip transitions from a particular shell
// simple helper functions for ReadBadnellAIData
inline void InitTransition(const TransitionProxy& t);
inline int irsl2ind(vector<t_BadnellLevel>& level, int irsl);

// UTA lines below this absorption oscillator strength value will be ignored -
// F+13 paper plotted f not gf
const realnum f_cutoff = 1.e-4f;

void atmdat_readin(void)
{
	long int i, 
	  iCase ,
	  ion, 
	  ipDens ,
	  ipISO ,
	  ipTemp ,
	  j, 
	  ig0, 
	  ig1, 
	  imax, 
	  nelem, 
	  nelec, 
	  magic1, 
	  magic2,
	  mol;

	char cha;
	char chS2[3];

	long ipZ;
	bool lgEOL;

	FILE *ioDATA;
	char chLine[FILENAME_PATH_LENGTH_2] , 
		chFilename[FILENAME_PATH_LENGTH_2];

	static bool lgFirstCall = true;

	DEBUG_ENTRY( "atmdat_readin()" );

	/* do nothing if not first call */
	if( !lgFirstCall )
	{ 
		/* do not do anything, but make sure that number of zones has not increased */
		bool lgTooBig = false;
		for( j=0; j < iterations.iter_malloc; j++ )
		{
			if( geometry.nend[j]>=struc.nzlim )
				lgTooBig = true;
		}
		if( lgTooBig )
		{
			fprintf(ioQQQ," This is the second or later calculation in a grid.\n");
			fprintf(ioQQQ," The number of zones has been increased beyond what it was on the first calculation.\n");
			fprintf(ioQQQ," This can\'t be done since space has already been allocated.\n");
			fprintf(ioQQQ," Have the first calculation do the largest number of zones so that an increase is not needed.\n");
			fprintf(ioQQQ," Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
		return;
	}

	lgFirstCall = false; /* do not reevaluate again */

	/* make sure that molecules have been initialized - this will fail
	 * if this routine is called before size of molecular network is known */
	if( !mole_global.num_total )
	{
		/* mole_global.num_comole_calc can't be zero */
		TotalInsanity();
	}

	/* create space for the structure variables */
	/* nzlim will be limit, and number allocated */
	/* >>chng 01 jul 28, define this var, do all following mallocs */
	for( j=0; j < iterations.iter_malloc; j++ )
	{
		struc.nzlim = MAX2( struc.nzlim , geometry.nend[j] );
	}

	/* sloppy, but add one extra for safely */
	++struc.nzlim;

	struc.coolstr = ((double*)MALLOC( (size_t)(struc.nzlim)*sizeof(double )));

	struc.heatstr = ((double*)MALLOC( (size_t)(struc.nzlim)*sizeof(double )));

	struc.testr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.volstr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.drad_x_fillfac = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.histr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.hiistr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.ednstr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.o3str = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.pressure = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.windv = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.AccelTotalOutward = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.AccelGravity = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.pres_radiation_lines_curr = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.GasPressure = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.hden = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.DenParticles = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.DenMass = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.drad = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.depth = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.depth_last = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.drad_last = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	struc.xLyman_depth = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	if (mole_global.num_calc != 0)
	{
		struc.molecules = ((realnum**)MALLOC( (size_t)(mole_global.num_calc)*sizeof(realnum* )));
	}
	else
	{
		struc.molecules = NULL;
	}

	for(mol=0;mol<mole_global.num_calc;mol++)
	{
		struc.molecules[mol] = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));
	}

	struc.H2_abund = ((realnum*)MALLOC( (size_t)(struc.nzlim)*sizeof(realnum )));

	/* create space for gas phase abundances array, first create space for the elements */
	struc.gas_phase = (realnum **)MALLOC(sizeof(realnum *)*(unsigned)LIMELM );
	/* last the zones */
	for( ipZ=0; ipZ< LIMELM;++ipZ )
	{
		struc.gas_phase[ipZ] = 
			(realnum*)MALLOC(sizeof(realnum)*(unsigned)(struc.nzlim) );
	}

	/* create space for struc.xIonDense array, first create space for the zones */
	struc.xIonDense = (realnum ***)MALLOC(sizeof(realnum **)*(unsigned)(LIMELM) );

	for( ipZ=0; ipZ<LIMELM;++ipZ )
	{
		/* space for the ionization stage */
		struc.xIonDense[ipZ] = 
			(realnum**)MALLOC(sizeof(realnum*)*(unsigned)(LIMELM+1) );
		/* now create diagonal of space for zones */
		for( ion=0; ion < (LIMELM+1); ++ion )
		{
			struc.xIonDense[ipZ][ion] = 
				(realnum*)MALLOC(sizeof(realnum)*(unsigned)(struc.nzlim) );
		}
	}

	struc.StatesElem = (realnum ****)MALLOC(sizeof(realnum ***)*(unsigned)(LIMELM) );
	for( nelem=ipHYDROGEN; nelem<LIMELM; ++nelem)
	{
		if( dense.lgElmtOn[nelem] )
		{
			struc.StatesElem[nelem] = (realnum***)MALLOC(sizeof(realnum**)*(unsigned)(nelem+1) );
			for( long ion=0; ion<nelem+1; ion++ )
			{
				long ipISO = nelem-ion;
				if( ipISO < NISO )
				{
					struc.StatesElem[nelem][ion] = (realnum**)MALLOC(sizeof(realnum*)*(unsigned)iso_sp[ipISO][nelem].numLevels_max);
					for( long level=0; level < iso_sp[ipISO][nelem].numLevels_max; ++level )
					{
						struc.StatesElem[nelem][ion][level] = 
							(realnum*)MALLOC(sizeof(realnum)*(unsigned)(struc.nzlim) );
					}
				}
				else
				{
					fixit();  // for now, point non-iso ions to NULL
					struc.StatesElem[nelem][ion] = NULL;
				}
			}
		}
	}

	/* some structure variables */
	for( i=0; i < struc.nzlim; i++ )
	{
		struc.testr[i] = 0.;
		struc.volstr[i] = 0.;
		struc.drad_x_fillfac[i] = 0.;
		struc.histr[i] = 0.;
		struc.hiistr[i] = 0.;
		struc.ednstr[i] = 0.;
		struc.o3str[i] = 0.;
		struc.heatstr[i] = 0.;
		struc.coolstr[i] = 0.;
		struc.pressure[i] = 0.;
		struc.pres_radiation_lines_curr[i] = 0.;
		struc.GasPressure[i] = 0.;
		struc.DenParticles[i] = 0.;
		struc.depth[i] = 0.;
		for(mol=0;mol<mole_global.num_calc;mol++)
		{
			struc.molecules[mol][i] = 0.;
		}
		struc.H2_abund[i] = 0.;
	}

	/* allocate space for some arrays used by dynamics routines, and zero out vars */
	DynaCreateArrays( );

	/*************************************************************
	 *                                                           *
	 * get the level 1 lines, with real collision data set       *
	 *                                                           *
	 *************************************************************/

	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_readin reading level1.dat\n");

	ioDATA = open_data( "level1.dat", "r" );

	/* first line is a version number and does not count */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_readin could not read first line of level1.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* count how many lines are in the file, ignoring all lines
	 * starting with '#' */
	nLevel1 = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* we want to count the lines that do not start with #
		 * since these contain data */
		if( chLine[0] != '#')
			++nLevel1;
	}

	/* now malloc the TauLines array */
	TauLines.resize(nLevel1+1);
	AllTransitions.push_back(TauLines);

	for( i=0; i<nLevel1+1; i++ )
	{
		TauLines[i].Junk();
		TauLines[i].AddLoState();
		TauLines[i].AddHiState();
		TauLines[i].AddLine2Stack();
	}

	/* now rewind the file so we can read it a second time*/
	if( fseek( ioDATA , 0 , SEEK_SET ) != 0 )
	{
		fprintf( ioQQQ, " atmdat_readin could not rewind level1.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " atmdat_readin could not read first line of level1.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* level 1 magic number */
	nelem = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	nelec = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	ion = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

	/* magic number
	 * the following is the set of numbers that appear at the start of level1.dat  */
	int nYr=12, nMonth=1, nDay=13;
	if( ( nelem != nYr ) || ( nelec != nMonth ) || ( ion != nDay ) )
	{
		fprintf( ioQQQ, 
			" atmdat_readin: the version of level1.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			" Please obtain the current version from the Cloudy web site.\n" );
		fprintf( ioQQQ, 
			" I expected to find the number %i %i %i and got %2.2li %2.2li %2.2li instead.\n" ,
			nYr, nMonth, nDay,
			nelem , nelec , ion );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* this starts at 1 not 0 since zero is reserved for the
	 * dummy line - we counted number of lines above */
	i = 1;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		if( i >= (nLevel1+1) )
			TotalInsanity();

		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			/* first two cols of line has chem element symbol, 
			 * try to match it with the list of names 
			 * there are no H lines in the level 1 set so zero
			 * is an invalid result */

			/* copy first two char into chS2 and null terminate */
			strncpy( chS2 , chLine , 2);
			chS2[2] = 0;

			ipZ = 0;
			for( j=0; j<LIMELM; ++j)
			{
				if( strcmp( elementnames.chElementSym[j], chS2 ) == 0 )
				{
					/* ipZ is the actual atomic number starting from 1 */
					ipZ = j + 1;
					break;
				}
			}

			/* this happens if no valid chemical element in first two cols on this line */
			if( ipZ == 0 )
			{
				fprintf( ioQQQ, 
					" atmdat_readin could not identify chem symbol on this level 1line:\n");
				fprintf( ioQQQ, "%s\n", chLine );
				fprintf( ioQQQ, "looking for this string==%2s==\n",chS2 );
				cdEXIT(EXIT_FAILURE);
			}

			/* now stuff them into the array, will convert over to proper units later */
			(*TauLines[i].Hi()).nelem() = (int)ipZ;
			/* start the scan early on the line -- the element label will not be
			 * picked up, first number will be the ion stage after the label, as in c  4 */
			j = 1;
			(*TauLines[i].Hi()).IonStg() = (int)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 1.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			(*TauLines[i].Lo()).nelem() = (*TauLines[i].Hi()).nelem();
			(*TauLines[i].Lo()).IonStg() = (*TauLines[i].Hi()).IonStg();

			TauLines[i].WLAng() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 2.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			TauLines[i].EnergyWN() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 3.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			(*TauLines[i].Lo()).g() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 4.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			(*TauLines[i].Hi()).g() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 5.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			/* gf is log if negative */
			TauLines[i].Emis().gf() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 6.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			if( TauLines[i].Emis().gf() < 0. )
				TauLines[i].Emis().gf() = (realnum)pow((realnum)10.f,TauLines[i].Emis().gf());

			/* Emis().Aul() is log if negative */
			TauLines[i].Emis().Aul() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 7.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}
			if( TauLines[i].Emis().Aul() < 0. )
				TauLines[i].Emis().Aul() = (realnum)pow((realnum)10.f,TauLines[i].Emis().Aul());

			TauLines[i].Emis().iRedisFun() = (int)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			if( lgEOL )
			{
				fprintf( ioQQQ, " There should have been a number on this level1 line 8.   Sorry.\n" );
				fprintf( ioQQQ, "string==%s==\n" ,chLine );
				cdEXIT(EXIT_FAILURE);
			}

			/* this is new in c94.01 and returns nothing (0.) for most lines */
			TauLines[i].Coll().col_str() = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);

			/* finally increment i */
			++i;
		}
	}

	/* now close the file */
	fclose(ioDATA);
	if( trace.lgTrace )
		fprintf( ioQQQ, " reading level1.dat OK\n" );


	/*************************************************************
	 *                                                           *
	 * get the level 2 line, opacity project, data set           *
	 *                                                           *
	 *************************************************************/

	/* nWindLine is initialized to the dimension of the vector when it is
	 * initialized in the definition at the start of this file.  
	 * it is set to -1 with the "no level2" command, which
	 * stops us from trying to establish this vector */
	if( nWindLine > 0 )
	{
		/* begin level2 level 2 wind block */
		/* open file with level 2 line data */

		/* create the TauLine2 emline array */
		TauLine2.resize(nWindLine);
		AllTransitions.push_back(TauLine2);
		cs1_flag_lev2 = ((realnum *)MALLOC( (size_t)nWindLine*sizeof(realnum )));

		/* first initialize entire array to dangerously large negative numbers */
		for( i=0; i< nWindLine; ++i )
		{
			/* >>chng 99 jul 16, from setting all t[] to flt_max, to call for
			 * following, each member of structure set to own type of impossible value */
			TauLine2[i].Junk();

			TauLine2[i].AddHiState();
			TauLine2[i].AddLoState();
			TauLine2[i].AddLine2Stack();
		}

		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading level2.dat\n");

		ioDATA = open_data( "level2.dat", "r" );

		/* get magic number off first line */
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " level2.dat error getting magic number\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* scan off three numbers, should be the yr mn dy of creation date */
		i = 1;
		nelem = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
		nelec = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
		ion = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);

		/* level2.dat magic number
		 * the following is the set of numbers that appear at the start of level2.dat */
		if( lgEOL || ( nelem != 9 ) || ( nelec != 11 ) || ( ion != 18 ) )
		{
			fprintf( ioQQQ, 
				" atmdat_readin: the version of level2.dat is not the current version.\n" );
			fprintf( ioQQQ, 
				" I expected to find the number 09 11 18 and got %2.2li %2.2li %2.2li instead.\n" ,
				nelem , nelec , ion );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			fprintf( ioQQQ, "Please obtain the correct version.\n");
			cdEXIT(EXIT_FAILURE);
		}

		/* now get the actual data */
		i = 0;
		while( i < nWindLine )
		{
			/* this must be double for sscanf to work below */
			double tt[7];
			if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
			{
				fprintf( ioQQQ, " level2.dat error getting line  %li\n", i );
				cdEXIT(EXIT_FAILURE);
			}
			/* option to skip any line starting with # */
			if( chLine[0]!='#' )
			{
				/*printf("string is %s\n",chLine );*/
				sscanf( chLine , "%lf %lf %lf %lf %lf %lf %lf " , 
					&tt[0] ,
					&tt[1] ,
					&tt[2] ,
					&tt[3] ,
					&tt[4] ,
					&tt[5] ,
					&tt[6] );
				/* these are readjusted into their final form in the structure 
				 * in routine lines_setup*/
				(*TauLine2[i].Hi()).nelem() = (int)tt[0];
				(*TauLine2[i].Hi()).IonStg() = (int)tt[1];
				(*TauLine2[i].Lo()).g() = (realnum)tt[2];
				(*TauLine2[i].Hi()).g() = (realnum)tt[3];
				TauLine2[i].Emis().gf() = (realnum)tt[4];
				TauLine2[i].EnergyWN() = (realnum)tt[5];
				cs1_flag_lev2[i] = (realnum)tt[6];
				++i;
			}
		}
		/* get magic number off last line */
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " level2.dat error getting last magic number\n" );
			cdEXIT(EXIT_FAILURE);
		}
		sscanf( chLine , "%ld" , &magic2 );
		if( 999 != magic2 )
		{
			fprintf( ioQQQ, " level2.dat ends will wrong magic number=%ld \n", 
			  magic2 );
			cdEXIT(EXIT_FAILURE);
		}
		fclose( ioDATA );
		if( trace.lgTrace )
			fprintf( ioQQQ," reading level2.dat OK\n");

		/* end of level2 block*/
	}

	/* the UTA line sets - trace will print summary of various data sources */

	/* reserve space for all data sets, no worries if too small though... */
	UTALines.reserve( 4000 );
	AllTransitions.push_back(UTALines);

	// version of element symbols in lower case and without spaces....
	const char* chElmSymLC[] =
		{ "h", "he", "li", "be", "b", "c", "n", "o", "f", "ne", "na", "mg", "al", "si", "p",
		  "s", "cl", "ar", "k", "ca", "sc", "ti", "v", "cr", "mn", "fe", "co", "ni", "cu", "zn" };

	// save cite for UTA sources, insure no double counting
	char chUTA_ref[LIMELM][LIMELM][5];
	for( long nelem=0; nelem < LIMELM; ++nelem )
		for( long ion=0; ion <= nelem; ++ion )
			strcpy( chUTA_ref[nelem][ion] , "" );

	/* first read in the Badnell data */
	for( ipISO=ipLI_LIKE; ipISO < ipAL_LIKE; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			// ion = 0 for neutral atom
			long ion = nelem - ipISO;
			strcat( chUTA_ref[nelem][ion] , "B" );

			bitset<IS_TOP> Skip;
			if( ipISO < ipNA_LIKE )
			{
				// construct file name
				ostringstream oss;
				// Badnell calls Li-like series He-like, etc...
				oss << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				// Badnell uses ion = 1 for neutral atom
				oss << chElmSymLC[nelem] << ion+1 << "ic1-2.dat";
				// now read the data...
				ReadBadnellAIData( oss.str(), nelem, ion, UTALines, Skip );
			}
			else
			{
				// from Na-like onwards both K-shell (ic1-3) and L-shell (ic2-3)
				// excitation are treated in two separate files
				ostringstream oss;
				oss << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				oss << chElmSymLC[nelem] << ion+1 << "ic1-3.dat";
				ReadBadnellAIData( oss.str(), nelem, ion, UTALines, Skip );

				// Kisielius L1, L2-shell data sets take precedence for Na, Mg, and Al-like iron
				if( ionbal.lgInnerShell_Kisielius && nelem == ipIRON &&
				    ipISO >= ipNA_LIKE && ipISO <= ipAL_LIKE )
				{
					Skip[IS_L1_SHELL] = 1;
					Skip[IS_L2_SHELL] = 1;
				}

				ostringstream oss2;
				oss2 << "UTA/nrb00_" << chElmSymLC[ipISO-1] << "_";
				oss2 << chElmSymLC[nelem] << ion+1 << "ic2-3.dat";
				ReadBadnellAIData( oss2.str(), nelem, ion, UTALines, Skip );
			}
		}
	}

	/* these are the statistical weights for the ground levels of all ions of iron */
	const realnum StatWeightGroundLevelIron[] =
		{ 9.f, 10.f, 9.f, 6.f, 1.f, 4.f, 5.f, 4.f, 1.f, 4.f, 5.f, 4.f, 1.f,
		  2.f, 1.f, 2.f, 1.f, 4.f, 5.f, 4.f, 1.f, 2.f, 1.f, 2.f, 1.f, 2.f };

	// blank line that will be pushed on the UTA line stack
	qList BlankStates(1);
	TransitionList BlankList("BlankList",&BlankStates,1);
	TransitionList::iterator BlankLine = BlankList.begin();
	(*BlankLine).Junk();

	/* next read in the Gu file */
	if( ionbal.lgInnerShell_Gu06 )
	{
		/* read the Gu et al. (2006) data
		 * >>refer	Fe	UTA	Gu, M. F., Holczer T., Behar E., & Kahn S. M. 2006, ApJ 641, 1227-1232 */
		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading UTA_Gu06.dat\n");

		FILE *ioGU06 = open_data( "UTA/UTA_Gu06.dat", "r" );

		nelem = 0;
		/* get up to magic number in Gu 06 file - there is a large header of 
		 * comments with the first data the magic number */
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioGU06 ) != NULL )
		{
			/* we want to break on first line that does not start with #
			 * since that contains data */
			if( chLine[0] != '#')
				break;
		}
		/* now get Gu magic number */
		sscanf( chLine, "%li %li %li", &nelem, &nelec, &ion );

		/* is magic number correct?
		 * the following is the set of numbers that appear at the start of UTA_Gu06.dat 2006 11 17 */
		if( nelem != 2007 || nelec != 1 || ion != 23 )
		{
			fprintf( ioQQQ, 
				" atmdat_readin: the version of UTA_Gu06.dat is not the current version.\n" );
			fprintf( ioQQQ, 
				" I expected to find the number 2007 1 23 and got %li %li %li instead.\n" ,
				nelem , nelec , ion );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}

		int nelemGu =-1, ionGu=-1;
		while( read_whole_line( chLine, (int)sizeof(chLine), ioGU06 ) != NULL )
		{
			if( chLine[0] != '#' )
			{
				long int i2;
				double WLAng, Aul, oscill, Aauto;

				sscanf( chLine, "%4li%5li%8lf%13lf%12lf",
					&ion, &i2, &WLAng, &Aul, &Aauto );
				sscanf( &chLine[54], "%13lf", &oscill );

				// avoid duplication of ions: anything upto and including the Mg-like
				// series is covered by the Badnell data set, Al-like is covered by
				// the Kisielius data set if it is turned on.
				ipISO = ipIRON - ion + 1;
				int ipThres = ionbal.lgInnerShell_Kisielius ? ipAL_LIKE : ipMG_LIKE;
				if( ipISO <= ipThres )
					continue;

				UTALines.push_back( *BlankLine );
				InitTransition( UTALines.back() );

				/* all these are iron, first number was ion stage with 0 the atom */
				(*UTALines.back().Hi()).nelem() = ipIRON+1;

				/* now do stage of ionization */
				ASSERT( ion > 0 && ion <= ipIRON );
				/* the ion stage - 1 is atom - this data file has different
				 * format from other two - others are 0 for atom */
				(*UTALines.back().Hi()).IonStg() = ion;
				if( ipIRON!=nelemGu || ion!=ionGu )
				{
					// one label per ion
					nelemGu = ipIRON;
					ionGu = ion;
					strcat( chUTA_ref[ipIRON][ion-1] , "G" );
				}

				/* these are the statistical weights 
				 * lower levels are not included in the original data file */
				if( strstr_s( chLine, "(J=1/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 2.f;
				else if( strstr_s( chLine, "(J=1)" ) != NULL )
					(*UTALines.back().Hi()).g() = 3.f;
				else if( strstr_s( chLine, "(J=3/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 4.f;
				else if( strstr_s( chLine, "(J=2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 5.f;
				else if( strstr_s( chLine, "(J=5/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 6.f;
				else if( strstr_s( chLine, "(J=3)" ) != NULL )
					(*UTALines.back().Hi()).g() = 7.f;
				else if( strstr_s( chLine, "(J=7/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 8.f;
				else if( strstr_s( chLine, "(J=4)" ) != NULL )
					(*UTALines.back().Hi()).g() = 9.f;
				else if( strstr_s( chLine, "(J=9/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 10.f;
				else if( strstr_s( chLine, "(J=5)" ) != NULL )
					(*UTALines.back().Hi()).g() = 11.f;
				else if( strstr_s( chLine, "(J=11/2)" ) != NULL )
					(*UTALines.back().Hi()).g() = 12.f;
				else
					TotalInsanity();
				(*UTALines.back().Lo()).g() = StatWeightGroundLevelIron[ion-1];

				/* wavelength in Angstroms */
				UTALines.back().WLAng() = (realnum)WLAng;
				UTALines.back().EnergyWN() = (realnum)(1e8/WLAng);

				/* store branching ratio for autoionization */
				double frac_ioniz = Aauto/(Aul + Aauto);
				ASSERT( frac_ioniz >= 0. &&  frac_ioniz <= 1. );
				UTALines.back().Emis().AutoIonizFrac() = (realnum)frac_ioniz;

				/* save gf scanned from line */
				UTALines.back().Emis().gf() = (*UTALines.back().Lo()).g() * (realnum)oscill;

				/* this is true spontaneous rate for doubly excited state to inner 
				 * shell UTA, and is used for pumping, and also relaxing back to inner shell */
				UTALines.back().Emis().Aul() =
					(realnum)eina( UTALines.back().Emis().gf(),
					UTALines.back().EnergyWN(), (*UTALines.back().Hi()).g() );

				ASSERT( fp_equal_tol( (realnum)Aul, UTALines.back().Emis().Aul(), 1.e-3f*(realnum)Aul ) );

				UTALines.back().Emis().iRedisFun() = ipPRD;

				UTALines.back().Emis().dampXvel() = (realnum)(
						(UTALines.back().Emis().Aul()+Aauto) /
						UTALines.back().EnergyWN()/PI4);
				ASSERT( UTALines.back().Emis().dampXvel()>0. );

				// ignore line if too weak
				if( UTALines.back().Emis().gf() < StatWeightGroundLevelIron[ion-1] * f_cutoff )
					UTALines.pop_back();
			}
		}

		fclose( ioGU06 );

		if( trace.lgTrace )
			fprintf( ioQQQ, " reading UTA_Gu06.dat OK\n" );
	}
	else
	{
		/* this branch, get the Behar 01 data */
		/* >>refer	Fe	UTA	Behar, E., Sako, M., & Kahn, S. M. 2001, ApJ, 563, 497-504 */
		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading UTA_Behar.dat\n");

		FILE *ioBEHAR = open_data( "UTA/UTA_Behar.dat", "r" );

		/* UTA_Behar.dat magic number */
		nelem = 0;
		if( read_whole_line( chLine, (int)sizeof(chLine), ioBEHAR ) != NULL )
			sscanf( chLine, "%li %li %li", &nelem, &nelec, &ion );

		/* magic number
		 * the following is the set of numbers that appear at the start of UTA_Behar.dat 2002 8 19 */
		if( nelem != 2002 || nelec != 10 || ion != 22 )
		{
			fprintf( ioQQQ, 
				" atmdat_readin: the version of UTA_Behar.dat is not the current version.\n" );
			fprintf( ioQQQ, 
				" I expected to find the number 2002 10 22 and got %li %li %li instead.\n" ,
				nelem , nelec , ion );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}

		/* read in the Behar 01 data */
		int nelemBehar=-1, ionBehar=-1;
		while( read_whole_line( chLine, (int)sizeof(chLine), ioBEHAR ) != NULL )
		{
			if( chLine[0] != '#' )
			{
				long int i1, i2, i3;
				double f1, f2, oscill;
				double frac_relax;

				sscanf( chLine ,"%li\t%li\t%li\t%lf\t%lf\t%lf\t%lf",
					&i1,&i2,&i3,&f1,&f2,&frac_relax,&oscill );

				// skip line if too weak
				if( oscill < f_cutoff )
					continue;

				// avoid duplication of ions: anything upto and including the Mg-like
				// series is covered by the Badnell data set, Al-like is covered by
				// the Kisielius data set if it is turned on.
				ipISO = ipIRON - i1;
				int ipThres = ionbal.lgInnerShell_Kisielius ? ipAL_LIKE : ipMG_LIKE;
				if( ipISO <= ipThres )
					continue;

				UTALines.push_back( *BlankLine );
				InitTransition( UTALines.back());

				/* all these are iron, first number was ion stage with 0 the atom */
				(*UTALines.back().Hi()).nelem() = ipIRON+1;

				/* now do stage of ionization */
				ASSERT( i1 >= 0 && i1 < ipIRON );
				/* NB - the plus one is because 1 will be subtracted when used,
				* in the original data file i1 is 0 for the atom */
				(*UTALines.back().Hi()).IonStg() = i1 + 1;
				if( ipIRON!=nelemBehar || i1!=ionBehar )
				{
					// one label per ion
					nelemBehar = ipIRON;
					ionBehar = i1;
					strcat( chUTA_ref[ipIRON][i1] , "b" );
				}

				/* wavelength in Angstroms */
				UTALines.back().WLAng() = (realnum)f1;
				UTALines.back().EnergyWN() = 1.e8f/(realnum)f1;

				(*UTALines.back().Lo()).g() = StatWeightGroundLevelIron[i1];
				/* now derive g_up by comparing Aul (is f2) and the oscillator strength
				 * this procedure will cause flu_test to be a factor g_up too small */
				double flu_test = GetGF( f2, UTALines.back().EnergyWN(), 1. )/(*UTALines.back().Lo()).g();
				/* so dividing the real oscillator strength by flu_test yields g_up... */
				(*UTALines.back().Hi()).g() = (realnum)nint( oscill/flu_test );

				ASSERT( i2 == 0 || fp_equal( (*UTALines.back().Hi()).g(), (realnum)i2 ) );
				ASSERT( i3 == 0 || fp_equal( (*UTALines.back().Lo()).g(), (realnum)i3 ) );

				/* this is true spontaneous rate for doubly excited state to inner shell UTA,
				* and is used for pumping, and also relaxing back to inner shell */
				UTALines.back().Emis().Aul() = (realnum)f2;

				UTALines.back().Emis().gf() = (*UTALines.back().Lo()).g() * (realnum)oscill;

				UTALines.back().Emis().iRedisFun() = ipPRD;

				UTALines.back().Emis().dampXvel() = (realnum)(
						(UTALines.back().Emis().Aul()/frac_relax) /
						UTALines.back().EnergyWN()/PI4);
				ASSERT( UTALines.back().Emis().dampXvel()>0. );

				/* store branching ratio for autoionization */
				ASSERT( frac_relax >= 0.f &&  frac_relax <= 1.f );
				UTALines.back().Emis().AutoIonizFrac() = 1.f-(realnum)frac_relax;
			}
		}

		fclose( ioBEHAR );

		if( trace.lgTrace )
			fprintf( ioQQQ, " reading UTA_Behar.dat OK\n" );
	}

	if( ionbal.lgInnerShell_Kisielius )
	{
		/* last read in the Romas Kisielius data
		 *>>refer	Fe	UTA	Kisielius, R., Hibbert, A.. Ferland, G. J., et al. 2003, MNRAS, 344, 696 */
		if( trace.lgTrace )
			fprintf( ioQQQ," atmdat_readin reading UTA_Kisielius.dat\n");

		FILE *ioROMAS = open_data( "UTA/UTA_Kisielius.dat", "r" );

		// magic number
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioROMAS ) != NULL )
		{
			// skip any #
			if( chLine[0] != '#')
				break;
		}
		sscanf( chLine, "%li %li %li", &nelem, &nelec, &ion );
		if( nelem != 11 || nelec != 8 || ion != 25 )
		{
			fprintf( ioQQQ,
				" atmdat_readin: the version of UTA_Kisielius.dat is not the current version.\n" );
			fprintf( ioQQQ,
				" I expected to find the number 11 8 25 and got %li %li %li instead.\n" ,
				nelem , nelec , ion );
			fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
			cdEXIT(EXIT_FAILURE);
		}

		long int nRomasUsed = 0 , nRomasTotal = 0;
		int nelemRomas=-1 , ionRomas=-1;
		FILE *ioROMASused;
		bool lgSaveRomasUsed = false;
		if( lgSaveRomasUsed )
		{
			if( (ioROMASused=fopen("RomasUsed.txt","w"))==NULL)
			{
				fprintf(ioQQQ,"could not open RomasUsed.txt\n");
				cdEXIT(EXIT_FAILURE);
			}
		}
		while( read_whole_line( chLine, (int)sizeof(chLine), ioROMAS ) != NULL )
		{
			/* only look at lines without '#' in first col */
			if( chLine[0] != '#' )
			{
				long int i1, i2, i3;
				double f1, f2, oscill;
				double frac_relax;

				++nRomasTotal;
				sscanf( chLine, "%li\t%li\t%li\t%lf\t%lf\t%lf\t%lf",
					&i1,&i2,&i3,&f1,&f2,&frac_relax,&oscill );

				// skip line if too weak
				// Fe+13 has 25 000 lines, most are vastly weak and do not add to integrated rate
				// the cutoff of 1e-4 reduces the number of Romas UTA lines from 84224 to 573
				if( oscill < f_cutoff )
					continue;

				if( lgSaveRomasUsed )
					fprintf(ioROMASused , "%s" , chLine);

				/* For Fe XIV both levels of the 2P* ground term are present in the data.
				 * following true, ignore the transitions from the excited level
				 * false, assume ground term populated by statistical weight if
				 * "fudge 0" also appears in input stream */
				const bool lgAllowSplitFe14 = false;
				if( lgAllowSplitFe14 || i2 == StatWeightGroundLevelIron[i1] )
				{
					++nRomasUsed;
					UTALines.push_back( *BlankLine );
					InitTransition( UTALines.back());

					/* all these are iron, first number was ion stage with 0 the atom */
					(*UTALines.back().Hi()).nelem() = ipIRON+1;

					/* now do stage of ionization */
					ASSERT( i1 >= 0 && i1 < ipIRON );
					/* NB - the plus one is because 1 will be subtracted when used,
					 * in the original data file i1 is 0 for the atom */
					(*UTALines.back().Hi()).IonStg() = i1 + 1;

					realnum facpop = 1.;
					// allow low or high density limit for level populations in ground term of Fe XIV
					if( i1 == 13 )
					{
						// Fe XIV - fudge -1 returns number of fudge parameters, 0 if not specified
						if( lgAllowSplitFe14 && fudge(-1) )
						{
							if( i2 == StatWeightGroundLevelIron[i1] )
								facpop= 0.333;
							else
								facpop = 0.667;
						}
						else
						{
							if( i2 == StatWeightGroundLevelIron[i1] )
								facpop= 1.;
							else
								facpop = 1e-10;
						}
					}
					if( ipIRON!=nelemRomas || i1!=ionRomas )
					{
						// one label per ion
						nelemRomas = ipIRON;
						ionRomas = i1;
						strcat( chUTA_ref[ipIRON][i1] , "K" );
					}

					/* these were the statistical weights */
					(*UTALines.back().Hi()).g() = (realnum)i3;
					(*UTALines.back().Lo()).g() = (realnum)i2;

					UTALines.back().WLAng() = (realnum)f1;
					UTALines.back().EnergyWN() = 1.e8f/(realnum)f1;
					// print transitions contributing to 15.5A UTA feature
#					if 0
					if( i1==13 && f1>15.35 && f1<15.55)
					{
						fprintf(ioQQQ,"DEBUG %li\t%.5f\t%.3e\n",i2, f1 , oscill * facpop);
					}
#					endif

					/* this is true spontaneous rate for doubly excited state to inner shell,
					 * and is used for pumping, and also relaxing back to inner shell */
					UTALines.back().Emis().gf() = (*UTALines.back().Lo()).g() * (realnum)oscill*facpop;
					UTALines.back().Emis().Aul() =
						(realnum)eina( UTALines.back().Emis().gf(),
							UTALines.back().EnergyWN(),
							(*UTALines.back().Hi()).g() );

					UTALines.back().Emis().iRedisFun() = ipPRD;

					/* store branching ratio for autoionization */
					ASSERT( frac_relax >= 0.f &&  frac_relax <= 1.f );
					UTALines.back().Emis().AutoIonizFrac() = 1.f-(realnum)frac_relax;

					// Romas data do not have autoionization rates so take typical
					// value of 1e15 s-1, suggested by Badnell OP data
					UTALines.back().Emis().dampXvel() = (realnum)(
							(UTALines.back().Emis().Aul()+1e15) /
							UTALines.back().EnergyWN()/PI4);
					ASSERT( UTALines.back().Emis().dampXvel()>0. );
					//fprintf(ioQQQ,"DEBUGGG %li %.3e\n", nRomasUsed, UTALines.back().Emis().dampXvel() );
				}
			}
		}

		fclose( ioROMAS );
		if( lgSaveRomasUsed )
			fclose( ioROMASused );

		if( trace.lgTrace )
			fprintf( ioQQQ, " reading UTA_Kisielius.dat OK,used %li lines from a total of %li\n" , nRomasUsed , nRomasTotal );
	}

	nUTA = (long)UTALines.size();

	/* option to dump UTA lines */
	if( trace.lgTrace )
	{
		fprintf(ioQQQ,"\nUTA data sources; B=Badnell 05; G==Gu 06, b=Behar, K=2011 paper\n");
		fprintf(ioQQQ," ion ");
		for( long ion=0; ion<=LIMELM; ++ion )
			fprintf(ioQQQ,"%4li",ion);
		fprintf(ioQQQ,"\n");
		for( long nelem=0; nelem<LIMELM; ++nelem )
		{
			fprintf(ioQQQ,"%4s ", elementnames.chElementSym[nelem] );
			for( long ion=0; ion<=nelem; ++ion )
			{
				fprintf(ioQQQ,"%4s",chUTA_ref[nelem][ion] );
			}
			fprintf(ioQQQ,"\n");
		}
		fprintf(ioQQQ," ion ");
		for( long ion=0; ion<=LIMELM; ++ion )
			fprintf(ioQQQ,"%4li",ion);
		fprintf(ioQQQ,"\n\n");
	}

	if( false )
	{
		for( i=0; i < nUTA; ++i )
			dprintf( ioQQQ, "%5ld %s %2ld wavl %7.3f glo %2g gup %2g Aul %.2e gf %.2e ai branch %.3f\n",
				 i,
				 elementnames.chElementSym[(*UTALines[i].Hi()).nelem()-1],
				 (*UTALines[i].Hi()).IonStg(),
				 UTALines[i].WLAng(),
				 (*UTALines[i].Lo()).g(),
				 (*UTALines[i].Hi()).g(),
				 UTALines[i].Emis().Aul(),
				 UTALines[i].Emis().gf(),
				 UTALines[i].Emis().AutoIonizFrac() );
		cdEXIT(EXIT_SUCCESS);
	}

	/* end UTA */

	/* read in data for the set of hyperfine structure lines, and allocate
	 * space for the transition HFLines[nHFLines] structure */
	HyperfineCreate();

	/* Make sure that if hybrid is on, then Stout/Chianti are on */
	if( atmdat.lgChiantiHybrid && !atmdat.lgChiantiOn)
	{
		TotalInsanity();
	}
	if( atmdat.lgStoutHybrid && !atmdat.lgStoutOn )
	{
		TotalInsanity();
	}

	/* read in atomic and molecular models from third-party databases */
	if( atmdat.lgLamdaOn || atmdat.lgChiantiOn || atmdat.lgStoutOn)
		database_readin();
	else
		nSpecies = 0;

	/* initialize the large block of level 1 real lines, and OP level 2 lines */
	lines_setup();

	/* mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat mewe_gbar.dat ========*/
	/* read in g-bar data taken from
	 *>>refer	all	gbar	Mewe, R., Gronenschild, E. H. B. M., van den Oord, G. H. J. 1985, A&AS, 62, 197 */
	/* open file with Mewe coefficients */

	if( trace.lgTrace )
		fprintf( ioQQQ," atmdat_readin reading mewe_gbar.dat\n");

	ioDATA = open_data( "mewe_gbar.dat", "r" );

	/* get magic number off first line */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " mewe_gbar.dat error getting magic number\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* check the number */
	sscanf( chLine , "%ld" , &magic1 );
	if( magic1 != 9101 )
	{
		fprintf( ioQQQ, " mewe_gbar.dat starts with wrong magic number=%ld \n", 
		  magic1 );
		cdEXIT(EXIT_FAILURE);
	}

	/* now get the actual data, indices are correct for c, in Fort went from 2 to 210 */
	for( i=1; i < 210; i++ )
	{
		if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " mewe_gbar.dat error getting line  %li\n", i );
			cdEXIT(EXIT_FAILURE);
		}
		/*printf("%s\n",chLine);*/
		double help[4];
		sscanf( chLine, "%lf %lf %lf %lf ", &help[0], &help[1], &help[2], &help[3] );
		for( int l=0; l < 4; ++l )
			MeweCoef.g[i][l] = (realnum)help[l];
	}

	/* get magic number off last line */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " mewe_gbar.dat error getting last magic number\n" );
		cdEXIT(EXIT_FAILURE);
	}

	sscanf( chLine , "%ld" , &magic2 );

	if( magic1 != magic2 )
	{
		fprintf( ioQQQ, " mewe_gbar.dat ends will wrong magic number=%ld \n", 
		  magic2 );
		cdEXIT(EXIT_FAILURE);
	}

	fclose( ioDATA );

	if( trace.lgTrace )
		fprintf( ioQQQ," reading mewe_gbar.dat OK \n");

	 /* This is what remains of the t_yield initialization
	  * this should not be in the constructor of t_yield
	  * since it initializes a different struct */
	for( nelem=0; nelem < LIMELM; nelem++ )
		for( ion=0; ion < LIMELM; ion++ )
			Heavy.nsShells[nelem][ion] = LONG_MAX;

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with #
	 * this loop goes from lithium to Zn */
	for( nelem=2; nelem < LIMELM; nelem++ )
	{
		/* nelem is on the shifted C scale, so 2 is Li */
		for( ion=0; ion <= nelem; ion++ )
		{
			/* number of bound electrons, = atomic number for neutral */
			nelec = nelem - ion + 1;
			/* one of dima's routines to determine the number of electrons
			 * for this species, nelem +1 to shift to physical number */
			/* subroutine atmdat_outer_shell(iz,in,imax,ig0,ig1)
			 * iz - atomic number from 1 to 30 (integer) 
			 * in - number of electrons from 1 to iz (integer)
			 * Output: imax - number of the outer shell
			 */
			atmdat_outer_shell(nelem+1,nelec,&imax,&ig0,&ig1);

			ASSERT( imax > 0 && imax <= 10 );

			/* nsShells[nelem][ion] is outer shell number for ion with nelec electrons
			 * on physics scale, with K shell being 1 */
			Heavy.nsShells[nelem][ion] = imax;
		}
	}

	/*************************************************************
	 *                                                           *
	 * get the Auger electron yield data set                     *
	 *                                                           *
	 *************************************************************/

	t_yield::Inst().init_yield();

	/****************************************************************
	 *                                                              *
	 * get the Hummer and Storey model case A, B results, these are *
	 * the two data files e1b.dat and e2b.dat, for H and He         *
	 *                                                              *
	 ****************************************************************/

	/* now get Hummer and Storey case b data, loop over H, He */
	/* >>chng 01 aug 08, generalized this to both case A and B, and all elements
	 * up to HS_NZ, now 8 */
	for( ipZ=0; ipZ<HS_NZ; ++ipZ )
	{
		/* don't do the minor elements, Li-B */
		if( ipZ>1 && ipZ<5 ) continue;

		for( iCase=0; iCase<2; ++iCase )
		{
			/* open file with case B data */
			/* >>chng 01 aug 08, add HS_ to start of file names to indicate Hummer Storey origin */
			/* create file name for this charge
			 * first character after e is charge number for this element,
			 * then follows whether this is case A or case B data */
			sprintf( chFilename, "HS_e%ld%c.dat", ipZ+1, ( iCase == 0 ) ? 'a' : 'b' );

			if( trace.lgTrace )
				fprintf( ioQQQ," atmdat_readin reading Hummer Storey emission file %s\n",chFilename );

			/* open the file */
			ioDATA = open_data( chFilename, "r" );

			/* read in the number of temperatures and densities*/
			i = fscanf( ioDATA, "%li %li ", 
				&atmdat.ntemp[iCase][ipZ], &atmdat.nDensity[iCase][ipZ] );

			/* check that ntemp and nDensity are below NHSDIM, 
			 * set to 15 in atmdat_HS_caseB.h */
			assert (atmdat.ntemp[iCase][ipZ] >0 && atmdat.ntemp[iCase][ipZ] <= NHSDIM );
			assert (atmdat.nDensity[iCase][ipZ] > 0 && atmdat.nDensity[iCase][ipZ] <= NHSDIM);

			/* loop reading in line emissivities for all temperatures*/
			for( ipTemp=0; ipTemp < atmdat.ntemp[iCase][ipZ]; ipTemp++ )
			{
				for( ipDens=0; ipDens < atmdat.nDensity[iCase][ipZ]; ipDens++ )
				{
					long int junk, junk2 , ne;
					i = fscanf( ioDATA, " %lf %li %lf %c %li %ld ", 
						&atmdat.Density[iCase][ipZ][ipDens], &junk , 
						&atmdat.ElecTemp[iCase][ipZ][ipTemp], &cha , &junk2 , 
						&atmdat.ncut[iCase][ipZ] );

					ne = atmdat.ncut[iCase][ipZ]*(atmdat.ncut[iCase][ipZ] - 1)/2;
					ASSERT( ne<=NLINEHS );
					for( j=0; j < ne; j++ )
					{
						i = fscanf( ioDATA, "%lf ",
							    &atmdat.Emiss[iCase][ipZ][ipTemp][ipDens][j] );
					}
				}
			}

			/*this is end of read-in loop */
			fclose(ioDATA); 
			if( trace.lgTrace )
				fprintf( ioQQQ," reading %s OK\n", chFilename );

#			if 0
			/* print list of densities and temperatures */
			for( ipDens=0; ipDens<atmdat.nDensity[iCase][ipZ]; ipDens++ )
			{
				fprintf(ioQQQ," %e,", atmdat.Density[iCase][ipZ][ipDens]);
			}
			fprintf(ioQQQ,"\n");
			for( ipTemp=0; ipTemp<atmdat.ntemp[iCase][ipZ]; ipTemp++ )
			{
				fprintf(ioQQQ," %e,", atmdat.ElecTemp[iCase][ipZ][ipTemp]);
			}
			fprintf(ioQQQ,"\n");
#			endif
		}
	}

	// read cross sections for neutral helium
	read_SH98_He1_cross_sections();
	// read cross sections for some he-like ions
	read_Helike_cross_sections();

	/* Verner's model FeII atom
	 * do not call if Netzer model used, or it Fe(2) is zero
	 * exception is when code is searching for first soln */
	/* read the atomic data from files */
	/* convert line arrays to internal form needed by code */
	FeIICreate();

	// set up spline coefficients for two-photon continuua
	atmdat_2phot_setSplineCoefs();

	return;
}

STATIC void read_SH98_He1_cross_sections(void)
{
	DEBUG_ENTRY( "read_SH98_He1_cross_sections()" );

	FILE *ioDATA;
	bool lgEOL;

	char chPath[FILENAME_PATH_LENGTH_2],
		chDirectory[FILENAME_PATH_LENGTH_2], 
		chLine[FILENAME_PATH_LENGTH_2];

	const int ipNUM_FILES = 10;

	char chFileNames[ipNUM_FILES][10] = {
		"p0202.3se",
		"p0202.3po",
		"p0202.3ge",
		"p0202.3fo",
		"p0202.3de",
		"p0202.1se",
		"p0202.1po",
		"p0202.1ge",
		"p0202.1fo",
		"p0202.1de" };

	HS_He1_Xsectn = ((double****)MALLOC( (size_t)(26)*sizeof(double***)));
	HS_He1_Energy = ((double****)MALLOC( (size_t)(26)*sizeof(double***)));

	// point these to NULL and use real quantum numbers rather than starting at 0
	HS_He1_Xsectn[0] = NULL;
	HS_He1_Energy[0] = NULL;

	for( long in = 1; in<=25; in++ )
	{
		// malloc n values of angular momentum, but not more than 5
		HS_He1_Xsectn[in] = ((double***)MALLOC( (size_t)(MIN2(5,in))*sizeof(double**)));
		HS_He1_Energy[in] = ((double***)MALLOC( (size_t)(MIN2(5,in))*sizeof(double**)));
		for( long il = 0; il<MIN2(5,in); il++ )
		{
			// malloc two values of spin
			HS_He1_Xsectn[in][il] = ((double**)MALLOC( (size_t)(2)*sizeof(double*)));
			HS_He1_Energy[in][il] = ((double**)MALLOC( (size_t)(2)*sizeof(double*)));
			HS_He1_Xsectn[in][il][0] = ((double*)MALLOC( (size_t)(NUM_HS98_DATA_POINTS)*sizeof(double)));
			HS_He1_Energy[in][il][0] = ((double*)MALLOC( (size_t)(NUM_HS98_DATA_POINTS)*sizeof(double)));
			HS_He1_Xsectn[in][il][1] = ((double*)MALLOC( (size_t)(NUM_HS98_DATA_POINTS)*sizeof(double)));
			HS_He1_Energy[in][il][1] = ((double*)MALLOC( (size_t)(NUM_HS98_DATA_POINTS)*sizeof(double)));
		}
	}

#	ifdef _WIN32
	strcpy( chDirectory, "sh98_he1\\pi\\" );
#	else
	strcpy( chDirectory, "sh98_he1/pi/" );
#	endif
	
	//HS_He1_data[25][<=4][2]

	for( long ipFile=0; ipFile<ipNUM_FILES; ipFile++ )
	{
		long S, L, index, N=0;
		long UNUSED P;

		strcpy( chPath, chDirectory );
		strcat( chPath, chFileNames[ipFile] );
		ioDATA = open_data( chPath, "r" );
		
		while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
		{
			long i=1, s;
			long i1, i2, i3;
			long numDataPoints;
			
			// first line (read above in while) is not needed except that "0 0 0" marks EOF
			i1 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			i2 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			i3 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			if( i1==0 && i2==0 && i3==0 )
				break;

			// don't need next two lines in each set
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
			
			i=1;
			// 4th line in each set identifies the quantum level
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
			S = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			L = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			P = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			index = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );

			//indices start with unity
			ASSERT( index >= 1 );

			// index is energy order in that series and is therefore related
			// to principal quantum number.  For triplet S must add one because
			// series starts at n=2. 
			if( L==0 && S==3 )
				N = L + index + 1;
			else
				N = L + index;

			// data go up to n=25
			ASSERT( N<=25 );

			if( S==1 )
				s=0;
			else if( S==3 )
				s=1;
			else
				TotalInsanity();

			i=1;
			// 5th line in each set contains the number of energies
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
			//first value is not needed
			FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			numDataPoints = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			// each set has exactly 811 data points, might as well assert this
			ASSERT( numDataPoints == NUM_HS98_DATA_POINTS );

			// don't need 6th line either
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );

			// now begin reading in lines
			// there must be exactly numDataPoints ordered pairs,
			// throw an assert for any deviation
			for( long k=0; k<NUM_HS98_DATA_POINTS; k++ )
			{
				i=1;
				read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
				HS_He1_Energy[N][L][s][k] = (double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
				HS_He1_Xsectn[N][L][s][k] = (double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			}
		}

		// we reached end of file, assert last quantum number was 25
		ASSERT( N==25 );

		fclose( ioDATA );
	}

	return;
}

STATIC void read_Helike_cross_sections(void)
{
	DEBUG_ENTRY( "read_Helike_cross_sections()" );

	FILE *ioDATA;
	bool lgEOL;

	char chLine[FILENAME_PATH_LENGTH_2];
	char chFileName[23] = "helike_pcs_topbase.dat";

	// the data only go as high as n=10	
	const int MaxN = 10;
	long last_i1=0;

	// data will be used up to calcium (nelem=19)
	// data exists up to iron but we will not use it because it has strong resonance
	// features, but is nonetheless very nearly hydrogenic
	OP_Helike_Xsectn = ((double*****)MALLOC( (size_t)(ipCALCIUM+1)*sizeof(double****)));
	OP_Helike_Energy = ((double*****)MALLOC( (size_t)(ipCALCIUM+1)*sizeof(double****)));
	OP_Helike_NumPts = ((long****)MALLOC( (size_t)(ipCALCIUM+1)*sizeof(long***)));

	// point these to NULL because we will never use them
	OP_Helike_Xsectn[ipHYDROGEN] = NULL;
	OP_Helike_Energy[ipHYDROGEN] = NULL;
	OP_Helike_NumPts[ipHYDROGEN] = 0;
	OP_Helike_Xsectn[ipHELIUM] = NULL;
	OP_Helike_Energy[ipHELIUM] = NULL;
	OP_Helike_NumPts[ipHELIUM] = 0; 

	for( long nelem = ipLITHIUM; nelem<= ipCALCIUM; nelem++ )
	{
		// malloc principal quantum number
		OP_Helike_Xsectn[nelem] = ((double****)MALLOC( (size_t)(MaxN+1)*sizeof(double***)));
		OP_Helike_Energy[nelem] = ((double****)MALLOC( (size_t)(MaxN+1)*sizeof(double***)));
		OP_Helike_NumPts[nelem] = ((long***)MALLOC( (size_t)(MaxN+1)*sizeof(long**)));

		for( long in = 1; in<=MaxN; in++ )
		{
			// malloc angular momentum
			OP_Helike_Xsectn[nelem][in] = ((double***)MALLOC( (size_t)(in)*sizeof(double**)));
			OP_Helike_Energy[nelem][in] = ((double***)MALLOC( (size_t)(in)*sizeof(double**)));
			OP_Helike_NumPts[nelem][in] = ((long**)MALLOC( (size_t)(in)*sizeof(long*)));
			for( long il = 0; il<in; il++ )
			{
				// malloc two values of spin
				OP_Helike_Xsectn[nelem][in][il] = ((double**)MALLOC( (size_t)(2)*sizeof(double*)));
				OP_Helike_Energy[nelem][in][il] = ((double**)MALLOC( (size_t)(2)*sizeof(double*)));
				OP_Helike_NumPts[nelem][in][il] = ((long*)MALLOC( (size_t)(2)*sizeof(long)));
				// point these to NULL for now, won't know how many we need until we begin parsing
				OP_Helike_Xsectn[nelem][in][il][0] = NULL;
				OP_Helike_Energy[nelem][in][il][0] = NULL;
				OP_Helike_NumPts[nelem][in][il][0] = 0;
				OP_Helike_Xsectn[nelem][in][il][1] = NULL;
				OP_Helike_Energy[nelem][in][il][1] = NULL;
				OP_Helike_NumPts[nelem][in][il][1] = 0;
			}
		}
	}

	ioDATA = open_data( chFileName, "r" );

	// Header looks like this:
	// ================================================
	//       I  NZ  NE  ISLP  ILV        E(RYD)      NP
	// ================================================
	//      1   3   2   100    1  -5.53159E+00      55

	// so skip the first three lines.
	for( long i=0; i<3; i++)
	{
		if(	read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
		{
			fprintf( ioQQQ,"PROBLEM corruption in TOPbase Helike pcs datafile.\nSorry\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		long i=1;
		long n, l, s;
		long i1, i2, i3, i4, i5, i7;
		double i6;
		
		i1 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i2 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i3 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i4 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i5 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i6 = (double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
		i7 = (long)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );

		if( lgEOL )
		{
			fprintf( ioQQQ,"PROBLEM corruption in TOPbase Helike pcs datafile.\nSorry\n" );
			cdEXIT(EXIT_FAILURE);
		}

		// this marks end of data.
		if( i1==i2 && i1==i3 && i1==i4 && i1==i5 && i1==i7 && i1==-1 )
			break;

		// first parameter is level index (overall, not just for series or charge)
		// check that it is as expected
		ASSERT( i1>0 && i1==(last_i1 + 1) && i1<=795 ); 
		last_i1 = i1;
		// second parameter is nuclear charge
		ASSERT( (i2-1)<=ipCALCIUM && (i2-1)>=ipLITHIUM );
		// third parameter must be 2 for helike.
		ASSERT( i3==2 );
		// fourth parameter is (2S+1)*100+L*10+P 
		ASSERT( i4>=100 && i4<400 );
		if( i4 >= 300 )
			s=1;
		else
			s=0;
		l = (i4 - (2*s+1)*100)/10;
		// data only goes up to l=2
		ASSERT( l<=2 );
		// fifth is index in the series, related to principal quantum number
		ASSERT( i5>=1 && i5<=10 );
		if( s==1 && l==0 )
			n = i5 + 1;
		else
			n = i5 + l;
		ASSERT( l<=MaxN );
		// sixth is threshhold energy, don't need but assert negative
		// \todo 3 save this and renorm cross-section with ratio of actual to recorded Eth?
		ASSERT( i6 < 0. );
		// seventh parameter is number of data points, can be zero, use to MALLOC otherwise
		OP_Helike_NumPts[i2-1][n][l][s] = i7;
		if( i7==0 )
			continue;
		
		ASSERT( i7 > 0 );
		ASSERT( OP_Helike_Xsectn[i2-1][n][l][s]==NULL );
		ASSERT( OP_Helike_Energy[i2-1][n][l][s]==NULL );
		OP_Helike_Xsectn[i2-1][n][l][s] = ((double*)MALLOC( (size_t)(i7)*sizeof(double)));
		OP_Helike_Energy[i2-1][n][l][s] = ((double*)MALLOC( (size_t)(i7)*sizeof(double)));

		// now begin reading in lines
		// there must be exactly i7 ordered pairs,
		for( long k=0; k<i7; k++ )
		{
			i=1;
			read_whole_line( chLine , (int)sizeof(chLine) , ioDATA );
			OP_Helike_Energy[i2-1][n][l][s][k] = (double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			OP_Helike_Xsectn[i2-1][n][l][s][k] = (double)FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );

			// make sure data is well-behaved
			if( k > 0 )
			{
				ASSERT( OP_Helike_Energy[i2-1][n][l][s][k] > OP_Helike_Energy[i2-1][n][l][s][k-1] );
				ASSERT( OP_Helike_Xsectn[i2-1][n][l][s][k] < OP_Helike_Xsectn[i2-1][n][l][s][k-1] );
			}

			// try to read one more item off the line and verify that lgEOL is true
			FFmtRead( chLine, &i, sizeof(chLine), &lgEOL );
			ASSERT( lgEOL );
		}
	}

	fclose( ioDATA );

	return;
}

t_yield::t_yield()
{
	 /* preset two arrays that will hold auger data 
	  * set to very large values so that code will blow
	  * if not set properly below */
	for( int nelem=0; nelem < LIMELM; nelem++ )
	{
		for( int ion=0; ion < LIMELM; ion++ )
		{
			for( int ns=0; ns < 7; ns++ )
			{
				n_elec_eject[nelem][ion][ns] = LONG_MAX;
				for( int nelec=0; nelec < 10; nelec++ )
				{
					frac_elec_eject[nelem][ion][ns][nelec] = FLT_MAX;
				}
			}
		}
	}

	lgKillAuger = false;
}

void t_yield::init_yield()
{
	char chLine[FILENAME_PATH_LENGTH_2];
	const char* chFilename;

	/* following is double for sscanf to work */
	double temp[15];

	DEBUG_ENTRY( "init_yield()" );

	/*************************************************************
	 *                                                           *
	 * get the Auger electron yield data set                     *
	 *                                                           *
	 *************************************************************/

	/* NB NB -- This test of Heavy.nsShells remains here to assure
	 * that it contains meaningful values since needed below, once
	 * t_Heavy is a Singleton, it can be removed !!! */
	ASSERT( Heavy.nsShells[2][0] > 0 );

	/* hydrogen and helium will not be done below, so set yields here*/
	n_elec_eject[ipHYDROGEN][0][0] = 1;
	n_elec_eject[ipHELIUM][0][0] = 1;
	n_elec_eject[ipHELIUM][1][0] = 1;

	frac_elec_eject[ipHYDROGEN][0][0][0] = 1;
	frac_elec_eject[ipHELIUM][0][0][0] = 1;
	frac_elec_eject[ipHELIUM][1][0][0] = 1;

	/* open file auger.dat, yield data file that came from
	 * >>refer	all	auger	Kaastra, J. S., & Mewe, R. 1993, A&AS, 97, 443-482 */
	chFilename = "mewe_nelectron.dat";

	if( trace.lgTrace )
		fprintf( ioQQQ, " init_yield reading %s\n", chFilename );

	FILE *ioDATA;

	/* open the file */
	ioDATA = open_data( chFilename, "r" );

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with #
	 * this loop goes from lithium to Zn */
	for( int nelem=2; nelem < LIMELM; nelem++ )
	{
		/* nelem is on the shifted C scale, so 2 is Li */
		for( int ion=0; ion <= nelem; ion++ )
		{
			for( int ns=0; ns < Heavy.nsShells[nelem][ion]; ns++ )
			{
				char ch1 = '#';
				/* the * is the old comment char, accept it, but really want # */
				while( ch1 == '#' || ch1 == '*' )
				{
					if( read_whole_line( chLine, (int)sizeof(chLine), ioDATA ) == NULL )
					{
						fprintf( ioQQQ, " %s error getting line %i\n", chFilename, ns );
						cdEXIT(EXIT_FAILURE);
					}
					ch1 = chLine[0];
				}
				/*printf("string is %s\n",chLine );*/
				sscanf( chLine, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
					&temp[0], &temp[1], &temp[2], &temp[3], &temp[4],
					&temp[5], &temp[6], &temp[7], &temp[8], &temp[9],
					&temp[10],&temp[11],&temp[12],&temp[13],&temp[14] );
				n_elec_eject[nelem][ion][ns] = (long int)temp[3];

				ASSERT( n_elec_eject[nelem][ion][ns] >= 0 && n_elec_eject[nelem][ion][ns] < 11 );

				/* can pop off up to 10 auger electrons, these are the probabilities*/
				for( int j=0; j < 10; j++ )
				{
					frac_elec_eject[nelem][ion][ns][j] = (realnum)temp[j+4];
					ASSERT( frac_elec_eject[nelem][ion][ns][j] >= 0. );
				}
			}
		}
		/* activate this print statement to get yield for k-shell of all atoms */
		/*fprintf(ioQQQ,"yyyield\t%li\t%.4e\n", nelem+1, frac_elec_eject[nelem][0][0][0] );*/
	}

	fclose( ioDATA );

	if( trace.lgTrace )
	{
		/* this is set with no auger command */
		if( lgKillAuger )
			fprintf( ioQQQ, " Auger yields will be killed.\n");
		fprintf( ioQQQ, " reading %s OK\n", chFilename );
	}

	/* open file mewe_fluor.dat, yield data file that came from
	 * >>refer	all	auger	Kaastra, J. S., & Mewe, R. 1993, 97, 443-482 */
	chFilename = "mewe_fluor.dat";

	if( trace.lgTrace )
		fprintf( ioQQQ, " init_yield reading %s\n", chFilename );

	/* open the file */
	ioDATA = open_data( chFilename, "r" );

	/* now read in all auger yields
	 * will do elements from li on up,
	 * skip any line starting with # */
	do 
	{
		if( read_whole_line( chLine, (int)sizeof(chLine), ioDATA ) == NULL )
		{
			fprintf( ioQQQ, " %s error getting line %i\n", chFilename, 0 );
			cdEXIT(EXIT_FAILURE);
		}
	}
	while( chLine[0] == '#' );

	bool lgEOL = false;

	nfl_lines = 0;
	do
	{
		const int NKM = 10;
		int nDima[NKM] = { 0, 1, 2, 2, 3, 4, 4, 5, 5, 6 };
		int nAuger;

		if( nfl_lines >= MEWE_FLUOR )
			TotalInsanity();

		/*printf("string is %s\n",chLine );*/
		sscanf( chLine, "%lf %lf %lf %lf %lf %lf %lf", 
			&temp[0], &temp[1], &temp[2], &temp[3], &temp[4], 
			&temp[5], &temp[6] );

		/* the atomic number, C is 5 */
		nfl_nelem[nfl_lines] = (int)temp[0]-1;
		ASSERT( nfl_nelem[nfl_lines] >= 0 && nfl_nelem[nfl_lines] < LIMELM );

		/* the ion stage for target, atom is 0 */
		nfl_ion[nfl_lines] = (int)temp[1]-1;
		ASSERT( nfl_ion[nfl_lines] >= 0 && nfl_ion[nfl_lines] <= nfl_nelem[nfl_lines]+1 );

		/* the target's shell */
		nfl_nshell[nfl_lines] = nDima[(long)temp[2]-1];
		ASSERT( nfl_nshell[nfl_lines] >= 0 && 
			/* nsShells is shell number, where K is 1 */
			nfl_nshell[nfl_lines] < Heavy.nsShells[nfl_nelem[nfl_lines]][nfl_ion[nfl_lines]]-1 );

		/* this is the number of Auger electrons ejected */
		nAuger = (int)temp[3];
		/* so this is the spectrum of the photons that are emitted */
		nfl_ion_emit[nfl_lines] = nfl_ion[nfl_lines] + nAuger + 1;
		/* must be gt 0 since at least photoelectron comes off */
		ASSERT( nfl_ion_emit[nfl_lines] > 0 && nfl_ion_emit[nfl_lines] <= nfl_nelem[nfl_lines]+1);

		/* this is the type of line as defined in their paper */
		nfl_nLine[nfl_lines] = (int)temp[4];

		/* energy in Ryd */
		fl_energy[nfl_lines] = (realnum)temp[5] / (realnum)EVRYD;
		ASSERT( fl_energy[nfl_lines] > 0. );

		/* fluor yield */
		fl_yield[nfl_lines] = (realnum)temp[6];
		/* NB cannot assert <=1 since data file has yields around 1.3 - 1.4 */
		ASSERT( fl_yield[nfl_lines] >= 0 );

		++nfl_lines;

		do
		{
			if( read_whole_line( chLine, (int)sizeof(chLine), ioDATA ) == NULL )
				lgEOL = true;
		} 
		while( chLine[0]=='#' && !lgEOL );
	} 
	while( !lgEOL );

	fclose( ioDATA );
	if( trace.lgTrace )
		fprintf( ioQQQ, " reading %s OK\n", chFilename );
}

// read autoionization data from Badnell data file 
STATIC void ReadBadnellAIData(const string& fnam,      // filename containing the Badnell data
			      long nelem,              // nelem is on C scale
			      long ion,                // ion is on C scale
			      TransitionList& UTA,     // UTA lines will be pushed on this stack
			      bitset<IS_TOP> Skip)     // option to skip transitions from a particular shell
{
	DEBUG_ENTRY( "ReadBadnellAIData()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," ReadBadnellAIData reading %s\n", fnam.c_str() );

	fstream ioDATA;
	open_data( ioDATA, fnam.c_str(), mode_r );

	string line;
	getline( ioDATA, line );
	ASSERT( line.substr(0,4) == "SEQ=" );
	getline( ioDATA, line );
	getline( ioDATA, line );
	// we don't need the parent level data, so we will skip it...
	ASSERT( line.substr(3,21) == "PARENT LEVEL INDEXING" );
	int nParent;
	istringstream iss( line.substr(65,4) );
	iss >> nParent;
	// data lines containing data for all levels of the parent ion will span nMulti lines
	int nMulti = (nParent+5)/6;
	for( int i=0; i < nParent+5; ++i )
		getline( ioDATA, line );

	// here starts the header for the level data we need
	ASSERT( line.substr(3,26) == "IC RESOLVED LEVEL INDEXING" );
	int nLevel;
	istringstream iss2( line.substr(63,6) );
	iss2 >> nLevel;
	// skip rest of header
	for( int i=0; i < 3; ++i )
		getline( ioDATA, line );

	// now get the level data
	vector<t_BadnellLevel> level( nLevel );
	for( int i=0; i < nLevel; ++i )
	{
		getline( ioDATA, line );
		istringstream iss3( line );
		int indx, irsl;
		iss3 >> indx >> irsl;
		level[indx-1].irsl = irsl;
		level[indx-1].config = line.substr(16,20);
		istringstream iss4( line.substr(37,1) );
		iss4 >> level[indx-1].S;
		istringstream iss5( line.substr(39,1) );
		iss5 >> level[indx-1].L;
		istringstream iss6( line.substr(41,4) );
		double J;
		iss6 >> J;
		level[indx-1].g = nint(2.*J + 1.);
		istringstream iss7( line.substr(46,11) );
		iss7 >> level[indx-1].energy;

		// which inner shell has been excited?
		level[indx-1].lgAutoIonizing = ( line[57] == '*' );
		if( level[indx-1].lgAutoIonizing )
		{
			if( level[indx-1].config.find( "1S1" ) != string::npos )
				level[indx-1].WhichShell = IS_K_SHELL;
			else if( level[indx-1].config.find( "2S1" ) != string::npos )
				level[indx-1].WhichShell = IS_L1_SHELL;
			else if( level[indx-1].config.find( "2P5" ) != string::npos )
				level[indx-1].WhichShell = IS_L2_SHELL;
			else
				TotalInsanity();
		}
		else
		{
			level[indx-1].WhichShell = IS_NONE;
		}
	}

	// levels are done, now move on to the lines
	// first search for start of the header
	while( getline( ioDATA, line ) )
	{
		if( line.find( "IRSL  IRSL" ) != string::npos )
			break;
	}
	// skip rest of the header
	for( int i=0; i < nMulti-1; ++i )
		getline( ioDATA, line );

	// blank line that will be pushed on the UTA line stack
	qList BlankStates(1);
	TransitionList BlankList("BlankList",&BlankStates,1);
	TransitionList::iterator BlankLine = BlankList.begin();
	(*BlankLine).Junk();

	// start reading the line data
	while( getline( ioDATA, line ) )
	{
		// have we reached the end of the line data?
		if( line.size() < 10 )
			break;

		// test if there is an autoionization rate on this line
		// if not, it only contains radiative data and we skip it...
		if( line.size() < 50 )
			continue;

		// there may be an asterisk here; wipe it out since we don't need it
		line[19] = ' ';

		int irsl_lo, irsl_hi, dum;
		double edif, Bij, Rji, Aai;
		istringstream iss8( line );
		// irsl_lo: index for lower level of transition
		// irsl_hi: index for upper level of transition
		// edif: energy difference between levels in Ryd
		// Bij: UPWARD Einstein A for transition irsl_lo -> irsl_hi
		// Rji: sum of Aji for all radiative transitions to lower levels
		// Aai: autoionization rate from level irsl_hi
		iss8 >> irsl_lo >> irsl_hi >> dum >> dum >> edif >> Bij >> Rji >> Aai;
		// ind_lo and ind_hi are on C scale
		int ind_lo = irsl2ind( level, irsl_lo );
		int ind_hi = irsl2ind( level, irsl_hi );
		ASSERT( level[ind_hi].lgAutoIonizing );
		// skip rest of partial autoionization rates
		for( int i=0; i < nMulti-1; ++i )
			getline( ioDATA, line );
		// skip this transition if it does not originate from ground level
		// or if the user requested to skip excitations from this inner shell
		if( ind_lo == 0 && !Skip[level[ind_hi].WhichShell] )
		{
			UTA.push_back( *BlankLine );
			InitTransition( UTA.back() );

			// t_emission has nelem and ion on fortran scale...
			(*UTA.back().Hi()).nelem() = nelem+1;
			(*UTA.back().Hi()).IonStg() = ion+1;

			(*UTA.back().Hi()).g() = (realnum)level[ind_hi].g;
			(*UTA.back().Lo()).g() = (realnum)level[ind_lo].g;

			double WavNum = edif*RYD_INF;

			/* wavelength in Angstroms */
			UTA.back().WLAng() = (realnum)(1e8/WavNum);
			UTA.back().EnergyWN() = (realnum)WavNum;

			/* store branching ratio for autoionization */
			double frac_ioniz = Aai/(Rji + Aai);
			ASSERT( frac_ioniz >= 0. &&  frac_ioniz <= 1. );
			UTA.back().Emis().AutoIonizFrac() = (realnum)frac_ioniz;

			/* this is true spontaneous rate for doubly excited state to ground
			 * and is used for pumping, and also relaxing back to inner shell */
			/* Badnell gives UPWARD transition rate Alu, an unusual notation,
			 * convert it here to the normal downward transition rate Aul */
			UTA.back().Emis().Aul() = (realnum)(Bij*(*UTA.back().Lo()).g()/(*UTA.back().Hi()).g());
			UTA.back().Emis().gf() =
				(realnum)GetGF( UTA.back().Emis().Aul(), UTA.back().EnergyWN(), (*UTA.back().Hi()).g() );

			UTA.back().Emis().iRedisFun() = ipPRD;

			UTA.back().Emis().dampXvel() = (realnum)((Rji+Aai)/UTA.back().EnergyWN()/PI4);

			ASSERT( UTA.back().Emis().dampXvel() > 0. );

			// remove this line if it is too weak
			if( UTA.back().Emis().gf() < level[ind_lo].g * f_cutoff )
				UTA.pop_back();
		}
	}

	// perform a sanity check on the tail of the file
	getline( ioDATA, line );
	ASSERT( line.substr(3,7) == "NRSLMX=" );

	ioDATA.close();

	if( trace.lgTrace )
		fprintf( ioQQQ, " reading %s OK\n", fnam.c_str() );
}

inline void InitTransition(const TransitionProxy& t)
{
	t.AddHiState();
	t.AddLoState();
	t.AddLine2Stack();
}

inline int irsl2ind(vector<t_BadnellLevel>& level, int irsl)
{
	for( unsigned int i=0; i < level.size(); ++i )
	{
		if( level[i].irsl == irsl )
			return (int)i;
	}
	// we should never get here...
	TotalInsanity();
}
