/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseCommands main command line parser, decode command, then call other routines to read */

/* CHANGES: (M. Salz 17.05.2013)
 *  - switch dense_tabden() to
      interpol_tabulated(), which is basically the same but not
      density specific, and does also temperature and velocity
      interpolation */

#include "cddefines.h"
#include "physconst.h"
#include "optimize.h"
#include "doppvel.h"
#include "stopcalc.h"
#include "abund.h"
#include "geometry.h"
#include "dense.h"
#include "plot.h"
#include "grid.h"
#include "rfield.h"
#include "grainvar.h"
#include "dynamics.h"
#include "magnetic.h"
#include "trace.h"
#include "atmdat.h"
#include "h2.h"
#include "mole.h"
#include "hmi.h"
#include "rt.h"
#include "thermal.h"
#include "opacity.h"
#include "atomfeii.h"
#include "called.h"
#include "wind.h"
#include "hextra.h"
#include "iterations.h"
#include "radius.h"
#include "input.h" 
#include "monitor_results.h"
#include "phycon.h"
#include "fudgec.h"
#include "version.h"
#include "conv.h"
#include "parse.h"
#include "cosmology.h"
#include "pressure.h"
#include "parser.h"
#include "dark_matter.h"

void ParseAperture(Parser &p);
void ParseAtom(Parser &p);
void ParseBremsstrahlung(Parser &p);
void ParseCExtra(Parser &p);
void ParseCMBOuter(Parser &p);
void ParseCosm(Parser &p);
void ParseCovering(Parser &p);
void ParseCylinder(Parser &p);
void ParseDarkMatter(Parser &p);
void ParseDielectronic(Parser &);
void ParseDiffuse(Parser &p);
void ParseDistance(Parser &p);
void ParseDoubleTau(Parser &);
void ParseEden(Parser &p);
void ParseEnergy(Parser &p);
void ParseFail(Parser &p);
void ParseFill(Parser &p);
void ParseF_nuSpecific(Parser &p);
void ParseForceTemperature(Parser &p);
void ParseFudge(Parser &p);
void ParsePGrains(Parser &);
void ParseGravity(Parser &p);
void ParseHeLike(Parser &);
void ParseHelp(Parser &);
void ParseHExtra(Parser &p);
void ParseConvHighT(Parser &);
void ParseHydrogen(Parser &);
void ParseInitCount(Parser &p);
void ParseIntensity(Parser &p);
void ParseIterations(Parser &p);
void ParseL_nu(Parser &p);
void ParseLaser(Parser &p);
void ParseLuminosity(Parser &p);
void ParseNeutrons(Parser &p);
void ParseNuF_nu(Parser &p);
void ParseNuL_nu(Parser &p);
void ParsePhi(Parser &p);
void ParseQH(Parser &p);
void ParseRoberto(Parser &);
void ParseSpecial(Parser &);
void ParseTauMin(Parser &p);
void ParseTitle(Parser &);
void ParseTolerance(Parser &);
void ParseVLaw(Parser &p);
void ParseTurbulence(Parser &p);

void ParseCommands(void)
{
	bool lgStop ,
		lgStop_not_enough_info;

	DEBUG_ENTRY( "ParseCommands()" );

	/* following says abundances are solar  */
	abund.lgAbnSolar = true;

	/* there are no plots desired yet */
	plotCom.lgPlotON = false;
	plotCom.nplot = 0;

	/* this flag remembers whether grains have ever been turned on.  It is passed
	 * to routine ParseAbundances - there grains will be turned on with commands
	 * such as abundances ism, unless grains were previously set 
	 * with a grains command.  this way abundances will not clobber explicitly set
	 * grains. */

	radius.Radius = 0.;
	radius.rinner = 0.;
	rfield.nShape = 0;

	/* will be set true if no buffering command entered */
	input.lgSetNoBuffering = false;

	/* initialize some code variables in case assert command used in input stream */
	InitMonitorResults();

	for( long int i=0; i < LIMSPC; i++ )
	{
		strcpy( rfield.chRSpec[i], "UNKN" );
	}
	optimize.nparm = 0;

	/* this is an option to turn on trace printout on the nth
	 * call from the optimizer */
	if( optimize.lgTrOpt )
	{
		/* nTrOpt was set with "optimize trace" command,
		 * is iteration to turn on trace */
		if( optimize.nTrOpt == optimize.nOptimiz )
		{
			trace.lgTrace = true;
			called.lgTalk = cpu.i().lgMPI_talk();
			trace.lgTrOvrd = true;
			fprintf( ioQQQ, " READR turns on trace from optimize option.\n" );
		}
	}

	/* say this is a beta version if we are talking and it is the truth */
	if( t_version::Inst().nBetaVer > 0 && called.lgTalk )
	{
		fprintf( ioQQQ, 
			 "\n                               This is a beta release of Cloudy, and is intended for testing only.\n" );
		fprintf( ioQQQ,
			 "Please help make Cloudy better by posing problems or suggestions on http://tech.groups.yahoo.com/group/cloudy_simulations/.\n\n" );
	}

	if( called.lgTalk )
	{
		/* this code prints pretty lines at top of output box */
		int indent = (int)((122 - strlen(t_version::Inst().chVersion))/2);
		fprintf( ioQQQ, "%*cCloudy %s\n", indent, ' ', t_version::Inst().chVersion );

		fprintf( ioQQQ, "%57cwww.nublado.org\n\n", ' ' );

		/* now print box and date of version, before printing commands */
		fprintf( ioQQQ, "%23c", ' ' );
		fprintf( ioQQQ, "**************************************");
		fprintf( ioQQQ, "%7.7s", t_version::Inst().chDate);
		fprintf( ioQQQ, "**************************************\n");

		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );
	}

	/* read in commands and print them   */

	/* following signals which read is in progress,
	 * 1 is forward read, of input commands
	 * -1 is reverse read from bottom of line array, of cloudy.ini file */
	input.iReadWay = 1;

	/* initialize array reader, this sub does nothing but set
	 * initial value of a variable, depending on value of iReadWay */
	input.init();

	static const CloudyCommand commands[] = {
		{"ABSO",ParseAbsMag},
		/* enter luminosity in absolute magnitudes, in reads2 */			
		{"AGE",ParseAge},
		/* enter age of cloud so we can check for time-steady reads2 */
		{"AGN ",ParseAgn},
		/* enter generic style AGN continuum, in reads2 */
		{"ABUN",ParseAbundancesNonSolar},
		{"APER",ParseAperture},
		{"ATOM",ParseAtom},
		{"BACK", ParseBackgrd},
		/* cosmic background, in parse_backgrd */
		{"BLAC", ParseBlackbody},
		/* black body, in reads2 */
		{"BREM", ParseBremsstrahlung},
		{"CASE",ParseCaseB},
		/* do Case A or Case B (usually Case B since most realistic) */
		{"CEXT",ParseCExtra},
		/* CMB command */
		{"CLUM",ParseFill},
		{"CMB", ParseCMBOuter},
		/* cosmic thermal background radiation, argument is redshift */
		/* if no number on line then (realnum)FFmtRead returns z=0; i.e., now */
		{"COMP", ParseCompile},
		/* compile ascii version of stellar atmosphere continua in volk */
		{"CONS",ParseConstant},
		/* constant temperature, pressure, density, or gas pressure
		 * in readsun */
		{"CORO",ParseCoronal},
		/* coronal equilibrium; set constant temperature to number on line
		 *  in readsun */
		{"COSMOLOGY",ParseCosmology},
		{"COSMI", ParseCosmicRays},
		{"COSM",ParseCosm}, // Backstop for ambiguity
		{"COVE",ParseCovering},
		{"CRAS",ParseCrashDo},
		/* any of several tests to cause the code to crash */
		{"CYLI",ParseCylinder},
		{"DARK",ParseDarkMatter},
		{"DIEL",ParseDielectronic},
		{"DIFF",ParseDiffuse},
		{"DIST",ParseDistance},
		{"DLAW",ParseDLaw},
		/* either use dense_fabden routine, or read in table of points */
		{"DOUB",ParseDoubleTau},
		/* double optical depth scale after each iteration */
		{"DRIV",ParseDrive},
		/* call one of several drivers, readsun */
		{"EDEN",ParseEden},
		/* option to add "extra" electrons, to test Compton limit
		 *  for very low T(star) - option is log of eden */
		{"ELEM",ParseElement},
		/* element command;
		 * scale or abundance options, to change abundance of specific element
		 * read option to change order of elements
		 * in reads2.f */
		{"ENER",ParseEnergy},			
		{"EXTI",ParseExtinguish},
		/* extinguish ionizing continuum by absorbing column AFTER
		 * setting luminosity or Q(H).  First number is the column
		 * density (log), second number is leakage (def=0%)
		 * last number is lowest energy (ryd), last two may be omitted
		 * from right to left 
		 * 
		 * extinction is actually done in extin, which is called by ContSetIntensity */
		{"FAIL",ParseFail},
		/* reset number of temp failures allowed, default=20 */
		{"FILL",ParseFill},
		{"FLUC",ParseFluc},
		/* rapid density fluctuations, in readsun */
		{"F(NU",ParseF_nuSpecific},
		/* this is the specific flux at nu
		 *  following says F(nu) not nuF(nu) */
		{"FORC",ParseForceTemperature},
		/* force temperature of first zone, but don't keep constant
		 * allow to then go to nearest equilibrium
		 *  log if < 10 */
		{"FUDG",ParseFudge},
		{"GLOB",ParseGlobule},
		/* globule with density increasing inward
		 * parameters are outer density, radius of globule, and density power */
		{"GRAI",ParseGrain},
		/* read parameters dealing with grains, in reads2 */
		{"PGRA",ParsePGrains},
		{"GRAV",ParseGravity},
		/* (self-)Gravity forces: Yago Ascasibar (UAM, Spring 2009) */
		{"GRID",ParseGrid},
		/* option to run grid of models by varying certain parameters
		 * in reads2 */
		{"HDEN",ParseHDEN},
		/* parse the hden command to set the hydrogen density, in reads2 */
		{"HELI",ParseHeLike},
		{"HELP",ParseHelp},
		{"HEXT",ParseHExtra},
		/* "extra" heating rate, so that first= log(erg(cm-3, s-1},
		 * second optional number is scale radius, so that HXTOT = TurbHeat*SEXP(DEPTH/SCALE)
		 * if missing then constant heating.
		 * third option is depth from shielded face, to mimic irradiation from both sides*/
		{"HIGH",ParseConvHighT},
		/* approach equilibrium from high te */
		{"HYDROGEN",ParseHydrogen},
		{"ILLU",ParseIlluminate},
		// illuminate command
		{"INIT",ParseInitCount},
		{"INTEN",ParseIntensity},
		{"INTER",ParseInterp},
		/* interpolate on input tables for continuum, set of power  laws used
		 * input ordered pairs nu( ryd or log(Hz},, log(fnu)
		 * additional lines begin CONTINUE
		 * first check that this is the one and only INTERP command
		 * in readsun */
		{"IONI",ParseIonParI},
		/* inter ionization parameter U=Q/12 R*R N C;
		 * defined per hydrogen, not per electron (as before)
		 * radius must also be entered if spherical, not needed if plane */
		{"ITER",ParseIterations},
		{"L(NU",ParseL_nu},
		{"LASE",ParseLaser},
		{"LUMI",ParseLuminosity},
		{"MAGN", ParseMagnet},
		/* parse the magnetic field command, routine in magnetic.c */
		{"MAP ",ParseMap},
		/* do cooling space map for specified zones, in reads2 */
		{"META",ParseMetal},
		/* read depletion for metals, all elements heavier than He
		 * in reads2 */
		{"MONI",ParseMonitorResults},
		/* monitor that code predicts certain results, in monitor_results.h */
		{"NEUT",ParseNeutrons},
		{"NO ", ParseDont},
		/* don't do something, in readsun */
		{"NORM",ParseNorm},
		/* normalize lines to this rather than h-b, sec number is scale factor */
		{"NUF(",ParseNuF_nu},
		{"NUL(",ParseNuL_nu},
		{"OPTI",ParseOptimize},
		/* option to optimize the model by varying certain parameters
		 * in reads2 */
		{"PHI(",ParsePhi},
		{"PLOT", ParsePlot},
		/* make plot of log(nu.f(nu}, vs log(nu) or opacity
		 * in file plot */
		{"POWE", ParsePowerlawContinuum},
		/* power law with cutoff, in reads2 */
		{"PRIN",ParsePrint},
		/* adjust print schedule, in readsun */
		{"PUNC",ParseSave},
		/* save something, in save */
		{"Q(H)",ParseQH},
		{"RATI",ParseRatio},
		/* enter a continuum luminosity as a ratio of
		 * nuFnu for this continuum relative to a previous continuum
		 * format; first number is ratio of second to first continuum
		 * second number is energy for this ratio
		 * if third number on line, then 2nd number is energy of
		 * first continuum, while 3rd number is energy of second continuum
		 * in reads2 */
		{"RADI", ParseRadius},
		/* log of inner and outer radii, default second=infinity,
		 * if R2<R1 then R2=R1+R2
		 * there is an optional keyword, "PARSEC" on the line
		 * to use PC as units, reads2 */
		{"ROBE",ParseRoberto},
		{"SAVE",ParseSave},
		{"SET ",ParseSet},
		/* set something, in reads2 */
		{"SPEC", ParseSpecial},
		{"SPHE", ParseSphere},
		/* compute a spherical model, diffuse field from other side in
		 * in reads2 */
		{"STAT",ParseState},
		/* either get or put the code's state as a file on disk */
		{"STOP",ParseStop},
		/* stop model at desired zone, temperature, column density or tau-912
		 * in readsun */
		{"TABL",ParseTable},
		/* interpolate on input tables for continuum, set of power  laws used
		 * input stored in big BLOCK data
		 * first check that this is the one and only INTERP command
		 * in readsun */
		{"TAUM",ParseTauMin},
		{"TEST",ParseTest},
		/* parse the test command and its options */
		{"TIME",ParseDynaTime},
		/* parse the time dependent command, in dynamics.c */
		{"TITL",ParseTitle},
		{"TLAW",	ParseTLaw},
		/* some temperature vs depth law */
		{"TOLE", ParseTolerance},
		{"TRAC", ParseTrace},
		/* turn on trace output, in reads2 */
		{"VLAW",ParseVLaw},
		{"TURB",ParseTurbulence},
		{"WIND",ParseDynaWind},
		/* NB - advection and wind commands are now a single command */
		/* parse the wind command, in dynamics.c */
		{"XI",ParseIonParX},
		{NULL,NULL}}; // {NULL,NULL} sentinel must come last

	Parser p(commands);

	p.m_nqh = 0;
	p.m_nInitFile=0;/* used to count number of init files read in */
	p.m_lgDSet = false;
	p.m_lgEOF = false;

	/* read until eof or blank line, then return control back to main program */

	while (p.getline())
	{
		/* line starting with blank is taken as end of commands */
		if ( p.last( ) )
			break;

		/* echo the line but only if it does not contain the keyword HIDE */
		p.echo( );

		/* check whether VARY is on line */
		if( p.nMatch("VARY") )
		{
			optimize.lgVarOn = true;
			if( optimize.nparm + 1 > LIMPAR )
			{
				fprintf( ioQQQ, " Too many VARY lines entered; the limit is%4ld\n", 
				  LIMPAR );
				cdEXIT(EXIT_FAILURE);
			}
		}
		else
		{
			optimize.lgVarOn = false;
		}

		if( p.isCommandComment() )
		{
			((void)0); // do nothing for comments
		}
		else if( p.isVar() )
		{
			p.doSetVar();
		} 
		else
		{
			long int i;
			for (i=0; commands[i].name != NULL; ++i)
			{
				if (p.Command(commands[i].name,commands[i].action))
					break;
			}
			if (commands[i].name == NULL)
			{
				p.CommandError();	// No command was recognized
			}
		}
	}

	/*------------------------------------------------------------------- */
	/* fall through - hit lgEOF or blank line */

	if( called.lgTalk )
	{
		fprintf( ioQQQ, "%23c*%81c*\n", ' ', ' ' );
		fprintf( ioQQQ, "%23c***********************************************************************************\n\n\n\n", ' ' );
	}

	/* this is an option to turn on trace printout on the nth
	 * call from the optimizer - only allow trace if
	 * this is the case and nOptimiz 1 below nTrOpt */
	if( optimize.lgTrOpt )
	{
		/* nTrOpt was set with "optimize trace" command,
		 * is iteration to turn on trace */
		if( optimize.nTrOpt != optimize.nOptimiz )
		{
			trace.lgTrace = false;
			/* following overrides turning on trace elsewhere */
			trace.lgTrOvrd = false;
		}
		else
		{
			trace.lgTrace = true;
			called.lgTalk = cpu.i().lgMPI_talk();
			trace.lgTrOvrd = true;
			fprintf( ioQQQ, " READR turns on trace from optimize option.\n" );
		}
	}

	/* set density from DLAW command, must be done here since it may depend on later input commands */
	if( strcmp(dense.chDenseLaw,"DLW1") == 0 )
	{
		dense.SetGasPhaseDensity( ipHYDROGEN, (realnum)dense_fabden(radius.Radius,radius.depth) );

		if( dense.gas_phase[ipHYDROGEN] <= 0. )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density set by DLAW must be > 0.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if( strcmp(dense.chDenseLaw,"DLW2") == 0 )
	{
		dense.SetGasPhaseDensity( ipHYDROGEN, (realnum)interpol_tabulated(
						radius.Radius, radius.depth,
						dense.lgTabDepth, dense.lgTabLinear,
						dense.tabrad, dense.tabval,
						dense.nvals ) );

		if( dense.gas_phase[ipHYDROGEN] <= 0. )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density set by DLAW must be > 0.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else if( strcmp(dense.chDenseLaw,"DLW3") == 0 )
	{
		dense.SetGasPhaseDensity( ipHYDROGEN, (realnum)dense_parametric_wind(radius.Radius) );

		if( dense.gas_phase[ipHYDROGEN] <= 0. )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density set by DLAW must be > 0.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	/* start checks on parameters set properly - this begins with same line saying start .. */

	/* lgStop_not_enough_info says that not enough info for model, so stop 
	 * set true in following tests if anything missing */
	lgStop_not_enough_info = false;
	lgStop = false;

	/* check whether hydrogen density has been set - this value was set to 0 in zero */
	if( dense.gas_phase[ipHYDROGEN] <= 0. )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Hydrogen density MUST be specified.\n" );
		lgStop_not_enough_info = true;
		lgStop = true;
		/* need to set something since used below - will abort
		 * since lgStop is set */
		dense.SetGasPhaseDensity( ipHYDROGEN, 1. );
	}

	/* the SAVE XSPEC command cannot be combined with negative increments on the GRID command */
	if( grid.lgSaveXspec && grid.lgNegativeIncrements )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " PROBLEM DISASTER The SAVE XSPEC command cannot be combined with negative grid increments.\n" );
			fprintf( ioQQQ, " PROBLEM DISASTER Please check your GRID commands.\n\n\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}

	if( geometry.covaper < 0.f || geometry.iEmissPower == 2 )
		geometry.covaper = geometry.covgeo;

	radius.rinner = radius.Radius;

	/* mass flux for wind model - used for mass conservation */
	wind.emdot = dense.gas_phase[ipHYDROGEN]*wind.windv0;

	/* set converge criteria - limit number of iterations and zones */
	if( iterations.lgConverge_set)
	{
		iterations.itermx = MIN2( iterations.itermx , iterations.lim_iter );
		for( long int j=0; j < iterations.iter_malloc; j++ )
		{
			geometry.nend[j] = MIN2( geometry.nend[j] , iterations.lim_zone );
		}
	}

	/* NSAVE is number of lines saved, =0 if no commands entered */
	if( input.nSave < 0 )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER No commands were entered - whats up?\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* iterate to convergence and wind models are mutually exclusive */
	if( wind.lgBallistic() && conv.lgAutoIt )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE PROBLEM Due to the nature of the Sobolev approximation, it makes no sense to converge a windy model.\n" );
			fprintf( ioQQQ, " NOTE Iterate to convergence is turned off\n\n\n" );
		}
		conv.lgAutoIt = false;
		iterations.itermx = 0;
	}

	/* iterate to convergence and case b do not make sense together */
	/* WJH 22 May 2004: unless we are doing i-front dynamics (negative wind.windv0) */
	if( opac.lgCaseB && conv.lgAutoIt && (wind.lgBallistic() || wind.lgStatic()) )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE Case B is an artificial test, it makes no sense to converge this model.\n" );
			fprintf( ioQQQ, " NOTE Iterate to convergence is turned off.\n\n\n" );
		}
		conv.lgAutoIt = false;
		iterations.itermx = 0;
	}

	/* specifying a density power and constant pressure makes no sense */
	if( dense.DensityPower!=0. && strcmp( dense.chDenseLaw, "CPRE" )==0 )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE Specifying both a density power law and constant pressure is impossible.\n" );
		}
		lgStop = true;
	}

	if( !rfield.lgIonizReevaluate && strcmp( dense.chDenseLaw, "CDEN" )!=0 )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE NO REEVALUATE IONIZATION can only be used with constant density.\n" );
			fprintf( ioQQQ, " NOTE Resetting to reevaluate ionization.\n\n" );
		}
		rfield.lgIonizReevaluate = true;
	}

	if( !rfield.lgOpacityReevaluate && strcmp( dense.chDenseLaw, "CDEN" )!=0 )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE NO REEVALUATE OPACITY can only be used with constant density.\n" );
			fprintf( ioQQQ, " NOTE Resetting to reevaluate opacity.\n\n" );
		}
		rfield.lgOpacityReevaluate = true;
	}

	/* check that a symmetry is specified if gravity from an external mass has been added */
	if( pressure.external_mass[0].size()>0 && pressure.gravity_symmetry==-1 )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE Gravity from an external mass has been added, but no symmetry (spherical/mid-plane) was specified.\n" );
			fprintf( ioQQQ, " NOTE It will be ignored.\n\n\n" );
		}
	}

	/* check if the combination of stopping column density and H density are physically plausible */
	double thickness = min( MIN3( StopCalc.tauend, StopCalc.colpls, StopCalc.colnut ),
			 MIN3( StopCalc.col_h2, StopCalc.col_h2_nut, StopCalc.HColStop ) );
	if( thickness < COLUMN_INIT )
	{
		/* a stop column density was specified - check on physical thickness this corresponds to */
		thickness /= (dense.gas_phase[ipHYDROGEN]*geometry.FillFac);
		/* don't complain if outer radius set small with `stop thickness' command. */
		if( thickness > 1e25 && radius.StopThickness[0] > 1e25 )
		{
			fprintf( ioQQQ, 
				"NOTE The specified column density and hydrogen density correspond to a thickness of %.2e cm.\n",
				thickness);
			fprintf( ioQQQ, 
				"NOTE This seems large to me.\n");
			fprintf(ioQQQ,"NOTE a very large radius may cause overflow.\n\n");
		}
	}

	if( gv.lgDColOn && thermal.ConstGrainTemp>0 && called.lgTalk )
	{
		/* warn if constant grain temperature but gas-grain thermal effects
		 * are still included */
		fprintf( ioQQQ, 
			"NOTE The grain temperatures are set to a constant value with the "
			"CONSTANT GRAIN TEMPERATURE command, but "
			"energy exchange \n");
		fprintf( ioQQQ, 
			"NOTE is still included.  The grain-gas heating-cooling will be incorrect.  "
			"Consider turning off gas-grain collisional energy\n");
		fprintf( ioQQQ, 
			"NOTE exchange with the NO GRAIN GAS COLLISIONAL ENERGY EXCHANGE command.\n\n\n");
	}

	if( !rfield.lgDoLineTrans && rfield.lgOpacityFine )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE NO LINE TRANSER set but fine opacities still computed.\n" );
			fprintf( ioQQQ, " NOTE Turning off fine opacities.\n\n" );
		}
		rfield.lgOpacityFine = false;
	}

	if( h2.lgEnabled && (!rfield.lgDoLineTrans || !rfield.lgOpacityFine) )
	{
		if( called.lgTalk )
		{
			fprintf( ioQQQ, " NOTE Large H2 molecule turned on but line transfer and fine opacities are not.\n" );
			fprintf( ioQQQ, " NOTE Turning on line transfer and fine opacities.\n\n" );
		}
		rfield.lgOpacityFine = true;
		rfield.lgDoLineTrans = true;
	}

	if( rfield.lgMustBlockHIon && !rfield.lgBlockHIon )
	{
		/* one of the input continua needs to have H-ionizing radiation
		 * blocked with extinguish command, but it was not done */
		fprintf( ioQQQ, "\n NOTE\n"
			" NOTE One of the incident continuum is a form used when no H-ionizing radiation is produced.\n" );
		fprintf( ioQQQ, " NOTE You must also include the EXTINGUISH command to make sure this is done.\n" );
		fprintf( ioQQQ, " NOTE The EXTINGUISH command was not included.\n" );
		fprintf( ioQQQ, " NOTE YOU MAY BE MAKING A BIG MISTAKE!!\n NOTE\n\n\n\n" );
	}

	/* if stop temp set below default then we are going into cold and possibly molecular
	 * gas - check some parameters in this case */
	if( called.lgTalk && (StopCalc.TempLoStopZone < phycon.TEMP_STOP_DEFAULT || 
		/* thermal.ConstTemp def is zero, set pos when constant temperature is set */
		(thermal.ConstTemp > 0. && thermal.ConstTemp < phycon.TEMP_STOP_DEFAULT ) ) )
	{
		/* print warning if temperature set below default but cosmic rays not turned on 
		 * do not print if molecules are off */
		if( (hextra.cryden == 0.) && !mole_global.lgNoMole )
		{
			fprintf( ioQQQ, "\n NOTE\n"
							" NOTE The simulation is going into neutral gas but cosmic rays are not included.\n" );
			fprintf( ioQQQ, " NOTE Ion-molecule chemistry will not occur without a source of ionization.\n" );
			fprintf( ioQQQ, " NOTE The chemistry network may collapse deep in molecular regions.\n" );
			fprintf( ioQQQ, " NOTE Consider adding galactic background cosmic rays with the COSMIC RAYS BACKGROUND command.\n" );
			fprintf( ioQQQ, " NOTE You may be making a BIG mistake.\n NOTE\n\n\n\n" );
		}
	}

	/* dense.gas_phase[ipHYDROGEN] is linear hydrogen density (cm-3) */
	/* test for hydrogen density properly set has already been done above */
	if( called.lgTalk && dense.gas_phase[ipHYDROGEN] < 1e-4 )
	{
		fprintf( ioQQQ, " NOTE Is the entered value of the hydrogen density (%.2e) reasonable?\n",
			dense.gas_phase[ipHYDROGEN]);
		fprintf( ioQQQ, " NOTE It seems pretty low to me.\n\n\n" );
	}
	else if( called.lgTalk && dense.gas_phase[ipHYDROGEN] > 1e15 )
	{
		fprintf( ioQQQ, " NOTE Is this value of the hydrogen density reasonable?\n" );
		fprintf( ioQQQ, " NOTE It seems pretty high to me.\n\n\n" );
	}

	/* is the model going to crash because of extreme density? */
	if( called.lgTalk && !lgStop && !lgStop_not_enough_info )
	{
		if( dense.gas_phase[ipHYDROGEN] < 1e-6 || dense.gas_phase[ipHYDROGEN] > 1e19 )
		{
			fprintf( ioQQQ, " NOTE Simulation may crash because of extreme "
				"density.  The value was %.2e\n\n" , 
				dense.gas_phase[ipHYDROGEN] );
		}
	}

	if( rfield.nShape == 0 )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER No incident radiation field was specified - "
			"at least put in the CMB.\n" );
		lgStop = true;
		lgStop_not_enough_info = true;
	}

	if( (p.m_nqh) == 0 )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Luminosity of continuum MUST be specified.\n" );
		lgStop = true;
		lgStop_not_enough_info = true;
	}

	/* Testing on 0 is safe since this it is initialized that way
	 * only print comment if continuum has been specified,
	 * if continuum not given then we are aborting anyway */
	if( radius.Radius == 0. && rfield.nShape> 0)
	{
		fprintf( ioQQQ, " PROBLEM DISASTER Starting radius MUST be specified.\n" );
		lgStop = true;
		lgStop_not_enough_info = true;
	}

	if( rfield.nShape != (p.m_nqh) )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER There were not the same number of continuum shapes and luminosities entered.\n" );
		lgStop = true;
	}

	/* we only want to do this test on the first call to the command
	 * parser - it will be called many more times but with no grid command
	 * during the actual grid calculation */
	static bool lgFirstPass = true;

	/* the NO VARY option sets this flag, and can be used to turn off
	 * the grid command, as well as the optimizer */
	if( optimize.lgNoVary && grid.lgGrid )
	{
		/* ignore grids */
		grid.lgGrid = false;
		optimize.nparm = 0;
		grid.nGridCommands = 0;
	}

	if( lgFirstPass && grid.lgGrid && (optimize.nparm!=grid.nGridCommands) )
	{
		/* number of grid vary options do match */
		fprintf( ioQQQ, " PROBLEM DISASTER The GRID command was entered "
			"but there were %li GRID commands and %li commands with a VARY option.\n" ,
			grid.nGridCommands , optimize.nparm);
		fprintf( ioQQQ, " There must be the same number of GRIDs and VARY.\n" );
		lgStop = true;
	}
	lgFirstPass = false;

	if( lgStop_not_enough_info )
	{
		fprintf( ioQQQ, " PROBLEM DISASTER I do not have enough information to do the simulation, I cannot go on.\n" );
		fprintf( ioQQQ, "\n\n Sorry.\n\n\n" );
		cdEXIT(EXIT_FAILURE);
	}

	{
		bool lgParserTest = false;
		if ( lgParserTest ) 
		{
			// Quit after parse phase, to allow quick verification that
			// there are no immediate syntax handling failures.
			fprintf(ioQQQ,"Parser phase PASSED\n");
			exit(0);
		}
	}

	if( lgStop )
	{
		cdEXIT(EXIT_FAILURE);
	}

	/* end checks on parameters set properly - this begins with same line saying start .. */
	return;
}

void ParseAbundancesNonSolar(Parser &p)
{
	/* chemical abundances, readsun, p.m_lgDSet is set true if grains command
	 * ever entered, so that abundances command does not override grains command */
	ParseAbundances(p);
	/* abundances no longer solar */
	abund.lgAbnSolar = false;
}

void ParseAperture(Parser &p)
{
	DEBUG_ENTRY("ParseAperture()");
	/* aperture command to simulate pencil beam or long slit */
	
	/* if the "BEAM" or "SLIT" option is specified then only part 
	 * of the geometry is observed, and intensities
	 * should not be weighted by r^2.  There are two limiting cases, SLIT in which
	 * the slit is longer than the diameter of the nebula and the contribution to the
	 * detected luminosity goes as r^1, and BEAM when the contribution is r^0, 
	 * the same as plane parallel 
	 * 
	 * default value of geometry.iEmissPower is 2 (set in zero.c) for full geometry  
	 */
	if( p.nMatch("SLIT") )
	{
		/* long slit is case where slit is longer than diameter, so emissivity contributes
		 * r^1 to the observed luminosity */
		geometry.iEmissPower = 1;
	}
	else if( p.nMatch("BEAM") )
	{
		/* pencil beam is case where we view the nebula through a narrow square
		 * centered on the central source, this gives r^0 dependence */
		geometry.iEmissPower = 0;
	}
	else if( p.nMatch("SIZE") )
	{
		/* set the aperture size: slit width or suface area of the pencil beam */
		/* units are arcsec for slit width, or arcsec^2 for pencil beam */
		geometry.size = realnum(p.FFmtRead());
		if( p.lgEOL() )
		{
			p.NoNumb("aperture size");
		}
		if( geometry.size <= 0.f )
		{
			fprintf( ioQQQ, " The aperture size must be positive.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		geometry.lgSizeSet = true;
	}
	else if( p.nMatch("COVE") )
	{
		/* set the aperture covering factor, see Hazy 1 for a discussion
		 * this is a dimensionless fraction between 0 and 1 */
		geometry.covaper = realnum(p.FFmtRead());
		if( p.lgEOL() )
		{
			p.NoNumb("aperture covering factor");
		}
		if( geometry.covaper <= 0.f || geometry.covaper > 1.f )
		{
			fprintf( ioQQQ, " The aperture covering factor must be > 0 and <= 1.  Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		fprintf( ioQQQ, " One of the keywords SLIT, BEAM, SIZE or COVEring factor must appear.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
void ParseAtom(Parser &p)
{
	DEBUG_ENTRY("ParseAtom()");
	/* accept both forms of feii */
	if( p.nMatch("FEII") || p.nMatch("FE II") )
	{
		/* parse the atom FeII command - routine in atom_feii.c */
		ParseAtomFeII(p);
	}
	
	else if( p.nMatch("H-LI") )
	{
		/* parse the atom h-like command */
		ParseAtomISO(ipH_LIKE, p);
	}
	
	else if( p.nMatch("HE-L") )
	{
		/* parse the atom he-like command */
		ParseAtomISO(ipHE_LIKE, p);
	}
	
	else if( p.nMatch("ROTO") || p.nMatch(" CO ") )
	{
		/* ATOM ROTOR same as ATOM CO */
		fprintf(ioQQQ," The old CO models no longer exist, and this command is no longer supported.\n" );
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	else if( p.nMatch(" H2 ") || p.nMatch(" HD ") )
	{
		ParseAtomH2(p);
	}

	/* enable models from CHIANTI database */
	else if (p.nMatch("CHIA"))
	{
		char chString_quotes_lowercase[INPUT_LINE_LENGTH];
		bool lgQuotesFound = true;
		if (p.GetQuote(chString_quotes_lowercase, false))
			lgQuotesFound = false;

		// option to specify different CloudyChianti.ini file, was initialized
		// with string CloudyChianti.ini
		if (lgQuotesFound == true)
			strcpy(atmdat.chCloudyChiantiFile, chString_quotes_lowercase);

		atmdat.lgChiantiOn = true;
		if (p.nMatch(" OFF"))
		{
			atmdat.lgChiantiOn = false;
			atmdat.lgChiantiHybrid = false;
		}

		// hybrid, chianti with OP for higher energy lines
		// Turn off hybrid, use Chianti only
		if (p.nMatch(" NO HYBR"))
			atmdat.lgChiantiHybrid = false;

		// Print which species are being used in output and # of levels
		if (p.nMatch(" PR"))
			atmdat.lgChiantiPrint = true;

		// Use Experimental energies exclusively. Default use experimental.
		if (p.nMatch(" THEO"))
			atmdat.lgChiantiExp = false;

		// Input the maximum number of Chianti levels to use
		if (p.nMatch(" LEV"))
		{
			if (p.nMatch(" MAX"))
			{
				atmdat.nChiantiMaxLevels = LONG_MAX;
				atmdat.nChiantiMaxLevelsFe = LONG_MAX;
			}
			else
			{
				atmdat.nChiantiMaxLevelsFe = (long)p.FFmtRead();
				atmdat.nChiantiMaxLevels = (long)p.FFmtRead();
				if (p.lgEOL())
				{
					p.NoNumb("two numbers, the maximum number of levels in Fe, and in other elements,");
				}
				if( atmdat.nChiantiMaxLevelsFe < 2 || atmdat.nChiantiMaxLevels < 2 )
				{
					fprintf(ioQQQ,
							" \nPROBLEM The maximum number of chianti levels should be two or greater.\n");
					fprintf(ioQQQ, " To turn off the Chianti data use \"Set Chianti off\" instead.\n");
					fprintf(ioQQQ, " See Hazy 1 for details.\n");
					cdEXIT( EXIT_FAILURE );
				}
			}
			atmdat.lgChiantiLevelsSet = true;
		}
	}

	/* enable models from STOUT database */
	else if (p.nMatch("STOUT"))
	{
		char chString_quotes_lowercase[INPUT_LINE_LENGTH];
		bool lgQuotesFound = true;
		if (p.GetQuote(chString_quotes_lowercase, false))
			lgQuotesFound = false;

		// option to specify different Stout.ini file, was initialized
		// with string Stout.ini
		if (lgQuotesFound == true)
			strcpy(atmdat.chStoutFile, chString_quotes_lowercase);

		atmdat.lgStoutOn = true;
		// Print which species are being used in output and # of levels
		if (p.nMatch(" PR"))
			atmdat.lgStoutPrint = true;

		if (p.nMatch(" OFF"))
		{
			atmdat.lgStoutOn = false;
			atmdat.lgStoutHybrid = false;
		}

		// hybrid, Stout with OP for higher energy lines
		// Turn off hybrid, use Stout only
		if (p.nMatch(" NO HYBR"))
			atmdat.lgStoutHybrid = false;

		// Input the maximum number of Stout levels to use
		if (p.nMatch(" LEV"))
		{
			if (p.nMatch(" MAX"))
			{
				atmdat.nStoutMaxLevels = LONG_MAX;
			}
			else
			{
				atmdat.nStoutMaxLevels = (long)p.FFmtRead();
				if (p.lgEOL())
				{
					p.NoNumb("the maximum number of Stout levels,");
				}
				if( atmdat.nStoutMaxLevels < 2 )
				{
					fprintf(ioQQQ,
							" \nPROBLEM The maximum number of Stout levels should be two or greater.\n");
					fprintf(ioQQQ, " To turn off the Stout data use \"Set Stout off\" instead.\n");
					fprintf(ioQQQ, " See Hazy 1 for details.\n");
					cdEXIT( EXIT_FAILURE );
				}
			}
		}

	}

	/* enable models from LAMDA (Leiden Atomic and Molecular Database) */
	else if (p.nMatch("LAMD"))
	{
		if (p.nMatch(" OFF"))
			atmdat.lgLamdaOn = false;
		else if (p.nMatch(" ON "))
			atmdat.lgLamdaOn = true;
		else
		{
			/* should not have happened ... */
			fprintf(ioQQQ, " There should have been an option on this SET LAMDA command.\n");
			fprintf(ioQQQ, " consult Hazy to find valid options.\n Sorry.\n");
			cdEXIT(EXIT_FAILURE);
		}
	}
	else
	{
		fprintf( ioQQQ, " I could not recognize a keyword on this atom command.\n");
		fprintf( ioQQQ, " The available keys are FeII, H-Like, He-like, rotor and H2.\n");
		fprintf( ioQQQ, " Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
void ParseBremsstrahlung(Parser &p)
{
	DEBUG_ENTRY("ParseBremsstrahlung()");
	/* bremsstrahlung continuum from central object */
	strcpy( rfield.chSpType[rfield.nShape], "BREMS" );
	rfield.slope[rfield.nShape] = 
		(realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("temperature");
	}
	
	/* temperature is interpreted as log if <= 10, or if keyword is present */
	if( rfield.slope[rfield.nShape] <= 10. || p.nMatch(" LOG") )
	{
		rfield.slope[rfield.nShape] =  
			pow(10.,rfield.slope[rfield.nShape]);
	}
	rfield.cutoff[rfield.nShape][0] = 0.;
	
	/* option for vary keyword */
	if( optimize.lgVarOn )
	{
		/* only one parameter */
		optimize.nvarxt[optimize.nparm] = 1;
		strcpy( optimize.chVarFmt[optimize.nparm], "BREMS, T=%f LOG" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* log of temp will be pointer */
		optimize.vparm[0][optimize.nparm] =  (realnum)log10(rfield.slope[rfield.nShape]);
		optimize.vincr[optimize.nparm] = 0.5;
		++optimize.nparm;
	}
	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
void ParseCExtra(Parser &p)
{
	/* add "extra" cooling, power law temp dependence */
	thermal.lgCExtraOn = true;
	thermal.CoolExtra = (realnum)pow(10.,p.FFmtRead());
	if( p.lgEOL() )
	{
		p.NoNumb("extra cooling");
	}
	thermal.cextpw = (realnum)p.FFmtRead();
}
void ParseCMBOuter(Parser &p)
{
	cosmology.redshift_start = (realnum)p.FFmtRead();
	cosmology.redshift_current = cosmology.redshift_start;
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	ParseCMB(cosmology.redshift_start,&(p.m_nqh));
	
	/* option to also set the hydrogen density at given redshift. */
	if( p.nMatch("DENS") && !p.lgEOL() )
	{
		char chStuff[INPUT_LINE_LENGTH];
		
		/* hydrogen density */
		sprintf( chStuff , "HDEN %.2e LINEAR", GetDensity( cosmology.redshift_start ) );
		p.setline(chStuff);
		p.set_point(4);
		ParseHDEN(p);
	}
}
void ParseCosm(Parser &)
{
	DEBUG_ENTRY("ParseCosm()");
	fprintf(ioQQQ,
			  "This command is now ambiguous -- please specify either COSMIC RAYS or COSMOLOGY.\nSorry.\n");
	cdEXIT(EXIT_FAILURE);
}
void ParseCovering(Parser &p)
{
	DEBUG_ENTRY("ParseCovering()");
	/* covering factor for gas */
	/* The geometric covering factor accounts for how much of 4\pi is
	 * covered by gas, and so linearly multiplies the predicted intensities */
	geometry.covgeo = (realnum)p.FFmtRead();
	
	if( p.lgEOL() )
	{
		p.NoNumb("covering factor");
	}
	
	/* if negative then log, convert to linear */
	if( geometry.covgeo <= 0. )
	{
		geometry.covgeo = (realnum)pow((realnum)10.f,geometry.covgeo);
	}
	
	if( geometry.covgeo > 1. )
	{
		fprintf( ioQQQ, " A covering factor greater than 1 makes no physical sense.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	/* radiative transfer covering factor will be equal to the geometric one */
	geometry.covrt = geometry.covgeo;
}
void ParseCylinder(Parser &p)
{
	/* height of cylinder in cm */
	radius.lgCylnOn = true;
	radius.CylindHigh = pow(10.,p.FFmtRead());
	if( p.lgEOL() )
	{
		p.NoNumb("height");
	}
}
void ParseDarkMatter(Parser &p)
{
	DEBUG_ENTRY( "ParseDarkMatter()" );

	if( p.nMatch(" NFW") )
	{
		/* specify an NFW profile */
		/* >>refer dark matter	Navarro, Frenk, & White, 1996, ApJ, 462, 563 */

		dark.r_200 = (realnum)p.getNumberCheckAlwaysLog("NFW r_200");
		/* Let r_s have a default corresponding to c_200 = 10. */
		dark.r_s = (realnum)p.getNumberDefaultAlwaysLog("NFW r_s", log10(dark.r_200)-1. );
		dark.lgNFW_Set = true;

		/* option for vary keyword */
		if( optimize.lgVarOn )
		{
			/* only one parameter */
			optimize.nvarxt[optimize.nparm] = 1;
			strcpy( optimize.chVarFmt[optimize.nparm], "DARK NFW %f" );
			/* pointer to where to write */
			optimize.nvfpnt[optimize.nparm] = input.nRead;
			/* log of temp will be pointer */
			optimize.vparm[0][optimize.nparm] =  (realnum)log10(dark.r_200);
			optimize.vincr[optimize.nparm] = 0.5;
			++optimize.nparm;
		}
	}
	else
	{
		fprintf( ioQQQ, " Did not recognize a valid option for this DARK command.\nSorry.\n\n" );
		cdEXIT(EXIT_FAILURE);
	}

	return;
}
void ParseDielectronic(Parser &)
{
	DEBUG_ENTRY("ParseDielectronic()");
	/* >>chng 05 dec 21, change dielectronic command to set dielectronic recombination command */
	fprintf( ioQQQ, " The DIELectronic command has been replaced with the SET DIELectronic recombination command.\n" );
	fprintf( ioQQQ, " Please have a look at Hazy.\n Sorry.\n\n" );
	cdEXIT(EXIT_FAILURE);
}
void ParseDiffuse(Parser &p)
{
	DEBUG_ENTRY("ParseDiffuse()");
	/* set method of transferring diffuse fields,
	 * default is outward only, cdDffTrns label "OU2", set in zero.c */
	if( p.nMatch(" OTS") )
	{
		if( p.nMatch("SIMP") )
		{
			/* this is simple ots, which never allows photons to escape */
			strcpy( rfield.chDffTrns, "OSS" );
		}
		else
		{
			/* this is regular ots, which allows photons to escape */
			strcpy( rfield.chDffTrns, "OTS" );
		}
		rfield.lgOutOnly = false;
	}
	else if( p.nMatch(" OUT") )
	{
		rfield.lgOutOnly = true;
		long int j = (long int)p.FFmtRead();
		if( p.lgEOL() )
		{
			/* this is the default set in zero */
			strcpy( rfield.chDffTrns, "OU2" );
		}
		else
		{
			if( j > 0 && j < 10 )
			{
				sprintf( rfield.chDffTrns, "OU%1ld", j );
			}
			else
			{
				fprintf( ioQQQ, " must be between 1 and 9 \n" );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	
	else
	{
		fprintf( ioQQQ, " There should have been OUTward or OTS on this line.  Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
}
void ParseDistance(Parser &p)
{
	/* distance to the object */
	radius.distance = p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("distance");
	}
			/* default is for quantity to be log of distance, linear keyword
			 * overrides this - is LINEAR is not on line then exp */
	if( !p.nMatch("LINE" ) )
	{
		radius.distance = pow(10., radius.distance );
	}
	/* default is radius in cm - if parsecs appears then must 
	 * convert to cm */
	if( p.nMatch("PARS") )
	{
		radius.distance *= PARSEC;
	}
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "DISTANCE %f LOG" );
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = realnum(log10(radius.distance));
		optimize.vincr[optimize.nparm] = 0.3f;
		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}
}
void ParseDoubleTau(Parser &)
{
	rt.DoubleTau = 2.;
}
void ParseEden(Parser &p)
{
	dense.EdenExtra = (realnum)pow(10.,p.FFmtRead());
	if( p.lgEOL() )
	{
		p.NoNumb("electron density");
	}
	/* warn that this model is meaningless */
	phycon.lgPhysOK = false;
}
void ParseEnergy(Parser &p)
{
	DEBUG_ENTRY("ParseEnergy()");
	/* energy density (degrees K) of this continuum source */
	if( (p.m_nqh) >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* silly, but calms down the lint */
	ASSERT( (p.m_nqh) < LIMSPC );
	strcpy( rfield.chRSpec[(p.m_nqh)], "SQCM" );
	realnum teset = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("energy density");
	}
	/* set radius to very large value if not already set */
	/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
	if( !radius.lgRadiusKnown )
	{
		/* set radius to large value */
		radius.Radius = pow(10.,radius.rdfalt);
	}
	
	/* convert temp to log, recognize linear option */
	if( !p.nMatch(" LOG") && (p.nMatch("LINE") || teset > 10.) )
	{
		/* option for linear temperature, must store log */
				teset = (realnum)log10(teset);
	}
	
	if( teset > 5. )
	{
		fprintf( ioQQQ, " This intensity may be too large.  The code may crash due to overflow.  Was log intended?\n" );
	}
	
	/* teset is now log of temp, now get log of total luminosity */
	strcpy( rfield.chSpNorm[(p.m_nqh)], "LUMI" );
	
	/* full range of continuum will be used */
	rfield.range[(p.m_nqh)][0] = rfield.emm;
	rfield.range[(p.m_nqh)][1] = rfield.egamry;
	rfield.totpow[(p.m_nqh)] = (4.*(teset) - 4.2464476 + 0.60206);
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "ENERGY DENSITY %f LOG" );
		if( rfield.lgTimeVary[(p.m_nqh)] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = teset;
		optimize.vincr[optimize.nparm] = 0.1f;
		optimize.nvarxt[optimize.nparm] = 1;
		++optimize.nparm;
	}
	
	/* finally increment number of continua */
	++(p.m_nqh);
}
void ParseFail(Parser &p)
{
	/* save previous value */
	long int j = conv.LimFail;
	conv.LimFail = (long int)p.FFmtRead();
	if( p.lgEOL() )
	{
				p.NoNumb("limit");
	}

	/* >>chng 01 mar 14, switch logic on maps */
	/* ' map' flag, make map when failure, default is no map,
	 * second check is so that no map does not trigger a map */
	if( p.nMatch(" MAP") && !p.nMatch(" NO ") )
	{
		conv.lgMap = true;
	}
	
	/* complain if failures was increased above default */
	if( conv.LimFail > j )
	{
		fprintf( ioQQQ, " This command should not be necessary.\n" );
		fprintf( ioQQQ, " Please show this input stream to Gary Ferland if this command is really needed for this simulation.\n" );
	}
}
void ParseFill(Parser &p)
{
	/* filling factor, power law exponent (default=1., 0.) */
	realnum a = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("filling factor");
	}
	
	if( a <= 0. || p.nMatch(" LOG") )
	{
		/* number less than or equal to 0, must have been entered as log */
		geometry.FillFac = (realnum)pow((realnum)10.f,a);
	}
	else
	{
		/* number greater than zero, must have been the real thing */
		geometry.FillFac = a;
	}
	if( geometry.FillFac > 1.0 )
	{
		if( called.lgTalk )
			fprintf( ioQQQ, " Filling factor > 1, reset to 1\n" );
		geometry.FillFac = 1.;
	}
	geometry.fiscal = geometry.FillFac;
	
	/* option to have power law dependence, filpow set to 0 in zerologic */
	geometry.filpow = (realnum)p.FFmtRead();
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "FILLING FACTOR= %f LOG power= %f" );
		
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)log10(geometry.FillFac);
		optimize.vincr[optimize.nparm] = 0.5f;
		
		/* power law dependence here, but cannot be varied */
		optimize.vparm[1][optimize.nparm] = geometry.filpow;
		optimize.nvarxt[optimize.nparm] = 2;
		
		/* do not allow filling factor to go positive */
		optimize.varang[optimize.nparm][0] = -FLT_MAX;
		optimize.varang[optimize.nparm][1] = 0.f;
		++optimize.nparm;
	}
}
void ParseF_nuSpecific(Parser &p)
{
	bool lgNu2 = false;
	/* in reads2 */
	ParseF_nu(p,"SQCM",lgNu2);
}
void ParseForceTemperature(Parser &p)
{
	thermal.ConstTemp = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("temperature");
	}
	
	if( p.nMatch(" LOG") || (thermal.ConstTemp <= 10. && 
									 !p.nMatch("LINE")) )
	{
		thermal.ConstTemp = (realnum)pow((realnum)10.f,thermal.ConstTemp);
	}
	
	/* low energy bounds of continuum array do not permit too-low a temp */
	if( thermal.ConstTemp < 3. )
	{
		fprintf( ioQQQ, " TE reset to 3K: entered number too small.\n" );
		thermal.ConstTemp = 3.;
	}
}	
void ParseFudge(Parser &p)
{
	/* enter a fudge factor */
	fudgec.nfudge = 0;
	for( long int j=0; j < NFUDGC; j++ )
	{
		fudgec.fudgea[j] = (realnum)p.FFmtRead();
		/* fudgec.nfudge is number of factors on the line */
		if( !p.lgEOL() )
			fudgec.nfudge = j+1;
	}
	if( fudgec.nfudge == 0 )
		p.NoNumb("fudge factor");
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		/* no luminosity options on vary */
		optimize.nvarxt[optimize.nparm] = 1;
		strcpy( optimize.chVarFmt[optimize.nparm], "FUDGE= %f" );
		/* enter remaining parameters as constants */
		char chHold[1000];
		for( long int j=1; j<fudgec.nfudge; ++j )
		{
			sprintf( chHold , " %f" , fudgec.fudgea[j] );
			strcat( optimize.chVarFmt[optimize.nparm] , chHold );
		}
		optimize.lgOptimizeAsLinear[optimize.nparm] = true;
		
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* fudge factor stored here  */
		optimize.vparm[0][optimize.nparm] = fudgec.fudgea[0];
		/* the increment in the first steps away from the original value */
		optimize.vincr[optimize.nparm] = abs(0.2f*fudgec.fudgea[0]);
		/* we have no clue what to use when initial estimate is 0... */
		if( optimize.vincr[optimize.nparm] == 0.f )
			optimize.vincr[optimize.nparm] = 1.f;
		++optimize.nparm;
	}
}
void ParsePGrains(Parser &)
{
	DEBUG_ENTRY("ParsePGrains()");
	fprintf(ioQQQ," Sorry, this command is obsolete, you can now use the normal GRAINS command.\n");
	cdEXIT(EXIT_FAILURE);
}
void ParseGravity(Parser &p)
{
	DEBUG_ENTRY("ParseGravity()");

	if( p.nMatch("EXTE") )
	{
		/* external mass: M/Msun or Sigma/(Msun/pc^2) depending on symmetry  */
		/*                if no number is read, FFmtRead returns 0           */
		/*                default is linear, unless "log" keyword is present */
		double M_i = p.FFmtRead();

		if( !p.lgEOL() && p.nMatch("LOG") ) 
			M_i = pow( 10., M_i );
		pressure.external_mass[0].push_back( M_i );

		/* extent of the mass distribution, in pc             */
		/* default is linear, unless "log" keyword is present */
		double x_i = p.FFmtRead();

		if( !p.lgEOL() && p.nMatch("LOG") ) 
			x_i = pow( 10., x_i );
		pressure.external_mass[1].push_back( x_i * PARSEC );

		/* exponential index of the mass distribution */
		pressure.external_mass[2].push_back( p.FFmtRead() );
	}
	else
	{
		/* choose spherical or mid-plane symmetry for the gas mass distribution */
		if( p.nMatch("SPHE") )
		{
			pressure.gravity_symmetry = 0;
		}
		else if( p.nMatch("PLAN") )
		{
			pressure.gravity_symmetry = 1;
		}
		else
		{
			fprintf( ioQQQ, " The symmetry of the gravitational mass must be specified explicitly. Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* external mass, proportional to the gas density      */
		/* default is linear, unless "log" keyword is present */
		pressure.self_mass_factor = p.FFmtRead();

		if( p.lgEOL() )
			pressure.self_mass_factor = 1.;
		else if( p.nMatch("LOG") )
		{
			pressure.self_mass_factor = pow( 10., pressure.self_mass_factor );
		}
	}
}
void ParseHeLike(Parser &)
{
	DEBUG_ENTRY("ParseHeLike()");
	fprintf(ioQQQ,"Sorry, this command is replaced with ATOM HE-LIKE\n");
	cdEXIT(EXIT_FAILURE);
}
void ParseHelp(Parser &p)
{
	DEBUG_ENTRY("ParseHelp()");
	p.help(ioQQQ);
}
void ParseHExtra(Parser &p)
{
	DEBUG_ENTRY(" ParseHExtra()");
	hextra.TurbHeat = (realnum)pow(10.,p.FFmtRead());
	if( p.lgEOL() )
		p.NoNumb("extra heating first parameter" );

	/* save initial value in case reset with time dependent command */
	hextra.TurbHeatSave = hextra.TurbHeat;
	
	/* remember which type of scale dependence we will use */
	const char *chHextraScale;
	/* these are optional numbers on depth or density option */
	realnum par1=0. , par2=0.;
	long int nPar;
	if( p.nMatch("DEPT") )
	{
		/* depend on depth */
		hextra.lgHextraDepth = true;
		chHextraScale = "DEPTH";
		/* optional scale radius */
		realnum a = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("depth");
		else
			hextra.turrad = (realnum)pow((realnum)10.f,a);
		
		/* depth from shielded face, to mimic illumination from both sides 
		 * may not be specified */
		realnum b = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			hextra.turback = 0.;
			nPar = 2;
		}
		else
		{
			hextra.turback = (realnum)pow((realnum)10.f,b);
			nPar = 3;
		}
		par1 = a;
		par2 = b;
	}
	else if( p.nMatch("DENS") )
	{
		/* depends on density */
		chHextraScale = "DENSITY";
		hextra.lgHextraDensity = true;
		
		/* the log of the density */
		realnum a = (realnum)p.FFmtRead();
		/* if number not entered then unity is used, heating 
		 * proportional to density */
		hextra.HextraScaleDensity = (realnum)pow((realnum)10.f,a);
		par1 = a;
		par2 = 0;
		nPar = 2;
	}
	else if( p.nMatch("SS") )
	{
		/* alpha disk model, specify alpha, mass of black hole, and distance */
		chHextraScale = "SS";
		hextra.lgHextraSS = true;
		
		/* the parameter alpha of alpha-disk model */
		hextra.HextraSSalpha = hextra.TurbHeat;
	
		/*  mass in solar masses of center black hole */
		realnum a = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("hextraSS Mass");
		hextra.HextraSS_M  = a*SOLAR_MASS;
		
		/* radius (cm) from center */
		realnum b = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("hextraSS radius");
		hextra.HextraSSradius = b;
		par1 = a;
		par2 = 0;
		nPar = 2;
	}
	else
	{
		/* constant heating */
		chHextraScale = "";
		nPar = 1;
	}
	
	/* option to have heating vary with time */
	if( p.nMatch( "TIME" ) )
		hextra.lgTurbHeatVaryTime = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		if( hextra.lgHextraSS )
		{
			fprintf(ioQQQ,"Sorry, HEXTRA SS command does not now support vary option.\n");
			cdEXIT(EXIT_FAILURE);
		}
		/* 1 to 3 options on hextra vary */
		optimize.nvarxt[optimize.nparm] = nPar;
		optimize.vparm[0][optimize.nparm] = log10(hextra.TurbHeat);
		optimize.vparm[1][optimize.nparm] = par1;
		optimize.vparm[2][optimize.nparm] = par2;
		
		/* the keyword LOG is not used above, but is checked elsewhere */
		strcpy( optimize.chVarFmt[optimize.nparm], "HEXTra %f LOG " );
		strcat( optimize.chVarFmt[optimize.nparm], chHextraScale );
		while( nPar > 1 )
		{
			/* add extra spots to write parameters */
			--nPar;
			strcat( optimize.chVarFmt[optimize.nparm], " %f" );
		}
		if( hextra.lgTurbHeatVaryTime )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* the increment in the first steps away from the original value */
		optimize.vincr[optimize.nparm] = 0.1f;
		++optimize.nparm;
	}
}

void ParseConvHighT(Parser &)
{
	thermal.lgTeHigh = true;
}
void ParseHydrogen(Parser &)
{
	DEBUG_ENTRY("ParseHydrogen()");
	fprintf(ioQQQ," Sorry, this command has been replaced with the ATOM H-LIKE command.\n");
	cdEXIT(EXIT_FAILURE);
}
void ParseInitCount(Parser &p)
{
	DEBUG_ENTRY("ParseInitCount()");
	/* read cloudy.ini initialization file
	 * need to read in file every time, since some vars
	 * get reset in zero - would be unsafe not to read in again */
	ParseInit(p);
	
	/* check that only one init file was in the input stream -
	 * we cannot now read more than one */
	++p.m_nInitFile;
	if( p.m_nInitFile > 1 )
	{
		fprintf( ioQQQ, 
					" This is the second init file, I can only handle one.\nSorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	/* we will put the data for the ini file at the end of the array of
	 * line images, this is the increment for stuffing the lines in - negative */
	input.iReadWay = -1;
	
	/* initialize array reader, this sub does nothing but set
	 * initial value of a variable, depending on value of iReadWay */
	input.init();
}
void ParseIntensity(Parser &p)
{
	DEBUG_ENTRY("ParseIntensity()");
	/* intensity of this continuum source */
	if( (p.m_nqh) >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* silly, but calms down the lint */
	ASSERT( (p.m_nqh) < LIMSPC );
	strcpy( rfield.chRSpec[(p.m_nqh)], "SQCM" );
	rfield.totpow[(p.m_nqh)] = p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("intensity");
	
	/* set radius to very large value if not already set */
	/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
	if( !radius.lgRadiusKnown )
	{
		/* set radius to large value */
		radius.Radius = pow(10.,radius.rdfalt);
	}
	if( p.nMatch("LINE") )
	{
		/* silly, but calms down the lint */
		ASSERT( (p.m_nqh) < LIMSPC );
		/* option for linear input parameter */
		rfield.totpow[(p.m_nqh)] = log10(rfield.totpow[(p.m_nqh)]);
	}
	strcpy( rfield.chSpNorm[(p.m_nqh)], "LUMI" );
	/* ParseRangeOption in readsun */
	ParseRangeOption(p);
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		/* create the format we will use to write the command */
		strcpy( optimize.chVarFmt[optimize.nparm], "INTENSITY %f LOG range %f %f" );
		if( rfield.lgTimeVary[(p.m_nqh)] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* array index for this command */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)rfield.totpow[(p.m_nqh)];
		optimize.vincr[optimize.nparm] = 0.5;
		/* range option, but cannot be varied */
		optimize.vparm[1][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][0]);
		optimize.vparm[2][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][1]);
		optimize.nvarxt[optimize.nparm] = 3;
		++optimize.nparm;
	}
	++(p.m_nqh);
}
void ParseIterations(Parser &p)
{
	/* iterate to get optical depths and diffuse field properly */
	iterations.itermx = (long int)p.FFmtRead() - 1;
	iterations.itermx = MAX2(iterations.itermx,1);
	/* >>chng 06 mar 19, malloc itrdim arrays 
	 * iterations.iter_malloc is -1 before space allocated and set to
	 * itermx if not previously set */
	if( iterations.itermx > iterations.iter_malloc - 1 )
	{
		long int iter_malloc_save = iterations.iter_malloc;
		/* add one for mismatch between array dim and iterations defn,
		 * another few for safety */
		iterations.iter_malloc = iterations.itermx+3;
		iterations.IterPrnt = (long int*)REALLOC(iterations.IterPrnt ,
															  (size_t)iterations.iter_malloc*sizeof(long int) );
		geometry.nend = (long int*)REALLOC(geometry.nend ,
													  (size_t)iterations.iter_malloc*sizeof(long int) );
		radius.StopThickness = (double*)REALLOC(radius.StopThickness ,
													(size_t)iterations.iter_malloc*sizeof(double) );
		/* >>chng 06 jun 07, need to set IterPrnt, thickness, and nend */
		for(long int j=iter_malloc_save; j<iterations.iter_malloc; ++j )
		{
			iterations.IterPrnt[j] = iterations.IterPrnt[iter_malloc_save-1];
			geometry.nend[j] = geometry.nend[iter_malloc_save-1];
			radius.StopThickness[j] = radius.StopThickness[iter_malloc_save-1];
		}
	}
	
	if( p.nMatch("CONV") )
	{
		/* option to keep iterating until it converges, checks on convergence
		 * are done in update, and checked again in prtcomment */
		conv.lgAutoIt = true;
		/* above would have been legitimate setting of ITERMX, set to default 10 */
		if( p.lgEOL() )
		{
			iterations.itermx = 10 - 1;
		}
		realnum a = (realnum)p.FFmtRead();
		/* change convergence criteria, preset in zero */
		if( !p.lgEOL() )
		{
			conv.autocv = a;
		}
	}
}	
void ParseL_nu(Parser &p)
{
	/* this is the specific luminosity at nu
	 * following says n(nu) not nuF(nu) */
	bool lgNu2 = false;
	/* in reads2 */
	ParseF_nu(p,"4 PI",lgNu2);
}
void ParseLaser(Parser &p)
{
	DEBUG_ENTRY("ParseLaser()");
	/* mostly one frequency (a laser) to check gamma's,*/
	/* say the continuum type */
	strcpy( rfield.chSpType[rfield.nShape], "LASER" );
	
	/* scan off the laser's energy ion Rydbergs */
	rfield.slope[rfield.nShape] = p.FFmtRead();
	
	/* negative energies are logs */
	if( rfield.slope[rfield.nShape] <= 0. )
	{
		rfield.slope[rfield.nShape] =
			pow(10.,rfield.slope[rfield.nShape]);
	}
	if( p.lgEOL() )
	{
		p.NoNumb("energy");
	}
	
	/* next number is optional relative width of laser */
	rfield.cutoff[rfield.nShape][0] = p.FFmtRead();
	if( p.lgEOL() )
	{
		/* default width is 0.05 */
		rfield.cutoff[rfield.nShape][0] = 0.05;
	}
	
	/* increment counter making sure we still have space enough */
	++rfield.nShape;
	if( rfield.nShape >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
				cdEXIT(EXIT_FAILURE);
	}
}
void ParseLuminosity(Parser &p)
{
	DEBUG_ENTRY("ParseLuminosity()");
	/* luminosity of this continuum source */
	if( (p.m_nqh) >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	strcpy( rfield.chRSpec[(p.m_nqh)], "4 PI" );
	rfield.totpow[(p.m_nqh)] = p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("luminosity");
	}
	if( p.nMatch("LINE") )
	{
		/* option for linear input parameter */
		rfield.totpow[(p.m_nqh)] = log10(rfield.totpow[(p.m_nqh)]);
	}
	
	strcpy( rfield.chSpNorm[(p.m_nqh)], "LUMI" );
	
	/* the solar option - number is log total solar luminosity */
	if( p.nMatch("SOLA") )
	{
		/* option to use log of total luminosity in solar units */
		rfield.range[(p.m_nqh)][0] = rfield.emm;
		rfield.range[(p.m_nqh)][1] = rfield.egamry;
		/* log of solar luminosity in erg s-1 */
		rfield.totpow[(p.m_nqh)] += 33.5827f;
	}
	else
	{
		/* check if range is present and parse it if it is - ParseRangeOption in readsun 
		 * this includes TOTAL keyword for total luminosity */
		ParseRangeOption(p);
	}
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "LUMINOSITY %f LOG range %f %f" );
		if( rfield.lgTimeVary[(p.m_nqh)] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)rfield.totpow[(p.m_nqh)];
		optimize.vincr[optimize.nparm] = 0.5;
		/* range option in, but cannot be varied */
		optimize.vparm[1][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][0]);
		optimize.vparm[2][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][1]);
		optimize.nvarxt[optimize.nparm] = 3;
		++optimize.nparm;
	}
	++(p.m_nqh);
}
void ParseNeutrons(Parser &p)
{
	/* heating and ionization due to fast neutrons */
	hextra.lgNeutrnHeatOn = true;
	
	/* first number is neutron luminosity
	 * relative to bolometric luminosity */
	hextra.frcneu = (realnum)p.FFmtRead();
	if( p.lgEOL() )
		p.NoNumb("neutron luminosity");
	
	/* save as log of fraction */
	if( hextra.frcneu > 0. )
		hextra.frcneu = (realnum)log10(hextra.frcneu);
	
	/* second number is efficiency */
	hextra.effneu = (realnum)p.FFmtRead();
	if( p.lgEOL() )
	{
		hextra.effneu = 1.0;
	}
	else
	{
		if( hextra.effneu <= 0. )
			hextra.effneu = (realnum)pow((realnum)10.f,hextra.effneu);
	}
}
void ParseNuF_nu(Parser &p)
{
	/* flux density of this continuum source, at optional frequency
	 *  expressed as product nu*f_nu */
	bool lgNu2 = true;
	/* in reads2 */
	ParseF_nu(p,"SQCM",lgNu2);
}
void ParseNuL_nu(Parser &p)
{
	/* specific luminosity density of this continuum source, at opt freq
	 * expressed as product nu*f_nu */
	bool lgNu2 = true;
	/* in reads2 */
	ParseF_nu(p,"4 PI",lgNu2);
}
void ParsePhi(Parser &p)
{
	DEBUG_ENTRY("ParsePhi()");
	/* enter phi(h), the number of h-ionizing photons/cm2 */
	if( (p.m_nqh) >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* silly, but calms down the lint */
	ASSERT( (p.m_nqh) < LIMSPC );
	strcpy( rfield.chRSpec[(p.m_nqh)], "SQCM" );
	strcpy( rfield.chSpNorm[(p.m_nqh)], "PHI " );
	rfield.totpow[(p.m_nqh)] = p.FFmtRead();
	if( p.lgEOL() )
	{
		p.NoNumb("number of h-ionizing photons");
	}
	/* set radius to very large value if not already set */
	/* >>chng 01 jul 24, from Radius == 0 to this, as per PvH comments */
	if( !radius.lgRadiusKnown )
	{
		/* set radius to large value */
		radius.Radius = pow(10.,radius.rdfalt);
	}
	/* check if continuum so intense that crash is likely */
	if( rfield.totpow[(p.m_nqh)] > 29. )
	{
		fprintf( ioQQQ, " Is the flux for this continuum correct?\n" );
		fprintf( ioQQQ, " It appears too bright to me.\n" );
	}
	/* ParseRangeOption in readsun xx parse_rangeoption*/
	ParseRangeOption(p);
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "phi(h) %f LOG range %f %f" );
		if( rfield.lgTimeVary[(p.m_nqh)] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)rfield.totpow[(p.m_nqh)];
		optimize.vincr[optimize.nparm] = 0.5;
		/* range option in, but cannot be varied */
		optimize.vparm[1][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][0]);
		optimize.vparm[2][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][1]);
		optimize.nvarxt[optimize.nparm] = 3;
		++optimize.nparm;
	}
	++(p.m_nqh);
}
void ParseQH(Parser &p)
{
	DEBUG_ENTRY("ParseQH()");
	/* log of number of ionizing photons */
	if( (p.m_nqh) >= LIMSPC )
	{
		/* too many continua were entered */
		fprintf( ioQQQ, " Too many continua entered; increase LIMSPC\n" );
		cdEXIT(EXIT_FAILURE);
	}
	
	/* silly, but calms down the lint */
	ASSERT( (p.m_nqh) < LIMSPC );
	strcpy( rfield.chRSpec[(p.m_nqh)], "4 PI" );
	strcpy( rfield.chSpNorm[(p.m_nqh)], "Q(H)" );
	rfield.totpow[(p.m_nqh)] = p.FFmtRead();
	if( rfield.totpow[(p.m_nqh)] > 100. && called.lgTalk )
	{
		fprintf( ioQQQ, " Is this reasonable?\n" );
	}
	if( p.lgEOL() )
	{
		p.NoNumb("number of ionizing photons");
	}
	/* ParseRangeOption in readsun */
	ParseRangeOption(p);
	
	/* >>chng 06 mar 22, add time option to vary only some continua with time */
	if( p.nMatch( "TIME" ) )
		rfield.lgTimeVary[(p.m_nqh)] = true;
	
	/* vary option */
	if( optimize.lgVarOn )
	{
		strcpy( optimize.chVarFmt[optimize.nparm], "Q(H) %f LOG range %f %f" );
		if( rfield.lgTimeVary[(p.m_nqh)] )
			strcat( optimize.chVarFmt[optimize.nparm], " TIME" );
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		optimize.vparm[0][optimize.nparm] = (realnum)rfield.totpow[(p.m_nqh)];
		optimize.vincr[optimize.nparm] = 0.5;
		/* range option in, but cannot be varied */
		optimize.vparm[1][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][0]);
		optimize.vparm[2][optimize.nparm] = (realnum)log10(rfield.range[(p.m_nqh)][1]);
		optimize.nvarxt[optimize.nparm] = 3;
		++optimize.nparm;
	}
	/* increment number of luminosity sources specified */
	++(p.m_nqh);
}
void ParseRoberto(Parser &)
{
	/* this is the Roberto Terlivich command, no telling if it still works */
	radius.dRadSign = -1.;
}
void ParseSpecial(Parser &)
{
	DEBUG_ENTRY("ParseSpecial()");
	/* special key, can do anything */
	cdEXIT(EXIT_FAILURE);
}
void ParseTauMin(Parser &p)
{
	/* taumin command minimum optical depths for lines dafault 1e-20 */
	opac.taumin = (realnum)pow(10.,p.FFmtRead());
	if( p.lgEOL() )
		p.NoNumb("minimum optical depth");	
}
void ParseTitle(Parser &p)
{
	/* read in title of model starting in col 5 -- prefer to get string
		in quotes, but for the moment if it's not present take what you
		can get */
	if (p.GetQuote(input.chTitle,false) != 0)
		strcpy( input.chTitle , p.getRawTail().c_str()+1 );
}	
void ParseTolerance(Parser &)
{
	DEBUG_ENTRY("ParseTolerance()");
	fprintf(ioQQQ,
			  "Sorry, this command has been replaced with the SET TEMPERATURE TOLERANCE command.\n");
	cdEXIT(EXIT_FAILURE);
}
void ParseVLaw(Parser &p)
{
	/* velocity power law as a function of radius. */
	DoppVel.TurbVelLaw = (realnum)p.FFmtRead();
	
	DoppVel.lgTurbLawOn = true;
	/** \todo 2 is there a need to keep this negative? */
	/* for now, insist upon non-positive power,
	 * so that velocity decreases with increasing radius. */
	ASSERT( DoppVel.TurbVelLaw <= 0.f );
}
void ParseTurbulence(Parser &p)
{
	DEBUG_ENTRY("ParseTurbulence()");
	string ExtraPars;

	if( p.nMatch("EQUIPART") )
	{
		/* turbulence equipartition option */
		DoppVel.lgTurbEquiMag = true;

		/* this is optional F parameter from equation 34 of 
		 *>>refer	pressure	turb	Heiles, C. & Troland, T.H. 2005, ApJ, 624, 773 
		 * turbulent energy density will be rho F v^2 / 2 */
		DoppVel.Heiles_Troland_F = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			/* this is the default value of 3 for isotropic turbulence */
			DoppVel.Heiles_Troland_F = 3.f;
		}
	}
	else
	{
		/* line has turbulent velocity in km/s */
		DoppVel.lgTurbEquiMag = false;
		DoppVel.TurbVel = (realnum)p.FFmtRead();
		if( p.lgEOL() )
			p.NoNumb("microturbulent velocity");	

		if( p.nMatch(" LOG") )
		{
			if( DoppVel.TurbVel > 32. )
			{
				fprintf( ioQQQ, "PROBLEM the log of the turbulence is "
					 "%.2e - I cannot handle a number this big.\n",
					 DoppVel.TurbVel );
				fprintf( ioQQQ, " The line image was\n" );
				p.PrintLine(ioQQQ);
				fprintf( ioQQQ, " Sorry.\n" );
				cdEXIT(EXIT_FAILURE);
			}
			DoppVel.TurbVel = (realnum)pow((realnum)10.f,DoppVel.TurbVel);
		}
		/* now convert from km/s to cm/s */
		DoppVel.TurbVel *= 1e5;

		if( DoppVel.TurbVel < 0. )
		{
			fprintf( ioQQQ, " PROBLEM: the turbulent velocity needs to be > 0, but this was entered: %e\n",
				 DoppVel.TurbVel );
			fprintf( ioQQQ, " Bailing out. Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}
		else if( DoppVel.TurbVel >= SPEEDLIGHT )
		{
			fprintf( ioQQQ, " PROBLEM: A turbulent velocity greater than speed of light is not allowed, this was entered: %e\n",
				 DoppVel.TurbVel );
			fprintf( ioQQQ, " Bailing out. Sorry.\n" );
			cdEXIT(EXIT_FAILURE);
		}

		/* this is optional F parameter from equation 34 of 
		 *>>refer	pressure	turb	Heiles, C. & Troland, T.H. 2005, ApJ, 624, 773 
		 * turbulent energy density will be rho F v^2 / 2 */
		DoppVel.Heiles_Troland_F = (realnum)p.FFmtRead();
		if( p.lgEOL() )
		{
			/* this is the default value of 3 for isotropic turbulence */
			DoppVel.Heiles_Troland_F = 3.f;
		}

		/* option to have turbulence be dissipative, keyword is dissipate,
		 * argument is log of scale length in cm */
		if( p.nMatch("DISS") )
		{
			DoppVel.DispScale = (realnum)pow(10., p.FFmtRead() );
			if( p.lgEOL() )
				p.NoNumb("turbulence dissipation scale");
			ExtraPars += " DISSIPATE %f";
		}
	}

	/* include turbulent pressure in equation of state?
	 * >>chng 06 mar 24, include turbulent pressure in gas equation of state by default,
	 * but not if NO PRESSURE occurs */
	if( p.nMatch(" NO ") && p.nMatch("PRES") )
	{
		DoppVel.lgTurb_pressure = false;
		ExtraPars += " NO PRESSURE";
	}
	else
	{
		/* default is to include turbulent pressure in equation of state */
		DoppVel.lgTurb_pressure = true;
	}

	/* vary option */
	if( optimize.lgVarOn && !p.nMatch("EQUIPART") )
	{
		/* only one parameter to vary */
		optimize.nvarxt[optimize.nparm] = 2;
		// the line image
		strcpy( optimize.chVarFmt[optimize.nparm], "TURBULENCE= %f LOG %f" );
		strcat( optimize.chVarFmt[optimize.nparm], ExtraPars.c_str() );
		
		/* pointer to where to write */
		optimize.nvfpnt[optimize.nparm] = input.nRead;
		/* turbulent velocity */
		optimize.vparm[0][optimize.nparm] = log10(DoppVel.TurbVel/1e5f);
		optimize.vparm[1][optimize.nparm] = DoppVel.Heiles_Troland_F;
		if( p.nMatch("DISS") )
		{
			optimize.nvarxt[optimize.nparm] = 3;
			optimize.vparm[2][optimize.nparm] = log10(DoppVel.DispScale);
		}
		optimize.vincr[optimize.nparm] = 0.1f;
		++optimize.nparm;
	}

	/* set initial turbulence */
	DoppVel.TurbVelZero = DoppVel.TurbVel;
}
