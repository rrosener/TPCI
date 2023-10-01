/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*cloudy the main routine, this IS Cloudy, return 0 normal exit, 1 error exit,
 * called by maincl when used as standalone program */
/*BadStart announce that things are so bad the calculation cannot even start */
#include "cddefines.h"
#include "save.h"
#include "noexec.h"
#include "lines.h"
#include "abund.h"
#include "continuum.h"
#include "warnings.h"
#include "atmdat.h"
#include "prt.h"
#include "conv.h"
#include "parse.h"
#include "dynamics.h"
#include "init.h"
#include "opacity.h"
#include "state.h"
#include "rt.h"
#include "monitor_results.h"
#include "zones.h"
#include "iterations.h"
#include "plot.h"
#include "radius.h"
#include "grid.h"
#include "cloudy.h"
#include "ionbal.h"
#include "dense.h"
#include "energy.h"
#include "called.h"

/*BadStart announce that things are so bad the calculation cannot even start */
STATIC void BadStart();

/* returns 1 if disaster strikes, 0 if everything appears ok */
bool cloudy()
{
	DEBUG_ENTRY( "cloudy()" );

	bool lgOK;

	/* 
	 * this is the main routine to drive the modules that compute a model
	 * when cloudy is used as a stand-alone program 
	 * the main program (maincl) calls cdInit then cdDrive
	 * this sub is called by cdDrive which returns upon exiting
	 *
	 * this routine uses the following variables:
	 *
	 * nzone 
	 * this is the zone number, and is incremented here
	 * is zero during search phase, 1 for first zone at illuminated face
	 * logical function iter_end_check returns a true condition if NZONE reaches
	 * NEND(ITER), the limit to the number of zones in this iteration
	 * nzone is totally controlled in this subroutine
	 *
	 * iteration 
	 * this is the iteration number, it is 1 for the first iteration
	 * and is incremented in this subroutine at the end of each iteration
	 *
	 * iterations.itermx
	 * this is the limit to the number of iterations, and is entered
	 * by the user
	 * This routine returns when iteration > iterations.itermx
	 */

	/* nzone is zero while initial search for conditions takes place 
	 * and is incremented to 1 for start of first zone */
	nzone = 0;
	fnzone = 0.; 

	/* iteration is iteration number, calculation is complete when iteration > iterations.itermx */
	iteration = 1;

	/* this initializes variables at the start of each simulation
	 * in a grid, before the parser is called - this must set any values
	 * that may be changed by the command parser */
	InitDefaultsPreparse();

	/* scan in and parse input data */
	ParseCommands();

	/* fix abundances of elements, in abundances.cpp */
	AbundancesSet();

	ASSERT(lgElemsConserved());

	/* one time creation of some variables */
	InitCoreloadPostparse();

	/* initialize vars that can only be done after parsing the commands 
	 * sets up CO network since this depends on which elements are on
	 * and whether chemistry is enabled */
	InitSimPostparse();

	/* ContCreateMesh calls fill to set up continuum energy mesh if first call, 
	 * otherwise reset to original mesh.
	 * This is AFTER ParseCommands so that
	 * path and mesh size can be set with commands before mesh is set */
	ContCreateMesh();

	/* create several data sets by allocating memory and reading in data files, 
	 * but only if this is first call */
	atmdat_readin();

	/* fix pointers to ionization edges and frequency array
	 * calls iso_create
	 * this routine only returns if this is a later call of code */
	ContCreatePointers();

	/* Badnell_rec_init This code is written by Terry Yun, 2005 *
	 * It reads dielectronic recombination rate coefficient fits into 3D arrays */
	Badnell_rec_init();

	ASSERT(lgElemsConserved());

	/* set continuum normalization, continuum itself, and ionization stage limits */
	ContSetIntensity();

	ASSERT(lgElemsConserved());

	/* print header */
	PrtHeader();

	/* this is an option to stop after initial setup */
	if( noexec.lgNoExec )
		return false;

	/* guess some outward optical depths, set inward optical depths, 
	 * also calls RT_line_all so escape probabilities ok before printout of trace */
	RT_tau_init();

	/* generate initial set of opacities, but only if this is the first call 
	 * in this core load
	 * grains only exist AFTER this routine exits */
	OpacityCreateAll();

	ASSERT(lgElemsConserved());

	/* this checks that various parts of the code still work properly */
	SanityCheck("begin");

	/* option to recover state from previous calculation */
	if( state.lgGet_state )
		state_get_put( "get" );

	ASSERT(lgElemsConserved());

	/* find the initial temperature, punt if initial conditions outside range of code,
	 * abort condition set by flag lgAbort */
	if( ConvInitSolution() )
	{
		// create line stacks for possible printout before bailing out
		LineStackCreate();
		BadStart();
		return true;
	}
	// create line stacks ...
	LineStackCreate();

	/* set thickness of first zone */
	radius_first();

	/* find thickness of next zone */
	if( radius_next() )
	{
		BadStart();
		return true;
	}

	/* set up some zone variables, correct continuum for sphericity, 
	 * after return, radius is equal to the inner radius, not outer radius
	 * of this zone */
	ZoneStart("init");

	/* print all abundances, gas phase and grains, in abundances.c */
	/* >>chng 06 mar 05, move AbundancesPrt() after ZoneStart("init") so that
	 * GrnVryDpth() gives correct values for the inner face of the cloud, PvH */
	AbundancesPrt();

	/* this is an option to stop after printing header only */
	if( prt.lgOnlyHead )
		return false;

	plot("FIRST");

	/* outer loop is over number of iterations
	 * >>chng 05 mar 12, chng from test on true to not aborting */
	while( !lgAbort )
	{
		IterStart();
		nzone = 0;
		fnzone = 0.;

		/* loop over zones across cloud for this iteration, 
		 * iter_end_check checks whether model is complete and this iteration is done
		 * returns true is current iteration is complete 
		 * calls PrtZone to print zone results */
		while( !iter_end_check() )
		{
			/* the zone number, 0 during search phase, first zone is 1 */
			++nzone;
			/* this is the zone number plus the number of calls to bottom solvers
			 * from top pressure solver, divided by 100 */
			fnzone = (double)nzone;

			/* use changes in opacity, temperature, ionization, to fix new dr for next zone */
			/* >>chng 03 oct 29, move radius_next up to here from below, so that
			 * precise correct zone thickness will be used in current zone, so that
			 * column density and thickness will be exact 
			 * abort condition is possible */
			if( radius_next() )
				break;

			/* following sets zone thickness, drad, to drnext */
			/* set variable dealing with new radius, in zones.c */
			ZoneStart("incr");

			/* converge the pressure-temperature-ionization solution for this zone 
			 * NB ignoring return value - should be ok (return 1 for abort) 
			 * we can't break here in case of abort since there is still housekeeping
			 * that must be done in following routines */
			ConvPresTempEdenIoniz();

			/* generate diffuse emission from this zone, add to outward & reflected beams */
			RT_diffuse();

			/* do work associated with incrementing this radius, 
			 * total attenuation and dilution of radiation fields takes place here,
			 * reflected continuum incremented here
			 * various mean and extremes of solution are also remembered here */
			radius_increment();

			/* attenuation of diffuse and beamed continua */
			RT_continuum();

			/* increment optical depths 
			 * final evaluation of line rates and cooling */
			RT_tau_inc();

			/* fill in emission line array, adds outward lines */
			/* >>chng 99 dec 29, moved to here from below RT_tau_inc, 
			 * lines adds lines to outward beam,
			 * and these are attenuated in radius_increment */
			lines();

			/* possibly save out some results from this zone */
			SaveDo("MIDL");

			/* do some things to finish out this zone */
			ZoneEnd();

			// default false due to slow time - set true with set check energy every zone
			// to allow per zone check of energy conservation, relatively slow
			// when chianti is fully enabled
			if( continuum.lgCheckEnergyEveryZone )
			{
				/* Ensure that energy has been conserved.  Complain and punt if not */
				if( !lgConserveEnergy() )
				{
					fprintf( ioQQQ, " PROBLEM DISASTER Energy was not conserved at zone %li\n", nzone );
					ShowMe();
					lgAbort = true;
				}
			}
		}
		/* end loop over zones */

		/* close out this iteration, in startenditer.c */
		IterEnd();

		/* print out some comments, generate warning and cautions*/
		PrtComment();

		/* save stuff only needed on completion of this iteration */
		SaveDo("LAST" );

		/* second call to plot routine, to complete plots for this iteration */
		plot("SECND");

		/* print out the results */
		PrtFinal();

		/*ConvIterCheck check whether model has converged or more iterations
		 * are needed - implements the iter to convergence command 
		 * acts by setting iterations.itermx to iteration if converged*/
		ConvIterCheck();

		/* option to save state */
		if( state.lgPut_state )
			state_get_put( "put" );

		/* >>chng 06 mar 03 move block to here, had been after PrtFinal,
		 * so that save state will save reset optical depths */
		/* this is the normal exit, occurs if we reached limit to number of iterations,
		 * or if code has set busted */
		/* >>chng 06 mar 18, add flag for time dependent simulation completed */
		if( iteration > iterations.itermx || lgAbort || dynamics.lgStatic_completed )
			break;

		/* reset limiting and initial optical depth variables */
		RT_tau_reset();

		/* increment iteration counter */
		++iteration;

		/* reinitialize some variables to initial conditions at previous first zone
		 * routine in startenditer.c */
		IterRestart();

		/* reset zone number to 0 - done here since this is only routine that sets nzone */
		nzone = 0;
		fnzone = 0.;

		ZoneStart("init");

		/* find new initial temperature, punt if initial conditions outside range of code,
		 * if ConvInitSolution() fails, lgAbort will always be true, so we check that */
		if( ConvInitSolution() )
			break;
	}

	CloseSaveFiles( false );

	/* this checks that various parts of the code worked properly */
	SanityCheck("final");

	if( called.lgTalk )
	{
		fprintf(ioQQQ,"---------------Convergence statistics---------------\n"); 
		fprintf(ioQQQ,"%10.3g mean iterations/state convergence\n",((double)conv.getCounter(CONV_BASE_LOOPS))/
			(MAX2((double)conv.getCounter(CONV_BASE_CALLS),1.0)));
		fprintf(ioQQQ,"%10.3g mean cx acceleration loops/iteration\n",((double)conv.getCounter(CONV_BASE_ACCELS))/
			(MAX2((double)conv.getCounter(CONV_BASE_LOOPS),1.0)));
		fprintf(ioQQQ,"%10.3g mean iso convergence loops/ion solve\n",((double)conv.getCounter(ISO_LOOPS))/
			(MAX2((double)conv.getCounter(ION_SOLVES),1.0)));
		fprintf(ioQQQ,"%10.3g mean steps/chemistry solve\n",((double)conv.getCounter(MOLE_SOLVE_STEPS))/
			(MAX2((double)conv.getCounter(MOLE_SOLVE),1.0)));
		fprintf(ioQQQ,"%10.3g mean step length searches/chemistry step\n",((double)conv.getCounter(NEWTON_LOOP))/
			(MAX2((double)conv.getCounter(NEWTON),1.0)));
		fprintf(ioQQQ,"----------------------------------------------------\n\n");
	}
	
	/* check whether any asserts were present and failed.  
	 * return is true if ok, false if not.  routine also checks
	 * number of warnings and returns false if not ok */
	lgOK = lgCheckMonitors(ioQQQ);

	if( lgOK && !warnings.lgWarngs && !lgAbort )
	{
		/* no failed asserts or warnings */
		return false;
	}
	else
	{
		/* there were failed asserts or warnings */
		return true;
	}
}

/*BadStart announce that things are so bad the calculation cannot even start */
STATIC void BadStart()
{
	char chLine[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "BadStart()" );

	/* initialize the line saver */
	wcnint();
	sprintf( warnings.chRgcln[0], "   Calculation stopped because initial conditions out of bounds." );
	sprintf( chLine, " W-Calculation could not begin." );
	warnin(chLine);

	if( grid.lgGrid )
	{
		/* possibly save out some results from this zone when doing grid
		 * since grid save files expect something at this grid point */
		SaveDo("MIDL");
		SaveDo("LAST" );
	}
	return;
}
