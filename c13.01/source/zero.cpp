/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*zero zero out or initialize variables, called by cdInit, but also by optimize_func during optimization,
 * this is called before any commands are parsed, called one time per model, at very start */
/*rfield_optac_zero zero out rfield arrays between certain limits */
#include "cddefines.h"
#include "physconst.h"
#include "iterations.h"
#include "hydrogenic.h"
#include "oxy.h"
#include "doppvel.h"
#include "dense.h"
#include "hextra.h"
#include "grains.h"
#include "magnetic.h"
#include "state.h"
#include "rt.h"
#include "he.h"
#include "struc.h"
#include "h2.h"
#include "co.h"
#include "coolheavy.h"
#include "lines.h"
#include "dynamics.h"
#include "carb.h"
#include "mean.h"
#include "atomfeii.h"
#include "iso.h"
#include "conv.h"
#include "geometry.h"
#include "timesc.h"
#include "peimbt.h"
#include "ionbal.h"
#include "continuum.h"
#include "atmdat.h"
#include "mole.h"
#include "ca.h"
#include "input.h"
#include "atoms.h"
#include "pressure.h"
#include "numderiv.h"
#include "colden.h"
#include "yield.h"
#include "hmi.h"
#include "rfield.h"
#include "abund.h"
#include "radius.h"
#include "opacity.h"
#include "secondaries.h"
#include "called.h"
#include "phycon.h"
#include "warnings.h"
#include "thermal.h"
#include "cooling.h"
#include "fe.h"
#include "hyperfine.h"
#include "init.h"
#include "dark_matter.h"

// //////////////////////////////////////////////////////////////////////////
//
//
// NB DO NOT ADD VARIABLES TO THIS FILE!  THE GOAL IS TO REMOVE THIS FILE
// initialization of variables should be done in one of the ini_*.cpp routines
//
//
// //////////////////////////////////////////////////////////////////////////

/* zero out or initialize variables, called by cdInit, but also by optimize_func 
 * during optimization, called before command parser, one time per model, 
 * in a grid one time per grid point (so called nGrid times), 
 * only one time in multi-iteration models */
void zero(void)
{
	//long int i;

	/* this is used to signify the first call to this routine.  At that
	 * stage some memory has not been allocated so must not initialize,
	 * set false at very end of this routine */
	static bool lgFirstCall = true;

	DEBUG_ENTRY( "zero()" );

	/* this routine is called exactly one time at the start of
	 * the calculation of a single model.  When the code is used as a subroutine
	 * this routine is called one time for each model.  It is called before
	 * the boundary conditions are read in, and is never called again
	 * during that calculation of the one model.  
	 * All default variables that must be initialized before a calculation starts
	 * must appear in the routine.  In a grid they are reset for each model
	 */

	/* parameters having to do with magnetic field */
	Magnetic_init();

	/* set all initial abundances */
	AbundancesZero();

	/* zero out parameters needed by large FeII atom */
	FeIIZero();

	/* zero out warnings, cautions, notes, etc */
	wcnint();

	/* this is number of iterations that have been malloced - we could 
	 * increase this if more iterations are needed */
	iterations.iter_malloc = 200;
	/* >>chng 06 jun 27, only malloc on first call - memory leak */
	if( lgFirstCall)
	{
		iterations.IterPrnt = (long int*)MALLOC( (size_t)iterations.iter_malloc*sizeof(long int) );
		geometry.nend = (long int*)MALLOC( (size_t)iterations.iter_malloc*sizeof(long int) );
		radius.StopThickness = (double*)MALLOC( (size_t)iterations.iter_malloc*sizeof(double) );
		radius.StopRadius = (double*)MALLOC( (size_t)iterations.iter_malloc*sizeof(double) );
	}
	for( long i=0; i < iterations.iter_malloc; i++ )
	{
		iterations.IterPrnt[i] = 10000;
	}
	iterations.itermx = 0;
	/* this implements set coverage command */
	iterations.lgConverge_set = false;
	iteration = 0;

	/* limits for highest and lowest stages of ionization in TrimStage */
	ionbal.trimhi = 1e-6;
	ionbal.lgTrimhiOn = true;
	ionbal.trimlo = 1e-10;

	hyperfine.lgLya_pump_21cm = true;

	/* variable to do with geometry */
	geometry.nprint = 1000;
	geometry.lgZoneSet = false;
	geometry.lgZoneTrp = false;
	geometry.lgEndDflt = true;

	/* some variables for saving the codes' state */
	state.lgGet_state = false;
	state.lgPut_state = false;
	state.lgState_print = false;

	/* this is default number of zones
	 * >>chng 96 jun 5, from 400 to 500 for thickest corners4 grid */
	/* >>chng 04 jan 30, from 600 to 800, code uses finer zoning today */
	/* >>chng 04 dec 24, from 800 to 1400, so that HII region - molecular cloud
	 * sims do not need set nend - all sims in test suite will run ok without set nend */
	geometry.nEndDflt = 1400;

	for( long i=0; i < iterations.iter_malloc; i++ )
	{
		geometry.nend[i] = geometry.nEndDflt;
		/*>>chng 03 nov 13, from 1e30 to 1e31, because default inner radius raised to 1e30 */
		radius.StopThickness[i] = 1e31;
		radius.StopRadius[i] = -1.;
	}

	geometry.fiscal = 1.;
	geometry.FillFac = 1.;
	geometry.filpow = 0.;

	/* default is open geometry, not sphere */
	geometry.lgSphere = false;
	/* the radiative transport covering factor */
	geometry.covrt = 0.;
	/* the geometric covering factor */
	geometry.covgeo = 1.;
	/* default is expanding when geometry set */
	geometry.lgStatic = false;
	/* option to tell code not to complain when geometry static done without iterating,
	 * set with (OK) option on geometry command */
	geometry.lgStaticNoIt = false;
	/* this is exponent for emissivity contributing to observed luminosity, r^2.
	 * set to 1 with aperture slit, to 0 with aperture beam command */
	geometry.iEmissPower = 2;

	/* this counts number of times ionize is called by PressureChange, in current zone
	 * these are reset here, so that we count from first zone not search */
	conv.nPres2Ioniz = 0;

	/* clear flag indicating the ionization convergence failures 
	 * have ever occurred in current zone
	conv.lgConvIonizThisZone = false; */

	conv.resetCounters();

	/* general abort flag */
	lgAbort = false;

	/* cooling tolerance heating tolerance - allowed error in heating -  cooling balance */
	/*conv.HeatCoolRelErrorAllowed = 0.02f;*/
	/* >>chng 04 sep 25, change te tolerance from 0.02 to 4x smaller, 0.005, drove instabilities
	 * in chemistry */
	conv.HeatCoolRelErrorAllowed = 0.005f;

	/* this is the default allowed relative error in the electron density */
	conv.EdenErrorAllowed = 1e-2;

	conv.IonizErrorAllowed = 1e-2;

	conv.dCmHdT = 0.;

	conv.LimFail = 20;
	conv.lgMap = false;

	/* this counts how many times ionize is called in this model after startr,
	 * and is flag used by ionize to understand it is being called the first time*/
	conv.nTotalIoniz = 0;
	/* these are variables to remember the biggest error in the
	 * electron density, and the zone in which it occurred */
	conv.BigEdenError = 0.;
	conv.AverEdenError = 0.;
	conv.BigHeatCoolError = 0.;
	conv.AverHeatCoolError = 0.;
	conv.BigPressError = 0.;
	conv.AverPressError = 0.;
	strcpy( conv.chSolverEden, "vWDB" );
	strcpy( conv.chSolverTemp, "vWDB" );
	strcpy( conv.chNotConverged, "none" );
	strcpy( conv.chConvEden, "none" );
	conv.resetConvIoniz();
	/* iterate to convergence flag */
	conv.lgAutoIt = false;
	/* convergence criteria */
	conv.autocv = 0.20f;
	conv.lgConvTemp = true;
	conv.lgConvPres = true;
	conv.lgConvEden = true;
	conv.lgUpdateCouplings = false;
	/* >>chng 04 jan 25, only set lgConvIoniz true where used in ConvXXX path */
	/*conv.lgConvIoniz = true;*/

	/* this option, use the new atmdat_rad_rec recombination rates */
	t_ADfA::Inst().set_version( PHFIT96 );

	/* age of the cloud, to check for time-steady */
	timesc.CloudAgeSet = -1.f;
	/* some timescale for CO and H2 */
	timesc.time_H2_Dest_longest = 0.;
	timesc.time_H2_Form_longest = 0.;
	/* remains neg if not evaluated */
	timesc.time_H2_Dest_here = -1.;
	timesc.time_H2_Form_here = 0.;

	timesc.BigCOMoleForm = 0.;

	timesc.TimeH21cm = 0.;
	timesc.sound_speed_isothermal = 0.;

	peimbt.tsqden = 1e7;

	/* CO related variables */
	co.codfrc = 0.;
	co.codtot = 0.;
	co.CODissHeat = 0.;

	NumDeriv.lgNumDeriv = false;

	/* index within the line in the line stack 
	 * default is Hbeta total - the third line in the stack
	 * 0th is a zero for sanity, 1st is unit, 2nd is a comment */
	/* >>chng 02 apr 22 from 2 to 3 since added unit at 1 */
	/* >>chng 06 mar 11, from 3 to -1 will now set to "H  1" 4861 */
	LineSave.ipNormWavL = -1;
	LineSave.WavLNorm = 4861.36f;
	LineSave.lgNormSet = false;
	LineSave.sig_figs = 4;

	/* the label for the normalization line */
	strcpy( LineSave.chNormLab, "    " );

	/* the scale factor for the normalization line */
	LineSave.ScaleNormLine = 1.;

	/* this is scale factor, reset with set resolution command, for setting
	 * the continuum resolution.  Setting to 0.1 will increase resolution by 10x.
	 * this multiplies the resolution contained in the continuum_mesh.ini file */
	continuum.ResolutionScaleFactor = 1.;

	continuum.lgCoStarInterpolationCaution = false;
	continuum.lgCon0 = false;

	/* upper limit to energies of inner shell opacities in Ryd
	 * this is 1 MeV by default */
	continuum.EnergyKshell = 7.35e4;

	/* free free heating, cooling, net */
	CoolHeavy.lgFreeOn = true;
	CoolHeavy.brems_cool_h = 0.;
	CoolHeavy.colmet = 0.;

	CoolHeavy.brems_cool_net = 0.;
	hydro.cintot = 0.;

	/* option to print emissivity instead of intensity/luminosity */
	hydro.lgHiPop2 = false;
	hydro.pop2mx = 0.;

	/* flag for Lya masing */
	hydro.HCollIonMax = 0.;

	/* type of hydrogen atom top off, options are " add" and "scal" 
	 * in versions 90 this was " add", but was "scal" in 91
	 * >>chng 99 jan 16, changed back to " add"*/
	/*strcpy( hydro.chHTopType, "scal" );*/
	strcpy( hydro.chHTopType, " add" );

	/* Lya excitation temperature, counter for hotter than gas */
	hydro.TexcLya = 0.;
	hydro.TLyaMax = 0.;
	hydro.nLyaHot = 0;

	/* option to kill damping wings of Lya */
	hydro.DampOnFac = 1.;

	/* is continuum pumping of H Lyman lines included?  yes, but turned off
	 * with atom h-like Lyman pumping off command */
	hydro.lgLymanPumping = true;

	/* multiplicative factor for all continuum pumping of H I Lyman lines,
	 * account for possible emission in the line */
	hydro.xLymanPumpingScaleFactor = 1.f;

	/* >>refer	abundance	D/H	Pettini, M., & Bowen, D.V., 2001, ApJ, 560, 41 */
	/* quoted error is +/- 0.35 */
	hydro.D2H_ratio = 1.65e-5;

	/* zero fractions of He0 destruction due to 23S */
	he.nzone = 0;
	he.frac_he0dest_23S = 0.;
	he.frac_he0dest_23S_photo = 0.;

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		/* option to disable continuum lowering */
		iso_ctrl.lgContinuumLoweringEnabled[ipISO] = true;

		/* flag set by compile he-like command, says to regenerate table of recombination coef */
		iso_ctrl.lgCompileRecomb[ipISO] = false;
		iso_ctrl.lgNoRecombInterp[ipISO] = false;

		/* how the gbar cs will be treated - set with atom he-like gbar command */
		/** \todo	2	change this to CS_new */
		iso_ctrl.lgCS_Vriens[ipISO] = true;
		iso_ctrl.lgCS_Vrinceanu[ipISO] = true;

		fixit(); /* make this the default for ipH_LIKE if not too slow.  */
		iso_ctrl.lgCS_Vrinceanu[ipH_LIKE] = false;

		iso_ctrl.lgCS_therm_ave[ipISO] = false;
		iso_ctrl.lgCS_None[ipISO] = false;
		/* when set try actually set to 1 or 2, depending on which fit is to be used,
		 * 1 is the broken power law fit */
		/* >>chng 02 dec 21, change to broken power law fit */
		iso_ctrl.nCS_new[ipISO] = 1;
		/* This flag says whether the density is high enough that helium is sufficiently l-mixed. */
		iso_ctrl.lgCritDensLMix[ipISO] = true;
		/* flag saying whether to include fine-structure mixing in spontaneous decays	
		 * set with ATOM HE-LIKE FSM command */
		iso_ctrl.lgFSM[ipISO] = 0;
		/* This is the flag saying whether to generate errors.  false means don't.	*/
		iso_ctrl.lgRandErrGen[ipISO] = false;
		/* this is the flag saying whether we should include excess recombination in the
		 * helike sequence.  Should only be off if testing effect of top off approximations. */
		iso_ctrl.lgTopoff[ipISO] = true;
		/* Dielectronic recombination for helike ions is on by default.	*/
		iso_ctrl.lgDielRecom[ipISO] = true;

		/* number of Lyman lines to include in opacities, this can be vastly larger
		 * than the number of actual levels in the model atom */
		iso_ctrl.nLyman[ipISO] = 100;
		iso_ctrl.nLyman_malloc[ipISO] = 100;

		/* controls whether l-mixing and collisional ionization included */
		iso_ctrl.lgColl_l_mixing[ipISO] = true;
		iso_ctrl.lgColl_excite[ipISO] = true;
		iso_ctrl.lgColl_ionize[ipISO] = true;
		iso_ctrl.lgLTE_levels[ipISO] = false;
		iso_ctrl.lgPrintNumberOfLevels = false;
	}

	/* Dielectronic recombination forming hydrogen-like ions does not exist. */
	iso_ctrl.lgDielRecom[ipH_LIKE] = false;

	/* smallest transition probability allowed */
	iso_ctrl.SmallA = 1e-30f;

	/* reset with SET IND2 command, turns on/off induced two photon */
	iso_ctrl.lgInd2nu_On = false;

	/* hydrogen redistribution functions */
	iso_ctrl.ipLyaRedist[ipH_LIKE] = ipPRD;
	iso_ctrl.ipResoRedist[ipH_LIKE] = ipCRD;
	iso_ctrl.ipSubRedist[ipH_LIKE] = ipCRDW;

	/* this is the upper level for each Lya, which uses the special ipLY_A */
	iso_ctrl.nLyaLevel[ipH_LIKE] = ipH2p;
	iso_ctrl.nLyaLevel[ipHE_LIKE] = ipHe2p1P;

	/* he-like redistribution functions */
	iso_ctrl.ipLyaRedist[ipHE_LIKE] = ipPRD;
	iso_ctrl.ipResoRedist[ipHE_LIKE] = ipCRD;
	iso_ctrl.ipSubRedist[ipHE_LIKE] = ipCRDW;

	iso_ctrl.lgPessimisticErrors = false;

	/* do not average collision strengths - evaluate at kT 
	 * set true with command SET COLLISION STRENGHTS AVERAGE */
	iso_ctrl.lgCollStrenThermAver = false;


	/**********************************************************************
	 * all parameters having to do with secondary ionization 
	 * by suprathermal electrons 
	 **********************************************************************/
	secondaries.SetCsupra = 0.;
	secondaries.lgCSetOn = false;
	secondaries.lgSecOFF = false;
	secondaries.SecHIonMax = 0.;

	secondaries.HeatEfficPrimary = 1.;
	secondaries.SecIon2PrimaryErg = 0.;
	secondaries.SecExcitLya2PrimaryErg = 0.;
	secondaries.x12tot = 0.;
	secondaries.sec2total = 0.;

	if( lgFirstCall )
	{
		/* malloc space for supra[nelem][ion] */
		secondaries.csupra = (realnum **)MALLOC( (unsigned)LIMELM*sizeof(realnum *) );
		secondaries.csupra_effic = (realnum **)MALLOC( (unsigned)LIMELM*sizeof(realnum *) );
		for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
		{
			secondaries.csupra[nelem] = (realnum *)MALLOC( (unsigned)(nelem+1)*sizeof(realnum) );
			secondaries.csupra_effic[nelem] = (realnum *)MALLOC( (unsigned)(nelem+1)*sizeof(realnum) );
		}
	}
	for( long nelem=ipHYDROGEN; nelem<LIMELM; ++nelem )
	{
		for( long ion=0; ion<nelem+1; ++ion )
		{
			/* secondary ionization rate for each species */
			secondaries.csupra[nelem][ion] = 0.;
			/* the rate of each species relative to H0 */
			secondaries.csupra_effic[nelem][ion] = 1.f;
		}
	}
	/* this scale factor is from table 10 of Tielens & Hollenbach 1985 */
	secondaries.csupra_effic[ipHELIUM][0] = 1.08f;

	/* on first call, these arrays do not exist, only zero here on 
	 * second and later calls, on first call, create them */
	if( lgFirstCall )
	{
		/* these will save bound electron recoil information data */
		ionbal.ipCompRecoil = 
			(long**)MALLOC(sizeof(long*)*(unsigned)LIMELM );
		ionbal.CompRecoilIonRate = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		ionbal.CompRecoilIonRateSave = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		ionbal.CompRecoilHeatRate = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		ionbal.CompRecoilHeatRateSave = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		ionbal.PhotoRate_Shell = 
			(double****)MALLOC(sizeof(double***)*(unsigned)LIMELM );
		ionbal.CollIonRate_Ground = 
			(double***)MALLOC(sizeof(double**)*(unsigned)LIMELM );
		ionbal.UTA_ionize_rate = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		ionbal.UTA_heat_rate = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );

		/* these are source and sink terms for heavy element ionization balance from the
		 * chemistry */
		mole.source = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		mole.sink = 
			(double**)MALLOC(sizeof(double*)*(unsigned)LIMELM );
		mole.xMoleChTrRate = 
			(realnum***)MALLOC(sizeof(realnum**)*(unsigned)LIMELM );

		/* space for ionization recombination arrays */
		ionbal.RateIoniz = (double ***)MALLOC(sizeof(double **)*(unsigned)LIMELM );
		ionbal.RateRecomTot = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		ionbal.RateRecomIso = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		ionbal.RR_rate_coef_used = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		ionbal.RR_Verner_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );

		/* rate coefficients [cm3 s-1] for Badnell DR recombination */
		ionbal.DR_Badnell_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		ionbal.RR_Badnell_rate_coef = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );
		ionbal.CX_recomb_rate_used = (double **)MALLOC(sizeof(double *)*(unsigned)LIMELM );

		/* create arrays for ions */
		for( long nelem=0; nelem<LIMELM; ++nelem )
		{
			ionbal.DR_Badnell_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.RR_Badnell_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.CX_recomb_rate_used[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );

			ionbal.RateIoniz[nelem] = (double **)MALLOC(sizeof(double *)*(unsigned)(nelem+1) );
			ionbal.RateRecomTot[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );					

			for( long ion=0; ion<nelem+1; ++ion )
			{
				ionbal.RateIoniz[nelem][ion] = (double *)MALLOC(sizeof(double )*(unsigned)(nelem+2) );
				for( long ion2=0; ion2<nelem+2; ++ion2 )
					ionbal.RateIoniz[nelem][ion][ion2] = 0.;
			}

			ionbal.RR_rate_coef_used[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.RR_Verner_rate_coef[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.UTA_ionize_rate[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.UTA_heat_rate[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.ipCompRecoil[nelem] = 
				(long*)MALLOC(sizeof(long)*(unsigned)(nelem+1) );
			ionbal.CompRecoilIonRate[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.CompRecoilIonRateSave[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.CompRecoilHeatRate[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.CompRecoilHeatRateSave[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+1) );
			ionbal.PhotoRate_Shell[nelem] = 
				(double***)MALLOC(sizeof(double**)*(unsigned)(nelem+1) );
			ionbal.CollIonRate_Ground[nelem] = 
				(double**)MALLOC(sizeof(double*)*(unsigned)(nelem+1) );
			/* chemistry source and sink terms for ionization ladders */
			mole.source[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+2) );
			mole.sink[nelem] = 
				(double*)MALLOC(sizeof(double)*(unsigned)(nelem+2) );
			mole.xMoleChTrRate[nelem] = 
				(realnum**)MALLOC(sizeof(realnum*)*(unsigned)(nelem+2) );
			for( long ion=0; ion<nelem+2; ++ion )
			{
				mole.xMoleChTrRate[nelem][ion] = 
					(realnum*)MALLOC(sizeof(realnum)*(unsigned)(nelem+2) );
			}
			ionbal.RateRecomIso[nelem] = (double *)MALLOC(sizeof(double)*(unsigned)(NISO) );
			for( long ipISO=0; ipISO<NISO; ++ipISO )
			{
				ionbal.RateRecomIso[nelem][ipISO] = 0.;	
			}

			for( long ion=0; ion<nelem+1; ++ion )
			{
				/* >>chng 03 aug 09, set these to impossible values */
				ionbal.RateRecomTot[nelem][ion] = -1.;
				ionbal.UTA_ionize_rate[nelem][ion] = -1.;
				ionbal.UTA_heat_rate[nelem][ion] = -1.;
				ionbal.ipCompRecoil[nelem][ion] = -99;
				ionbal.CompRecoilIonRate[nelem][ion] = -1.;
				ionbal.CompRecoilIonRateSave[nelem][ion] = -1.;
				ionbal.CompRecoilHeatRate[nelem][ion] = -1.;
				ionbal.CompRecoilHeatRateSave[nelem][ion] = -1.;

				/* finish mallocing space */
				ionbal.PhotoRate_Shell[nelem][ion] = 
					(double**)MALLOC(sizeof(double*)*(unsigned)NSHELLS );
				ionbal.CollIonRate_Ground[nelem][ion] = 
					(double*)MALLOC(sizeof(double)*(unsigned)2 );
				for( long ns=0; ns<NSHELLS; ++ns )
				{
					ionbal.PhotoRate_Shell[nelem][ion][ns] = 
						(double*)MALLOC(sizeof(double)*(unsigned)3 );
				}

				/* now set to impossible values */
				ionbal.ipCompRecoil[nelem][ion] = -100000;
				ionbal.DR_Badnell_rate_coef[nelem][ion] = 0.;
				ionbal.RR_Badnell_rate_coef[nelem][ion] = 0.;
			}

			set_NaN( ionbal.RR_rate_coef_used[nelem], nelem+1 );
			set_NaN( ionbal.RR_Verner_rate_coef[nelem], nelem+1 );
			set_NaN( ionbal.CX_recomb_rate_used[nelem], nelem+1 );
		}
	}

	/* now zero out these arrays */
	for( long nelem=0; nelem< LIMELM; ++nelem )
	{
		for( long ion=0; ion<nelem+1; ++ion )
		{

			ionbal.CompRecoilHeatRate[nelem][ion] = 0.;
			ionbal.CompRecoilIonRate[nelem][ion] = 0.;
			ionbal.UTA_ionize_rate[nelem][ion] = 0.;
			ionbal.UTA_heat_rate[nelem][ion] = 0.;
			ionbal.CollIonRate_Ground[nelem][ion][0] = 0.;
			ionbal.CollIonRate_Ground[nelem][ion][1] = 0.;
			ionbal.RateRecomTot[nelem][ion] = 0.;
			for( long ns=0; ns < NSHELLS; ++ns )
			{
				/* must be zero since ion routines use these when
				 * not yet defined */
				ionbal.PhotoRate_Shell[nelem][ion][ns][0] = 0.;
				ionbal.PhotoRate_Shell[nelem][ion][ns][1] = 0.;
				ionbal.PhotoRate_Shell[nelem][ion][ns][2] = 0.;
			}
		}
		/* these have one more ion than above */
		for( long ion=0; ion<nelem+2; ++ion )
		{
			/* zero out the source and sink arrays */
			mole.source[nelem][ion] = 0.;
			mole.sink[nelem][ion] = 0.;
			for( long ion2=0; ion2<nelem+2; ++ion2 )
			{
				mole.xMoleChTrRate[nelem][ion][ion2] = 0.;
			}
		}
	}

	ionbal.lgPhotoIoniz_On = true;
	ionbal.lgCompRecoil = true;

	/* these three adjust the treatment of UTA ionization */
	ionbal.lgInnerShellLine_on = true;
	ionbal.lgInnerShell_Kisielius = true;
	ionbal.lgInnerShell_Gu06 = true;

	/* default condition is burgess suppressed, Nussbaumer and Storey not */
	ionbal.lgSupDie[0] = true;
	ionbal.lgSupDie[1] = false;

	ionbal.lgNoCota = false;
	for( long nelem = 0; nelem < LIMELM; ++nelem )
	{
		ionbal.CotaRate[nelem] = 0.;
	}
	ionbal.ilt = 0;
	ionbal.iltln = 0;
	ionbal.ilthn = 0;
	ionbal.ihthn = 0;
	ionbal.ifail = 0;
	ionbal.lgGrainIonRecom = true;

	/* option to print recombination coefficient then exit */
	ionbal.lgRecom_Badnell_print = false;
	ionbal.guess_noise = 0.;

	/**********************************************************************
	 * these are options to print errors to special window,
	 * set with print errors command, 
	 * output will go to standard error
	 * defined in cdInit 
	 **********************************************************************/
	lgPrnErr = false;
	ioPrnErr = stderr;

	/* main arrays to save ionization fractions*/
	dense.zero();
	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		dense.SetGasPhaseDensity( nelem, 0. );
		for( long ion=0; ion < LIMELM+1; ion++ )
		{
			dense.xIonDense[nelem][ion] = 0.;
		}
	}
	dense.xMassTotal = 0.;

	/* this is the simple Fred Hamann FeII atom */
	t_fe2ovr_la::Inst().zero_opacity();

	/* zero out volume and column density save arrays */
	mean.MeanZero();

	/* zero out heating rates */
	HeatZero();

	/* some parameters dealing with calcium */
	ca.Ca2RmLya = 0.;
	ca.popca2ex = 0.;
	ca.Ca3d = 0.;
	ca.Ca4p = 0.;
	ca.dstCala = 0.;

	/* C12/C13 isotope ratio, sets the ratio of C12O16 to C13O16 and 
	 * C13/C12 1909 */
	// Negative says use default in external file, unless command appears.
	co.C12_C13_isotope_ratio = -1.;

	/* this is the default allowed relative error in the pressure */
	conv.PressureErrorAllowed = 0.01f;

	conv.MaxFractionalDensityStepPerIteration = 0.03;

	/* default error in total gas-phase density of each element, including molecules */
	conv.GasPhaseAbundErrorAllowed = 1e-5f;

	/* this is abort option set with SET PRESIONIZ command */
	conv.limPres2Ioniz = 3000;

	conv.nTeFail = 0;
	conv.nTotalFailures = 0;
	conv.nPreFail = 0;
	conv.failmx = 0.;
	conv.nIonFail = 0;
	conv.nPopFail = 0;
	conv.nNeFail = 0;
	conv.nGrainFail = 0;
	conv.dCmHdT = 0.;

	/* some titles and line images */
	for( long i=0; i<74; ++i)
	{
		input.chTitle[i] = ' ';
	}
	input.chTitle[75] = '\0';

	/* velocity field information */
	/* the turbulent velocity at illuminated face, internally in cm/s,
	 * but entered with turbulence command in km/s */
	DoppVel.TurbVel = 0.;
	/* is a turbulent gradient imposed?  Default is no. */
	DoppVel.lgTurbLawOn = false;
	/* the log of the turbulence gradient power law. Default is zero. */
	DoppVel.TurbVelLaw = 0.;
	/* the parameter F in eq 34 of
	 *>>refer	pressure	turb	Heiles, C. & Troland, T.H. 2005, 624, 773 */
	DoppVel.Heiles_Troland_F = 0.;
	/* is TurbVel included in pressure? - can be done two ways, with the velocity
	 * being set of with equipartition - true when TurbVel set if not equipartition,
	 * false with NO PRESSURE option on turbulence command */
	DoppVel.lgTurb_pressure = true;
	/* The scale in cm over which the turbulence is dissipated.  Normally 0,
	 * only set if dissipate keyword appears on turbulence command */
	DoppVel.DispScale = 0.;
	/* equipartition option on turbulence command, to set turbulence from B */
	DoppVel.lgTurbEquiMag = false;

	/* pressure related variables */

	pressure.RhoGravity_dark = 0.;
	pressure.RhoGravity_self = 0.;
	pressure.RhoGravity_external = 0.;
	pressure.RhoGravity = 0.;
	pressure.IntegRhoGravity = 0.;
	pressure.gravity_symmetry = -1;
	pressure.self_mass_factor = 1.;

	pressure.PresRamCurr = 0.;
	pressure.pres_radiation_lines_curr = 0.;
	pressure.lgPradCap = false;
	pressure.lgPradDen = false;
	pressure.lgLineRadPresOn = true;
	/* normally abort when radiation pressure exceeds gas pressure in const pres mod,
	 * this is option to keep going, set with NO ABORT on constant pressure command */
	pressure.lgRadPresAbortOK = true;
	/* Ditto for whether to stop at sonic point, this gets set to false
	 * for some of the dynamics pressure modes (strongd, shock, antishock)*/
	pressure.lgSonicPointAbortOK = true;
	/* this flag will say we hit the sonic point */
	pressure.lgSonicPoint = false;
	/* True when no physical solution for desired pressure in strong D fronts */
	pressure.lgStrongDLimbo = false; 

	pressure.RadBetaMax = 0.;
	pressure.ipPradMax_nzone = -1;
	pressure.PresMax = 0.;

	/* initial and current pressures */
	pressure.PresTotlInit = 0.;
	pressure.PresTotlCurr = 0.;

	/* zero out some dynamics variables */
	DynaZero();

	phycon.lgPhysOK = true;
	/* largest relative changes in Te, ne, H+, H2, and CO in structure
	 * this is computed as part of prtcomment so does not exist when code not talking,
	 * set to zero in zero and still zero if prtcomment not called */
	phycon.BigJumpTe = 0.;
	phycon.BigJumpne = 0.;
	phycon.BigJumpH2 = 0.;
	phycon.BigJumpCO = 0.;

	dense.xNucleiTotal = 1.;
	/* WJH */
	dense.xMassDensity0 = -1.0f; 

	// needed for TempChange to work but arrays needed for EdenChange to
	// work are not yet defined
	dense.eden = 1.;

	/* now set physical conditions array 
	 * following will force updating all temperature - density variables */
	TempChange( 1e4 , true);

	/* this is a scale factor that changes the n(H0)*1.7e-4 that is added to the
	* electron density to account for collisions with atomic H.  it is an order of
	* magnitude guess, so this command provides ability to see whether it affects results */
	dense.HCorrFac = 1.f;

	dark.lgNFW_Set = false;
	dark.r_200 = 0.;
	dark.r_s = 0.;

	atoms.nNegOI = 0;
	for( long i=0; i< N_OI_LEVELS; ++i )
	{
		atoms.popoi[i] = 0.;
	}
	atoms.popmg2 = 0.;

	/* do we want to save negative opacities */
	opac.lgNegOpacIO = false;

	opac.otsmin = 0.;

	/* this flag says to use the static opacities,
	 * only evaluate them at start of each new zone.
	 * when set false with 
	 * no static opacities
	 * command, always reevaluate them */
	opac.lgOpacStatic =  true;

	/* set true in radinc if negative opacities ever occur */
	opac.lgOpacNeg = false;

	/* can turn of scattering opacities for some geometries */
	opac.lgScatON = true;

	/* variables having to do with compiling and/or using the
	 * ancillary file of stored opacities */
	opac.lgCompileOpac = false;
	/* "no file opacity" command sets following var false, says not to use file */
	/** \todo	2	file opacities are disabled for now - reinstate this when arrays settle down */
	opac.lgUseFileOpac = false;

	dynamics.Upstream_density = 0.;
	/* effects of fast neutrons */
	hextra.frcneu = 0.;
	hextra.effneu = 1.;
	hextra.totneu = 0.;
	hextra.lgNeutrnHeatOn = false;
	hextra.CrsSecNeutron = 4e-26;

	opac.stimax[0] = 0.;
	opac.stimax[1] = 0.;


	hmi.H2_total = 0.;
	hmi.H2_total_f = 0.f;
	hmi.HD_total = 0.;
	hmi.H2_frac_abund_set = 0.;
	hmi.hmihet = 0.;
	hmi.h2plus_heat = 0.;
	hmi.HeatH2Dish_used = 0.;
	hmi.HeatH2Dexc_used = 0.;
	hmi.HeatH2Dish_TH85 = 0.;
	hmi.HeatH2Dexc_TH85 = 0.;
	hmi.UV_Cont_rel2_Draine_DB96_face = 0.;
	hmi.UV_Cont_rel2_Draine_DB96_depth = 0.;
	hmi.UV_Cont_rel2_Habing_TH85_face = 0.;
	hmi.UV_Cont_rel2_Habing_TH85_depth = 0.;
	hmi.HeatH2DexcMax = 0.;
	hmi.CoolH2DexcMax = 0.;
	hmi.hmitot = 0.;
	hmi.H2Opacity = 0.;
	hmi.hmidep = 1.;
	hmi.h2dep = 1.;
	hmi.h2pdep = 1.;
	hmi.h3pdep = 1.;

	/* option to kill effects of H2 in CO chemistry - set with
	 * set Leiden hack h2* off */
	hmi.lgLeiden_Keep_ipMH2s = true;
	hmi.lgLeidenCRHack = true;

	/* flag to turn off molecular network */
	mole_global.lgNoMole = false;
	mole_global.lgNoHeavyMole = false;
	/* capture of molecules onto grain surfaces - formation of ices
	 * flag says to include this process - turned off with the
	 * command NO GRAIN MOLECULES */
	mole_global.lgGrain_mole_deplete = true;
	/* flag saying that H2O water destruction rate went to zero */
	mole_global.lgH2Ozer = false;
	/* option to turn on the UMIST rates, naturally this will be 1, set to zero
	   with the set UMIST rates command */
	mole_global.lgLeidenHack = false;
	/* option to use diffuse cloud chemical rates from Table 8 of
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	mole_global.lgFederman = true;
	/* option to use effective temperature as defined in
	 * >> refer Zsargo, J. & Federman, S. R. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	mole_global.lgNonEquilChem = false;
	/** option to set proton elimination rates to zero
	 * >>refer	CO	chemistry	Huntress, W. T., 1977, ApJS, 33, 495
	 * By default, this is false - changed with set chemistry command */
	mole_global.lgProtElim = true;
	/** option to not include neutrals in the non-equilibrium scheme
	 * >> refer Federman, S. R. & Zsargo, J. 2003, ApJ, 589, 319
	 * By default, this is false - changed with set chemistry command */
	mole_global.lgNeutrals = true;
	/* option to use H2 continuum dissociation cross sections computed by P.C. Stancil
	 * By default, this is true - changed with "set H2 continuum dissociation xxx" command
	 * options are "Stancil" or "AD69" */
	mole_global.lgStancil = false;
	// all isotopes are currently disabled by default
	mole_global.lgTreatIsotopes.resize( LIMELM );
	fill( mole_global.lgTreatIsotopes.begin(), mole_global.lgTreatIsotopes.end(), false );

	/* this says which estimate of the rate of the Solomon process to use,
	 * default is Tielens & Hollenbach 1985a, changed with 
	 * set h2 Solomon command, options are TH85 and BD96,
	 * the second for the Bertoldi & Draine rates - they
	 * differ by 1 dex.   when large H2 turned on this is ignored */
	/* the Tielens & Hollenbach 1985 treatment */
	hmi.chH2_small_model_type = 'T';
	/* the improved H2 formalism given by 
	*>>refer	H2	dissoc	Burton, M.G., Hollenbach, D.J. & Tielens, A.G.G.M 
	>>refcon	1990, ApJ, 365, 620 */
	hmi.chH2_small_model_type = 'H';
	/* the Bertoldi & Draine 1996 treatment */
	/* >>chng 03 nov 15, change default to BD96 */
	hmi.chH2_small_model_type = 'B';
	/* >>chng 05 dec 08, use the Elwert et al. approximations as the default */
	hmi.chH2_small_model_type = 'E';

	/* set NaN */
	set_NaN( hmi.HeatH2Dish_used ); 
	set_NaN( hmi.HeatH2Dish_TH85 );
	set_NaN( hmi.HeatH2Dish_BD96 );
	set_NaN( hmi.HeatH2Dish_BHT90 );
	set_NaN( hmi.HeatH2Dish_ELWERT );

	/** HeatH2Dexc_used is heating due to collisional deexcitation of vib-excited 
	* H2 actually used */
	set_NaN( hmi.HeatH2Dexc_used );
	set_NaN( hmi.HeatH2Dexc_TH85 );
	set_NaN( hmi.HeatH2Dexc_BD96 );
	set_NaN( hmi.HeatH2Dexc_BHT90 );
	set_NaN( hmi.HeatH2Dexc_ELWERT );

	/** these are derivative wrt temp for collisional processes within X */
	set_NaN( hmi.deriv_HeatH2Dexc_used );
	set_NaN( hmi.deriv_HeatH2Dexc_TH85 );
	set_NaN( hmi.deriv_HeatH2Dexc_BD96 );
	set_NaN( hmi.deriv_HeatH2Dexc_BHT90 );
	set_NaN( hmi.deriv_HeatH2Dexc_ELWERT );

	set_NaN( hmi.H2_Solomon_dissoc_rate_used_H2g );
	set_NaN( hmi.H2_Solomon_dissoc_rate_TH85_H2g );
	set_NaN( hmi.H2_Solomon_dissoc_rate_BHT90_H2g );
	set_NaN( hmi.H2_Solomon_dissoc_rate_BD96_H2g );
	set_NaN( hmi.H2_Solomon_dissoc_rate_ELWERT_H2g );

	set_NaN( hmi.H2_Solomon_dissoc_rate_used_H2s );
	set_NaN( hmi.H2_Solomon_dissoc_rate_TH85_H2s );
	set_NaN( hmi.H2_Solomon_dissoc_rate_BHT90_H2s );
	set_NaN( hmi.H2_Solomon_dissoc_rate_BD96_H2s );
	set_NaN( hmi.H2_Solomon_dissoc_rate_ELWERT_H2s );

	/** the Solomon process rate H2 dissociates into X continuum - actually used */
	/**set_NaN( hmi.H2_Solomon_dissoc_rate_used );*/
	/** H2 + hnu => 2H from TH85 */
	/** H2 + hnu => 2H actually used */
	set_NaN( hmi.H2_photodissoc_used_H2g );
	set_NaN( hmi.H2_photodissoc_used_H2s );
	set_NaN( hmi.H2_photodissoc_ELWERT_H2g );
	set_NaN( hmi.H2_photodissoc_ELWERT_H2s );
	set_NaN( hmi.H2_photodissoc_TH85 );
	set_NaN( hmi.H2_photodissoc_BHT90 );

	/* default grain formation pumping - Takahashi 2001 */
	hmi.chGrainFormPump = 'T';

	/* set which approximation for Jura rate - Cazaux & Tielens 
	 * >>refer	H2	form	Cazaux, S., & Tielens, A.G.G.M., 2002, ApJ, 575, L29 */
	hmi.chJura = 'C';

	/* scale factor to multiply Jura rate, set Jura rate command */
	hmi.ScaleJura = 1.f;

	/* binding energy for change in H2 population while on grain surface,
	 * set with "set h2 Tad" command */
	hmi.Tad = 800.;

	hmi.lgH2_Thermal_BigH2 = true;
	hmi.lgH2_Chemistry_BigH2 = true;

	/* zero out some column densities */
	for( long i=0; i < NCOLD; i++ )
	{
		colden.colden[i] = 0.;
	}
	colden.He123S = 0.;
	colden.coldenH2_ov_vel = 0.;

	/* F=0 and F=1 column densities of H0*/
	colden.H0_21cm_upper =0;
	colden.H0_21cm_lower =0;

	for( long i=0; i < 5; i++ )
	{
		colden.C2Pops[i] = 0.;
		colden.C2Colden[i] = 0.;
		/* pops and column density for SiII atom */
		colden.Si2Pops[i] = 0.;
		colden.Si2Colden[i] = 0.;
	}
	for( long i=0; i < 3; i++ )
	{
		/* pops and column density for CI atom */
		colden.C1Pops[i] = 0.;
		colden.C1Colden[i] = 0.;
		/* pops and column density for OI atom */
		colden.O1Pops[i] = 0.;
		colden.O1Colden[i] = 0.;
		/* pops and column density for CIII atom */
		colden.C3Pops[i] = 0.;
	}
	for( long i=0; i < 4; i++ )
	{
		/* pops and column density for CIII atom */
		colden.C3Colden[i] = 0.;
	}

	/* variables to do with Jeans mass and radius */
	colden.TotMassColl = 0.;
	colden.tmas = 0.;
	colden.wmas = 0.;
	colden.rjnmin = FLT_MAX;
	colden.ajmmin = FLT_MAX;

	/* variables dealing with the radius */
	radius.rinner = 0.;
	radius.distance = 0.;
	radius.Radius = 0.;
	radius.Radius_mid_zone = 0.;
	radius.depth = DEPTH_OFFSET;
	radius.depth_mid_zone = DEPTH_OFFSET/2.;
	radius.depth_x_fillfac = 0.;
	radius.lgRadiusKnown = false;
	radius.drad = 0.;
	radius.drad_mid_zone = 0.;
	radius.r1r0sq = 1.;
	/* this is changed with the roberto command, to go from out to in */
	radius.dRadSign = 1.;

	/* RDFALT is log of default starting radius (cm) */
	/* >>chng 03 nov 12, from 25 to 30 for Lya clouds */
	/*radius.rdfalt = 25.;*/
	radius.rdfalt = 30.;

	/* set default cylinder thickness */
	radius.CylindHigh = 1e35f;
	radius.lgCylnOn = false;

	radius.drad_x_fillfac = 1.;
	radius.darea_x_fillfac = 1.;
	radius.dVeffVol = 1.;
	radius.dVeffAper = 1.;
	radius.drNext = 1.;
	radius.dRNeff = 1.;
	radius.lgdR2Small = false;

	radius.sdrmin = SMALLFLOAT;
	radius.lgSdrminRel = false;
	radius.sdrmax = 1e30;
	radius.lgSdrmaxRel = false;
	radius.lgSMinON = false;
	radius.lgDrMnOn = true;
	radius.lgFixed = false;
	radius.sdrmin_rel_depth = 1e-5;

	radius.lgDrMinUsed = false;

	rfield.lgHabing = false;

	/* flag to turn off Lya ots */
	rfield.lgLyaOTS = true;
	/* HeII rec and Lya ots */
	rfield.lgHeIIOTS = true;
	rfield.lgKillOTSLine = false;
	rfield.lgKillOutLine = false;
	rfield.lgKillOutCont = false;

	/* rfield.DiffPumpOn is unity unless process disabled by setting to 1
	 * with no diffuse line pumping command */
	rfield.DiffPumpOn = 1.;

	/* 02 jun 13, by Ryan...added this line	*/
	rfield.lgCompileGauntFF = false;

	/* >>chng 03 nov 28, add option to not do line transfer */
	rfield.lgDoLineTrans = true;

	/* flag saying whether to constantly reevaluated opacities -
	 * set false with no opacity reevaluate command */
	rfield.lgOpacityReevaluate = true;

	/* flag saying whether to constantly reevaluated ionization -
	 * set false with no ionization reevaluate command */
	rfield.lgIonizReevaluate = true;
	/* this element is default for choosing line width */
	rfield.fine_opac_nelem = ipIRON;
	/* there will be this many resolution elements in one FWHM for this element,
	 * at the lowest temperature to be considered */
	rfield.fine_opac_nresolv = 1;
	/* continuum scale factor for case of time varying continuum */
	rfield.time_continuum_scale = 1.;
	/* will fine optical depths be punched? */
	rfield.lgSaveOpacityFine = false;

	/* first is set true if one of the incident continua needs to have
	 * H-ionizing radiation blocked.  Second is set true is it is blocked
	 * with extinguish command - want both true if first is true */
	rfield.lgMustBlockHIon = false;
	rfield.lgBlockHIon = false;

	/* reset some variable related to cooling */
	CoolZero();

	thermal.lgCNegChk = true;
	thermal.CoolHeatMax = 0.;
	thermal.wlCoolHeatMax = 0;
	thermal.totcol = 0.;
	thermal.heatl = 0.;
	thermal.coolheat = 0.;
	thermal.lgCExtraOn = false;
	thermal.CoolExtra = 0.;
	thermal.ctot = 1.;

	thermal.htot = 1.;
	thermal.power = 0.;
	thermal.FreeFreeTotHeat = 0.;
	thermal.char_tran_cool = 0.;
	thermal.char_tran_heat = 0.;

	fnzone = 0.;
	nzone = 0;
	/* save initial condition for talk in case PRINT LAST used */
	called.lgTalkSave = called.lgTalk;

	oxy.poiii2 = 0.;
	oxy.poiii3 = 0.;
	oxy.poiexc = 0.;

	oxy.d5007r = 0.;
	oxy.d5007t = 0.;
	oxy.d4363 = 0.;
	oxy.d6300 = 0.;

	atmdat.nsbig = 0;

	/***************************************************
	 * charge transfer ionization and recombination 
	 ***************************************************/
	/* HCharHeatMax, HCharCoolMax are largest fractions of heating in cur zone
	 * or cooling due to ct */
	atmdat.HCharHeatMax = 0.;
	atmdat.HCharCoolMax = 0.;

	for ( long nelem=0; nelem < t_atmdat::NCX; ++nelem)
	{
		atmdat.CharExcIonTotal[nelem] = 0.;
		atmdat.CharExcRecTotal[nelem] = 0.;
	}
	atmdat.HIonFrac = 0.;
	atmdat.HIonFracMax = 0.;
	/* option to turn off all charge transfer, turned off with no charge transfer command */
	atmdat.lgCTOn = true;

	/* flag saying that charge transfer heating should be included,
	 * turned off with no CTHeat commmand */
	atmdat.HCharHeatOn = 1.;
	for( long nelem1=0; nelem1 < t_atmdat::NCX; ++nelem1)
	{
		for( long nelem=0; nelem< LIMELM; ++nelem )
		{
			for( long ion=0; ion<LIMELM; ++ion )
			{
				atmdat.CharExcIonOf[nelem1][nelem][ion] = 0.;
				atmdat.CharExcRecTo[nelem1][nelem][ion] = 0.;
			}
		}
	}

	/* >>chng 97 jan 6, from 0 to 8.5e-10*q as per Alex Dalgarno e-mail
	 * >>chng 97 feb 6, from 8.5e-10*q 1.92e-9x as per Alex Dalgarno e-mail */
	atmdat.HCTAlex = 1.92e-9;

	for( long nelem=0; nelem < LIMELM; nelem++ )
	{
		/* these are depletion scale factors */
		abund.depset[nelem] = 1.;
		/*begin sanity check */
		if( abund.depset[nelem] == 0. )
		{
			fprintf( ioQQQ, " ZERO finds insane abundance or depletion.\n" );
			fprintf( ioQQQ, " atomic number=%6ld abundance=%10.2e depletion=%10.2e\n", 
			  nelem, abund.solar[nelem], abund.depset[nelem] );
			ShowMe();
			cdEXIT(EXIT_FAILURE);
		}
		/*end sanity check */
	}

	/* typical ISM depletion factors, subjective mean of Cowie and Songaila
	 * and Jenkins 
	 * */
	abund.Depletion[0] = 1.;
	abund.Depletion[1] = 1.;
	abund.Depletion[2] = .16f;
	abund.Depletion[3] = .6f;
	abund.Depletion[4] = .13f;
	abund.Depletion[5] = 0.4f;
	abund.Depletion[6] = 1.0f;
	abund.Depletion[7] = 0.6f;
	abund.Depletion[8] = .3f;
	abund.Depletion[9] = 1.f;
	abund.Depletion[10] = 0.2f;
	abund.Depletion[11] = 0.2f;
	abund.Depletion[12] = 0.01f;
	abund.Depletion[13] = 0.03f;
	abund.Depletion[14] = .25f;
	abund.Depletion[15] = 1.0f;
	abund.Depletion[16] = 0.4f;
	abund.Depletion[17] = 1.0f;
	abund.Depletion[18] = .3f;
	abund.Depletion[19] = 1e-4f;
	abund.Depletion[20] = 5e-3f;
	abund.Depletion[21] = 8e-3f;
	abund.Depletion[22] = 6e-3f;
	abund.Depletion[23] = 6e-3f;
	abund.Depletion[24] = 5e-2f;
	abund.Depletion[25] = 0.01f;
	abund.Depletion[26] = 0.01f;
	abund.Depletion[27] = 0.01f;
	abund.Depletion[28] = .1f;
	abund.Depletion[29] = .25f;

	abund.lgDepln = false;
	abund.ScaleMetals = 1.;

	/* this tells the code to use standard Auger yields */
	t_yield::Inst().reset_yield();

	rt.dTauMase = 0.;
	rt.lgMaserCapHit = false;
	rt.lgMaserSetDR = false;

	rt.DoubleTau = 1.;
	rt.lgFstOn = true;
	rt.lgElecScatEscape = true;

	/* there was a call to TestCode */
	lgTestCodeCalled = false;
	/* test code enabled with set test command */
	lgTestCodeEnabled = false;

	/* zero out some grain variables */
	GrainZero();

	/* this is flag saying whether this is very first call,
	 * a time when space has not been allocated */
	lgFirstCall = false;
	return;
}
