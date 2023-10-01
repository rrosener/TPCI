/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*OpacityCreateAll compute initial set of opacities for all species */
/*OpacityCreate1Element generate ionic subshell opacities by calling t_ADfA::Inst().phfit */
/*opacity_more_memory allocate more memory for opacity stack */
/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
/*hmiopc derive total H- H minus opacity */
/*rayleh compute Rayleigh scattering cross section for Lya */
/*OpacityValenceRescale routine to rescale non-OP valence shell cross sections */
/******************************************************************************
 *NB NB NB  NB NB NB NB NB NB NB  NB NB NB NB
 * everything set here must be written to the opacity store files
 *
 ****************************************************************************** */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "continuum.h"
#include "iso.h"
#include "hydrogenic.h"
#include "oxy.h"
#include "trace.h"
#include "heavy.h"
#include "rfield.h"
#include "hmi.h"
#include "atmdat.h"
#include "save.h"
#include "grains.h"
#include "thirdparty.h"
#include "hydro_bauman.h"
#include "opacity.h"
#include "helike_recom.h"
#include "taulines.h"
#include "h2.h"
#include "h2_priv.h"
#include "ipoint.h"
#include "mole.h"

static const int NCSH2P = 10;

/* limit to number of opacity cells available in the opacity stack*/
static long int ndimOpacityStack = 2600000L;

/*OpacityCreate1Element generate opacities for entire element by calling t_ADfA::Inst().phfit */
STATIC void OpacityCreate1Element(long int nelem);

/*opacity_more_memory allocate more memory for opacity stack */
STATIC void opacity_more_memory(void);

/*hmiopc derive total H- H minus opacity */
STATIC double hmiopc(double freq);

/*rayleh compute Rayleigh scattering cross section for Lya */
STATIC double rayleh(double ener);

/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
STATIC double Opacity_iso_photo_cs( double energy , long ipISO , long nelem , long index );

/*OpacityCreateReilMan generate photoionization cross sections from Reilman and Manson points */
STATIC void OpacityCreateReilMan(long int low, 
  long int ihi, 
  const realnum cross[], 
  long int ncross, 
  long int *ipop, 
  const char *chLabl);

static bool lgRealloc = false;

/*OpacityCreatePowerLaw generate array of cross sections using a simple power law fit */
STATIC void OpacityCreatePowerLaw(
	/* lower energy limit on continuum mesh */
	long int ilo, 
	/* upper energy limit on continuum mesh */
	long int ihi, 
	/* threshold cross section */
	double cross, 
	/* power law index */
	double s, 
	/* pointer to opacity offset where this starts */
	long int *ip);

/*ofit compute cross sections for all shells of atomic oxygen */
STATIC double ofit(double e, 
	  realnum opart[]);

/*OpacityValenceRescale routine to rescale non-OP valence shell cross sections for atom */
STATIC void OpacityValenceRescale(
	/* element number on C scale */
	long int nelem ,
	/* scale factor, must be >= 0. */
	double scale )
{

	long int ion , nshell , low , ihi , ipop , ip;

	DEBUG_ENTRY( "OpacityValenceRescale()" );

	/* return if element is not turned on 
	 * >>chng 05 oct 19, this had not been done, so low in the opacity offset below was
	 * not set, and opacity index was negative - only problem when K turned off */
	if( !dense.lgElmtOn[nelem] )
	{
		return;
	}

	/* scale better be >= 0. */
	ASSERT( scale >= 0. );

	ion = 0;
	/* this is valence shell on C scale */
	nshell = Heavy.nsShells[nelem][ion] - 1;

	/* set lower and upper limits to this range */
	low = opac.ipElement[nelem][ion][nshell][0];
	ihi = opac.ipElement[nelem][ion][nshell][1];
	ipop = opac.ipElement[nelem][ion][nshell][2];

	/* loop over energy range of this shell */
	for( ip=low-1; ip < ihi; ip++ )
	{
		opac.OpacStack[ip-low+ipop] *= scale;
	}
	return;
}

void OpacityCreateAll(void)
{
	long int i, 
	  ipISO ,
	  need ,
	  nelem;

	realnum opart[7];

	double crs, 
	  dx,
	  eps, 
	  thom, 
	  thres, 
	  x;

	/* >>chng 02 may 29, change to lgOpacMalloced */
	/* remember whether opacities have ever been evaluated in this coreload
	static bool lgOpEval = false; */

	/* fits to cross section for photo dist of H_2^+ */
	static const realnum csh2p[NCSH2P]={6.75f,0.24f,8.68f,2.5f,10.54f,7.1f,12.46f,
	  6.0f,14.28f,2.7f};

	DEBUG_ENTRY( "OpacityCreateAll()" );

	/* H2+ h2plus h2+ */

	/* make and print dust opacities
	 * fill in dstab and dstsc, totals, zero if no dust
	 * may be different if different grains are turned on */
	GrainsInit();

	/* flag lgOpacMalloced says whether opacity stack has been generated
	 * only do this one time per core load  */
	/* >>chng 02 may 29, from lgOpEval to lgOpacMalloced */
	if( lgOpacMalloced )
	{
		/* this is not the first time code called */
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " OpacityCreateAll called but NOT evaluated since already done.\n" );
		}
		return;
	}

	lgOpacMalloced = true;

	/* create the space for the opacity stack */
	opac.OpacStack = (double*)MALLOC((size_t)ndimOpacityStack*sizeof(double));

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " OpacityCreateAll called, evaluating.\n" );
	}

	/* zero out opac since this array sometimes addressed before OpacityAddTotal called */
	for( i=0; i < rfield.nupper; i++ )
	{
		opac.opacity_abs[i] = 0.;
	}

	/* nOpacTot is number of opacity cells in OpacStack filled so far by opacity generating routines */
	opac.nOpacTot = 0;

	/* photoionization of h, he-like iso electronic sequences */
	for( ipISO=ipH_LIKE; ipISO<=ipHE_LIKE; ++ipISO )
	{
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] )
			{
				long int nupper;

				/* this is the opacity offset in the general purpose pointer array */
				/* indices are type, shell. ion, element, so this is the inner shell,
				 * NB - this works for H and He, but the second index should be 1 for Li */
				opac.ipElement[nelem][nelem-ipISO][0][2] = opac.nOpacTot + 1;
	
				fixit(); // opacities really need to be owned by states or species and stored in STL containers
					// so we don't have to mess with remembering what the upper and lower limits are

				// all iso states go to high-energy limit of code 
				nupper = rfield.nupper;
				for( long index=0; index < iso_sp[ipISO][nelem].numLevels_max; index++ )
				{
					/* this is array index to the opacity offset */
					iso_sp[ipISO][nelem].fb[index].ipOpac = opac.nOpacTot + 1;

					/* first make sure that first energy point is at least near the limit */
					/* >>chng 01 sep 23, increased factor form 0.98 to 0.94, needed since cells now go
					 * so far into radio, where resolution is poor */
					ASSERT( rfield.AnuOrg[iso_sp[ipISO][nelem].fb[index].ipIsoLevNIonCon-1] > 0.94f * 
						iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd );

					/* number of cells we will need to do this level */
					need = nupper - iso_sp[ipISO][nelem].fb[index].ipIsoLevNIonCon + 1;
					ASSERT( need > 0 );

					if( opac.nOpacTot + need > ndimOpacityStack )
						opacity_more_memory();

					for( i=iso_sp[ipISO][nelem].fb[index].ipIsoLevNIonCon-1; i < nupper; i++ )
					{
						double crs = 
							Opacity_iso_photo_cs( rfield.AnuOrg[i] , ipISO , nelem , index );
						opac.OpacStack[i-iso_sp[ipISO][nelem].fb[index].ipIsoLevNIonCon+iso_sp[ipISO][nelem].fb[index].ipOpac] = crs; 
						if( index==iso_sp[ipISO][nelem].numLevels_max-1 )
							iso_sp[ipISO][nelem].HighestLevelOpacStack.push_back( crs );
					}

					opac.nOpacTot += need;
				}
			}
		}
	}

	/* H2 continuum dissociation opacity */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( (*diatom)->lgEnabled && mole_global.lgStancil )
		{
			for( vector< diss_tran >::iterator tran = (*diatom)->Diss_Trans.begin(); tran != (*diatom)->Diss_Trans.end(); ++tran )
			{
				/* choose to integrate from 0.1 to 4 Ryd, data only extends from 0.7 to ~2 Ryd */
				long lower_limit = ipoint(tran->energies[0]);
 				long upper_limit = ipoint(tran->energies.back());
 				upper_limit = MIN2( upper_limit, rfield.nflux-1 );
 				long ipCont_Diss = opac.nOpacTot + 1;
 				long num_points = 0;
 
 				/* >>chng 02 may 08, by Ryan.  Added this and other checks for allotted memory.	*/
 				if( opac.nOpacTot + upper_limit - lower_limit + 1 > ndimOpacityStack )
 					opacity_more_memory();
 
				for(i = lower_limit; i <= upper_limit; ++i)
				{
					opac.OpacStack[ipCont_Diss + num_points - 1] = 
						MolDissocCrossSection( *tran, rfield.anu[i] );
					++num_points;
				}
				opac.nOpacTot += num_points;
			}
		}
	}
	
	/* This check will get us through Klein-Nishina below.	*/
	if( opac.nOpacTot + iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon + rfield.nupper > ndimOpacityStack )
		opacity_more_memory();

	/* Lyman alpha damping wings - Rayleigh scattering */
	opac.ipRayScat = opac.nOpacTot + 1;
	for( i=0; i < iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon; i++ )
	{
		opac.OpacStack[i-1+opac.ipRayScat] = rayleh(rfield.AnuOrg[i]);
	}
	opac.nOpacTot += iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon;

	/* ==============================================================
	 * this block of code defines the electron scattering cross section
	 * for all energies */

	/* assume Thomson scattering up to ipCKshell, 20.6 Ryd=0.3 keV */
	thom = 6.65e-25;
	opac.iopcom = opac.nOpacTot + 1;
	for( i=0; i < opac.ipCKshell; i++ )
	{
		opac.OpacStack[i-1+opac.iopcom] = thom;
		/*fprintf(ioQQQ,"%.3e\t%.3e\n", 
			rfield.AnuOrg[i]*EVRYD , opac.OpacStack[i-1+opac.iopcom] );*/
	}

	/* Klein-Nishina from eqn 7.5, 
	 * >>refer	Klein-Nishina	cs	Rybicki and Lightman */
	for( i=opac.ipCKshell; i < rfield.nupper; i++ )
	{
		dx = rfield.AnuOrg[i]/3.7573e4;

		opac.OpacStack[i-1+opac.iopcom] = thom*3.e0/4.e0*((1.e0 + 
		  dx)/POW3(dx)*(2.e0*dx*(1.e0 + dx)/(1.e0 + 2.e0*dx) - log(1.e0+
		  2.e0*dx)) + 1.e0/2.e0/dx*log(1.e0+2.e0*dx) - (1.e0 + 3.e0*
		  dx)/POW3(1.e0 + 2.e0*dx));
		/*fprintf(ioQQQ,"%.3e\t%.3e\n", 
			rfield.AnuOrg[i]*EVRYD , opac.OpacStack[i-1+opac.iopcom] );*/
	}
	opac.nOpacTot += rfield.nupper - 1 + 1;

	/* ============================================================== */

	/* This check will get us through "H- hminus H minus bound-free opacity" below.	*/
	/* >>chng 02 may 08, by Ryan.  Added this and other checks for allotted memory.	*/
	if( opac.nOpacTot + 3*rfield.nupper - opac.ippr + iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon - hmi.iphmin + 2 > ndimOpacityStack )
		opacity_more_memory();

	/* pair production */
	opac.ioppr = opac.nOpacTot + 1;
	for( i=opac.ippr-1; i < rfield.nupper; i++ )
	{
		/* pair production heating rate for unscreened H + He
		 * fit to figure 41 of Experimental Nuclear Physics,
		 * Vol1, E.Segre, ed */

		x = rfield.AnuOrg[i]/7.512e4*2.;

		opac.OpacStack[i-opac.ippr+opac.ioppr] = 5.793e-28*
		  POW2((-0.46737118 + x*(0.349255416 + x*0.002179893))/(1. + 
		  x*(0.130471301 + x*0.000524906)));
		/*fprintf(ioQQQ,"%.3e\t%.3e\n", 
			rfield.AnuOrg[i]*EVRYD , opac.OpacStack[i-opac.ippr+opac.ioppr] );*/
	}
	opac.nOpacTot += rfield.nupper - opac.ippr + 1;

	/* brems (free-free) opacity */
	opac.ipBrems = opac.nOpacTot + 1;

	for( i=0; i < rfield.nupper; i++ )
	{
		/* missing factor of 1E-20 to avoid underflow
		 * free free opacity needs g(ff)*(1-exp(hn/kT))/SQRT(T)*1E-20 */
		opac.OpacStack[i-1+opac.ipBrems] = 
			/*(realnum)(1.03680e-18/POW3(rfield.AnuOrg[i]));*/
			/* >>chng 00 jun 05, descale by 1e10 so that underflow at high-energy
			 * end does not occur */
			1.03680e-8/POW3(rfield.AnuOrg[i]);
	}
	opac.nOpacTot += rfield.nupper - 1 + 1;

	opac.iphmra = opac.nOpacTot + 1;
	for( i=0; i < rfield.nupper; i++ )
	{
		/* following is ratio of h minus to neut h bremss opacity */
		opac.OpacStack[i-1+opac.iphmra] = 0.1175*rfield.anusqr[i];
	}
	opac.nOpacTot += rfield.nupper - 1 + 1;

	opac.iphmop = opac.nOpacTot + 1;
	for( i=hmi.iphmin-1; i < iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon; i++ )
	{
		/* H- hminus H minus bound-free opacity */
		opac.OpacStack[i-hmi.iphmin+opac.iphmop] = 
			hmiopc(rfield.AnuOrg[i]);
	}
	opac.nOpacTot += iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon - hmi.iphmin + 1;

	/* ============================================================== */

	/* This check will get us through "H2 photoionization cross section" below.	*/
	/* >>chng 07 oct 10, by Ryan.  Added this check for allotted memory.	*/
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
	{
		if( opac.nOpacTot + rfield.nupper - (*diatom)->ip_photo_opac_thresh + 1 > ndimOpacityStack )
			opacity_more_memory();

		(*diatom)->ip_photo_opac_offset = opac.nOpacTot + 1;
		opac.nOpacTot += (*diatom)->OpacityCreate( opac.OpacStack );
	}

	/* H2+ H2P h2plus photoabsorption
	 * cs from 
	 * >>refer	H2+	photodissoc	Buckingham, R.A., Reid, S., & Spence, R. 1952, MNRAS 112, 382, 0 K temp */
	OpacityCreateReilMan(opac.ih2pnt[0],opac.ih2pnt[1],csh2p,NCSH2P,&opac.ih2pof,
	  "H2+ ");

	/* This check will get us through "HeI singlets neutral helium ground" below.	*/
	/* >>chng 02 may 08, by Ryan.  Added this and other checks for allotted memory.	*/
	if( opac.nOpacTot + rfield.nupper - iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon + 1 > ndimOpacityStack )
		opacity_more_memory();

	/* HeI singlets neutral helium ground */
	opac.iophe1[0] = opac.nOpacTot + 1;
	opac.ipElement[ipHELIUM][0][0][2] = opac.iophe1[0];
	for( i=iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon-1; i < rfield.nupper; i++ )
	{
		crs = t_ADfA::Inst().phfit(2,2,1,rfield.AnuOrg[i]*EVRYD);
		opac.OpacStack[i-iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon+opac.iophe1[0]] = 
			crs*1e-18;
	}
	opac.nOpacTot += rfield.nupper - iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ipIsoLevNIonCon + 1;

	/* these are opacity offset points that would be defined in OpacityCreate1Element,
	 * but this routine will not be called for H and He
	 * generate all heavy element opacities, everything heavier than He,
	 * nelem is on the C scale, so Li is 2 */
	/*>>chng 99 jan 27, do not reevaluate hydrogenic opacity below */
	for( nelem=2; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			OpacityCreate1Element(nelem);
		}
	}

	/* option to rescale atoms of some elements that were not done by opacity project
	 * the valence shell - two arguments - element number on C scale, and scale factor */
	/*>>chng 05 sep 26, fudge factor to get atomic K fraction along well defined line of sight
	 * to be observed value - this is ratio of cross sections, actual value is very uncertain since
	 * differences betweeen Verner & opacity project are huge */
	OpacityValenceRescale( ipPOTASSIUM , 5. );

	/* now add on some special cases, where exicted states, etc, come in */
	/* Nitrogen
	 * >>refer	n1	photo	Henry, R., ApJ 161, 1153.
	 * photoionization of excited level of N+ */
	OpacityCreatePowerLaw(opac.in1[0],opac.in1[1],9e-18,1.75,&opac.in1[2]);

	/* atomic Oxygen
	 * only do this if 1996 Verner results are used */
	if( dense.lgElmtOn[ipOXYGEN] && t_ADfA::Inst().get_version() == PHFIT96 )
	{
		/* This check will get us through this loop.	*/
		/* >>chng 02 may 08, by Ryan.  Added this and other checks for allotted memory.	*/
		if( opac.nOpacTot + opac.ipElement[ipOXYGEN][0][2][1] - opac.ipElement[ipOXYGEN][0][2][0] + 1 > ndimOpacityStack )
			opacity_more_memory();

		/* integrate over energy range of the valence shell of atomic oxygen*/
		for( i=opac.ipElement[ipOXYGEN][0][2][0]-1; i < opac.ipElement[ipOXYGEN][0][2][1]; i++ )
		{
			/* call special routine to evaluate partial cross section for OI shells */
			eps = rfield.AnuOrg[i]*EVRYD;
			crs = ofit(eps,opart);

			/* this will be total cs of all processes leaving shell 3 */
			crs = opart[0];
			for( long n=1; n < 6; n++ )
			{
				/* add up table of cross sections */
				crs += opart[n];
			}
			/* convert to cgs and overwrite cross sections set by OpacityCreate1Element */
			crs *= 1e-18;
			opac.OpacStack[i-opac.ipElement[ipOXYGEN][0][2][0]+opac.ipElement[ipOXYGEN][0][2][2]] = crs;
		}
		/* >>chng 02 may 09 - this was a significant error */
		/* >>chng 02 may 08, by Ryan.  This loop did not update total slots filled.	*/
		opac.nOpacTot += opac.ipElement[ipOXYGEN][0][2][1] - opac.ipElement[ipOXYGEN][0][2][0] + 1;
	}

	/* Henry nubmers for 1S excit state of OI, OP data very sparse */
	OpacityCreatePowerLaw(opac.ipo1exc[0],opac.ipo1exc[1],4.64e-18,0.,&opac.ipo1exc[2]);

	/* photoionization of excited level of O2+ 1D making 5007
	 * fit to TopBase Opacity Project cs */
	OpacityCreatePowerLaw(opac.ipo3exc[0],opac.ipo3exc[1],3.8e-18,0.,&opac.ipo3exc[2]);

	/* photoionization of excited level of O2+ 1S making 4363 */
	OpacityCreatePowerLaw(opac.ipo3exc3[0],opac.ipo3exc3[1],5.5e-18,0.01,
	  &opac.ipo3exc3[2]);

	/* This check will get us through the next two steps.	*/
	/* >>chng 02 may 08, by Ryan.  Added this and other checks for allotted memory.	*/
	if( opac.nOpacTot + iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon - oxy.i2d + 1 
		+ iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - opac.ipmgex + 1 > ndimOpacityStack )
		opacity_more_memory();

	/* photoionization to excited states of O+ */
	opac.iopo2d = opac.nOpacTot + 1;
	thres = rfield.AnuOrg[oxy.i2d-1];
	for( i=oxy.i2d-1; i < iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon; i++ )
	{
		crs = 3.85e-18*(4.4*pow(rfield.AnuOrg[i]/thres,-1.5) - 3.38*
		  pow(rfield.AnuOrg[i]/thres,-2.5));

		opac.OpacStack[i-oxy.i2d+opac.iopo2d] = crs;
	}
	opac.nOpacTot += iso_sp[ipH_LIKE][ipHELIUM].fb[ipH1s].ipIsoLevNIonCon - oxy.i2d + 1;

	/* magnesium
	 * photoionization of excited level of Mg+
	 * fit to opacity project data Dima got */
	opac.ipOpMgEx = opac.nOpacTot + 1;
	for( i=opac.ipmgex-1; i < iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon; i++ )
	{
		opac.OpacStack[i-opac.ipmgex+opac.ipOpMgEx] = 
			(0.2602325880970085 + 
		  445.8558249365131*exp(-rfield.AnuOrg[i]/0.1009243952792674))*
		  1e-18;
	}
	opac.nOpacTot += iso_sp[ipH_LIKE][ipHYDROGEN].fb[ipH1s].ipIsoLevNIonCon - opac.ipmgex + 1;

	ASSERT( opac.nOpacTot < ndimOpacityStack );

	/* Calcium
	 * excited states of Ca+ */
	OpacityCreatePowerLaw(opac.ica2ex[0],opac.ica2ex[1],4e-18,1.,&opac.ica2op);

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, 
			" OpacityCreateAll return OK, number of opacity cells used in OPSC= %ld and OPSV is dimensioned %ld\n", 
		  opac.nOpacTot, ndimOpacityStack );
	}

	/* option to compile opacities into file for later use 
	 * this is executed if the 'compile opacities' command is entered */
	if( opac.lgCompileOpac )
	{
		fprintf( ioQQQ, "The COMPILE OPACITIES command is currently not supported\n" );
		cdEXIT(EXIT_FAILURE);
	}

	if( lgRealloc )
		fprintf(ioQQQ,
		" Please consider revising ndimOpacityStack in OpacityCreateAll, a total of %li cells were needed.\n\n" , opac.nOpacTot);
	return;
}
/*OpacityCreatePowerLaw generate array of cross sections using a simple power law fit */
STATIC void OpacityCreatePowerLaw(
	/* lower energy limit on continuum mesh */
	long int ilo, 
	/* upper energy limit on continuum mesh */
	long int ihi, 
	/* threshold cross section */
	double cross, 
	/* power law index */
	double s, 
	/* pointer to opacity offset where this starts */
	long int *ip)
{
	long int i;
	double thres;

	DEBUG_ENTRY( "OpacityCreatePowerLaw()" );

	/* non-positive cross section is unphysical */
	ASSERT( cross > 0. );

	/* place in the opacity stack where we will stuff cross sections */
	*ip = opac.nOpacTot + 1;
	ASSERT( *ip > 0 );
	ASSERT( ilo > 0 );
	thres = rfield.anu[ilo-1];

	if( opac.nOpacTot + ihi - ilo + 1 > ndimOpacityStack )
		opacity_more_memory();

	for( i=ilo-1; i < ihi; i++ )
	{
		opac.OpacStack[i-ilo+*ip] = cross*pow(rfield.anu[i]/thres,-s);
	}

	opac.nOpacTot += ihi - ilo + 1;
	return;
}

/*OpacityCreateReilMan generate photoionization cross sections from Reilman and Manson points */
STATIC void OpacityCreateReilMan(long int low, 
  long int ihi, 
  const realnum cross[], 
  long int ncross, 
  long int *ipop, 
  const char *chLabl)
{
	long int i, 
	  ics, 
	  j, 
	  ncr;

	const int NOP = 100;
	realnum cs[NOP], 
	  en[NOP], 
	  slope;

	DEBUG_ENTRY( "OpacityCreateReilMan()" );

	/* this is the opacity entering routine designed for
	 * the Reilman and Manson tables.  It works with incident
	 * photon energy (entered in eV) and cross sections in megabarns
	 * */
	*ipop = opac.nOpacTot + 1;
	ASSERT( *ipop > 0 );

	ncr = ncross/2;
	if( ncr > NOP )
	{
		fprintf( ioQQQ, " Too many opacities were entered into OpacityCreateReilMan.  Increase the value of NOP.\n" );
		fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
		cdEXIT(EXIT_FAILURE);
	}

	/* the array CROSS has ordered pairs of elements.
	 * the first is the energy in eV (not Ryd)
	 * and the second is the cross section in megabarns */
	for( i=0; i < ncr; i++ )
	{
		en[i] = cross[i*2]/13.6f;
		cs[i] = cross[(i+1)*2-1]*1e-18f;
	}

	ASSERT( low>0 );
	if( en[0] > rfield.anu[low-1] )
	{
		fprintf( ioQQQ, 
			" OpacityCreateReilMan: The entered opacity energy bandwidth is not large enough (low fail).\n" );
		fprintf( ioQQQ, 
			" The desired energy (Ryd) was%12.5eeV and the lowest entered in the array was%12.5e eV\n", 
		  rfield.anu[low-1]*EVRYD, en[0]*EVRYD );

		fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
		fprintf( ioQQQ, " The original energy (eV) and cross section (mb) arrays follow:\n" );
		fprintf( ioQQQ, " " );

		for( i=0; i < ncross; i++ )
		{
			fprintf( ioQQQ, "%11.4e", cross[i] );
		}

		fprintf( ioQQQ, "\n" );
		cdEXIT(EXIT_FAILURE);
	}

	slope = (cs[1] - cs[0])/(en[1] - en[0]);
	ics = 1;

	if( opac.nOpacTot + ihi - low + 1 > ndimOpacityStack )
	  opacity_more_memory();

	/* now fill in the opacities using linear interpolation */
	for( i=low-1; i < ihi; i++ )
	{
		if( rfield.anu[i] > en[ics-1] && rfield.anu[i] <= en[ics] )
		{
			opac.OpacStack[i-low+*ipop] = cs[ics-1] + slope*(rfield.anu[i] - 
			  en[ics-1]);
		}

		else
		{
			ics += 1;
			if( ics + 1 > ncr )
			{
				fprintf( ioQQQ, " OpacityCreateReilMan: The entered opacity energy bandwidth is not large enough (high fail).\n" );
				fprintf( ioQQQ, " The entered energy was %10.2eeV and the highest in the array was %10.2eeV\n", 
				  rfield.anu[i]*13.6, en[ncr-1]*13.6 );
				fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl
				   );
				fprintf( ioQQQ, " The lowest energy enterd in the array was%10.2e eV\n", 
				  en[0]*13.65 );
				fprintf( ioQQQ, " The highest energy ever needed would be%10.2eeV\n", 
				  rfield.anu[ihi-1]*13.6 );
				fprintf( ioQQQ, " The lowest energy needed was%10.2eeV\n", 
				  rfield.anu[low-1]*13.6 );
				cdEXIT(EXIT_FAILURE);
			}

			slope = (cs[ics] - cs[ics-1])/(en[ics] - en[ics-1]);
			if( rfield.anu[i] > en[ics-1] && rfield.anu[i] <= en[ics] )
			{
				opac.OpacStack[i-low+*ipop] = cs[ics-1] + slope*(rfield.anu[i] - 
				  en[ics-1]);
			}
			else
			{
				ASSERT( i > 0);
				fprintf( ioQQQ, " Internal logical error in OpacityCreateReilMan.\n" );
				fprintf( ioQQQ, " The desired energy (%10.2eeV), I=%5ld, is not within the next energy bound%10.2e%10.2e\n", 
				  rfield.anu[i]*13.6, i, en[ics-1], en[ics] );

				fprintf( ioQQQ, " The previous energy (eV) was%10.2e\n", 
				  rfield.anu[i-1]*13.6 );

				fprintf( ioQQQ, " Here comes the energy array.  ICS=%4ld\n", 
				  ics );

				for( j=0; j < ncr; j++ )
				{
					fprintf( ioQQQ, "%10.2e", en[j] );
				}
				fprintf( ioQQQ, "\n" );

				fprintf( ioQQQ, " chLabl was %4.4s\n", chLabl );
				cdEXIT(EXIT_FAILURE);
			}
		}
	}
	/* >>chng 02 may 09, this was a significant logcal error */
	/* >>chng 02 may 08, by Ryan.  This routine did not update the total slots filled.	*/
	opac.nOpacTot += ihi - low + 1;
	return;
}


/*ofit compute cross sections for all shells of atomic oxygen */
STATIC double ofit(double e, 
	  realnum opart[])
{
	double otot,
	  q, 
	  x;

	static const double y[7][5] = {
		{8.915,3995.,3.242,10.44,0.0},
		{11.31,1498.,5.27,7.319,0.0},
		{10.5,1.059e05,1.263,13.04,0.0},
		{19.49,48.47,8.806,5.983,0.0},
		{50.,4.244e04,0.1913,7.012,4.454e-02},
		{110.5,0.1588,148.3,-3.38,3.589e-02},
		{177.4,32.37,381.2,1.083,0.0}
	};
	static const double eth[7]={13.62,16.94,18.79,28.48,50.,110.5,538.};
	static const long l[7]={1,1,1,0,1,1,0};

	DEBUG_ENTRY( "ofit()" );
	/*compute cross sections for all shells of atomic oxygen
	 * Photoionization of OI
	 * Input parameter:   e - photon energy, eV
	 * Output parameters: otot - total photoionization cross section, Mb
	 *  opart(1) - 2p-shell photoionization, goes to 4So
	 *  opart(2) - 2p-shell photoionization, goes to 2Do
	 *  opart(3) - 2p-shell photoionization, goes to 2Po
	 *  opart(4) - 2s-shell photoionization
	 *  opart(5) - double photoionization, goes to O++
	 *  opart(6) - triple photoionization, goes to O+++
	 *  opart(7) - 1s-shell photoionization */

	otot = 0.0;

	for( int i=0; i < 7; i++ )
	{
		opart[i] = 0.0;
	}

	for( int i=0; i < 7; i++ )
	{
		if( e >= eth[i] )
		{
			// this assert is trivially true, but it helps PGCC
			ASSERT( i < 7 );
			q = 5.5 - 0.5*y[i][3] + l[i];

			x = e/y[i][0];

			opart[i] = (realnum)(y[i][1]*(POW2(x - 1.0) + POW2(y[i][4]))/
			  pow(x,q)/pow(1.0 + sqrt(x/y[i][2]),y[i][3]));

			otot += opart[i];
		}
	}
	return(otot);
}

/******************************************************************************/

/*OpacityCreate1Element generate ionic subshell opacities by calling t_ADfA::Inst().phfit */
STATIC void OpacityCreate1Element(
		  /* atomic number on the C scale, lowest ever called will be Li=2 */
		  long int nelem)
{
	long int ihi, 
	  ip, 
	  ipop, 
	  limit, 
	  low, 
	  need, 
	  nelec, 
	  ion, 
	  nshell;
	double cs; 
	double energy;

	DEBUG_ENTRY( "OpacityCreate1Element()" );

	/* confirm range of validity of atomic number, Li=2 should be the lightest */
	ASSERT( nelem >= 2 );
	ASSERT( nelem < LIMELM );

	/*>>chng 99 jan 27, no longer redo hydrogenic opacity here */
	/* for( ion=0; ion <= nelem; ion++ )*/
	for( ion=0; ion < nelem; ion++ )
	{

		/* will be used for a sanity check on number of hits in a cell*/
		for( ip=0; ip < rfield.nupper; ip++ )
		{
			opac.opacity_abs[ip] = 0.;
		}

		/* number of bound electrons */
		nelec = nelem+1 - ion;

		/* loop over all shells, from innermost K shell to valence */
		for( nshell=0; nshell < Heavy.nsShells[nelem][ion]; nshell++ )
		{
			/* this is array index for start of this shell within large opacity stack */
			opac.ipElement[nelem][ion][nshell][2] = opac.nOpacTot +  1;

			/* this is continuum upper limit to array index for energy range of this shell */
			limit = opac.ipElement[nelem][ion][nshell][1];

			/* this is number of cells in continuum needed to store opacity */
			need = limit - opac.ipElement[nelem][ion][nshell][0] + 1;

			/* check that opac will have enough frequeny cells */
			if( opac.nOpacTot + need > ndimOpacityStack )
				opacity_more_memory();

			/* set lower and upper limits to this range */
			low = opac.ipElement[nelem][ion][nshell][0];
			ihi = opac.ipElement[nelem][ion][nshell][1];
			ipop = opac.ipElement[nelem][ion][nshell][2];

			/* make sure indices are within correct bounds,
			 * mainly check on logic for detecting missing shells */
			ASSERT( low <= ihi || low<5 );

			/* loop over energy range of this shell */
			for( ip=low-1; ip < ihi; ip++ )
			{
				/* photo energy MAX so that we never eval below threshold */
				energy = MAX2(rfield.AnuOrg[ip]*EVRYD , 
					t_ADfA::Inst().ph1(nshell,nelec-1,nelem,0));

				/* the cross section in mega barns */
				cs = t_ADfA::Inst().phfit(nelem+1,nelec,nshell+1,energy);
				/* cannot assert that cs is positive since, at edge of shell,
				 * energy might be slightly below threshold and hence zero,
				 * due to finite size of continuum bins */
				opac.OpacStack[ip-low+ipop] = cs*1e-18;

				/* add this to total opacity, which we will confirm to be greater than zero below */
				opac.opacity_abs[ip] += cs;
			}

			opac.nOpacTot += ihi - low + 1;

			/* save pointers option */
			if( save.lgPunPoint )
			{
				fprintf( save.ipPoint, "%3ld%3ld%3ld%10.2e%10.2e%10.2e%10.2e\n", 
				  nelem, ion, nshell, rfield.anu[low-1], rfield.anu[ihi-1], 
				  opac.OpacStack[ipop-1], opac.OpacStack[ihi-low+ipop-1] );
			}
		}

		ASSERT( Heavy.nsShells[nelem][ion] >= 1 );
		/*confirm that total opacity is greater than zero  */
		for( ip=opac.ipElement[nelem][ion][Heavy.nsShells[nelem][ion]-1][0]-1; 
			ip < continuum.KshellLimit; ip++ )
		{
			ASSERT( opac.opacity_abs[ip] > 0. );
		}

	}
	return;
}

/*opacity_more_memory allocate more memory for opacity stack */
STATIC void opacity_more_memory(void)
{

	DEBUG_ENTRY( "opacity_more_memory()" );

	/* double size */
	ndimOpacityStack *= 2;
	opac.OpacStack = (double *)REALLOC(  opac.OpacStack , (size_t)ndimOpacityStack*sizeof(double) );
	fprintf( ioQQQ, " NOTE OpacityCreate1Element needed more opacity cells than ndimOpacityStack,  please consider increasing it.\n" );
	fprintf( ioQQQ, " NOTE OpacityCreate1Element doubled memory allocation to %li.\n",ndimOpacityStack );
	lgRealloc = true;
	return;
}

/*Opacity_iso_photo_cs returns photoionization cross section for isoelectronic sequences */
STATIC double Opacity_iso_photo_cs( 
		/* photon energy ryd */
		double EgammaRyd , 
		/* iso sequence */
		long ipISO , 
		/* charge, 0 for H */
		long nelem , 
		/* index */
		long index )
{
	double crs=-DBL_MAX;

	DEBUG_ENTRY( "Opacity_iso_photo_cs()" );

	if( ipISO==ipH_LIKE )
	{
		if( index==0 )
		{
			/* this is the ground state, use Dima's routine, which works in eV
			 * and returns megabarns */
			double EgammaEV = MAX2(EgammaRyd*(realnum)EVRYD , t_ADfA::Inst().ph1(0,0,nelem,0));
			crs = t_ADfA::Inst().phfit(nelem+1,1,1,EgammaEV)* 1e-18;
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else if( index < iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max )
		{
			/* photon energy relative to threshold */
			double photon = MAX2( EgammaRyd/iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, 1. + FLT_EPSILON*2. );

			crs = H_photo_cs( photon , N_(index), L_(index), nelem+1 );
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else if( N_(index) <= NHYDRO_MAX_LEVEL )
		{
			/* for first cell, depending on the current resolution of the energy mesh,
			 * the center of the first cell can be below the ionization limit of the
			 * level.  do not let the energy fall below this limit */
			/* This will make sure that we don't call epsilon below threshold,
			 * the factor 1.001 was chosen so that t_ADfA::Inst().hpfit, which works
			 * in terms of Dima's Rydberg constant, is not tripped below threshold */
			EgammaRyd = MAX2( EgammaRyd , iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd*1.001f );

			crs = t_ADfA::Inst().hpfit(nelem+1,N_(index),EgammaRyd*EVRYD);
			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
		else
		{
			/* photon energy relative to threshold */
			double photon = MAX2( EgammaRyd/iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, 1. + FLT_EPSILON*2. );

			/* cross section for collapsed level should be 
			 * roughly equal to cross-section for yrast level,
			 * so third parameter is n - 1. */
			crs = H_photo_cs( photon , N_(index), N_(index)-1, nelem+1 );

			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
	}
	else if( ipISO==ipHE_LIKE )
	{
		EgammaRyd = MAX2( EgammaRyd , iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd);
		/* this would be a collapsed level */
		if( index >= iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max )
		{
			long int nup = iso_sp[ipHE_LIKE][nelem].n_HighestResolved_max + index + 1 -
				(iso_sp[ipHE_LIKE][nelem].numLevels_max - iso_sp[ipHE_LIKE][nelem].nCollapsed_max);

			/* this is a collapsed level - this is hydrogenic routine and
			 * first he-like energy may not agree exactly with threshold for H */
			crs = t_ADfA::Inst().hpfit(nelem,nup ,EgammaRyd*EVRYD);
			/* make sure cross section is reasonable if away from threshold */
			ASSERT( 
				(EgammaRyd < iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd*1.02) ||
				(crs > 0. && crs < 1e-10) );
		}
		else
		{
			long n = N_(index);
			long l = L_(index);
			long S = S_(index);
			/* He_cross_section returns cross section (cm^-2), 
			 * given EgammaRyd, the photon energy in Ryd,
			 * quantum numbers n, l, and S,
			 * nelem is charge, equal to 1 for Helium,
			 * this is a wrapper for cross_section */
			crs = He_cross_section( EgammaRyd, iso_sp[ipISO][nelem].fb[index].xIsoLevNIonRyd, n, l, S, nelem );

			/* make sure cross section is reasonable */
			ASSERT( crs > 0. && crs < 1e-10 );
		}
	}
	else
		TotalInsanity();
	return(crs);
}

/*hmiopc derive total H- H minus opacity */
static const int NCRS = 33;

STATIC double hmiopc(double freq)
{
	double energy, 
	  hmiopc_v, 
	  x, 
	  y;
	static double y2[NCRS];
	static double crs[NCRS]={0.,0.124,0.398,0.708,1.054,1.437,1.805,
	  2.176,2.518,2.842,3.126,3.377,3.580,3.741,3.851,3.913,3.925,
	  3.887,3.805,3.676,3.511,3.306,3.071,2.810,2.523,2.219,1.898,
	  1.567,1.233,.912,.629,.39,.19};
	static double ener[NCRS]={0.,0.001459,0.003296,0.005256,0.007351,
	  0.009595,0.01201,0.01460,0.01741,0.02044,0.02375,0.02735,0.03129,
	  0.03563,0.04043,0.04576,0.05171,0.05841,0.06601,0.07469,0.08470,
	  0.09638,0.1102,0.1268,0.1470,0.1723,0.2049,0.2483,0.3090,0.4001,
	  0.5520,0.8557,1.7669};
	static bool lgFirst = true;

	DEBUG_ENTRY( "hmiopc()" );

	/* bound free cross section (x10**-17 cm^2) from Doughty et al
	 * 1966, MNRAS 132, 255; good agreement with Wishart MNRAS 187, 59p. */

	/* photoelectron energy, add 0.05552 to get incoming energy (Ryd) */


	if( lgFirst )
	{
		/* set up coefficients for spline */
		spline(ener,crs,NCRS,2e31,2e31,y2);
		lgFirst = false;
	}

	energy = freq - 0.05552;
	if( energy < ener[0] || energy > ener[NCRS-1] )
	{
		hmiopc_v = 0.;
	}
	else
	{
		x = energy;
		splint(ener,crs,y2,NCRS,x,&y);
		hmiopc_v = y*1e-17;
	}
	return( hmiopc_v );
}

/*rayleh compute Rayleigh scattering cross section for Lya */
STATIC double rayleh(double ener)
{
	double rayleh_v;

	DEBUG_ENTRY( "rayleh()" );

	/** \todo	2	update to astro-ph/0308073, Lee, H-W, ApJ in press */
	/* do hydrogen Rayleigh scattering cross sections;
	 * fits to 
	 *>>refer	Ly	scattering	Gavrila, M., 1967, Physical Review 163, 147
	 * and Mihalas radiative damping
	 *
	 * >>chng 96 aug 15, changed logic to do more terms for each part of
	 * rayleigh scattering
	 * if( ener.lt.0.05 ) then
	 *  rayleh = 8.41e-25 * ener**4 * DampOnFac
	 * */
	if( ener < 0.05 )
	{
		rayleh_v = (8.41e-25*powi(ener,4) + 3.37e-24*powi(ener,6))*
		  hydro.DampOnFac;
	}

	else if( ener < 0.646 )
	{
		rayleh_v = (8.41e-25*powi(ener,4) + 3.37e-24*powi(ener,6) + 
		  4.71e-22*powi(ener,14))*hydro.DampOnFac;
	}

	else if( ener >= 0.646 && ener < 1.0 )
	{
		rayleh_v = fabs(0.74959-ener);
		rayleh_v = 1.788e5/POW2(FR1RYD*MAX2(0.001,rayleh_v));
		/*  typical energy between Ly-a and Ly-beta */
		rayleh_v = MAX2(rayleh_v,1e-24)*hydro.DampOnFac;
	}

	else
	{
		rayleh_v = 0.;
	}
	return( rayleh_v );
}
