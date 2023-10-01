/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SanityCheck, check that various parts of the code still work, called by Cloudy after continuum
 * and optical depth arrays are set up, but before initial temperature and ionization */
#include "cddefines.h"
#include "physconst.h"
#include "thirdparty.h"
#include "dense.h"
#include "elementnames.h"
#include "continuum.h"
#include "helike_recom.h"
#include "rfield.h"
#include "taulines.h"
#include "hypho.h"
#include "iso.h"
#include "opacity.h"
#include "hydro_bauman.h"
#include "hydrogenic.h"
#include "heavy.h"
#include "trace.h"
#include "cloudy.h"

namespace {
	class my_Integrand: public std::unary_function<double, double>
	{
	public:
		double operator() (double x)
		{
			return sin(x);
		}
		// the constructor is needed to avoid warnings by the Sun Studio compiler
		my_Integrand() {}
	};

	// create a version of sin with C++ linkage
	// this is done because we need a C++ pointer to this routine below
	double mySin(double x)
	{
		return sin(x);
	}
}

/* NB - this routine must not change any global variables - any that are changed as part of
 * a test must be reset, so that the code retains state */

STATIC void SanityCheckBegin(void );
STATIC void SanityCheckFinal(void );
/* chJob is either "begin" or "final" 
 * "begin is before code starts up
 * "final is after model is complete */
void SanityCheck( const char *chJob )
{
	DEBUG_ENTRY( "SanityCheck()" );

	if( strcmp(chJob,"begin") == 0 )
	{
		SanityCheckBegin();
	}
	else if( strcmp(chJob,"final") == 0 )
	{
		SanityCheckFinal();
	}
	else
	{
		fprintf(ioQQQ,"SanityCheck called with insane argument.\n");
		cdEXIT(EXIT_FAILURE);
	}
}

STATIC void SanityCheckFinal(void )
{
	/* PrtComment also has some ending checks on sanity */
}

STATIC void SanityCheckBegin(void )
{
	bool lgOK=true;
	int lgFlag;// error return for spsort, 0 success, >=1 for errors
	int32 ner, ipiv[3];
	long	i , 
		j , 
		nelem , 
		ion ,
		nshells;
	double *A;

	/* this will be charge to the 4th power */
	double Aul ,
		error,
		Z4, gaunt;

	long n, logu, loggamma2;

	const int NDIM = 10;
	double x , ans1 , ans2  , xMatrix[NDIM][NDIM] , yVector[NDIM] ,
		rcond;
	realnum *fvector;
	long int *ipvector;

	DEBUG_ENTRY( "SanityCheck()" );

	/*********************************************************
	 *                                                       *
	 * confirm that various part of cloudy still work        *
	 *                                                       *
	 *********************************************************/

	/* if this is no longer true at end, we have a problem */
	lgOK = true;

	/*********************************************************
	 *                                                       *
	 * check that all the Lyas As are ok                     *
	 *                                                       *
	 *********************************************************/
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		/* this element may be turned off */
		if( dense.lgElmtOn[nelem] )
		{ 
			/* H_Einstein_A( n, l, np, lp, iz ) - all are on physics scale */
			Aul = H_Einstein_A( 2, 1, 1, 0, nelem+1 );
			/*fprintf(ioQQQ,"%li\t%.4e\n", nelem+1, iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Aul() );*/
			if( fabs(Aul - iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Aul() ) /Aul > 0.01 )
			{
				fprintf(ioQQQ," SanityCheck found insane H-like As.\n");
				lgOK = false;
			}
		}
	}

	/*********************************************************
	 *                                                       *
	 * check that gaunt factors are good                     *
	 *                                                       *
	 *********************************************************/
	/* Uncommenting each of the four print statements here
	 * will produce a nice table comparable to Sutherland 98, Table 2.	*/
	/* fprintf(ioQQQ,"u\t-4\t-3\t-2\t-1\t0\t1\t2\t3\t4\n");*/
	for( logu=-4; logu<=4; logu++)
	{
		/*fprintf(ioQQQ,"%li\t", logu);*/
		for(loggamma2=-4; loggamma2<=4; loggamma2++)
		{ 
			double SutherlandGff[9][9]=
			{	{5.5243, 5.5213, 5.4983, 5.3780, 5.0090, 4.4354, 3.8317, 3.2472, 2.7008},
				{4.2581, 4.2577, 4.2403, 4.1307, 3.7816, 3.2436, 2.7008, 2.2126, 1.8041},
				{3.0048, 3.0125, 3.0152, 2.9434, 2.6560, 2.2131, 1.8071, 1.4933, 1.2771},
				{1.8153, 1.8367, 1.8880, 1.9243, 1.7825, 1.5088, 1.2886, 1.1507, 1.0747},
				{0.8531, 0.8815, 0.9698, 1.1699, 1.2939, 1.1988, 1.1033, 1.0501, 1.0237},
				{0.3101, 0.3283, 0.3900, 0.5894, 0.9725, 1.1284, 1.0825, 1.0419, 1.0202},
				{0.1007, 0.1080, 0.1335, 0.2281, 0.5171, 0.9561, 1.1065, 1.0693, 1.0355},
				{0.0320, 0.0344, 0.0432, 0.0772, 0.1997, 0.5146, 0.9548, 1.1042, 1.0680},
				{0.0101, 0.0109, 0.0138, 0.0249, 0.0675, 0.1987, 0.5146, 0.9547, 1.1040}};

			gaunt = cont_gaunt_calc( TE1RYD/pow(10.,(double)loggamma2), 1., pow(10.,(double)(logu-loggamma2)) );
			error = fabs( gaunt - SutherlandGff[logu+4][loggamma2+4] ) /gaunt;
			/*fprintf(ioQQQ,"%1.3f\t", gaunt);*/
			if( error>0.11 || ( loggamma2<2 && error>0.05 ) )
			{
				fprintf(ioQQQ," SanityCheck found insane gff. log(u) %li, log(gamma2) %li, error %.3e\n",
					logu, loggamma2, error );
				lgOK = false;
			}
		}
		/*fprintf(ioQQQ,"\n");*/
	}

	/*********************************************************
	 *                                                       *
	 * check some transition probabililties for he-like ions *
	 *                                                       *
	 *********************************************************/
	for( nelem=1; nelem<LIMELM; ++nelem )
	{
		/* the helike 9-1 transition, A(3^3P to 2^3S) */
		double as[]={
		 /* updated with Johnson values */
		 0.       ,9.47e+006 ,3.44e+008 ,1.74e+009 ,5.51e+009 ,1.34e+010 ,
		2.79e+010 ,5.32E+010 ,8.81e+010 ,1.46E+011 ,2.15e+011 ,3.15e+011 ,
		4.46e+011 ,6.39E+011 ,8.26e+011 ,1.09e+012 ,1.41e+012 ,1.86E+012 ,
		2.26e+012 ,2.80e+012 ,3.44e+012 ,4.18e+012 ,5.04e+012 ,6.02e+012 ,
		7.14e+012 ,8.40e+012 ,9.83e+012 ,1.14e+013 ,1.32e+013 ,1.52e+013
		};

		if( iso_sp[ipHE_LIKE][nelem].numLevels_max > 8 && dense.lgElmtOn[nelem])
		{	
			/* following used to print current values of As */
			if( fabs( as[nelem] - iso_sp[ipHE_LIKE][nelem].trans(9,1).Emis().Aul() ) /as[nelem] > 0.025 )
			{
				fprintf(ioQQQ,
					" SanityCheck found insane He-like As: expected, nelem=%li found: %.2e %.2e.\n",
					nelem,
					as[nelem] , 
					iso_sp[ipHE_LIKE][nelem].trans(9,1).Emis().Aul() );
				lgOK = false;
			}
		}
	}

	/* only do this test if case b is not in effect */
	if( !opac.lgCaseB )
	{

		for( i = 0; i <=110; i++ )
		{
			double DrakeTotalAuls[111] = {
				-1.0000E+00, -1.0000E+00, -1.0000E+00, 1.02160E+07,
				1.02160E+07, 1.02160E+07, 1.80090E+09, 2.78530E+07,
				1.82990E+07, 1.05480E+07, 7.07210E+07, 6.37210E+07,
				5.79960E+08, 1.60330E+07, 1.13640E+07, 7.21900E+06,
				3.11920E+07, 2.69830E+07, 1.38380E+07, 1.38330E+07,
				2.52270E+08, 9.20720E+06, 6.82220E+06, 4.56010E+06,
				1.64120E+07, 1.39290E+07, 7.16030E+06, 7.15560E+06,
				4.25840E+06, 4.25830E+06, 1.31150E+08, 5.62960E+06,
				4.29430E+06, 2.95570E+06, 9.66980E+06, 8.12340E+06,
				4.19010E+06, 4.18650E+06, 2.48120E+06, 2.48120E+06,
				1.64590E+06, 1.64590E+06, 7.65750E+07, 3.65330E+06,
				2.84420E+06, 1.99470E+06, 6.16640E+06, 5.14950E+06,
				2.66460E+06, 2.66200E+06, 1.57560E+06, 1.57560E+06,
				1.04170E+06, 1.04170E+06, 7.41210E+05, 7.41210E+05,
				4.84990E+07, 2.49130E+06, 1.96890E+06, 1.39900E+06,
				4.16900E+06, 3.46850E+06, 1.79980E+06, 1.79790E+06,
				1.06410E+06, 1.06410E+06, 7.02480E+05, 7.02480E+05,
				4.98460E+05, 4.98460E+05, -1.0000E+00, -1.0000E+00,
				3.26190E+07, 1.76920E+06, 1.41440E+06, 1.01460E+06,
				2.94830E+06, 2.44680E+06, 1.27280E+06, 1.27140E+06,
				7.52800E+05, 7.52790E+05, 4.96740E+05, 4.96740E+05,
				3.51970E+05, 3.51970E+05, -1.0000E+00, -1.0000E+00,
				-1.0000E+00, -1.0000E+00, 2.29740E+07, 1.29900E+06,
				1.04800E+06, 7.57160E+05, 2.16090E+06, 1.79030E+06,
				9.33210E+05, 9.32120E+05, 5.52310E+05, 5.52310E+05,
				3.64460E+05, 3.64460E+05, 2.58070E+05, 2.58070E+05,
				-1.0000E+00, -1.0000E+00, -1.0000E+00, -1.0000E+00,
				-1.0000E+00, -1.0000E+00, 1.67840E+07};

			if( DrakeTotalAuls[i] > 0. && 
				i < iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max - iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max)
			{
				if( fabs( DrakeTotalAuls[i] - (1./iso_sp[ipHE_LIKE][ipHELIUM].st[i].lifetime()) ) /DrakeTotalAuls[i] > 0.001 )
				{
					fprintf(ioQQQ,
						" SanityCheck found helium lifetime outside 0.1 pct of Drake values: index, expected, found: %li %.4e %.4e\n",
						i,
						DrakeTotalAuls[i], 
						(1./iso_sp[ipHE_LIKE][ipHELIUM].st[i].lifetime()) );
					lgOK = false;
				}
			}
		}
	}

	/*********************************************************
	 *                                                       *
	 * check the threshold photoionization cs for He I       *
	 *                                                       *
	 *********************************************************/
	if( dense.lgElmtOn[ipHELIUM] )
	{
		/* HeI photoionization cross sections at threshold for lowest 20 levels */
		const int NHE1CS = 20;
		double he1cs[NHE1CS] = 
		{
			5.480E-18 , 9.253E-18 , 1.598E-17 , 1.598E-17 , 1.598E-17 , 1.348E-17 , 
			8.025E-18 , 1.449E-17 , 2.852E-17 , 1.848E-17 , 1.813E-17 , 2.699E-17 , 
			1.077E-17 , 2.038E-17 , 4.159E-17 , 3.670E-17 , 3.575E-17 , 1.900E-17 , 
			1.900E-17 , 4.175E-17 
		};

		/* loop over levels and check on photo cross section */
		j = MIN2( NHE1CS+1 , iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max -iso_sp[ipHE_LIKE][ipHELIUM].nCollapsed_max );
		for( n=1; n<j; ++n )
		{
			/* above list of levels does not include the ground */
			i = iso_sp[ipHE_LIKE][ipHELIUM].fb[n].ipOpac;
			ASSERT( i>0 );
			/*fprintf(ioQQQ,"%li\t%lin", n , i );*/
			/* >>chng 02 apr 10, from 0.01 to 0.02, values stored
			 * where taken from calc at low contin resolution, when continuum
			 * resolution changed this changes too */
			/*fprintf(ioQQQ,"%li %.2e\n", n,( he1cs[n-1] - opac.OpacStack[i - 1] ) /he1cs[n-1] );*/
			/* >>chng 02 jul 16, limt from 0.02 to 0.04, so that "set resolution 4" will work */
			/* >>chng 04 may 18, levels 10 and 11 are about 12% off - because of energy binning, chng from 0.08 to 0.15 */ 
			if( fabs( he1cs[n-1] - opac.OpacStack[i - 1] ) /he1cs[n-1] > 0.15 )
			{
				fprintf(ioQQQ,
					" SanityCheck found insane HeI photo cs: expected, n=%li found: %.3e %.3e.\n",
					n,
					he1cs[n-1] , 
					opac.OpacStack[i - 1]);
				fprintf(ioQQQ,
					" n=%li, l=%li, s=%li\n",
					iso_sp[ipHE_LIKE][ipHELIUM].st[n].n() ,
					iso_sp[ipHE_LIKE][ipHELIUM].st[n].l() ,
					iso_sp[ipHE_LIKE][ipHELIUM].st[n].S());
				lgOK = false;
			}
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		long nelem = ipISO;
		/* Check for agreement between on-the-fly and interpolation calculations
		 * of recombination, but only if interpolation is turned on. */
		if( !iso_ctrl.lgNoRecombInterp[ipISO] )
		{
			/* check the recombination coefficients for ground state */
			error = fabs( iso_recomb_check( ipISO, nelem , 0 , 7500. ) );
			if( error > 0.01 )
			{
				fprintf(ioQQQ,
					" SanityCheck found insane1 %s %s recom coef: expected, n=%i error: %.2e \n",
					iso_ctrl.chISO[ipISO],
					elementnames.chElementSym[nelem],
					0,
					error );
				lgOK = false;
			}

			/* check the recombination coefficients for ground state of the root of each iso sequence */
			error = fabs( iso_recomb_check( ipISO, nelem , 1 , 12500. ) );
			if( error > 0.01 )
			{
				fprintf(ioQQQ,
					" SanityCheck found insane2 %s %s recom coef: expected, n=%i error: %.2e \n",
					iso_ctrl.chISO[ipISO],
					elementnames.chElementSym[nelem],
					1,
					error );
				lgOK = false;
			}
		}
	}

	/*********************************************************
	 *                                                       *
	 * check out the sorting routine                         *
	 *                                                       *
	 *********************************************************/

	const int NSORT = 100 ;

	fvector = (realnum *)MALLOC((NSORT)*sizeof(realnum) );

	ipvector = (long *)MALLOC((NSORT)*sizeof(long int) );

	nelem = 1;
	/* make up some unsorted values */
	for( i=0; i<NSORT; ++i )
	{
		nelem *= -1;
		fvector[i] = (realnum)nelem * ((realnum)NSORT-i);
	}

	/*spsort netlib routine to sort array returning sorted indices */
	spsort(fvector, 
		   NSORT, 
		  ipvector, 
		  /* flag saying what to do - 1 sorts into increasing order, not changing
		   * the original routine */
		  1, 
		  &lgFlag);

	if( lgFlag ) lgOK = false;

	for( i=1; i<NSORT; ++i )
	{
		/*fprintf(ioQQQ," %li %li %.0f\n", 
			i, ipvector[i],fvector[ipvector[i]] );*/
		if( fvector[ipvector[i]] <= fvector[ipvector[i-1]] )
		{
			fprintf(ioQQQ," SanityCheck found insane sort\n");
			lgOK = false;
		}
	}

	free( fvector );
	free( ipvector);

#	if 0
	ttemp = (realnum)sqrt(phycon.te);
	/* check that the temperatures make sense */
	if( fabs(ttemp - phycon.sqrte )/ttemp > 1e-5 )
	{
		fprintf(ioQQQ , "SanityCheck finds insane te %e sqrt te %e sqrte %e dif %e\n",
			phycon.te , 
			sqrt(phycon.te) , 
			phycon.sqrte , 
			fabs(sqrt(phycon.te) - phycon.sqrte ) );
		lgOK = false;
	}
#	endif

	/*********************************************************
	 *                                                       *
	 * confirm that widflx and anu arrays correspond         *
	 * to one another                                        *
	 *                                                       *
	 *********************************************************/

#	if 0
	/* this check on widflx can't be used since some sharpling curved continua, like laser,
	 * totally fail due to non-linear nature of widflx and anu relationship */
#	if !defined(NDEBUG)
	x = 0.;
	for( i=1; i<rfield.nupper-1; ++i )
	{
		if( fabs( ((rfield.anu[i+1]-rfield.anu[i]) + (rfield.anu[i]-rfield.anu[i-1])) /rfield.widflx[i] /2.-1.) > 0.02 )
		{
			ans1 = fabs( ((rfield.anu[i+1]-rfield.anu[i]) + (rfield.anu[i]-rfield.anu[i-1])) /rfield.widflx[i] /2.-1.);
			fprintf(ioQQQ," SanityCheck found insane widflx anu[i+1]=%e anu[i]=%e widflx=%e delta=%e rel err %e\n",
			rfield.anu[i+1] , rfield.anu[i] , rfield.widflx[i] , rfield.anu[i+1] -rfield.anu[i] , ans1 );
			lgOK = false;
			x = MAX2( ans1 , x);
		}
		/* problems when at energy where resolution of grid changes dramatically */
		/* this is resolution at current energy */
		ans1 = rfield.widflx[i] / rfield.anu[i];
		if( (rfield.anu[i]+rfield.widflx[i]/2.)*(1.-ans1/10.) > rfield.anu[i+1] - rfield.widflx[i+1]/2.) 
		{
			fprintf(ioQQQ," SanityCheck found insane overlap1 widflx %e %e %e %e %e %e\n",
			rfield.anu[i] , rfield.widflx[i], rfield.anu[i] + rfield.widflx[i]/2. , rfield.anu[i+1], 
			rfield.widflx[i+1], rfield.anu[i+1] -rfield.widflx[i+1]/2. );
			lgOK = false;
		}
		if( !lgOK )
		{
			fprintf(ioQQQ," big error was %e\n", x);
		}
	}
#	endif
#	endif


	/*********************************************************
	 *                                                       *
	 * confirm that hydrogen einstein As are still valid     *
	 *                                                       *
	 *********************************************************/
	for( nelem=0; nelem<2; ++nelem )
	{
		/* this element may be turned off */
		if( dense.lgElmtOn[nelem] )
		{ 
			/*Z4 = (double)(POW2(nelem+1)*POW2(nelem+1));*/
			/* form charge to the 4th power */
			Z4 = (double)(nelem+1);
			Z4 *= Z4;
			Z4 *= Z4;
			/* H Lya */
			ans1 = iso_sp[ipH_LIKE][nelem].trans(ipH2p,ipH1s).Emis().Aul();
			ans2 = 6.265e8*Z4;
			if( fabs(ans1-ans2)/ans2 > 1e-3 )
			{
				fprintf(ioQQQ , "SanityCheck finds insane A for H Lya %g %g nelem=%li\n",
					ans1 , ans2 , nelem );
				lgOK = false;
			}
		}
	}

	/* check that hydrogenic branching ratios add up to unity */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			int ipHi, ipLo;
			for( ipHi=4; ipHi< iso_sp[ipH_LIKE][nelem].numLevels_max-iso_sp[ipH_LIKE][nelem].nCollapsed_max; ++ipHi )
			{
				double sum = 0.;
				for( ipLo=0; ipLo<ipHi; ++ipLo )
				{
					sum += iso_sp[ipH_LIKE][nelem].BranchRatio[ipHi][ipLo];
				}
				if( fabs(sum-1.)>0.01 ) 
				{
					fprintf(ioQQQ , 
						"SanityCheck H branching ratio sum not unity for nelem=%li upper level=%i sum=%.3e\n",
						nelem, ipHi, sum );
					lgOK = false;
				}
			}
		}
	}

	/* check that hydrogenic lifetimes are sane (compare inverse sum of As with closed form of lifetime) */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		if( dense.lgElmtOn[nelem] )
		{
			int ipHi, ipLo;
			for( ipHi=1; ipHi< iso_sp[ipH_LIKE][nelem].numLevels_max-iso_sp[ipH_LIKE][nelem].nCollapsed_max; ++ipHi )
			{
				double inverse_sum = 0.;
				double sum = 0.;
				long ipISO = ipH_LIKE;
				
				/* we do not have an accurate closed form for l=0 lifetimes.  Everything else should be very accurate. */
				if( L_(ipHi)==0 )
					continue;

				for( ipLo=0; ipLo<ipHi; ++ipLo )
				{
					sum += iso_sp[ipH_LIKE][nelem].trans(ipHi,ipLo).Emis().Aul();
				}

				inverse_sum = 1./sum;

				double lifetime = iso_state_lifetime( ipH_LIKE, nelem, N_(ipHi), L_(ipHi) );
				double percent_error = (1. - inverse_sum/lifetime)*100;
				/* The closed form seems to consistently yield lifetimes roughly 0.2% smaller than the inverse sum
				 * Transition probabilities go like energy cubed.  The cube of the finite mass rydberg is about 0.16% less than unity.  
				 * is that the difference? */
				/* for now enforce to better than 0.5% */
				if( fabs(percent_error) > 0.5 ) 
				{
					fprintf(ioQQQ , 
						"SanityCheck hydrogenic lifetime not consistent with closed form for nelem=%2li upper level=%2i inverse-sum= %.3e lifetime= %.3e error= %.2f%%\n",
						nelem, ipHi, inverse_sum, lifetime, percent_error );
					lgOK = false;
				}
			}
		}
	}


	/* check photo cross sections for H */
	long ipISO = ipH_LIKE;
	nelem = 0;
	for( long n=1; n <= iso_sp[ipISO][nelem].n_HighestResolved_max; ++n )
	{
		// This sanity check is of little value, as the cross-section is calculated 
		// with the same routine here as when populating the opacity stack.  
		// However, the cell mid-point energy used in populating OpacStack can be 
		// greater than the energy implied by the relative energy below. 
		// For very Rydberg, very yrare levels, the cross-section can fall off fast 
		// enough over that energy difference to fail this check.
		// Will reported a sim that failed the check for some levels with n >= 71. 
		// So just don't bother checking above 60.    
		if( n > 60 )
			break;

		double rel_photon_energy = 1. + FLT_EPSILON*2.;
		for( long l=0; l < n; l++ )
		{
			double cs;
			long index = iso_sp[ipISO][nelem].QuantumNumbers2Index[n][l][2];

			cs = H_photo_cs( rel_photon_energy, n, l, nelem+1 );

			error = fabs(cs - opac.OpacStack[iso_sp[ipISO][nelem].fb[index].ipOpac-1] )/
				( (cs + opac.OpacStack[iso_sp[ipISO][nelem].fb[index].ipOpac-1] )/2.);

			if( error > 0.05 )
			{
				fprintf(ioQQQ , "SanityCheck finds insane H photo cs, n,l = %li, %li, expected, found = %e, %e, error = %e\n",
					n, l, opac.OpacStack[iso_sp[ipISO][nelem].fb[index].ipOpac-1], cs, error );
				lgOK = false;
			}
		}
	}
	
	/*********************************************************
	 *                                                       *
	 * confirm that Gaussian integration routines still work *
	 *                                                       *
	 *********************************************************/
	ASSERT( fp_equal( qg32( 0., PI, mySin ) , 2. ) );

	/* And again with the new structure */
	my_Integrand func;
	Integrator<my_Integrand,Gaussian32> mySine;
	ASSERT( fp_equal( mySine.sum( 0., PI, func ), 2. ) );

	/*********************************************************
	 *                                                       *
	 * confirm that exponential integral routines still work *
	 *                                                       *
	 *********************************************************/

	/* check that first and second exponential integrals are ok,
	 * step through range of values, beginning with following */
	x = 1e-3;
	do
	{
		/* check that fast e1 routine is ok */
		ans1 = ee1(x);
		ans2 = expn( 1 , x );
		if( fabs(ans1-ans2)/(ans1+ans2) > 1e-6 )
		{
			fprintf(ioQQQ , "SanityCheck finds insane E1 %g %g %g\n",
				x , ans1 , ans2 );
			lgOK = false;
		}

		/* check that e2 is ok */
		ans1 = e2(x);
		ans2 = expn( 2 , x );
		if( fabs(ans1-ans2)/(ans1+ans2) > 1e-6 )
		{
			fprintf(ioQQQ , "SanityCheck finds insane E2 %g %g %g\n",
				x , ans1 , ans2 );
			lgOK = false;
		}

		/* now increment x */
		x *= 2.;
		/* following limit set by sexp returning zero, used in ee1 */
	} while( x < 64. );

	/*********************************************************
	 *                                                       *
	 * confirm that matrix inversion routine still works     *
	 *                                                       *
	 *********************************************************/

	/* these are the answer, chosen to get xvec 1,2,3 */
	yVector[0] = 1.;
	yVector[1] = 3.;
	yVector[2] = 3.;

	/* zero out the main matrix */
	for(i=0;i<3;++i)
	{
		for( j=0;j<3;++j )
		{
			xMatrix[i][j] = 0.;
		}
	}

	/* remember that order is column, row, alphabetical order, rc */
	xMatrix[0][0] = 1.;
	xMatrix[0][1] = 1.;
	xMatrix[1][1] = 1.;
	xMatrix[2][2] = 1.;

	/* this is the default matrix solver */
	/* this test is the 1-d matrix with 2-d macro simulation */
	/* LDA is right dimension of matrix */

	/* MALLOC space for the  1-d array */
	A = (double*)MALLOC( sizeof(double)*NDIM*NDIM );

	/* copy over the main matrix */
	for( i=0; i < 3; ++i )
	{
		for( j=0; j < 3; ++j )
		{
			A[i*NDIM+j] = xMatrix[i][j];
		}
	}

	ner = 0;

	/*void DGETRF(long,long,double*,long,long[],long*);*/
 	/*void DGETRF(int,int,double*,int,int[],int*);*/
  	getrf_wrapper(3, 3, A, NDIM, ipiv, &ner);
	if( ner != 0 )
	{
		fprintf( ioQQQ, " SanityCheck DGETRF error\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* usage DGETRS, 'N' = no transpose
		* order of matrix,
		* number of cols in bvec, =1
		* array
		* leading dim of array */
 	/*void DGETRS(char,int,int,double*,int,int[],double*,int,int*);*/
	getrs_wrapper('N', 3, 1, A, NDIM, ipiv, yVector, 3, &ner);

	if( ner != 0 )
	{
		fprintf( ioQQQ, " SanityCheck DGETRS error\n" );
		cdEXIT(EXIT_FAILURE);
	}
	/* release the vector */
	free( A );

	/* now check on validity of the solution, demand that this
	 * simple problem have come within a few epsilons of the
	 * correct answer */

	/* find largest deviation */
	rcond = 0.;
	for(i=0;i<3;++i)
	{
		x = fabs( yVector[i]-i-1.);
		rcond = MAX2( rcond, x );
		/*printf(" %g ", yVector[i]);*/
	}

	if( rcond>DBL_EPSILON)
	{
		fprintf(ioQQQ,
			"SanityCheck found too large a deviation in matrix solver = %g \n", 
			rcond);
		/* set flag saying that things are not ok */
		lgOK = false;
	}
	/* end matrix inversion check */

	/* these pointers were set to INT_MIN in ContCreatePointers,
	 * then set to valid numbers in ipShells and OpacityCreate1Element
	 * this checks that all values have been properly filled */
	for( nelem=0; nelem<LIMELM; ++nelem )
	{
		/* must reset state of code after tests performed, remember state here */
		realnum xIonF[NISO][LIMELM];
		double hbn[NISO][LIMELM], hn[NISO][LIMELM];

		if( dense.lgElmtOn[nelem] )
		{
			/* set these abundances so that opacities can be checked below */
			hbn[ipH_LIKE][nelem] = iso_sp[ipH_LIKE][nelem].fb[0].DepartCoef;
			hn[ipH_LIKE][nelem] = iso_sp[ipH_LIKE][nelem].st[0].Pop();
			xIonF[ipH_LIKE][nelem] = dense.xIonDense[nelem][nelem+1-ipH_LIKE];

			iso_sp[ipH_LIKE][nelem].fb[0].DepartCoef = 0.;
			iso_sp[ipH_LIKE][nelem].st[0].Pop() = 1.;
			dense.xIonDense[nelem][nelem+1-ipH_LIKE] = 1.;

			if( nelem > ipHYDROGEN )
			{

				hbn[ipHE_LIKE][nelem] = iso_sp[ipHE_LIKE][nelem].fb[0].DepartCoef;
				hn[ipHE_LIKE][nelem] = iso_sp[ipHE_LIKE][nelem].st[0].Pop();
				xIonF[ipHE_LIKE][nelem] = dense.xIonDense[nelem][nelem+1-ipHE_LIKE];

				/* this does not exist for hydrogen itself */
				iso_sp[ipHE_LIKE][nelem].fb[0].DepartCoef = 0.;
				iso_sp[ipHE_LIKE][nelem].st[0].Pop() = 1.;
				dense.xIonDense[nelem][nelem+1-ipHE_LIKE] = 1.;
			}

			for( ion=0; ion<=nelem; ++ion )
			{
				/* loop over all shells that are defined */
				for( nshells=0; nshells<Heavy.nsShells[nelem][ion]; ++nshells )
				{
					for( j=0; j<3; ++j )
					{
						/* >>chng 00 apr 05, array index is on fortran scale so must be
						 * >= 1.  This test had been <0, correct for C.  Caught by Peter van Hoof */
						if( opac.ipElement[nelem][ion][nshells][j] <=0 )
						{
							/* this is not possible */
							fprintf(ioQQQ,
								"SanityCheck found insane ipElement for nelem=%li ion=%li nshells=%li j=%li \n", 
								nelem , ion , nshells, j );
							fprintf(ioQQQ,
								"value was %li  \n", opac.ipElement[nelem][ion][nshells][j] );
							/* set flag saying that things are not ok */
							lgOK = false;
						}
					}
				}

				if( nelem > 1 )
				{
					vector<realnum> saveion;
					saveion.resize(nelem+2);
					/* check that photoionization cross sections are ok */
					for( j=0; j <= (nelem + 1); j++ )
					{
						saveion[j] = dense.xIonDense[nelem][j];
						dense.xIonDense[nelem][j] = 0.;
					}

					dense.xIonDense[nelem][ion] = 1.;

					OpacityZero();
					opac.lgRedoStatic = true;

					/* generate opacity with standard routine - this is the one
					 * called in OpacityAddTotal to make opacities in usual calculations */
					OpacityAdd1Element(nelem);

					/* this starts one beyond energy of threshold since cs may be zero there */
					for( j=Heavy.ipHeavy[nelem][ion]; j < MIN2(rfield.nflux,continuum.KshellLimit); j++ )
					{
						if( opac.opacity_abs[j]+opac.OpacStatic[j] < FLT_MIN )
						{
							/* this is not possible */
							fprintf(ioQQQ,
								"SanityCheck found non-positive photo cs for nelem=%li ion=%li \n", 
								nelem , ion );
							fprintf(ioQQQ,
								"value was %.2e + %.2e nelem %li ion %li at energy %.2e\n", 
								opac.opacity_abs[j] ,
								opac.OpacStatic[j] ,
								nelem , 
								ion , 
								rfield.anu[j]);
							/* set flag saying that things are not ok */
							lgOK = false;
							break;/**/
						}
					}
					/* reset the ionization distribution */
					for( j=0; j <= (nelem + 1); j++ )
					{
						dense.xIonDense[nelem][j] = saveion[j];
					}

				}
			}
			iso_sp[ipH_LIKE][nelem].fb[ipH1s].DepartCoef = hbn[ipH_LIKE][nelem];
			iso_sp[ipH_LIKE][nelem].st[ipH1s].Pop() = hn[ipH_LIKE][nelem];
			dense.xIonDense[nelem][nelem+1-ipH_LIKE] = xIonF[ipH_LIKE][nelem];

			if( nelem > ipHYDROGEN )
			{
				iso_sp[ipHE_LIKE][nelem].fb[ipHe1s1S].DepartCoef = hbn[ipHE_LIKE][nelem];
				iso_sp[ipHE_LIKE][nelem].st[ipHe1s1S].Pop() = hn[ipHE_LIKE][nelem];
				dense.xIonDense[nelem][nelem+1-ipHE_LIKE] = xIonF[ipHE_LIKE][nelem];
			}
		}
	}


	/*********************************************************
	 *                                                       *
	 * everything is done, all checks make, did we pass them?*
	 *                                                       *
	 *********************************************************/

	if( lgOK )
	{
		/*return if ok */
		if( trace.lgTrace )
		{
			fprintf( ioQQQ, " SanityCheck returns OK\n");
		}
		return;
	}

	else
	{
		/* stop since problem encountered, lgEOF set false */
		fprintf(ioQQQ , "SanityCheck finds insanity so exiting\n");
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}
}
