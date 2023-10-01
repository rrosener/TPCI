/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "conv.h"
#include "dense.h"
#include "heavy.h"
#include "helike_cs.h"
#include "hydroeinsta.h"
#include "hydrogenic.h"
#include "hydro_vs_rates.h"
#include "ionbal.h"
#include "iso.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h"
#include "rfield.h"
#include "secondaries.h"
#include "trace.h"
#include "taulines.h"

/* These are masses relative to the proton mass of the electron, proton, he+, and alpha particle. */
static double ColliderMass[4] = {ELECTRON_MASS/PROTON_MASS, 1.0, 4.0, 4.0};

void iso_collisional_ionization( long ipISO, long nelem )
{
	ASSERT( ipISO < NISO );

	DEBUG_ENTRY( "iso_collisional_ionization()" );
	
	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* collisional ionization from ground */
	sp->fb[0].ColIoniz = iso_ctrl.lgColl_ionize[ipISO] *
		t_ADfA::Inst().coll_ion_wrapper( nelem, nelem-ipISO, phycon.te );

	iso_put_error(ipISO,nelem,sp->numLevels_max,0,IPCOLLIS,0.20f,0.20f);
	
	for( long ipHi=1; ipHi<sp->numLevels_max; ipHi++ )
	{
		if( nelem == ipISO )
		{
			/* use routine from Vriens and Smeets (1981). */
			/* >>refer	iso neutral	col.ion.	Vriens, L., & Smeets, A.H.M. 1980, Phys Rev A 22, 940 */
			sp->fb[ipHi].ColIoniz = hydro_vs_ioniz( sp->fb[ipHi].xIsoLevNIonRyd, phycon.te );
		}
		else
		{
			/* ions */
			/* use hydrogenic ionization rates for ions
			 * >>refer	iso ions	col.ion.	Allen 1973, Astro. Quan. for low Te.
			 * >>refer	iso ions	col.ion.	Sampson and Zhang 1988, ApJ, 335, 516 for High Te.
			 * */
			sp->fb[ipHi].ColIoniz = 
				Hion_coll_ioniz_ratecoef( ipISO, nelem, N_(ipHi), sp->fb[ipHi].xIsoLevNIonRyd, phycon.te  );
		}
		
		// iso_ctrl.lgColl_ionize is option to turn off collisions, "atom XX-like collis off" command
		sp->fb[ipHi].ColIoniz *= iso_ctrl.lgColl_ionize[ipISO];
 	
		iso_put_error(ipISO,nelem,sp->numLevels_max,ipHi,IPCOLLIS,0.20f,0.20f);
	}
	
	/* Here we arbitrarily scale the highest level ionization to account for the fact
	 * that, if the atom is not full size, this level should be interacting with higher
	 * levels and not just the continuum.  We did add on collisional excitation terms instead
	 * but this caused a problem at low temperatures because the collisional ionization was 
	 * a sum of terms with different Boltzmann factors, while PopLTE had just one Boltzmann
	 * factor.  The result was a collisional recombination that had residual exponentials of
	 * the form exp(x/kT), which blew up at small T.	*/
	if( 0 && !sp->lgLevelsLowered )
	{
		sp->fb[sp->numLevels_max-1].ColIoniz *= 100.;
	}
	
	return;
}

void iso_suprathermal( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_suprathermal()" );

	/* check that we were called with valid parameters */
	ASSERT( ipISO < NISO );
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );

	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/***********************************************************************
	 *                                                                     *
	 * get secondary excitation by suprathermal electrons                  *
	 *                                                                     *
	 ***********************************************************************/

	for( long i=1; i < sp->numLevels_max; i++ )
	{
		if( sp->trans(i,0).ipCont() > 0 )
		{
			/* get secondaries for all iso lines by scaling LyA 
			 * excitation by ratio of cross section (oscillator strength/energy) 
			 * Born approximation or plane-wave approximation based on
			 *>>refer	HI	excitation	Shemansky, D.E., et al., 1985, ApJ, 296, 774 */
			sp->trans(i,0).Coll().rate_lu_nontherm_set() = secondaries.x12tot *
				(sp->trans(i,0).Emis().gf()/
				sp->trans(i,0).EnergyWN()) /
				(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).Emis().gf()/
				iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).EnergyWN()) *
				iso_ctrl.lgColl_excite[ipISO];;
		}
		else
			sp->trans(i,0).Coll().rate_lu_nontherm_set() = 0.;
	}

	return;
}

/*============================*/
/* evaluate collisional rates */
void iso_collide( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collide()" );

	/* this array stores the last temperature at which collision data were evaluated for
	 * each species of the isoelectronic sequence. */
	static double TeUsed[NISO][LIMELM]={
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0} };

	/* check that we were called with valid parameters */
	ASSERT( ipISO < NISO );
	ASSERT( nelem >= ipISO );
	ASSERT( nelem < LIMELM );
	
	t_iso_sp* sp = &iso_sp[ipISO][nelem];

	/* skip most of this routine if temperature has not changed,
	 * the check on conv.nTotalIoniz is to make sure that we redo this
	 * on the very first call in a grid calc - it is 0 on the first call */
	if( fp_equal( TeUsed[ipISO][nelem], phycon.te ) && conv.nTotalIoniz )
	{
		ASSERT( sp->trans( iso_ctrl.nLyaLevel[ipISO], 0 ).Coll().ColUL( colliders ) >= 0. );

		if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
		{
			fprintf( ioQQQ, 
				"     iso_collide called %s nelem %li - no reeval Boltz fac, LTE dens\n",
				iso_ctrl.chISO[ipISO], nelem );
		}
	}
	else
	{
		TeUsed[ipISO][nelem] = phycon.te;

		if( trace.lgTrace  && (trace.lgHBug||trace.lgHeBug) )
		{
			fprintf( ioQQQ, 
				"     iso_collide called %s nelem %li - will reeval Boltz fac, LTE dens\n",
				iso_ctrl.chISO[ipISO], nelem );
		}
		
		/**********************************************************
		 *                                                        *
		 * Boltzmann factors for all levels,  and                 *
		 * collisional ionization and excitation                  *
		 *                                                        *
		 **********************************************************/

		/* HION_LTE_POP	is planck^2 / (2 pi m_e k ), raised to 3/2 next */
		double factor = HION_LTE_POP*dense.AtomicWeight[nelem]/
			(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);

		/* term in () is stat weight of electron * ion */
		double ConvLTEPOP = pow(factor,1.5)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;

		sp->lgPopLTE_OK = true;

		// this is the maximum value of iso.PopLTE (units cm^3) that corresponds
		// to the minimum positive density values.  A smaller density will be
		// regarded as zero, and the product PopLTE*n_e*n_Z+ will also be zero.
#define MAX_POP_LTE	(MAX_DENSITY/dense.density_low_limit/dense.density_low_limit)

		/* fully define Boltzmann factors to continuum for model levels */
		for( long ipLo=0; ipLo<sp->numLevels_max; ipLo++ )
		{
			/* this Boltzmann factor is exp( +ioniz energy / Te ) */
			sp->st[ipLo].Boltzmann() = 
				dsexp(sp->st[ipLo].energy().Ryd()/phycon.te_ryd);
			sp->st[ipLo].ConBoltz() =
				dsexp(sp->fb[ipLo].xIsoLevNIonRyd/phycon.te_ryd);

			/***************************************
			 *                                     *
			 * LTE abundances for all levels       *
			 *                                     *
			 ***************************************/

			if( sp->st[ipLo].ConBoltz() > SMALLDOUBLE )
			{
				/* LTE population of given level. */
				sp->fb[ipLo].PopLTE = 
					sp->st[ipLo].g() / sp->st[ipLo].ConBoltz() * ConvLTEPOP;
				ASSERT( sp->fb[ipLo].PopLTE < BIGDOUBLE );
			}
			else
			{
				sp->fb[ipLo].PopLTE = 0.;
			}

			sp->fb[ipLo].PopLTE = MIN2( sp->fb[ipLo].PopLTE, MAX_POP_LTE );

			/* now check for any zeros - if present then matrix cannot be used */
			if( sp->fb[ipLo].PopLTE <= 0. )
			{
				sp->lgPopLTE_OK = false;
			}
		}

		iso_collisional_ionization( ipISO, nelem );

		/***********************************************************
		 *                                                         *
		 * collisional deexcitation for all lines in iso sequence  *
		 *                                                         *
		 ***********************************************************/

		if( iso_ctrl.lgColl_excite[ipISO] )
		{
			for( long ipHi=1; ipHi<sp->numLevels_max; ipHi++ )
			{
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{
					for( long ipCollider = ipELECTRON; ipCollider <= ipALPHA; ipCollider++ )
					{
						double cs_temp = 0.;

						if( N_(ipHi) == N_(ipLo) && !iso_ctrl.lgColl_l_mixing[ipISO] )
							cs_temp = 0.;
						else if( N_(ipHi)-N_(ipLo) > 2 && ipCollider > ipELECTRON )
							cs_temp = 0.;
						else if( ipISO == ipH_LIKE )
							cs_temp =  HydroCSInterp( nelem , ipHi , ipLo, ipCollider );
						else if( ipISO == ipHE_LIKE )
							cs_temp = HeCSInterp( nelem , ipHi , ipLo, ipCollider );
						else
							TotalInsanity();

						cs_temp *=  iso_ctrl.lgColl_excite[ipISO];
					
						if( opac.lgCaseB_HummerStorey && N_(ipHi)==N_(ipLo)+1 && abs(L_(ipHi)-L_(ipLo))==1 )
						{
							double Aul = HydroEinstA( N_(ipLo), N_(ipHi) );
							cs_temp *= ( sp->trans(ipHi,ipLo).Emis().gf() / sp->st[ipLo].g() ) /
								( GetGF(Aul, sp->trans(ipHi,ipLo).EnergyWN(), 2.*N_(ipHi)*N_(ipHi))/(2.*N_(ipLo)*N_(ipLo)) );
						}
	
						// store electron collision strength in generic collision strength
						if( ipCollider == ipELECTRON )
							sp->trans(ipHi,ipLo).Coll().col_str() = (realnum) cs_temp;

						double reduced_mass_collider_system = dense.AtomicWeight[nelem]*ColliderMass[ipCollider]/
							(dense.AtomicWeight[nelem]+ColliderMass[ipCollider])*ATOMIC_MASS_UNIT;

						if( ipCollider == ipELECTRON )
							reduced_mass_collider_system  = ELECTRON_MASS;

						double rateCoef = cs_temp *
								pow(ELECTRON_MASS/reduced_mass_collider_system, 1.5) * COLL_CONST/
								( phycon.sqrte * (double)sp->st[ipHi].g() );
					
						if( !opac.lgCaseB_HummerStorey )
						{	
							if( ipISO == ipH_LIKE )
							{
								if( N_(ipHi) > sp->n_HighestResolved_max &&
									N_(ipLo) <= sp->n_HighestResolved_max )
								{
									rateCoef *= (8./3.)*(log(double(N_(ipHi)))+2.);
								}
							}
							else if( ipISO == ipHE_LIKE )
							{
								fixit();
								/* This is intended to be a trick to get the correct collisional excitation from
								 * collapsed levels to resolved levels.  Is it needed or does the stat weight used
								 * above handle it automatically? If it is needed, is this correct? */
								if( N_(ipHi) > sp->n_HighestResolved_max &&
									N_(ipLo) <= sp->n_HighestResolved_max )
								{
									rateCoef *= (2./3.)*(log(double(N_(ipHi)))+2.);
								}
							}
						}
	
						sp->trans(ipHi,ipLo).Coll().rate_coef_ul_set()[ipCollider] =
								(realnum) rateCoef;
					}
					
					/* check for sanity */
					ASSERT( sp->trans(ipHi,ipLo).Coll().rate_coef_ul()[ipELECTRON] >= 0. );

					if( N_(ipHi) <= 5 && N_(ipLo) <= 2 )
						iso_put_error( ipISO, nelem, ipHi, ipLo, IPCOLLIS, 0.10f, 0.30f );
					else
						iso_put_error( ipISO, nelem, ipHi, ipLo, IPCOLLIS, 0.20f, 0.30f );
				}
			}
		}
		
		if( (trace.lgTrace && trace.lgIsoTraceFull[ipISO]) && (nelem == trace.ipIsoTrace[ipISO]) )
		{
			fprintf( ioQQQ, "     iso_collide: %s Z=%li de-excitation rates coefficients\n", iso_ctrl.chISO[ipISO], nelem + 1 );
			long upper_limit = sp->numLevels_local;
			for( long ipHi=1; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ, " %li\t", ipHi );
				for( long ipLo=0; ipLo < ipHi; ipLo++ )
				{
					fprintf( ioQQQ,PrintEfmt("%10.3e", 
							sp->trans(ipHi,ipLo).Coll().ColUL( colliders ) / dense.EdenHCorr ));
				}
				fprintf( ioQQQ, "\n" );
			}

			fprintf( ioQQQ, "     iso_collide: %s Z=%li collisional ionization coefficients\n", iso_ctrl.chISO[ipISO], nelem + 1 );
			for( long ipHi=0; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  sp->fb[ipHi].ColIoniz ));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "     iso_collide: %s Z=%li continuum boltzmann factor\n", iso_ctrl.chISO[ipISO], nelem + 1 );
			for( long ipHi=0; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  sp->st[ipHi].ConBoltz() ));
			}
			fprintf( ioQQQ, "\n" );

			fprintf( ioQQQ, "     iso_collide: %s Z=%li continuum boltzmann factor\n", iso_ctrl.chISO[ipISO], nelem + 1 );
			for( long ipHi=0; ipHi < upper_limit; ipHi++ )
			{
				fprintf( ioQQQ,PrintEfmt("%10.3e",  sp->fb[ipHi].PopLTE ));
			}
			fprintf( ioQQQ, "\n" );
		}

		/* the case b hummer and storey option,
		 * this kills collisional excitation and ionization from n=1 and n=2 */
		if( opac.lgCaseB_HummerStorey )
		{
			for( long ipLo=0; ipLo<sp->numLevels_max-1; ipLo++ )
			{
				if( N_(ipLo)>=3 )
					break;
				
				sp->fb[ipLo].ColIoniz = 0.;

				for( long ipHi=ipLo+1; ipHi<sp->numLevels_max; ipHi++ )
				{
					/* don't disable 2-2 collisions */
					if( N_(ipLo)==2 && N_(ipHi)==2 )
						continue;

					sp->trans(ipHi,ipLo).Coll().col_str() = 0.;
					for( long k=0; k<ipNCOLLIDER; ++k )
						sp->trans(ipHi,ipLo).Coll().rate_coef_ul_set()[k] = 0.;
				}
			}
		}
	}

	iso_suprathermal( ipISO, nelem );

	/* this must be reevaluated every time since eden can change when Te does not */
	/* save into main array - collisional ionization by thermal electrons */
	ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0] = 
		sp->fb[0].ColIoniz*dense.EdenHCorr;

	/* cooling due to collisional ionization, which only includes thermal electrons */
	ionbal.CollIonRate_Ground[nelem][nelem-ipISO][1] = 
		ionbal.CollIonRate_Ground[nelem][nelem-ipISO][0]*
		rfield.anu[Heavy.ipHeavy[nelem][nelem-ipISO]-1]*EN1RYD;

	return;
}
