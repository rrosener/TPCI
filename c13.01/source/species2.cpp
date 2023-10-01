/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "atmdat.h"
#include "phycon.h"
#include "taulines.h"
#include "mole.h"
#include "mole_priv.h"
#include "atoms.h"
#include "string.h"
#include "thirdparty.h"
#include "dense.h"
#include "conv.h"
#include "h2.h"
#include "physconst.h"
#include "secondaries.h"
#include "thermal.h"
#include "cooling.h"
#include "lines_service.h"
#include "atmdat.h"

STATIC double LeidenCollRate(long, long, const TransitionProxy& ,double);
STATIC double StoutCollRate(long ipSpecies, long ipCollider, const TransitionProxy&, double ftemp);
STATIC double ChiantiCollRate(long ipSpecies, long ipCollider, const TransitionProxy&, double ftemp);

static const bool DEBUGSTATE = false;

static double *g, *ex, *pops, *depart, *source, *sink;
static double **AulEscp, **col_str, **AulDest, **AulPump, **CollRate;

/*Solving for the level populations*/

void dBase_solve(void )
{
	realnum abund;
	DEBUG_ENTRY( "dBase_solve()" );
	static bool lgFirstPass = true;
	static long maxNumLevels = 1;
	double totalHeating = 0.;

	if( nSpecies==0 )
		return;

	if( lgFirstPass )
	{
		for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
			maxNumLevels = MAX2( maxNumLevels, dBaseSpecies[ipSpecies].numLevels_max );

		g = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		ex = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		pops = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		depart = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		source = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));
		sink = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));

		AulEscp = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		col_str = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		AulDest = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 
		AulPump = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *));  
		CollRate = (double **)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double *)); 

		for( long j=0; j< maxNumLevels; j++ )
		{
			AulEscp[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			col_str[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			AulDest[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
			AulPump[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double));  
			CollRate[j] = (double *)MALLOC( (unsigned long)(maxNumLevels)*sizeof(double)); 
		}

		lgFirstPass = false;
	}

	// zero all of these values
	memset( g, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( ex, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( pops, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( depart, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( source, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	memset( sink, 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	for( long j=0; j< maxNumLevels; j++ )
	{
		memset( AulEscp[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( col_str[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( AulDest[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( AulPump[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
		memset( CollRate[j], 0, (unsigned long)(maxNumLevels)*sizeof(double) );
	}


	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		const char *spName = dBaseSpecies[ipSpecies].chLabel;
		double cooltl, coolder;
		int nNegPop;
		bool lgZeroPop, lgDeBug = false;

		dBaseSpecies[ipSpecies].CoolTotal = 0.;

#if 0
		//limit for now to small number of levels
		dBaseSpecies[ipSpecies].numLevels_local = MIN2( dBaseSpecies[ipSpecies].numLevels_local, 10 );
#endif

		/* first find current density (cm-3) of species */
		if( dBaseSpecies[ipSpecies].lgMolecular )
		{
			molezone *SpeciesCurrent;
			/** \todo	0	this pointer should be cached one time, and the species
			 * removed from the list if it is not computed */
			if( (SpeciesCurrent = findspecieslocal(dBaseSpecies[ipSpecies].chLabel)) == null_molezone )
			{
				/* did not find the species - print warning for now */
				if( !conv.nTotalIoniz )
					fprintf(ioQQQ," PROBLEM dBase_solve did not find molecular species %li\n",ipSpecies);
			}
			abund = (realnum)SpeciesCurrent->den;
		}
		else
		{
			/* an atom or ion */
			ASSERT( dBaseStates[ipSpecies][0].nelem()<=LIMELM && dBaseStates[ipSpecies][0].IonStg()<=dBaseStates[ipSpecies][0].nelem()+1 );
			abund = dense.xIonDense[ dBaseStates[ipSpecies][0].nelem()-1 ][ dBaseStates[ipSpecies][0].IonStg()-1 ];
		}

		abund *= dBaseSpecies[ipSpecies].fracType * dBaseSpecies[ipSpecies].fracIsotopologue;

		// initialization at start of each iteration
		if( conv.nTotalIoniz == 0)
			dBaseSpecies[ipSpecies].lgActive = true;

		bool lgMakeInActive = (abund <= 1e-20 * dense.xNucleiTotal);
		if( lgMakeInActive && dBaseSpecies[ipSpecies].lgActive )
		{
			// zero out populations and intensities, if previously not set
			dBaseStates[ipSpecies][0].Pop() = 0.;
			for(long ipHi = 1; ipHi < dBaseSpecies[ipSpecies].numLevels_max; ipHi++ )
			{	
				dBaseStates[ipSpecies][ipHi].Pop() = 0.;
# ifdef USE_NLTE7
				dBaseStates[ipSpecies][ipHi].DestCollBB() = 0.;
				dBaseStates[ipSpecies][ipHi].DestPhotoBB() = 0.;
# endif
			}
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				(*tr).Emis().phots() = 0.;
				(*tr).Emis().xIntensity() = 0.;
				(*tr).Coll().col_str() = 0.;
				(*tr).Coll().cool() = 0.;
				(*tr).Coll().heat() = 0.;
			}
			dBaseSpecies[ipSpecies].lgActive = false;
		}

		if( !lgMakeInActive )
			dBaseSpecies[ipSpecies].lgActive = true;

		if( !dBaseSpecies[ipSpecies].lgActive )
			continue;

		// we always hit search phase first, reset number of levels
		if( conv.lgSearch )
			dBaseSpecies[ipSpecies].numLevels_local = dBaseSpecies[ipSpecies].numLevels_max;
		for( long ipLo = 0; ipLo < dBaseSpecies[ipSpecies].numLevels_local; ipLo++ )
		{
			/* statistical weights & Excitation Energies*/
			g[ipLo] = dBaseStates[ipSpecies][ipLo].g() ;
			// parts of the code assert that ground is at zero energy - this is
			// not true for the stored molecular data - so rescale to zero
			ex[ipLo] = dBaseStates[ipSpecies][ipLo].energy().WN() - 
					  dBaseStates[ipSpecies][0].energy().WN();
			/* zero some rates */	
			source[ipLo] = 0.;
			sink[ipLo] = 0.;
		}

		// non-zero was due to roundoff errors on 32-bit
		if( ex[0] <= dBaseStates[ipSpecies][0].energy().WN()* 10. *DBL_EPSILON )
			ex[0] = 0.;
		else
			TotalInsanity();

		for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
			{
				if (ipHi > ipLo)
				{
					AulEscp[ipHi][ipLo] = SMALLFLOAT;
					AulDest[ipHi][ipLo] = SMALLFLOAT;
					AulPump[ipLo][ipHi] = SMALLFLOAT;
				}
				else
				{
					AulEscp[ipHi][ipLo] = 0.;
					AulDest[ipHi][ipLo] = 0.;
					AulPump[ipLo][ipHi] = 0.;
				}
			}
		}

		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
				continue;
			int ipLo = (*tr).ipLo();
			AulEscp[ipHi][ipLo] = (*tr).Emis().Aul()*
				((*tr).Emis().Pesc() + (*tr).Emis().Pelec_esc());
			AulDest[ipHi][ipLo] = (*tr).Emis().Aul()*(*tr).Emis().Pdest();
			AulPump[ipLo][ipHi] = (*tr).Emis().pump();
		}

		/*Setting all the collision strengths and collision rate to zero*/
		for( long ipHi= 0; ipHi<dBaseSpecies[ipSpecies].numLevels_local; ipHi++)
		{
			for( long ipLo= 0; ipLo<dBaseSpecies[ipSpecies].numLevels_local; ipLo++)
			{
				col_str[ipHi][ipLo] = 0.;
				CollRate[ipHi][ipLo] = 0.;
			}
		}

		/*Setting all the collision strengths and collision rate to zero*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			CollisionZero( (*tr).Coll() );
			for( long k=0; k<ipNCOLLIDER; ++k )
				(*tr).Coll().rate_coef_ul_set()[k] = 0.f;
		}

		/* update the collision rates */
		/* molecule */
		if( dBaseSpecies[ipSpecies].lgLAMDA )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					/*using the collision rate coefficients directly*/
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
						(realnum)LeidenCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
			}
		}
		/* Chianti */
		else if( dBaseTrans[ipSpecies].chLabel() == "Chianti" )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				if ((*tr).ipHi() >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
						(realnum)ChiantiCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
			}
		}
		/* Stout */
		else if( dBaseTrans[ipSpecies].chLabel() == "Stout" )
		{
			for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
				 tr != dBaseTrans[ipSpecies].end(); ++tr)
			{
				if ((*tr).ipHi() >= dBaseSpecies[ipSpecies].numLevels_local)
					continue;
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					(*tr).Coll().rate_coef_ul_set()[intCollNo] =
							(realnum)StoutCollRate(ipSpecies, intCollNo, *tr, phycon.te);
				}
			}
		}
		else
			TotalInsanity();
		
		/* guess some missing data */
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			int ipLo = (*tr).ipLo();
			const CollisionProxy &coll_temp = (*tr).Coll();

				/* make educated guesses for some missing data */
			if( dBaseSpecies[ipSpecies].lgMolecular )
			{
				ASSERT( dBaseSpecies[ipSpecies].lgLAMDA );
				/*The collision rate coefficients for helium should not be present and that for molecular hydrogen should be present*/
				if( AtmolCollRateCoeff[ipSpecies][ipATOM_HE].temps.size() == 0 &&
					AtmolCollRateCoeff[ipSpecies][ipH2].temps.size() != 0 )
				{
					coll_temp.rate_coef_ul_set()[ipATOM_HE] = 0.7f * coll_temp.rate_coef_ul()[ipH2];
				}

				/* Put in something for hydrogen collisions if not in database */
				if( AtmolCollRateCoeff[ipSpecies][ipATOM_H].temps.size() == 0 )
				{
					if( AtmolCollRateCoeff[ipSpecies][ipATOM_HE].temps.size() != 0 ) //He0
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 2.0f * coll_temp.rate_coef_ul()[ipATOM_HE];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2_ORTHO].temps.size() != 0 ) //ortho-H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2_ORTHO];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2_PARA].temps.size() != 0 ) //para-H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2_PARA];
					}
					else if( AtmolCollRateCoeff[ipSpecies][ipH2].temps.size() != 0 ) // total H2
					{
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1.4f * coll_temp.rate_coef_ul()[ipH2];
					}
					else
						coll_temp.rate_coef_ul_set()[ipATOM_H] = 1e-13f * (realnum)g[ipLo];
				}

				/* Put in something for proton collisions if not in database */
				if( AtmolCollRateCoeff[ipSpecies][ipPROTON].temps.size() == 0 )
				{
					if( AtmolCollRateCoeff[ipSpecies][ipHE_PLUS].temps.size() != 0 ) //He+
					{
						coll_temp.rate_coef_ul_set()[ipPROTON] = 2.0f * coll_temp.rate_coef_ul()[ipHE_PLUS];
					}
					else
						coll_temp.rate_coef_ul_set()[ipPROTON] = 1e-13f * (realnum)g[ipLo];

				}
				
#if	0
				/* if nothing else has been done, just put a small rate coefficient in */
				for( long intCollNo=0; intCollNo<ipNCOLLIDER; intCollNo++)
				{
					if( coll_temp.rate_coef_ul()[intCollNo] == 0. )
						coll_temp.rate_coef_ul_set()[intCollNo] = 1e-13;
				}
#endif
			}
			else
			{
				/* test for transitions without collision data */
				if( (*tr).Coll().rate_coef_ul_set()[ipELECTRON] == 0. )
				{
					/* Specific transitions of Fe 3,4,and 5 have data that are not in the Chianti/Kurucz files*/
					if( strcmp(dBaseSpecies[ipSpecies].chLabel,"Fe 3") == 0 && ipHi < 14 )
					{
						coll_temp.col_str() = Fe3_cs(ipLo,ipHi);
					}
					else if( strcmp(dBaseSpecies[ipSpecies].chLabel,"Fe 4") == 0 && ipHi < 12 )
					{
						coll_temp.col_str() = Fe4_cs(ipLo,ipHi);
					}
					else if( strcmp(dBaseSpecies[ipSpecies].chLabel,"Fe 5") == 0 && ipHi < 14 )
					{
						coll_temp.col_str() = Fe5_cs(ipLo,ipHi);
					}
					else if( atmdat.lgGbarOn )
					{
						/* All other transitions should use gbar if enabled */
						MakeCS(*tr);
					}
					else
					{
						//If gbar is off, use the default collision strength value.
						coll_temp.col_str() = atmdat.collstrDefault;
					}

					coll_temp.rate_coef_ul_set()[ipELECTRON] = (COLL_CONST*coll_temp.col_str())/
							(g[ipHi]*phycon.sqrte);
				}
			}
		}

		/*Updating the CollRate*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			int ipLo = (*tr).ipLo();
			CollRate[ipHi][ipLo] = (*tr).Coll().ColUL( colliders );
			
			CollRate[ipLo][ipHi] = CollRate[ipHi][ipLo] * g[ipHi] / g[ipLo] *
				dsexp( (*tr).EnergyK() / phycon.te );
			
			/* now add in excitations resulting from cosmic ray secondaries */
			// \todo 2 add branch to do forbidden transitions	
			// this g-bar only works for permitted lines
			if( (*tr).ipCont() > 0 )
			{
				/* get secondaries for all permitted lines by scaling LyA 
				 * excitation by ratio of cross section (oscillator strength/energy) 
				 * Born approximation or plane-wave approximation */
				(*tr).Coll().rate_lu_nontherm_set() = secondaries.x12tot *
					((*tr).Emis().gf()/(*tr).EnergyWN()) /
					(iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).Emis().gf()/
					iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,0).EnergyWN());
					
				CollRate[ipLo][ipHi] += (*tr).Coll().rate_lu_nontherm();
				CollRate[ipHi][ipLo] += (*tr).Coll().rate_lu_nontherm() * (*(*tr).Lo()).g() / (*(*tr).Hi()).g();
			}
		}

		multi_arr<double,2> Cool(dBaseSpecies[ipSpecies].numLevels_local, dBaseSpecies[ipSpecies].numLevels_local);

		/* solve the n-level atom */
		atom_levelN(
			/* dBaseSpecies[ipSpecies].numLevels_local is the number of levels to compute*/ 
			dBaseSpecies[ipSpecies].numLevels_local, 
			/* ABUND is total abundance of species, used for nth equation
			 * if balance equations are homogeneous */
			abund, 
			/* g(dBaseSpecies[ipSpecies].numLevels_local) is statistical weight of levels */
			g, 
			/* EX(dBaseSpecies[ipSpecies].numLevels_local) is excitation potential of levels, deg K or wavenumbers
			 * 0 for lowest level, all are energy rel to ground NOT d(ENER) */
			ex, 
			/* this is 'K' for ex[] as Kelvin deg, is 'w' for wavenumbers */
			'w',
			/* populations [cm-3] of each level as deduced here */
			pops, 
			/* departure coefficient, derived below */
			depart,
			/* net transition rate, A * esc prob, s-1 */
			&AulEscp, 
			/* col str from up to low */
			&col_str, 
			/* AulDest[ihi][ilo] is destruction rate, trans from ihi to ilo, A * dest prob,
			 * asserts confirm that [ihi][ilo] is zero */
			&AulDest, 
			/* AulPump[ilo][ihi] is pumping rate from lower to upper level (s^-1), (hi,lo) must be zero  */
			&AulPump, 
			/* collision rates (s^-1), evaluated here and returned for cooling by calling function,
			 * unless following flag is true.  If true then calling function has already filled
			 * in these rates.  CollRate[ipSpecies][j] is rate from ipSpecies to j */
			&CollRate,
			/* this is an additional creation rate from continuum, normally zero, units cm-3 s-1 */
			source,
			/* this is an additional destruction rate to continuum, normally zero, units s-1 */
			sink,
			// flag saying whether CollRate already done (true), or we need to do it here (false),
			true,
			/* total cooling and its derivative, set here but nothing done with it*/
			&cooltl, 
			&coolder, 
			/* string used to identify calling program in case of error */
			spName, 
			/* nNegPop flag indicating what we have done
			 * positive if negative populations occurred
			 * zero if normal calculation done
			 * negative if too cold (for some atoms other routine will be called in this case) */
			&nNegPop,
			/* true if populations are zero, either due to zero abundance of very low temperature */
			&lgZeroPop,
			/* option to print debug information */
			lgDeBug,
			/* option to do the molecule in LTE */
			dBaseSpecies[ipSpecies].lgLTE,
			/* cooling per line */
			&Cool);

		if( nNegPop > 0 )
		{
			/* negative populations occurred */
			fprintf(ioQQQ," PROBLEM in dBase_solve, atom_levelN returned negative population .\n");
			cdEXIT( EXIT_FAILURE );
		}

		// highest levels may have no population
		while( (pops[dBaseSpecies[ipSpecies].numLevels_local-1]<=0 ) &&
			(dBaseSpecies[ipSpecies].numLevels_local > 1) )
				--dBaseSpecies[ipSpecies].numLevels_local;

		for( long j=0;j< dBaseSpecies[ipSpecies].numLevels_local; j++ )
			dBaseStates[ipSpecies][j].Pop() = MAX2(pops[j],SMALLFLOAT);

		for( long j=dBaseSpecies[ipSpecies].numLevels_local;
			j< dBaseSpecies[ipSpecies].numLevels_max; j++ )
				dBaseStates[ipSpecies][j].Pop() = 0.;

# ifdef USE_NLTE7
		// do sums of totals out (destruction) and into (creation) level
		for( int lvl = 0; lvl < dBaseSpecies[ipSpecies].numLevels_local; lvl++)
		{
			for( int lvl2 = 0; lvl2 < dBaseSpecies[ipSpecies].numLevels_local; lvl2++)
			{
				if( lvl != lvl2 )
				{
					dBaseStates[ipSpecies][lvl].DestCollBB() += CollRate[lvl][lvl2];
					dBaseStates[ipSpecies][lvl].CreatCollBB() += CollRate[lvl2][lvl]*pops[lvl2];
					if( lvl2 < lvl)
						dBaseStates[ipSpecies][lvl].DestPhotoBB() += AulDest[lvl][lvl2];
				}
			}
		}
# endif

		/*Atmol  line*/
		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			int ipLo = (*tr).ipLo();
			(*tr).Coll().cool() = max(Cool[ipHi][ipLo],0.);
			(*tr).Coll().heat() = max(-Cool[ipHi][ipLo],0.);
			
			if ( (*tr).ipCont() > 0 )
			{
			
				(*tr).Emis().phots() =	(*tr).Emis().Aul() * (*(*tr).Hi()).Pop() *
						((*tr).Emis().Pesc() + (*tr).Emis().Pelec_esc());

				(*tr).Emis().xIntensity() = (*tr).Emis().phots() * (*tr).EnergyErg();

					/* population of lower level rel to ion, corrected for stim em */
				(*tr).Emis().PopOpc() = (*(*tr).Lo()).Pop() - (*(*tr).Hi()).Pop()*
					(*(*tr).Lo()).g()/(*(*tr).Hi()).g();

				/* it's possible for this sum to be zero.  Set ratio to zero in that case. */
				if( CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi] > 0. )
				{
					(*tr).Emis().ColOvTot() = CollRate[ipLo][ipHi]/
						(CollRate[ipLo][ipHi]+AulPump[ipLo][ipHi]);
				}
				else
					(*tr).Emis().ColOvTot() = 0.;
				
				// this is only used for save line printout.  Maybe colliders may be involved, but
				// this simple approximation of a "collision strength" should be good enough for
				// the purposes of that printout.
				(*tr).Coll().col_str() = (realnum)( (*tr).Coll().rate_coef_ul()[ipELECTRON] *
					(g[ipHi]*phycon.sqrte)/COLL_CONST);
			}
			else
			{
				(*tr).Coll().col_str() = 0.;
			}
		}

		for(TransitionList::iterator tr = dBaseTrans[ipSpecies].begin();
			 tr != dBaseTrans[ipSpecies].end(); ++tr)
		{
			int ipHi = (*tr).ipHi();
			if (ipHi < dBaseSpecies[ipSpecies].numLevels_local)
				continue;
			EmLineZero( (*tr).Emis() );		
		}

		dBaseSpecies[ipSpecies].CoolTotal = cooltl;
		CoolAdd( dBaseSpecies[ipSpecies].chLabel, 0., max(0.,dBaseSpecies[ipSpecies].CoolTotal) );
		if( dBaseSpecies[ipSpecies].lgMolecular )
			thermal.elementcool[LIMELM] += max(0.,dBaseSpecies[ipSpecies].CoolTotal);
		else
			thermal.elementcool[dBaseStates[ipSpecies][0].nelem()-1] += max(0.,dBaseSpecies[ipSpecies].CoolTotal);
		totalHeating += max(0., -dBaseSpecies[ipSpecies].CoolTotal);
		thermal.dCooldT += coolder;

		/* option to print departure coefficients */
		{
			enum {DEBUG_LOC=false};

			if( DEBUG_LOC )
			{
				fprintf( ioQQQ, " Departure coefficients for species %li\n", ipSpecies );
				for( long j=0; j< dBaseSpecies[ipSpecies].numLevels_local; j++ )
				{
					fprintf( ioQQQ, " level %li \t Depar Coef %e\n", j, depart[j] );
				}
			}
		}
	}

	// total heating for all dBase dBaseSpecies	
	thermal.heating[0][27] = totalHeating;

	return;
}

/*Leiden*/
STATIC double LeidenCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY("LeidenCollRate()");
	double ret_collrate = InterpCollRate( AtmolCollRateCoeff[ipSpecies][ipCollider], tr.ipHi(), tr.ipLo(), ftemp);
	return ret_collrate;
}

/*STOUT*/
STATIC double StoutCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY("StoutCollRate()");

	double rate = 0.;
	// deexcitation rate or collision strength?
	bool lgIsRate = StoutCollData[ipSpecies][tr.ipHi()][tr.ipLo()][ipCollider].lgIsRate;
	int n = StoutCollData[ipSpecies][tr.ipHi()][tr.ipLo()][ipCollider].ntemps;
	if( n < 2)
		return 0.;

	double *x = (double*)MALLOC(n*sizeof(double));
	double *y = (double*)MALLOC(n*sizeof(double));

	double fupsilon = 0.;
	for(int j = 0; j < n; j ++)
	{
		x[j] = StoutCollData[ipSpecies][tr.ipHi()][tr.ipLo()][ipCollider].temps[j];
		y[j] = StoutCollData[ipSpecies][tr.ipHi()][tr.ipLo()][ipCollider].collstrs[j];
		ASSERT( x[j] > 0. && y[j] > 0.);
	}
	//If the temperature is above or below the temperature range, use the CS from the closest temperature.
	//Otherwise, do the linear interpolation.
	if( ftemp < x[0] )
	{
		fupsilon = y[0];
	}
	else if( ftemp > x[n-1] )
	{
		fupsilon = y[n-1];
	}
	else
	{
		fupsilon = linint(x,y,n,ftemp);
	}

	free(x);
	free(y);
	ASSERT(fupsilon > 0);

	/* We can deal with derexcitation rates and collision strengths currently */
	if( lgIsRate )
	{
		rate = fupsilon;
	}
	else
	{
		/* convert the collision strength to a collision rate coefficient */
		/* This formula converting collision strength to collision rate coefficient works fine for the electrons*/
		/* For any other collider the mass would be different*/
		rate = (COLL_CONST*fupsilon)/((*tr.Hi()).g()*sqrt(ftemp));
	}

	return rate;
}

/*CHIANTI*/
STATIC double ChiantiCollRate(long ipSpecies, long ipCollider, const TransitionProxy& tr, double ftemp)
{
	DEBUG_ENTRY("ChiantiCollRate()");

	double rate = 0.;
	double fupsilon = CHIANTI_Upsilon(ipSpecies, ipCollider, tr.ipHi(), tr.ipLo(), ftemp);

	/* NB NB - if proton colliders, the upsilons returned here are actually already rate coefficients. */
	/* these are designated by a collider index and a transition type */
	if( ipCollider == ipPROTON && AtmolCollSplines[ipSpecies][tr.ipHi()][tr.ipLo()][ipCollider].intTranType == 6 )
	{
		rate = fupsilon;
	}
	else if( ipCollider == ipELECTRON )
	{
		/* convert the collision strength to a collision rate coefficient */
		/*This formula converting collision strength to collision rate coefficient works fine for the electrons*/
		/*For any other collider the mass would be different*/
		rate = (COLL_CONST*fupsilon)/((*tr.Hi()).g()*sqrt(ftemp));
	}
	else
		rate = 0.;

	return rate;
}

double CHIANTI_Upsilon(long ipSpecies, long ipCollider, long ipHi, long ipLo, double ftemp)
{
	double fdeltae,fscalingparam,fkte,fxt,fsups,fups;
	int intxs,inttype,intsplinepts;

	DEBUG_ENTRY("CHIANTI_Upsilon()");

	if( AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline == NULL )
	{
		return 0.;
	}

	intsplinepts = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].nSplinePts;
	inttype = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].intTranType;
	fdeltae = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].EnergyDiff;
	fscalingparam = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].ScalingParam;

	fkte = ftemp/fdeltae/1.57888e5;

	/*Way the temperature is scaled*/
	/*Burgess&Tully 1992:Paper gives only types 1 to 4*/
	/*Found that the expressions were the same for 5 & 6 from the associated routine DESCALE_ALL*/
	/*What about 7,8&9?*/
	if( inttype ==1 || inttype==4 )
	{
		fxt = 1-(log(fscalingparam)/(log(fkte+fscalingparam)));
	}
	else if(inttype  == 2 || inttype == 3||inttype == 5 || inttype == 6)
	{
		fxt = fkte/(fkte+fscalingparam);
	}
	else
		TotalInsanity();

	double xs[9],*spl,*spl2;
	/*Creating spline points array*/
	if(intsplinepts == 5)
	{
		spl = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline;
		for(intxs=0;intxs<5;intxs++)
		{
			xs[intxs] = 0.25*intxs;
			if(DEBUGSTATE)
			{
				printf("The xs and spl values are %f and %f \n",xs[intxs],spl[intxs]);
				getchar();
			}
		}
	}
	else if(intsplinepts == 9)
	{
		spl = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].collspline;	
		for( intxs=0; intxs<9; intxs++ )
		{
			xs[intxs] = 0.125*intxs;
			if(DEBUGSTATE)
			{
				printf("The xs and spl values are %f and %f \n",xs[intxs],spl[intxs]);
				getchar();
			}
		}
	}
	else
	{
		TotalInsanity();
	}

	/*Finding the second derivative*/
	spl2 = AtmolCollSplines[ipSpecies][ipHi][ipLo][ipCollider].SplineSecDer;

	if(DEBUGSTATE)
	{
		printf("\n");
		for(intxs=0;intxs<intsplinepts;intxs++)
		{
			printf("The %d value of 2nd derivative is %f \n",intxs+1,spl2[intxs]);
		}
	}

	/*Extracting out the value*/
	//splint(xs,spl,spl2,intsplinepts,fxt,&fsups);

	fsups = linint( xs, spl, intsplinepts, fxt);

	/*Finding upsilon*/
	if(inttype == 1)
	{
		fups = fsups*log(fkte+exp(1.0));
	}
	else if(inttype == 2)
	{
		fups = fsups;
	}
	else if(inttype == 3)
	{
		fups = fsups/(fkte+1.0) ;
	}
	else if(inttype == 4)
	{
		fups = fsups*log(fkte+fscalingparam) ;
	}
	else if(inttype == 5)
	{
		fups = fsups/fkte ;
	}
	else if(inttype == 6)
	{
		fups = pow(10.0,fsups) ;
	}
	else
	{
		TotalInsanity();
	}

	if( fups < 0. ) 
	{
		fprintf( ioQQQ," WARNING: Negative upsilon in species %s, collider %li, indices %4li %4li, Te = %e.\n",
				dBaseSpecies[ipSpecies].chLabel, ipCollider, ipHi, ipLo, ftemp );
		fups = 0.;
	}
	ASSERT(fups>=0);
	return(fups);
}

