/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*iso_create create data for hydrogen and helium, 1 per coreload, called by ContCreatePointers 
 * in turn called after commands parsed */
/*iso_zero zero data for hydrogen and helium */
#include "cddefines.h"
#include "atmdat.h"
#include "dense.h"
#include "elementnames.h"
#include "helike.h"
#include "helike_einsta.h"
#include "hydro_bauman.h"
#include "hydrogenic.h"
#include "hydroeinsta.h"
#include "iso.h"
#include "lines_service.h"
#include "opacity.h"
#include "phycon.h"
#include "physconst.h"
#include "secondaries.h"
#include "taulines.h"
#include "thirdparty.h"

/*iso_zero zero data for hydrogen and helium */
STATIC void iso_zero(void);

/* allocate memory for iso sequence structures */
STATIC void iso_allocate(void);

/* define levels of iso sequences and assign quantum numbers to those levels */
STATIC void iso_assign_quantum_numbers(void);

STATIC void FillExtraLymanLine( const TransitionList::iterator& t, long ipISO, long nelem, long nHi );

STATIC void iso_satellite( void );

char chL[21]={'S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U','V','W','X','Y','Z'};

void iso_create(void)
{
	long int ipHi, 
		ipLo;

	static int nCalled = 0;

	double HIonPoten;

	DEBUG_ENTRY( "iso_create()" );

	/* > 1 if not first call, then just zero arrays out */
	if( nCalled > 0 )
	{
		iso_zero();
		return;
	}

	/* this is first call, increment the nCalled counterso never do this again */
	++nCalled;

	/* these are the statistical weights of the ions */
	iso_ctrl.stat_ion[ipH_LIKE] = 1.f;
	iso_ctrl.stat_ion[ipHE_LIKE] = 2.f;

	/* this routine allocates all the memory
	 * needed for the iso sequence structures */
	iso_allocate();

	/* loop over iso sequences and assign quantum numbers to all levels */
	iso_assign_quantum_numbers();

	/* this is a dummy line, junk it too. */
	(*TauDummy).Junk();
	(*TauDummy).AddHiState();
	(*TauDummy).AddLoState();
	(*TauDummy).AddLine2Stack();

	/********************************************/
	/**********  Line and level energies ********/
	/********************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		/* main hydrogenic arrays, fill with sane values */
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* Dima's array has ionization potentials in eV, but not on same
				 * scale as cloudy itself*/
				/* extra factor accounts for this */
				HIonPoten = t_ADfA::Inst().ph1(0,0,nelem,0)/EVRYD* 0.9998787;
				ASSERT(HIonPoten > 0.);

				double EnergyRydGround = 0.;
				/* go from ground to the highest level */
				for( ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					double EnergyWN, EnergyRyd;

					if( ipISO == ipH_LIKE )
					{
						EnergyRyd = HIonPoten/POW2((double)N_(ipHi));
					}
					else if( ipISO == ipHE_LIKE )
					{
						EnergyRyd = helike_energy( nelem, ipHi ) * WAVNRYD;
					}
					else
					{
						/* Other iso sequences don't exist yet. */
						TotalInsanity();
					}

					/* >>chng 02 feb 09, change test to >= 0 since we now use 0 for 2s-2p */
					ASSERT(EnergyRyd >= 0.);

					iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd = EnergyRyd;
					if (ipHi == 0)
						EnergyRydGround = EnergyRyd;
					iso_sp[ipISO][nelem].st[ipHi].energy().set(EnergyRydGround-EnergyRyd);

					/* now loop from ground to level below ipHi */
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						EnergyWN = RYD_INF * (iso_sp[ipISO][nelem].fb[ipLo].xIsoLevNIonRyd -
							iso_sp[ipISO][nelem].fb[ipHi].xIsoLevNIonRyd);

						/* This is the minimum line energy we will allow.  */
						/* \todo 2 wire this to lowest energy of code. */
						if( EnergyWN==0 && ipISO==ipHE_LIKE )
							EnergyWN = 0.0001;

						if( EnergyWN < 0. )
							EnergyWN = -1.0 * EnergyWN;

						/* transition energy in various units: */
						iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() = (realnum)EnergyWN;

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() >= 0.);
						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyErg() >= 0.);
						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyK() >= 0.);

						/** \todo 2 this should be changed for j-resolved levels */
						if( N_(ipLo)==N_(ipHi) && ipISO==ipH_LIKE )
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng() = 0.;
						}
						else
						{
							/* make following an air wavelength */
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng() = 
								(realnum)(1.0e8/
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN()/
								RefIndex( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN()));
							ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).WLAng() > 0.);
						}

					}
				}

				/* fill the extra Lyman lines */
				for( ipHi=2; ipHi < iso_ctrl.nLyman_malloc[ipISO]; ipHi++ )
				{
					FillExtraLymanLine( ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi], ipISO, nelem, ipHi );
				}
			}
		}
	}

	/***************************************************************/
	/***** Set up recombination tables for later interpolation *****/
	/***************************************************************/
	/* NB - the above is all we need if we are compiling recombination tables. */
	iso_recomb_malloc();
	iso_recomb_setup( ipH_LIKE );
	iso_recomb_setup( ipHE_LIKE );
	iso_recomb_auxiliary_free();

	/* set up helium collision tables */
	HeCollidSetup();

	/***********************************************************************************/
	/**********  Transition Probabilities, Redistribution Functions, Opacitites ********/
	/***********************************************************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		if( ipISO == ipH_LIKE )
		{
			/* do nothing here */
		}
		else if( ipISO == ipHE_LIKE )
		{
			/* This routine reads in transition probabilities from a file. */ 
			HelikeTransProbSetup();
		}
		else
		{
			TotalInsanity();
		}

		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < (iso_sp[ipISO][nelem].numLevels_max - 1); ipLo++ )
				{
					for( ipHi=ipLo + 1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
					{
						realnum Aul;

						/* transition prob, EinstA uses current H atom indices */
						if( ipISO == ipH_LIKE )
						{
							Aul = hydro_transprob( nelem, ipHi, ipLo );
						}
						else if( ipISO == ipHE_LIKE )
						{
							Aul = helike_transprob(nelem, ipHi, ipLo);
						}
						else
						{
							TotalInsanity();
						}

						if( Aul <= iso_ctrl.SmallA )
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipEmis() = -1;
						else
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).AddLine2Stack();

						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() = Aul;

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() > 0.);

						if( ipLo == 0 && ipHi == iso_ctrl.nLyaLevel[ipISO] )
						{
							long redis = iso_ctrl.ipLyaRedist[ipISO];
							// H LyA has a special redistribution function 
							if( ipISO==ipH_LIKE && nelem==ipHYDROGEN )
								redis = ipLY_A;
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = redis;
						}
						else if( ipLo == 0 )
						{
							/* these are rest of Lyman lines, 
							 * complete redistribution, doppler core only, K2 core, default ipCRD */
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = iso_ctrl.ipResoRedist[ipISO];
						}
						else
						{
							/* all lines coming from excited states, default is complete
							 * redis with wings, ipCRDW*/
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().iRedisFun() = iso_ctrl.ipSubRedist[ipISO];
						}

						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA ||
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() <= 0.)
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() = 0.;
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() = 0.;
						}
						else
						{
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() = 
								(realnum)(GetGF(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul(),
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN(),
								iso_sp[ipISO][nelem].st[ipHi].g()));
							ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf() > 0.);

							/* derive the abs coef, call to function is gf, wl (A), g_low */
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() = 
								(realnum)(abscf(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().gf(),
								iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN(),
								iso_sp[ipISO][nelem].st[ipLo].g()));
							ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().opacity() > 0.);
						}
					}
				}
			}
		}
	}

	/************************************************/
	/**********  Fine Structure Mixing - FSM ********/
	/************************************************/
	if( iso_ctrl.lgFSM[ipHE_LIKE] )
	{
		/* set some special optical depth values */
		for( long nelem=ipHE_LIKE; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipHi=ipHe2s3S; ipHi<iso_sp[ipHE_LIKE][nelem].numLevels_max; ipHi++ )
				{
					for( ipLo=ipHe1s1S; ipLo<ipHi; ipLo++ )
					{
						DoFSMixing( nelem, ipLo, ipHi );
					}
				}
			}
		}
	}

	/* following comes out very slightly off, correct here */
	iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3s,ipH2s).WLAng() = 1640.f;
	iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3s,ipH2p).WLAng() = 1640.f;
	if( iso_sp[ipH_LIKE][ipHELIUM].n_HighestResolved_max >=3 )
	{
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3p,ipH2s).WLAng() = 1640.f;
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3p,ipH2p).WLAng() = 1640.f;
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3d,ipH2s).WLAng() = 1640.f;
		iso_sp[ipH_LIKE][ipHELIUM].trans(ipH3d,ipH2p).WLAng() = 1640.f;
	}

	/****************************************************/
	/**********  lifetimes and damping constants ********/
	/****************************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* define these for H and He always */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				/* these are not defined and must never be used */
				iso_sp[ipISO][nelem].st[0].lifetime() = -FLT_MAX;

				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					iso_sp[ipISO][nelem].st[ipHi].lifetime() = SMALLFLOAT;
					
					for( ipLo=0; ipLo < ipHi; ipLo++ )  
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						iso_sp[ipISO][nelem].st[ipHi].lifetime() += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();
					}

					/* sum of A's was just stuffed, now invert for lifetime. */
					iso_sp[ipISO][nelem].st[ipHi].lifetime() = 1./iso_sp[ipISO][nelem].st[ipHi].lifetime();

					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() <= 0. )
							continue;

						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
							continue;

						iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel() = (realnum)( 
							(1.f/iso_sp[ipISO][nelem].st[ipHi].lifetime())/
							PI4/iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN());

						ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel()> 0.);
					}
				}
			}
		}
	}

	/* zero out some line information */
	iso_zero();

	/* loop over iso sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem == ipISO || dense.lgElmtOn[nelem] ) 
			{
				/* calculate cascade probabilities, branching ratios, and associated errors. */
				iso_cascade( ipISO, nelem);
			}
		}
	}

	iso_satellite();

	for( long nelem=ipHYDROGEN; nelem < LIMELM; ++nelem )
		iso_satellite_update( nelem );

	/***************************************/
	/**********  Stark Broadening **********/
	/***************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( long ipHi= 1; ipHi < iso_sp[ipISO][nelem].numLevels_max; ++ipHi )
				{
					for( long ipLo=0; ipLo < ipHi; ++ipLo )
					{
						iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk = 0.;
						iso_sp[ipISO][nelem].ex[ipHi][ipLo].pestrk_up = 0.;
					}
				}
			}
		}
	}

	return;
}

/* ============================================================================== */
STATIC void iso_zero(void)
{
	DEBUG_ENTRY( "iso_zero()" );

	hydro.HLineWidth = 0.;

	/****************************************************/
	/**********  initialize some variables **********/
	/****************************************************/
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( long ipHi=0; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					iso_sp[ipISO][nelem].st[ipHi].Pop() = 0.;
					iso_sp[ipISO][nelem].fb[ipHi].Reset();
				}
				if (ipISO <= nelem) 
					iso_sp[ipISO][nelem].st[0].Pop() =  
						dense.xIonDense[nelem][nelem-ipISO]; 
			}

			if( nelem < 2 )
			{			
				iso_collapsed_bnl_set( ipISO, nelem );
				//iso_collapsed_bnl_print( ipISO, nelem );
				iso_collapsed_Aul_update( ipISO, nelem );
				iso_collapsed_lifetimes_update( ipISO, nelem );
			}
		}
	}

	/* ground state of H and He is different since totally determine
	 * their own opacities */
	iso_sp[ipH_LIKE][ipHYDROGEN].fb[0].ConOpacRatio = 1e-5;
	iso_sp[ipH_LIKE][ipHELIUM].fb[0].ConOpacRatio = 1e-5;
	iso_sp[ipHE_LIKE][ipHELIUM].fb[0].ConOpacRatio = 1e-5;

	return;
}

STATIC void iso_allocate(void)
{

	DEBUG_ENTRY( "iso_allocate()" );

	/* the hydrogen and helium like iso-sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				t_iso_sp *sp = &iso_sp[ipISO][nelem];
				sp->numLevels_malloc = sp->numLevels_max;

				ASSERT( sp->numLevels_max > 0 );
				ASSERT( iso_ctrl.nLyman_malloc[ipISO] == iso_ctrl.nLyman[ipISO] );

				sp->CachedAs.reserve( MAX2(1, sp->nCollapsed_max) );

				sp->ipTrans.reserve( sp->numLevels_max );
				sp->ex.reserve( sp->numLevels_max );
				sp->CascadeProb.reserve( sp->numLevels_max );
				sp->BranchRatio.reserve( sp->numLevels_max );
				//sp->st.resize( sp->numLevels_max );
				sp->fb.resize( sp->numLevels_max );

				for( long i = 0; i < sp->nCollapsed_max; ++i )
				{
					sp->CachedAs.reserve( i, sp->numLevels_max - sp->nCollapsed_max );
					for( long i1 = 0; i1 < sp->numLevels_max - sp->nCollapsed_max; ++i1 )
					{
						/* allocate two spaces delta L +/- 1 */
						sp->CachedAs.reserve( i, i1, 2 );
					}
				}

				sp->QuantumNumbers2Index.reserve( sp->n_HighestResolved_max + sp->nCollapsed_max + 1 );

				for( long i = 1; i <= sp->n_HighestResolved_max + sp->nCollapsed_max; ++i )
				{
					/* Allocate proper number of angular momentum quantum number.	*/
					sp->QuantumNumbers2Index.reserve( i, i );

					for( long i1 = 0; i1 < i; ++i1 )
					{
						/* This may have to change for other iso sequences. */
						ASSERT( NISO == 2 );
						/* Allocate 4 spaces for multiplicity.	H-like will be accessed with "2" for doublet,
						 * He-like will be accessed via "1" for singlet or "3" for triplet.  "0" will not be used. */
						sp->QuantumNumbers2Index.reserve( i, i1, 4 );
					}
				}
				
				for( long n=1; n < sp->numLevels_max; ++n )
				{
					sp->ipTrans.reserve( n, n );
				}

				for( long n=0; n < sp->numLevels_max; ++n )
				{
					sp->ex.reserve( n, sp->numLevels_max );
					sp->CascadeProb.reserve( n, sp->numLevels_max );
					sp->BranchRatio.reserve( n, sp->numLevels_max );
				}

				sp->ipTrans.alloc();
				sp->ex.alloc();
				sp->CascadeProb.alloc();
				sp->BranchRatio.alloc();

				sp->CachedAs.alloc();
				sp->QuantumNumbers2Index.alloc();
				sp->bnl_effective.alloc( sp->QuantumNumbers2Index.clone() );
				sp->QuantumNumbers2Index.invalidate();
				sp->bnl_effective.invalidate();
			}
		}
	}

	ipSatelliteLines.reserve( NISO );
	ipExtraLymanLines.reserve( NISO );

	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		ipSatelliteLines.reserve( ipISO, LIMELM );
		ipExtraLymanLines.reserve( ipISO, LIMELM );

		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				ASSERT( iso_sp[ipISO][nelem].numLevels_max > 0 );

				ipSatelliteLines.reserve( ipISO, nelem, iso_sp[ipISO][nelem].numLevels_max );
				ipExtraLymanLines.reserve( ipISO, nelem, iso_ctrl.nLyman_malloc[ipISO] );
			}
		}
	}

	ipSatelliteLines.alloc();
	ipExtraLymanLines.alloc();

	Transitions.resize(NISO);
	SatelliteLines.resize(NISO);
	ExtraLymanLines.resize(NISO);
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		Transitions[ipISO].reserve(LIMELM);
		SatelliteLines[ipISO].reserve(LIMELM);
		ExtraLymanLines[ipISO].reserve(LIMELM);
		for( long nelem=0; nelem < ipISO; ++nelem )
		{
			Transitions[ipISO].push_back(
				TransitionList("Insanity",&AnonStates));
			SatelliteLines[ipISO].push_back(
				TransitionList("Insanity",&AnonStates));
			ExtraLymanLines[ipISO].push_back(
				TransitionList("Insanity",&AnonStates));
		}
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				Transitions[ipISO].push_back(
					TransitionList("Isosequence",&iso_sp[ipISO][nelem].st));
				SatelliteLines[ipISO].push_back(
					TransitionList("SatelliteLines",&iso_sp[ipISO][nelem].st));
				ExtraLymanLines[ipISO].push_back(
					TransitionList("ExtraLymanLines",&iso_sp[ipISO][nelem].st));
			}
			else
			{
				Transitions[ipISO].push_back(
					TransitionList("Insanity",&AnonStates));
				SatelliteLines[ipISO].push_back(
					TransitionList("Insanity",&AnonStates));
				ExtraLymanLines[ipISO].push_back(
					TransitionList("Insanity",&AnonStates));
			}
		}
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			/* only grab core for elements that are turned on */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					SatelliteLines[ipISO][nelem].resize( iso_sp[ipISO][nelem].numLevels_max );
					AllTransitions.push_back(SatelliteLines[ipISO][nelem]);
					unsigned int nLine = 0;
					for( long ipLo=0; ipLo<iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
					{
						/* Upper level is continuum, use a generic state
						 * lower level is the same as the index. */
						ipSatelliteLines[ipISO][nelem][ipLo] = nLine;
						SatelliteLines[ipISO][nelem][nLine].Junk();
						long ipHi = iso_sp[ipISO][nelem].numLevels_max;
						SatelliteLines[ipISO][nelem][nLine].setHi(ipHi);
						SatelliteLines[ipISO][nelem][nLine].setLo(ipLo);
						SatelliteLines[ipISO][nelem][nLine].AddLine2Stack();
						++nLine;
					}
					ASSERT(SatelliteLines[ipISO][nelem].size() == nLine);
				}

				//iso_sp[ipISO][nelem].tr.resize( iso_sp[ipISO][nelem].ipTrans.size() );
				//iso_sp[ipISO][nelem].tr.states() = &iso_sp[ipISO][nelem].st;
				Transitions[ipISO][nelem].resize( iso_sp[ipISO][nelem].ipTrans.size() );
				AllTransitions.push_back(Transitions[ipISO][nelem]);
				unsigned int nTransition=0;
				for( long ipHi=1; ipHi<iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
				{
					for( long ipLo=0; ipLo < ipHi; ipLo++ )
					{
						/* set ENTIRE array to impossible values, in case of bad pointer */
						iso_sp[ipISO][nelem].ipTrans[ipHi][ipLo] = nTransition;
						Transitions[ipISO][nelem][nTransition].Junk();
						Transitions[ipISO][nelem][nTransition].setHi(ipHi);
						Transitions[ipISO][nelem][nTransition].setLo(ipLo);
						++nTransition;
					}
				}
				ASSERT(Transitions[ipISO][nelem].size() == nTransition);
				iso_sp[ipISO][nelem].tr = &Transitions[ipISO][nelem];

				/* junk the extra Lyman lines */
				AllTransitions.push_back(ExtraLymanLines[ipISO][nelem]);
				ExtraLymanLines[ipISO][nelem].resize(iso_ctrl.nLyman_malloc[ipISO]-2);
				ExtraLymanLines[ipISO][nelem].states() = &iso_sp[ipISO][nelem].st;
				unsigned int nExtraLyman = 0;
				for( long ipHi=2; ipHi < iso_ctrl.nLyman_malloc[ipISO]; ipHi++ )
				{
					ipExtraLymanLines[ipISO][nelem][ipHi] = nExtraLyman;
					ExtraLymanLines[ipISO][nelem][nExtraLyman].Junk();
					long ipHi_offset = iso_sp[ipISO][nelem].numLevels_max + ipHi - 2;
					if( iso_ctrl.lgDielRecom[ipISO] )
						ipHi_offset += 1;
					ExtraLymanLines[ipISO][nelem][nExtraLyman].setHi(ipHi_offset);
					/* lower level is just ground state of the ion */
					ExtraLymanLines[ipISO][nelem][nExtraLyman].setLo(0);
					ExtraLymanLines[ipISO][nelem][nExtraLyman].AddLine2Stack();
					++nExtraLyman;
				}
				ASSERT(ExtraLymanLines[ipISO][nelem].size() == nExtraLyman);
			}
		}
	}

	// associate line and level stacks with species
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; ++nelem )
		{
			if( dense.lgElmtOn[nelem] )
			{
				long ion = nelem - ipISO;
				ASSERT( ion >= 0 && ion <= nelem );
				char chLabel[6] = {'\0'}, chTemp[4] = {'\0'};
				sprintf( chLabel, "%s", elementnames.chElementSym[nelem] );
				if( chLabel[1]==' ' )
					chLabel[1] = '\0';
				else
					chLabel[2] = '\0';
				if( ion==1 )
					sprintf( chTemp, "+" );
				else if( ion>0 )
					sprintf( chTemp, "+%li", ion );
				strcat( chLabel, chTemp );

				molecule *spmole = findspecies(chLabel);
				ASSERT( spmole != null_mole );
				mole.species[ spmole->index ].levels = &iso_sp[ipISO][nelem].st;
				mole.species[ spmole->index ].lines = &Transitions[ipISO][nelem];
			}
		}
	}

	return;
}

STATIC void iso_assign_quantum_numbers(void)
{
	long int
	  ipLo,
	  level,
	  i,
	  in,
	  il,
	  is,
	  ij;

	DEBUG_ENTRY( "iso_assign_quantum_numbers()" );

	for( long nelem=ipHYDROGEN; nelem < LIMELM; nelem++ )
	{
		long ipISO = ipH_LIKE;
		/* only check elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* 2 for doublet */
			is = ipDOUBLET;

			/* this loop is over quantum number n */
			for( in = 1L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					iso_sp[ipISO][nelem].st[i].n() = in;
					iso_sp[ipISO][nelem].st[i].S() = is;
					iso_sp[ipISO][nelem].st[i].l() = il;
					iso_sp[ipISO][nelem].st[i].j() = -1;
					iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = i;
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso_sp[ipISO][nelem].n_HighestResolved_max + 1;
			for( level = i; level< iso_sp[ipISO][nelem].numLevels_max; ++level)
			{
				iso_sp[ipISO][nelem].st[level].n() = in;
				iso_sp[ipISO][nelem].st[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].j() = -1;
				/* Point every l to same index for collapsed levels.	*/
				for( il = 0; il < in; ++il )
				{
					iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = level;
				}
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso_sp[ipISO][nelem].numLevels_max );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1) ) );

			/* Verify states and QuantumNumbers2Index agree in all cases	*/
			for( in = 2L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].n() == in );
					if( in <= iso_sp[ipISO][nelem].n_HighestResolved_max )
					{
						/* Must only check these for resolved levels...
						 * collapsed levels have pointers for l and s that will blow if used.	*/
						ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].l() == il );
						ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].S() == is );
					}
				}
			}
		}
	}

	/* then do he-like */
	for( long nelem=ipHELIUM; nelem < LIMELM; nelem++ )
	{
		long ipISO = ipHE_LIKE;
		/* only check elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem] )
		{
			i = 0;

			/* this loop is over quantum number n */
			for( in = 1L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					for( is = 3L; is >= 1L; is -= 2 )
					{
						/* All levels except singlet P follow the ordering scheme:	*/
						/*	lower l's have lower energy	*/
						/* 	triplets have lower energy	*/
						if( (il == 1L) && (is == 1L) )
							continue;
						/* n = 1 has no triplet, of course.	*/
						if( (in == 1L) && (is == 3L) )
							continue;

						/* singlets */
						if( is == 1 )
						{
							iso_sp[ipISO][nelem].st[i].n() = in;
							iso_sp[ipISO][nelem].st[i].S() = is;
							iso_sp[ipISO][nelem].st[i].l() = il;
							/* this is not a typo, J=L for singlets.  */
							iso_sp[ipISO][nelem].st[i].j() = il;
							iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = i;
							++i;
						}
						/* 2 triplet P is j-resolved */
						else if( (in == 2) && (il == 1) && (is == 3) )
						{
							ij = 0;
							do 
							{
								iso_sp[ipISO][nelem].st[i].n() = in;
								iso_sp[ipISO][nelem].st[i].S() = is;
								iso_sp[ipISO][nelem].st[i].l() = il;
								iso_sp[ipISO][nelem].st[i].j() = ij;
								iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = i;
								++i;
								++ij;
								/* repeat this for the separate j-levels within 2^3P. */
							}	while ( ij < 3 );
						}
						else
						{
							iso_sp[ipISO][nelem].st[i].n() = in;
							iso_sp[ipISO][nelem].st[i].S() = is;
							iso_sp[ipISO][nelem].st[i].l() = il;
							iso_sp[ipISO][nelem].st[i].j() = -1L;
							iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = i;
							++i;
						}
					}
				}
				/*	Insert singlet P at the end of every sequence for a given n.	*/
				if( in > 1L )
				{
					iso_sp[ipISO][nelem].st[i].n() = in;
					iso_sp[ipISO][nelem].st[i].S() = 1L;
					iso_sp[ipISO][nelem].st[i].l() = 1L;
					iso_sp[ipISO][nelem].st[i].j() = 1L;
					iso_sp[ipISO][nelem].QuantumNumbers2Index[in][1][1] = i;
					++i;
				}
			}
			/* now do the collapsed levels */
			in = iso_sp[ipISO][nelem].n_HighestResolved_max + 1;
			for( level = i; level< iso_sp[ipISO][nelem].numLevels_max; ++level)
			{
				iso_sp[ipISO][nelem].st[level].n() = in;
				iso_sp[ipISO][nelem].st[level].S() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].l() = -LONG_MAX;
				iso_sp[ipISO][nelem].st[level].j() = -1;
				/* Point every l and s to same index for collapsed levels.	*/
				for( il = 0; il < in; ++il )
				{
					for( is = 1; is <= 3; is += 2 )
					{
						iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] = level;
					}
				}
				++in;
			}
			--in;

			/* confirm that we did not overrun the array */
			ASSERT( i <= iso_sp[ipISO][nelem].numLevels_max );

			/* confirm that n is positive and not greater than the max n. */
			ASSERT( (in > 0) && (in < (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1) ) );

			/* Verify states and QuantumNumbers2Index agree in all cases	*/
			for( in = 2L; in <= iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max; ++in )
			{
				for( il = 0L; il < in; ++il )
				{
					for( is = 3L; is >= 1; is -= 2 )
					{
						/* Ground state is not triplicate.	*/
						if( (in == 1L) && (is == 3L) )
							continue;

						ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].n() == in );
						if( in <= iso_sp[ipISO][nelem].n_HighestResolved_max )
						{
							/* Must only check these for resolved levels...
							 * collapsed levels have pointers for l and s that will blow if used.	*/
							ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].l() == il );
							ASSERT( iso_sp[ipISO][nelem].st[ iso_sp[ipISO][nelem].QuantumNumbers2Index[in][il][is] ].S() == is );
						}
					}
				}
			}
		}
	}

	for( long ipISO=ipH_LIKE; ipISO<NISO; ipISO++ )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* must always do helium even if turned off */
			if( nelem < 2 || dense.lgElmtOn[nelem] )
			{
				for( ipLo=ipH1s; ipLo < iso_sp[ipISO][nelem].numLevels_max; ipLo++ )
				{
					iso_sp[ipISO][nelem].st[ipLo].nelem() = (int)(nelem+1);
					iso_sp[ipISO][nelem].st[ipLo].IonStg() = (int)(nelem+1-ipISO);

					if( iso_sp[ipISO][nelem].st[ipLo].j() >= 0 )
					{
						iso_sp[ipISO][nelem].st[ipLo].g() = 2.f*iso_sp[ipISO][nelem].st[ipLo].j()+1.f;
					}
					else if( iso_sp[ipISO][nelem].st[ipLo].l() >= 0 )
					{
						iso_sp[ipISO][nelem].st[ipLo].g() = (2.f*iso_sp[ipISO][nelem].st[ipLo].l()+1.f) *
							iso_sp[ipISO][nelem].st[ipLo].S();
					}
					else
					{
						if( ipISO == ipH_LIKE )
							iso_sp[ipISO][nelem].st[ipLo].g() = 2.f*(realnum)POW2( iso_sp[ipISO][nelem].st[ipLo].n() );
						else if( ipISO == ipHE_LIKE )
							iso_sp[ipISO][nelem].st[ipLo].g() = 4.f*(realnum)POW2( iso_sp[ipISO][nelem].st[ipLo].n() );
						else
						{
							/* replace this with correct thing if more sequences are added. */
							TotalInsanity();
						}
					}
					char chConfiguration[11] = "          ";
					long nCharactersWritten = 0;

					ASSERT( iso_sp[ipISO][nelem].st[ipLo].n() < 1000 );

					/* include j only if defined.  */
					if( iso_sp[ipISO][nelem].st[ipLo].n() > iso_sp[ipISO][nelem].n_HighestResolved_max )
					{
						nCharactersWritten = sprintf( chConfiguration, "n=%3li", 
							iso_sp[ipISO][nelem].st[ipLo].n() );
					}
					else if( iso_sp[ipISO][nelem].st[ipLo].j() > 0 )
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c_%li", 
							iso_sp[ipISO][nelem].st[ipLo].n(), 
							iso_sp[ipISO][nelem].st[ipLo].S(),
							chL[ MIN2( 20, iso_sp[ipISO][nelem].st[ipLo].l() ) ],
							iso_sp[ipISO][nelem].st[ipLo].j() );
					}
					else
					{
						nCharactersWritten = sprintf( chConfiguration, "%3li^%li%c", 
							iso_sp[ipISO][nelem].st[ipLo].n(), 
							iso_sp[ipISO][nelem].st[ipLo].S(),
							chL[ MIN2( 20, iso_sp[ipISO][nelem].st[ipLo].l()) ] );
					}

					ASSERT( nCharactersWritten <= 10 );
					chConfiguration[10] = '\0';

					strncpy( iso_sp[ipISO][nelem].st[ipLo].chConfig(), chConfiguration, 10 );
				}
			}
		}
	}
	return;
}

#if defined(__ICC) && defined(__i386)
#pragma optimization_level 1
#endif
STATIC void FillExtraLymanLine( const TransitionList::iterator& t, long ipISO, long nelem, long nHi )
{
	double Enerwn, Aul;

	DEBUG_ENTRY( "FillExtraLymanLine()" );

	/* atomic number or charge and stage: */
	(*(*t).Hi()).nelem() = (int)(nelem+1);
	(*(*t).Hi()).IonStg() = (int)(nelem+1-ipISO);
	
	(*(*t).Hi()).n() = nHi;

	/* statistical weight is same as statistical weight of corresponding LyA. */
	(*(*t).Hi()).g() = iso_sp[ipISO][nelem].st[iso_ctrl.nLyaLevel[ipISO]].g();

	/* energies */
	Enerwn = iso_sp[ipISO][nelem].fb[0].xIsoLevNIonRyd * RYD_INF * (  1. - 1./POW2((double)nHi) );

	/* transition energy in various units:*/
	(*t).EnergyWN() = (realnum)(Enerwn);
	(*t).WLAng() = (realnum)(1.0e8/ Enerwn/ RefIndex(Enerwn));
	(*(*t).Hi()).energy().set( Enerwn, "cm^-1" );

	if( ipISO == ipH_LIKE )
	{
		Aul = H_Einstein_A( nHi, 1, 1, 0, nelem+1 );
	}
	else
	{
		if( nelem == ipHELIUM )
		{
			/* A simple fit for the calculation of Helium lyman Aul's.	*/
			Aul = (1.508e10) / pow((double)nHi,2.975);
		}
		else 
		{
			/* Fit to values given in 
			 * >>refer	He-like	As	Johnson, W.R., Savukov, I.M., Safronova, U.I., & 
			 * >>refercon	Dalgarno, A., 2002, ApJS 141, 543J	*/
			/* originally astro.ph. 0201454  */
			Aul = 1.375E10 * pow((double)nelem, 3.9) / pow((double)nHi,3.1);
		}
	}

	(*t).Emis().Aul() = (realnum)Aul;

	(*(*t).Hi()).lifetime() = iso_state_lifetime( ipISO, nelem, nHi, 1 );

	(*t).Emis().dampXvel() = (realnum)( 1.f / (*(*t).Hi()).lifetime() / PI4 / (*t).EnergyWN() );

	(*t).Emis().iRedisFun() = iso_ctrl.ipResoRedist[ipISO];

	(*t).Emis().gf() = (realnum)(GetGF((*t).Emis().Aul(),	(*t).EnergyWN(), (*(*t).Hi()).g()));

	/* derive the abs coef, call to function is Emis().gf(), wl (A), g_low */
	(*t).Emis().opacity() = (realnum)(abscf((*t).Emis().gf(), (*t).EnergyWN(), (*(*t).Lo()).g()));

	/* create array indices that will blow up */
	(*t).ipCont() = INT_MIN;
	(*t).Emis().ipFine() = INT_MIN;

	{
		/* option to print particulars of some line when called
			* a prettier print statement is near where chSpin is defined below
			* search for "pretty print" */
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC )
		{
			fprintf(ioQQQ,"%li\t%li\t%.2e\t%.2e\n",
				nelem+1,
				nHi,
				(*t).Emis().Aul() , 	
				(*t).Emis().opacity()
				);
		}
	}
	return;
}

/* calculate radiative lifetime of an individual iso state */
double iso_state_lifetime( long ipISO, long nelem, long n, long l )
{
	/* >>refer hydro lifetimes	Horbatsch, M. W., Horbatsch, M. and Hessels, E. A. 2005, JPhysB, 38, 1765 */

	double tau, t0, eps2;
	/* mass of electron */
	double m = ELECTRON_MASS;
	/* nuclear mass */
	double M = (double)dense.AtomicWeight[nelem] * ATOMIC_MASS_UNIT;
	double mu = (m*M)/(M+m);
	long z = 1;
	long Z = nelem + 1 - ipISO;
	
	DEBUG_ENTRY( "iso_state_lifetime()" );

	/* this should not be used for l=0 per the Horbatsch et al. paper */
	ASSERT( l > 0 );

	eps2 = 1. - ( l*l + l + 8./47. - (l+1.)/69./n ) / POW2( (double)n );

	t0 = 3. * H_BAR * pow( (double)n, 5.) / 
		( 2. * POW4( (double)( z * Z ) ) * pow( FINE_STRUCTURE, 5. ) * mu * POW2( SPEEDLIGHT ) ) *
		POW2( (m + M)/(Z*m + z*M) );

	tau = t0 * ( 1. - eps2 ) / 
		( 1. + 19./88.*( (1./eps2 - 1.) * log( 1. - eps2 ) + 1. - 
		0.5 * eps2 - 0.025 * eps2 * eps2 ) );

	if( ipISO == ipHE_LIKE )
	{
		/* iso_state_lifetime is not spin specific, must exclude helike triplet here. */
		tau /= 3.;
		/* this is also necessary to correct the helike lifetimes */
		tau *= 1.1722 * pow( (double)nelem, 0.1 );
	}

	/* would probably need a new lifetime algorithm for any other iso sequences. */
	ASSERT( ipISO <= ipHE_LIKE );
	ASSERT( tau > 0. );

	return tau;
}

/* calculate cascade probabilities, branching ratios, and associated errors. */
void iso_cascade( long ipISO, long nelem )
{
	/* The sum of all A's coming out of a given n,
	 * Below we assert a monotonic trend. */
	double *SumAPerN;

	long int i, j, ipLo, ipHi;

	DEBUG_ENTRY( "iso_cascade()" );

	SumAPerN = ((double*)MALLOC( (size_t)(iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1 )*sizeof(double )));
	memset( SumAPerN, 0, (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max + 1 )*sizeof(double ) );

	/* Initialize some ground state stuff, easier here than in loops.	*/
	iso_sp[ipISO][nelem].CascadeProb[0][0] = 1.;
	if( iso_ctrl.lgRandErrGen[ipISO] )
	{
		iso_sp[ipISO][nelem].fb[0].SigmaAtot = 0.;
		iso_sp[ipISO][nelem].ex[0][0].SigmaCascadeProb = 0.;
	}

	/***************************************************************************/
	/****** Cascade probabilities, Branching ratios, and associated errors *****/
	/***************************************************************************/
	for( ipHi=1; ipHi<iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		double SumAs = 0.;

		/** Cascade probabilities are as defined in Robbins 68,
		 * generalized here for cascade probability for any iso sequence.	
		 * >>refer He triplets	Robbins, R.R. 1968, ApJ 151, 497R	
		 * >>refer He triplets	Robbins, R.R. 1968a, ApJ 151, 511R	*/

		/* initialize variables. */
		iso_sp[ipISO][nelem].CascadeProb[ipHi][ipHi] = 1.;
		iso_sp[ipISO][nelem].CascadeProb[ipHi][0] = 0.;
		iso_sp[ipISO][nelem].BranchRatio[ipHi][0] = 0.;

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot = 0.;
			iso_sp[ipISO][nelem].ex[ipHi][ipHi].SigmaCascadeProb = 0.;
		}

		long ipLoStart = 0;
		if( opac.lgCaseB && L_(ipHi)==1 && (ipISO==ipH_LIKE || S_(ipHi)==1) )
			ipLoStart = 1;

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			SumAs += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();
		}

		for( ipLo=ipLoStart; ipLo<ipHi; ipLo++ )
		{
			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
			{
				iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] = 0.;
				iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] = 0.;
				continue;
			}

			iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] = 0.;
			iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] = 
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() / SumAs;

			ASSERT( iso_sp[ipISO][nelem].BranchRatio[ipHi][ipLo] <= 1.0000001 );

			SumAPerN[N_(ipHi)] += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();

			/* there are some negative energy transitions, where the order
			 * has changed, but these are not optically allowed, these are
			 * same n, different L, forbidden transitions */
			ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() > 0. ||
				iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA );

			if( iso_ctrl.lgRandErrGen[ipISO] )
			{
				ASSERT( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[IPRAD] >= 0. );
				/* Uncertainties in A's are added in quadrature, square root is taken below. */
				iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot += 
					pow( iso_sp[ipISO][nelem].ex[ipHi][ipLo].Error[IPRAD] * 
					(double)iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul(), 2. );
			}
		}

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			/* Uncertainties in A's are added in quadrature above, square root taken here. */
			iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot = sqrt( iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot );
		}

		/* cascade probabilities */
		for( i=0; i<ipHi; i++ )
		{
			for( ipLo=0; ipLo<=i; ipLo++ )
			{
				iso_sp[ipISO][nelem].CascadeProb[ipHi][ipLo] += iso_sp[ipISO][nelem].BranchRatio[ipHi][i] * iso_sp[ipISO][nelem].CascadeProb[i][ipLo];
			}
		}

		if( iso_ctrl.lgRandErrGen[ipISO] )
		{
			for( ipLo=0; ipLo<ipHi; ipLo++ )
			{
				double SigmaCul = 0.;
				for( i=ipLo; i<ipHi; i++ )
				{
					if( iso_sp[ipISO][nelem].trans(ipHi,i).Emis().Aul() > iso_ctrl.SmallA )
					{
						/* Uncertainties in A's and cascade probabilities */
						double SigmaA = iso_sp[ipISO][nelem].ex[ipHi][i].Error[IPRAD] * 
							iso_sp[ipISO][nelem].trans(ipHi,i).Emis().Aul();
						SigmaCul += 
							pow(SigmaA*iso_sp[ipISO][nelem].CascadeProb[i][ipLo]*iso_sp[ipISO][nelem].st[ipHi].lifetime(), 2.) +
							pow(iso_sp[ipISO][nelem].fb[ipHi].SigmaAtot*iso_sp[ipISO][nelem].BranchRatio[ipHi][i]*
							iso_sp[ipISO][nelem].CascadeProb[i][ipLo]*iso_sp[ipISO][nelem].st[ipHi].lifetime(), 2.) +
							pow(iso_sp[ipISO][nelem].ex[i][ipLo].SigmaCascadeProb*iso_sp[ipISO][nelem].BranchRatio[ipHi][i], 2.);
					}
				}
				SigmaCul = sqrt(SigmaCul);
				iso_sp[ipISO][nelem].ex[ipHi][ipLo].SigmaCascadeProb = SigmaCul;
			}
		}
	}

	/************************************************************************/
	/*** Allowed decay conversion probabilities. See Robbins68b, Table 1. ***/
	/************************************************************************/
	{
		enum {DEBUG_LOC=false};

		if( DEBUG_LOC && (nelem == ipHELIUM) && (ipISO==ipHE_LIKE) )
		{
			/* To output Bm(n,l; ipLo), set ipLo, hi_l, and hi_s accordingly.	*/
			long int hi_l,hi_s;
			double Bm;

			/* these must be set for following output to make sense
			 * as is, a dangerous bit of code - set NaN for safety */
			hi_s = -100000;
			hi_l = -100000;
			ipLo = -100000;
			/* tripS to 2^3P	*/
			//hi_l = 0, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripD to 2^3P	*/
			//hi_l = 2, hi_s = 3, ipLo = ipHe2p3P0;

			/* tripP to 2^3S	*/
			//hi_l = 1, hi_s = 3, ipLo = ipHe2s3S;	

			ASSERT( hi_l != iso_sp[ipISO][nelem].st[ipLo].l() );

			fprintf(ioQQQ,"Bm(n,%ld,%ld;%ld)\n",hi_l,hi_s,ipLo);
			fprintf(ioQQQ,"m\t2\t\t3\t\t4\t\t5\t\t6\n");

			for( ipHi=ipHe2p3P2; ipHi<iso_sp[ipISO][nelem].numLevels_max-iso_sp[ipISO][nelem].nCollapsed_max; ipHi++ )
			{
				/* Pick out excitations from metastable 2tripS to ntripP.	*/
				if( (iso_sp[ipISO][nelem].st[ipHi].l() == 1) && (iso_sp[ipISO][nelem].st[ipHi].S() == 3) )
				{
					fprintf(ioQQQ,"\n%ld\t",iso_sp[ipISO][nelem].st[ipHi].n());
					j = 0;
					Bm = 0;
					for( i = ipLo; i<=ipHi; i++)
					{
						if( (iso_sp[ipISO][nelem].st[i].l() == hi_l) && (iso_sp[ipISO][nelem].st[i].S() == hi_s)  )
						{
							if( (ipLo == ipHe2p3P0) && (i > ipHe2p3P2) )
							{
								Bm += iso_sp[ipISO][nelem].CascadeProb[ipHi][i] * ( iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P0] + 
									iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P1] + iso_sp[ipISO][nelem].BranchRatio[i][ipHe2p3P2] );
							}
							else
								Bm += iso_sp[ipISO][nelem].CascadeProb[ipHi][i] * iso_sp[ipISO][nelem].BranchRatio[i][ipLo];

							if( (i == ipHe2p3P0) || (i == ipHe2p3P1) || (i == ipHe2p3P2) )
							{
								j++;
								if(j == 3)
								{
									fprintf(ioQQQ,"%2.4e\t",Bm);
									Bm = 0;
								}
							}
							else
							{
								fprintf(ioQQQ,"%2.4e\t",Bm);
								Bm = 0;
							}
						}
					}
				}
			}
			fprintf(ioQQQ,"\n\n");
		}
	}

	/******************************************************/
	/***  Lifetimes should increase monotonically with  ***/
	/***  increasing n...Make sure the A's decrease.    ***/
	/******************************************************/
	for( i=2; i < iso_sp[ipISO][nelem].n_HighestResolved_max; ++i)
	{
		ASSERT( (SumAPerN[i] > SumAPerN[i+1]) || opac.lgCaseB );
	}

	{
		enum {DEBUG_LOC=false};
		if( DEBUG_LOC /* && (ipISO == ipH_LIKE) && (nelem == ipHYDROGEN) */)
		{
			for( i = 2; i<= (iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max); ++i)
			{
				fprintf(ioQQQ,"n %ld\t lifetime %.4e\n", i, 1./SumAPerN[i]);
			}
		}
	}

	free( SumAPerN );

	return;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void iso_satellite( void )
{
	DEBUG_ENTRY( "iso_satellite()" );

	for( long ipISO = ipHE_LIKE; ipISO < NISO; ipISO++ )
	{
		for( long nelem = ipISO; nelem < LIMELM; nelem++ )
		{
			if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
			{
				for( long i=0; i<iso_sp[ipISO][nelem].numLevels_max; i++ )
				{
					char chLab[5]="    ";

					TransitionList::iterator tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][i];
					(*tr).Zero();
			
					/* Make approximation that all levels have energy of H-like 2s level */
					/* Lines to 1s2s have roughly energy of parent Ly-alpha.  So lines to 1snL will have an energy
					 * smaller by the difference between nL and 2s energies.  Therefore, the following has
					 * energy of parent Ly-alpha MINUS the difference between daughter level and daughter n=2 level. */
					(*tr).WLAng() = (realnum)(RYDLAM/
						((iso_sp[ipISO-1][nelem].fb[0].xIsoLevNIonRyd - iso_sp[ipISO-1][nelem].fb[1].xIsoLevNIonRyd) -
						 (iso_sp[ipISO][nelem].fb[1].xIsoLevNIonRyd- iso_sp[ipISO][nelem].fb[i].xIsoLevNIonRyd)) );

					(*tr).EnergyWN() = 1.e8f / 
						(*tr).WLAng();

					/* generate label for this ion */
					sprintf( chLab, "%2s%2ld",elementnames.chElementSym[nelem], nelem+1-ipISO );

					(*tr).Emis().iRedisFun() = ipCRDW;
					/* this is not the usual nelem, is it atomic not C scale. */
					(*(*tr).Hi()).nelem() = nelem + 1;
					(*(*tr).Hi()).IonStg() = nelem + 1 - ipISO;
					fixit(); /* what should the stat weight of the upper level be? For now say 2. */
					(*(*tr).Hi()).g() = 2.f;
					// The lower level must already be initialized.  
					ASSERT( (*(*tr).Lo()).g() == iso_sp[ipISO][nelem].st[i].g() );
					//(*(*tr).Lo()).g() = iso_sp[ipISO][nelem].st[i].g();
					(*tr).Emis().PopOpc() = 
						(*(*tr).Lo()).Pop();

					(*tr).Emis().pump() = 0.;

				}
			}
		}
	}

	return;
}

void iso_satellite_update( long nelem )
{
	double ConBoltz, LTE_pop=SMALLFLOAT+FLT_EPSILON, factor1, ConvLTEPOP;
	
	DEBUG_ENTRY( "iso_satellite_update()" );

	for( long ipISO = ipHE_LIKE; ipISO < MIN2(NISO,nelem+1); ipISO++ )
	{
		if( dense.lgElmtOn[nelem] && iso_ctrl.lgDielRecom[ipISO] ) 
		{
			for( long i=0; i<iso_sp[ipISO][nelem].numLevels_max; i++ )
			{
				double dr_rate = iso_sp[ipISO][nelem].fb[i].DielecRecomb * iso_ctrl.lgDielRecom[ipISO];
				
				TransitionList::iterator tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][i];
				(*tr).Emis().phots() = 
					dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO];
				
				(*tr).Emis().xIntensity() = 
					(*tr).Emis().phots() * 
					ERG1CM * (*tr).EnergyWN();
				
				/* We set line intensity above using a rate, but here we need a transition probability.  
				 * We can obtain this by dividing dr_rate by the population of the autoionizing level.
				 * We assume this level is in statistical equilibrium. */
				factor1 = HION_LTE_POP*dense.AtomicWeight[nelem]/
					(dense.AtomicWeight[nelem]+ELECTRON_MASS/ATOMIC_MASS_UNIT);
				
				/* term in () is stat weight of electron * ion */
				ConvLTEPOP = pow(factor1,1.5)/(2.*iso_ctrl.stat_ion[ipISO])/phycon.te32;
				
				/* This Boltzmann factor is exp( +ioniz energy / Te ).  For simplicity, we make
				 * the fair approximation that all of the autoionizing levels have an energy
				 * equal to the parents n=2. */
				ConBoltz = dsexp(iso_sp[ipISO-1][nelem].fb[1].xIsoLevNIonRyd/phycon.te_ryd);
				
				if( ConBoltz >= SMALLDOUBLE )
				{
					/* The energy used to calculate ConBoltz above
					 * should be negative since this is above the continuum, but 
					 * to be safe we calculate ConBoltz with a positive energy above
					 * and multiply by it here instead of dividing.  */
					LTE_pop = (*(*tr).Hi()).g() * ConBoltz * ConvLTEPOP;
				}
				
				LTE_pop = max( LTE_pop, 1e-30f );
				
				/* Now the transition probability is simply dr_rate/LTE_pop. */
				(*tr).Emis().Aul() = (realnum)(dr_rate/LTE_pop);
				(*tr).Emis().Aul() = 
					max( iso_ctrl.SmallA, (*tr).Emis().Aul() );
				
				(*tr).Emis().gf() = (realnum)GetGF(
					(*tr).Emis().Aul(),
					(*tr).EnergyWN(),
					(*(*tr).Hi()).g());
				
				(*tr).Emis().gf() = 
					max( 1e-20f, (*tr).Emis().gf() );
				
				(*(*tr).Hi()).Pop() = LTE_pop * dense.xIonDense[nelem][nelem+1-ipISO] * dense.eden;
				
				(*tr).Emis().PopOpc() = 
					(*(*tr).Lo()).Pop() - 
					(*(*tr).Hi()).Pop() * 
					(*(*tr).Lo()).g()/(*(*tr).Hi()).g();
				
				(*tr).Emis().opacity() = 
					(realnum)(abscf((*tr).Emis().gf(),
										 (*tr).EnergyWN(),
										 (*(*tr).Lo()).g()));
				
				/* a typical transition probability is of order 1e10 s-1 */
				double lifetime = 1e-10;
				
				(*tr).Emis().dampXvel() = (realnum)( 
					(1.f/lifetime)/PI4/(*tr).EnergyWN());
			}
		}
	}
	
	return;
}

long iso_get_total_num_levels( long ipISO, long nmaxResolved, long numCollapsed )
{
	DEBUG_ENTRY( "iso_get_total_num_levels()" );

	long tot_num_levels;

	/* return the number of levels up to and including nmaxResolved PLUS 
	 * the number (numCollapsed) of collapsed n-levels */		

	if( ipISO == ipH_LIKE )
	{
		tot_num_levels = (long)( nmaxResolved * 0.5 *( nmaxResolved + 1 ) ) + numCollapsed;
	}
	else if( ipISO == ipHE_LIKE )
	{
		tot_num_levels = nmaxResolved*nmaxResolved + nmaxResolved + 1 + numCollapsed;
	}
	else
		TotalInsanity();

	return tot_num_levels;
}

void iso_update_num_levels( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_update_num_levels()" );

	/* This is the minimum resolved nmax. */
	ASSERT( iso_sp[ipISO][nelem].n_HighestResolved_max >= 3 );

	iso_sp[ipISO][nelem].numLevels_max = 
		iso_get_total_num_levels( ipISO, iso_sp[ipISO][nelem].n_HighestResolved_max, iso_sp[ipISO][nelem].nCollapsed_max );

	if( iso_sp[ipISO][nelem].numLevels_max > iso_sp[ipISO][nelem].numLevels_malloc )
	{
		fprintf( ioQQQ, "The number of levels for ipISO %li, nelem %li, has been increased since the initial coreload.\n",
			ipISO, nelem );
		fprintf( ioQQQ, "This cannot be done.\n" );
		cdEXIT(EXIT_FAILURE);
	}

	/* set local copies to the max values */
	iso_sp[ipISO][nelem].numLevels_local = iso_sp[ipISO][nelem].numLevels_max;
	iso_sp[ipISO][nelem].nCollapsed_local = iso_sp[ipISO][nelem].nCollapsed_max;
	iso_sp[ipISO][nelem].n_HighestResolved_local = iso_sp[ipISO][nelem].n_HighestResolved_max;

	/* find the largest number of levels in any element in all iso sequences
	 * we will allocate one matrix for ionization solver, and just use a piece of that memory
	 * for smaller models. */
	max_num_levels = MAX2( max_num_levels, iso_sp[ipISO][nelem].numLevels_max);

	return;
}

void iso_collapsed_bnl_set( long ipISO, long nelem )
{

	DEBUG_ENTRY( "iso_collapsed_bnl_set()" );

	double bnl_array[4][3][4][10] = {
		{
			/* H */
			{
				{6.13E-01,	2.56E-01,	1.51E-01,	2.74E-01,	3.98E-01,	4.98E-01,	5.71E-01,	6.33E-01,	7.28E-01,	9.59E-01},
				{1.31E+00,	5.17E-01,	2.76E-01,	4.47E-01,	5.87E-01,	6.82E-01,	7.44E-01,	8.05E-01,	9.30E-01,	1.27E+00},
				{1.94E+00,	7.32E-01,	3.63E-01,	5.48E-01,	6.83E-01,	7.66E-01,	8.19E-01,	8.80E-01,	1.02E+00,	1.43E+00},
				{2.53E+00,	9.15E-01,	4.28E-01,	6.16E-01,	7.42E-01,	8.13E-01,	8.60E-01,	9.22E-01,	1.08E+00,	1.56E+00}
			},
			{
				{5.63E-01,	2.65E-01,	1.55E-01,	2.76E-01,	3.91E-01,	4.75E-01,	5.24E-01,	5.45E-01,	5.51E-01,	5.53E-01},
				{1.21E+00,	5.30E-01,	2.81E-01,	4.48E-01,	5.80E-01,	6.62E-01,	7.05E-01,	7.24E-01,	7.36E-01,	7.46E-01},
				{1.81E+00,	7.46E-01,	3.68E-01,	5.49E-01,	6.78E-01,	7.51E-01,	7.88E-01,	8.09E-01,	8.26E-01,	8.43E-01},
				{2.38E+00,	9.27E-01,	4.33E-01,	6.17E-01,	7.38E-01,	8.05E-01,	8.40E-01,	8.65E-01,	8.92E-01,	9.22E-01}
			},
			{
				{2.97E-01,	2.76E-01,	2.41E-01,	3.04E-01,	3.66E-01,	4.10E-01,	4.35E-01,	4.48E-01,	4.52E-01,	4.53E-01},
				{5.63E-01,	5.04E-01,	3.92E-01,	4.67E-01,	5.39E-01,	5.85E-01,	6.10E-01,	6.20E-01,	6.23E-01,	6.23E-01},
				{7.93E-01,	6.90E-01,	4.94E-01,	5.65E-01,	6.36E-01,	6.79E-01,	7.00E-01,	7.09E-01,	7.11E-01,	7.11E-01},
				{1.04E+00,	8.66E-01,	5.62E-01,	6.31E-01,	7.01E-01,	7.43E-01,	7.63E-01,	7.71E-01,	7.73E-01,	7.73E-01}
			}
		},
		{
			/* He+ */
			{
				{6.70E-02,	2.93E-02,	1.94E-02,	4.20E-02,	7.40E-02,	1.12E-01,	1.51E-01,	1.86E-01,	2.26E-01,	3.84E-01},
				{2.39E-01,	1.03E-01,	6.52E-02,	1.31E-01,	2.11E-01,	2.91E-01,	3.61E-01,	4.17E-01,	4.85E-01,	8.00E-01},
				{4.26E-01,	1.80E-01,	1.10E-01,	2.09E-01,	3.18E-01,	4.15E-01,	4.93E-01,	5.54E-01,	6.34E-01,	1.04E+00},
				{6.11E-01,	2.55E-01,	1.51E-01,	2.74E-01,	3.99E-01,	5.02E-01,	5.80E-01,	6.41E-01,	7.30E-01,	1.21E+00}
			},
			{
				{6.79E-02,	3.00E-02,	2.00E-02,	4.30E-02,	7.48E-02,	1.11E-01,	1.44E-01,	1.70E-01,	1.87E-01,	1.96E-01},
				{2.40E-01,	1.04E-01,	6.62E-02,	1.32E-01,	2.11E-01,	2.87E-01,	3.51E-01,	3.98E-01,	4.32E-01,	4.58E-01},
				{4.26E-01,	1.81E-01,	1.11E-01,	2.10E-01,	3.17E-01,	4.12E-01,	4.89E-01,	5.53E-01,	6.14E-01,	6.84E-01},
				{6.12E-01,	2.55E-01,	1.51E-01,	2.73E-01,	3.97E-01,	4.98E-01,	5.77E-01,	6.51E-01,	7.82E-01,	1.18E+00}
			},
			{
				{4.98E-02,	3.47E-02,	2.31E-02,	4.54E-02,	7.14E-02,	9.37E-02,	1.08E-01,	1.13E-01,	1.13E-01,	1.11E-01},
				{1.75E-01,	1.16E-01,	7.36E-02,	1.36E-01,	2.01E-01,	2.50E-01,	2.76E-01,	2.84E-01,	2.81E-01,	2.77E-01},
				{3.38E-01,	1.97E-01,	1.18E-01,	2.13E-01,	3.06E-01,	3.72E-01,	4.06E-01,	4.15E-01,	4.10E-01,	4.04E-01},
				{6.01E-01,	2.60E-01,	1.53E-01,	2.76E-01,	3.95E-01,	4.87E-01,	5.45E-01,	5.76E-01,	5.93E-01,	6.05E-01}
			}
		},
		{
			/* He singlets */
			{
				{1.77E-01,	3.59E-01,	1.54E-01,	2.75E-01,	3.98E-01,	4.94E-01,	5.51E-01,	5.68E-01,	5.46E-01,	4.97E-01},
				{4.09E-01,	7.23E-01,	2.83E-01,	4.48E-01,	5.89E-01,	6.78E-01,	7.22E-01,	7.30E-01,	7.07E-01,	6.65E-01},
				{6.40E-01,	1.02E+00,	3.74E-01,	5.49E-01,	6.85E-01,	7.63E-01,	7.98E-01,	8.03E-01,	7.84E-01,	7.53E-01},
				{8.70E-01,	1.28E+00,	4.42E-01,	6.17E-01,	7.44E-01,	8.13E-01,	8.42E-01,	8.46E-01,	8.34E-01,	8.13E-01}
			},
			{
				{1.78E-01,	3.62E-01,	1.55E-01,	2.73E-01,	3.91E-01,	4.73E-01,	5.10E-01,	5.04E-01,	4.70E-01,	4.32E-01},
				{4.08E-01,	7.26E-01,	2.83E-01,	4.45E-01,	5.79E-01,	6.54E-01,	6.78E-01,	6.64E-01,	6.30E-01,	5.98E-01},
				{6.37E-01,	1.03E+00,	3.73E-01,	5.46E-01,	6.75E-01,	7.40E-01,	7.57E-01,	7.43E-01,	7.15E-01,	6.92E-01},
				{8.65E-01,	1.28E+00,	4.41E-01,	6.14E-01,	7.35E-01,	7.92E-01,	8.05E-01,	7.95E-01,	7.74E-01,	7.59E-01}
			},
			{
				{2.07E-01,	3.73E-01,	1.73E-01,	2.85E-01,	4.03E-01,	4.76E-01,	5.06E-01,	5.03E-01,	4.84E-01,	4.63E-01},
				{4.32E-01,	7.13E-01,	3.06E-01,	4.54E-01,	5.81E-01,	6.44E-01,	6.59E-01,	6.49E-01,	6.28E-01,	6.11E-01},
				{6.40E-01,	9.85E-01,	3.98E-01,	5.53E-01,	6.74E-01,	7.27E-01,	7.36E-01,	7.26E-01,	7.10E-01,	6.98E-01},
				{8.38E-01,	1.21E+00,	4.67E-01,	6.20E-01,	7.34E-01,	7.79E-01,	7.87E-01,	7.79E-01,	7.69E-01,	7.63E-01}
			}
		},
		{
			/* He triplets */
			{
				{9.31E-02,	3.96E-01,	1.36E-01,	2.74E-01,	3.99E-01,	4.95E-01,	5.52E-01,	5.70E-01,	5.48E-01,	4.96E-01},
				{2.25E-01,	8.46E-01,	2.49E-01,	4.46E-01,	5.89E-01,	6.79E-01,	7.23E-01,	7.31E-01,	7.08E-01,	6.64E-01},
				{3.59E-01,	1.24E+00,	3.30E-01,	5.47E-01,	6.85E-01,	7.63E-01,	7.98E-01,	8.04E-01,	7.85E-01,	7.53E-01},
				{4.93E-01,	1.60E+00,	3.91E-01,	6.15E-01,	7.44E-01,	8.13E-01,	8.42E-01,	8.47E-01,	8.35E-01,	8.12E-01}
			},
			{
				{9.32E-02,	3.99E-01,	1.35E-01,	2.72E-01,	3.91E-01,	4.75E-01,	5.14E-01,	5.09E-01,	4.73E-01,	4.31E-01},
				{2.25E-01,	8.49E-01,	2.49E-01,	4.44E-01,	5.79E-01,	6.56E-01,	6.81E-01,	6.68E-01,	6.31E-01,	5.96E-01},
				{3.58E-01,	1.25E+00,	3.29E-01,	5.44E-01,	6.76E-01,	7.42E-01,	7.60E-01,	7.46E-01,	7.16E-01,	6.91E-01},
				{4.92E-01,	1.60E+00,	3.90E-01,	6.12E-01,	7.36E-01,	7.93E-01,	8.07E-01,	7.97E-01,	7.74E-01,	7.58E-01}
			},
			{
				{1.13E-01,	4.21E-01,	1.47E-01,	2.83E-01,	4.04E-01,	4.80E-01,	5.13E-01,	5.12E-01,	4.93E-01,	4.71E-01},
				{2.52E-01,	8.56E-01,	2.61E-01,	4.50E-01,	5.82E-01,	6.48E-01,	6.66E-01,	6.56E-01,	6.35E-01,	6.16E-01},
				{3.85E-01,	1.23E+00,	3.41E-01,	5.49E-01,	6.75E-01,	7.30E-01,	7.41E-01,	7.31E-01,	7.15E-01,	7.02E-01},
				{5.14E-01,	1.56E+00,	4.01E-01,	6.15E-01,	7.34E-01,	7.82E-01,	7.90E-01,	7.83E-01,	7.72E-01,	7.65E-01}
			}
		}
	};

	double temps[4] = {5000., 10000., 15000., 20000. };
	double log_dens[3] = {2., 4., 6.};
	long ipTe, ipDens;

	ASSERT( nelem <= 1 );

	/* find temperature in tabulated values.  */
	ipTe = hunt_bisect( temps, 4, phycon.te );			
	ipDens = hunt_bisect( log_dens, 3, log10(dense.eden) );			

	ASSERT( (ipTe >=0) && (ipTe < 3)  );
	ASSERT( (ipDens >=0) && (ipDens < 2)  );

	for( long nHi=iso_sp[ipISO][nelem].n_HighestResolved_max+1; nHi<=iso_sp[ipISO][nelem].n_HighestResolved_max+iso_sp[ipISO][nelem].nCollapsed_max; nHi++ )
	{
		for( long lHi=0; lHi<nHi; lHi++ )
		{
			for( long sHi=1; sHi<4; sHi++ )
			{
				if( ipISO == ipH_LIKE && sHi != 2 )
					continue;
				else if( ipISO == ipHE_LIKE && sHi != 1 && sHi != 3 )
					continue;

				double bnl_at_lo_den, bnl_at_hi_den, bnl;
				double bnl_max, bnl_min, temp, dens;

				long ipL = MIN2(9,lHi);
				long ip1;

				if( nelem==ipHYDROGEN )
					ip1 = 0;
				else if( nelem==ipHELIUM )
				{
					if( ipISO==ipH_LIKE )
						ip1 = 1;
					else if( ipISO==ipHE_LIKE )
					{
						if( sHi==1 )
							ip1 = 2;
						else if( sHi==3 )
							ip1 = 3;
						else
							TotalInsanity();
					}
					else
						TotalInsanity();
				}
				else
					TotalInsanity();

				temp = MAX2( temps[0], phycon.te );
				temp = MIN2( temps[3], temp );

				dens = MAX2( log_dens[0], log10(dense.eden) );
				dens = MIN2( log_dens[2], dens );

				/* Calculate the answer...must interpolate on two variables.
				 * First interpolate on T, at both the lower and upper densities.
				 * Then interpolate between these results for the right density.	*/
			
				if( temp < temps[0] && dens < log_dens[0] )
					bnl = bnl_array[ip1][0][0][ipL];
				else if( temp < temps[0] && dens >= log_dens[2] )
					bnl = bnl_array[ip1][2][0][ipL]; 
				else if( temp >= temps[3] && dens < log_dens[0] )
					bnl = bnl_array[ip1][0][3][ipL]; 
				else if( temp >= temps[3] && dens >= log_dens[2] )
					bnl = bnl_array[ip1][2][3][ipL];
				else
				{
					bnl_at_lo_den = ( temp - temps[ipTe]) / (temps[ipTe+1] - temps[ipTe]) *
						(bnl_array[ip1][ipDens][ipTe+1][ipL] - bnl_array[ip1][ipDens][ipTe][ipL]) + bnl_array[ip1][ipDens][ipTe][ipL];

					bnl_at_hi_den = ( temp - temps[ipTe]) / (temps[ipTe+1] - temps[ipTe]) *
						(bnl_array[ip1][ipDens+1][ipTe+1][ipL] - bnl_array[ip1][ipDens+1][ipTe][ipL]) + bnl_array[ip1][ipDens+1][ipTe][ipL];

					bnl = ( dens - log_dens[ipDens]) / (log_dens[ipDens+1] - log_dens[ipDens]) * 
						(bnl_at_hi_den - bnl_at_lo_den) + bnl_at_lo_den;
				}

				/** these are just sanity checks, the interpolated value should be between values at interpolation points */
				{
					bnl_max = MAX4( bnl_array[ip1][ipDens][ipTe+1][ipL], bnl_array[ip1][ipDens+1][ipTe+1][ipL],
						bnl_array[ip1][ipDens][ipTe][ipL], bnl_array[ip1][ipDens+1][ipTe][ipL] );
					ASSERT( bnl <= bnl_max );

					bnl_min = MIN4( bnl_array[ip1][ipDens][ipTe+1][ipL], bnl_array[ip1][ipDens+1][ipTe+1][ipL],
						bnl_array[ip1][ipDens][ipTe][ipL], bnl_array[ip1][ipDens+1][ipTe][ipL] );
					ASSERT( bnl >= bnl_min );
				}

				iso_sp[ipISO][nelem].bnl_effective[nHi][lHi][sHi] = bnl;

				ASSERT( iso_sp[ipISO][nelem].bnl_effective[nHi][lHi][sHi]  > 0. );
			}
		}
	}

	return;
}


void iso_collapsed_bnl_print( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_bnl_print()" );

	for( long is = 1; is<=3; ++is)
	{
		if( ipISO == ipH_LIKE && is != 2 )
			continue;
		else if( ipISO == ipHE_LIKE && is != 1 && is != 3 )
			continue;

		char chSpin[3][9]= {"singlets", "doublets", "triplets"};

		/* give element number and spin */
		fprintf(ioQQQ," %s %s  %s bnl\n",
			iso_ctrl.chISO[ipISO],
			elementnames.chElementSym[nelem],
			chSpin[is-1]);

		/* header with the l states */
		fprintf(ioQQQ," n\\l=>    ");
		for( long i =0; i < iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max; ++i)
		{
			fprintf(ioQQQ,"%2ld         ",i);
		}
		fprintf(ioQQQ,"\n");

		/* loop over prin quant numbers, one per line, with l across */
		for( long in = 1; in <= iso_sp[ipISO][nelem].n_HighestResolved_max + iso_sp[ipISO][nelem].nCollapsed_max; ++in)
		{
			if( is==3 && in==1 )
				continue;

			fprintf(ioQQQ," %2ld      ",in);

			for( long il = 0; il < in; ++il)
			{
				fprintf( ioQQQ, "%9.3e ", iso_sp[ipISO][nelem].bnl_effective[in][il][is] );
			}
			fprintf(ioQQQ,"\n");
		}
	}

	return;
}

void iso_collapsed_Aul_update( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_Aul_update()" );

	long ipFirstCollapsed = iso_sp[ipISO][nelem].numLevels_max - iso_sp[ipISO][nelem].nCollapsed_max;

	for( long ipLo=0; ipLo<ipFirstCollapsed; ipLo++ )
	{
		long spin = iso_sp[ipISO][nelem].st[ipLo].S();

		/* calculate effective Aul's from collapsed levels */
		for( long nHi=iso_sp[ipISO][nelem].n_HighestResolved_max+1; nHi<=iso_sp[ipISO][nelem].n_HighestResolved_max+iso_sp[ipISO][nelem].nCollapsed_max; nHi++ )
		{
			realnum Auls[2] = {
				iso_sp[ipISO][nelem].CachedAs[ nHi-iso_sp[ipISO][nelem].n_HighestResolved_max-1 ][ ipLo ][0],
				iso_sp[ipISO][nelem].CachedAs[ nHi-iso_sp[ipISO][nelem].n_HighestResolved_max-1 ][ ipLo ][1] };

			realnum EffectiveAul = 
				Auls[0]*spin*(2.f*(L_(ipLo)+1.f)+1.f)*(realnum)iso_sp[ipISO][nelem].bnl_effective[nHi][ L_(ipLo)+1 ][spin];
			
			/* this is for n,L-1 -> n',L
			 * make sure L-1 exists. */
			if( L_(ipLo) > 0 )
			{
				EffectiveAul += 
					Auls[1]*spin*(2.f*(L_(ipLo)-1.f)+1.f)*(realnum)iso_sp[ipISO][nelem].bnl_effective[nHi][ L_(ipLo)-1 ][spin];
			}

			if( ipISO==ipH_LIKE )
				EffectiveAul /= (2.f*nHi*nHi);
			else if( ipISO==ipHE_LIKE )
				EffectiveAul /= (4.f*nHi*nHi);
			else
				TotalInsanity();

			long ipHi = iso_sp[ipISO][nelem].QuantumNumbers2Index[nHi][ L_(ipLo)+1 ][spin];

			/* FINALLY, put the effective A in the proper Emis structure. */
			iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() = EffectiveAul;

			ASSERT( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() > 0. );
		}
	}

	return;
}

void iso_collapsed_lifetimes_update( long ipISO, long nelem )
{
	DEBUG_ENTRY( "iso_collapsed_Aul_update()" );

	for( long ipHi=iso_sp[ipISO][nelem].numLevels_max- iso_sp[ipISO][nelem].nCollapsed_max; ipHi < iso_sp[ipISO][nelem].numLevels_max; ipHi++ )
	{
		iso_sp[ipISO][nelem].st[ipHi].lifetime() = SMALLFLOAT;

		for( long ipLo=0; ipLo < ipHi; ipLo++ )  
		{
			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
				continue;

			iso_sp[ipISO][nelem].st[ipHi].lifetime() += iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul();
		}

		/* sum of A's was just stuffed, now invert for lifetime. */
		iso_sp[ipISO][nelem].st[ipHi].lifetime() = 1./iso_sp[ipISO][nelem].st[ipHi].lifetime();

		for( long ipLo=0; ipLo < ipHi; ipLo++ )
		{
			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN() <= 0. )
				continue;

			if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().Aul() <= iso_ctrl.SmallA )
				continue;

			iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel() = (realnum)( 
				(1.f/iso_sp[ipISO][nelem].st[ipHi].lifetime())/
				PI4/iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyWN());

			ASSERT(iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().dampXvel()> 0.);
		}
	}

	return;
}

#if	0
STATIC void Prt_AGN_table( void )
{
	/* the designation of the levels, chLevel[n][string] */
	multi_arr<char,2> chLevel(max_num_levels,10);

	/* create spectroscopic designation of labels */
	for( long ipLo=0; ipLo < iso_sp[ipISO][ipISO].numLevels_max-iso_sp[ipISO][ipISO].nCollapsed_max; ++ipLo )
	{
		long nelem = ipISO;
		sprintf( &chLevel.front(ipLo) , "%li %li%c", N_(ipLo), S_(ipLo), chL[MIN2(20,L_(ipLo))] );
	}

	/* option to print cs data for AGN */
	/* create spectroscopic designation of labels */
	{
		/* option to print particulars of some line when called */
		enum {AGN=false};
		if( AGN )
		{
#			define NTEMP 6
			double te[NTEMP]={6000.,8000.,10000.,15000.,20000.,25000. };
			double telog[NTEMP] ,
				cs ,
				ratecoef;
			long nelem = ipHELIUM;
			fprintf( ioQQQ,"trans");
			for( long i=0; i < NTEMP; ++i )
			{
				telog[i] = log10( te[i] );
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			for( long i=0; i < NTEMP; ++i )
			{
				fprintf( ioQQQ,"\t%.3e",te[i]);
			}
			fprintf(ioQQQ,"\n");

			for( long ipHi=ipHe2s3S; ipHi< iso_sp[ipHE_LIKE][ipHELIUM].numLevels_max; ++ipHi )
			{
				for( long ipLo=ipHe1s1S; ipLo < ipHi; ++ipLo )
				{

					/* deltaN = 0 transitions may be wrong because 
					 * COLL_CONST below is only correct for electron colliders */
					if( N_(ipHi) == N_(ipLo) )
						continue;

					/* print the designations of the lower and upper levels */
					fprintf( ioQQQ,"%s - %s",
						 &chLevel.front(ipLo) , &chLevel.front(ipHi) );

					/* print the interpolated collision strengths */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						/* print cs */
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						fprintf(ioQQQ,"\t%.2e", cs );
					}

					/* print the rate coefficients */
					for( long i=0; i < NTEMP; ++i )
					{
						phycon.alogte = telog[i];
						phycon.te = pow(10.,telog[i] );
						tfidle(false);
						cs = HeCSInterp( nelem , ipHi , ipLo, ipELECTRON );
						/* collisional deexcitation rate */
						ratecoef = cs/sqrt(phycon.te)*COLL_CONST/iso_sp[ipHE_LIKE][nelem].st[ipLo].g() *
							sexp( iso_sp[ipHE_LIKE][nelem].trans(ipHi,ipLo).EnergyK() / phycon.te );
						fprintf(ioQQQ,"\t%.2e", ratecoef );
					}
					fprintf(ioQQQ,"\n");
				}
			}
			cdEXIT(EXIT_FAILURE);
		}
	}

	return;
}
#endif
