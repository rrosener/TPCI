/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*lines_helium put He-like iso sequence into line intensity stack */
/*TempInterp interpolates on a grid of values to produce predicted value at current Te.*/
#include "cddefines.h"
#include "dense.h"
#include "prt.h"
#include "helike.h"
#include "iso.h"
#include "atmdat.h"
#include "lines.h"
#include "lines_service.h"
#include "phycon.h"
#include "physconst.h"
#include "taulines.h"
#include "thirdparty.h"
#include "trace.h"

#define NUMTEMPS	21
#define NUMDENS		14		
typedef struct 
{
	/* index for upper and lower levels of line */
	long int ipHi;
	long int ipLo;

	char label[5];

} stdLines;

STATIC void GetStandardHeLines(void);
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te );
STATIC void DoSatelliteLines( long nelem );

static bool lgFirstRun = true;
static double CaBDensities[NUMDENS];
static double CaBTemps[NUMTEMPS];
static long NumLines;
static double ****CaBIntensity;
static stdLines **CaBLines;

void lines_helium(void)
{
	long ipISO = ipHE_LIKE;
	long int i, nelem, ipHi, ipLo;
	char chLabel[5]="    ";

	double 
	  sum,
	  photons_3889_plus_7065 = 0.;

	DEBUG_ENTRY( "lines_helium()" );

	if( trace.lgTrace )
		fprintf( ioQQQ, "   prt_lines_helium called\n" );

	// this can be changed with the atom levels command but must be at
	// least 3.
	ASSERT( iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 3 );

	i = StuffComment( "He-like iso-sequence" );
	linadd( 0., (realnum)i , "####", 'i',
		" start He-like iso sequence");

	linadd(MAX2(0.,iso_sp[ipHE_LIKE][ipHELIUM].xLineTotCool),506,"Clin",'c',
		"  total collisional cooling due to all HeI lines ");

	linadd(MAX2(0.,-iso_sp[ipHE_LIKE][ipHELIUM].xLineTotCool),506,"Hlin",'h'	,
		"  total collisional heating due to all HeI lines ");

	/* read in Case A and B lines from data file	*/
	if( lgFirstRun )
	{
		GetStandardHeLines();
		lgFirstRun = false;
	}

	/* this is the main printout, where line intensities are entered into the stack */
	for( nelem=ipISO; nelem < LIMELM; nelem++ )
	{
		if( dense.lgElmtOn[nelem] )
		{
			t_iso_sp* sp = &iso_sp[ipHE_LIKE][nelem];

			ASSERT( sp->n_HighestResolved_max >= 3 );

			if( nelem == ipHELIUM )
			{
				for( i=0; i< NumLines; i++ )
				{
					ipHi = CaBLines[nelem][i].ipHi;
					ipLo = CaBLines[nelem][i].ipLo;
					double intens_at_Te[NUMDENS];
					for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
						intens_at_Te[ipDens] = TempInterp2( CaBTemps, CaBIntensity[nelem][i][ipDens], NUMTEMPS, phycon.te );
					double intens = linint( CaBDensities, intens_at_Te, NUMDENS, log10(dense.eden) );
					intens = pow( 10., intens ) * dense.xIonDense[nelem][nelem+1-ipISO]*dense.eden;
					ASSERT( intens >= 0. );
					linadd( intens, atmdat.CaseBWlHeI[i], CaBLines[nelem][i].label,'i', "Case B intensity ");
				}
			}

			// add two-photon details here
			if( LineSave.ipass == 0 )
			{
				/* chIonLbl is function that generates a null terminated 4 char string, of form "C  2" 
				 * the result, chLable, is only used when ipass == 0, can be undefined otherwise */
				chIonLbl(chLabel, nelem+1, nelem+1-ipISO);
			}
			for( vector<two_photon>::iterator tnu = sp->TwoNu.begin(); tnu != sp->TwoNu.end(); ++tnu )
			{
				fixit(); // This was multiplied by Pesc when treated as a line, now what?  Only used for printout?
				fixit(); // below should be 'i' instead of 'r' ?
				linadd(	tnu->AulTotal * tnu->E2nu * EN1RYD * (*tnu->Pop), 
					 0, chLabel, 'r',
					" two photon continuum ");

				linadd(	tnu->induc_dn * tnu->E2nu * EN1RYD * (*tnu->Pop), 
					22, chLabel ,'i',
					" induced two photon emission ");
			}

			/* here we will create an entry for the three lines 
			 * coming from 2 3P to 1 1S - the entry called TOTL will
			 * appear before the lines of the multiplet */
			sum = 0.;
			for( i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
			{
				if( sp->trans(i,ipHe1s1S).ipCont() <= 0 )
					continue;

				sum += 
					sp->trans(i,ipHe1s1S).Emis().Aul()*
					sp->st[i].Pop()*
					(sp->trans(i,ipHe1s1S).Emis().Pesc() +
					sp->trans(i,ipHe1s1S).Emis().Pelec_esc() ) *
					sp->trans(i,ipHe1s1S).EnergyErg();
			}

			linadd(sum,sp->trans(ipHe2p3P1,ipHe1s1S).WLAng(),"TOTL",'i' ,
				" total emission in He-like intercombination lines from 2p3P to ground ");

			/* set number of levels we want to print, first is default,
			 * only print real levels, second is set with "print line
			 * iso collapsed" command */
			long int nLoop  = sp->numLevels_max - sp->nCollapsed_max;
			if( prt.lgPrnIsoCollapsed )
				nLoop  = sp->numLevels_max;

			/* now do real permitted lines */
			/* NB NB - low and high must be in this order so that all balmer, paschen,
			 * etc series line up correctly in final printout */
			/* >>chng 01 jun 13, bring 23P lines back together */
			for( ipLo=0; ipLo < ipHe2p3P0; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* >>chng 01 may 30, do not add fake he-like lines (majority) to line stack */
					/* >>chng 01 dec 11, use variable for smallest A */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 ) 
						continue;

					/* recombine fine-structure lines since the energies are
					 * not resolved anyway.	*/
					if( iso_ctrl.lgFSM[ipISO] && ( abs(sp->st[ipHi].l() -
						sp->st[ipLo].l())==1 )
						&& (sp->st[ipLo].l()>1) 
						&& (sp->st[ipHi].l()>1) 
						&& ( sp->st[ipHi].n() ==
						sp->st[ipLo].n() ) )
					{
						/* skip if both singlets. */
						if( (sp->st[ipHi].S()==1) 
							&& (sp->st[ipLo].S()==1) )
						{
							continue;
						}
						else if( (sp->st[ipHi].S()==3) 
							&& (sp->st[ipLo].S()==3) )
						{

							/* singlet to singlet*/
							sp->trans(ipHi+1,ipLo+1).Emis().phots() = 
								sp->trans(ipHi,ipLo+1).Emis().Aul()*
								sp->st[ipHi].Pop()*
								(sp->trans(ipHi,ipLo+1).Emis().Pesc() +
								sp->trans(ipHi,ipLo+1).Emis().Pelec_esc() ) +
								sp->trans(ipHi+1,ipLo+1).Emis().Aul()*
								sp->st[ipHi+1].Pop()*
								(sp->trans(ipHi+1,ipLo+1).Emis().Pesc() +
								sp->trans(ipHi+1,ipLo+1).Emis().Pelec_esc() );

							sp->trans(ipHi+1,ipLo+1).Emis().xIntensity() = 
								sp->trans(ipHi+1,ipLo+1).Emis().phots() *
								sp->trans(ipHi+1,ipLo+1).EnergyErg();

							PutLine(sp->trans(ipHi+1,ipLo+1), " ");

							/* triplet to triplet */
							sp->trans(ipHi,ipLo).Emis().phots() = 
								sp->trans(ipHi,ipLo).Emis().Aul()*
								sp->st[ipHi].Pop()*
								(sp->trans(ipHi,ipLo).Emis().Pesc() +
								sp->trans(ipHi,ipLo).Emis().Pelec_esc() ) +
								sp->trans(ipHi+1,ipLo).Emis().Aul()*
								sp->st[ipHi+1].Pop()*
								(sp->trans(ipHi+1,ipLo).Emis().Pesc() +
								sp->trans(ipHi+1,ipLo).Emis().Pelec_esc() );

							sp->trans(ipHi,ipLo).Emis().xIntensity() = 
								sp->trans(ipHi,ipLo).Emis().phots() *
								sp->trans(ipHi,ipLo).EnergyErg();

							PutLine(sp->trans(ipHi,ipLo), " ");
						}
					}

					else if( ipLo==ipHe2s3S && ipHi == ipHe2p3P0 )
					{
						/* here we will create an entry for the three lines 
						 * coming from 2 3P to 2 3S - the entry called TOTL will
						 * appear before the lines of the multiplet 
						 * for He I this is 10830 */

						realnum av_wl = 0.;
						sum = 0.;
						for( i=ipHe2p3P0; i <= ipHe2p3P2; i++ )
						{
							sum += 
								sp->trans(i,ipLo).Emis().Aul()*
								sp->st[i].Pop()*
								(sp->trans(i,ipLo).Emis().Pesc() +
								sp->trans(i,ipLo).Emis().Pelec_esc() )*
								sp->trans(i,ipLo).EnergyErg();
							av_wl += sp->trans(i,ipLo).WLAng();
						}
						av_wl /= 3.;
#						if 0
						{
#						include "elementnames.h"
#						include "prt.h"
						fprintf(ioQQQ,"DEBUG 2P - 2S multiplet wl %s ",
							elementnames.chElementSym[nelem] );
						prt_wl( ioQQQ , av_wl );
						fprintf(ioQQQ,"\n" );
						}
#						endif

						linadd(sum,av_wl,"TOTL",'i',
							"total emission in He-like lines, use average of three line wavelengths " );

						/* also add this with the regular label, so it is correctly picked up by assert case-b command */
						linadd(sum,av_wl,chLabel,'i',
							"total emission in He-like lines, use average of three line wavelengths " );

						/*>>chng 05 sep 8, added the following so that the component
						 * from ipHe2p3P0 is printed, in addition to the total. */

						/* find number of photons escaping cloud */
						sp->trans(ipHi,ipLo).Emis().phots() = 
							sp->trans(ipHi,ipLo).Emis().Aul()*
							sp->st[ipHi].Pop()*
							(sp->trans(ipHi,ipLo).Emis().Pesc() +
							sp->trans(ipHi,ipLo).Emis().Pelec_esc() );

						/* now find line intensity */
						sp->trans(ipHi,ipLo).Emis().xIntensity() = 
							sp->trans(ipHi,ipLo).Emis().phots()*
							sp->trans(ipHi,ipLo).EnergyErg();

						if( iso_ctrl.lgRandErrGen[ipISO] )
						{
							sp->trans(ipHi,ipLo).Emis().phots() *=
								sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
							sp->trans(ipHi,ipLo).Emis().xIntensity() *= 
								sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
						}
						PutLine(sp->trans(ipHi,ipLo), " ");
					}

					else
					{
						/* find number of photons escaping cloud */
						sp->trans(ipHi,ipLo).Emis().phots() = 
							sp->trans(ipHi,ipLo).Emis().Aul()*
							sp->st[ipHi].Pop()*
							(sp->trans(ipHi,ipLo).Emis().Pesc() +
							sp->trans(ipHi,ipLo).Emis().Pelec_esc() );

						/* now find line intensity */
						/* >>chng 01 jan 15, put cast double to force double evaluation */
						sp->trans(ipHi,ipLo).Emis().xIntensity() = 
							sp->trans(ipHi,ipLo).Emis().phots()*
							sp->trans(ipHi,ipLo).EnergyErg();

						if( iso_ctrl.lgRandErrGen[ipISO] )
						{
							sp->trans(ipHi,ipLo).Emis().phots() *=
								sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
							sp->trans(ipHi,ipLo).Emis().xIntensity() *= 
								sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
						}

						if( abs( L_(ipHi) - L_(ipLo) ) != 1 )
						{
							EnterTheseLast.push_back( ipHi );
							continue;
						}

						/* 
						fprintf(ioQQQ,"1 loop %li %li %.1f\n", ipLo,ipHi, 
							sp->trans(ipHi,ipLo).WLAng() ); */
						PutLine(sp->trans(ipHi,ipLo),
							"total intensity of He-like line");
						{
							/* option to print particulars of some line when called
							 * a prettier print statement is near where chSpin is defined below*/
							enum {DEBUG_LOC=false};
							if( DEBUG_LOC )
							{
								if( nelem==1 && ipLo==0 && ipHi==1 )
								{
									fprintf(ioQQQ,"he1 626 %.2e %.2e \n", 
										sp->trans(ipHi,ipLo).Emis().TauIn(),
										sp->trans(ipHi,ipLo).Emis().TauTot()
										);
								}
							}
						}
					}
				}
		
				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
					PutLine( sp->trans(*it,ipLo),
						"predicted line, all processes included");
			}

			/* this sum will bring together the three lines going to J levels within 23P */
			for( ipHi=ipHe2p3P2+1; ipHi < nLoop; ipHi++ )
			{
				double sumcool , sumheat ,
					save , savecool , saveheat;

				sum = 0;
				sumcool = 0.;
				sumheat = 0.;
				for( ipLo=ipHe2p3P0; ipLo <= ipHe2p3P2; ++ipLo )
				{
					if( sp->trans(ipHi,ipLo).ipCont() <= 0 )
						continue;

					/* find number of photons escaping cloud */
					sp->trans(ipHi,ipLo).Emis().phots() = 
						sp->trans(ipHi,ipLo).Emis().Aul()*
						sp->st[ipHi].Pop()*
						(sp->trans(ipHi,ipLo).Emis().Pesc() +
						sp->trans(ipHi,ipLo).Emis().Pelec_esc() );

					/* now find line intensity */
					/* >>chng 01 jan 15, put cast double to force double evaluation */
					sp->trans(ipHi,ipLo).Emis().xIntensity() = 
						sp->trans(ipHi,ipLo).Emis().phots()*
						sp->trans(ipHi,ipLo).EnergyErg();

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						sp->trans(ipHi,ipLo).Emis().phots() *=
							sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
						sp->trans(ipHi,ipLo).Emis().xIntensity() *= 
							sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
					}

					sumcool += sp->trans(ipHi,ipLo).Coll().cool();
					sumheat += sp->trans(ipHi,ipLo).Coll().heat();
					sum += sp->trans(ipHi,ipLo).Emis().xIntensity();
				}

				/* skip non-radiative lines */
				if( sp->trans(ipHi,ipHe2p3P2).ipCont() < 1 ) 
					continue;

				/* this will enter .xIntensity() into the line stack */
				save = sp->trans(ipHi,ipHe2p3P2).Emis().xIntensity();
				savecool = sp->trans(ipHi,ipHe2p3P2).Coll().cool();
				saveheat = sp->trans(ipHi,ipHe2p3P2).Coll().heat();

				sp->trans(ipHi,ipHe2p3P2).Emis().xIntensity() = sum;
				sp->trans(ipHi,ipHe2p3P2).Coll().cool() = sumcool;
				sp->trans(ipHi,ipHe2p3P2).Coll().heat() = sumheat;

				/*fprintf(ioQQQ,"2 loop %li %li %.1f\n", ipHe2p3P2,ipHi, 
					sp->trans(ipHi,ipHe2p3P2).WLAng() );*/
				PutLine(sp->trans(ipHi,ipHe2p3P2),
					"predicted line, all processes included");

				sp->trans(ipHi,ipHe2p3P2).Emis().xIntensity() = save;
				sp->trans(ipHi,ipHe2p3P2).Coll().cool() = savecool;
				sp->trans(ipHi,ipHe2p3P2).Coll().heat() = saveheat;
			}
			for( ipLo=ipHe2p3P2+1; ipLo < nLoop-1; ipLo++ )
			{
				vector<long> EnterTheseLast;
				for( ipHi=ipLo+1; ipHi < nLoop; ipHi++ )
				{
					/* skip non-radiative lines */
					if( sp->trans(ipHi,ipLo).ipCont() < 1 ) 
						continue;

					/* find number of photons escaping cloud */
					sp->trans(ipHi,ipLo).Emis().phots() = 
						sp->trans(ipHi,ipLo).Emis().Aul()*
						sp->st[ipHi].Pop()*
						(sp->trans(ipHi,ipLo).Emis().Pesc() +
						sp->trans(ipHi,ipLo).Emis().Pelec_esc() );

					/* now find line intensity */
					/* >>chng 01 jan 15, put cast double to force double evaluation */
					sp->trans(ipHi,ipLo).Emis().xIntensity() = 
						sp->trans(ipHi,ipLo).Emis().phots()*
						sp->trans(ipHi,ipLo).EnergyErg();

					if( iso_ctrl.lgRandErrGen[ipISO] )
					{
						sp->trans(ipHi,ipLo).Emis().phots() *=
							sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
						sp->trans(ipHi,ipLo).Emis().xIntensity() *= 
							sp->ex[ipHi][ipLo].ErrorFactor[IPRAD];
					}
						
					if( abs( L_(ipHi) - L_(ipLo) ) != 1 )
					{
						EnterTheseLast.push_back( ipHi );
						continue;
					}
				
					/* this will enter .xIntensity() into the line stack */
					PutLine(sp->trans(ipHi,ipLo),
						"predicted line, all processes included");
				}

				// enter these lines last because they are generally weaker quadrupole transitions.
				for( vector<long>::iterator it = EnterTheseLast.begin(); it != EnterTheseLast.end(); it++ )
					PutLine(sp->trans(*it,ipLo),
						"predicted line, all processes included");
			}

			/* Now put the satellite lines in */
			if( iso_ctrl.lgDielRecom[ipISO] )
				DoSatelliteLines(nelem);
		}
	}

	if( iso_sp[ipHE_LIKE][ipHELIUM].n_HighestResolved_max >= 4 &&
		( iso_sp[ipH_LIKE][ipHYDROGEN].n_HighestResolved_max + iso_sp[ipH_LIKE][ipHYDROGEN].nCollapsed_max ) >=8 )
	{
		t_iso_sp* sp = &iso_sp[ipHE_LIKE][ipHELIUM];
		const long ipHe4s3S = iso_sp[ipHE_LIKE][ipHELIUM].QuantumNumbers2Index[4][0][3];
		const long ipHe4p3P = iso_sp[ipHE_LIKE][ipHELIUM].QuantumNumbers2Index[4][1][3];

		/* For info only, add the total photon flux in 3889 and 7065,
		* and 3188, 4713, and 5876. */
		photons_3889_plus_7065 =
			/* to 2p3P2 */
			sp->trans(ipHe3s3S,ipHe2p3P2).Emis().xIntensity()/
			sp->trans(ipHe3s3S,ipHe2p3P2).EnergyErg() +
			sp->trans(ipHe3d3D,ipHe2p3P2).Emis().xIntensity()/
			sp->trans(ipHe3d3D,ipHe2p3P2).EnergyErg() +
			sp->trans(ipHe4s3S,ipHe2p3P2).Emis().xIntensity()/
			sp->trans(ipHe4s3S,ipHe2p3P2).EnergyErg() +
			/* to 2p3P1 */
			sp->trans(ipHe3s3S,ipHe2p3P1).Emis().xIntensity()/
			sp->trans(ipHe3s3S,ipHe2p3P1).EnergyErg() +
			sp->trans(ipHe3d3D,ipHe2p3P1).Emis().xIntensity()/
			sp->trans(ipHe3d3D,ipHe2p3P1).EnergyErg() +
			sp->trans(ipHe4s3S,ipHe2p3P1).Emis().xIntensity()/
			sp->trans(ipHe4s3S,ipHe2p3P1).EnergyErg() +
			/* to 2p3P0 */
			sp->trans(ipHe3s3S,ipHe2p3P0).Emis().xIntensity()/
			sp->trans(ipHe3s3S,ipHe2p3P0).EnergyErg() +
			sp->trans(ipHe3d3D,ipHe2p3P0).Emis().xIntensity()/
			sp->trans(ipHe3d3D,ipHe2p3P0).EnergyErg() +
			sp->trans(ipHe4s3S,ipHe2p3P0).Emis().xIntensity()/
			sp->trans(ipHe4s3S,ipHe2p3P0).EnergyErg() +
			/* to 2s3S */
			sp->trans(ipHe3p3P,ipHe2s3S).Emis().xIntensity()/
			sp->trans(ipHe3p3P,ipHe2s3S).EnergyErg() +
			sp->trans(ipHe4p3P,ipHe2s3S).Emis().xIntensity()/
			sp->trans(ipHe4p3P,ipHe2s3S).EnergyErg();

		long upperIndexofH8 = iso_sp[ipH_LIKE][ipHYDROGEN].QuantumNumbers2Index[8][1][2];

		/* Add in photon flux of H8 3889 */
		photons_3889_plus_7065 += 
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(upperIndexofH8,1).Emis().xIntensity()/
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(upperIndexofH8,1).EnergyErg();

		/* now multiply by ergs of normalization line, so that relative flux of
		* this line will be ratio of photon fluxes. */
		if( LineSave.WavLNorm > 0 )
			photons_3889_plus_7065 *= (ERG1CM*1.e8)/LineSave.WavLNorm;
		linadd( photons_3889_plus_7065, 3889., "Pho+", 'i',
			"photon sum given in Porter et al. 2007 (astro-ph/0611579)");
	}

	/* ====================================================
	 * end helium
	 * ====================================================*/

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, "   lines_helium returns\n" );
	}
	return;
}


STATIC void GetStandardHeLines(void)
{
	FILE *ioDATA;
	bool lgEOL, lgHIT;
	long i, i1, i2;

#	define chLine_LENGTH 1000
	char chLine[chLine_LENGTH];

	DEBUG_ENTRY( "GetStandardHeLines()" );

	if( trace.lgTrace )
		fprintf( ioQQQ," lines_helium opening he1_case_b.dat\n");

	ioDATA = open_data( "he1_case_b.dat", "r" );

	/* check that magic number is ok */
	if( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) == NULL )
	{
		fprintf( ioQQQ, " lines_helium could not read first line of he1_case_b.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}
	i = 1;
	/* first number is magic number, second is number of lines in file	*/
	i1 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	i2 = (long)FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
	NumLines = i2;

	/* the following is to check the numbers that appear at the start of he1_case_b.dat */
	if( i1 !=CASEBMAGIC )
	{
		fprintf( ioQQQ, 
			" lines_helium: the version of he1_case_ab.dat is not the current version.\n" );
		fprintf( ioQQQ, 
			" lines_helium: I expected to find the number %i and got %li instead.\n" ,
			CASEBMAGIC, i1 );
		fprintf( ioQQQ, "Here is the line image:\n==%s==\n", chLine );
		cdEXIT(EXIT_FAILURE);
	}

	/* get the array of densities */
	lgHIT = false;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			lgHIT = true;
			long j = 0;
			lgEOL = false;
			i = 1;
			while( !lgEOL && j < NUMDENS)
			{
				CaBDensities[j] = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
				++j;
			}
			break;
		}
	}

	/* get the array of temperatures */
	lgHIT = false;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		/* only look at lines without '#' in first col */
		if( chLine[0] != '#')
		{
			lgHIT = true;
			long j = 0;
			lgEOL = false;
			i = 1;
			while( !lgEOL && j < NUMTEMPS)
			{
				CaBTemps[j] = FFmtRead(chLine,&i,sizeof(chLine),&lgEOL);
				++j;
			}
			break;
		}
	}

	if( !lgHIT )
	{
		fprintf( ioQQQ, " lines_helium could not find line of temperatures in he1_case_ab.dat.\n");
		cdEXIT(EXIT_FAILURE);
	}

	/* create space for array of values, if not already done */
	CaBIntensity = (double ****)MALLOC(sizeof(double ***)*(unsigned)LIMELM );
	/* create space for array of values, if not already done */
	CaBLines = (stdLines **)MALLOC(sizeof(stdLines *)*(unsigned)LIMELM );

	for( long nelem=ipHELIUM; nelem<LIMELM; ++nelem )
	{
		/** \todo	2	- this structure is currently only used for helium itself...
		 * stuff numbers in for other elements, or drop the [nelem] dimension off
		 * of CaBLines	*/
		if( nelem != ipHELIUM )
			continue;

		/* only grab core for elements that are turned on */
		if( nelem == ipHELIUM || dense.lgElmtOn[nelem])
		{
			/* create space for array of values, if not already done */
			CaBIntensity[nelem] = (double ***)MALLOC(sizeof(double **)*(unsigned)(i2) );
			CaBLines[nelem] = (stdLines *)MALLOC(sizeof(stdLines )*(unsigned)(i2) );

			/* avoid allocating 0 bytes, some OS return NULL pointer, PvH 
			CaBIntensity[nelem][0] = NULL;*/
			for( long j = 0; j < i2; ++j )
			{
				CaBIntensity[nelem][j] = (double **)MALLOC(sizeof(double*)*(unsigned)NUMDENS );

				CaBLines[nelem][j].ipHi = -1;
				CaBLines[nelem][j].ipLo = -1;
				strcpy( CaBLines[nelem][j].label , "    " );

				for( long k = 0; k < NUMDENS; ++k )
				{
					CaBIntensity[nelem][j][k] = (double *)MALLOC(sizeof(double)*(unsigned)NUMTEMPS );
					for( long l = 0; l < NUMTEMPS; ++l )
					{
						CaBIntensity[nelem][j][k][l] = 0.;
					}
				}
			}
		}
	}

	/* now read in the case B line data */
	lgHIT = false;
	long nelem = ipHELIUM;
	int line = 0;
	while( read_whole_line( chLine , (int)sizeof(chLine) , ioDATA ) != NULL )
	{
		char *chTemp = NULL;

		/* only look at lines without '#' in first col */
		if( (chLine[0] == ' ') || (chLine[0]=='\n') )
			break;
		if( chLine[0] != '#')
		{
			lgHIT = true;

			/* get lower and upper level index */
			long j = 1;
			// the first number is the wavelength, which is only used if the
			// model atom is too small to include this transition
			realnum wl = (realnum)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			long ipLo = (long)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			long ipHi = (long)FFmtRead(chLine,&j,sizeof(chLine),&lgEOL);
			CaBLines[nelem][line].ipLo = ipLo;
			CaBLines[nelem][line].ipHi = ipHi;

			ASSERT( CaBLines[nelem][line].ipLo < CaBLines[nelem][line].ipHi );

			strncpy( CaBLines[nelem][line].label, "Ca B" , 4 );
			CaBLines[nelem][line].label[4] = 0;

			t_iso_sp* sp = &iso_sp[ipHE_LIKE][nelem];
			if( ipHi < sp->numLevels_max - sp->nCollapsed_max )
				atmdat.CaseBWlHeI.push_back( sp->trans(ipHi,ipLo).WLAng() );
			else
				atmdat.CaseBWlHeI.push_back( wl );
		
			for( long ipDens = 0; ipDens < NUMDENS; ++ipDens )
			{
				if( read_whole_line( chLine, (int)sizeof(chLine) , ioDATA ) == NULL )
				{
					fprintf( ioQQQ, " lines_helium could not scan case B lines, current indices: %li %li\n",
						CaBLines[nelem][line].ipHi,
						CaBLines[nelem][line].ipLo );
					cdEXIT(EXIT_FAILURE);
				}

				chTemp = chLine;
				j = 1;
				long den = (long)FFmtRead(chTemp,&j,sizeof(chTemp),&lgEOL);
				ASSERT( den == ipDens + 1 );
			
				for( long ipTe = 0; ipTe < NUMTEMPS; ++ipTe )
				{
					double b;
					if( (chTemp = strchr_s( chTemp, '\t' )) == NULL )
					{
						fprintf( ioQQQ, " lines_helium could not scan case B lines, current indices: %li %li\n",
							CaBLines[nelem][line].ipHi,
							CaBLines[nelem][line].ipLo );
						cdEXIT(EXIT_FAILURE);
					}
					++chTemp;
					sscanf( chTemp, "%le" , &b );
					CaBIntensity[nelem][line][ipDens][ipTe] = b;
				}
			}
			line++;
		}
	}

	ASSERT( line == NumLines );
	ASSERT( atmdat.CaseBWlHeI.size() == (unsigned)line );

	/* close the data file */
	fclose( ioDATA );
	return;
}

/** \todo	there is a virtually identical routine in helike_recom.cpp -> combine */
STATIC double TempInterp2( double* TempArray , double* ValueArray, long NumElements, double Te )
{
	long int ipTe=-1;
	double rate = 0.;
	long i0;

	DEBUG_ENTRY( "TempInterp2()" );

	/* te totally unknown */
	if( Te <= TempArray[0] )
	{
		return ValueArray[0] + log10( TempArray[0] ) - log10( Te );
	}
	else if( Te >= TempArray[NumElements-1] )
	{
		return ValueArray[NumElements-1];
	}

	/* now search for temperature */
	ipTe = hunt_bisect( TempArray, NumElements, Te );			

	ASSERT( (ipTe >=0) && (ipTe < NumElements-1)  );

	/* Do a four-point interpolation */
	const int ORDER = 3; /* order of the fitting polynomial */
	i0 = max(min(ipTe-ORDER/2,NumElements-ORDER-1),0);
	rate = lagrange( &TempArray[i0], &ValueArray[i0], ORDER+1, Te );

	return rate;
}

/** \todo	2	say where these come from	*/	
/* For double-ionization discussions, see Lindsay, Rejoub, & Stebbings 2002	*/
/* Also read Itza-Ortiz, Godunov, Wang, and McGuire 2001.	*/
STATIC void DoSatelliteLines( long nelem )
{
	long ipISO = ipHE_LIKE;
	
	DEBUG_ENTRY( "DoSatelliteLines()" );

	ASSERT( iso_ctrl.lgDielRecom[ipISO] && dense.lgElmtOn[nelem] );

	for( long i=0; i < iso_sp[ipISO][nelem].numLevels_max; i++ )
	{
		double dr_rate = iso_sp[ipISO][nelem].fb[i].DielecRecomb;
		TransitionProxy tr = SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][i]];

		tr.Emis().phots() = 
			dr_rate * dense.eden * dense.xIonDense[nelem][nelem+1-ipISO];
		
		tr.Emis().xIntensity() = 
			tr.Emis().phots() * ERG1CM * tr.EnergyWN();
		tr.Emis().pump() = 0.;

		PutLine( tr, "iso satellite line" );
	}

	return;
}
