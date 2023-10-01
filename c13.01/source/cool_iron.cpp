/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*CoolIron compute iron cooling */
/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion */
/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
/*Fe2_cooling compute cooling due to FeII emission */
/*Fe3_cs return collision strengths for low level Fe3 lines */
/*Fe4_cs return collision strengths for low level Fe4 lines */
/*Fe5_cs return collision strengths for low level Fe5 lines */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "coolheavy.h"
#include "taulines.h"
#include "phycon.h"
#include "iso.h"
#include "conv.h"
#include "cooling.h"
#include "trace.h"
#include "hydrogenic.h"
#include "ligbar.h"
#include "cooling.h"
#include "thermal.h"
#include "lines_service.h"
#include "atoms.h"
#include "atomfeii.h"
#include "fe.h"

/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion */
STATIC void Fe3Lev14(void);

/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
STATIC void Fe4Lev12(void);

/*Fe2_cooling compute cooling due to FeII emission */
STATIC void Fe2_cooling( void )
{
	const int NLFE2 = 6;
	long int i, j;
	int nNegPop;

	static double **AulPump,
		**CollRate,
		**AulEscp,
		**col_str ,
		**AulDest, 
		*depart,
		*pops,
		*destroy,
		*create;

	static bool lgFirst=true;
	bool lgZeroPop;

	/* stat weights for Fred's 6 level model FeII atom */
	static double gFe2[NLFE2]={1.,1.,1.,1.,1.,1.};
	/* excitation energies (Kelvin) for Fred's 6 level model FeII atom */
	static double ex[NLFE2]={0.,3.32e4,5.68e4,6.95e4,1.15e5,1.31e5};

	/* these are used to only evaluated FeII one time per temperature, zone, 
	 * and abundance */
	static double TUsed = 0.; 
	static double AbunUsed = 0.;
	/* remember which zone last evaluated with */
	static long int nZUsed=-1,
		/* make sure at least two calls per zone */
		nCall=0;

	DEBUG_ENTRY( "Fe2_cooling()" );

	/* return if nothing to do */
	if( dense.xIonDense[ipIRON][1] == 0. )
	{
		/* zero abundance, do nothing */
		/* heating cooling and derivative from large atom */
		FeII.Fe2_large_cool = 0.;
		FeII.Fe2_large_heat = 0.;
		FeII.ddT_Fe2_large_cool = 0.;

		/* cooling and derivative from simple UV atom */
		FeII.Fe2_UVsimp_cool = 0.;
		FeII.ddT_Fe2_UVsimp_cool = 0.;

		/* now zero out intensities of all FeII lines and level populations */
		FeIIIntenZero();
		return;
	}

	/* evaluate FeII model atom, by default only the lowest 16 levels
	 * but number of levels can be increased with atom feii command */

	/* totally reevaluate FeII atom if new zone, or cooling is significant
	 * and temperature changed, we are in search phase,
	 * lgSlow option set true with atom FeII slow, forces constant
	 * evaluation of atom */
	if( FeII.lgSlow || 
		conv.lgFirstSweepThisZone || conv.lgLastSweepThisZone ||
	    /* check whether things have changed on later calls */
	    ( !fp_equal( phycon.te, TUsed ) && fabs(FeII.Fe2_large_cool/thermal.ctot) > 0.002 &&  
	      fabs(dense.xIonDense[ipIRON][1]-AbunUsed)/SDIV(AbunUsed) > 0.002 ) ||
	    ( !fp_equal( phycon.te, TUsed ) && fabs(FeII.Fe2_large_cool/thermal.ctot) > 0.01) )
	{

		if( conv.lgFirstSweepThisZone )
			/* first call this zone set nCall to zero*/
			nCall = 0;
		else
			/* not first call, increment, check above to make sure at least
			 * two evaluations */
			++nCall;

		/* option to trace convergence and FeII calls */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, "        CoolIron5 calling FeIILevelPops since ");
			if( conv.lgFirstSweepThisZone )
			{
				fprintf( ioQQQ, 
					"first sweep this zone." );
			}
			else if( ! fp_equal( phycon.te, TUsed ) )
			{
				fprintf( ioQQQ, 
					"temperature changed, old new are %g %g, nCall %li ", 
					TUsed, phycon.te , nCall);
			}
			else if( nzone != nZUsed )
			{
				fprintf( ioQQQ, 
					"new zone, nCall %li ", nCall );
			}
			else if( FeII.lgSlow )
			{
				fprintf( ioQQQ, 
					"FeII.lgSlow set  %li", nCall );
			}
			else if( conv.lgSearch )
			{
				fprintf( ioQQQ, 
					" in search phase  %li", nCall );
			}
			else if( nCall < 2 )
			{
				fprintf( ioQQQ, 
					"not second nCall %li " , nCall );
			}
			else if( ! fp_equal( phycon.te, TUsed ) && FeII.Fe2_large_cool/thermal.ctot > 0.001 )
			{
				fprintf( ioQQQ, 
					"temp or cooling changed, new are %g %g nCall %li ", 
					phycon.te, FeII.Fe2_large_cool, nCall );
			}
			else
			{
				fprintf(ioQQQ, "????");
			}
			fprintf(ioQQQ, "\n");
		}

		/* remember parameters for current conditions */
		TUsed = phycon.te;
		AbunUsed = dense.xIonDense[ipIRON][1];
		nZUsed = nzone;

		/* this print turned on with atom FeII print command */
		if( FeII.lgPrint )
		{
			fprintf(ioQQQ,
				" FeIILevelPops called zone %4li te %5f abun %10e c(fe/tot):%6f nCall %li\n", 
				nzone,phycon.te,AbunUsed,FeII.Fe2_large_cool/thermal.ctot,nCall);
		}

		/* this solves the multi-level problem, 
		 * sets FeII.Fe2_large_cool, 
				FeII.Fe2_large_heat, & 
				FeII.ddT_Fe2_large_cool 
				but does nothing with them */

		// line RT update, followed by level population solver
		FeII_RT_Make();
		FeIILevelPops();
		{
			enum{DEBUG_LOC=false};
			if( DEBUG_LOC && iteration > 1 && nzone >=4 )
			{
				fprintf(ioQQQ,"DEBUG1\t%li\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", 
					nzone,
					phycon.te,
					dense.gas_phase[ipHYDROGEN],
					dense.eden,
					FeII.Fe2_large_cool , 
					FeII.ddT_Fe2_large_cool ,
					FeII.Fe2_large_cool/dense.eden/dense.gas_phase[ipHYDROGEN] , 
					thermal.ctot );
			}
		}

		if( trace.nTrConvg >= 5 || FeII.lgPrint)
		{
			/* spacing needed to get proper trace convergence output */
			fprintf( ioQQQ, "           FeIILevelPops5 returned cool=%.2e heat=%.2e derivative=%.2e\n",
					FeII.Fe2_large_cool,FeII.Fe2_large_heat ,FeII.ddT_Fe2_large_cool);
		}

	}
	else if( ! fp_equal( dense.xIonDense[ipIRON][1], AbunUsed ) )
	{
		realnum ratio;
		/* this branch, same zone and temperature, but small change in abundance, so just
		 * rescale cooling and derivative by this change.  assumption is that very small changes
		 * in abundance occurs as ots rates damp out */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, 
				"       CoolIron rescaling FeIILevelPops since small change, CFe2=%.2e CTOT=%.2e\n",
				FeII.Fe2_large_cool,thermal.ctot);
		}
		ratio = dense.xIonDense[ipIRON][1]/AbunUsed;
		FeII.Fe2_large_cool *= ratio;
		FeII.ddT_Fe2_large_cool *= ratio;
		FeII.Fe2_large_heat *= ratio;
		AbunUsed = dense.xIonDense[ipIRON][1];
	}
	else
	{
		/* this is case where temp is unchanged, so heating and cooling same too */
		if( trace.nTrConvg >= 5 )
		{
			fprintf( ioQQQ, "       CoolIron NOT calling FeIILevelPops\n");
		}
	}

	/* now update total cooling and its derivative 
	 * all paths flow through here */
	/* FeII.Fecool = FeII.Fe2_large_cool; */
	CoolAdd("Fe 2",0,MAX2(0.,FeII.Fe2_large_cool));

	/* add negative cooling to heating stack */
	thermal.heating[25][27] = MAX2(0.,FeII.Fe2_large_heat);

	/* counts as heating derivative if negative cooling */
	if( FeII.Fe2_large_cool > 0. )
	{
		/* >>chng 01 mar 16, add factor of 3 due to conv problems after changing damper */
		thermal.dCooldT += 3.*FeII.ddT_Fe2_large_cool;
	}

	if( trace.lgTrace && trace.lgCoolTr )
	{
		fprintf( ioQQQ, " Large FeII returns te, cooling, dc=%11.3e%11.3e%11.3e\n", 
		  phycon.te, FeII.Fe2_large_cool, FeII.ddT_Fe2_large_cool );
	}

	/* >>chng 05 nov 29, still do simple UV atom if only ground term is done */
	if( !FeII.lgFeIILargeOn )
	{

		/* following treatment of Fe II follows
		 * >>refer	fe2	model	Wills, B.J., Wills, D., Netzer, H. 1985, ApJ, 288, 143
		 * all elements are used, and must be set to zero if zero */

		/* set up space for simple  model of UV FeII emission */
		if( lgFirst )
		{
			/* will never do this again in this core load */
			lgFirst = false;
			/* allocate the 1D arrays*/
			pops = (double *)MALLOC( sizeof(double)*(NLFE2) );
			create = (double *)MALLOC( sizeof(double)*(NLFE2) );
			destroy = (double *)MALLOC( sizeof(double)*(NLFE2) );
			depart = (double *)MALLOC( sizeof(double)*(NLFE2) );
			/* create space for the 2D arrays */
			AulPump = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			CollRate = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			AulDest = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			AulEscp = ((double **)MALLOC((NLFE2)*sizeof(double *)));
			col_str = ((double **)MALLOC((NLFE2)*sizeof(double *)));

			for( i=0; i < NLFE2; ++i )
			{
				AulPump[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				CollRate[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				AulDest[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				AulEscp[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
				col_str[i] = ((double *)MALLOC((NLFE2)*sizeof(double )));
			}
		}

		/*zero out all arrays, then check that upper diagonal remains zero below */
		for( i=0; i < NLFE2; i++ )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE2; j++ )
			{
				/*data[j][i] = 0.;*/
				col_str[j][i] = 0.;
				AulEscp[j][i] = 0.;
				AulDest[j][i] = 0.;
				AulPump[j][i] = 0.;
			}
		}

		/* now put in real data for lines */
		AulEscp[1][0] = 1.;
		AulEscp[2][0] = ( TauLines[ipTuv3].Emis().Pesc() + TauLines[ipTuv3].Emis().Pelec_esc())*TauLines[ipTuv3].Emis().Aul();
		AulDest[2][0] = TauLines[ipTuv3].Emis().Pdest()*TauLines[ipTuv3].Emis().Aul();
		AulPump[0][2] = TauLines[ipTuv3].Emis().pump();

		AulEscp[5][0] = (TauLines[ipTFe16].Emis().Pesc() + TauLines[ipTFe16].Emis().Pelec_esc())*TauLines[ipTFe16].Emis().Aul();
		AulDest[5][0] = TauLines[ipTFe16].Emis().Pdest()*TauLines[ipTFe16].Emis().Aul();
		/* continuum pumping of n=6 */
		AulPump[0][5] = TauLines[ipTFe16].Emis().pump();
		/* Ly-alpha pumping */

		double PumpLyaFeII = iso_sp[ipH_LIKE][ipHYDROGEN].st[ipH2p].Pop()*
			iso_sp[ipH_LIKE][ipHYDROGEN].trans(ipH2p,ipH1s).Emis().Aul()*
			hydro.dstfe2lya/SDIV(dense.xIonDense[ipIRON][1]);

		AulPump[0][5] += PumpLyaFeII;

		AulEscp[2][1] = (TauLines[ipTr48].Emis().Pesc() + TauLines[ipTr48].Emis().Pelec_esc())*TauLines[ipTr48].Emis().Aul();
		AulDest[2][1] = TauLines[ipTr48].Emis().Pdest()*TauLines[ipTr48].Emis().Aul();
		AulPump[1][2] = TauLines[ipTr48].Emis().pump();

		AulEscp[5][1] = (TauLines[ipTFe26].Emis().Pesc() + TauLines[ipTFe26].Emis().Pelec_esc())*TauLines[ipTFe26].Emis().Aul();
		AulDest[5][1] = TauLines[ipTFe26].Emis().Pdest()*TauLines[ipTFe26].Emis().Aul();
		AulPump[1][5] = TauLines[ipTFe26].Emis().pump();

		AulEscp[3][2] = (TauLines[ipTFe34].Emis().Pesc() + TauLines[ipTFe34].Emis().Pelec_esc())*TauLines[ipTFe34].Emis().Aul();
		AulDest[3][2] = TauLines[ipTFe34].Emis().Pdest()*TauLines[ipTFe34].Emis().Aul();
		AulPump[2][3] = TauLines[ipTFe34].Emis().pump();

		AulEscp[4][2] = (TauLines[ipTFe35].Emis().Pesc() + TauLines[ipTFe35].Emis().Pelec_esc())*TauLines[ipTFe35].Emis().Aul();
		AulDest[4][2] = TauLines[ipTFe35].Emis().Pdest()*TauLines[ipTFe35].Emis().Aul();
		AulPump[2][4] = TauLines[ipTFe35].Emis().pump();

		AulEscp[5][3] = (TauLines[ipTFe46].Emis().Pesc() + TauLines[ipTFe46].Emis().Pelec_esc())*TauLines[ipTFe46].Emis().Aul();
		AulDest[5][3] = TauLines[ipTFe46].Emis().Pdest()*TauLines[ipTFe46].Emis().Aul();
		AulPump[3][5] = TauLines[ipTFe46].Emis().pump();

		AulEscp[5][4] = (TauLines[ipTFe56].Emis().Pesc() + TauLines[ipTFe56].Emis().Pelec_esc())*TauLines[ipTFe56].Emis().Aul();
		AulDest[5][4] = TauLines[ipTFe56].Emis().Pdest()*TauLines[ipTFe56].Emis().Aul();
		AulPump[4][5] = TauLines[ipTFe56].Emis().pump();

		/* these are collision strengths */
		col_str[1][0] = 1.;
		col_str[2][0] = 12.;
		col_str[3][0] = 1.;
		col_str[4][0] = 1.;
		col_str[5][0] = 12.;
		col_str[2][1] = 6.;
		col_str[3][1] = 1.;
		col_str[4][1] = 1.;
		col_str[5][1] = 12.;
		col_str[3][2] = 6.;
		col_str[4][2] = 12.;
		col_str[5][2] = 1.;
		col_str[4][3] = 1.;
		col_str[5][3] = 12.;
		col_str[5][4] = 6.;

		/*void atom_levelN(long,long,realnum,double[],double[],double[],double*,
		double*,double*,long*,realnum*,realnum*,STRING,int*);*/
		atom_levelN(NLFE2,
			dense.xIonDense[ipIRON][1],
			gFe2,
			ex,
			'K',
			pops,
			depart,
			&AulEscp ,
			&col_str,
			&AulDest,
			&AulPump,
			&CollRate,
			create,
			destroy,
			false,/* say atom_levelN should evaluate coll rates from cs */
			/*&&ipdest,*/
			&FeII.Fe2_UVsimp_cool,
			&FeII.ddT_Fe2_UVsimp_cool,
			"FeII",
			&nNegPop,
			&lgZeroPop,
			false );

		/* nNegPop positive if negative pops occurred, negative if too cold */
		if( nNegPop > 0 /*negative if too cold - that is not negative and is OK */ )
		{
			fprintf(ioQQQ," PROBLEM, atom_levelN returned negative population for simple UV FeII.\n");
		}

		/* add heating - cooling by this process */;
		CoolAdd("Fe 2",0,MAX2(0.,FeII.Fe2_UVsimp_cool));
		thermal.heating[25][27] = MAX2(0.,-FeII.Fe2_UVsimp_cool);
		thermal.dCooldT += FeII.ddT_Fe2_UVsimp_cool;

		/* LIMLEVELN is the dim of the PopLevels vector */
		ASSERT( NLFE2 <= LIMLEVELN );
		for( i=0; i < NLFE2; ++i )
		{
			atoms.PopLevels[i] = pops[i];
			atoms.DepLTELevels[i] = depart[i];
		}

		(*TauLines[ipTuv3].Lo()).Pop() = pops[0];
		(*TauLines[ipTuv3].Hi()).Pop() = pops[2];
		TauLines[ipTuv3].Emis().PopOpc() = (pops[0] - pops[2]);
		TauLines[ipTuv3].Emis().phots() = pops[2]*AulEscp[2][0];
		TauLines[ipTuv3].Emis().xIntensity() = 
			TauLines[ipTuv3].Emis().phots()*TauLines[ipTuv3].EnergyErg();

		(*TauLines[ipTr48].Lo()).Pop() = pops[1];
		(*TauLines[ipTr48].Hi()).Pop() = pops[2];
		TauLines[ipTr48].Emis().PopOpc() = (pops[1] - pops[2]);
		TauLines[ipTr48].Emis().phots() = pops[2]*AulEscp[2][1];
		TauLines[ipTr48].Emis().xIntensity() = 
			TauLines[ipTr48].Emis().phots()*TauLines[ipTr48].EnergyErg();

		FeII.for7 = pops[1]*AulEscp[1][0]*4.65e-12;

		(*TauLines[ipTFe16].Lo()).Pop() = pops[0];
		(*TauLines[ipTFe16].Hi()).Pop() = pops[5];
		/* Lyman alpha optical depths are not known on first iteration,
		 * inward optical depths used, so line trapping overestimated,
		 * this can cause artificial maser in FeII - prevent by not
		 * including stimulated emission correction on first iteration */
		TauLines[ipTFe16].Emis().PopOpc() = (pops[0] - pops[5]);
		TauLines[ipTFe16].Emis().phots() = pops[5]*AulEscp[5][0];
		TauLines[ipTFe16].Emis().xIntensity() = 
			TauLines[ipTFe16].Emis().phots()*TauLines[ipTFe16].EnergyErg();

		(*TauLines[ipTFe26].Lo()).Pop() = pops[1];
		(*TauLines[ipTFe26].Hi()).Pop() = pops[5];
		TauLines[ipTFe26].Emis().PopOpc() = (pops[1] - pops[5]);
		TauLines[ipTFe26].Emis().phots() = pops[5]*AulEscp[5][1];
		TauLines[ipTFe26].Emis().xIntensity() = 
			TauLines[ipTFe26].Emis().phots()*TauLines[ipTFe26].EnergyErg();

		(*TauLines[ipTFe34].Lo()).Pop() = pops[2];
		(*TauLines[ipTFe34].Hi()).Pop() = pops[3];
		TauLines[ipTFe34].Emis().PopOpc() = (pops[2] - pops[3]);
		TauLines[ipTFe34].Emis().phots() = pops[3]*AulEscp[3][2];
		TauLines[ipTFe34].Emis().xIntensity() = 
			TauLines[ipTFe34].Emis().phots()*TauLines[ipTFe34].EnergyErg();

		(*TauLines[ipTFe35].Lo()).Pop() = pops[2];
		(*TauLines[ipTFe35].Hi()).Pop() = pops[4];
		TauLines[ipTFe35].Emis().PopOpc() = (pops[2] - pops[4]);
		TauLines[ipTFe35].Emis().phots() = pops[4]*AulEscp[4][2];
		TauLines[ipTFe35].Emis().xIntensity() = 
			TauLines[ipTFe35].Emis().phots()*TauLines[ipTFe35].EnergyErg();

		(*TauLines[ipTFe46].Lo()).Pop() = pops[3];
		(*TauLines[ipTFe46].Hi()).Pop() = pops[5];
		TauLines[ipTFe46].Emis().PopOpc() = (pops[3] - pops[5]);
		TauLines[ipTFe46].Emis().phots() = pops[5]*AulEscp[5][3];
		TauLines[ipTFe46].Emis().xIntensity() = 
			TauLines[ipTFe46].Emis().phots()*TauLines[ipTFe46].EnergyErg();

		(*TauLines[ipTFe56].Lo()).Pop() = pops[4];
		(*TauLines[ipTFe56].Hi()).Pop() = pops[5];
		TauLines[ipTFe56].Emis().PopOpc() = (pops[4] - pops[5]);
		TauLines[ipTFe56].Emis().phots() = pops[5]*AulEscp[5][4];
		TauLines[ipTFe56].Emis().xIntensity() = 
			TauLines[ipTFe56].Emis().phots()*TauLines[ipTFe56].EnergyErg();

		/* Jack's funny FeII lines, data from 
		 * >>refer	fe2	energy	Johansson, S., Brage, T., Leckrone, D.S., Nave, G. &
		 * >>refercon Wahlgren, G.M. 1995, ApJ 446, 361 */
		PutCS(10.,TauLines[ipT191]);
		atom_level2(TauLines[ipT191]);
	}

	{
		/*@-redef@*/
		enum{DEBUG_LOC=false};
		/*@+redef@*/
		if( DEBUG_LOC && iteration > 1 && nzone >=4 )
		{
			fprintf(ioQQQ,"DEBUG2\t%.2e\t%.2e\t%.2e\n",
				phycon.te,
				FeII.Fe2_large_cool , 
				FeII.Fe2_UVsimp_cool );
		}
	}

	return;

}

/*CoolIron - calculation total cooling due to Fe */
void CoolIron(void)
{
	realnum rate;

	DEBUG_ENTRY( "CoolIron()" );

	/* cooling by FeI 24m, 34.2 m */
	/* >>chng 03 nov 15, add these lines */
	/** \todo	2	- ground term is actually a fix level system, the vectors are
	 * created, with pointers ipFe1_54m , ipFe1_111m, must add collision date, use
	 * larger model atom */
	/*>>refer	Fe1	cs	Hollenbach & McKee 89 */
	/* the 24.0 micron line */
	rate = (realnum)(1.2e-7 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*8.0e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		8.0e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	LineConvRate2CS( TauLines[ipFe1_24m]  , rate );

	rate = (realnum)(9.3e-8 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*5.3e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		5.3e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	LineConvRate2CS( TauLines[ipFe1_35m]  , rate );

	rate = (realnum)(1.2e-7 * dense.eden + 
		/* >>chng 05 jul 05, eden to cdsqte */
		/*6.9e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]) / dense.eden);*/
		6.9e-10*pow((phycon.te/100.), 0.17 )*dense.xIonDense[ipHYDROGEN][0]);
	(*(*TauDummy).Hi()).g() = (*TauLines[ipFe1_35m].Hi()).g();
	LineConvRate2CS( *TauDummy  , rate );
	/* this says that line is a dummy, not real one */
	(*(*TauDummy).Hi()).g() = 0.;

	atom_level3(TauLines[ipFe1_24m],TauLines[ipFe1_35m],*TauDummy);

	/* series of FeI lines from Dima Verner's list, each 2-lev atom
	 *
	 * Fe I 3884 */
	MakeCS(TauLines[ipFeI3884]);
	atom_level2(TauLines[ipFeI3884]);

	/* Fe I 3729 */
	MakeCS(TauLines[ipFeI3729]);
	atom_level2(TauLines[ipFeI3729]);

	/* Fe I 3457 */
	MakeCS(TauLines[ipFeI3457]);
	atom_level2(TauLines[ipFeI3457]);

	/* Fe I 3021 */
	MakeCS(TauLines[ipFeI3021]);
	atom_level2(TauLines[ipFeI3021]);

	/* Fe I 2966 */
	MakeCS(TauLines[ipFeI2966]);
	atom_level2(TauLines[ipFeI2966]);

	/* >>chng 05 dec 03, move Fe2 FeII Fe II cooling into separate routine */
	Fe2_cooling();

	/* lump 3p and 3f together; cs=
	 * >>refer	fe3	as	Garstang, R.H., Robb, W.D., Rountree, S.P. 1978, ApJ, 222, 384
	 * A from
	 * >>refer	fe3	as	Garstang, R.H., 1957, Vistas in Astronomy, 1, 268
	 * FE III 5270, is 20.9% of total 
	 * >>chng 05 feb 18, Kevin Blagrave email
	 * average wavelength is 4823 with statistical weight averaging of upper energy level,
	 * as per , change 5th number from 2.948 to 2.984, also photon energy
	 * from 3.78 to 4.12 */

	/* >>chng 05 dec 16, FeIII update by Kevin Blagrave */
	/* FeIII 1122 entire multiplet - atomic data=A's Dima, CS = guess */
	PutCS(25.,TauLines[ipT1122]);
	atom_level2(TauLines[ipT1122]);

	/* call 14 level atom */
	Fe3Lev14();

	/* call 12 level atom */
	Fe4Lev12();

	/* FE V 3892 + 3839, data from Shields */
	CoolHeavy.c3892 = atom_pop2(7.4,25.,5.,0.6,3.7e4,dense.xIonDense[ipIRON][4])*
	  5.11e-12;
	CoolAdd("Fe 5",3892,CoolHeavy.c3892);

	return;
}

/*Fe4Lev12 compute populations and cooling due to 12 level Fe IV ion */
STATIC void Fe4Lev12(void)
{
	const int NLFE4 = 12;
	bool lgZeroPop;
	int nNegPop;
	long int i, 
	  j;
	static bool lgFirst=true;

	double dfe4dt;

	/*static long int **ipdest; */
	static double 
		**AulEscp,
		**col_str,
		**AulDest, 
		depart[NLFE4],
		pops[NLFE4], 
		destroy[NLFE4], 
		create[NLFE4],
		**CollRate,
		**AulPump;

	static const double Fe4A[NLFE4][NLFE4] = {
		{0.,0.,0.,1.e-5,0.,1.368,.89,0.,1.3e-3,1.8e-4,.056,.028},
		{0.,0.,2.6e-8,0.,0.,0.,0.,0.,1.7e-7,0.,0.,0.},
		{0.,0.,0.,0.,3.5e-7,6.4e-10,0.,0.,6.315e-4,0.,6.7e-7,0.},
		{0.,0.,0.,0.,1.1e-6,6.8e-5,8.6e-6,3.4e-10,7.6e-5,1.e-7,5.8e-4,2.8e-4},
		{0.,0.,0.,0.,0.,1.5e-5,1.3e-9,0.,7.6e-4,0.,1.1e-6,6.0e-7},
		{0.,0.,0.,0.,0.,0.,1.1e-5,1.2e-13,.038,9.9e-7,.022,.018},
		{0.,0.,0.,0.,0.,0.,0.,3.7e-5,2.9e-6,.034,3.5e-3,.039},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,.058,3.1e-6,1.4e-3},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.3e-4,3.1e-14},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.9e-19,1.0e-5},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.3e-7},
		{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}
	};

	static const double gfe4[NLFE4]={6.,12.,10.,6.,8.,6.,4.,2.,8.,2.,6.,4.};

	/* excitation energies in Kelvin 
	static double ex[NLFE4]={0.,46395.8,46464.,46475.6,46483.,50725.,
	  50838.,50945.,55796.,55966.,56021.,56025.};*/
	/*>>refer	Fe3	energies	version 3 NIST Atomic Spectra Database */
	/*>>chng 05 dec 17, from Kelvin above to excitation energies in wn */
	static const double excit_wn[NLFE4]={0.,32245.5,32292.8,32301.2,32305.7,35253.8,
	  35333.3,35406.6,38779.4,38896.7,38935.1,38938.2};

	DEBUG_ENTRY( "Fe4Lev12()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		AulPump = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		AulEscp = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE4)*sizeof(double *)));
		for( i=0; i < NLFE4; ++i )
		{
			AulPump[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			AulEscp[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE4)*sizeof(double )));
		}
	}

	/* bail if no Fe+3 */
	if( dense.xIonDense[ipIRON][3] <= 0. )
	{
		fe.Fe4CoolTot = 0.;
		fe.fe40401 = 0.;
		fe.fe42836 = 0.;
		fe.fe42829 = 0.;
		fe.fe42567 = 0.;
		fe.fe41207 = 0.;
		fe.fe41206 = 0.;
		fe.fe41106 = 0.;
		fe.fe41007 = 0.;
		fe.fe41008 = 0.;
		fe.fe40906 = 0.;
		CoolAdd("Fe 4",0,0.);

		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE4 <= LIMLEVELN);
		for( i=0; i < NLFE4; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}
	/* number of levels in model ion */

	/* these are in wavenumbers
	 * data excit_wn/ 0., 32245.5, 32293., 32301.2, 
	 *  1  32306., 35254., 35333., 35407., 38779., 38897., 38935.,
	 *  2  38938./ 
	 * excitation energies in Kelvin */

	/* A's are from Garstang, R.H., MNRAS 118, 572 (1958).
	 * each set is for a lower level indicated by second element in array,
	 * index runs over upper level
	 * A's are saved into arrays as data(up,lo) */

	/* collision strengths from Berrington and Pelan  Ast Ap S 114, 367.
	 * order is cs(low,up) */

	/* all elements are used, and must be set to zero if zero */
	for( i=0; i < NLFE4; i++ )
	{
		create[i] = 0.;
		destroy[i] = 0.;
		for( j=0; j < NLFE4; j++ )
		{
			col_str[j][i] = 0.;
			AulEscp[j][i] = 0.;
			AulDest[j][i] = 0.;
			AulPump[j][i] = 0.;
		}
	}

	/* fill in Einstein As and collision strengths */
	for( long ipHi=1; ipHi < NLFE4; ipHi++ )
	{
		for( long ipLo=0; ipLo < ipHi; ipLo++ )
		{
			AulEscp[ipHi][ipLo] = Fe4A[ipLo][ipHi];
			col_str[ipHi][ipLo] = Fe4_cs(ipLo, ipHi);
		}
	}

	/* leveln itself is well-protected against zero abundances,
	 * low temperatures */

	atom_levelN(NLFE4,
		dense.xIonDense[ipIRON][3],
		gfe4,
		excit_wn,
		'w',
		pops,
		depart,
		&AulEscp ,
		&col_str ,
		&AulDest,
		&AulPump,
		&CollRate,
		create,
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,
		&fe.Fe4CoolTot,
		&dfe4dt,
		"FeIV",
		/* nNegPop positive if negative pops occured, negative if too cold */
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE4 <= LIMLEVELN );
	for( i=0; i < NLFE4; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )
	{
		fprintf( ioQQQ, " fe4levl2 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	CoolAdd("Fe 4",0,fe.Fe4CoolTot);

	thermal.dCooldT += dfe4dt;

	/* three transitions hst can observe
	 * 4 -1 3095.9A and 5 -1 3095.9A */
	fe.fe40401 = (pops[3]*Fe4A[0][3]*(excit_wn[3] - excit_wn[0]) + 
		pops[4]*Fe4A[0][4]*(excit_wn[4] - excit_wn[0]) )*T1CM*BOLTZMANN;

	fe.fe42836 = pops[5]*Fe4A[0][5]*(excit_wn[5] - excit_wn[0])*T1CM*BOLTZMANN;

	fe.fe42829 = pops[6]*Fe4A[0][6]*(excit_wn[5] - excit_wn[0])*T1CM*BOLTZMANN;

	fe.fe42567 = (pops[10]*Fe4A[0][10]*(excit_wn[10] - excit_wn[0]) + 
		pops[11]*Fe4A[0][11]*(excit_wn[10] - excit_wn[0]))*T1CM*BOLTZMANN;

	fe.fe41207 = pops[11]*Fe4A[6][11]*(excit_wn[11] - excit_wn[6])*T1CM*BOLTZMANN;
	fe.fe41206 = pops[11]*Fe4A[5][11]*(excit_wn[11] - excit_wn[5])*T1CM*BOLTZMANN;
	fe.fe41106 = pops[10]*Fe4A[5][10]*(excit_wn[10] - excit_wn[5])*T1CM*BOLTZMANN;
	fe.fe41007 = pops[9]*Fe4A[6][9]*(excit_wn[9] - excit_wn[6])*T1CM*BOLTZMANN;
	fe.fe41008 = pops[9]*Fe4A[7][9]*(excit_wn[9] - excit_wn[7])*T1CM*BOLTZMANN;
	fe.fe40906 = pops[8]*Fe4A[5][8]*(excit_wn[8] - excit_wn[5])*T1CM*BOLTZMANN;
	return;
}

/*Fe3Lev14 compute populations and cooling due to 14 level Fe III ion 
 * >>chng 05 dec 17, code provided by Kevin Blagrave and described in
 *>>refer	Fe3	model	Blagrave, K.P.M., Martin, P.G. & Baldwin, J.A. 
 *>>refercon 2006, ApJ, 644, 1006B (astro-ph 0603167) */
STATIC void Fe3Lev14(void)
{
	bool lgZeroPop;
	int nNegPop;
	long int i,
		j;
	static bool lgFirst=true;

	double dfe3dt;

	long int ihi , ilo;
	static double
		*depart,
		*pops,
		*destroy,
		*create ,
		**AulDest,
		**CollRate,
		**AulPump,
		**AulNet,
		**col_str;

	/* statistical weights for the n levels */
	static double gfe3[NLFE3]={9.,7.,5.,3.,1.,5.,13.,11.,9.,3.,1.,9.,7.,5.};

	/*refer Fe3	energies	NIST version 3 Atomic Spectra Database */
	/* from smallest to largest energy (not in multiplet groupings) */

	/* energy in wavenumbers */
	static double excit_wn[NLFE3]={
		0.0    ,   436.2,   738.9,   932.4,  1027.3,
		19404.8, 20051.1, 20300.8, 20481.9, 20688.4,
		21208.5, 21462.2, 21699.9, 21857.2 };

	DEBUG_ENTRY( "Fe3Lev14()" );

	if( lgFirst )
	{
		/* will never do this again */
		lgFirst = false;
		/* allocate the 1D arrays*/
		/* create space for the 2D arrays */
		depart = ((double *)MALLOC((NLFE3)*sizeof(double)));
		pops = ((double *)MALLOC((NLFE3)*sizeof(double)));
		destroy = ((double *)MALLOC((NLFE3)*sizeof(double)));
		create = ((double *)MALLOC((NLFE3)*sizeof(double)));
		/* now the 2-d arrays */
		fe.Fe3_wl = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		fe.Fe3_emiss = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulNet = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		col_str = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulPump = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		CollRate = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		AulDest = ((double **)MALLOC((NLFE3)*sizeof(double *)));
		for( i=0; i < NLFE3; ++i )
		{
			fe.Fe3_wl[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			fe.Fe3_emiss[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulNet[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			col_str[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulPump[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			CollRate[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
			AulDest[i] = ((double *)MALLOC((NLFE3)*sizeof(double )));
		}

		/* set some to constant values after zeroing out */
		for( i=0; i < NLFE3; ++i )
		{
			create[i] = 0.;
			destroy[i] = 0.;
			for( j=0; j < NLFE3; ++j )
			{
				AulNet[i][j] = 0.;
				col_str[i][j] = 0.;
				CollRate[i][j] = 0.;
				AulDest[i][j] = 0.;
				AulPump[i][j] = 0.;
				fe.Fe3_wl[i][j] = 0.;
				fe.Fe3_emiss[i][j] = 0.;
			}
		}
		/* calculates wavelengths of transitions */
		/* dividing by RefIndex converts the vacuum wavelength to air wavelength */
		for( ihi=1; ihi < NLFE3; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe3_wl[ihi][ilo] = 1e8/(excit_wn[ihi]-excit_wn[ilo]) / 
					RefIndex( (excit_wn[ihi]-excit_wn[ilo]) );
			}
		}

		/* assume FeIII is optically thin - just use As as net escape */
		/*>>refer	Fe3	as	Quinet, P., 1996, A&AS, 116, 573      */
		AulNet[1][0] = 2.8e-3;
		AulNet[7][0] = 4.9e-6;
		AulNet[8][0] = 5.7e-3;
		AulNet[11][0] = 4.5e-1;
		AulNet[12][0] = 4.2e-2;

		AulNet[2][1] = 1.8e-3;
		AulNet[5][1] = 4.2e-1;
		AulNet[8][1] = 1.0e-3;
		AulNet[11][1] = 8.4e-2;
		AulNet[12][1] = 2.5e-1;
		AulNet[13][1] = 2.7e-2;

		AulNet[3][2] = 7.0e-4;
		AulNet[5][2] = 5.1e-5;
		AulNet[9][2] = 5.4e-1;
		AulNet[12][2] = 8.5e-2;
		AulNet[13][2] = 9.8e-2;

		AulNet[4][3] = 1.4e-4;
		AulNet[5][3] = 3.9e-2;
		AulNet[9][3] = 4.1e-5;
		AulNet[10][3] = 7.0e-1;
		AulNet[13][3] = 4.7e-2;

		AulNet[9][4] = 9.3e-2;

		AulNet[9][5] = 4.7e-2;
		AulNet[12][5] = 2.5e-6;
		AulNet[13][5] = 1.7e-5;

		AulNet[7][6] = 2.7e-4;

		AulNet[8][7] = 1.2e-4;
		AulNet[11][7] = 6.6e-4;

		AulNet[11][8] = 1.6e-3;
		AulNet[12][8] = 7.8e-4;

		AulNet[10][9] = 8.4e-3;
		AulNet[13][9] = 2.8e-7;

		AulNet[12][11] = 3.0e-4;

		AulNet[13][12] = 1.4e-4;

		for( int ipHi = 1; ipHi < NLFE3; ipHi++)
		{
			for( int ipLo = 0; ipLo < ipHi; ipLo++)
			{
				col_str[ipHi][ipLo] = Fe3_cs(ipLo,ipHi);
			}
		}
	}

	/* bail if no ions */
	if( dense.xIonDense[ipIRON][2] <= 0. )
	{
		CoolAdd("Fe 3",0,0.);

		fe.Fe3CoolTot = 0.;   
		for( ihi=1; ihi < NLFE3; ++ihi )
		{
			for( ilo=0; ilo < ihi; ++ilo )
			{
				fe.Fe3_emiss[ihi][ilo] = 0.;
			}
		}
		/* level populations */
		/* LIMLEVELN is the dimension of the atoms vectors */
		ASSERT( NLFE3 <= LIMLEVELN);
		for( i=0; i < NLFE3; i++ )
		{
			atoms.PopLevels[i] = 0.;
			atoms.DepLTELevels[i] = 1.;
		}
		return;
	}

	/* nNegPop positive if negative pops occurred, negative if too cold */
	atom_levelN(
		/* number of levels */
		NLFE3,
		/* the abundance of the ion, cm-3 */
		dense.xIonDense[ipIRON][2],
		/* the statistical weights */
		gfe3,
		/* the excitation energies */
		excit_wn,
		'w',
		/* the derived populations - cm-3 */
		pops,
		/* the derived departure coefficients */
		depart,
		/* the net emission rate, Aul * escape prob */
		&AulNet ,
		/* the collision strengths */
		&col_str ,
		/* A * destruction prob */
		&AulDest,
		/* pumping rate */
		&AulPump,
		/* collision rate, s-1, must defined if no collision strengths */
		&CollRate,
		/* creation vector */
		create,
		/* destruction vector */
		destroy,
		/* say atom_levelN should evaluate coll rates from cs */
		false,   
		&fe.Fe3CoolTot,
		&dfe3dt,
		"Fe 3",
		&nNegPop,
		&lgZeroPop,
		false );

	/* LIMLEVELN is the dim of the PopLevels vector */
	ASSERT( NLFE3 <= LIMLEVELN );
	for( i=0; i < NLFE3; ++i )
	{
		atoms.PopLevels[i] = pops[i];
		atoms.DepLTELevels[i] = depart[i];
	}

	if( nNegPop > 0 )
	{
		fprintf( ioQQQ, " Fe3Lev14 found negative populations\n" );
		ShowMe();
		cdEXIT(EXIT_FAILURE);
	}

	/* add cooling then its derivative */
	CoolAdd("Fe 3",0,fe.Fe3CoolTot);
	/* derivative of cooling */
	thermal.dCooldT += dfe3dt;

	/* find emission in each line */
	for( ihi=1; ihi < NLFE3; ++ihi )
	{
		for( ilo=0; ilo < ihi; ++ilo )
		{
			/* emission in these lines */
			fe.Fe3_emiss[ihi][ilo] = pops[ihi]*AulNet[ihi][ilo]*(excit_wn[ihi] - excit_wn[ilo])*T1CM*BOLTZMANN;
		}
	}
	return;
}

double Fe3_cs(long ipLo, long ipHi)
{
	static double col_str[14][14];

	static double lgOneTimeMustInit=true;
	if( lgOneTimeMustInit )
	{
		lgOneTimeMustInit = false;
		/* collision strengths at log T 4 */
		/** \todo	2	put in temperature dependence */
		/*>>refer	Fe3	cs	Zhang, H.  1996, 119, 523 */
		col_str[1][0] = 2.92;
		col_str[2][0] = 1.24;
		col_str[3][0] = 0.595;
		col_str[4][0] = 0.180;
		col_str[5][0] = 0.580;
		col_str[6][0] = 1.34;
		col_str[7][0] = 0.489;
		col_str[8][0] = 0.0926;
		col_str[9][0] = 0.165;
		col_str[10][0] = 0.0213;
		col_str[11][0] = 1.07;
		col_str[12][0] = 0.435;
		col_str[13][0] = 0.157;
	
		col_str[2][1] = 2.06;
		col_str[3][1] = 0.799;
		col_str[4][1] = 0.225;
		col_str[5][1] = 0.335;
		col_str[6][1] = 0.555;
		col_str[7][1] = 0.609;
		col_str[8][1] = 0.367;
		col_str[9][1] = 0.195;
		col_str[10][1] = 0.0698;
		col_str[11][1] = 0.538;
		col_str[12][1] = 0.484;
		col_str[13][1] = 0.285;
	
		col_str[3][2] = 1.29;
		col_str[4][2] = 0.312;
		col_str[5][2] = 0.173;
		col_str[6][2] = 0.178;
		col_str[7][2] = 0.430;
		col_str[8][2] = 0.486;
		col_str[9][2] = 0.179;
		col_str[10][2] = 0.0741;
		col_str[11][2] = 0.249;
		col_str[12][2] = 0.362;
		col_str[13][2] = 0.324;
	
		col_str[4][3] = 0.493;
		col_str[5][3] = 0.0767;
		col_str[6][3] = 0.0348;
		col_str[7][3] = 0.223;
		col_str[8][3] = 0.401;
		col_str[9][3] = 0.126;
		col_str[10][3] = 0.0528;
		col_str[11][3] = 0.101;
		col_str[12][3] = 0.207;
		col_str[13][3] = 0.253;
	
		col_str[5][4] = 0.0211;
		col_str[6][4] = 0.00122;
		col_str[7][4] = 0.0653;
		col_str[8][4] = 0.154;
		col_str[9][4] = 0.0453;
		col_str[10][4] = 0.0189;
		col_str[11][4] = 0.0265;
		col_str[12][4] = 0.0654;
		col_str[13][4] = 0.0950;
	
		col_str[6][5] = 0.403;
		col_str[7][5] = 0.213;
		col_str[8][5] = 0.0939;
		col_str[9][5] = 1.10;
		col_str[10][5] = 0.282;
		col_str[11][5] = 0.942;
		col_str[12][5] = 0.768;
		col_str[13][5] = 0.579;
	
		col_str[7][6] = 2.84; /* 10-9 */
		col_str[8][6] = 0.379; /* 11-9 */
		col_str[9][6] = 0.0876;  /* 7-9 */
		col_str[10][6] = 0.00807; /* 8-9 */
		col_str[11][6] = 1.85; /* 12-9 */
		col_str[12][6] = 0.667; /* 13-9 */
		col_str[13][6] = 0.0905; /* 14-9 */
	
		col_str[8][7] = 3.07; /* 11-10 */
		col_str[9][7] = 0.167;   /* 7-10 */
		col_str[10][7] = 0.0526;  /* 8-10 */
		col_str[11][7] = 0.814; /* 12-10 */
		col_str[12][7] = 0.837; /* 13-10 */
		col_str[13][7] = 0.626; /* 14-10 */
	
		col_str[9][8] = 0.181; /* 7-11 */
		col_str[10][8] = 0.0854; /* 8-11 */
		col_str[11][8] = 0.180; /* 12-11 */
		col_str[12][8] = 0.778; /* 13-11 */
		col_str[13][8] = 0.941; /* 14-11 */
	
		col_str[10][9] = 0.377; /* 8-7 */
		col_str[11][9] = 0.603; /* 12-7 */
		col_str[12][9] = 0.472; /* 13-7 */
		col_str[13][9] = 0.302; /* 14-7 */
	
		col_str[11][10] = 0.216; /* 12-8 */
		col_str[12][10] = 0.137; /* 13-8 */
		col_str[13][10] = 0.106; /* 14-8 */
	
		col_str[12][11] = 1.25;
		col_str[13][11] = 0.292;
	
		col_str[13][12] = 1.10;
	}

	ASSERT( ipHi > ipLo );
	double CollisionStrength = col_str[ipHi][ipLo];
	ASSERT( CollisionStrength >0. );

	return( CollisionStrength );
}

double Fe4_cs(long ipLo, long ipHi)
{
	const int NLFE4=12;

	/* collision strengths from Berrington and Pelan  Ast Ap S 114, 367.
	 * order is cs(low,up) */
	static const double Fe4CS[NLFE4][NLFE4] = {
			{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.98,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.8167,3.72,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.49,0.0475,0.330,0.,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.6533,0.473,2.26,1.64,0.,0.,0.,0.,0.,0.,0.,0.},
			{0.45,0.686,0.446,0.106,0.254,0.,0.,0.,0.,0.,0.,0.},
			{0.30,0.392,0.152,0.269,0.199,0.605,0.,0.,0.,0.,0.,0.},
			{0.15,0.0207,0.190,0.0857,0.166,0.195,0.327,0.,0.,0.,0.,0.},
			{0.512,1.23,0.733,0.174,0.398,0.623,0.335,0.102,0.,0.,0.,0.},
			{0.128,0.0583,0.185,0.200,0.188,0.0835,0.127,0.0498,0.0787,0.,0.,0.},
			{0.384,0.578,0.534,0.363,0.417,0.396,0.210,0.171,0.810,0.101,0.,0.},
			{0.256,0.234,0.306,0.318,0.403,0.209,0.195,0.112,0.195,0.458,0.727,0.}
	};

	ASSERT( ipHi > ipLo );
	double CollisionStrength = Fe4CS[ipHi][ipLo];
	ASSERT( CollisionStrength >0. );

	return( CollisionStrength );
}

double Fe5_cs(long ipLo, long ipHi)
{
	const int NLFE5 = 14;
	static double col_str[NLFE5][NLFE5];
	/*>>refer	Fe5	cs	Shields ApJ 219, 559. */


	static double lgOneTimeMustInit=true;
	if( lgOneTimeMustInit )
	{
		lgOneTimeMustInit = false;
		for( int i = 0;i < NLFE5;i++)
		{
			for( int j = 0;j < NLFE5;j++)
			{
				col_str[i][j] = 1.;
			}
		}

		col_str[10][3] = 1.4; //3896A
		col_str[7][2] = 1.1; //4072A
		col_str[13][4] = 3.7; //3892A
		col_str[12][3] = 3.7; //3839A
		col_str[11][2] = 2.0; //3795A
	}

	ASSERT( ipHi > ipLo );
	double CollisionStrength = col_str[ipHi][ipLo];
	ASSERT( CollisionStrength >0. );

	return( CollisionStrength );
}
