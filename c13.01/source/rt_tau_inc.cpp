/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_tau_inc increment optical depths once per zone, called after radius_increment */
#include "cddefines.h"
#include "taulines.h"
#include "iso.h"
#include "rfield.h"
#include "trace.h"
#include "dense.h"
#include "hyperfine.h"
#include "wind.h"
#include "prt.h"
#include "conv.h"
#include "h2.h"
#include "mole.h"
#include "hmi.h"
#include "opacity.h"
#include "cooling.h"
#include "thermal.h"
#include "radius.h"
#include "atomfeii.h"
#include "rt.h"
#include "doppvel.h"
#include "mole.h"

/*RT_tau_inc increment optical depths once per zone, called after radius_increment */
void RT_tau_inc(void)
{

	long int i, 
	  ipHi, 
	  ipLo;

	DEBUG_ENTRY( "RT_tau_inc()" );

	if( trace.lgTrace )
	{
		fprintf( ioQQQ, " RT_tau_inc called.\n" );
	}

	/* call RT_line_all one last time in this zone, to get fine opacities defined */
	ASSERT( !conv.lgFirstSweepThisZone );
	conv.lgLastSweepThisZone = true;
	RT_line_all( );
	/* this may have updated some escape/destruction rates - force update
	 * to all cooling lines */
	CoolEvaluate( &thermal.ctot );

	if( nzone <=1 )
	{
		opac.telec = (realnum)(radius.drad_x_fillfac*dense.eden*6.65e-25);
		opac.thmin = (realnum)(radius.drad_x_fillfac*findspecieslocal("H-")->den*3.9e-17*
			(1. - rfield.ContBoltz[hmi.iphmin-1]/ hmi.hmidep));
	}
	else
	{
		opac.telec += (realnum)(radius.drad_x_fillfac*dense.eden*6.65e-25);
		opac.thmin += (realnum)(radius.drad_x_fillfac*findspecieslocal("H-")->den*3.9e-17*
			(1. - rfield.ContBoltz[hmi.iphmin-1]/ hmi.hmidep));
	}

	/* prevent maser runaway */
	rt.dTauMase = 0;
	rt.mas_species = 0;
	rt.mas_ion = 0;
	rt.mas_hi = 0;
	rt.mas_lo = 0;

	/* all lines in iso sequences */
	for( long ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		for( long nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			/* this is the parent ion, for HI lines, is 1, 
			 * for element He is 1 for He-like (HeI) and 2 for H-like (HeII) */
			int ion = nelem+1-ipISO;
			/* do not evaluate in case where trivial parent ion */
			if( ion <=dense.IonHigh[nelem] && dense.xIonDense[nelem][ion] > dense.density_low_limit )
			{
				if( iso_ctrl.lgDielRecom[ipISO] )
				{
					// SatelliteLines are indexed by lower level
					for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local; ipLo++ )
					{
						RT_line_one_tauinc(SatelliteLines[ipISO][nelem][ipSatelliteLines[ipISO][nelem][ipLo]], ipISO, nelem, -1, ipLo, 
							GetDopplerWidth(dense.AtomicWeight[nelem]) );
					}
				}

				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
				{
					for( ipLo=0; ipLo < ipHi; ipLo++ )
					{
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() <= 0 )
							continue;

						/* actually do the work */
						RT_line_one_tauinc(iso_sp[ipISO][nelem].trans(ipHi,ipLo), ipISO, nelem, ipHi, ipLo, 
							GetDopplerWidth(dense.AtomicWeight[nelem]) );
					}
				}
				ipLo = 0;
				/* these are the extra Lyman lines */
				for( ipHi=iso_sp[ipISO][nelem].st[iso_sp[ipISO][nelem].numLevels_local-1].n()+1; ipHi < iso_ctrl.nLyman[ipISO]; ipHi++ )
				{
					TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi];
					(*tr).Emis().PopOpc() = iso_sp[ipISO][nelem].st[0].Pop();

					/* actually do the work */
					RT_line_one_tauinc(*tr, -1 ,ipISO, nelem, ipHi,
						GetDopplerWidth(dense.AtomicWeight[nelem]) );
				}
			}
		}
	}

	/* increment optical depths for all heavy element lines
	 * same routine does wind and static,
	 * does not start from 0 since first line is dummy */
	for( i=1; i <= nLevel1; i++ )
	{
		RT_line_one_tauinc(TauLines[i], -2, -2, -2, i, GetDopplerWidth(dense.AtomicWeight[(*TauLines[i].Hi()).nelem()-1]) );
	}

	/* all lines in cooling with g-bar */
	for( i=0; i < nWindLine; i++ )
	{
		/* do not include H-like or He-like in the level two lines since
		 * these are already counted in iso sequences */
		if( (*TauLine2[i].Hi()).IonStg() < (*TauLine2[i].Hi()).nelem()+1-NISO )
		{
			RT_line_one_tauinc(TauLine2[i], -3, -3, -3, i, GetDopplerWidth(dense.AtomicWeight[(*TauLine2[i].Hi()).nelem()-1]) );
		}
	}

	/* the block of inner shell lines */
	for( i=0; i < nUTA; i++ )
	{
		/* populations have not been set */
		UTALines[i].Emis().PopOpc() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
		(*UTALines[i].Lo()).Pop() = dense.xIonDense[(*UTALines[i].Hi()).nelem()-1][(*UTALines[i].Hi()).IonStg()-1];
		(*UTALines[i].Hi()).Pop() = 0.;
		RT_line_one_tauinc(UTALines[i], -4 , -4 , -4 , i, GetDopplerWidth(dense.AtomicWeight[(*UTALines[i].Hi()).nelem()-1]) );
	}

	/* all hyper fine structure lines  */
	for( i=0; i < nHFLines; i++ )
	{
		/* remember current gas-phase abundances */
		realnum save = dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1];

		/* bail if no abundance */
		if( save<=0. ) continue;

		/* set gas-phase abundance to total times isotope ratio */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] *= hyperfine.HFLabundance[i];

		RT_line_one_tauinc(HFLines[i] , -5 , -5 , -5 , i, GetDopplerWidth(dense.AtomicWeight[(*HFLines[i].Hi()).nelem()-1]) );

		/* put the correct gas-phase abundance back in the array */
		dense.xIonDense[(*HFLines[i].Hi()).nelem()-1][(*HFLines[i].Hi()).IonStg()-1] = save;
	}

	/* do large FeII atom if this is enabled */
	FeII_RT_TauInc();

	/* increment optical depth for the H2 molecule */
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		(*diatom)->H2_RT_tau_inc();

	/* database Lines*/
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			realnum DopplerWidth = GetDopplerWidth( dBaseSpecies[ipSpecies].fmolweight );
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{	
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
					continue;
				int ipLo = (*tr).ipLo();

				RT_line_one_tauinc( *tr, -10, ipSpecies, ipHi, ipLo, DopplerWidth );
			}
		}
	}
	
	/* following is for static atmosphere */
	if( wind.lgStatic() )
	{
		/* iron fe feii fe2  - overlapping feii lines */
		t_fe2ovr_la::Inst().tau_inc();
	}

	if( trace.lgTrace && trace.lgOptcBug )
	{
		fprintf( ioQQQ, " RT_tau_inc updated optical depths:\n" );
		prtmet();
	}

	if( trace.lgTrace )
		fprintf( ioQQQ, " RT_tau_inc returns.\n" );

	return;
}
