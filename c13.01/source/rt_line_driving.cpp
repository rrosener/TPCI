/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*RT_line_driving derive radiative acceleration due to line absorption of incident continuum,
 * return value is line radiative acceleration */
#include "cddefines.h"
#include "physconst.h"
#include "iso.h"
#include "dense.h"
#include "taulines.h"
#include "h2.h"
#include "atomfeii.h"
#include "rt.h"

/*RT_line_driving derive radiative acceleration due to line absorption of incident continuum,
 * return value is line radiative acceleration */
double RT_line_driving(void)
{
	long int i, 
	  ipHi, 
	  nelem, 
	  ipLo,
	  ipISO;

	double AllHeavy, 
	  AllRest, 
	  OneLine, 
	  fe2drive, 
	  forlin_v, 
	  h2drive,
	  accel_iso[NISO];

	/* following used for debugging */
	/* double 
	  RestMax, 
	  HeavMax, 
	  hydromax;
	  long int 
	  ipRestMax, 
	  ihmax; */

	DEBUG_ENTRY( "RT_line_driving()" );

	/* this function finds the total rate the gas absorbs energy
	 * this result is divided by the calling routine to find the
	 * momentum absorbed by the gas, and eventually the radiative acceleration
	 *
	 * the energy absorbed by the line is
	 * Abundance * energy * A *(g_up/g_lo) * occnum * escape prob
	 * where occnum is the photon occupation number, and the g's are
	 * the ratios of statistical weights */

	/* total energy absorbed in this zone, per cubic cm
	 * do hydrogen first */

	for( ipISO=ipH_LIKE; ipISO<NISO; ++ipISO )
	{
		accel_iso[ipISO] = 0;
		for( nelem=ipISO; nelem < LIMELM; nelem++ )
		{
			if( (dense.IonHigh[nelem] >= nelem + 1-ipISO)  )
			{
				for( ipHi=1; ipHi < iso_sp[ipISO][nelem].numLevels_local; ipHi++ )
				{
					/* do not put in highest level since its not real */
					for( ipLo=0; ipLo < ipHi - 1; ipLo++ )
					{
						/* do not include bogus lines */
						if( iso_sp[ipISO][nelem].trans(ipHi,ipLo).ipCont() > 0 )
						{
							OneLine = iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().pump()*
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).EnergyErg()*
							iso_sp[ipISO][nelem].trans(ipHi,ipLo).Emis().PopOpc();

							accel_iso[ipISO] += OneLine;
						}
					}
				}
			
				if( iso_ctrl.lgDielRecom[ipISO] )
				{	
					// SatelliteLines are indexed by lower level, summed over satellite levels
					for( ipLo=0; ipLo < iso_sp[ipISO][nelem].numLevels_local; ipLo++ )
					{
						/* do not include bogus lines */
						TransitionList::iterator tr = SatelliteLines[ipISO][nelem].begin()+ipSatelliteLines[ipISO][nelem][ipLo];
						if((*tr).ipCont() > 0 )
						{
							OneLine = (*tr).Emis().pump()*
								(*tr).EnergyErg()*
								(*tr).Emis().PopOpc();

							accel_iso[ipISO] += OneLine;
						}
					}
				}

				for( ipHi=iso_sp[ipISO][nelem].st[iso_sp[ipISO][nelem].numLevels_local-1].n()+1; ipHi < iso_ctrl.nLyman[ipISO]; ipHi++ )
				{
					/* do not include bogus lines */
					TransitionList::iterator tr = ExtraLymanLines[ipISO][nelem].begin()+ipExtraLymanLines[ipISO][nelem][ipHi];
					if( (*tr).ipCont() > 0 )
					{
						OneLine = (*tr).Emis().pump()*
							(*tr).EnergyErg()*
							(*tr).Emis().PopOpc();

						accel_iso[ipISO] += OneLine;
					}

				}
			}
		}
	}

	/* all heavy element lines in calculation of cooling 
	 * these are the level 1 lines */
	AllHeavy = 0.;
	for( i=1; i <= nLevel1; i++ )
	{
		OneLine = 
		  TauLines[i].Emis().pump()*
		  TauLines[i].EnergyErg()*
		  TauLines[i].Emis().PopOpc();
		AllHeavy += OneLine;
	}

	/* all heavy element lines treated with g-bar 
	 * these are the level 2 lines, f should be ok */
	AllRest = 0.;
	for( i=0; i < nWindLine; i++ )
	{
		OneLine = 
			TauLine2[i].Emis().pump()*
			TauLine2[i].EnergyErg()*
			TauLine2[i].Emis().PopOpc();
		AllRest += OneLine;
	}
	for( i=0; i < nUTA; i++ )
	{
		OneLine = 
			UTALines[i].Emis().pump()*
			UTALines[i].EnergyErg()*
			UTALines[i].Emis().PopOpc();
		AllRest += OneLine;
	}
	for( i=0; i < nHFLines; i++ )
	{
		OneLine = 
			HFLines[i].Emis().pump()*
			HFLines[i].EnergyErg()*
			HFLines[i].Emis().PopOpc();
		AllRest += OneLine;
	}
        
	for( long ipSpecies=0; ipSpecies<nSpecies; ipSpecies++ )
	{
		if( dBaseSpecies[ipSpecies].lgActive )
		{
			for (TransitionList::iterator tr=dBaseTrans[ipSpecies].begin(); 
				  tr != dBaseTrans[ipSpecies].end(); ++tr)
			{	
				int ipHi = (*tr).ipHi();
				if (ipHi >= dBaseSpecies[ipSpecies].numLevels_local || (*tr).ipCont() <= 0)
					continue;
				OneLine = (*tr).EnergyErg()*(*tr).Emis().pump()*(*tr).Emis().PopOpc();
				AllRest += OneLine;
			}
		}
	}

	/* the H2 molecule */
	h2drive = 0.;
	for( diatom_iter diatom = diatoms.begin(); diatom != diatoms.end(); ++diatom )
		h2drive += (*diatom)->H2_Accel();

	/* The large model FeII atom */
	fe2drive = 0.;
	FeIIAccel(&fe2drive);

	forlin_v = AllHeavy + accel_iso[ipH_LIKE] + accel_iso[ipHE_LIKE] + 
		fe2drive + h2drive + AllRest;

	/*fprintf(ioQQQ," wind te %e %e %e %e %e\n", 	
		AllHeavy , HydroAccel , fe2drive , he1l , AllRest );*/
	return( forlin_v );
}
