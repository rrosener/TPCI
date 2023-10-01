/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/* mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
/* H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
 * the full sim, called by H2_LevelPops*/
/* H2_Solomon_rate find rates between H2s and H2g and other levels,
 * for eventual use in the chemistry */
/* gs_rate evaluate rates between ground and star states */
/* H2_He_coll_init initialize H2 - He collision data set
 * H2_He_coll interpolate on h2 - He collision data set to return rate at temp*/
#include "cddefines.h" 
#include "phycon.h" 
#include "hmi.h" 
#include "taulines.h" 
#include "h2.h" 
#include "h2_priv.h" 
#include "rfield.h"
#include "thermal.h"
#include "gammas.h"
#include "ionbal.h"

/*H2_Solomon_rate find rates between H2s and H2g and other levels,
 * for eventual use in the chemistry */
void diatomics::H2_Solomon_rate( void )
{
	DEBUG_ENTRY( "H2_Solomon_rate()" );

	/* iElecLo will always be X in this routine */

	/* find rate (s-1) h2 dissociation into X continuum by Solomon process and
	 * assign to the TH85 g and s states 
	 * these will go back into the chemistry network */

	/* rates [s-1] for dissociation from s or g, into electronic excited states  
	 * followed by dissociation */
	Solomon_dissoc_rate_g = 0.;
	Solomon_dissoc_rate_s = 0.;

	/* these are used in a print statement - are they needed? */
	Solomon_elec_decay_g = 0.;
	Solomon_elec_decay_s = 0.;

	/* at this point we have already evaluated the sum of the radiative rates out
	 * of the electronic excited states - this is H2_rad_rate_out[electronic][vib][rot] 
	 * and this includes both decays into the continuum and bound states of X */

	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		qList::iterator Hi = (*tr).Hi() ;	
		if( (*Hi).n() < 1 ) continue;
		double factor = (double)H2_dissprob[(*Hi).n()][(*Hi).v()][(*Hi).J()]/H2_rad_rate_out[(*Hi).n()][(*Hi).v()][(*Hi).J()];
		/* this is the rate [cm-3 s-1] that mole goes from 
		 * lower level into electronic excited states then 
		 * into continuum */
		double rate_up_cont = (*(*tr).Lo()).Pop() * (*tr).Emis().pump() * factor;

		/* rate electronic state decays into H2g */
		double elec_decay = (*(*tr).Hi()).Pop() * (*tr).Emis().Aul() *
			((*tr).Emis().Pesc() + (*tr).Emis().Pdest() + (*tr).Emis().Pelec_esc());
		if( (*(*tr).Lo()).energy().WN() > ENERGY_H2_STAR && hmi.lgLeiden_Keep_ipMH2s )
		{
			/* this is H2g up to excited then to continuum - 
		 	 * cm-3 s-1 at this point */
			Solomon_dissoc_rate_s += rate_up_cont;
			/* rate electronic state decays into H2g */
			Solomon_elec_decay_s += elec_decay;
		}
		else
		{
			/* this is H2g up to excited then to continuum - 
			 * cm-3 s-1 at this point */
			Solomon_dissoc_rate_g += rate_up_cont;
			/* rate electronic state decays into H2g */
			Solomon_elec_decay_g += elec_decay;
		}
	}

	double H2_sum_excit_elec_den = GetExcitedElecDensity();

	/* at this point units of Solomon_elec_decay_g, H2s are cm-3 s-1 
	 * since populations are included -
	 * div by pops to get actual dissocation rate, s-1 */
	if( *dense_total > SMALLFLOAT )
	{
		Solomon_elec_decay_g /= SDIV( H2_sum_excit_elec_den );
		Solomon_elec_decay_s /= SDIV( H2_sum_excit_elec_den );

		/* will be used for H2s-> H + H */
		Solomon_dissoc_rate_s /= SDIV(H2_den_s);

		/* will be used for H2g-> H + H */
		Solomon_dissoc_rate_g /= SDIV(H2_den_g);

	}
	else
	{
		Solomon_dissoc_rate_s = 0.; 
		Solomon_dissoc_rate_g = 0.;	

	}
	/*fprintf(ioQQQ,"DEBUG H2 new %.2e %.2e %.2e %.2e \n",
		Solomon_elec_decay_g ,
		Solomon_elec_decay_s ,
		Solomon_dissoc_rate_s,
		Solomon_dissoc_rate_g );*/

	return;
}

/* gs_rate evaluate rate between ground and star states */
double diatomics::gs_rate( void )
{
	DEBUG_ENTRY( "diatomics::gs_rate()" );

	/* rate g goes to s */
	double ground_to_star_rate = 0.;

	/* loop over all levels in H2g */
	for( long ipLoX=0; ipLoX < nEner_H2_ground; ++ipLoX )
	{
		/* now find all pumps up to electronic excited states */
		/* sum over all electronic states, finding dissociation rate */
		for( long iElecHi=1; iElecHi<n_elec_states; ++iElecHi )
		{
			for( long iVibHi=0; iVibHi<=nVib_hi[iElecHi]; ++iVibHi )
			{
				for( long iRotHi=Jlowest[iElecHi]; iRotHi<=nRot_hi[iElecHi][iVibHi]; ++iRotHi )
				{
					long ipHi = ipEnergySort[iElecHi][iVibHi][iRotHi];
					if( lgH2_radiative[ipHi][ipLoX] )
					{
						/* this is the rate [cm-3 s-1] that mole goes from 
						 * lower level into electronic excited states then 
						 * into continuum */
						double rate_up_cont = 
							states[ipLoX].Pop() *
							trans[ ipTransitionSort[ipHi][ipLoX] ].Emis().pump();
						double decay_star = H2_rad_rate_out[iElecHi][iVibHi][iRotHi] - H2_dissprob[iElecHi][iVibHi][iRotHi];
						/* loop over all other levels in H2g, subtracting 
						 * their rate - remainder is rate into star, this is 
						 * usually only a few levels */
						for( long ipOther=0; ipOther < nEner_H2_ground; ++ipOther )
						{
							if( lgH2_radiative[ipHi][ipOther] ) 
							{
								EmissionProxy em = trans[ ipTransitionSort[ipHi][ipOther] ].Emis();
								decay_star -= em.Aul() * ( em.Pesc() + em.Pdest() + em.Pelec_esc() );
							}
						}
						/* MAX because may underflow to negative numbers is rates very large 
						 * this is fraction that returns to H2s */
						decay_star = MAX2(0., decay_star)/SDIV(H2_rad_rate_out[iElecHi][iVibHi][iRotHi]);
						ground_to_star_rate += rate_up_cont*decay_star;

					}/* end if line exists */
				}/* end loop rot electronic excited */
			}/* end loop vib electronic excited */
		}/* end loop electronic electronic excited */
	}

	/* at this point units are cm-3 s-1 - convert to rate s-1 */
	ground_to_star_rate /= SDIV( H2_den_g );
	return ground_to_star_rate;
}

long diatomics::OpacityCreate( double *stack )
{
	DEBUG_ENTRY( "diatomics::OpacityCreate()" );

	ASSERT( photoion_opacity_fun != NULL );

	for( long i=ip_photo_opac_thresh-1; i < rfield.nupper; i++ )
	{
		/* >>chng 05 nov 24, add H2 photoionization cross section */
		//opac.OpacStack[i-ip_photo_opac_thresh + ip_photo_opac_offset] =
		stack[ i - ip_photo_opac_thresh + ip_photo_opac_offset ] =
			photoion_opacity_fun( rfield.AnuOrg[i] );
			//Yan_H2_CS(rfield.AnuOrg[i]);
	}
         
	//opac.nOpacTot += rfield.nupper - ip_photo_opac_thresh + 1;
	return rfield.nupper - ip_photo_opac_thresh + 1;
}

/* H2_zero_pops_too_low - zero out some H2 variables if we decide not to compute
 * the full sim, called by H2_LevelPops*/
void diatomics::H2_zero_pops_too_low( void )
{
	DEBUG_ENTRY( "H2_zero_pops_too_low()" );

	// zero populations
	fill_n( pops_per_elec, N_ELEC, 0. );
	pops_per_vib.zero();
	for( qList::iterator st = states.begin(); st != states.end(); ++st )
	{
		double pop = H2_populations_LTE[ (*st).n() ][ (*st).v() ][ (*st).J() ] * (*dense_total);
		H2_old_populations[ (*st).n() ][ (*st).v() ][ (*st).J() ] = pop;
		(*st).Pop() = pop;
	}


	for( TransitionList::iterator tr = trans.begin(); tr != rad_end; ++tr )
	{
		/* population of lower level with correction for stim emission */
		(*tr).Emis().PopOpc() = (*(*tr).Lo()).Pop() - (*(*tr).Hi()).Pop() * (*(*tr).Lo()).g() / (*(*tr).Hi()).g();
								
		/* following two heat exchange excitation, deexcitation */
		(*tr).Coll().cool() = 0.;
		(*tr).Coll().heat() = 0.;

		/* intensity of line */
		(*tr).Emis().xIntensity() = 0.;
		(*tr).Emis().phots() = 0.;
		(*tr).Emis().ots() = 0.;
	}

	photodissoc_BigH2_H2s = 0.;
	photodissoc_BigH2_H2g = 0.;
	HeatDiss = 0.;
	HeatDexc = 0.;
	HeatDexc_deriv = 0.;
	Solomon_dissoc_rate_g = 0.;
	Solomon_dissoc_rate_s = 0.;
	return;
}

/*mole_H2_LTE sets Boltzmann factors and LTE unit population of large H2 molecular */
void diatomics::mole_H2_LTE( void )
{
	DEBUG_ENTRY( "mole_H2_LTE()" );

	/* do we need to update the Boltzmann factors and unit LTE populations? */
	if( ! fp_equal( phycon.te, TeUsedBoltz ) )
	{
		double part_fun = 0.;
		TeUsedBoltz = phycon.te;
		/* loop over all levels setting H2_Boltzmann and deriving partition function */
		for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			long iVib = (*st).v();
			long iRot = (*st).J();
			H2_Boltzmann[iElec][iVib][iRot] = 
				/* energy is relative to lowest level in the molecule, v=0, J=0,
				 * so Boltzmann factor is relative to this level */
				sexp( (*st).energy().K() / phycon.te );
			/* sum the partition function - Boltzmann factor times statistical weight */
			part_fun += H2_Boltzmann[iElec][iVib][iRot] * (*st).g();
			ASSERT( part_fun > 0 );
		}
		/* have partition function, set H2_populations_LTE (populations for unit H2 density) */
		for( qList::const_iterator st = states.begin(); st != states.end(); ++st )
		{
			long iElec = (*st).n();
			long iVib = (*st).v();
			long iRot = (*st).J();
			/* these are the H2 LTE populations for a unit H2 density -
			 * these populations will sum up to unity */
			H2_populations_LTE[iElec][iVib][iRot] = 
				H2_Boltzmann[iElec][iVib][iRot] * 
				(*st).g() / part_fun;
		}
		if( nTRACE >= n_trace_full ) 
			fprintf(ioQQQ,
			"mole_H2_LTE set H2_Boltzmann factors, T=%.2f, partition function is %.2f\n",
			phycon.te,
			part_fun);
	}

	return;
}

/*H2_Reset called by IterRestart to reset variables that are needed after an iteration */
void diatomics::H2_Reset( void )
{

	DEBUG_ENTRY( "H2_Reset()" );

	if(nTRACE) 
		fprintf(ioQQQ,
		"\n***************H2_Reset called, resetting nCall_this_iteration, zone %.2f iteration %li\n", 
		fnzone,
		iteration );

	/* number of times large molecules evaluated in this iteration,
	 * is false if never evaluated, on next evaluation will start with LTE populations */
	nCall_this_iteration = 0;

	/* these remember the largest and smallest factors needed to
	 * renormalize the H2 chemistry */
	renorm_max = 1.;
	renorm_min = 1.;

	/* counters used by H2_itrzn to find number of calls of h2 per zone */
	nH2_pops = 0;
	nH2_zone = 0;
	/* this is used to establish zone number for evaluation of number of levels in matrix */
	nzone_nlevel_set = 0;

	nzoneAsEval = -1;
	iterationAsEval = -1; 

	TeUsedColl = -1;
	TeUsedBoltz = -1;

	lgEvaluated = false;

	/* zero out array used to save emission line intensities */
	H2_SaveLine.zero();

	/* this was option to print all electronic states in the main printout - but
	 * number of electronic states was not known at initialization so set to -1,
	 * will set properly now */
	if( nElecLevelOutput < 1 )
		nElecLevelOutput = n_elec_states;

	return;
}

/*Yan_H2_CS - cross sections for the photoionization of H2 
 * may be represented analytically in the paper 
 *>>refer	H2	photo cs	Yan, M.; Sadeghpour, H. R.; Dalgarno, A. 1998, ApJ, 496, 1044
 * This is revised version given in ERRATUM 2001, ApJ, 559, 1194 
 * return value is photo cs in cm-2 */
double Yan_H2_CS( double energy_ryd /* photon energy in ryd */)
{

	double energy_keV; /* keV */
	double cross_section; /* barns */
	double x; /* x = E/15.4 */
	double xsqrt , x15 , x2;
	double energy = energy_ryd * EVRYD;

	DEBUG_ENTRY( "Yan_H2_CS()" );

	/* energy relative to threshold */
	x = energy / 15.4;
	energy_keV = energy/1000.0;
	xsqrt = sqrt(x);
	x15 = x * xsqrt;
	x2 = x*x;

	if( energy < 15.4 )
	{
		cross_section = 0.;
	}

	else if(energy >= 15.4 && energy < 18.)
	{
		cross_section = 1e7 * (1 - 197.448/xsqrt + 438.823/x - 260.481/x15 + 17.915/x2);
		/* this equation is obviously negative for x = 1 */
		cross_section = MAX2( 0. , cross_section );
	}

	else if(energy >= 18. && energy <= 30.)
	{
		cross_section = (-145.528 +351.394*xsqrt - 274.294*x + 74.320*x15)/pow(energy_keV,3.5);
	}

	else if(energy > 30. && energy <= 85.)
	{
		cross_section = (65.304 - 91.762*xsqrt + 51.778*x - 9.364*x15)/pow(energy_keV,3.5);
	}

	/* the high-energy tail */
	else 
	{
		cross_section = 45.57*(1 - 2.003/xsqrt - 4.806/x + 50.577/x15 - 171.044/x2 + 231.608/(xsqrt*x2) - 81.885/(x*x2))/pow(energy_keV,3.5);
	}

	return( cross_section * 1e-24 );
}

void diatomics::CalcPhotoionizationRate(void)
{
	DEBUG_ENTRY( "diatomics::CalcPhotoionizationRate()" );

	/* must reevaluate During search phase */
	if( ( nzone_eval!=nzone || iteration_evaluated!=iteration ) || !nzone )
	{
		t_phoHeat photoHeat;
		/* generally not important, do one time per zone */
		photoionize_rate = 
			GammaK( ip_photo_opac_thresh,
				rfield.nupper,
				ip_photo_opac_offset,1.,
				&photoHeat)*
			ionbal.lgPhotoIoniz_On +
			/* Compton recoil ionization - we include this in the H2 photoionization
			 * rate but not the heating rate - factor of two since assume 2H
			 * is same as two H0 at such high energies */
			2.*ionbal.CompRecoilIonRate[ipHYDROGEN][0];
			
		/* photo heating - this has units s-1 - needs H2 density
		 * to become vol heat rate */
		photo_heat_soft = photoHeat.HeatLowEnr * ionbal.lgPhotoIoniz_On;
		photo_heat_hard = photoHeat.HeatHiEnr * ionbal.lgPhotoIoniz_On;
		nzone_eval = nzone;
		iteration_evaluated = iteration;
	}

	return;
}

