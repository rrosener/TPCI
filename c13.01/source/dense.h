/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DENSE_H_
#define DENSE_H_

/* dense.h density related variables */

/* CHANGES: (M. Salz 14.05.2013)
 *  - add the logic variable lgTabLinear to trigger
 *    linear interpolation
 *  - changed variable names for the sake of conformaty
 *    with temperatur and velocity tabulated values
 *  - lgDLWDepth --> lgTabDepth
 *  - frad, fhden --> tabrad, tabval*/

/**dense_fabden called by dlaw command, returns density for any density law 
 \param radius
 \param depth
*/
double dense_fabden(double radius, 
  double depth);

/**dense_tabden interpolate on table of points for density with dlaw table command, by K Volk 
\param r0
\param depth
*/
double dense_tabden(double r0, 
  double depth);

/*dense_parametric_wind called by dlaw command, returns density for any density law */
double dense_parametric_wind(double rad);

const int LIMTABDLAW = 500;

/** variables dealing with pressure across model */
struct t_dense
{
	t_dense() {
		/* list of atomic weights, mass in AMU, used for thermal line widths */
		/* >>refer	all	atomic weight	Coplen, T.B. 2001, J. Phys. Chem REf Data, 30, 701 */
		AtomicWeight[0] = 1.00794f;
		AtomicWeight[1] = 4.0026022f;
		AtomicWeight[2] = 6.9412f;
		AtomicWeight[3] = 9.0121823f;
		AtomicWeight[4] = 10.8117f;
		AtomicWeight[5] = 12.01078f;
		AtomicWeight[6] = 14.00672f;
		AtomicWeight[7] = 15.99943f;
		AtomicWeight[8] = 18.9984032f;
		AtomicWeight[9] = 20.17976f;
		AtomicWeight[10] = 22.989770f;
		AtomicWeight[11] = 24.30506f;
		AtomicWeight[12] = 26.9815382f;
		AtomicWeight[13] = 28.08553f;
		AtomicWeight[14] = 30.9737612f;
		AtomicWeight[15] = 32.0655f;
		AtomicWeight[16] = 35.4532f;
		AtomicWeight[17] = 39.9481f;
		AtomicWeight[18] = 39.09831f;
		AtomicWeight[19] = 40.0784f;
		AtomicWeight[20] = 44.9559108f;
		AtomicWeight[21] = 47.8671f;
		AtomicWeight[22] = 50.94151f;
		AtomicWeight[23] = 51.99616f;
		AtomicWeight[24] = 54.9380499f;
		AtomicWeight[25] = 55.8472f;
		AtomicWeight[26] = 58.9332009f;
		AtomicWeight[27] = 58.69342f;
		AtomicWeight[28] = 63.5463f;
		AtomicWeight[29] = 65.392f;
		
		zero();
	}

	/** dense.gas_phase is the total gas phase abundances,
	 * including anything within molecules,
	 * but not including grains */
	realnum gas_phase[LIMELM];
	void SetGasPhaseDensity( const long nelem, const realnum density );

	/** vector of atomic weights for all elements, set in zerologic */
	realnum AtomicWeight[LIMELM];

private:
	/** dense.xMolecules density of elements locked in molecules, 
	 * this is included in gas_phase */
	realnum m_xMolecules[LIMELM];

public:
	realnum xMolecules(long nelem)
	{
		return m_xMolecules[nelem];
	}
	void updateXMolecules();
	void zero();

	/**xMassDensity grams per cc */
	realnum xMassDensity;
			
	/** WJH: fiducial value that corresponds to hden set in init file, this is 
	 * used for setting the mass-flux in dynamic models */
	realnum xMassDensity0;
		
	/** total number of particles per cubic centimeter */
	realnum pden;
		
	/** mean AMU per particle */
	realnum wmole;
		
	/** total number of nuclei, set in PressureTotal */
	realnum xNucleiTotal;
		
	/** total mass in grams PER 4 pi rinner^2 */
	realnum xMassTotal;

	/** this is scale factor that multiplies the correction factor for neutral
	 * hydrogen collisions, def 1, changed with set command */
	realnum HCorrFac;

	/** indices for lowest stage of ionization of the elements on C scale, 
	 * lowest is 0 for atom, -1 if element turned off, 
	 * the first stage of ionization with positive abundance is [IonLow] where 0 is the atom,
	 * the highest stage of ionization with positive abundance is [IonHigh], 
	 * NB NB so loops should be  
	 * ion=IonLow, ion<=IonHigh */
	long int IonLow[LIMELM+1];
	long int IonHigh[LIMELM+1];

	/** dense.xIonDense[nelem][i] is density of ith ionization stage (cm^-3),
	 * [nelem][0] is atom, [][1]) the first ion
	 * nelem = 0 for H, 1 for he, etc */
	double xIonDense[LIMELM][LIMELM+1];

	/** which ions have chianti enabled? */
	bool lgIonChiantiOn[LIMELM][LIMELM+1];

	/** which ions have stout enabled? */
	bool lgIonStoutOn[LIMELM][LIMELM+1];

	/** Maximum wavenumber in chianti **/
	double maxWN[LIMELM][LIMELM+1];

	/** this is lower limit to abundance of element that will be include
	 * in the calculation, default is zero, set with command
	 * ELEMENT LIMIT OFF XXX command */
	realnum AbundanceLimit;

	/** array of logical variables saying whether an element is enable (true)
	 * or disabled (false).  It is set totally true in zero
	 * and is set false with the "element off" command.  
	 * In SetAbundances if can be reset so that an element that was disabled
	 * on the first model in a core load is not later enabled */
	bool lgElmtOn[LIMELM];

	/** will we solve for ionization (false) or specify it with element ionization cmnd true */
	bool lgSetIoniz[LIMELM];

	/** dense.SetIoniz the ionization fractions that are set when lgSetIoniz set true,
	 * gas phase abundance is this times total aabundance
	 * Ionization fraction for [nelem][ion] */
	realnum SetIoniz[LIMELM][LIMELM+1];

	/** label describing the density law for current calculation
	 * 'DLW2' is dense_tabden interpolated table */
	char chDenseLaw[5];

	/* this says keep initial density constant, 
	 * so pressure from iter to iter not really const */
	bool lgDenseInitConstant;

	// let initial pressure vary with time
	bool lgPressureVaryTime;

	//  required number is timescale for time variation
	double PressureVaryTimeTimescale;
	//  optional number is index for time variation
	double PressureVaryTimeIndex;

	/** parameters set by the dlaw command, used by dense_fabden (maybe) */
	double DensityLaw[10];

	/** options on set atomic data command */
	bool lgAsChoose[LIMELM][LIMELM];
	bool lgCSChoose[LIMELM][LIMELM];

	/**frad is log radius in cm, fhden is log hden*/
	// realnum frad[LIMTABDLAW];
	// realnum fhden[LIMTABDLAW];
	realnum tabrad[LIMTABDLAW];
	realnum tabval[LIMTABDLAW];

	/**number of values in above table */
	long int nvals;

	/**lg is true if depth, false if radius to be used*/
	bool lgTabDepth;

	/**lg is true if linear, false if log interpolation to be used*/
	bool lgTabLinear;
	
	/** electron density, units cm-3  */
	double eden;

	/** max and min eden over this iteration */
	double EdenMax , EdenMin;

	/** lowest allowed density for any ion = if density falls below
	 * this then set to zero in ion_trim */
	double density_low_limit;

	/**zone where bad electron density was detected */
	long int nzEdenBad;

	/**EdenSet electron density set with set eden command*/
	realnum EdenSet;

	/** extra electron density, set with eden command */
	realnum EdenExtra;

	/** option to set electron fraction, n_e/n_H */
	realnum EdenFraction;

	/** square root of electron density, set in tfidle */
	double SqrtEden;

	/** EdenHCorr is eden + hi * 1.7e-4, includes correction for H atom collisions,
	 * units cm-3 */
	double EdenHCorr;

	realnum EdenHCorr_f;

	/** this is the true eden as set in eden_sum, we will try to converge eden to this */
	double EdenTrue;

	/** fraction of electron density due to ions rather than molecules and grains */
	double eden_from_metals;

	/**flags set when bad electron density is detected */
	bool lgEdenBad;

	/** edensqte is eden/sqrte */
	double edensqte; 

	/** cdsqte is eden/sqrte times 8.629e-6 - this multiplies the collision strength to 
	 * produce the deexcitation rate coefficient, s-1 
	 * 8.629e-6 is COLL_CONST	*/
	double cdsqte;

	/** parameters dealing with hydrogen density scaling as power of radius
	 *DensityPower is power */
	realnum DensityPower;
	realnum rscale;
	realnum den0;

	/** set true when density fluctuations are turned on */
	bool lgDenFlucOn;

	/** set false when fluctuations are over col den rather than radius, set with column options
	 * on fluctuations command*/
	bool lgDenFlucRadius;

	/** parameters for the density fluctuations command */
	realnum flong;
	realnum cfirst; 
	realnum csecnd;
	realnum flcPhase;

};
extern t_dense dense;

// test whether nuclei conservation holds
bool lgElemsConserved(void);

#include "iso.h"

// test whether the sum of a set of states agrees with the ion population to the requested tolerance.
void lgStatesConserved ( long nelem, long ionStage, qList states,  long numStates, realnum err_tol, long loop_ion );

void SumDensities( void );

void ScaleAllDensities( const realnum factor );
void ScaleIonDensities( const long nelem, const realnum factor );
bool AbundChange( );

realnum scalingDensity(void);
realnum scalingZoneDensity(long i);

#endif /* DENSE_H_ */
