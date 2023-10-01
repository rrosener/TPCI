/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CONV_H_
#define CONV_H_

/**\file conv.h 
 * this is the series of routines that converge the pressure, temperature, 
 * electron density, and ionization, for a zone.  Ideally, only the top
 * routine, ConvPresTempEdenIoniz, should be public 
 */

/**ConvIterCheck check whether model has converged or whether more iterations
 * are needed - implements the iter to converg comnd */
void ConvIterCheck( void );

/**ConvInitSolution drive search for initial solution at illuminated face,
 * called by cloudy, returns 1 if outside conditions appropriate for cloudy  */
bool ConvInitSolution();

/**ConvPresTempEdenIoniz solve for current pressure, calls PressureChange, ConvTempEdenIoniz,
 * called by cloudy 
 * returns 0 if ok, 1 if disaster */

int ConvPresTempEdenIoniz(void);

/**ConvTempEdenIoniz determine  temperature, called by ConPresTempEdenIoniz,
 * calls ConvEdenIoniz to get electron density and ionization 
 * returns 0 if ok, 1 if disaster */
int ConvTempEdenIoniz(void);

/**ConvEdenIoniz called by ConvTempEdenIoniz and ConvInitIonize, 
 * it calls ConvIoniz and converges the electron density 
 * \return returns 1 if error condition, 0 if ok */
int ConvEdenIoniz(void);
 
/**ConvIoniz called by ConvEdenIonz, it calls ConvBase until converged 
 * \return returns 1 if error condition, 0 if ok */
int ConvIoniz(void);

/**ConvFail handle convergence failure 
   \param chMode[] chMode is one of "pres", "eden", "ioni", "pops", "grai", "temp"
   \param chDetail string giving details about the convergence failure
*/
void ConvFail(
	/* chMode is one of "pres", "eden", "ioni", "pops", "grai", "temp" */
	const char chMode[],
	/* chDetail - string giving details about the convergence failure */
	const char chDetail[] );

/**ConvBase main routine to drive ionization solution for all species, find total opacity
 * called by ConvIoniz/
 *lgConverg check whether ionization of element nelem has converged 
 *\param loopi this tells how many times ConvBase has been called by ConvIoniz while
 * trying to converge electron density == 0 on first call - allows
 * logic in ConvBase to check for ots oscillations
 */
int ConvBase(long loopi);

/**eden_sum sum free electron density over all species, sets variable erredn.EdenTrue
 * called by ConvEdenIoniz which actually controls the electron density updates 
 * returns 0 if all is ok, 1 if need to abort calc */
int eden_sum(void);

/**EdenChange - changes electron density and dependent variables */
void EdenChange( double EdenNew );

// Variables monitoring progress of convergence
enum counter_type
{
	MOLE_SOLVE, 
	MOLE_SOLVE_STEPS, 
	NEWTON, 
	NEWTON_LOOP,
	CONV_BASE_CALLS,
	CONV_BASE_LOOPS,
	CONV_BASE_ACCELS,
	ION_SOLVES,
	ISO_LOOPS,
	EDEN_CHANGES,
	TEMP_CHANGES,
	PRES_CHANGES,
	NTYPES
};

/**
 * the variables that deal with the convergence of the model 
 */
struct t_conv {

	/** this says why the electron density did not converge,  */
	char chConvEden[INPUT_LINE_LENGTH];

private:
	/** this says why the ionization did not converge, reasons can be a large
	 * change in the level of ionization, or in the heating */
	char m_chConvIoniz[INPUT_LINE_LENGTH];

public:
	void resetConvIoniz()
	{
		m_lgConvIoniz = true;
		strncpy(m_chConvIoniz, "NONE!!!!!", INPUT_LINE_LENGTH-1);
		m_BadConvIoniz[0] = 0.0;
		m_BadConvIoniz[1] = 0.0;
	}
	void setConvIonizFail(const char* reason, double oldval, double newval)
	{
		m_lgConvIoniz = false;
		strncpy(m_chConvIoniz, reason, INPUT_LINE_LENGTH-1);
		m_chConvIoniz[INPUT_LINE_LENGTH-1] = '\0';
		m_BadConvIoniz[0] = oldval;
		m_BadConvIoniz[1] = newval;
	}
	bool lgConvIoniz() const
	{
		return m_lgConvIoniz;
	}
	const char *chConvIoniz() const
	{
		return m_chConvIoniz;
	}
	double convIonizOldVal() const
	{
		return m_BadConvIoniz[0];
	}
	double convIonizNewVal() const
	{
		return m_BadConvIoniz[1];
	}


	/** this gives the reason the model was declared not converged
	* when 'iter to convergence' command is given */
	char chNotConverged[INPUT_LINE_LENGTH];

private:
	/** this flag is used in ConvPresTempEdenIoniz to check that ionization has converged */
	bool m_lgConvIoniz;

public:
	/** this flag is used in ConvPresTempEdenIoniz to check that populations have converged */
	bool lgConvPops;

private:
	/** when the lgConvIoniz flag is set false, the old and new numbers,
	 * the reason for the lack of convergence, should be set to following */
	double m_BadConvIoniz[2];

public:
	/** this will count the number of ionizations in one call from ConvPresTempEdenIoniz*/
	long int nPres2Ioniz;

	/** first sweep through solvers in this zone. also true is search phase */
	bool lgFirstSweepThisZone;
	/** last sweep through solvers in this zone, so update fine opacities */
	bool lgLastSweepThisZone;

	/** a limit to the above, in case where one zone takes forever to not converge,
	 * usually very large, set with SET PRESIONIZ command */
	long int limPres2Ioniz;

	/** counts the number of calls to conv base after iteration starts
	 * can be used to determine very first pass through an iteration
	 * reset to zero at start of each iteration in a simulation */
	long int nTotalIoniz;

	/** the same counter but set to zero after the initial solution is
	 * converged, so this is a measure of the number of calls within
	 * true zones */
	long int nTotalIoniz_start;

	/** conv.lgSearch is true during initial temp-ion search phase,
	 * false after first zone established */
	bool lgSearch;

	/** remember the average electron density error */
	realnum AverEdenError;

	/** remember the biggest and average heating-cooling error */
	realnum BigHeatCoolError;
	realnum AverHeatCoolError;

	/** remember the biggest and average pressure error */
	realnum BigPressError;
	realnum AverPressError;

	/** flag set in ConvBase, saying whether ionization stage is trimmed down */
	bool lgIonStageTrimed;

	/** this is true if ots rates are oscillating, and this is why ionization
	 * is not converged */
	bool lgOscilOTS;

	/** true if temperature is converged, false if not */
	bool lgConvTemp;

	/** true if pressure is converged, false if not */
	bool lgConvPres;

	/** true when the electron density has converged */
	bool lgConvEden;

	/** total number of all falures, used to trigger abort */
	long int nTotalFailures;

	/**nTeFail number of temperature failures*/
	long int nTeFail;

	/**failmx is largest relative error in heating cooling match*/
	realnum failmx;

	/**nPreFail is number of pressure failures*/
	long int nPreFail;

	/**nNeFail is number of electron density failures*/
	long int nNeFail;

	/** remember the biggest electron density error as test of convergence quality */
	realnum BigEdenError;

	/**nIonFail is number of ionization failures*/
	long int nIonFail;

	/**nIonFail is number of level population failures */
	long int nPopFail;

	/** number of grain ionization balance failures */
	long int nGrainFail;

	/** number of chemistry solution failures */
	long int nChemFail;

	/** LimFail is limit to number of te failures, set with "failures" cmnd */
	long int LimFail;

	/** lgMap is option to map failures */
	bool lgMap;

	/** zones where converge problems occurred */
	long int ifailz[12];

	/** which electron density solver to use? 
	 * set with set eden solver command, simple and new */
	char chSolverEden[20];

	/** which temperature density solver to use? 
	 * set with set eden solver command, default and brent */
	char chSolverTemp[20];

	/** flag saying that calculation stopped for bad reason
	 * mostly set in lgEndfun */
	bool lgBadStop;

	/** says "interate to convergence" command given */
	bool lgAutoIt;

	/** says update inter-couplings after every ion solution */
	bool lgUpdateCouplings;

	/** a convergence criteria */
	realnum autocv;

	/** this is relative error in the electron density we want
	 * set in zero to 0.01
	 * reset with set eden error command */
	double EdenErrorAllowed;

	/** this is relative error in the pressure,
	 * initialized to 0.02 in 
	 * reset with set pressure error command */
	realnum PressureErrorAllowed;

	double MaxFractionalDensityStepPerIteration;

	/** allowed error in heating - cooling balance, set with TOLERANCE command
	 * default set in zerologic */
	realnum HeatCoolRelErrorAllowed;

	realnum IonizErrorAllowed;

	/** allowed error in total gas-phase density of each element, including molecules
	 * change with "SET DENSITY TOLERANCE" command
	 * default set in zero */
	realnum GasPhaseAbundErrorAllowed;

	/** numerical estimate of d(cooling-heating)/dT */
	double dCmHdT;

	/** 1-sigma uncertainty in numerical estimate of d(cooling-heating)/dT */
	double sigma_dCmHdT;

	/** these are used to retain the density pressure history in current zone,
	 * can be output with save pressure history */
	vector<double> hist_pres_density, hist_pres_current, hist_pres_error;
	long int hist_pres_nzone;

	/** these are used to retain the temp/heat/cooling history in current zone,
	 * can be output with save temperature history */
	vector<double> hist_temp_temp, hist_temp_heat, hist_temp_cool;
	long int hist_temp_nzone;

private:
	// Variables monitoring progress of convergence
	long m_counters[NTYPES];
	long m_counters_zone[NTYPES];
public:
	void incrementCounter( const counter_type type )
	{
		++m_counters[type];
		++m_counters_zone[type];
	}
	void resetCounters()
	{
		for( long i=0; i<NTYPES; ++i )
			m_counters[i] = 0;
	}
	void resetCountersZone()
	{
		for( long i=0; i<NTYPES; ++i )
			m_counters_zone[i] = 0;
	}
	long getCounter( const long type )
	{
		return m_counters[type];
	}
	long getCounterZone( const long type )
	{
		return m_counters_zone[type];
	}
};


extern t_conv conv;

#endif /* CONV_H_ */
