/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

/* CHANGES: (M. Salz 21.05.2013)
 *  - define routine cdEDEN_depth( double )
 *    for obtaining electron density structure 
 *    after the simulation
 *  - cdDenPart_depth
 *    cdDenMass_depth
 *    cdWindVel_depth
 *    cdCooling_depth
 *    cdHeating_depth 
 *    cdRadAcce_depth*/

#ifndef CDDRIVE_H_
#define CDDRIVE_H_

/**\file cddrive.h
 \verbatim
  this file contains definitions for the suite of routines that
  allow the code to be driven as a subroutine.
  These routines set up model parameters, 
  control the execution of Cloudy, and obtain results once complete
  these are the only "public" routines, and only these should
  be accessed when controlling Cloudy
 
  DRIVING CLOUDY FROM A FORTRAN PROGRAM:
  This should not be too hard - the recommended approach is to use
  the cfortran.h file described at http://www-zeus.desy.de/~burow/cfortran/ 

  A note on return error conditions:

  Many functions return an int to indicate success or failure.
  I try to follow standard Unix/C convention.
  A return value of zero usually indicates that the routine was successful,
  and a non-zero value (not always 1) indicates failure.  This is conventional
  in both C and Unix.  So the way to call Cloudy and test for success is
 
  if( cdDdrive() )
  {
 	   printf(" Cloudy failed.\n");
  }
 
  Although I try to follow this, there ARE exceptions.
 \endverbatim
 */

/**
 * cdInit 
 * This routine must be called before any of the others - 
 * it reinitializes many variables, and must be called before any
 * of the other routines.  In a large grid of calculations it must be repeatedly called
 * before the start of the new calculation and after all results have
 * been obtained from the previous model */
void cdInit();

/**
 * cdTalk 
 * tells the code whether or not to produce any of its normal output, 
 * If the argument is true (or if it is not called at all) it produces 
 * output, produces no output if it is false */
void cdTalk(bool);

/**
 * cdOutput
 * This tells the code where to send output.  The arguments are as
 * for the stdio.h fopen call, but the resulting file pointer is checked 
 * for validity.  All further log output will go to this file.  
 * If filename = "", output is switched to stdout (and mode is ignored).
 * If this routine is not called then all output will go to 
 * stdout, the standard c output */
void cdOutput( const char* filename = "", const char *mode = "w" );

/**
 * cdInput
 * This tells the code where to get input.  The arguments are as
 * for the stdio.h fopen call, but the resulting file pointer is checked 
 * for validity.  All further input will come from this file.  
 * If filename = "", input is switched to stdin (and mode is ignored).
 * If this routine is not called then all input will come from 
 * stdin, the standard c input */
void cdInput( const char* filename = "", const char *mode = "r" );

/** 
 * returns depth structure of previous model 
 * \param cdDepth[]
*/
void cdDepth_depth( double cdDepth[] );

/**
 * cdnZone 
 * returns number of zones */
long int cdnZone();

/**
 * cdB21cm
 * returns B as measured by 21 cm 
 * assumes tangled field weighted by n(H0)/T */
double cdB21cm();

/**
 * cdRead 
 * This sends commands to the code.  The normal set of commands
 * described in Part I of Hazy must be entered into a null-terminated
 * string.  These strings are then fed to Cloudy with this command.  The function
 * returns the number of commands that can still be entered before the command
 * stack is full.  The code will stop if you try to continue giving it commands
 * after the command has returned zero. This return value is the opposite of the
 * standard - a non-zero return is normal 
*/
int cdRead( const char* );

/** 
 * cdPrtWL print line wavelengths in Angstroms in the standard format - 
 * a wrapper for prt_wl which must be kept parallel with sprt_wl
 * both of those live in pdt.c 
 *\param [out] *ioOUT output file handle
 *\param [in] the emission line wavelength
*/
void cdPrtWL( FILE *io , realnum wl );

/** debugLine provides a debugging hook into the main line array 
* loops over whole array and finds every line that matches length,
* the wavelength, the argument to the function
* put breakpoint inside if test 
* return value is number of matches, also prints all matches
*\param [in] the emission line wavelength
*\param [out] the number of matches
*/
long debugLine( realnum wavelength );

/**
 * cdNoExec 
 * This provides option to have the code prepare the initial conditions for a model,
 * but not actually try to compute the model.  I use this when setting up a large
 * grid so that I can quickly run through the full grid as a check that the commands
 * are entered properly and the parameters I am going to vary do so properly.
 * The command is then commented out when the grid is properly set up. */
void cdNoExec();

/**
 * cdDrive 
 * This command actually computes a model.
 * It returns 0 if the calculation was successful, and 1 if an error
 * condition was encountered */
int cdDrive();


/* The next two routines confirm that the previous calculation was ok 
 * or produce a list of error conditions */

/**
 * cdErrors
 * After the calculation is completed, a summary of all error messages can be
 * be generated by calling this routine.  The argument is the output file 
 *\param [out] *ioOUT output file
 */
void cdErrors(FILE* );

 /**
 * cdNwcns
 * This command returns the number of warnings, cautions, notes, surprises, 
 * assorted types of failures found the last computed model
   
  \param *lgAbort abort status, if non-zero then big problems happened
  \param *NumberWarnings the number of warnings
  \param *NumberCautions the number of cautions  
  \param *NumberNotes the number of notes   
  \param *NumberSurprises the number of surprises
  \param *NumberTempFailures the number of temperature convergence failures
  \param *NumberPresFailures the number of pressure convergence failures
  \param *NumberIonFailures the number of ionization convergence failures 
  \param *NumberNeFailures the number of electron density convergence failures
 */ 
void cdNwcns(
  bool *lgAbort ,
  long int *NumberWarnings, 
  long int *NumberCautions, 
  long int *NumberNotes, 
  long int *NumberSurprises, 
  long int *NumberTempFailures, 
  long int *NumberPresFailures,
  long int *NumberIonFailures, 
  long int *NumberNeFailures );

/** This prints the reason why the model stopped, and the model geometry, on 
 * the io file pointed to by the file handle */
void cdReasonGeo(FILE*);

/**
 * These routines produce lists of warnings, cautions, notes, surprises
 * These are the routines that are called by cdErrors.  Normally
 * cdErrors and not these routines would be called by the user. */
void cdWarnings(FILE*);
 /**
  * produces list of cautions*/
void cdCautions(FILE*);
 /** produces list of surprises*/
void cdSurprises(FILE*);
 /** produces list of Notes */
void cdNotes(FILE*);

/***********************************************************
 *
 * The next routines examine the predictions of the previous model
 *
 ***********************************************************/

/** 
 * cdLine
 * This routine finds the predicted intensity of any line in the standard output.
 *
 * 
 * \param *chLabel 1st parameter is the 4-char label + null terminated label, as it appears in the output.
 * \param wavelength 2nd parameter is the float wavelength in Angstroms, not how it appears in printout.
 *   The first four digits must agree with the number in the printout, but the units must be Angstroms.
 * 3rd parameter is the predicted intensity relative to the normalization line.
 * 4th par is the log of the predicted luminosity or intensity of the line (ergs).
 * \param *relint 5th is pointer to relative intensity, a double that is returned
 * \param *absint 6th is pointer to log of absolute intensity
 * \param lEmergent - emergent or intrinsic intensity
 *
 * \return return value:
 * The routine returns an index (>0) of the array element within stack if it finds the line, 
 * It returns the negative of the total number of lines if it could not find the line.
 * (this is a debugging aid)
 * note that this returns a long int since there are LOTS of lines
 * this also IS NOT the standard C convention for success or failure */
 
long int cdLine(
	const char *chLabel, 
	realnum wavelength, 
	double *relint, 
	double *absint);

long int cdLine(
	const char *chLabel, 
	realnum wavelength, 
	double *relint, 
	double *absint,
	// 0 is intrinsic,
	// 1 emergent
	// 2 is intrinsic cumulative,
	// 3 emergent cumulative
	int LineType );


 /**cdLine_ip get the predicted line intensity, using index for line in stack 
 \param ipLine
 \param *relint linear intensity relative to normalization line
 \param *absint log of luminosity or intensity of line
 \param lgEmergent - intrinsic or emergent intensity
 */ 
void cdLine_ip(long int ipLine, 
	  double *relint, 
	  double *absint ,
		// 0 is intrinsic,
		// 1 emergent
		// 2 is intrinsic cumulative,
		// 3 emergent cumulative
		int LineType );
void cdLine_ip(long int ipLine, 
			   double *relint, 
			   double *absint );

/**
 \verbatim
 * cdColm 
 * This obtains the column density of a species in the previously computed model.
 * The first parameter is a 4 character + NULL terminated string which gives
 * the first 4 char of the element name as spelled by Cloudy, either upper or lower case.
 * The second parameter is the stage of ionization, 1 for atom, 2 for first ion, etc; 0 is special.
 *
 * examples: 
 * column density of atomic carbon
 * cdColm( "carb" , 1 , &col1 );
 *
 * doubly ionized helium
 * cdColm( "heli" , 3 , &col3 );
 *
 * molecular hydrogen
 * cdColm("H2  " , 0 , &col2 );
 *
 * If the ion stage is zero then the routine will check the first label
 * for the values "H2  ", "OH  ", "CO  " and "CII* *", 
 * and will return the H2, OH, CO or CII* column density in this case
 *
 * The column density [cm-2] is returned as the third parameter in all cases
 *
 * The function returns 0 if it found the species, 1 if it failed 
 \endverbatim
*/

int cdColm(const char*, long, double* );

/**cdH2_colden return column density in H2, returns -1 if cannot find state,
 * header is in cdDrive, source in h2.c */

double cdH2_colden( long iVib , long iRot );

/**
 * cdEmis
 * This routine finds the local emissivity for any line.
 * The first argument is the 4 character (null terminated) label as it appears in the 
 * standard output. 
 * The second argument is float wavelength as it appears in the standard output.
 * The emissivity (erg /cm^3 /s) is returned as the last parameter.
 * cdEms returns the index of the line in the stack if it was found, 
 * the negative of the total number of lines in the stack if it could not find the line 
 \param *chLabel 4 char null terminated string label
 \param wavelength line wavelength
 \param *emiss the vol emissivity of this line in last computed zone
 \param lgEmergent intrinsic or emergent emissivities
 \return return value will be index of line within stack
*/
// this version returns intrinsic
long int cdEmis(
	char *chLabel,
	realnum wavelength, 
	double *emiss );

long int cdEmis(
	char *chLabel,
	realnum wavelength, 
	double *emiss ,
	bool lgEmergent );


/** cdEms_ip obtain the local emissivity for a line with known index
 \param ipLine index of the line in the stack
 \param *emiss the vol emissivity of this line in last computed zone 
 \param lgEmergent intrinsic or emergent emissivities
*/
void cdEmis_ip(
	long int ipLine, 
	double *emiss ,
	bool lgEmergent);

/**
 * cdCooling_last
 * The returns the total cooling (erg cm^-3 s^-1) for the last computed zone */
double cdCooling_last();

/**
 * cdCooling_depth
 * The returns the total cooling (erg cm^-3 s^-1) structure 
 * \param Cool_struc[]
*/
void cdCooling_depth( double Cool_struc[] );

/**
 * cdHeating_last
 * returns the total heating (erg cm^-3 s^-1) for the last computed zone */
double cdHeating_last();

/**
 * cdHeating_depth
 * returns the total heating (erg cm^-3 s^-1) structue 
 * \param Heat_struc[]
*/
void cdHeating_depth( double Heat_struc[] );

/**cdDenPart_depth return particle density struc. of previous model */
void cdDenPart_depth( double DenPart[] );

/**cdDenMass_depth return mass density struc. of previous model */
void cdDenMass_depth( double DenMass[] );

/**cdWindVel_depth return velocity struc. of previous model */
void cdWindVel_depth( double WindVel[] );

/**cdRadAcce_depth 
 * return total radiative acceleration struc. 
 * of previous model (outward direction) 
 */
void cdRadAcce_depth( double RadAccel[] );

/**cdEDEN_last return electron density of last zone */
double cdEDEN_last();

/** 
 * returns electron density structure of previous model 
 * \param cdEDEN[]
*/
void cdEDEN_depth( double cdEDEN[] );

 /** cdPressure_last
 * This returns the pressure and its constituents for the last computed zone. 
 \param  *TotalPressure total pressure, all forms
 \param  *GasPressure gas pressure 
 \param  *RadiationPressure radiation pressure
 */ 
void cdPressure_last(
	double *TotalPressure,
	double *GasPressure,
	double *RadiationPressure);

/**
 * cdPressure_depth
 * This returns the pressure and its constituents for the last iteration. 
 * space was allocated in the calling routine for the vectors - 
 * before calling this, cdnZone should have been called to get the number of
 * zones, then space allocated for the arrays
 \param  TotalPressure[] total pressure, all forms
 \param  GasPressure[] gas pressure
 \param  RadiationPressure[] radiation pressure
 */ 
void cdPressure_depth(
	double TotalPressure[],
	double GasPressure[],
	double RadiationPressure[]);

/**
 * cdTemp_last
 * returns the temperature of the last zone on last iteration */
double cdTemp_last();

/**
 \verbatim
 * cdIonFrac
 * This returns the ionization fraction for any element included in the calculation. 
 * The first parameter is 4 char null terminated string giving the first 4 letters of
 * element name as spelled by Cloudy.  
 * The second parameter is an integer giving the ionization stage, 
 * 1 for atom, 2 for first ion, etc.
 * The third parameter returns the predicted ionization fraction of that ion stage.
 * The last parameter is an 8 character + null string that says either "volume" or "radius",
 * to specify whether the average should be weighted by volume or radius.
 * The return value is 0 if the routine could find the species and
 * non-zero if it failed to find the element 
\endverbatim 
\param *chLabel four char string, null terminated, giving the element name 
\param IonStage IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number
\param *fracin  will be fractional ionization 
\param *chWeight how to weight the average, must be "VOLUME" or "RADIUS" 
\param lgDensity if true then weighting also has electron density, if false then only volume or radius 
 */ 
int cdIonFrac(
	const char *chLabel, 
	long int IonStage, 
	double *fracin, 
	const char *chWeight ,
	bool lgDensity );

/**
 * cdVersion
 * The argument is a string with at least 8 characters that will receive a null terminated
 * string with the version number of the code. */
void cdVersion(char chString[] );

/**
 * cdDate
 * The argument is a string with at least 8 char that will receive a null terminated
 * string with the date of the current version of the code. */
void cdDate(char chString[] );

/* The following pairs of routines can keep track of the execution time for one model -
 * cdSetExecTime called first (in cdInit, not by the user) to initialize timer.
 * When cdExecTime is called it will return the elapsed time in seconds
 * since cdInit called cdSetExecTime*/

/** normally called by cdInit, this routine sets initial variables for times */
void cdSetExecTime();

/** cdExecTime returns the elapsed time cpu time (sec) that has elapsed 
 * since cdInit called cdSetExecTime.*/
double cdExecTime();

/**
 \verbatim
 * cdGetLineList will read in a list of emission line labels and wavelengths
 * from a file.  I use it for generating LOC grids.  
 * Two files (cdGetLineList and cdGetLineList) are included in the main data 
 * distribution and have list of strong emission lines for high and low density gas.
 * other files can be created by the user.
 *
 * The first argument is the name of the file to read.
 * It it is void ("") then the routine will open LineList_BLR.dat 
 *
 * The next two arguments are references to vectors holding the
 * list of labels and wavelengths.  The routine will allocate the
 * needed space, but the vectors are defined in the calling routine.  
 * in the calling routine the two variable should be declared like this:
 * vector<char*> chLabels;
 * vector<realnum> wavelength;
 * They would appear as follows in the call to the routine:
 * chGetLineList("", chLabels , wavelength );
 *
 * cdGetLineList returns the number of lines it found in the file if it was successful,
 * and -1 if it could not open the file.
 * 
 \endverbatim
*/

long int cdGetLineList(
	const char chFile[],
	vector<char*>& chLabels,
	vector<realnum>& wl);

/** 
 * cdTimescales returns the longest thermal, recombination, and H2 formation 
 * timescales that occurred in the previous model
 \param  *TTherm the thermal cooling timescale
 \param  *THRecom the hydrogen recombination timescale 
 \param  *TH2 the H2 formation timescale
 */ 
void cdTimescales(
	double *TTherm , 
	double *THRecom , 
	double *TH2 );

/* ******************************************************************
 *
 * next part deals with FeII bands.  There are two types, the tabulated
 * band that are defined in FeII_bands.ini, and the psuedo-continuum bins
 * that are generatedby the code in FeIIContCreate.
 * nFeIIConBins is number of continuum bins in FeII_Cont 
 * nFeIIBands is number of bands in FeII_bands.ini, and are saved in FeII_Bands
 * the bands are created by hand and the entries in FeII_bands.ini are
 * meant to be created by a person  */

/* the declarations for the next four are in FeIILevelPops.c */
/** this is the number of bands read in from FeII_bands.ini */
extern long int nFeIIBands;

/** number of bands in continuum array */
extern long int nFeIIConBins;

/* band wavelength, lower and upper bounds, in vacuum Angstroms */
/** FeII.bands[n][3], where n is the number of bands in fe2bands.dat
 * these bands are defined in fe2bands.dat and read in at startup
 * of calculation */
extern realnum **FeII_Bands; 

/* continuum wavelengths, lower and upper bounds, in vacuum Angstroms
 * third is integrated intensity */
/** FeII_Cont[n][3], where n is the number of cells needed
 * these bands are defined in cdGetFeIIBands */
extern realnum **FeII_Cont; 


/** 
 * this routine returns the spectrum needed for Keith Arnaud's XSPEC
 * X-Ray analysis code.  It should be called after cdDrive has successfully
 * computed a model.  The calling routine must ensure that the  vectors
 * have enough space to store the resulting spectrum, 
 * given the bounds and energy resolution 
 \param Option 
 	\verbatim 
	 option - the type of spectrum to be returned
	 1			the incident continuum 4\pi nuJ_nu, , erg cm-2 s-1
	
	 2			the attenuated incident continuum, same units
	 3			the reflected continuum, same units
	
	 4			diffuse continuous emission outward direction
	 5			diffuse continuous emission, reflected
	
	 6			collisional+recombination lines, outward
	 7			collisional+recombination lines, reflected
	
				all lines and continuum emitted by the cloud assume full coverage of 
				continuum source
	\endverbatim 
 \param nEnergy the number of cells + 1
 \param ReturnedSpectrum[] the returned spectrum, same size is two energy spectra (see option), returns nEnergy -1 pts
*/
void cdSPEC( 
	int Option ,
    long int nEnergy ,
    double ReturnedSpectrum[] );


 /**  
	\param Option 
	\verbatim the type of spectrum to be returned
	 1			the incident continuum 4\pi nuJ_nu, , erg cm-2 s-1
	
	 2			the attenuated incident continuum, same units
	 3			the reflected continuum, same units
	
	 4			diffuse emission outward direction
	 5			diffuse emission, reflected
	
	 6			collisional+recombination lines, outward
	 7			collisional+recombination lines, reflected
	
				all lines and continuum emitted by the cloud assume full coverage of 
				continuum source
	\endverbatim
	 
	\param nEnergy the number of cells + 1
	\param ipLoEnergy
	\param ipHiEnergy
	\param ReturnedSpectrum[] the returned spectrum, same size is two energy spectra (see option), returns nEnergy -1 pts
 */ 
void cdSPEC2( 
	int Option ,
    long int nEnergy ,
	long int ipLoEnergy,
	long int ipHiEnergy,
    realnum ReturnedSpectrum[] );

/** cdTemp \verbatim
 * This routine finds the mean electron temperature for any ionization stage 
 * It returns 0 if it could find the species, 1 if it could not find the species.
 * The first argument is a null terminated 4 char string that gives the element
 * name as spelled by Cloudy.  
 * The second argument is ion stage, 1 for atom, 2 for first ion, etc
 * This third argument will be returned as result,
 * Last parameter is either "VOLUME" or "RADIUS" to give weighting 
 *
 * if the ion stage is zero then the element label will have a special meaning.
 * The string "21CM" is will return the 21 cm temperature.
 * The string "H2  " will return the temperature weighted wrt the H2 density \endverbatim
 \param *chLabel four char string, null terminated, giving the element name
 \param IonStage IonStage is ionization stage, 1 for atom, up to N+1 where N is atomic number
 \param *TeMean will be temperature
 \param *chWeight how to weight the average, must be "VOLUME" or "RADIUS"
 */ 
int cdTemp(
	const char *chLabel, 
	long int IonStage, 
	double *TeMean, 
	const char *chWeight );

/** cdPrintCommands( FILE *)
 * This routine prints all input commands into file whose handle is the argument 
 * \param *ioOUT [out] output file handle
 */
void cdPrintCommands( FILE * );

/** wrapper to close all save files */
void cdClosePunchFiles();

 /**cdH2_Line returns 1 if we found the line, 
  * or false==0 if we did not find the line because ortho-para transition
  * or upper level has lower energy than lower level 
  * NB - this is in mole_h2_io.c  
  \param iElecHi indices for the upper level 
  \param iVibHi indices for the upper level 
  \param iRotHi indices for the upper level 
  \param iElecLo indices for lower level
  \param iVibLo indices for lower level
  \param iRotLo indices for lower level
  \param *relint linear intensity relative to normalization line
  \param *absint log of luminosity or intensity of line 
 */ 
long int cdH2_Line(
	  /* indices for the upper level */
	  long int iElecHi, 
	  long int iVibHi ,
	  long int iRotHi ,
	  /* indices for lower level */
	  long int iElecLo, 
	  long int iVibLo ,
	  long int iRotLo ,
	  /* linear intensity relative to normalization line*/
	  double *relint, 
	  /* log of luminosity or intensity of line */
	  double *absint );

/* none of the following are generally needed */

/** this is the value that will be set true when cdInit is called.  
 * Other routines will check that this is true when they are called, 
 * to verify that cdInit was called first.  Definition is in cdInit.cpp */
extern bool lgcdInitCalled;

#endif /* CDDRIVE_H_ */
