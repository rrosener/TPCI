/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

/* CHANGES: (M. Salz 17.05.2013)
 *  - introduce the function
      interpol_tabulated() here, because radius values
      are necessary for the interpolation */

#ifndef RADIUS_H_
#define RADIUS_H_

/* radius.h */

/**radius_next use adaptive logic to find next zone thickness 
 * return 0 if ok, 1 for abort */
int radius_next(void);

/**radius_first derive thickness of first zone */
void radius_first(void);

/**radius_increment do work associated with geometry increments of this zone */
void radius_increment(void);

/**interpol_tabulated interpolates on table of points for density/temparature or velocity
\param r0
\param depth
\param lgDepth
\param lgLinear
\param *tbrad
\param *tbval
\param numvals
*/
double interpol_tabulated(double r0, double depth,
	bool lgDepth, bool lgLinear,
	realnum* tbrad, realnum* tbval,
	long int numvals );

struct t_radius {
	double 
		/** the inner radius in cm */
		rinner, 

		/** the outer radius of the current zone */
		Radius, 

		/** the radius, to center of last or current zone */
		Radius_mid_zone,

		/** the thickness of the current zone */
		drad, 

		/** the distance between middle of previous zone and middle of this zone */
		drad_mid_zone,

		/** the depth, the distance from the outer edge of
		 * current zone to the illuminated face */
		depth, 

		/** depth from illuminated face to center of last or current zone */
		depth_mid_zone,

		/** an estimate of the depth to the shielded face */
		Depth2Go,

		/** ratio of square of outer edge of current zone to 
		 * radius of illuminated face of cloud - note continuum
		* is relative to outer edge after ZoneDone is called too*/
		r1r0sq,

		/** total physical thickness of modeled region, (NOT OUTER RADIUS)
		 * this can set set as a stopping criteria, but if
		 * not set is 1e30 before first iteration.  At end of each iteration, thickness
		 * is set to total depth from illuminated face to outer edge */
		*StopThickness/**[ITR DIM]*/;

	/** stopping radius for the model, set with STOP RADIUS command */
	double *StopRadius/**[ITR DIM]*/;

	/** next dr, as set in nextdr */
	double drNext;

	/** the distance to the object from Earth, 
	 * set with the distance command */
	double distance;

	/** sign of dr for going in or out, 1 (usually) or -1 */
	double dRadSign;

	/** drad_x_fillfac is drad * filling factor */
	double drad_x_fillfac;

	/** integrated dReff, integral of depth times filling factor */
	double depth_x_fillfac;

	/** darea_x_fillfac is 2pi * radius * drad * filling factor */
	double darea_x_fillfac;

	/** dVeff is effec vol relative to inner radius,
	 * this version is not affected by the APERTURE command */
	double dVeffVol;

	/** dVeff is effec vol relative to inner radius
	 * this version is affected by the APERTURE SLIT | BEAM commands
	 * it should ONLY be used for quantities observed through the aperture
	 * if the APERTURE command is not used, dVeffAper and dVeffVol are identical */
	double dVeffAper;

	/** dRNeff is next dr effective radius */
	double dRNeff;

	/** dVolOutwrd, dVolReflec, outward and reflected effective vols
	 * used to get outward and reflected beams,
	 * these include only the vol of the current shell times the covering
	 * factor, and a number between 0 and 1 that is the fraction of the beam that goes
	 * out or is reflected.  this is determined by the rt covering factor */
	double dVolOutwrd; 
	double dVolReflec;

	/** Beam vars are related to lines where inward and outward fracs known */
	/** BeamInIn inward part of inwardly directed beam, 0 if sphere */
	double BeamInIn;

	/** BeamInOut outward part of inwardly directed beam, 0 if not sphere */
	double BeamInOut;

	/** BeamOutOut outward part of outwardly directed beam */
	double BeamOutOut;

	/** flag saying that zone thickness became too small, likely because
	 * of an uncontrolled oscillation */
	bool lgdR2Small;

	/** this says whether radius has been set - if true then can do luminosities,
	 * if false then only intensities */
	bool lgRadiusKnown;

	/** lgCylnOn set true when cylinder command given
	 * cylind is half height in centimeters */
	double CylindHigh;
	bool lgCylnOn;

	/** default inner radius when none set, log r =25 in scalar */
	double rdfalt;

	/** variables that deal with the globule command,
	 * glbden, the density */
	realnum glbden, 
	  /** the radius for the globule command */
	  glbrad, 
	  /** the globule power */
	  glbpow, 
	  glbdst;

	/** flag to turn off dr checking in dextdr when globule command entered	*/
	bool lgDrMnOn;

	/** lgPredLumin flag set true if intensities entered into \f$4\pi\f$ st */
	bool lgPredLumin;

	/** log of4 pi r_inner^2, 0 if intensities are printed, 
	 * but is log of 4pi r_o^2 if any luminosity commands are entered */
	realnum pirsq;

	/** additive factor to convert stored line intensities within code 
	 * into a final desired unit, luminosity, flux at Earth, or surface brightness */
	double Conv2PrtInten; 

	/** these are 1e-30 and 1e30 by default, and are set with "set dr" cmnds
	 * used as one of a pair of limits to how big or small zones get
	 * set dr command forces constant dr by setting both to same number */
	double sdrmin;
	double sdrmax;
	double lgFixed;
	// minimum dr relative to depth into cloud.
	double sdrmin_rel_depth;

	/** false, then sdrmin/sdrmax are limits to step size in cm,
         *  true, then are relative fraction of the current radius */
	bool lgSdrminRel;
	bool lgSdrmaxRel;

	/**lgSMinON is flag saying that set drmin has been enteed*/
	bool lgSMinON;

	/** this flag controlled in radius_first and says whether this option
	 * caused the first zone to have larger than optimal thickness */
	bool lgDR2Big;

	/** fraction of initial thickness, set in firstdr
	* do not let dr get smaller than this */
	/** NB - drMinimum not used in code - delete? */
	realnum drMinimum;

	/** min and max dr found in previous iteration */
	double dr_min_last_iter;
	double dr_max_last_iter;

	/** set true is calculations stops because zone thickness gets too small */
	bool lgDrMinUsed;

	/** fractional change used in nextdr */
	realnum drChange;

	/** the Stromgren thickness */
	realnum thickness_stromgren;

	};

extern t_radius radius;


#endif /* RADIUS_H_ */
