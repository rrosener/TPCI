/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CO_H_
#define CO_H_

/* co.h */

struct t_co {

	/** CODissHeat is CO Photodissociation heating */
	realnum CODissHeat, 
	  /**< largest fraction of total heating */
	  codfrc, 
	  /**< total heating integrated over cloud */
	  codtot;

	/** C12/C13 isotope ratio, sets the ratio of C12O16 to C13O16 and 
	* C13/C12 1909, initialized in zero.c  */
	realnum C12_C13_isotope_ratio;

	double xDens13CO;
};

extern t_co co;

#endif /* CO_H_ */
