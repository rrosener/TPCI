/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CARB_H_
#define CARB_H_

/** pmp2s.h */
struct t_carb {
	double c8727,
  	  c9850;

	/** correction for depopulation of excited state of CI */
	realnum r9850;

	};
extern t_carb carb;


#endif /* CARB_H_ */
