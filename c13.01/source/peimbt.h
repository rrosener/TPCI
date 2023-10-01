/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PEIMBT_H_
#define PEIMBT_H_

/** peimbt.h */
struct t_peimbt {
	realnum toiiir, 
	  tbac, 
	  t2hstr, 
	  tohyox, 
	  t2hyox, 
	  tbcthn, 
	  t2o3str, 
	  tsqden; /**<tsqden is the highest hydrogen density to do print peimbert analysis*/
	};
extern t_peimbt peimbt;
#endif /* PEIMBT_H_ */
