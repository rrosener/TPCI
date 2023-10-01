/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef DATE_H_
#define DATE_H_

/* this file is created every morning (US East Coast) by the perl script
 * dateit.pl which lives in the date directory off the
 * current main directory.  It works with the version of date.f
 * that lives in that directory to create a version of date.h
 * containing the current date.  This is then copied to the source directory.  
 *
 * This is the GMT year, 2011 is 111 */
#define YEAR	113
/* month, January is 0, December is 11 */
#define	MONTH	1
/* day is correct */
#define	DAY	20

#endif /* DATE_H_ */
