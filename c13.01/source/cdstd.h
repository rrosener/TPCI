/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef CDSTD_H_
#define CDSTD_H_

// cdstd.h: define macros to select library API version.

// This must be included before all library #includes.
// Typically this done as part of cddefines.h, only required
// independently when cddefines.h is not first include.

// We *require* only POSIX1990.
// See e.g. Rochkind, Advanced UNIX Programming for more details.

// There appears to be a bug in the POSIX implementation in FreeBSD that was
// imported into Apple Darwin (at least version 11.4.2) where defining _POSIX_SOURCE
// will cause compilation errors as a result of missing type definitions (e.g., u_int)
// when including certain system header files. See e.g. this report:
// http://lists.freebsd.org/pipermail/freebsd-bugs/2011-April/044049.html
#if !defined(__APPLE__) && !defined(__FreeBSD__)
#define _POSIX_SOURCE
#endif

#endif // CDSTD_H_
