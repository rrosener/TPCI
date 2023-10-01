/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "cddefines.h"

int main ()
{
	ioQQQ = stdout;
	return UnitTest::RunAllTests();
}
