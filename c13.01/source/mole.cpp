/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "mole.h"

t_mole_global mole_global;
t_mole_local mole;

/*=================================================================*/
/*mole_Init called from cdInit to initialize CO routines */
void t_mole_global::init(void)
{
	static bool lgmole_Init_called=false;

	DEBUG_ENTRY( "Mole::init()" );

	/* prevent memory leaks */
	/* \todo	this is a temporary fix for PR14. We should improve the overall design
	 * of this code to prevent valid pointers being overwritten in a second call to mole_Init */
	if( lgmole_Init_called )
	{
		return;
	}

	/* say that we have been called */
	lgmole_Init_called = true;

	make_species();
	mole_make_list();
	mole_make_groups();
	mole.species.resize( mole_global.num_total );

	return;
}

void t_mole_global::zero(void)
{
	static bool lgFirstCall = true;
	static long int num_total_MALLOC=-1;
	
	DEBUG_ENTRY("t_mole_global::zero()");

	if( lgFirstCall )
	{
		lgFirstCall = false;
		num_total_MALLOC = mole_global.num_total;
	}
	else if( mole_global.num_total>num_total_MALLOC )
	{
		/* number of species has increased since last time - this can't happen
		 * tsuite / programs / comp4 has 95 first time, 98 second time */
		fprintf(ioQQQ,"DISASTER - the number of species in the CO network has increased.  This is not allowed.\n");
		fprintf(ioQQQ,"This could happen if an element was initially turned off or grains not included, then the element or grains was included.  There are not allowed.\n");
		fprintf(ioQQQ,"Sorry.\n");
		cdEXIT(EXIT_FAILURE);
	}

	for( long i=0; i<mole_global.num_total; ++i ) 
	{
		mole.species[i].zero();
	}
}
