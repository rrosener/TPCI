/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#include "cddefines.h"
#include "quantumstate.h"

/*StateJunk set all elements of transition struc to dangerous values */
void Junk(qStateProxy st )
{

	DEBUG_ENTRY( "qStateProxy::Junk()" );

	st.chLabel()[0] = '\0';

	/** statistical weight [dimensionless] */
	st.g() = -FLT_MAX;

	/** column density of state [cm-2] */
	st.ColDen() = -FLT_MAX;

	/** population of state [cm-3] */
	st.Pop() = -FLT_MAX;

	 /** ion stage of element, 1 for atom, 2 ion, etc */
	st.IonStg() = -10000;

	 /** atomic number of element, 1 for H, 2 for He, etc */
	st.nelem() = -10000;

	st.n()=st.l()=st.S()=st.j()=st.v()=st.J()=-1;

	return;
}

/*StateZero zeros out the structure */
void Zero( qStateProxy st)
{

	DEBUG_ENTRY( "qStateProxy::Zero()" );

	/** population of state [cm-3] */
	st.Pop() = 0.;
	return;
}
