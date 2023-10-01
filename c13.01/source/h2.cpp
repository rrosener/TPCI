/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "h2.h"
#include "h2_priv.h"
#include "hmi.h"
	
vector<diatomics*> diatoms;

diatomics h2("h2", 4100., &hmi.H2_total, Yan_H2_CS);
diatomics hd("hd", 4100., &hmi.HD_total, Yan_H2_CS);

void diatoms_init( void )
{
	DEBUG_ENTRY( "diatoms_init()" );

	diatoms.clear();
	diatoms.push_back( &h2 );
	diatoms.push_back( &hd );

	// H
	h2.coll_source[0].magic = 110416;
	h2.coll_source[0].filename = "coll_rates_H_07.dat";
	// He
	h2.coll_source[1].magic = 110416;
	h2.coll_source[1].filename = "coll_rates_He_ORNL.dat";
	// H2 ortho
	h2.coll_source[2].magic = 110416;
	h2.coll_source[2].filename = "coll_rates_H2ortho_ORNL.dat";
	// H2 para	
	h2.coll_source[3].magic = 110416;
	h2.coll_source[3].filename = "coll_rates_H2para_ORNL.dat";
	// proton
	h2.coll_source[4].magic = 110416;
	h2.coll_source[4].filename = "coll_rates_Hp.dat";

	// H
	hd.coll_source[0].magic = 110416;
	hd.coll_source[0].filename = "coll_rates_H_07.dat";
	// He
	hd.coll_source[1].magic = 110416;
	hd.coll_source[1].filename = "coll_rates_He_LeBourlot.dat";
	// H2 ortho
	hd.coll_source[2].magic = 110416;
	hd.coll_source[2].filename = "coll_rates_H2ortho_LeBourlot.dat";
	// H2 para	
	hd.coll_source[3].magic = 110416;
	hd.coll_source[3].filename = "coll_rates_H2para_LeBourlot.dat";
	// proton
	hd.coll_source[4].magic = 110416;
	hd.coll_source[4].filename = "coll_rates_Hp.dat";

	return;
}

