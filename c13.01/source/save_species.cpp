/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*SaveSpecies generate output for the save species command */
#include "cddefines.h"
#include "opacity.h"
#include "taulines.h"
#include "radius.h"
#include "phycon.h"
#include "save.h"
#include "mole.h"

/* save results for one particular species */
STATIC void SaveSpeciesOne( const size_t species_index, const char chKey[],
		FILE *ioPUN, long int ipPun, size_t maxLevels );

/*SaveSpecies generate output for the save species command */
void SaveSpecies(
		FILE* ioPUN,
		long int ipPun)
{
	DEBUG_ENTRY( "SaveSpecies()" );

	if( strcmp( save.chSaveArgs[ipPun], "LABE" )==0 )
	{
		if( save.lgPunHeader[ipPun] )
		{
			/* save list of species labels */
			fprintf( ioPUN, "#Species labels\n" );
			save.lgPunHeader[ipPun] = false;
			for( size_t i=0; i<mole_global.list.size(); ++i )
			{
				molecule *spg = &(*mole_global.list[i]);
				fprintf( ioPUN, "%s\n", spg->label.c_str() );
			}
		}
		return;
	}

	else if( strcmp(save.chSaveArgs[ipPun],"LEVL") == 0 )
	{
		/* number of levels active in this zone */
		if( save.lgPunHeader[ipPun] )
		{
			/* save list of species labels */
			fprintf( ioPUN, "#Species\tnumber of levels\n" );
			save.lgPunHeader[ipPun] = false;
		}
		for( size_t i=0; i<mole_global.list.size(); ++i )
		{
			molecule *spg = &(*mole_global.list[i]);
			molezone *sp = &mole.species[i];
			fprintf( ioPUN, "%s", spg->label.c_str() );
			if( sp->levels == NULL )
				fprintf( ioPUN, "\t%4lu\n", 0UL );
			else
				fprintf( ioPUN, "\t%4lu\n", (unsigned long)sp->levels->size() );
		}
		return;
	}

	/* remaining options are column densities, populations, and energies
	 * first branch; save results for all species if "" */
	if( strcmp( save.chSaveSpecies[ipPun], "" ) == 0 )
	{
		// max number of levels, for header print
		size_t mostLevels = 0;
		for( size_t i=0; i<mole_global.list.size(); ++i )
		{
			molezone *sp = &mole.species[i];
			if( sp->levels != NULL )
				mostLevels = MAX2(mostLevels, sp->levels->size() );
		}
		ASSERT( mostLevels > 1 );
		ASSERT( mostLevels < 10000 );

		// loop over species
		for( size_t i=0; i<mole_global.list.size(); ++i )
			SaveSpeciesOne( i, save.chSaveArgs[ipPun], ioPUN, ipPun, mostLevels );
	}
	else
	{
		const molecule *saveSpeciesGlobal = findspecies(save.chSaveSpecies[ipPun]);
		const molezone *saveSpecies = findspecieslocal(save.chSaveSpecies[ipPun]);

	  	if( saveSpecies == null_molezone )
		{
			fprintf( ioQQQ,"Could not find species %s, so SAVE SPECIES LABELS to get a list of all species."
					"\nSorry.\n", save.chSaveSpecies[ipPun] );
			cdEXIT(EXIT_FAILURE);
		}
	
		size_t numLevels = 0;
		if( saveSpecies->levels != NULL )
			numLevels = saveSpecies->levels->size();	
		SaveSpeciesOne( saveSpeciesGlobal->index, save.chSaveArgs[ipPun], ioPUN, ipPun, numLevels );
	}

	return;
}

/* print 0.000e+00 as simply 0 */
STATIC void PrintShortZero( FILE *ioPUN , double arg )
{
	DEBUG_ENTRY( "PrintShortZero()" );
	if( arg==0. )
		fprintf(ioPUN,"\t0");
	else
		fprintf(ioPUN,"\t%.3e", arg);

}

/* save results for one particular species */
STATIC void SaveSpeciesOne( const size_t species_index, const char chKey[],
		FILE *ioPUN, long int ipPun, size_t maxLevels )
{
	DEBUG_ENTRY( "SaveSpeciesOne()" );

	molecule *spg = &(*mole_global.list[species_index]);
	molezone *sp = &mole.species[species_index];

	if( spg == null_mole || sp == null_molezone )
		return;
	
	// one time print of energy levels
	if( strcmp( chKey, "ENER" )==0 )
	{
		if( save.lgPunHeader[ipPun] )
		{
			save.lgPunHeader[ipPun] = false;

			fprintf( ioPUN, "#species energies");
			for( size_t i = 0; i < maxLevels; ++i )
			{
				fprintf( ioPUN, "\t%lu", (unsigned long)i );
			}
			fprintf( ioPUN, "\n");
		}

		fprintf( ioPUN, "%s", spg->label.c_str() );
		if( sp->levels == NULL || sp->levels->size() == 0 )
		{
			fprintf( ioPUN, "\t%.6e", 0. );
		}
		else
		{
			for( qList::const_iterator st = sp->levels->begin(); st != sp->levels->end(); ++st )
			{
				ASSERT( (*st).g() > 0.f );
				fprintf( ioPUN, "\t%.6e", AnuUnit( (*st).energy().Ryd() ) );
			}
		}
		fprintf( ioPUN, "\n");
		return;
	}

	if( strcmp( chKey, "POPU" )==0 )
	{
		if( save.lgPunHeader[ipPun] )
		{
			fprintf( ioPUN, "#depth [cm] species populations [cm-3]");

			for( size_t i = 0; i < maxLevels; ++i )
			{
				fprintf( ioPUN, "\t%lu", (unsigned long)i );
			}
			fprintf( ioPUN, "\n");
			save.lgPunHeader[ipPun] = false;
		}

		fprintf( ioPUN, "%.5e", radius.depth_mid_zone );
		fprintf( ioPUN, "\t%s", spg->label.c_str() );

		if( sp->levels == NULL || sp->levels->size() == 0 )
		{
			PrintShortZero( ioPUN, sp->den );
		}
		else
		{
			bool lgZeroHit = false;
			// loop over levels for this species
			for( qList::const_iterator st = sp->levels->begin(); st != sp->levels->end(); ++st )
			{
				// don't print high levels which can have 0 abundance at low temperature
				if( !lgZeroHit )
					PrintShortZero( ioPUN, (*st).Pop() );
				if( (*st).Pop() == 0.)
					lgZeroHit = true;
			}
		}
		fprintf( ioPUN, "\n");
	}
	else if( strcmp( chKey, "COLU" )==0 )
	{

		if( save.lgPunHeader[ipPun] )
		{
			fprintf( ioPUN, "#species column density [cm-2]");

			for( size_t i = 0; i < maxLevels; ++i )
			{
				fprintf( ioPUN, "\t%lu", (unsigned long)i );
			}
			fprintf( ioPUN, "\n");
			save.lgPunHeader[ipPun] = false;
		}

		fprintf( ioPUN, "%s", spg->label.c_str() );

		if( sp->levels == NULL || sp->levels->size() == 0 )
		{
			PrintShortZero( ioPUN, sp->column );
		}
		else
		{
			// loop over levels
			bool lgZeroHit = false;
			for( qList::const_iterator st = sp->levels->begin(); st != sp->levels->end(); ++st )
			{
				// don't print high levels which can have 0 abundance at low temperature
				if( !lgZeroHit )
					PrintShortZero( ioPUN, (*st).ColDen() );
				if( (*st).ColDen() == 0.)
					lgZeroHit = true;
			}
		}
		fprintf( ioPUN, "\n");
	}
	
	return;
}
