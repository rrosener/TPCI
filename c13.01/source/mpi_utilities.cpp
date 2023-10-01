/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "mpi_utilities.h"
#include "save.h"
#include "dynamics.h"
#include "grid.h"

#ifndef MPI_ENABLED

namespace MPI
{
	t_MPI COMM_WORLD;
}

#endif

// NB NB this routine cannot throw any exceptions as it is executed outside
// the try{} block -- this includes mechanisms like ASSERT and cdEXIT!
void load_balance::init(int nJobs)
{
	if( nJobs <= 0 )
		return;

	bool lgMPI = cpu.i().lgMPI();

	p_jobs.resize( nJobs );
	if( lgMPI )
		p_ptr = MPI::COMM_WORLD.Get_rank();
	else
		p_ptr = 0;
	// the master rank will now set up a random sequence for the jobs
	// this way we hope to get statistical load balancing of the ranks
	if( p_ptr == 0 )
	{
		for( int i=0; i < nJobs; ++i )
			p_jobs[i] = i;

		if ( lgMPI )
		{
			// This may or may not seed the random number generator used by
			// random_shuffle. There is no portable C++ interface to do this :-(
			srand( unsigned( time(NULL) ) );
			random_shuffle( p_jobs.begin(), p_jobs.end() );
		}
	}
	// now broadcast the random sequence to the other ranks...
	if( lgMPI )
		MPI::COMM_WORLD.Bcast( &p_jobs[0], nJobs, MPI::type(p_jobs[0]), 0 );
}

STATIC void check_grid_file( const string& fnam, int j, int ipPun );

/** process_output: concatenate output files produced in MPI grid run */
void process_output()
{
	DEBUG_ENTRY( "process_output()" );

	// NOTE: when this routine is called all file handles have already been closed

	try
	{
		string main_input = save.chRedirectPrefix + ".in";
		string main_output = save.chRedirectPrefix + ".out";

		// first process main output files
		FILE* main_output_handle = open_data( main_output.c_str(), "a", AS_LOCAL_ONLY );
		for( int j=0; j < grid.totNumModels; ++j )
		{
			string Base = GridPointPrefix(j) + save.chRedirectPrefix;
			string in_name = Base + ".in";
			remove( in_name.c_str() );
			string out_name = Base + ".out";
			append_file( main_output_handle, out_name.c_str() );
			remove( out_name.c_str() );
		}
		fclose( main_output_handle );

		fstream main_input_handle;
		open_data( main_input_handle, main_input.c_str(), mode_r, AS_LOCAL_ONLY );
		string line;

		int ipPun = 0;

		while( getline( main_input_handle, line ) )
		{
			string caps_line;
			// create all caps version
			string::const_iterator p = line.begin();
			while( p != line.end() )
				caps_line.push_back( toupper(*p++) );
			if( caps_line.compare( 0, 4, "SAVE" ) == 0 || caps_line.compare( 0, 4, "PUNC" ) == 0 )
			{
				ASSERT( ipPun < save.nsave );
				string fnam = save.chFilenamePrefix;
				string::size_type p = line.find( '"' );
				fnam += line.substr( ++p );
				fnam.erase( fnam.find( '"' ) );
				// first do a minimal check on the validity of the save files
				for( int j=0; j < grid.totNumModels; ++j )
					check_grid_file( fnam, j, ipPun );
				// open in binary mode in case we are writing a FITS file
				FILE *dest = open_data( fnam.c_str(), "ab", AS_LOCAL_ONLY_TRY );
				if( dest != NULL )
				{
					if( save.lgSaveToSeparateFiles[ipPun] )
					{
						// keep the save files for each grid point separate
						// the main save file contains the save header
						// salvage it by prepending it to the first save file
						// this gives the same behavior as in non-MPI runs
						string gridnam = GridPointPrefix(0) + fnam;
						append_file( dest, gridnam.c_str() );
						fclose( dest );
						dest = NULL;
						// this will overwrite the old file gridnam
						rename( fnam.c_str(), gridnam.c_str() );
					}
					else
					{
						// concatenate the save files for each grid point
						for( int j=0; j < grid.totNumModels; ++j )
						{
							string gridnam = GridPointPrefix(j) + fnam;
							append_file( dest, gridnam.c_str() );
							remove( gridnam.c_str() );
						}
					}
					if( caps_line.find( "XSPE", 4 ) != string::npos )
					{
						// dest points to an empty file, so generate the complete FITS file now
						ASSERT( save.FITStype[ipPun] >= 0 &&
							save.FITStype[ipPun] < NUM_OUTPUT_TYPES );
						saveFITSfile( dest, save.FITStype[ipPun] );
						fseek( dest, 0, SEEK_END );
						ASSERT( ftell(dest)%2880 == 0 );
					}
					if( dest != NULL )
					{
						fclose( dest );
					}
				}
				else
				{
					fprintf( ioQQQ, "PROBLEM - could not open file %s\n", fnam.c_str() );
				}
				++ipPun;
			}
		}
	}
	catch( ... )
	{
		fprintf( ioQQQ, "PROBLEM - an internal error occurred while post-processing the grid output\n" );
	}
}

/** check_grid_file: check whether the save file is present and is terminated by a GRID_DELIMIT string */
STATIC void check_grid_file( const string& fnam, int j, int ipPun )
{
	DEBUG_ENTRY( "check_grid_file()" );

	// these are binary files, don't touch them...
	if( save.lgFITS[ipPun] )
		return;

	bool lgForceNoDelimiter = false;
	// in these cases there should not be a GRID_DELIMIT string...
	if( !save.lgHashEndIter[ipPun] || !save.lg_separate_iterations[ipPun] ||
	    dynamics.lgTimeDependentStatic || strcmp( save.chHashString , "TIME_DEP" ) == 0 )
		lgForceNoDelimiter = true;

	bool lgAppendDelimiter = true;
	bool lgAppendNewline = false;
	string gridnam = GridPointPrefix(j) + fnam;
	fstream str;
	open_data( str, gridnam.c_str(), mode_r, AS_LOCAL_ONLY_TRY );
	if( str.is_open() )
	{
		str.seekg( 0, ios_base::end );
		if( str.good() && str.tellg() > 0 )
		{
			// check if the file ends in a newline
			str.seekg( -1, ios_base::cur );
			char chr;
			str.get( chr );
			lgAppendNewline = ( chr != '\n' );
			// check if the GRID_DELIMIT string is present
			string line;
			str.seekg( 0, ios_base::beg );
			while( getline( str, line ) )
			{
				if( line.find( "GRID_DELIMIT" ) != string::npos )
					lgAppendDelimiter = false;
			}
		}
		str.close();
	}
	if( lgForceNoDelimiter )
		lgAppendDelimiter = false;
	if( lgAppendNewline || lgAppendDelimiter )
	{
		open_data( str, gridnam.c_str(), mode_a, AS_LOCAL_ONLY_TRY );
		if( str.is_open() )
		{
			if( lgAppendNewline )
				str << endl;
			if( lgAppendDelimiter )
			{
				str << save.chHashString << " GRID_DELIMIT -- grid";
				str << setfill( '0' ) << setw(9) << j << endl;
			}
			str.close();
		}
	}
}

/** append_file: append output produced on file <source> to open file descriptor <dest> */
void append_file( FILE *dest, const char *source )
{
	DEBUG_ENTRY( "append_file()" );

	FILE *src = open_data( source, "rb", AS_LOCAL_ONLY_TRY );
	if( src == NULL )
		return;

	// limited testing shows that using a 4 KiB buffer should
	// give performance that is at least very close to optimal
	// tests were done by concatenating 10 copies of a 62.7 MiB file
	const size_t BUF_SIZE = 4096;
	char buf[BUF_SIZE];

	while( ! feof(src) )
	{
		size_t nb = fread( buf, sizeof(char), BUF_SIZE, src );
		fwrite( buf, sizeof(char), nb, dest );
	}
	fclose(src);
	return;
}
