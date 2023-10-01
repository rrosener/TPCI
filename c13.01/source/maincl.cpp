/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*main program that reads input and calls cloudy to compute a single model, or
 * try to optimize an observed model.  Routine returns 0 if model is ok, 
 * and 1 if problems occurred. */
#include "cddefines.h"
#include "cddrive.h"
#include "input.h"
#include "prt.h"
#include "save.h"
#include "called.h"
#include "monitor_results.h"
#include "mpi_utilities.h"
#include "grid.h"

exit_type cdMain( int argc, const char* argv[] );

inline void print_delimiter(long nOptimiz)
{
	fprintf( ioQQQ, " ************************************************** GRID_DELIMIT" );
	if( nOptimiz >= 0 )
		fprintf( ioQQQ, " -- grid%9.9ld", nOptimiz );
	fprintf( ioQQQ, "\n" );
}

/** main: this is a wrapper around cdMain. It takes care of the MPI stuff
 * for non-MPI runs, this should do nothing more than call cdMain and exit. */
int main( int argc, char *argv[] )
{
	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "main()" );

	if( cpu.i().lgMPI() )
	{
		MPI::Init( argc, argv );

		cpu.i().set_nCPU( MPI::COMM_WORLD.Get_size() );
		cpu.i().set_nRANK( MPI::COMM_WORLD.Get_rank() );

		// MPI::Init() will have overwritten our signal handlers
		// so we need to set them again....
		cpu.i().set_signal_handlers();
	}

	// this will generate input files for each grid point,
	// or execute the Phymir run, whichever is appropriate
	exit_status = cdMain( argc, (const char**)argv );

	// wait for writing of input files to finish
	if( cpu.i().lgMPI() )
		MPI::COMM_WORLD.Barrier();

	// process the individual grid points
	if( grid.lgGrid )
	{
		// this was set to true after we wrote the last input script
		grid.lgGridDone = false;

		// from now on each rank will run its own model
		cpu.i().set_MPISingleRankMode( true );

		load_balance lb( grid.totNumModels );

		// Each MPI rank will get jobs assigned by lb and execute them.
		// If there are no jobs left, lb.next_job() will return -1.
		while( ( optimize.nOptimiz = lb.next_job() ) >= 0 )
		{
			const char** mpi_argv = new const char*[argc+2];

			string jobName = GridPointPrefix( optimize.nOptimiz );
			for( int i=0; i < argc; ++i )
				mpi_argv[i] = argv[i];
			mpi_argv[argc] = "-g";
			mpi_argv[argc+1] = jobName.c_str();

			exit_type retval = cdMain( argc+2, mpi_argv );

			exit_status = max( retval, exit_status );
			delete[] mpi_argv;

			++grid.seqNum;
		}

		lb.finalize();

		// gather the spectral information from all ranks for saving FITS files
		grid.lgGridDone = true;
		optimize.nOptimiz = grid.totNumModels;
		GridGatherInCloudy();

		// and concatenate the output
		if( cpu.i().lgMaster() )
			process_output();
	}

	// remove empty output files from slave ranks
	if( cpu.i().lgMPI() && cpu.i().lgMaster() )
	{
		for( long n=1; n < cpu.i().nCPU(); ++n )
		{
			ostringstream oss;
			oss << ".err" << setfill('0') << setw(2) << n;
			string slave_output = save.chRedirectPrefix + oss.str();
			FILE *io = open_data( slave_output.c_str(), "a", AS_LOCAL_ONLY );
			bool lgEmpty = ( ftell(io) == 0 );
			fclose( io );
			if( lgEmpty )
				remove( slave_output.c_str() );
		}
	}

	if( cpu.i().lgMPI() )
		MPI::Finalize();

	return exit_status;
}

/** cdMain: this is the main entry point for Cloudy */
exit_type cdMain( int argc, const char* argv[] )
{
	/* these will be used to count number of various problems */
	long int NumberWarnings, 
	  NumberCautions, 
	  NumberNotes, 
	  NumberSurprises, 
	  NumberTempFailures, 
	  NumberPresFailures,
	  NumberIonFailures, 
	  NumberNeFailures;

	bool lgAbort_exit,
	  lgEarly_exit=true,
	  lgFileIO;
	/* number of lines we can still read in */
	int nread=0;

	int i;
	const char *s, 
	  *prefix = "",
	  *gprefix = "", // grid prefix
	  *pprefix = "", // save prefix
	  *rprefix = ""; // redirect prefix
	string infile("");

	/* the length of the following vector will be the longest line image
	 * the code will be able to read here.  Cloudy itself will ignore anything 
	 * beyond INPUT_LINE_LENGTH, and checks that no information exists beyond it. 
	 * The code will stop if the input line is longer than INPUT_LINE_LENGTH
	 * since extra characters would become a new command line due to buffer overrun */
	char chLine[INPUT_LINE_LENGTH];

	/* indicates that a command line flag to redirect I/O has been used */
	lgFileIO = false;

	exit_type exit_status = ES_SUCCESS;

	DEBUG_ENTRY( "cdMain()" );

	try {
		/* 
		 * Handle argument input -- written generally, but presently handles
		 * only `-p prefix' or `-pprefix' to set input file as `prefix.in',
		 * output file as `prefix.out' and the save prefix.
		 */
		for( i=1; i < argc; i++ ) 
		{
			s = argv[i];
			if( *s != '-' || s[1] == '\0' ) 
			{
				if (infile != "")
				{
					fprintf( ioQQQ, "%s: only one input file argument allowed\n", argv[0] );
					cdEXIT(ES_FAILURE);
				}
				infile = s;
				if( infile.find( cpu.i().chDirSeparator() ) != string::npos )
				{
					fprintf( ioQQQ, "%s %s: read/write from subdirectories is not supported\n",
								argv[0], infile.c_str() );
					cdEXIT(ES_FAILURE);
				}
				FILE *fp = fopen(infile.c_str(),"r");
				if (!fp)
				{
					fprintf( ioQQQ, "%s: input filename `%s' not found\n", argv[0], infile.c_str() );
					cdEXIT(ES_FAILURE);
				}
				fclose(fp);
				size_t suffindex = infile.find_last_of(".");
				if (suffindex != string::npos)
					infile = infile.substr(0,suffindex);
				pprefix = rprefix = infile.c_str();
				lgFileIO = true;
			}
			else
			{
				while( s != NULL && *(++s) )
				{
					exit_type exit = ES_SUCCESS;
					switch( *s ) 
					{
					case 'a':
						cpu.i().setAssertAbort( true );
						break;
					case 'g':
					case 'p':
					case 'r':
						if( s[1] != '\0' )
						{
							prefix = s+1;
						}					
						else
						{
							if( ++i == argc || argv[i][0] == '-' )
							{
								fprintf( ioQQQ, "%s: no argument given for -%c flag\n",
									 argv[0], *s );
								cdEXIT(ES_FAILURE);
							}
							prefix = argv[i];
							if( strchr(prefix, cpu.i().chDirSeparator()) != NULL )
							{
								fprintf( ioQQQ, "%s -%c %s: writing in subdirectories is not supported\n",
									 argv[0], *s, prefix );
								cdEXIT(ES_FAILURE);
							}
						}
						if( *s == 'g' )
							gprefix = prefix;
						else if( *s == 'p' )
						{
							pprefix = prefix;
							rprefix = prefix;
						}
						else if( *s == 'r' )
						{
							// make sure we erase the effects of a possible earlier -p flag
							pprefix = "";
							rprefix = prefix;
						}
						else
							TotalInsanity();
						s = NULL;
						lgFileIO = true;
						break;
					default:
						fprintf( ioQQQ, "%s: argument %d, `%s': flag -%c not understood\n",
							 argv[0], i, argv[i], *s );
						exit = ES_FAILURE;
					case 'h':
						fprintf( ioQQQ, "\nSupported flags are:\n\n" );
						fprintf( ioQQQ, "-p example\n" );
						fprintf( ioQQQ, "    Cloudy reads the input from the file example.in\n" );
						fprintf( ioQQQ, "    and writes the output to the file example.out.\n" );
						fprintf( ioQQQ, "    Additionally, all file names in SAVE commands are\n" );
						fprintf( ioQQQ, "    prepended with the string \"example\", e.g. the\n" );
						fprintf( ioQQQ, "    output of SAVE DR \".dr\" will be in example.dr.\n" );
						fprintf( ioQQQ, "-r example\n" );
						fprintf( ioQQQ, "    This does the same as the -p switch, except that\n" );
						fprintf( ioQQQ, "    the names used in SAVE commands are not altered.\n" );
						fprintf( ioQQQ, "-g example\n" );
						fprintf( ioQQQ, "    This switch is used internally in MPI grid runs.\n" );
						fprintf( ioQQQ, "    Normal users should not use this switch.\n" );
						fprintf( ioQQQ, "-a\n" );
						fprintf( ioQQQ, "    This switch is used in debugging. It causes the\n" );
						fprintf( ioQQQ, "    code to crash rather than exit gracefully after\n" );
						fprintf( ioQQQ, "    a failed assert. This flag is deprecated.\n" );
						fprintf( ioQQQ, "-h\n" );
						fprintf( ioQQQ, "    Print this message.\n" );
						cdEXIT(exit);
					}
				}
			}
		}

		/* initialize the code for this run */
		cdInit();

		save.chGridPrefix = gprefix;
		save.chFilenamePrefix = pprefix;
		save.chRedirectPrefix = rprefix;

		/* following should be set true to print to file instead of std output */
		if( lgFileIO )
		{
			string Base = save.chGridPrefix + save.chRedirectPrefix;
			string InName = Base + ".in";
			string OutName;
			if( cpu.i().lgMPI_talk() )
				OutName = Base + ".out";
			else
			{
				ostringstream oss;
				oss << ".err" << setfill('0') << setw(2) << cpu.i().nRANK();
				OutName = Base + oss.str();
			}
			cdInput( InName.c_str(), "r" );
			cdOutput( OutName.c_str(), "w" );
		}

		if( optimize.nOptimiz == 0 && called.lgTalk && cpu.i().lgMPISingleRankMode() )
			print_delimiter(-1);
			
		nread = 1;
		/* keep reading input lines until end of file */
		while( read_whole_line(chLine, (int)sizeof(chLine), ioStdin)!= NULL )
		{
			char *chChar;
			/* when running from command prompt, last line can be \n or 
			 * \r (this can happen with gcc under cygwin in windows) 
			 * or space, so check for each here, and break if present,
			 * check on chLine[23] is for case where we are reading cloudy output */
			/*fprintf(ioQQQ,"DEBUG char0 %i\n", (int)chLine[0] );*/
			if( chLine[0] == '\n' || chLine[0] == '\r' || 
			    ( chLine[0] == ' ' && chLine[23] != '*' ) ||
			    (strncmp( chLine , "***" , 3 )==0 ) )
				break;
			
			/* read entire line of input - lgAbort set true if lines too long */
			if( !lgInputComment(chLine) )
			{

				if( (chChar = strchr_s(chLine , '\"' ) ) ==NULL )
				{
					/* check for underscore, probably meant as a space */
					while( (chChar = strchr_s(chLine , '_' ) ) !=NULL )
					{
						*chChar = ' ';
						input.lgUnderscoreFound = true;
					}
				}

				/* change _, [, and ] to space if no filename occurs on line 
				 * >>chng 06 sep 04 use routine to check for comments 
				 * do not remove _, [, and ] in comments */
				/* check for left or right bracket, probably meant as a space */
				while( (chChar = strchr_s(chLine , '[' ) ) !=NULL )
				{
					*chChar = ' ';
					input.lgBracketFound = true;
				}

				while( (chChar = strchr_s(chLine , ']' ) ) !=NULL )
				{
					*chChar = ' ';
					input.lgBracketFound = true;
				}
			}

			/* this is trick so that cloudy input can read cloudy output */
			/* are first 25 char of input string equal to output? */
			if( strncmp(chLine,"                       * ",25) == 0 )
			{
				/* reading cloudy output, send in shifted input */
				nread = cdRead( chLine+25 );
			}
			else
			{
				/* stuff the command line into the internal stack */
				nread = cdRead( chLine  );
			}
		}

		if( lgAbort )
		{
			/* input parser hit something REALLY bad */
			cdEXIT(ES_CLOUDY_ABORT);
		}

		if( nread <= 0 )
			fprintf(ioQQQ," Warning: limit to number of lines exceeded, %i\n", nread);

		// optimize.lgVaryOn catches both optimizer and grid runs
		if( ( cpu.i().lgMPI() || optimize.lgVaryOn ) && save.chRedirectPrefix.empty() )
		{
			if( cpu.i().lgMaster() )
			{
				if( cpu.i().lgMPI() )
					fprintf( ioQQQ, " Please use the style \"mpirun -n np /path/to/cloudy.exe -r input\" when doing grid\n"
						"or optimizer runs.  See http://trac.nublado.org/wiki/RunCode for more information.\n" );
				else
					fprintf( ioQQQ, " Please use the style \"/path/to/cloudy.exe -r input\" when doing grid\n"
						"or optimizer runs.  See http://trac.nublado.org/wiki/RunCode for more information.\n" );
			}
			// stop the grid from being executed any further
			grid.lgGrid = false;
			cdEXIT(ES_FAILURE);
		}

		/* actually call the code.  This routine figures out whether the code will do
		 * a single model or be used to optimize on a spectrum, by looking for the
		 * keyword VARY on command lines.  It will call routine cloudy if no vary commands
		 * occur, and lgOptimize_do if VARY does occur.  
		 * cdDrive returns 0 if calculation is ok, 1 if problems happened */
		if( cdDrive() )
			exit_status = ES_FAILURE;

		/* the last line of output will contain some interesting information about the model*/
		cdNwcns(
			/* abort status, this better be false, 0 */
			&lgAbort_exit,
			/* the number of warnings, cautions, notes, and surprises */
			&NumberWarnings, 
			&NumberCautions, 
			&NumberNotes, 
			&NumberSurprises, 
			/* the number of temperature convergence failures */
			&NumberTempFailures, 
			/* the number of pressure convergence failures */
			&NumberPresFailures,
			/* the number of ionization convergence failures */
			&NumberIonFailures, 
			/* the number of electron density convergence failures */
			&NumberNeFailures );

		ostringstream finalMsg;

		finalMsg << " Cloudy ends: " << nzone << " zone";
		if( nzone > 1 )
			finalMsg << "s";

		finalMsg << ", " << iteration << " iteration";
		if( iteration > 1 )
			finalMsg << "s";

		if( lgAbort_exit )
			finalMsg << ", ABORT DISASTER PROBLEM";

		if( NumberWarnings > 0 )
		{
			finalMsg << ", " << NumberWarnings << " warning";
			if( NumberWarnings > 1 )
				finalMsg << "s";
			/* this indicates error */
			exit_status = ES_FAILURE;
		}

		if( NumberCautions > 0 )
		{
			finalMsg << ", " << NumberCautions << " caution";
			if( NumberCautions > 1 )
				finalMsg << "s";
		}

		/* this flag was set in lgCheckMonitors*/
		if( !lgMonitorsOK )
		{
			finalMsg << ", ";
			/* some botches were three sigma */
			if( lgBigBotch  )
				finalMsg << "BIG ";
			finalMsg << "BOTCHED MONITORS!!!";
			/* this indicates error */
			exit_status = ES_FAILURE;
		}

		if( NumberTempFailures+NumberPresFailures+NumberIonFailures+NumberNeFailures > 0 )
		{
			finalMsg << ". Failures: " << NumberTempFailures << " thermal, ";
			finalMsg << NumberPresFailures << " pressure, ";
			finalMsg << NumberIonFailures << " ionization, ";
			finalMsg << NumberNeFailures << " electron density";
		}

		/* NB DO NOT CHANGE ANY ASPECT OF THE FOLLOWING STRINGS - THEY ARE USED TO RECORD
		 * EXEC TIME BY A PERL SCRIPT */
		if( prt.lgPrintTime )
		{
			/* print execution time [s] by default,
			 * need spaces around number so that logging perl script picks up correct number 
			 * ir_extime.pl script will delete through "ExecTime(s)" and remainder of line must be number */
			finalMsg << ".  ExecTime(s) " << fixed << setprecision(2) << cdExecTime();
		}

		if( called.lgTalk )
			fprintf( ioQQQ, "%s\n", finalMsg.str().c_str() );

		lgEarly_exit = false;

		/* cdDrive returned 1 if something bad happened, and 0 if everything is ok.  We will
		 * return 0 if everything is ok, and a non-zero error code if something bad happened.*/
		cdEXIT(exit_status);
	}
	catch( bad_alloc )
	{
		fprintf( ioQQQ, " DISASTER - A memory allocation has failed. Most likely your computer "
			 "ran out of memory.\n Try monitoring the memory use of your run. Bailing out...\n" );
		exit_status = ES_BAD_ALLOC;
	}
	catch( out_of_range& e )
	{
		fprintf( ioQQQ, " DISASTER - An out_of_range exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_OUT_OF_RANGE;
	}
	catch( bad_assert& e )
	{
		MyAssert( e.file(), e.line() , e.comment() );
		exit_status = ES_BAD_ASSERT;
	}
#ifdef CATCH_SIGNAL
	catch( bad_signal& e )
	{
		if( ioQQQ != NULL )
		{
			if( e.sig() == SIGINT || e.sig() == SIGQUIT )
			{
				fprintf( ioQQQ, " User interrupt request. Bailing out...\n" );
				exit_status = ES_USER_INTERRUPT;
			}
			else if( e.sig() == SIGTERM )
			{
				fprintf( ioQQQ, " Termination request. Bailing out...\n" );
				exit_status = ES_TERMINATION_REQUEST;
			}
			else if( e.sig() == SIGILL )
			{
				fprintf( ioQQQ, " DISASTER - An illegal instruction was found. Bailing out...\n" );
				exit_status = ES_ILLEGAL_INSTRUCTION;
			}
			else if( e.sig() == SIGFPE )
			{
				fprintf( ioQQQ, " DISASTER - A floating point exception occurred. Bailing out...\n" );
				exit_status = ES_FP_EXCEPTION;
			}
			else if( e.sig() == SIGSEGV )
			{
				fprintf( ioQQQ, " DISASTER - A segmentation violation occurred. Bailing out...\n" );
				exit_status = ES_SEGFAULT;
			}
#			ifdef SIGBUS
			else if( e.sig() == SIGBUS )
			{
				fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
				exit_status = ES_BUS_ERROR;
			}
#			endif
			else
			{
				fprintf( ioQQQ, " DISASTER - A signal %d was caught. Bailing out...\n", e.sig() );
				exit_status = ES_UNKNOWN_SIGNAL;
			}

		}
	}
#endif
	catch( cloudy_exit& e )
	{
		if( called.lgTalk )
		{
			ostringstream oss;
			oss << " [Stop in " << e.routine();
			oss << " at " << e.file() << ":" << e.line();
			if( e.exit_status() == 0 )
				oss << ", Cloudy exited OK]";
			else
				oss << ", something went wrong]";
			fprintf( ioQQQ, "%s\n", oss.str().c_str() );
		}
		exit_status = e.exit_status();
		if( exit_status == ES_FAILURE && !lgEarly_exit )
		{
			// try to make the error code more descriptive
			// when there is no early exit in the code, then these 3
			// should be the only reasons for a non-zero exit code
			if( NumberWarnings > 0 )
				exit_status = ES_WARNINGS;
			if( !lgMonitorsOK )
				exit_status = ES_BOTCHES;
			if( lgAbort_exit )
				exit_status = ES_CLOUDY_ABORT;
		}
	}
	catch( std::exception& e )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
			 e.what() );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}
	// generic catch-all in case we forget any specific exception above... so this MUST be the last one.
	catch( ... )
	{
		fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
		exit_status = ES_UNKNOWN_EXCEPTION;
	}

	if( called.lgTalk && cpu.i().lgMPISingleRankMode() )
		print_delimiter(optimize.nOptimiz);

	cdPrepareExit(exit_status);

	return exit_status;
}
