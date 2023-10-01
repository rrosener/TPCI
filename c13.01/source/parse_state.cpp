/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
/*ParseState save or recover previous state of the code */
#include "cddefines.h"
#include "parse.h"
#include "state.h"
#include "parser.h"

/*ParseState save or recover previous state of the code */
void ParseState(Parser &p)
{
	char chFilename[INPUT_LINE_LENGTH];

	DEBUG_ENTRY( "ParseState()" );

	/* 
	 * get file name for this save output.
	 * GetQuote does the following -
	 * first copy original version of file name into chLabel, 
	 * string does include null termination.
	 * set filename in OrgCard and second parameter to spaces so 
	 * that not picked up below as keyword
	 * last parameter says to abort if no quote found 
	 */
	p.GetQuote( chFilename , true );

	/* option to print all contents of arrays - BIG PRINTOUT! */
	if( p.nMatch("PRIN") )
		state.lgState_print = true;

	if( p.nMatch(" GET") )
	{
#		if 0
		state.ioGET_STATE = open_data( chFilename, "rb", AS_LOCAL_ONLY );
#		endif
		state.lgGet_state = true;
		strcpy( state.chGetFilename , chFilename );
	}
	else if( p.nMatch(" PUT") )
	{
#		if 0
		state.ioPUT_STATE = open_data( chFilename , "wb", AS_LOCAL_ONLY );
#		endif
		state.lgPut_state = true;
		strcpy( state.chPutFilename , chFilename );
		/* look for keyword ALL - says want to save state for all iterations,
		 * default is last iteration */
		if( p.nMatch(" ALL") )
		{
			state.lgPutAll = true;
		}
		else
		{
			state.lgPutAll = false;
		}
	}

	else
	{
		fprintf( ioQQQ, " The STATE command has two keywords, GET and PUT.  One must appear - I did not see it.\n Sorry.\n" );
		cdEXIT(EXIT_FAILURE);
	}
	return;
}
