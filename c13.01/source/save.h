/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef SAVE_H_
#define SAVE_H_

#include "energy.h"
#include "h2_priv.h"

/* save.h */
static const long LIMPUN = 100L;
static const long MAX_HEADER_SIZE = 20000L;

/** magic version number for the transmitted continuum output file */
static const long VERSION_TRNCON = 20100901L;

/**SaveDo produce save output during calculation 
\param chTime chTime is null terminated 4 char string, either "MIDL" or "LAST" 
*/
void SaveDo(
	const char *chTime); 

/** save one line, called by SaveLineStuff 
\param t
\param io
\param xLimit
\param index
\param DopplerWidth
*/
void Save1Line( 
	const TransitionProxy & t , 
	FILE * io , 
	realnum xLimit , 
	long index,
	realnum DopplerWidth);

/**SaveLineData punches selected line data for all lines transferred in code 
\param io
*/
NORETURN void SaveLineData(
	FILE * io);

/**save_opacity save total opacity in any element, save opacity command 
\param io
\param np
*/
void save_opacity(
	FILE * io, 
	long int np);

/**SaveSpecial generate output for the save special command 
\param io
\param chTime
*/
void SaveSpecial(
	FILE* io , 
  const char *chTime);

/**SaveSpecies generate output for the save species command 
\param ioPUN io unit for save
\param ipPun arrayin save array
*/
void SaveSpecies(
		FILE* ioPUN,
		long int ipPun );

/**Save1LineData data for save one line 
\param t
\param io
\param lgCS_2 this flag says whether collision strength should be saveed - should be false
	for multi level atoms since sums are not done properly 
\param lgPrint true print header line
*/
void Save1LineData( 
	const TransitionProxy & t , 
	FILE * io,
	 bool lgCS_2 ,
	 bool &lgPrint);

/**save_line do the save output 
\param ip the file we will write to 
\param chDo
\param intrinsic or emergent emission
*/
void save_line(
	FILE * ip,
  const char *chDo, 
  bool lgEmergent);

/**save_average parse save average command, or actually do the save output 
\param ipPun - array index for file for final save output 
*/
void save_average( 
	long int ipPun);

/** save_colden parse save column density command, or actually do the 
save lines output 
\param ip the file we will write to 
*/
void save_colden(
  /* the file we will write to */
  FILE * ioPUN );

/**Save_Line_RT parse the save line rt command - read in a set of lines 
\param ip the file we will write to
\param *chJob
*/
void Save_Line_RT(
	FILE * ip);

/** Save spectra to a FITS compatible file. 
\param io
\param option
*/
void saveFITSfile( 
	FILE* io, 
	int option );

/**SaveHeat save contributors to local heating, with save heat command, called by save_do 
\param io
*/
void SaveHeat(FILE* io);

/**CoolSave save coolants, called by save_do 
\param io
*/
void CoolSave(FILE * io, char chJob[]);

/**SaveGrid implements the SAVE GRID command
\param pnunit : file descriptor to write to
\param status : exit status of the grid point
*/
void SaveGrid(FILE* pnunit, exit_type status);

struct t_save {

	t_save()
	{
		for( long i=0; i < LIMPUN; ++i )
		{
			nLineList[i] = -1;
			lgFITS[i] = false;
			FITStype[i] = -1;
			nAverageList[i] = -1;
		}
	}

	~t_save()
	{
		for( long i=0; i < LIMPUN; ++i )
		{
			SaveLineListFree(i);
			SaveAverageFree(i);
		}
	}

	void SaveLineListFree(long i)
	{
		for( unsigned j=0; j < chLineListLabel[i].size(); ++j )
			delete[] chLineListLabel[i][j];
		chLineListLabel[i].clear();
		wlLineList[i].clear();
	}
	void SaveAverageFree(long i)
	{
		for( unsigned j=0; j < chAverageType[i].size(); ++j )
			delete[] chAverageType[i][j];
		chAverageType[i].clear();
		for( unsigned j=0; j < chAverageSpeciesLabel[i].size(); ++j )
			delete[] chAverageSpeciesLabel[i][j];
		chAverageSpeciesLabel[i].clear();
	}

	/** following are for save LineList option
	 * nLineList is number of em lines, -1 if not defined */
	long nLineList[LIMPUN];
	/** chLineListLabel is label for line list */
	vector<char*> chLineListLabel[LIMPUN];
	/** wlLineList is set of emission lines for LineList */
	vector<realnum> wlLineList[LIMPUN];
	/** flag saying whether to take ratio (true) of pairs */
	bool lgLineListRatio[LIMPUN];

	/** following are for save averages option
	* nAverageList is number of averages, -1 if not defined */
	long nAverageList[LIMPUN];
	/** chAverageType is label for type of average */
	vector<char*> chAverageType[LIMPUN];
	/** chAverageSpeciesLabel is label for species */
	vector<char*> chAverageSpeciesLabel[LIMPUN];
	/** nAverageIonList is set of ions for averages */
	vector<int> nAverageIonList[LIMPUN];
	/** nAverage2ndPar is set of second parameters for averages */
	vector<int> nAverage2ndPar[LIMPUN];

	/** this is the file where we will direct the output */
	FILE *ipPnunit[LIMPUN]; 

	/** option to append instead of overwrite */
	bool lgNoClobber[LIMPUN];

	/** these are for some special cases, same purpose as previous no clobber */
	bool lgPunConv_noclobber;
	bool lgDROn_noclobber;
	bool lgPunPoint_noclobber;
	bool lgioRecom_noclobber;
	bool lgQHSaveFile_noclobber;
	bool lgTraceConvergeBase_noclobber;
	bool lgSaveGrid_noclobber;

	/** is this a real save command, or one of the similar options like
	 * save dr, which is not done in save files */
	bool lgRealSave[LIMPUN]; 

	/** emergent or intrinsic line intensities */
	bool lgEmergent[LIMPUN];

	/** cumulative line intensities for time dependent sim */
	bool lgCumulative[LIMPUN];

	/** number of save commands entered */
	long int nsave;

	/**chSave - what is it we want to save? set in GetPunch, used in DoPunch */
	char chSave[LIMPUN][5];
	diatomics *whichDiatomToPrint[LIMPUN];

	/** which opacity to save out  */
	char chOpcTyp[LIMPUN][5];

	char chHeader[LIMPUN][MAX_HEADER_SIZE];

	const char *chNONSENSE;

	/** this flag tells us whether to save results of a grid to separate files
	 * for each grid point or all to the same file.  Different for different 
	 * save commands */
	bool lgSaveToSeparateFiles[LIMPUN];

	/** option to not insert end-of-iteration separator - used for save files
	 * that create one line per iteration */
	bool lg_separate_iterations[LIMPUN];

	/** flag saying whether we should save headers.  Used in grid punches so that
	 * the header only gets punched once. */
	bool lgPunHeader[LIMPUN];

	/** flag saying whether any save continuum commands were entered
	 * set true in parsecontinuum when save continuum entered,
	 * used in PrtComment to warn if continuum punched with no iterations */
	bool lgPunContinuum;

	/** punarg is set of optional arguments for the save command */
	realnum punarg[LIMPUN][3];

	/** punarg is set of optional arguments for the save command */
	string optname[LIMPUN];

	/** implement save every option - lgSaveEveryZone true if want to save
	 *every zone, nSaveEveryZone is number of zones to save */
	bool lgSaveEveryZone[LIMPUN];
	long int nSaveEveryZone[LIMPUN];

	/** set of optional arguments for save command, but as a string */
	char chSaveArgs[LIMPUN][5];

	char chSaveSpecies[LIMPUN][CHARS_SPECIES];

	/** lg flag lgPunLstIter for this save option,
	 do we only want to save on last iteration? */
	bool lgPunLstIter[LIMPUN];

	/** flag saying that this save file is in FITS format */
	bool lgFITS[LIMPUN];

	/** which FITS type is in this file */
	int FITStype[LIMPUN];

	/** option to say whether any FITS output should be punched, 
	 * initialized to false, but turned to true after last grid exec.  */
	bool lgPunchFits;

	/**chConPunEnr - units of continuum in save output */
	const char *chConPunEnr[LIMPUN];

	/** this global variable is index of save command loop in dopunch */
	long int ipConPun;

	/** should hash marks be printed after every iteration?
	 * default is yes, set no with no hash option on save command */
	bool lgHashEndIter[LIMPUN];

	/** this is the hash string, normally a set of hash marks, can be reset with
	 * set save hash command */
	char chHashString[INPUT_LINE_LENGTH];

	/** flush file after every iteration */
	bool lgFLUSH;

	/** this is a prefix that will be used at the start of all file names
	 * when doing an MPI grid run, normally an empty string, set with -g flag */
	string chGridPrefix;

	/** this is a prefix that will be set at the start of all save file names
	 * normally an empty string, set with PUNCH PREFIX command or -p flag */
	string chFilenamePrefix;

	/** this is the prefix that will be set at the start of the input and output file
	 * normally an empty string, set with -p or -r flag */
	string chRedirectPrefix;

	/** this is the name of the main output file we are writing to
	 * it is set by calling cdOutput() */
	string chOutputFile;

	/** set with save line intensities and save results commands,
	 * says whether results arrays produced by routine PunResults1Line should
	 * be column or array */
	char chPunRltType[7];

	/** option to save out pointers with save pointers command
	 * ipPoint is save file handle, lgPunPoint says whether we will do it
	 */
	FILE* ipPoint;
	bool lgPunPoint;

	/** unit number, and flag, for saving reason for continued iterations */
	bool lgPunConv;
	FILE* ipPunConv;

	/** these control saving choice of dr - this is not really a save command
	 * ipDRout is io unit, lgDROn says saving dr logic has been set,
	 * and lgDRPLst says to save the last iteration */
	FILE * ipDRout;
	bool lgDROn, 
	  lgDRPLst,
	  lgDRHash;

	/* set true save convergence base */
	bool lgTraceConvergeBase,
		lgTraceConvergeBaseHash;
	FILE * ipTraceConvergeBase;

	/** option to save recombination coefficients to external file */
	FILE* ioRecom;
	bool lgioRecom;

	/** option to save line intensities every for every zone
	 * logical variable says whether LinEvery was set */
	long int LinEvery;
	bool lgLinEvery;

	/** set skip sets this variable, which says how many cells to skip in save */
	long int ncSaveSkip;

	/** threshold for faintest cooling or heating to be punched
	 * default is set to 0.05 in scalar, can be reset with 'set weakheatcool' */
	realnum WeakHeatCool;

	/** contrast factor for lines to continuum in save output
	* default is 1 (gives correct line intensities) and changed with
	* set width command (enters width in km/sec) */
	realnum Resolution;
	// same thing but for absorption lines - default is unity, set to
	// SaveLWidth if ABSORPTION keyword occurs on save line width command
	realnum ResolutionAbs;

	char chSpeciesDominantRates[LIMPUN][CHARS_SPECIES];

	// the frequency at which the continuum volume emissivity should be saved
	Energy emisfreq[LIMPUN];
	long ipEmisFreq[LIMPUN];

};

extern t_save save;

#endif /* SAVE_H_ */
