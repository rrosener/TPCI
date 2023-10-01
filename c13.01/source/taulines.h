/* This file is part of Cloudy and is copyright (C)1978-2013 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef TAULINES_H_
#define TAULINES_H_

#include "transition.h"

extern bool lgStatesAdded;
extern bool lgLinesAdded;
extern qList AnonStates;

extern char **chSpecies;
extern species *dBaseSpecies;
extern vector<qList > dBaseStates;
extern vector< multi_arr<int,2> > ipdBaseTrans;
extern vector<TransitionList> dBaseTrans;
extern multi_arr<CollRateCoeffArray,2> AtmolCollRateCoeff;
extern CollSplinesArray ****AtmolCollSplines;
extern StoutColls ****StoutCollData;
extern long int nSpecies;
/*************************/
void database_readin( void );
void dBase_solve(void );

/**
 * this is a dummy optical depth array for non-existant lines
 */
extern TransitionProxy::iterator TauDummy;

/**>>chng 99 sep 14, comments into level1.dat, count number
 * of lines then MALLOC the space
 * NB must MALLOC nlines + 1 since dummy line is first one */
extern TransitionList TauLines;

/** this is the set of extra lines,
 * ExtraLymanLines[ipISO][ipZ][n]*/
extern multi_arr<int,3> ipExtraLymanLines;
extern vector<vector<TransitionList> > ExtraLymanLines;

/** the set of inner shell lines */
extern long int nUTA;
extern TransitionList UTALines;

/** this is the number of level 1 lines, and is set in atmdat_readin
 * by counter number of data lines in level1.dat */
extern long int nLevel1;
/**extern transition TauLines[NTAULINES+1];*/

/** these are the public parts of the hyperfine structure line transfer info
 * data gathered from hyperfine.dat using routines in hyperfine.c
 * the structure containing the hfs line information */
/* abundances of these isotopes relative to main species are in hyperfine.h */
extern TransitionList HFLines;
/** the number of lines */
extern long int nHFLines;

/**
 * main line arrays for hydrogenic ions<br><br>
 *
 * first dimension is atomic number<br>
 * second dim is upper level<br>
 * third dim is lower level<br>
 * nta dim is set of pointers for quantities within line transfer arrays<br>
 * in the forc translation, the upper level was too low by 1, since the<br>
 * fortran was starting at 1.  the lower dim was not changed by translation<br>
 * since it started from ip1s = 0 <br>
 * any place where the third dim has -1 is probably a remnant from forc and is wrong<br>
 */

//extern vector<vector<multi_arr<int,2> > > ipTransitions;
extern vector<vector<TransitionList> > Transitions;

extern vector<TransitionList> AllTransitions;
extern void checkTransitionListOfLists(vector<TransitionList>&); 

extern multi_arr<int,2> ipFe2LevN;
extern TransitionList Fe2LevN;

/** lines forming from doubly excited states */
extern multi_arr<int,3> ipSatelliteLines; /* [ipISO][nelem][level] */
extern vector<vector<TransitionList> > SatelliteLines; /* [ipISO][nelem][level] */

/** this will be set true once space is allocated for the HydroLines array.
 * from then on any HYDROGENIC LEVELS command will be ignored, this is
 * set to false in cddefines.c */
extern bool lgHydroMalloc;

// number of direct excitation routes in [N I] fluorescence
const int NI_NDP = 9;

/** the following pointers to lines within the level 1 stack are defined in
 * atmdat_readin where they are initially set to a large value, then reset
 * to point to the correct line within the stack 
 *
 * NB NB NB - lines must be entered both here and in atmdat_readin.cpp 
 * valid values are given in atmdat_lines_setup.cpp */
extern long ipT1656 , ipT9830 , ipT8727 , ipT1335 , 
	ipT1909 ,ipT977 , ipT1550 , ipT1548 , ipT386 , ipT310 , ipc31175 , ipT291 , ipT280 ,
	ipT274 , ipT270 , ipT312 , ipT610 , ipT370 , ipT157 , ipT1085 , 
	ipT990 , ipT1486 , ipT765 , ipT1243 , ipT1239 , ipT374g , ipT374x , ipT1200 ,
	ipT2140 , ipT671 , ipT315 , ipT324 , ipT333 , ipT209 , ipT122 , ipT205 ,
	ipT57 , ipT6300 , ipT6363 , ipT5577 , ipT834 , ipT1661 , ipT1666 , ipT835 ,
	ipT789 , ipT630 , ipT1304 , ipSi10_606 , ipT1039 , ipT8446 , ipT4368 , ipTOI13 ,
	ipTOI11 , ipTOI29 , ipTOI46 , ipTO1025 , ipT304 , ipT1214 , ipT150 , ipT146 ,
	ipT63 , ipTO88 , ipT52 , ipT26 , ipT1032 , ipT1037 , ipT770 , ipT780 , ipxNe0676 , ipT895 , ipT88 , ipTNe13 ,
	ipTNe36 , ipTNe16 , ipTNe14 , ipTNe24 , ipT5895 , ipfsNa373 , ipfsNa490 , ipfsNa421 ,
	ipxNa6143 , ipxNa6862 , ipxNa0746 , ipMgI2853 , ipMgI2026 , ipT2796 , ipT2804 ,
	ipT705 , ipT4561 , ipxMg51325 , ipxMg52417 , ipxMg52855 , ipxMg71190 , ipxMg72261 ,
	ipxMg72569 , ipxMg08303 , ipTMg610 , ipTMg625 , ipT58 , ipTMg4 , ipTMg14 , ipTMg6 ,
	ipfsMg790 , ipfsMg755 , ipAlI3957 , ipAlI3090 , ipT1855 , ipT1863 , ipT2670 ,
	ipAl529 , ipAl6366 , ipAl6912 , ipAl8575 , ipAl8370 , ipAl09204 , ipT639 ,
	ipTAl550 , ipTAl568 , ipTAl48 , ipSii2518 , ipSii2215 , ipT1808 ,
	ipT1207 , ipT1895 , ipT1394 , ipT1403 , ipT1527 , ipT1305 , ipT1260 , ipSi619 ,
	ipSi10143 , ipTSi499 , ipTSi521 , ipTSi41 , ipTSi35 , ipTSi25 , ipTSi65 ,
	ipTSi3 , ipTSi4 , ipP0260 , ipP0233 , ipP0318 , ipP713 ,
	ipT1256 , ipT1194 , ipTS1720 , ipT1198 , ipT786 ,
	ipTS19 , ipTS34 , ipTS11 , ipfsCl214 , ipfsCl233 ,
	ipCl04203 , ipCl04117 , ipCl973 ,
	ipTAr7 , ipTAr9 , ipTAr22 , ipTAr13 , ipTAr8 , ipAr06453 ,
	ipKI7745 , ipxK03462 , ipxK04598 , ipxK04154 ,
	ipxK07319 , ipCaI4228 , ipT3934 , ipT3969 , ipT8498 , ipT8542 ,
	ipT8662 , ipT7291 , ipT7324 , ipTCa3 ,
	ipSc05231 , ipSc13264 , ipTi06172 , ipVa07130 , ipVa15172 ,
	ipFeI3884 , ipFeI3729 , ipFeI3457 , ipFeI3021 , ipFeI2966 , ipTuv3 ,
	ipTr48 , ipTFe16 , ipTFe26 , ipTFe34 , ipTFe35 , ipTFe46 , ipTFe56 , ipT1122 ,
	ipT191 , ipCo11527;
/** NB NB NB - lines must be entered both here and in atmdat_lines_setup.cpp where they
 * are actually defined and initialized!!  */
extern long ipS4_1405,ipS4_1398,ipS4_1424,ipS4_1417,ipS4_1407,
	ipO4_1400,ipO4_1397,ipO4_1407,ipO4_1405,ipO4_1401,
	ipN3_1749,ipN3_1747,ipN3_1754,ipN3_1752,ipN3_1751,
	ipC2_2325,ipC2_2324,ipC2_2329,ipC2_2328,ipC2_2327,
	ipSi2_2334,ipSi2_2329,ipSi2_2350,ipSi2_2344,ipSi2_2336,
	ipS1_25m ,ipS1_56m, ipCl1_11m , ipFe1_24m, ipFe1_35m , ipFe1_54m , ipFe1_111m,
	ipNi1_7m , ipNi1_11m , ipSi1_130m , ipSi1_68m , ipNI_pumpDirect[NI_NDP],
	ipNI_pumpIndirect;

/* NB NB NB - lines must be entered both here and in atmdat_readin where they
 * are actually defined and initialized!!  */

/* all of Dima's level 2 lines */

/**number of level 2 lines, dim for WindLine array */
const int NWINDDIM = 6744;

/** this is set to 0 with no atom_level2 command, normally
 * equal to NWINDDIM, definition is in cddefines.c */
extern long	nWindLine;

/* these are the level two lines themselves */
/** pointers to element and ion, TauLine2[line number][pointer within vector] */
extern TransitionList TauLine2;

extern realnum *cs1_flag_lev2;

#endif /* TAULINES_H_ */
