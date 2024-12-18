/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief TPCI functions
  
  This is the main file of "The PLUTO CLOUDY Interface"
  (TPCI) and contains all the functions.
  
  Minor changes were applied to the source codes of PLUTO and CLOUDY.
  
  PLUTO:
  http://plutocode.ph.unito.it/
  
  CLOUDY:
  http://www.nublado.org/
  
  Copyright 2014 Michael Salz
  
  \author M. Salz (msalz@hs.uni-hamburg.de)
  \date   Jun 6, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "cddefines.h"
#include "cddrive.h"

#ifdef PARALLEL
  #include <mpi.h>
#endif

#define REALNUM_DEFINED YES
  extern "C" {
    #include "pluto.h"
  }
#undef REALNUM_DEFINED

#include "params.h"


/*! \name Inverse domain loop
    - goes inversely through the x1 domain 
    - includes one boundary point (x1_beg)
*/
/**@{ */
#define INV_IDOM_LOOP(i)  for ((i) = IEND; (i) >= IBEG-1; (i)--)
/**@} */

#define USE_CLOUDY YES
#define USE_ADVEC NO
#define CLOUDY_PRINT_FREQ  10
#define CLOUDY_CONVERGE NO

#define CHANGE_FAKTOR     0.1
#define FRAC_COOL_TIMESTEP  0.1

int counter = 0;

int CallCloudy(Data *d, Grid *grid, int Cl_ncalls, double x1_dom_len, int Pl_k, int Pl_j, int koff, int joff, int lg_last_step);
void CloudyInputScript(Data *d, Grid *grid, int Cl_ncalls, double x1_dom_len, int Pl_k, int Pl_j, int koff, int joff, int lg_last_step);
void CloudyGetResults( Grid *grid, int Pl_k, int Pl_j );
void MapCloudytoPLUTO( Grid *grid, double ***Pl_val, int Pl_k, int Pl_j, 
                       double *Cl_depth, double *Cl_val, int Cl_nzone );
void RadiativeHeating(Data *d);
void RadiativeTimestep(Data *d,  Time_Step *Dts, int lg_last_step);

int CloudyRadSolve(Data *d, Time_Step *Dts, Grid *grid, int restart, int lg_last_step)
/*!
 * This is the interface to Cloudy.
 * 
 * - write data if convergance mode
 * - check if call to Cloudy necessary
 *   -> initialize + start Cloudy
 *   -> retrieve heating/cooling + ionization
 * - apply heating/cooling
 *
 * \param  d      pointer to PLUTO Data structure;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 * \param  restart  integer : YES/NO
 *
 * \return An integer giving success / failure of the Cloudy run.
 *
 *********************************************************************** */
{
  //  return 0;
    
  //printf("counter %d\n", counter);
  counter++;
  if ( counter < 1000 ) {
    //return 0;
  }

  int i, j, k;
  int joff, koff;
  int Cl_success = 1;
  bool lg_solve_rad = false;
  static bool lg_first_call = true;
  static bool lg_second_call = true;
  static int Cl_ncalls;
  double fac;
  
  double x1_dom_len;
  x1_dom_len = 0.9995*(grid[IDIR].x[IEND] - grid[IDIR].x[IBEG-1])*g_unitLength;
  
  static double ***last_dn, ***last_pr;
  
  double ***mean_mol;
  mean_mol = GetUserVar("U_MEAN_MOL");
  
  
  if (last_dn == NULL){
    last_dn  = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    last_pr  = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
    DOM_LOOP(k,j,i){ // initialize the arrays
      last_dn[k][j][i]  = d->Vc[DN][k][j][i];
      last_pr[k][j][i]  = d->Vc[PR][k][j][i];
    }
    if( restart == YES ){
      // no good solution for restart:
      // set Cl_ncalls and comment the print1 and QUIT lines
      Cl_ncalls = 100;
//       print1 ("! You must set the restart number for Cloudy!\n");
//       QUIT_PLUTO(1);
    }
    else{
      Cl_ncalls = 1;
    }
  }
  
  /* ------------------------------------------
     mean mol. weight must be initialized,
     used to compute temperature in tlaw
     ------------------------------------------ */
  if ( lg_first_call ){
    KDOM_LOOP(k){
      JDOM_LOOP(j){
        INV_IDOM_LOOP(i){
          mean_mol[k][j][i] = mu;
        };
      }
    }
  }
  
  /* ------------------------------------------
       check that the x-dir is not decomposed
       and Cloudy can compute the complete
       domain
     ------------------------------------------ */
  
  #ifdef PARALLEL
   int  pardims[DIMENSIONS];
   AL_Get_parallel_dim(SZ, pardims);
   if ( pardims[IDIR] ){
     print1 ((char *)"\n! X direction may not be decomposed when Cloudy is used !\n");
     print1 ((char *)"! Please use the command line option -no-x1par          !\n\n");
     QUIT_PLUTO(1);
   }
  #endif
  
  /* ------------------------------------------
     for parallel computations the global
     j and k indices are needed to save the
     Cloudy output
     ------------------------------------------ */
  
  #ifdef PARALLEL 
    int lbeg[DIMENSIONS], lend[DIMENSIONS], lghosts[DIMENSIONS];
    AL_Get_bounds(SZ, lbeg, lend, lghosts, AL_C_INDEXES);
    #if   ( DIMENSIONS == 2 )
      joff = lbeg[JDIR]-JBEG;
      koff = 0;
    #elif ( DIMENSIONS == 3 )
      joff = lbeg[JDIR]-JBEG;
      koff = lbeg[KDIR]-KBEG;
    #endif
  #else
    joff = 0;
    koff = 0;
  #endif
  

  /* ------------------------------------------------------
      Check if Cloudy must be called. True if:
      a) convergence:
          not yet
      b) evolve:
         - if the 
               i) density
              ii) pressure
             iii) velocity (NOT YET)
           at some point in the domain is
           changed by more than CHANGE_FAKTOR
     ------------------------------------------------------ */
  #if ( CLOUDY_CONVERGE )
    
  #else
    double max_fac_dn = 0;
    double max_fac_pr = 0;
    DOM_LOOP(k,j,i){
      // density
      fac = abs(last_dn[k][j][i]-d->Vc[DN][k][j][i])/MAX(last_dn[k][j][i],d->Vc[DN][k][j][i]);
      if (fac > max_fac_dn)
	max_fac_dn = fac;
	    
      if ( fac >= CHANGE_FAKTOR ) {
        lg_solve_rad = true;
      }
      // pressure
      fac = abs(last_pr[k][j][i]-d->Vc[PR][k][j][i])/MAX(last_pr[k][j][i],d->Vc[PR][k][j][i]);
      if ( fac >= CHANGE_FAKTOR ) {
        lg_solve_rad = true;
      }
      
      if (fac > max_fac_pr)
        max_fac_pr = fac;
    };

    if (counter % 1000 != 0) {
      lg_solve_rad = false;                  
    }
    //printf("max_fac_dn=%e, max_fac_pr=%e\n", max_fac_dn, max_fac_pr);
  #endif
  
  #ifdef PARALLEL
   int lgSC1 = lg_solve_rad;
   int lgSC2 = 0;
   MPI_Barrier ( MPI_COMM_WORLD );
   MPI_Allreduce ( &lgSC1, &lgSC2, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD );
   lg_solve_rad = lgSC2;
  #endif
  
  if ( lg_solve_rad || lg_first_call || lg_last_step){
    
    /* ------------------------------------------------------
        Call Cloudy and return success = 0
        - this is a loop over 2nd and 3rd dimensions, since
          Cloudy computes along the 1st dimension
       ------------------------------------------------------ */

    print1 ("> Cloudy: Solving Irradiation - file #%d \n", Cl_ncalls);
    
    KDOM_LOOP(k){
      JDOM_LOOP(j){
        Cl_success = CallCloudy(d, grid, Cl_ncalls, x1_dom_len, k, j, koff, joff, lg_last_step);
        if ( Cl_success != 0 ) { 
          print1 ("\n! PROBLEM DISASTER in Cloudy -> Cannot continue\n\n");
          QUIT_PLUTO(1);
        }
      }
    }
    #ifdef PARALLEL
     MPI_Barrier (MPI_COMM_WORLD);  // all irradiation slices should be finished
    #endif
    
    /* ------------------------------------------------------
        save the state for check on change/convergence
       ------------------------------------------------------ */
    
    DOM_LOOP(k,j,i){
      last_dn[k][j][i]  = d->Vc[DN][k][j][i];
      last_pr[k][j][i]  = d->Vc[PR][k][j][i];
    };
    
    Cl_ncalls ++;
    if( !lg_first_call ){
      lg_second_call = false;
    }
    lg_first_call = false;
  }
  
  /* ------------------------------------------------------
      Apply the radiative heating/cooling
     ------------------------------------------------------ */
  
  RadiativeHeating(d);
  
  /* ------------------------------------------------------
      Check the timestep
     ------------------------------------------------------ */
  
  RadiativeTimestep(d, Dts, lg_last_step);
  
  return Cl_success;
}


void RadiativeHeating(Data *d)
/*!
 * Apply the radiatvie heating/cooling
 * 
 * The userdef variable U_RAD_HEAT contains the
 * net heating-cooling (erg cm-3 s-1) computed in Cloudy.
 * This is applied in every hydro step.
 *
 * \param  d  pointer to PLUTO Data structure;
 * 
 *********************************************************************** */
{
  int k, j, i;
  double unitErg;
  unitErg = g_unitDensity*pow(g_unitVelocity,3)/g_unitLength;
  
  double ***rad_heat;
  rad_heat = GetUserVar("U_RAD_HEAT");
  
  DOM_LOOP(k,j,i){
    d->Vc[PR][k][j][i] += rad_heat[k][j][i]*(g_gamma-1)*g_dt /unitErg;
  };
  
}

void RadiativeTimestep(Data *d,  Time_Step *Dts, int lg_last_step)
/*!
 * Control timestepping via Cooling/Heating rate
 * 
 * The userdef variable U_RAD_HEAT contains the
 * net heating-cooling (erg cm-3 s-1) computed in Cloudy.
 * Energy may not change by more than 10 percent (TBD)
 *
 * \param  d  pointer to PLUTO Data structure;
 * \param  Dts    pointer to time Step structure;
 * 
 *********************************************************************** */
{
  int k, j, i;
  double dtcool, coolheat;
  double unitErg;
  unitErg = g_unitDensity*pow(g_unitVelocity,3)/g_unitLength;
  
  double ***rad_heat;
  rad_heat = GetUserVar("U_RAD_HEAT");
  
  dtcool = 1.e38;
  DOM_LOOP(k,j,i){
    coolheat = fabs(rad_heat[k][j][i]/unitErg);
    if (coolheat > 0.0){
      dtcool = MIN(dtcool, FRAC_COOL_TIMESTEP*d->Vc[PR][k][j][i]/(g_gamma-1)/coolheat);
    }
  };
  
  Dts->dt_cool = dtcool;
  
  /* ------------------------------------------
      print warning if advection should be
      included
     ------------------------------------------ */
  #if ( !USE_ADVEC )
    if ( g_stepNumber%1000 == 0 || lg_last_step){
      int lg_advec_rec, lg_advec_mol;
      double ***hd_time, ***hrec_time, ***hmol_time;
      hd_time = GetUserVar("U_HD_TIME");
      hrec_time = GetUserVar("U_HREC_TIME");
      hmol_time = GetUserVar("U_HMOL_TIME");
    
      lg_advec_rec = 0;
      lg_advec_mol = 0;
      DOM_LOOP(k,j,i){
         if( hrec_time[k][j][i]>0.0 && hrec_time[k][j][i]>hd_time[k][j][i]){
           lg_advec_rec = 1;
         };
         if( hmol_time[k][j][i]>0.0 && hmol_time[k][j][i]>hd_time[k][j][i]){
           lg_advec_mol = 1;
         };
      };
      
      if( lg_advec_rec ){
        print1 ("! WARNING: H recombination timescale is longer than the HD timescale.\n");
        print1 ("!          Please consider using advection!\n");
      }
      if( lg_advec_mol ){
        print1 ("! WARNING: Molecular timescales could be longer than the HD timescale.\n");
        print1 ("!          Please consider using advection!\n");
      }
    }
  #endif
}


int CallCloudy(Data *d, Grid *grid, int Cl_ncalls, double x1_dom_len, int Pl_k, int Pl_j, int koff, int joff, int lg_last_step)
/*!
 * Initialize + start Cloudy
 * (developed from Cloudy template)
 * 
 * - creates the Cloudy input script
 * - executes Cloudy
 * - retrieves the results and saves them
 * - catches all possible errors
 *
 * \param  d      pointer to PLUTO Data structure;
 * \param  Dts    pointer to time Step structure;
 * \param  grid   pointer to grid structure.
 *
 * \return An integer giving success / failure of the Cloudy run.
 *
 *********************************************************************** */
{
  exit_type exit_status = ES_SUCCESS;

  DEBUG_ENTRY( "main()" );

  try {

    bool Cl_lgAbort;
    long Cl_nw , Cl_nc , Cl_nn , Cl_ns , Cl_nte , Cl_npe , Cl_nione, Cl_neden;
    int i;
    
    /* ------------------------------------------------------
        initialize Cloudy
       ------------------------------------------------------ */
    
    cdInit();
    
    /* ------------------------------------------------------
        generate the input script
       ------------------------------------------------------ */
    
    CloudyInputScript(d, grid, Cl_ncalls, x1_dom_len, Pl_k, Pl_j, koff, joff, lg_last_step);
    
    /* ------------------------------------------------------
        execute the input script from above
       ------------------------------------------------------ */
    printf("I'm here3\n");  
    if( cdDrive() )
    {
      exit_status = ES_FAILURE;
    }
    printf("I'm here3.2\n");  
    /* ------------------------------------------------------
        retrieve the error messages
       ------------------------------------------------------ */    
    
    cdNwcns( &Cl_lgAbort , &Cl_nw , &Cl_nc , &Cl_nn , &Cl_ns , &Cl_nte , &Cl_npe , &Cl_nione, &Cl_neden );
    printf("I'm here3.5\n");  
    /* ------------------------------------------------------
        if Cloudy did not abort, get the results
       ------------------------------------------------------ */ 
    
    if( !Cl_lgAbort )
    {
      printf("I'm here4\n");  
      CloudyGetResults( grid, Pl_k, Pl_j );
      printf("I'm here5\n");  
    }
    
    cdEXIT(exit_status);
  }
  // here we catch all the possible exceptions that the code can throw
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
#      ifdef SIGBUS
      else if( e.sig() == SIGBUS )
      {
        fprintf( ioQQQ, " DISASTER - A bus error occurred. Bailing out...\n" );
        exit_status = ES_BUS_ERROR;
      }
#      endif
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
    if( ioQQQ != NULL )
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
  }
  catch( std::exception& e )
  {
    fprintf( ioQQQ, " DISASTER - An unknown exception was caught, what() = %s. Bailing out...\n",
       e.what() );
    exit_status = ES_UNKNOWN_EXCEPTION;
  }
  // generic catch-all in case we forget any specific exception above... 
  // so this MUST be the last one.
  catch( ... )
  {
    fprintf( ioQQQ, " DISASTER - An unknown exception was caught. Bailing out...\n" );
    exit_status = ES_UNKNOWN_EXCEPTION;
  }

  cdPrepareExit(exit_status);

  return exit_status;
}



void CloudyInputScript(Data *d, Grid *grid, int Cl_ncalls, double x1_dom_len, int Pl_k, int Pl_j, int koff, int joff, int lg_last_step)
/*!
 * Create the input script
 *
 *********************************************************************** */
{
  long nleft;
  char chLine [50];
  int i;
  double x1, val, dxmin;
  double ***mean_mol, ***rad_heat;
  mean_mol = GetUserVar("U_MEAN_MOL");
  rad_heat = GetUserVar("U_RAD_HEAT");
  
  /* ****************** IRRADIATION SED ********************* */
  nleft = cdRead("CMB");
  nleft = cdRead( "cosmic rays background" );
  
  nleft = cdRead( "init \"spectra.ini\"" );
      
  /* ************* GEOMETRY AND DENSITY STRUCTURE ********** */
  nleft = cdRead("radius 2.3e11 linear");
  
  sprintf( chLine , "stop depth %10.4e linear", x1_dom_len);
  //printf("Limit %s\n", chLine);
  nleft = cdRead( chLine );
  
  /* **** PASS DENSITY STRUCTURE FROM PLUTO TO CLOUDY ***** */
  //nleft = cdRead( "print off hide" );
  nleft = cdRead( "dlaw table depth linear" );
  INV_IDOM_LOOP(i){
    x1  = (grid[IDIR].x[IEND] - grid[IDIR].x[i])*g_unitLength;
    val = d->Vc[DN][Pl_k][Pl_j][i]*g_unitDensity/(CONST_amu*mu) * hydrogen_frac;
    //printf("dlaw %e %e\n", x1, val);
    // HOW DO I GET THE HYDROGEN DENSITY AUTOMATICALLY ????
    // Solar: 1.427  ,  ISM: 1.426  ,  H + He: 1.408  ,  H: 1.008
    sprintf( chLine , "%11.5e %11.5e", x1, val);
    nleft = cdRead( chLine );
  };
  nleft = cdRead( "end of dlaw" );
  nleft = cdRead( "print on" );


  /* ** PASS TEMPERATURE STRUCTURE FROM PLUTO TO CLOUDY **** */
  //nleft = cdRead( "print off hide" );
  nleft = cdRead( "tlaw table depth linear" );
  INV_IDOM_LOOP(i){
    x1  = (grid[IDIR].x[IEND] - grid[IDIR].x[i])*g_unitLength;
    val = KELVIN *mean_mol[Pl_k][Pl_j][i] *d->Vc[PRS][Pl_k][Pl_j][i]/d->Vc[RHO][Pl_k][Pl_j][i];
    sprintf( chLine , "%11.5e %11.5e", x1, val);
    nleft = cdRead( chLine );
    //printf("tlaw %e %s\n", grid[IDIR].x[IEND], chLine);
  };
  nleft = cdRead( "end of tlaw" );
  nleft = cdRead( "print on" );

  
  /* *********** PASS THE VELOCITY STRUCTURE ************** */
  #if ( USE_ADVEC )
  //nleft = cdRead( "print off hide" );
    nleft = cdRead( "wind advection table depth linear" );
    INV_IDOM_LOOP(i){
      x1   = (grid[IDIR].x[IEND] - grid[IDIR].x[i])*g_unitLength;
      val  = (-1.0)*d->Vc[VX][Pl_k][Pl_j][i]*g_unitVelocity;
      //sprintf( chLine , "%11.5e %11.5e", x1, val);
      if( val < 0.0 ){
        sprintf( chLine , "%11.5e %11.5e", x1, val);
      }
      else{
        sprintf( chLine , "%11.5e %11.5e", x1, -1.e-10);
      }
      //printf("velocity: %s\n", chLine);
      nleft = cdRead( chLine );
    };
    nleft = cdRead( "end of velocity table" );
    nleft = cdRead( "print on" );
    nleft = cdRead( "iterate 150" );
    nleft = cdRead( "set dynamics advection length fraction 0.01" );
  #else
    nleft = cdRead( "iterate 2" );
  #endif


  /* ************ SPEED UP / ITERATE ********************** */ 
//   nleft = cdRead( "stop zone 1" );
//   nleft = cdRead( "atom h-like levels small" ); 
//   nleft = cdRead( "atom he-like levels small" ); 
//   nleft = cdRead( "no level2" );
     nleft = cdRead( "no molecules" );
//   nleft = cdRead( "no opacity reevaluation" );
//   nleft = cdRead( "no ionization reevaluation" );
//   nleft = cdRead( "no fine opacities" );
//   nleft = cdRead( "no line transfer" );
  
  /* *************** PHYSICAL STUFF *********************** */
  nleft = cdRead( "element limit off -2" );
  nleft = cdRead("metals 0 linear");
//   nleft = cdRead( "element limit off -2.0" );
  nleft = cdRead( "stop temperature linear 5 K" );
  nleft = cdRead( "turbulence 1 km/sec no pressure" );
  nleft = cdRead( "double optical depth" );   
  nleft = cdRead( "abundances GASS10 no grains" );
//   nleft = cdRead( "elements read" );
//   nleft = cdRead( "helium" );
//   nleft = cdRead( "carbon" );
//   nleft = cdRead( "nitrogen" );
//   nleft = cdRead( "oxygen" );
// //   nleft = cdRead( "neon" );
// //   nleft = cdRead( "magnesium" );
//   nleft = cdRead( "silicon" );
// //   nleft = cdRead( "sulphur" );
// //   nleft = cdRead( "argon" );
//   nleft = cdRead( "end  of elements" );
// 
//   nleft = cdRead( "element neon off " );
//   nleft = cdRead( "element magnesium off " );
//   nleft = cdRead( "element sulphur off " );
//   nleft = cdRead( "element argon off " );
//   nleft = cdRead( "element Lithium off " );
//   nleft = cdRead( "element Beryllium off" );
//   nleft = cdRead( "element Boron  off " );
//   nleft = cdRead( "element Fluorine  off " );
//   nleft = cdRead( "element Phosphor off" );
//   nleft = cdRead( "element Chlorine off" );
//   nleft = cdRead( "element Potassium off" );
//   nleft = cdRead( "element sodium off" );
//   nleft = cdRead( "element Aluminium off" );
//   nleft = cdRead( "element calcium off" );
//   nleft = cdRead( "element Scandium  off" );
//   nleft = cdRead( "element Titanium off" );
//   nleft = cdRead( "element Vanadium off" );
//   nleft = cdRead( "element Chromium off" );
//   nleft = cdRead( "element Manganese off" );
//   nleft = cdRead( "element Cobalt off" );
//   nleft = cdRead( "element Iron off" );
//   nleft = cdRead( "element Copper off" );
//   nleft = cdRead( "element Nickel off" );
//   nleft = cdRead( "element Zinc  off" );
  
  /* ******************* OUTPUT *************************** */
  
  nleft = cdRead( "print short" );
  nleft = cdRead( "print line faint -2 log" );
  printf("I'm here2\n");  
  cdTalk ( false );
  cdOutput( "cloudy.out", "a");
  if (Cl_ncalls == 1){cdOutput( "cloudy.out");}
  
  //if ( Cl_ncalls%CLOUDY_PRINT_FREQ == 0 || g_stepNumber == 0 || lg_last_step){ // if not print every Cloudy call
  if (true) {
    cdTalk ( true );
    #if ( DIMENSIONS == 1 )
      sprintf( chLine , "cl_data.%04d.out", Cl_ncalls);
      cdOutput( chLine);
      sprintf( chLine , "set save prefix \"cl_data.%04d.\"", Cl_ncalls);
      nleft = cdRead( chLine );
    #elif ( DIMENSIONS == 2 )
      sprintf( chLine , "cl_data.%04d.%02ld.out", Cl_ncalls, (Pl_j-JBEG+joff));
      cdOutput( chLine);
      sprintf( chLine , "set save prefix \"cl_data.%04d.%02ld.\"", Cl_ncalls, (Pl_j-JBEG+joff));
      nleft = cdRead( chLine );
    #else
      sprintf( chLine , "cl_data.%04d.%02ld.%02ld.out", Cl_ncalls, (Pl_j-JBEG+joff), (Pl_k-KBEG+koff));
      cdOutput( chLine);
      sprintf( chLine , "set save prefix \"cl_data.%04d.%02ld.%02ld.\"", Cl_ncalls, (Pl_j-JBEG+joff), (Pl_k-KBEG+koff));
      nleft = cdRead( chLine );
    #endif
    nleft = cdRead( "set save hash \"\"" );
//     nleft = cdRead( "save overview \"over.tab\" last" );
    nleft = cdRead( "save overview \"over.tab\" " );
    nleft = cdRead( "save pressure \"pres.tab\" last" );
//     nleft = cdRead( "save wind \"wind.tab\" last" );
    nleft = cdRead( "save wind \"wind.tab\" " );
    nleft = cdRead( "save dynamics advection \"dyna.tab\" last" );
    nleft = cdRead( "save continuum \"continuum.tab\" last units Angstrom" );
    // nleft = cdRead( "save hydrogen conditions \"H_cond.tab\" last" );
    nleft = cdRead( "save cooling \"cool.tab\" last" );
    // nleft = cdRead( "save heating \"heat.tab\" last" );
    nleft = cdRead( "save ages \"ages.tab\" last" );
    nleft = cdRead( "save species populations \"pops.tab\" last" );
    nleft = cdRead( "save species energies \"energies.tab\" last" );
    printf("I'm here3\n");  
  }
}

void CloudyGetResults( Grid *grid ,int Pl_k, int Pl_j )
/*!
 * Retrieve the results from the computation
 * 
 * - creates arrays to retieve the results from the last Cloudy
 *   computation
 * - calls the cdget-functions
 * - calls the maping function, which also saves the result in
 *   userdefined variables.
 *
 * \param [in] grid     pointer to grid structure.
 * \param [in] int      k-indice of current Cloudy run
 * \param [in] int      j-indice of current Cloudy run
 * 
 *********************************************************************** */
{
  int i;
  int Cl_nzone;
  static double *Cl_depth;
  static double *Cl_numden, *Cl_massden, *Cl_meanmol;    // do I need static???
  static double *Cl_cooling, *Cl_heating, *Cl_radheat;
  static double *Cl_radaccel, *Cl_heateff, *Cl_eden;
  double ***mean_mol, ***rad_heat, ***rad_accel, ***heat_eff, ***eden;
  double aux_heat;
  mean_mol = GetUserVar("U_MEAN_MOL");
  rad_heat = GetUserVar("U_RAD_HEAT");
  rad_accel = GetUserVar("U_RAD_ACCEL");
  heat_eff = GetUserVar("U_HEAT_EFF");
  eden = GetUserVar("U_EDEN");
  
  Cl_nzone = cdnZone();
 
  /* ------------------------------------------
      these vectors must be allocated each time
      because the number of zones in Cloudy
      changes
     ------------------------------------------ */
 
  Cl_depth     = ARRAY_1D(Cl_nzone, double);
  Cl_numden    = ARRAY_1D(Cl_nzone, double);
  Cl_massden   = ARRAY_1D(Cl_nzone, double);
  Cl_meanmol   = ARRAY_1D(Cl_nzone, double);
  Cl_cooling   = ARRAY_1D(Cl_nzone, double);
  Cl_heating   = ARRAY_1D(Cl_nzone, double);
  Cl_radheat   = ARRAY_1D(Cl_nzone, double);
  Cl_radaccel  = ARRAY_1D(Cl_nzone, double);
  Cl_heateff   = ARRAY_1D(Cl_nzone, double);
  Cl_eden      = ARRAY_1D(Cl_nzone, double);
  
  /* ------------------------------------------
      get results from last Cloudy run
     ------------------------------------------ */

  cdDepth_depth(Cl_depth);
  cdDenPart_depth(Cl_numden);
  cdDenMass_depth(Cl_massden);
  cdCooling_depth(Cl_cooling);
  cdHeating_depth(Cl_heating);
  cdRadAcce_depth(Cl_radaccel);
  cdEDEN_depth(Cl_eden);
  
  /* ------------------------------------------
      only pass difference of rad. heating/cooling
      if larger than 0.005*heating/cooling
      otherwise it is just noise
      - value from HeatCoolRelErrorAllowed in 
        Cloudy
     ------------------------------------------ */
  for ( i = 0; i < Cl_nzone; i++){
    Cl_meanmol[i] = Cl_massden[i]/(CONST_amu*Cl_numden[i]);
    aux_heat = fabs(Cl_heating[i] - Cl_cooling[i])/MAX(MAX(fabs(Cl_heating[i]),fabs(Cl_cooling[i])),1.e-30);
    Cl_radheat[i] = (aux_heat > 0.005 ? Cl_heating[i]-Cl_cooling[i]:0.0);
    Cl_heateff[i] = (aux_heat > 0.005 ? (Cl_heating[i] - Cl_cooling[i])/Cl_heating[i]:0.0);
    Cl_radaccel[i] *= -1;
  }
   
  /* ------------------------------------------
      linear interpolation onto userdef
      variables
      - first boundary point must also be
        assigned for mean mol. (needed to pass
        temp. from PLUTO to Cloudy
     ------------------------------------------ */
  
  MapCloudytoPLUTO( grid, mean_mol, Pl_k, Pl_j,
                    Cl_depth, Cl_meanmol, Cl_nzone );
  
  mean_mol[Pl_k][Pl_j][IBEG-1] = mean_mol[Pl_k][Pl_j][IBEG];
  
  MapCloudytoPLUTO( grid, rad_heat, Pl_k, Pl_j,
                    Cl_depth, Cl_radheat, Cl_nzone );
  
  MapCloudytoPLUTO( grid, rad_accel, Pl_k, Pl_j,
                    Cl_depth, Cl_radaccel, Cl_nzone );
  
  MapCloudytoPLUTO( grid, heat_eff, Pl_k, Pl_j,
                    Cl_depth, Cl_heateff, Cl_nzone );
  
  MapCloudytoPLUTO( grid, eden, Pl_k, Pl_j,
                    Cl_depth, Cl_eden, Cl_nzone );
  
  FreeArray1D(Cl_depth);
  FreeArray1D(Cl_numden);
  FreeArray1D(Cl_massden);
  FreeArray1D(Cl_meanmol);
  FreeArray1D(Cl_cooling);
  FreeArray1D(Cl_heating);
  FreeArray1D(Cl_radheat);
  FreeArray1D(Cl_radaccel);
  FreeArray1D(Cl_heateff);
  FreeArray1D(Cl_eden);
}

void MapCloudytoPLUTO( Grid *grid, double ***Pl_val, int Pl_k, int Pl_j, 
                       double *Cl_depth, double *Cl_val, int Cl_nzone )
/*!
 * Interpolate Cloudy results onto the PLUTO grid
 *
 * Interpolates one result (eg., radiative heating) onto a
 * userdef variable. The result is then available throughout the
 * code.   
 * 
 * \param [in] grid     pointer to grid structure.
 * \param [in] double   pointer userdef variable.
 * \param [in] int      k-indice of current Cloudy run
 * \param [in] int      j-indice of current Cloudy run
 * \param [in] double   pointer Cloudy depth structure.
 * \param [in] double   pointer Cloudy result structure.
 * \param [in] int      number of zone in Cloudy run
 *
 *********************************************************************** */
{
  int i,ilow, ihigh, imid;
  double x1, dlt_x;
  
  IDOM_LOOP(i){
    
    x1 = (grid[IDIR].x[IEND] - grid[IDIR].x[i])*g_unitLength;
    
    if (x1 > Cl_depth[Cl_nzone-1]){
      Pl_val[Pl_k][Pl_j][i] = Cl_val[Cl_nzone-1];
    }
    else if (x1 < Cl_depth[0]){
      Pl_val[Pl_k][Pl_j][i] = Cl_val[0];
    }
    else{
      /* *** TABLE LOOKUP *** */
      ilow = 0;
      ihigh = Cl_nzone - 1;
      while (ilow != (ihigh-1)){
        imid = (ilow+ihigh)/2;
        if (x1 <= Cl_depth[imid]){
          ihigh = imid;
        }else if (x1 > Cl_depth[imid]){
          ilow = imid;
        }
      }
      /* *** INTERPOLATE *** */
      dlt_x = (x1 - Cl_depth[ilow])/(Cl_depth[ihigh] - Cl_depth[ilow]);
      Pl_val[Pl_k][Pl_j][i] = Cl_val[ilow] + dlt_x*(Cl_val[ihigh] - Cl_val[ilow]);
    }
  } 
}







