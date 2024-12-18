#include "pluto.h"
#include "params.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *  CHANGES, M.Salz:
 *  - do not change the predifined user variables because TPCI
 *    uses them to save e.g. the cooling/heating distribution
 *  - you can add additional userdefined variables just as the
 *    PLUTO userguide describes
 *
 ***************************************************************** */
{
  int i, j, k;  
  double ***T, ***mean_mol, ***time;
  double ***eden, ***hd_time, ***hrec_time, ***hmol_time;
  double yearsec = 3600.*24.*365.25;
  
  T        = GetUserVar("U_TEMP");
  time     = GetUserVar("U_TIME");
  mean_mol = GetUserVar("U_MEAN_MOL");
  
  eden = GetUserVar("U_EDEN");
  hd_time = GetUserVar("U_HD_TIME");
  hrec_time = GetUserVar("U_HREC_TIME");
  hmol_time = GetUserVar("U_HMOL_TIME");
  
  if (g_stepNumber == 0){   
  /* ------------------------------------------
     init of mean mol. weight and rad. heating
     is not important here
     --> saves first PLUTO file correct
     ------------------------------------------ */
    double ***rad_heat;
    double ***rad_accel;
    double ***heat_eff;
    rad_heat = GetUserVar("U_RAD_HEAT");
    rad_accel = GetUserVar("U_RAD_ACCEL");
    heat_eff = GetUserVar("U_HEAT_EFF");
    DOM_LOOP(k,j,i){
      mean_mol[k][j][i] = mu;
      rad_heat[k][j][i] = 0.0;
      rad_accel[k][j][i] = 0.0;
      heat_eff[k][j][i] = 1.0;
      eden[k][j][i] = 1.0;
    }
  }
  
  DOM_LOOP(k,j,i){
    T[k][j][i] = KELVIN*mean_mol[k][j][i] *d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]; 
    hd_time[k][j][i] = grid[IDIR].dx[i]*g_unitLength/(( fabs(d->Vc[VX1][k][j][i]) < 1.e-10 ? 1.e-10:fabs(d->Vc[VX1][k][j][i]))*g_unitVelocity);
    // formulas are from the hazy userguide of Cloudy
    hrec_time[k][j][i] = 7.6*pow(T[k][j][i]*1.e-4,0.8)/(eden[k][j][i]*1.e-4)*yearsec;
    hmol_time[k][j][i] = 0.3*pow(T[k][j][i]*1.e-3,-0.8)/(eden[k][j][i]*1.e-9)*yearsec;
    time[k][j][i] = g_time;
  }
  
  
  
}
/* ************************************************************* */
void ChangeDumpVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;

}





