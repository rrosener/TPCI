//initial values

#define Rp (2 *  CONST_Rearth)
#define Mp (8.3 * CONST_Mearth)
#define mu 1.23
#define hydrogen_frac 0.92

#define init_temperature 1950
#define Ms (0.905 * CONST_Msun)
#define semimajor_axis (0.01544 * CONST_au)
#define n_surface 1e14
#define P_surface (n_surface * CONST_kB * init_temperature)
#define rho_surface (P_surface * (mu * CONST_mH) / (CONST_kB * init_temperature))
#define beta (CONST_G * Mp * (mu * CONST_mH) / (Rp * CONST_kB * init_temperature))
