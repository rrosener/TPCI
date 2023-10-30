//initial values

#define planet_radius (15.6 *  CONST_Rearth)
#define planet_mass (232 * CONST_Mearth)
#define mu 2.3
#define hydrogen_frac 0.91

#define init_temperature 1450
#define Ms (1.069 * CONST_Msun)
#define semimajor_axis (0.04634 * CONST_au)
#define n_surface 1e14
#define P_surface (n_surface * CONST_kB * init_temperature)
#define rho_surface (P_surface * (mu * CONST_mH) / (CONST_kB * init_temperature))
#define beta (CONST_G * planet_mass * (mu * CONST_mH) / (planet_radius * CONST_kB * init_temperature))
