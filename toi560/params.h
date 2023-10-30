//initial values

#define planet_radius (2.9 *  CONST_Rearth)
#define planet_mass (11 * CONST_Mearth)
#define mu 1.23
#define hydrogen_frac 0.92

#define init_temperature 740
#define Ms (0.75 * CONST_Msun)
#define semimajor_axis (0.0596 * CONST_au)
#define n_surface 1e14
#define P_surface (n_surface * CONST_kB * init_temperature)
#define rho_surface (P_surface * (mu * CONST_mH) / (CONST_kB * init_temperature))
#define beta (CONST_G * planet_mass * (mu * CONST_mH) / (planet_radius * CONST_kB * init_temperature))
