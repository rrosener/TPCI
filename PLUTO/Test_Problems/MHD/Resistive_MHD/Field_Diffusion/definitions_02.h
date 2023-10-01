#define  PHYSICS                 MHD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     3

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         DIV_CLEANING
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           SUPER_TIME_STEPPING
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  ETAX               0
#define  ETAY               1
#define  ETAZ               2

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     YES
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANLEER_LIM
#define  ASSIGN_VECTOR_POTENTIAL  NO
#define  UPDATE_VECTOR_POTENTIAL  NO
#define  PRIMITIVE_HANCOCK     YES
#define  STS_nu                0.01
