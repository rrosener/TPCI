#define  PHYSICS                 MHD
#define  DIMENSIONS              2
#define  COMPONENTS              2
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     6

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  MHD_FORMULATION         CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD        NO
#define  RESISTIVE_MHD           NO
#define  THERMAL_CONDUCTION      EXPLICIT
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- pointers to user-def parameters -- */

#define  BMAG               0
#define  THETA              1
#define  T_IN               2
#define  T_OUT              3
#define  RHO_IN             4
#define  RHO_OUT            5

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         YES
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANALBADA_LIM
#define  CT_EMF_AVERAGE        UCT_HLL
#define  CT_EN_CORRECTION      NO
#define  ASSIGN_VECTOR_POTENTIAL  YES
#define  UPDATE_VECTOR_POTENTIAL  NO
