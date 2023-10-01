#define  PHYSICS                 RMHD
#define  DIMENSIONS              1
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  INTERPOLATION           LINEAR
#define  TIME_STEPPING           HANCOCK
#define  DIMENSIONAL_SPLITTING   YES
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     16

/* -- physics dependent declarations -- */

#define    EOS                     IDEAL
#define    ENTROPY_SWITCH          NO
#define    MHD_FORMULATION         NONE

/* -- pointers to user-def parameters -- */

#define  RHO_LEFT           0
#define  VX_LEFT            1
#define  VY_LEFT            2
#define  VZ_LEFT            3
#define  BY_LEFT            4
#define  BZ_LEFT            5
#define  PR_LEFT            6
#define  RHO_RIGHT          7
#define  VX_RIGHT           8
#define  VY_RIGHT           9
#define  VZ_RIGHT           10
#define  BY_RIGHT           11
#define  BZ_RIGHT           12
#define  PR_RIGHT           13
#define  BX_CONST           14
#define  GAMMA_EOS          15

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING     NO
#define  WARNING_MESSAGES      YES
#define  PRINT_TO_FILE         NO
#define  INTERNAL_BOUNDARY     NO
#define  SHOCK_FLATTENING      NO
#define  ARTIFICIAL_VISCOSITY  NO
#define  CHAR_LIMITING         NO
#define  LIMITER               VANLEER_LIM
#define  ASSIGN_VECTOR_POTENTIAL  NO
#define  UPDATE_VECTOR_POTENTIAL  NO
#define  PRIMITIVE_HANCOCK  NO
