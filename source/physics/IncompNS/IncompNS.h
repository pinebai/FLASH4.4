!!***h* source/physics/IncompNS/IncompNSMain/IncompNS.h
!!
!! This is the internal header file for the 
!! Incompressible Navier Stokes Module
!!
!!***

#define DIR_X 1
#define DIR_Y 2
#define DIR_Z 3

#define INJECTION_PROL 0
#define DIVPRES_PROL   1
#define LINEAR_PROL    101
#define QUADRATIC_PROL 102

#define AB2_SCHM    2
#define AB2_SCHM_V  21
#define RK3_SCHM    3

#define INS_INTSCHM_MULTISTEP 100
#define INS_INTSCHM_RK        200

!!#define EXPONENTIAL_WBREF_RAMP 1
