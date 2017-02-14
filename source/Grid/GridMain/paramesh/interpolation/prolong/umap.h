! Any source file including this header file should include "Flash.h" before it. - KW

#ifdef FIXEDBLOCKSIZE
       PARAMETER ( mvx_m   = MAXCELLS)
#else
! This file should not be included by any source files when not using fixed block sizes,
! since this file should be only included when using a PARAMESH Grid implementation! - KW
#endif

       PARAMETER ( mui_m   = NUNK_VARS)


        PARAMETER ( mvxu_m  = mvx_m * mui_m )
        PARAMETER ( mvxmu_m = mvx_m * mvxu_m )
