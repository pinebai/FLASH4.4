!!****if* source/physics/Hydro/HydroMain/unsplit_rad/multiTemp/hy_uhd_MultiTempData
!!
!! NAME
!!   
!!  hy_uhd_MultiTempData 
!!
!!
!! SYNOPSIS
!!
!!  use hy_uhd_MultiTempData
!!
!!
!! DESCRIPTION
!!
!!***


module hy_uhd_MultiTempData

#include "Flash.h"
#include "constants.h"

  !!*****Runtime parameters*****
  integer, save :: hy_eosMode
!!$  integer, save :: hy_ppmEintFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEnerFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEintCFluxConstructionMeth = 0
!!$  integer, save :: hy_ppmEnerCFluxConstructionMeth = 0

  !!*****End Runtime parameters*****

  ! Based on Runtime parameters, these lower bounds may be used by
  ! Hydro as well as Eos implementations. - KW
  real, save :: hy_smallEion = 0.0, &
                hy_smallEele = 0.0, &
                hy_smallErad = 0.0

  !!*****End Directly Derived from Runtime parameters*****

  
  !!*****Constants database
  real,    save :: hy_eMass, hy_pMass, hy_eMassInUAmu


  ! For testing ways to advect components and handle shock heating
  
  integer, save :: hy_3Ttry_B, &
                   hy_3Ttry_E, &
                   hy_3Ttry_F, &
                   hy_3Ttry_G

  real, save :: hy_3Ttry_D
  integer, save :: hy_3Ttry_B_rad

  real, save :: hy_mtFactorB0, hy_mtFactorBrad0
  real, save :: hy_mtFactorB1, hy_mtFactorBrad1
  real, save :: hy_mtFactorB2, hy_mtFactorBrad2

  real, save :: hy_mtScaleWork, hy_mtScaleAccel
  real, save :: hy_mtScaleLorentz

  real,parameter :: hy_mtPradScaleFactor = 1.0 !hardwired for now
#ifdef FLLM_VAR
  integer,parameter :: hy_mtPradScaleVar = FLLM_VAR
#else
  integer,parameter :: hy_mtPradScaleVar = -1
#endif

  real, save :: hy_mtPresRatLambda3Min

  !! Variables for controlling radiation flux limiter smoothing
  logical, parameter :: hy_smoothFluxLim0 = .TRUE.
  integer, parameter :: hy_smoothFluxLim1 = 1 ! choose 0 or 1
  integer, save :: hy_smoothIterations
  integer, save :: hy_smoothMethod
  real,    save :: hy_smoothCoeff
  logical, save :: hy_useMinSmoothVarVal, hy_useMaxSmoothVarVal
  real,    save :: hy_minSmoothVarVal, hy_maxSmoothVarVal

end module hy_uhd_MultiTempData
