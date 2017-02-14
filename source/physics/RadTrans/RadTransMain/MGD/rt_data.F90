!!****if* source/physics/RadTrans/RadTransMain/MGD/rt_data
!!
!!  NAME 
!!    rt_data
!!
!!  SYNOPSIS
!!    use rt_data
!!
!!  DESCRIPTION 
!!    Stores data for MGD
!!
!!***

#include "Flash.h"

module rt_data
  implicit none
  
  logical, save :: rt_useMGD ! Flag to indiciate that MGD is in use

  ! The energy group boundaries for multigroup diffusion (ergs):
  real, save, allocatable :: rt_mgdBounds(:)

  ! Boundary conditions for multigroup diffusion calculation
  integer, save, allocatable :: rt_mgdDomainBC(:,:)

  ! The radiation energies to use in each group for Dirichlet
  ! boundaries
  real, save, allocatable :: rt_mgdBcVals(:,:)

  ! Flux limiter options:
  integer, save :: rt_mgdFlMode ! Indicates type of flux limiter
  real, save :: rt_mgdFlCoef ! The MGD flux limiter coefficient

  ! Actual number of energy groups in the simulation. This number must
  ! be less or equal to MGD_NGROUPS*meshCopyCount
  integer, save :: rt_mgdNumGroups

  ! Implicitness factor, 0.0 to 1.0
  real, save :: rt_mgdthetaImplct
  real, save :: rt_mgdthetaC=1.0  !additionaly for some solvers, implicitness for absorption coupling
  real, save :: rt_mgdthetaD=0.0  !additionaly for some solvers, implicitness for emission coupling

  logical, save :: rt_timeGroups
  logical, save :: rt_groupBarrier

  logical, save :: rt_tightIonCoupling=.FALSE.   !additionaly for some solvers - assume tight
                                                 !thermal coupling of electrons to ions?

  logical, save :: rt_computeDt

  real, save :: rt_precomputedDt
  real, save :: rt_precomputedMinLoc(5)
  real, save :: rt_tempChangeRelTol

  integer, parameter :: rt_gcMaskSize=NUNK_VARS
  logical,dimension(rt_gcMaskSize),save :: rt_gcMask

end module rt_data
