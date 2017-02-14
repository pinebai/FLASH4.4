!!****f* source/physics/sourceTerms/Stir/Stir_computeDt
!!
!! NAME
!!
!!  Stir_computeDt
!!
!! SYNOPSIS
!!
!!  Stir_computeDt(integer(IN) :: blockID,
!!                 )
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(OUT)   :: dt_stir, 
!!                 real(OUT)   :: dt_minloc(5)) 
!!
!!
!! DESCRIPTION
!!
!!  compute a stiring timestep limiter
!!
!!   The timestep limiter would be:
!!
!!
!!             dt     =  min(cfl, sqrt(dx/a)
!!               stir                   
!!
!!             where "a" is the acceleration field 
!!
!!
!! ARGUMENTS
!!
!!  blockID       --  local block ID
!!  
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_stir       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!!
!!
!! SEE ALSO
!!
!!  Driver_computeDt
!!
!!***

subroutine Stir_computeDt(blockID,                 & 
                           blkLimits,blkLimitsGC,        &
                           solnData,                     &
                           dt_stir, dt_minloc)

! this version is the stub and should therefore not do anything.
! return value of dt_stir is set hugely large so it is not limiting.

#include "constants.h"
#include "Flash.h"


  implicit none


  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_stir
  integer, intent(INOUT)  :: dt_minloc(5)

  return

end subroutine Stir_computeDt


