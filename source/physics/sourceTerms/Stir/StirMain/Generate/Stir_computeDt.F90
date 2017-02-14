!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/Stir_computeDt
!!
!! NAME
!!
!!  Stir_computeDt
!!
!! SYNOPSIS
!!
!!  Stir_computeDt(integer(IN) :: blockID,
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(OUT)   :: dt_stir, 
!!                 real(OUT)   :: dt_minloc(5)) 
!!
!!
!! DESCRIPTION
!!
!!  compute a stirring timestep limiter
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
! return value of dt_burn is set hugely large so it is not limiting.
  use Stir_data, ONLY : st_useStir,st_computeDt, st_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getDeltas
#include "constants.h"
#include "Flash.h"


  implicit none


  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_stir
  integer, intent(INOUT)  :: dt_minloc(5)
  real, dimension(MDIM)   :: delta
  real :: dt_temp
  integer :: i,j,k
!!===================================================================

  ! initialize the timestep from this block to some obscenely high number

  if ((.not. st_useStir).or.(.not.st_computeDt))  return

  dt_temp = HUGE(0.0)
  call Grid_getDeltas(blockID,delta)

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
           

#ifdef ACCX_VAR
           dt_temp = min(dt_temp,abs(delta(IAXIS)/solnData(ACCX_VAR,i,j,k)))
#endif

#ifdef ACCY_VAR
           dt_temp = min(dt_temp,abs(delta(JAXIS)/solnData(ACCY_VAR,i,j,k)))
#endif

#ifdef ACCZ_VAR
           dt_temp = min(dt_temp,abs(delta(KAXIS)/solnData(ACCZ_VAR,i,j,k)))
#endif
           
           if (dt_temp < dt_stir*dt_stir) then
              dt_stir = sqrt(dt_temp)
              dt_minloc(1) = i
              dt_minloc(2) = j
              dt_minloc(3) = k
              dt_minloc(4) = blockID
              dt_minloc(5) = st_meshMe
           end if
           
        enddo
     enddo
  enddo
  if(dt_stir <= 0.0) call Driver_abortFlash("[Stir]: computed dt is not positive! Aborting!")

  return
  
end subroutine Stir_computeDt


