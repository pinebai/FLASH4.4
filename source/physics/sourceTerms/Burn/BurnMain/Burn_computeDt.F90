!!****if* source/physics/sourceTerms/Burn/BurnMain/Burn_computeDt
!!
!! NAME
!!
!!  Burn_computeDt
!!
!! SYNOPSIS
!!
!!  Burn_computeDt(integer(IN) :: blockID,
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(INOUT) :: dt_burn, 
!!                 real(INOUT) :: dt_minloc(5)) 
!!
!!
!! DESCRIPTION
!!
!!  compute a burning timestep limiter, by trying to force the energy
!!  generation from burning to be smaller than the internal energy
!!  in a zone.
!!
!!   The timestep limiter would be:
!!
!!                                      eint
!!             dt     =  enucDtFactor * -----
!!               burn                   enuc
!!
!!  enuc is energy/volume/s, so the time factor is already in there, and we
!!  are actually doing
!!
!!             
!!                                      eint
!!             dt     =  enucDtFactor * -----    * dt
!!               burn                   enuc*dt
!!
!!  enuc*dt is the amount of energy / volume deposited in a zone by burning. 
!!  eint is the internal energy / volume in that zone.  If enuc*dt is 2x
!!  eint, then we want a timestep that is half the size.  
!!
!!  enucDtFactor is a prefactor to scaling the timestep.  In general, we aim
!!  for enuc*dt < enucDtFactor * eint.  For good coupling between the hydro
!!  and the burner, enucDtFactor should be < 1.
!!
!!
!! ARGUMENTS
!!
!!  blockID       --  local block ID
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_burn       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!!
!! PARAMETERS
!!
!!  enucDtFactor    A parameter, such that enuc*dt < enucDtFactor * eint,
!!                  that is, the energy release from burning divided by
!!                  the internal energy in that zone is < enucDtFactor.
!!
!! SEE ALSO
!!
!!  Driver_computeDt
!!
!! NOTE
!!  
!!  On some platforms the use of HUGE may cause problem, and may
!!  need to be replaced by a hardcoded number
!!
!!
!!***

!!REORDER(4): solnData

subroutine Burn_computeDt(blockID,  &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dt_burn, dt_minloc)

  use Burn_data, ONLY: bn_enucDtFactor, bn_useBurn, bn_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"

  !! arguments
  integer, intent(IN)   :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_burn
  integer, intent(INOUT)  :: dt_minloc(5)

  !! local variables
  real              :: dt_temp, dt_tempInv
  integer           :: temploc(5)
  integer           :: i, j, k

  real, PARAMETER :: SMALL = TINY(1.0)
  real :: eint_zone, energyRatioInv

!!===================================================================

  ! initialize the timestep from this block to some obscenely high number

  if (.not. bn_useBurn)  return

  dt_temp = HUGE(0.0)
  dt_tempInv = SMALL

  ! loop over all of the zones and compute the minimum eint/enuc
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef EINT_VAR
           ! compute the internal energy in the zone
           eint_zone = solnData(EINT_VAR,i,j,k) 
#else
           eint_zone = solnData(ENER_VAR,i,j,k) - &
                0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                solnData(VELY_VAR,i,j,k)**2 + &
                solnData(VELZ_VAR,i,j,k)**2)
#endif

           ! compute the ratio.  Note, it is the absolute value that matters.
           ! Also prevent a divide by zero by first computing and comparing
           ! the inverse of what we want, and then only (un)invert that inverse
           ! if it is a reasonable number.
           energyRatioInv = abs(solnData(ENUC_VAR,i,j,k)) / eint_zone

           if (energyRatioInv > dt_tempInv) then
              dt_tempInv = energyRatioInv
              dt_temp = 1.0 / energyRatioInv
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = bn_meshMe
           endif

        enddo
     enddo
  enddo


  ! Set the timestep from this block.
  ! A little bit of trickery to avoid multiplying HUGE by something that is > 1. - KW
  dt_temp = min( dt_temp, HUGE(0.0)/max(1.0,bn_enucDtFactor) )
  dt_temp = bn_enucDtFactor*dt_temp

  if (dt_temp < dt_burn) then
     dt_burn = dt_temp
     dt_minloc = temploc
  endif

  if(dt_burn <= 0.0) call Driver_abortFlash("[Burn]: computed dt is not positive! Aborting!")

  return

end subroutine Burn_computeDt
