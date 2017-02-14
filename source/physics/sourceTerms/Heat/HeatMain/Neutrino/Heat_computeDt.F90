!!****if* source/physics/sourceTerms/Heat/HeatMain/Neutrino/Heat_computeDt
!!
!! NAME
!!  
!!  Heat_computeDt
!!
!!
!! SYNOPSIS
!! 
!!  Heat_computeDt ( integer(IN) : blockID, 
!!                   
!!                  real(IN):  x(:), 
!!                  real(IN): dx(:), 
!!                  real(IN): uxgrid(:),
!!                  real(IN): y(:), 
!!                  real(IN): dy(:), 
!!                  real(IN): uygrid(:), 
!!                  real(IN):  z(:), 
!!                  real(IN): dz(:), 
!!                  real(IN): uzgrid(:), 
!!                  real,pointer :  solnData(:,:,:,:),   
!!                  real,(INOUT):   dt_check, 
!!                  integer(INOUT): dt_minloc(:) )
!!  
!! DESCRIPTION
!!
!!  Computes the timestep limiter for heating source term solver.
!! 
!!
!!
!! ARGUMENTS
!!
!!  blockID        local block ID
!!  
!!  x, y, z         three, directional coordinates
!!  d*              deltas in each {*=x, y z} directions
!!  u*grid          velocity of grid expansion in {*=x, y z} directions
!!  solnData        the physical, solution data from grid
!!  dt_check        variable to hold timestep constraint
!!  dt_minloc(5)    array to hold limiting zone info:  zone indices
!!
!!***

subroutine Heat_computeDt (blockID, &
                              x, dx, uxgrid, &
                              y, dy, uygrid, &
                              z, dz, uzgrid, &
                              blkLimits,blkLimitsGC,        &
                              solnData,   &
                              dt_check, dt_minloc )
#include "Flash.h"
#include "constants.h"

  use Heat_data, ONLY : useHeat, ht_meshMe, ht_heatTimeFac
                        
  implicit none

  integer, intent(IN) :: blockID
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dt_check
  integer,INTENT(INOUT)    :: dt_minloc(5)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif

 !! local variables
  real              :: dt_temp, dt_tempInv
  integer           :: temploc(5)
  integer           :: i, j, k

  real, PARAMETER :: SMALL = TINY(1.0)
  real :: eint_zone, energyRatioInv


  ! initialize the timestep from this block to some obscenely high number

  if (.not. useHeat)  return

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
           energyRatioInv = abs(solnData(DELE_VAR,i,j,k)) / eint_zone

           if (energyRatioInv > dt_tempInv) then
              dt_tempInv = energyRatioInv
              dt_temp = 1.0 / energyRatioInv
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = ht_meshMe
           endif

        enddo
     enddo
  enddo


  ! Set the timestep from this block.
  ! A little bit of trickery to avoid multiplying HUGE by something that is > 1. - KW
  dt_temp = min( dt_temp, HUGE(0.0)/max(1.0,ht_heatTimeFac) )
  dt_temp = ht_heatTimeFac*dt_temp

  if (dt_temp < dt_check) then
     dt_check = dt_temp
     dt_minloc = temploc
  endif

  if(dt_check <= 0.0) then
     print *, 'dt=', dt_check
     call Driver_abortFlash("[Heat]: computed dt is not positive! Aborting!")
  endif

  return
end subroutine Heat_computeDt

