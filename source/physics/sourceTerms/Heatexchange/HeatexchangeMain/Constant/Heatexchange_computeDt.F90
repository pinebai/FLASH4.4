!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Constant/Heatexchange_computeDt
!!
!! NAME
!!
!!  Heatexchange_computeDt
!!
!! SYNOPSIS
!!
!!  Heatexchange_computeDt(integer(IN) :: blockID,
!!                 )
!!                 integer(IN) :: blkLimits(2,MDIM)
!!                 integer(IN) :: blkLimitsGC(2,MDIM)
!!                 real,pointer::  solnData(:,:,:,:),   
!!                 real(OUT)   :: dt_heatXchg, 
!!                 real(OUT)   :: dt_minloc(5)) 
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
!!                                      eintN
!!             dt     =  hx_dtFactor * ------
!!               burn                   qNdot
!!
!!  qNdot is energy/mass/s, so the time factor is already in there, and we
!!  are doing
!!
!!             
!!                                      eintN
!!             dt     =  hx_dtFactor * ------   * dt
!!               burn                   qNdot*dt
!!
!!  qNdot*dt is the amount of energy / mass deposited in a zone by burning. 
!!  eintN is the internal energy / mass in that zone.  If qNdot*dt is 2x
!!  eintN, then we want a timestep that is half the size.  
!!
!!  hx_dtFactor is a prefactor to scaling the timestep.  In general, we aim
!!  for qNdot*dt < hx_dtFactor * eintN.  For good coupling between the hydro
!!  and the burner, hx_dtFactor should be < 1.
!!
!!
!! ARGUMENTS
!!
!!  blockID       --  local block ID
!!  
!!  blkLimits     --  the indices for the interior endpoints of the block
!!  blkLimitsGC   --  the indices for endpoints including the guardcells
!!  solnData      --  the physical, solution data from grid
!!  dt_heatXchg       --  variable to hold timestep constraint
!!  dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!                    (i,j,k), block ID, PE number
!!
!!
!! PARAMETERS
!!
!!  hx_dtFactor    A parameter, such that qNdot*dt < hx_dtFactor * eint,
!!                  that is, the energy release from burning divided by
!!                  the internal energy in that zone is < hx_dtFactor.
!!
!! SEE ALSO
!!
!!  Driver_computeDt
!!
!!
!!***

subroutine Heatexchange_computeDt(blockID,   &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dt_heatXchg, dt_minloc)

  use Heatexchange_data, ONLY:hx_meshMe, hx_dtFactor, hx_useHeatexchange, &
       hx_smallE, hx_relTol, &
       hx_c12, hx_c13, hx_c23
  
  implicit none

#include "constants.h"
#include "Flash.h"

  !! arguments
  integer, intent(IN)   :: blockID 
  integer, intent(IN),dimension(2,MDIM)::blkLimits,blkLimitsGC
  real, pointer           :: solnData(:,:,:,:) 
  real, intent(INOUT)     :: dt_heatXchg
  integer, intent(INOUT)  :: dt_minloc(5)

  !! local variables
  real              :: dt_temp, dt_tempInv
  real                       :: temp1,temp2,temp3, eint1,eint2,eint3
  real                       :: t12diff, t13diff, t23diff
  real                       :: q1dot, q2dot, q3dot
  real                       :: bsEint1,bsEint2,bsEint3
  real :: energyLogRate, energy1LogRate, energy2LogRate, energy3LogRate
  real              :: denom
  integer           :: temploc(5)
  integer           :: i, j, k, g

  real, PARAMETER :: SMALL = TINY(1.0)
  real, PARAMETER :: smallt = 1e4
  real :: eint_zone
!!===================================================================

  ! initialize the timestep from this block to some obscenely high number

  if (.not. hx_useHeatexchange)  return

  dt_temp = HUGE(0.0)
  dt_tempInv = SMALL

  ! loop over all of the zones and compute the minimum eintN/qNdot
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           temp1    = solnData(TION_VAR,i,j,k)
           temp2    = solnData(TELE_VAR,i,j,k)
           temp3    = solnData(TRAD_VAR,i,j,k)
           t12diff  = temp1 - temp2
           t13diff  = temp1 - temp3
           t23diff  = temp2 - temp3

#ifdef EINT_VAR
           ! compute the internal energy in the zone
           eint_zone = solnData(EINT_VAR,i,j,k) 
#else
           eint_zone = solnData(ENER_VAR,i,j,k) - &
                0.5*(solnData(VELX_VAR,i,j,k)**2 + &
                solnData(VELY_VAR,i,j,k)**2 + &
                solnData(VELZ_VAR,i,j,k)**2)
#endif


           eint1   = solnData(EION_VAR,i,j,k)
           eint2   = solnData(EELE_VAR,i,j,k)
           eint3   = solnData(ERAD_VAR,i,j,k)
!!$           eint1 = max(eint1,SMALL)
!!$           eint2 = max(eint2,SMALL)
!!$           eint3 = max(eint3,SMALL)

           bsEint1 = 1; bsEint2 = 1; bsEint3 = 1
           q3dot =   bsEint1 * hx_c13 * t13diff  +  bsEint2 * hx_c23 * t23diff
           q2dot =   bsEint1 * hx_c12 * t12diff  -  bsEint3 * hx_c23 * t23diff
           q1dot = -(bsEint2 * hx_c12 * t12diff  +  bsEint3 * hx_c13 * t13diff)

           ! compute the ratio.  Note, it is the absolute value that matters.
           ! Also prevent a divide by zero by first computing and comparing
           ! the inverse of what we want, and then only (un)invert that inverse
           ! if it is a reasonable number.
           if (temp1>0.0) then
              denom = max(min(abs(t12diff),abs(t13diff))*eint1/temp1,hx_relTol*eint1)
           else
              denom = eint1
           end if
           if (eint1 > 0.0 .AND. q1dot .LE. 0.0) then
              energy1LogRate = abs(q1dot) / denom
           else if (eint1 > 0.0 .AND. q1dot > 0.0) then
              energy1LogRate = abs(q1dot) / denom ! * 0.000001e-20
           else if (q1dot == 0) then
              energy1LogRate = 0.0
           else
              energy1LogRate = HUGE(0.0)
              print*,'#1:q1dot,eint1,energy1LogRate:',q1dot,eint1,energy1LogRate
           end if

           if (temp2>0.0) then
              denom = max(min(abs(t12diff),abs(t23diff))*eint2/temp2,hx_relTol*eint2)
           else
              denom = eint2
           end if
           if (eint2 > 0.0 .AND. q2dot .LE. 0.0) then
              energy2LogRate = abs(q2dot) / denom
           else if (eint2 > 0.0 .AND. q2dot > 0.0) then
              energy2LogRate = abs(q2dot) / denom ! * 0.000001e-20
           else if (q2dot == 0) then
              energy2LogRate = 0.0
           else
              energy2LogRate = HUGE(0.0)
              print*,'#2:q2dot,eint2,energy2LogRate:',q2dot,eint2,energy2LogRate
           end if

           if (temp3>0.0) then
              denom = max(min(abs(t13diff),abs(t23diff))*eint3/temp3,hx_relTol*eint3)
           else
              denom = eint3
           end if
!!$           print*,'#3:q3dot,eint3:',q3dot,eint3
           if (eint3 > 0.0 .AND. q3dot .LE. 0.0) then
              energy3LogRate = abs(q3dot) / denom
           else if (eint3 > 0.0 .AND. q3dot > 0.0) then
              energy3LogRate = abs(q3dot) / denom ! * 0.000001we-20
           else if (q3dot == 0) then
              energy3LogRate = 0.0
           else if (eint3 == 0.0 .AND. q3dot > 0.0) then
              energy3LogRate = abs(q3dot) / hx_smallE ! * 0.000001e-20
           else
              energy3LogRate = HUGE(0.0)
           end if
           

!!$           if(i==blkLimits(LOW,IAXIS))print*,'elR:',energy1LogRate, energy2LogRate, energy3LogRate
           energyLogRate = max(energy1LogRate, energy2LogRate, energy3LogRate)

           if (energyLogRate > dt_tempInv) then
999           format(12(1x,G14.3))
!!$              print 999,temp1,eint1,q1dot,energy1LogRate, temp2,eint2,q2dot,energy2LogRate, temp3,eint3,q3dot,energy3LogRate
              dt_tempInv = energyLogRate
              dt_temp = 1.0 / energyLogRate
              temploc(1) = i
              temploc(2) = j
              temploc(3) = k
              temploc(4) = blockID
              temploc(5) = hx_meshMe
           endif

        enddo
     enddo
  enddo


  ! Set the timestep from this block.
  ! A little bit of trickery to avoid multiplying HUGE by something that is > 1. - KW
  dt_temp = min( dt_temp, HUGE(0.0)/max(1.0,hx_dtFactor) )
  dt_temp = hx_dtFactor*dt_temp

  if (dt_temp < dt_heatXchg) then
     dt_heatXchg = dt_temp
     dt_minloc = temploc
  endif

  return

end subroutine Heatexchange_computeDt
