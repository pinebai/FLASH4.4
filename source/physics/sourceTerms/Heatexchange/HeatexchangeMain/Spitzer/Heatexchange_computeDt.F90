!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Spitzer/Heatexchange_computeDt
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
       hx_smallE, hx_relTol, hx_ieTimeCoef, hx_navo
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

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
  real :: eqtime

  ! EOS data to compute cvele:
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, dimension(EOS_NUM) :: eos_arr
  real :: massfrac(NSPECIES)
  real :: cvele, cvion
  real :: abar, zbar

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
           t12diff  = temp1 - temp2

           eint1   = solnData(EION_VAR,i,j,k)
           eint2   = solnData(EELE_VAR,i,j,k)

           ! Compute abar and zbar:
           call Eos_getAbarZbar(solnData(:,i,j,k),abar=abar,zbar=zbar)
           eos_arr(EOS_DENS) = solnData(DENS_VAR, i,j,k)
           eos_arr(EOS_TEMPION) = temp1
           eos_arr(EOS_TEMP)    = temp2
           eos_arr(EOS_TEMPELE) = temp2
           eos_arr(EOS_TEMPRAD) = solnData(TRAD_VAR,i,j,k)
           eos_arr(EOS_ABAR) = abar
           eos_arr(EOS_ZBAR) = zbar

           ! Get the electron specific heat:
           mask = .false.
           mask(EOS_CVELE) = .true.
           mask(EOS_CVION) = .true.
           mask(EOS_DET)   = .true.
           massfrac = solnData(SPECIES_BEGIN:SPECIES_BEGIN+NSPECIES-1,i,j,k)
           call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
           cvele = eos_arr(EOS_CVELE)
           cvion = eos_arr(EOS_CVION)

           call hx_ieEquilTime(&
                zbar, &
                abar, &
                temp2, temp1, &
                solnData(DENS_VAR, i,j,k)*hx_navo/abar, &
                eqtime)

           eqtime = eqtime * hx_ieTimeCoef

           q2dot = cvele / eqtime * t12diff
           q1dot = -q2dot

           ! compute the ratio.  Note, it is the absolute value that matters.
           ! Also prevent a divide by zero by first computing and comparing
           ! the inverse of what we want, and then only (un)invert that inverse
           ! if it is a reasonable number.
           if (temp1>0.0) then
              denom = max(abs(t12diff)*eint1/temp1,hx_relTol*eint1)
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
              denom = max(abs(t12diff)*eint2/temp2,hx_relTol*eint2)
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
           

!!$           if(i==blkLimits(LOW,IAXIS))print*,'elR:',energy1LogRate, energy2LogRate, energy3LogRate
           energyLogRate = max(energy1LogRate, energy2LogRate)

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
  dt_temp = min( dt_temp, HUGE(0.0)/max(1.0,(1.0+epsilon(1.0))*hx_dtFactor) )
  dt_temp = hx_dtFactor*dt_temp

  if (dt_temp < dt_heatXchg) then
     dt_heatXchg = dt_temp
     dt_minloc = temploc
  endif

  return

end subroutine Heatexchange_computeDt
