!!****if* source/physics/Eos/EosNuclear/eos_lowdens
!!
!! NAME
!!
!!  eos_lowdens
!!
!! SYNOPSIS
!!
!!  call eos_lowdens(real, INTENT(IN)  :: xdens,
!!                   real, INTENT(INOUT)  :: xtemp,
!!                   real, INTENT(IN)  :: xabar,
!!                   real, INTENT(IN)  :: xzbar,
!!                   real, INTENT(INOUT)  :: xener,
!!                   real, INTENT(OUT)  :: xpres,
!!                   real, INTENT(OUT)  :: xentr,
!!                   real, INTENT(OUT)  :: xgamc,
!!                   integer, INTENT(IN)  :: mode,
!!                   real, INTENT(IN)  :: precision)
!!
!! DESCRIPTION
!! eos for low density 
!!
!! ARGUMENTS
!!
!!   xdens : density
!!
!!   xtemp : temperature 
!!
!!   xabar : abar
!!
!!   xzbar : zbar
!!
!!   xener : energy
!!
!!   xpres : pressure
!!
!!   xentr : entropy
!!
!!   xgamc : gamc
!!
!!   mode : Eos mode
!!
!!   precision : precision
!!
!!
!!
!!***

subroutine eos_lowdens(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,mode,precision)

  use Driver_interface, only : Driver_abortFlash
  use Eos_data, ONLY : eos_smallt

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  ! Arguments
  integer, INTENT(IN) :: mode
  real, INTENT(IN) :: precision
  real, INTENT(IN) :: xDens, xAbar, xZbar
  real, INTENT(INOUT) :: xTemp, xEner
  real, INTENT(OUT) :: xEntr, xPres, xGamc

  integer, parameter :: maxNewton = 20
  integer :: k
  real :: ewant, temp0, err, tempNew
  real :: xdedt

  if (mode==MODE_DENS_TEMP) then

     call idealGas_rad(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,xdedt)

  elseif (mode==MODE_DENS_EI) then

     ewant = xEner
     temp0 = xTemp
     err = 0.0

     call idealGas_rad(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,xdedt)

     !check if we already have the right temp
     err = abs((xEner - ewant)/ewant)
     if (err <= precision) goto 70

     tempNew = xTemp - (xEner - ewant)/xdedt
     ! Don't allow the temperature to change by more than an order of magnitude 
     ! in a single iteration
     if (tempNew > 10.e0*temp0) tempNew = 10.e0*temp0
     if (tempNew < 0.1e0*temp0) tempNew = 0.1e0*temp0

     ! Check if we are freezing, if so set the temperature to smallt, and adjust 
     ! the error so we don't wait for this one
     if (tempNew .LT. eos_smallt) then
        tempNew = eos_smallt
        err     = 0.1*precision  ! wtf?
     endif

     xTemp = tempNew

     do k = 2, maxNewton
        if (err <= precision) goto 70

        call idealGas_rad(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,xdedt)

        tempNew = xTemp - (xEner - ewant)/xdedt
        ! Don't allow the temperature to change by more than an order of magnitude 
        ! in a single iteration
        if (tempNew > 10.e0*xTemp) tempNew = 10.e0*xTemp
        if (tempNew < 0.1e0*xTemp) tempNew = 0.1e0*xTemp

        ! Check if we are freezing, if so set the temperature to smallt, and adjust 
        ! the error so we don't wait for this one
        if (tempNew .LT. eos_smallt) then
           tempNew = eos_smallt
           err     = 0.1*precision  ! wtf?
        endif

        xTemp = tempNew

        err = abs((xEner - ewant)/ewant)

     enddo

     ! Too many iterations!
     print *, ' '
     print *, 'Newton-Raphson failed in routine Nuclear Eos'
     print *, '(e and rho as input):'
     print *, ' '
     print *, 'too many iterations'
     print *, ' '
     print *, ' temp = ', temp0, xTemp
     print *, ' dens = ', xDens
     print *, ' pres = ', xPres
     print *, ' ener = ', ewant, xEner
     print *, ' dedt = ', xdedt
     print *, '   Ye = ', xZbar / xAbar, xAbar
     print *, ' diff = ', (xEner-ewant)/xdedt
     print *, '  err = ', err


     call Driver_abortFlash('[Eos] Error: too many iterations in Nuclear Eos')

70   continue
     ! Crank it through one more time
     call idealGas_rad(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,xdedt)

  else
     call Driver_abortFlash('[Eos] Error: unknown input mode in routine Eos')
  endif

  return
end subroutine eos_lowdens


subroutine idealGas_rad(xDens,xTemp,xAbar,xZbar,xEner,xPres,xEntr,xGamc,xdedt)

  use eosmodule, ONLY : avo, kb_erg, sioncon, arad

  implicit none

  real, INTENT(IN)  :: xDens, xTemp, xAbar, xZbar
  real, INTENT(OUT) :: xEner, xPres, xEntr, xGamc, xdedt

  real :: deni, tempi, kavoy
  real :: ytot1, xni, x1, x2, y0, z0
  real :: pion, eion, sion
  real :: dpiondd, dpiondt, deiondd, deiondt
  real :: prad, erad, srad
  real :: dpraddd, dpraddt, deraddd, deraddt
  real :: xdpdrhoe, dpresdt
  real :: Tneg, radfac
  real :: presi, chit, chid, x7

  ytot1 = 1.0/xAbar
  kavoy = kb_erg*avo*ytot1
  deni = 1.0/xDens
  tempi = 1.0/xTemp

  !! ion section:
  xni     = avo * ytot1 * xDens

  pion    = xni * kb_erg * xTemp
  dpiondd = avo * ytot1 * kb_erg * xTemp
  dpiondt = xni * kb_erg

  eion    = 1.5e0 * pion * deni
  deiondd = (1.5e0 * dpiondd - eion) * deni
  deiondt = 1.5e0 * xni * kb_erg * deni
  
  !!  sackur-tetrode equation for the ion entropy of 
  !!  a single ideal gas characterized by abar
  x2      = xAbar*xAbar*sqrt(xAbar) * deni / avo
  y0      = sioncon * xTemp
  z0      = x2 * y0 * sqrt(y0)
  sion    = (pion/xDens + eion) * tempi + kavoy*log(z0)

  !!  radiation section:
  Tneg = (3.0 * kb_erg * xDens / (100. * xAbar * arad))**(1./3.)

  if (xDens >= 1.0e-9 .or. xTemp <= Tneg) then
     radfac = 1.0
  else
     radfac = exp((Tneg - xTemp)/Tneg)
  endif

  prad    = radfac * arad * xTemp * xTemp * xTemp * xTemp / 3.0
  dpraddt = 4.0e0 * prad * tempi
  dpraddd = 0.0e0

  x1      = prad * deni
  erad    = 3.0e0 * x1
  deraddd = -erad * deni
  deraddt = 4.0e0 * erad * tempi

  srad    = (x1 + erad)*tempi

  !! Now sum the components

  xPres    = pion + prad
  xEner    = eion + erad
  xEntr    = sion + srad
  ! Convert entropy to kB/baryon  
  xEntr = xEntr / (kb_erg*avo)

  xdpdrhoe = dpiondd + dpraddd
  xdedt    = deiondt + deraddt
  dpresdt  = dpiondt + dpraddt
 
  !!  form gamma_1
  presi = 1.0e0/xPres
  chit  = xTemp*presi * dpresdt
  chid  = xdpdrhoe * xDens * presi
  x7    = xPres * deni * chit/(xTemp * xdedt)
  xGamc  = chit*x7 + chid

  if (xdedt <= 2.) then
     print *, 'crap'
     stop
  endif

  return
end subroutine idealGas_rad
