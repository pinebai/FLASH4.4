!!****if* source/Simulation/SimulationMain/RTFlame/sim_hse_step
!!
!! NAME
!!
!!  sim_hse_step
!!
!! SYNOPSIS
!!
!!  call sim_hse_step(real, dimension(:), intent(INOUT)  :: dens,
!!                    real, dimension(:), intent(INOUT)  :: temp,
!!                    real, dimension(:), intent(IN)  :: ye,
!!                    real, dimension(:), intent(IN)  :: sumy,
!!                    integer, intent(IN)  :: n,
!!                    real, intent(IN)  :: inputg,
!!                    real, intent(IN)  :: delta,
!!                    integer, intent(IN)  :: direction,
!!                    integer, intent(IN)  :: order,
!!                    integer, intent(IN)  :: mode)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   dens : 
!!
!!   temp : 
!!
!!   ye : 
!!
!!   sumy : 
!!
!!   n : 
!!
!!   inputg : 
!!
!!   delta : 
!!
!!   direction : 
!!
!!   order : 
!!
!!   mode : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! Dean Townsley 2007

! routine to put one zone in HSE with its neighbor(s)
subroutine sim_hse_step(dens, temp, ye, sumy, n, inputg, delta, direction, order, mode)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_interface, ONLY: Eos

  implicit none
#include "Eos.h"
#include "constants.h"


  real, dimension(:), intent(IN)    :: ye, sumy
  real, dimension(:), intent(INOUT) :: dens, temp
  integer, intent(IN)               :: n               ! index to update
  real, intent(IN)                  :: inputg, delta
  integer, intent(IN)               :: direction, order, mode
! direction HSE_FORWARD or HSE_BACKWARD
!           detirmines whether we are deriving dens(n) from dens(n-2) and dens(n-1)
!           or from dens(n+1) and (n+2)
! mode :   HSE_CONSTENTR chooses an adiabatic gradient
!          HSE_CONSTTEMP chooses constant temperature from reference zone
!          HSE_SETTEMP   uses the temperature already set in the new zone

  real    :: densm1, dens0, pres0, temp0
  real    :: densp1, tempp1, sumyp1, yep1, presp1
  real    :: localg

  real, dimension(EOS_NUM) :: eosData
  real    :: error
  integer :: iter, fromn
  integer :: max_iter = 20

  real    :: hse_tol = 1e-6
  real    :: f, dfdd, newdensp1, dp_hse

  logical :: mask(EOS_VARS+1:EOS_NUM)

  real    :: dtdp_ad0

!==========================================================

  ! first select inputs and sign of gravity based on direction chosen
  if (direction == HSE_FORWARD) then
     if (order == 2) densm1 = dens(n-2)
     dens0  = dens(n-1)
     temp0  = temp(n-1)
     localg    = inputg
     fromn = n-1
  else
     if (order == 2) densm1 = dens(n+2)
     dens0  = dens(n+1)
     temp0  = temp(n+1)
     localg    = -inputg
     fromn = n+1
  endif
  densp1 = dens(n)
  if (mode==HSE_SETTEMP) tempp1 = temp(n)
  if (mode==HSE_CONSTTEMP) tempp1 = temp(fromn)
  sumyp1 = sumy(n)
  yep1   = ye(n)
  ! set pres0 from prperties in reference zone
  ! also set up the adiabatic derivative from this zone if we need it
  eosData(EOS_DENS) = dens(fromn)
  eosData(EOS_TEMP) = temp(fromn)
  eosData(EOS_ABAR) = 1.0/sumy(fromn)
  eosData(EOS_ZBAR) = ye(fromn)*eosData(EOS_ABAR)
  mask(:) = .false.
  if (mode == HSE_CONSTENTR) then
     mask(EOS_DPT) = .true.
     mask(EOS_DET) = .true.
  endif
  call Eos(MODE_DENS_TEMP, 1, eosData, mask=mask)
  pres0 = eosData(EOS_PRES)
  if (mode == HSE_CONSTENTR) then
     dtdp_ad0 = eosData(EOS_TEMP)/pres0 * eosData(EOS_DPT)/dens0/eosData(EOS_DET)/eosData(EOS_GAMC)
  endif

  ! we are solving eq (42) in Zingale etal 2002 for densp1 = rho_{+1}
  !    P_{+1} - P_{0}  =  g*delta/12 * (5 rho_{+1} + 8 rho_{0} - rho_{-1})
  !  were P_{+1} = P(rho_{+1}, T)

  ! initial things that are constant during iteration
  eosData(EOS_ABAR) = 1.0/sumyp1
  eosData(EOS_ZBAR) = yep1*eosData(EOS_ABAR)
  ! initialize things that change
  ! densp1 was initialized above from input
  error = 2*hse_tol
  iter = 0

  mask(:) = .false.
  mask(EOS_DPD) = .true.
  do while (error > hse_tol .and. iter < max_iter)
     if (order == 1) then
        dp_hse = localg*delta*0.5*(densp1+dens0)
     else if (order == 2) then
        dp_hse = localg*delta/12.0*(5*densp1+8*dens0-densm1)
     endif
     if (mode==HSE_CONSTENTR) then
        tempp1 = temp0 + dtdp_ad0 * dp_hse
     endif
     eosData(EOS_DENS) = densp1
     eosData(EOS_TEMP) = tempp1
     call Eos(MODE_DENS_TEMP, 1, eosData,mask=mask)
     presp1 = eosData(EOS_PRES)
     if (order == 1) then
        f = presp1 - pres0 - dp_hse
        dfdd = eosData(EOS_DPD) - localg*delta*0.5
     else if (order == 2) then
        f = presp1 - pres0 - dp_hse
        dfdd = eosData(EOS_DPD) - localg*delta/12.0*5
     endif
     newdensp1 = densp1 - f/dfdd
     if (newdensp1 < 0.1*densp1) newdensp1 = 0.1*densp1
     if (newdensp1 > 10*densp1) newdensp1 = 10*densp1
     error = abs(densp1-newdensp1)*2.0/(densp1+newdensp1)
     densp1 = newdensp1
     iter = iter+1
  enddo

  ! handle non-convergence
  if (iter >= max_iter) then
     write (6,*) 'HSE did not converge, dens, temp = ', densp1, tempp1
     call Driver_abortFlash("HSE did not converge")
  endif

  ! last update to temperature if adiabatic gradient is seleceted
  if (mode==HSE_CONSTENTR) then
     if (order == 1) then
        dp_hse = localg*delta*0.5*(densp1+dens0)
     else if (order == 2) then
        dp_hse = localg*delta/12.0*(5*densp1+8*dens0-densm1)
     endif
     tempp1 = temp0 + dtdp_ad0 * dp_hse
  endif

  dens(n) = densp1
  temp(n) = tempp1

  return

end subroutine
