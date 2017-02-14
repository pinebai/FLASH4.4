!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/gr_hseStep
!!
!! NAME
!!
!!  gr_hseStep
!!
!! SYNOPSIS
!!
!!  call gr_hseStep(real(INOUT)  :: dens(:),
!!                  real(INOUT)  :: temp(:),
!!                  real(IN)  :: ye(:),
!!                  real(IN)  :: sumy(:),
!!                  integer(IN)  :: n,
!!                  real(IN)  :: inputg,
!!                  real(IN)  :: delta,
!!                  integer(IN)  :: direction,
!!                  integer(IN)  :: order,
!!                  integer(IN)  :: mode,
!!         OPTIONAL,real(IN)     :: massFrac(:))
!!
!! DESCRIPTION
!!
!!   This routine puts a cell indicated by n in hydrastatic equilibrium (HSE)
!!   with a neighboring cell (or cells if order=2) by adjusting its density
!!   and possibly temperature.
!!   
!!   HSE is given by
!!
!!      d P /d x = rho * g
!!
!!   where P is the pressure, rho is the density and g is the gravity.
!!   Gravity is assumed spatially uniform and given by the argument inputg.
!!   Here for order=2 this is quantized following
!!   Zingale etal. 2002, ApJS, 143, 539 by solving their eq (42)
!!
!!     P_{+1} - P_{0}  =  g*delta/12 * (5 rho_{+1} + 8 rho_{0} - rho_{-1})
!!
!!   where P_{+1} = P(rho_{+1}, T_{+1}).  This form is consistent with the
!!   second order reconstruction done by PPM, but is still imperfect since
!!   the resulting riemann problem is not checked to give zero flux.
!!   For order=1 (necessary for starting things out) we use
!!
!!     P_{+1} - P_{0}  =  g*delta * (rho_{+1} + rho_{0})/2
!!
!!
!!   The quantities being adjusted are determined by
!!   the "mode" argument, which can be one of the following
!!
!!    HSE_CONSTTEMP   The temperature of the reference neighbor cell is
!!                    copied to the cell being worked on and the density is
!!                    adjusted to obtain HSE.
!!    HSE_SETTEMP     The temperature in the array temp(:) for the cell being
!!                    worked on is assumed to have already been set and the
!!                    density is adjusted to obtain HSE.
!!    HSE_CONSTENTR   Both the temperature and density of the cell being
!!                    worked on are adjusted so that the entropy of this cell
!!                    matches that of the reference neigbor and it is in HSE.
!!                    (This is a common condition for convection zones.)
!!
!!   The direction argument determines which of the cell's neigbors are used
!!   as the reference for the HSE.
!!    HSE_FORWARD     cell n is placed in HSE with cell n-1 (and n-2)
!!    HSE_BACKWARD    cell n is placed in HSE with cell n+1 (and n+2)
!!
!!   The directional convention for gravity 
!!
!!        |--------|
!!        |        |            ^
!!        |  n+1   |  -         |
!!        |        |  |         |  positive inputg
!!        |--------|  delta
!!        |        |  |
!!        |   n    |  -
!!        |        |
!!        |--------|
!!        |        |
!!        |  n-1   |
!!        |        |
!!
!! NOTES
!!   WARNING!! some weirdness here...
!!   NOTE on composition handling (sumy, ye, and massFrac):
!!     The mode of composition handling is *not determined here*!!!
!!     If the EOS is setup to use species it will use massFrac;
!!     otherwise, if the Ye,SumY based EOS is used, they will be used.
!!
!!     If the sumy(:) and ye(:) arrays are being used, this composition
!!     information can be non-uniform.  If massFrac is being used, both the
!!     reference cell(s) and the cell being put in HSE will all have the same
!!     composition.
!!
!!   TODO clean up behavior described in NOTE
!!          i.e. nominally make it so that massFrac can also be non-uniform
!!
!! ARGUMENTS
!!
!!   dens : 1D array of cell densities
!!                       dens(n) is used as a guess value
!!                       dens(n) is set to HSE value on return
!!
!!   temp : 1D array of cell temperatures
!!                       temp(n) may be modified depending on mode
!!
!!   ye :   1D array of cell Ye values ( number of electrons per baryon)
!!                       see NOTE above
!!
!!   sumy : 1D array of cell sum Y_i (number of ions per baryon)
!!                       see NOTE above
!!
!!   n :         index in arrays of element to be put in HSE
!!
!!   inputg :    gravity (see above for sign convention)
!!
!!   delta :     distance between neighboring cell centers (assumed uniform)
!!
!!   direction : determines which neighbor (above or below) to use as
!!               reference.  See constants defined above
!!
!!   order :     order numerical difference (1 or 2)
!!
!!   mode :      determines thermodynamic constraint to use for HSE:
!!               one of HSE_CONSTTEMP or HSE_CONSTENTR for an isothermal or
!!               isentropic boundary, respectively, or HSE_SETTEMP to use
!!               the temperature already set.
!!
!!   massFrac :  species abundance set (applied to all cells) for
!!               when composition is species based
!!
!! HISTORY
!!    2007       HSE BC Implementation based on Zingale/Dursi/etc.2002 - Dean Townsley
!!    Dec-2009   Renamed to integrate into more general code           - Klaus Weide
!!
!!***


! routine to put one zone in HSE with its neighbor(s)
subroutine gr_hseStep(dens, temp, ye, sumy, n, inputg, delta, direction, order, mode, massFrac)

  use gr_bcHseData, ONLY : HSE_FORWARD, HSE_BACKWARD, HSE_CONSTENTR, HSE_CONSTTEMP, HSE_SETTEMP
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
!           determines whether we are deriving dens(n) from dens(n-2) and dens(n-1)
!           or from dens(n+1) and (n+2)
! mode :   HSE_CONSTENTR chooses an adiabatic gradient
!          HSE_CONSTTEMP chooses constant temperature from reference zone
!          HSE_SETTEMP   uses the temperature already set in the new zone
  real, optional, intent(IN)        :: massFrac(:)

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
  ! set pres0 from properties in reference zone
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
  call Eos(MODE_DENS_TEMP, 1, eosData,massFrac,mask=mask)
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
     call Eos(MODE_DENS_TEMP, 1, eosData,massFrac,mask=mask)
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

end subroutine gr_hseStep
