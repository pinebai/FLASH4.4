!!****if* source/physics/sourceTerms/Flame/FlameEffects/BurnParametric/Flame_rhJumpReactive
!!
!! NAME
!!
!!  Flame_rhJumpReactive
!!
!! SYNOPSIS
!!
!!  call Flame_rhJumpReactive(real, dimension(EOS_NUM)(inout) :: eosdata_u,
!!                            real(in) :: qbar_u,
!!                            real, dimension(EOS_NUM)(out) :: eosdata_b,
!!                            real(out) :: qbar_b,
!!                            integer(in) :: eos_mode)
!!
!! DESCRIPTION
!!
!! calculate post-flame NSE state given unburned material with properties
!! eosData_u and qbar_u
!! eosData_u is updated once with Eos() called with eos_mode, allowing other
!! combinations of parameters other than density and temperature of the unburned
!! material to be specified
!!
!! Dean Townsley 2007,2008
!!
!! note that this routine assumes that eosData_u(EOS_ABAR) and eosData_u(EOS_ZBAR)
!! are already set appropriately.  i.e. this is a non-species-based routine
!!
!! ARGUMENTS
!!
!!   eosdata_u : 
!!
!!   qbar_u : 
!!
!!   eosdata_b : 
!!
!!   qbar_b : 
!!
!!   eos_mode : 
!!
!!
!!***


subroutine Flame_rhJumpReactive(eosData_u, qbar_u, eosData_b, qbar_b, eos_mode)

  use Eos_interface, ONLY : Eos
  use fl_effData, ONLY : fl_effsmlrho, fl_effEosTol

  use NSE_interface, ONLY: NSE_finalAtDens 
  
  implicit none

#include "constants.h"
#include "Eos.h"

  real, dimension(EOS_NUM), intent(inout) :: eosData_u
  real,    intent(in)                     :: qbar_u
  real, dimension(EOS_NUM), intent(out)   :: eosData_b
  real,    intent(out)                    :: qbar_b
  integer, intent(in)                     :: eos_mode

  integer, parameter :: max_newton = 50

  real, dimension(EOS_NUM)   :: eosData
  real   ::  ye, pres_u, hmq_u
  real   ::  dens_n, emq, pres_n, qbar, sumyi, tempguess, edot, yedot
  real   ::  error, dd, dpdd, f, dfdd, dens_n_old

  integer :: niters

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask 

  mask = .false.

  ! calculate thermodynimic info about initial state according to mode argument
  call Eos(eos_mode,1,eosData_u,mask=mask)
  pres_u = eosData_u(EOS_PRES)
  ye = eosData_u(EOS_ZBAR)/eosData_u(EOS_ABAR)
  hmq_u = eosData_u(EOS_EINT) + eosData_u(EOS_PRES)/eosData_u(EOS_DENS) - 9.6485e17*qbar_u

  ! A bit unusual here
  !  we will use the (dens,emq) table to construct the endpoint of
  ! a constant pressure burn.  Although this is what is stored in
  ! the (pres,hmq) table, due to interpolation accuracy and the fact
  ! that fully burned cells use the (dens,emq) table, we will instead use it
  ! to construct the fully burned state so that there is no "jolt" in
  ! fully burned cells at the beginning of the simulation.

  ! we solve for the density such that
  !  pres = pres_u
  ! with
  !  emq = hmq_u - pres_u / dens

  ! our initial guess
  dens_n = eosData_u(EOS_DENS)

  error = 2.0* fl_effEosTol
  niters = 0
  
  do while ( (error > fl_effEosTol) .and. (niters < max_newton) )

     ! get the nse state information
     emq = hmq_u - pres_u/dens_n
     call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n, emq)
     ! and pressure
     eosData(EOS_DENS) = dens_n
     eosData(EOS_TEMP) = tempguess
     eosData(EOS_EINT) = emq + 9.6485e17*qbar
     eosData(EOS_ABAR) = 1.e0/sumyi
     eosData(EOS_ZBAR) = ye * eosData(EOS_ABAR)
     call Eos(MODE_DENS_EI,1, eosData, mask=mask)

     pres_n = eosData(EOS_PRES)

     ! derivative
     dd = 1.0e-7*dens_n
     ! get the nse state information
     emq = hmq_u - pres_u/(dens_n+dd)
     call NSE_finalAtDens(qbar, sumyi, tempguess, edot, yedot, ye, dens_n+dd, emq)

     ! and pressure
     eosData(EOS_DENS) = dens_n + dd
     eosData(EOS_TEMP) = tempguess
     eosData(EOS_EINT) = emq + 9.6485e17*qbar
     eosData(EOS_ABAR) = 1.e0/sumyi
     eosData(EOS_ZBAR) = ye * eosData(EOS_ABAR)
     call Eos(MODE_DENS_EI,1, eosData, mask=mask)

     dpdd = (eosData(EOS_PRES)-pres_n)/dd

     ! function we are zeroing and deriv
     f = pres_n - pres_u
     dfdd = dpdd

     dens_n_old = dens_n
     dens_n = dens_n - f/dfdd

     if (dens_n .lt. fl_effsmlrho) then
        write (6,*) 'small density in nseJump'
        dens_n = 0.5*dens_n_old
     endif

     error = abs( (dens_n-dens_n_old)/dens_n )

     niters = niters + 1

  enddo
 
  if (niters >= max_newton) then
     write(6,*) 'exceeded number of newton steps in nsejump'
     write(6,*) 'ended with dens ', dens_n
  endif

  ! now fill output information
  emq = hmq_u - pres_u/dens_n
  call NSE_finalAtDens(qbar_b, sumyi, tempguess, edot, yedot, ye, dens_n, emq)

  ! and pressure
  eosData_b(EOS_DENS) = dens_n
  eosData_b(EOS_TEMP) = tempguess
  eosData_b(EOS_EINT) = emq + 9.6485e17*qbar_b
  eosData_b(EOS_ABAR) = 1.e0/sumyi
  eosData_b(EOS_ZBAR) = ye * eosData_b(EOS_ABAR)
  call Eos(MODE_DENS_EI,1, eosData_b, mask=mask)

  return
end subroutine Flame_rhJumpReactive

