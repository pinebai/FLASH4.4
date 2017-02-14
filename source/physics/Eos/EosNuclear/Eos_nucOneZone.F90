!!****if* source/physics/Eos/EosNuclear/Eos_nucOneZone
!!
!! NAME
!!
!!  Eos_nucOneZone
!!
!! SYNOPSIS
!!
!!  call Eos_nucOneZone(real(INOUT)  :: xdens,
!!                      real(INOUT)  :: xtemp,
!!                      real(IN)  :: xye,
!!                      real(INOUT)  :: xener,
!!                      real(INOUT)  :: xpres,
!!                      real(INOUT)  :: xentr,
!!                      real(OUT)  :: xVar,
!!                      integer(IN)  :: mode)
!!
!! DESCRIPTION
!!
!!  Routine for cranking the EOS on a single zone, in any old
!!  mode you like.  Calls kernel routine nuc_eos_full.  Can also
!!  be used to recover values from the table that are not provided
!!  by nuc_eos_short without doing a full EOS solve.
!!
!! ARGUMENTS
!!
!!   mode : Eos mode
!!
!! NOTES
!!      Parts of this subunit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in the kernel 
!!      directory are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      subunit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  
!!      stellarcollapse.org/equationofstate.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M. 2013, ApJ, 765, 29
!!
!!
!!***

subroutine Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdedt,xCs2,xXp,xXn,xXa,xXh,xVar,varID,mode)

#include "Flash.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_abortFlash
  use eosmodule, ONLY : precision, temp_mev_to_kelvin, e_zeroPoint, alltables

  implicit none 

  real, intent(INOUT) :: xDens
  real, intent(IN)    :: xYe
  real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
  integer, intent(IN) :: mode, varID
  real, intent(OUT) :: xXp, xXn, xXa,xXh,xdedt,xCs2,xVar

  interface
     subroutine nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar,xmu_e,xmu_n,xmu_p, &
          xmuhat,keytemp,keyerr,rfeps)
       implicit none
       real*8, intent(in)    :: xye
       real*8, intent(inout) :: xrho,xtemp,xenr,xent
       real*8, intent(out)   :: xprs,xcs2,xdedt
       real*8, intent(out)   :: xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp
       real*8, intent(out)   :: xabar,xzbar,xmu_e,xmu_n,xmu_p,xmuhat
       real*8, intent(in)    :: rfeps
       integer, intent(in)   :: keytemp
       integer, intent(out)  :: keyerr
     end subroutine nuc_eos_full
  end interface

  real :: xAbar,xZbar,xMu_e,xMu_n,xMu_p,xMuhat
  
  real, parameter :: KtoMev = 1./temp_mev_to_kelvin

  integer :: xMode, err
  real :: xdpderho, xGamc
  real :: xdpdrhoe,xGame

  real :: d1,d2,d3
  real :: lr, lt

  ! index var mapping:
  !  0 -> full eos call
  !  1 -> logpress
  !  2 -> logenergy
  !  3 -> entropy
  !  4 -> munu
  !  5 -> cs2
  !  6 -> dedT
  !  7 -> dpdrhoe
  !  8 -> dpderho
  !  9 -> muhat
  ! 10 -> mu_e
  ! 11 -> mu_p
  ! 12 -> mu_n
  ! 13 -> xa
  ! 14 -> xh
  ! 15 -> xn
  ! 16 -> xp
  ! 17 -> abar
  ! 18 -> zbar
  ! 19 -> gamma
  ! 20 -> mu_nu

  select case(mode)
  case(MODE_DENS_EI)
     xMode = 0
  case(MODE_DENS_TEMP)
     xMode = 1
  case(MODE_DENS_ENTR)
     xMode = 2
  case default
     call Driver_abortFlash('[Eos] Error: unsupported mode for Nuclear Eos')
  end select

  xTemp = xTemp * KtoMev
  xEner = xEner - e_zeroPoint

  if (varID == 0) then
     call nuc_eos_full(xDens,xTemp,xYe,xEner,xPres,xEntr,xCs2,xdedt,&
          xdpderho,xdpdrhoe,xXa,xXh,xXn,xXp,xAbar,xZbar,xMu_e,xMu_n,xMu_p, &
          xMuhat,xMode,err,precision)
  elseif (varID == 20) then
     ! this mode should only be used when grabbing variables for an already
     ! consistent thermodynamics.  Don't call this if thermo state has changed!
     lr = log10(xDens)
     lt = log10(xTemp)
     call findthis(lr,lt,xYe,xMu_e,alltables(:,:,:,10),d1,d2,d3)
     call findthis(lr,lt,xYe,xMu_p,alltables(:,:,:,11),d1,d2,d3)
     call findthis(lr,lt,xYe,xMu_n,alltables(:,:,:,12),d1,d2,d3)

     xVar = xMu_e - xMu_n + xMu_p

  else
     lr = log10(xDens)
     lt = log10(xTemp)
     call findthis(lr,lt,xYe,xVar,alltables(:,:,:,varID),d1,d2,d3)
  endif

!  if (err /= 0) then
!     call Driver_abortFlash('[Eos] Error in Eos_nucOneZone')
!  endif

  xTemp = xTemp * temp_mev_to_kelvin
  xEner = xEner + e_zeroPoint

  return
end subroutine Eos_nucOneZone
