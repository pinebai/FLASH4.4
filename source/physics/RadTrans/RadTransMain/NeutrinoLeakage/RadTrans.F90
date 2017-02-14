!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/RadTrans
!!
!!  NAME 
!!
!!  RadTrans
!!
!!  SYNOPSIS
!!
!!  call RadTrans( integer(IN) :: nblk,
!!                 integer(IN) :: blklst(nblk),
!!                 real(IN)    :: dt, 
!!       optional, integer(IN) :: pass)
!!
!!  DESCRIPTION 
!!      This routine computes neutrino source terms using the Rosswog
!!      leakage approach.  Kernel implementations are from GR1D.
!!
!! ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!   dt     : The time step
!!   pass   : reverses solve direction
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine RadTrans(nblk, blklst, dt, pass)
  use RadTrans_data, ONLY : rt_useRadTrans
  use Driver_interface, ONLY : Driver_getSimTime, Driver_getNStep
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Eos_nucInterface, ONLY : Eos_nucDetectBounce
  use rt_data, ONLY : rt_rayData, rt_tauRuff, rt_chiRoss, &
       rt_heatFlx, rt_heatErms, rt_heatEave, &
       rt_leakNumRad, rt_leakNumRays, rt_nstepStart, rt_nstep, &
       rt_bounceTime, rt_reducedTime, rt_reducedSteps

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  real,    intent(in) :: dt
  integer, intent(in), optional :: pass

  real :: t_bounce, time
  logical, save :: postBounce = .false.
  logical :: doSample

  if (.not. rt_useRadTrans) return
  call Timers_start("RadTrans")

  if (.not. postBounce) call Eos_nucDetectBounce(postBounce,bounceTime=rt_bounceTime)
  ! Only compute leakage post-bounce:
  if (.not. postBounce) then
     call Timers_stop("RadTrans")
     return
  end if

  call Driver_getSimTime(time)
  call Driver_getNStep(rt_nstep)
  ! Decide whether or not to compute leakage on this time step
  doSample = .TRUE.
  if (time - rt_bounceTime > rt_reducedTime) then
     if (mod(rt_nstep-rt_nstepStart, rt_reducedSteps) /= 0) then
        doSample = .FALSE.
     end if
  end if


  ! First, loop over blocks and fill in local parts of the leakage arrays
  if (doSample) then
     call Timers_start("ray sampling")
     call rt_sampleRays(nblk,blklst)
     call Timers_stop("ray sampling")

     ! Now every rank has the entire leakage array, now compute taus ray-by-ray
     call Timers_start("ray-by-ray tau")
     call rt_calcTau
     call Timers_stop("ray-by-ray tau")
  end if

  ! Now we must inerpolate back to hydro grid and apply source terms
  ! This is done EVERY time step, regardless of the truth of doSample
  call Timers_start("local leakage")
  call rt_calcLeak(nblk,blklst,dt)
  call Timers_stop("local leakage")

  call Timers_stop("RadTrans")
          
end subroutine RadTrans

