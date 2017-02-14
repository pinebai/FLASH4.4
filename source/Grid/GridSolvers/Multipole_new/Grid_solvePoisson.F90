!!****if* source/Grid/GridSolvers/Multipole_new/Grid_solvePoisson
!!
!! NAME
!!
!!  Grid_solvePoisson
!!
!! SYNOPSIS
!!
!!  Grid_solvePoisson (integer(in)    :: iSoln,
!!                     integer(in)    :: iSrc, 
!!                     integer(6)(in) :: bcTypes,
!!                     real(2,6)(in)  :: bcValues,
!!                     real(inout)    :: poisFact)
!!
!! DESCRIPTION
!!
!!  Poisson solver routine.  This module implements the multipole
!!  summation method for isolated boundary problems.  Periodic
!!  problems are not supported; the multipole method should be
!!  used only for approximately spherical matter distributions.
!!
!! ARGUMENTS
!!
!!  iSoln          : index to variable containing potential
!!  iSrc           : index to variable containing density
!!  bcTypes        : boundary types along various faces,
!!                   only used in verifying that they are isolated
!!  bcValues       : the values to boundary conditions, currently not used
!!  poisFact       : factor to be used in calculation
!!
!!***

subroutine Grid_solvePoisson (iSoln,                   &
                              iSrc,                    &
                              bcTypes,                 &
                              bcValues,                &
                              poisFact)

  use Grid_interface,    ONLY : GRID_PDE_BND_ISOLATED,         &
                                Grid_getListOfBlocks

  use Driver_interface,  ONLY : Driver_abortFlash

  use Timers_interface,  ONLY : Timers_start,                  &
                                Timers_stop

  use gr_mpoleInterface, ONLY : gr_mpoleCenterOfExpansion,     &
                                gr_mpoleRadialSampling,        &
                                gr_mpoleAllocateRadialArrays,  &
                                gr_mpoleSetRadialBinData,      &
                                gr_mpolePrintRadialInfo,       &
                                gr_mpoleMoments,               &
                                gr_mpoleCollectMoments,        &
                                gr_mpolePotentials,            &
                                gr_mpoleDumpMoments,           &
                                gr_mpoleDeallocateRadialArrays

  use gr_mpoleData,      ONLY : gr_mpoleBlockCount,            &
                                gr_mpoleBlockList,             &
                                gr_mpoleMomentsDump,           &
                                gr_mpoleMultiThreading,        &
                                gr_mpoleRadialInfoPrint

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(in)    :: iSoln, iSrc
  integer, intent(in)    :: bcTypes (6)
  real,    intent(in)    :: bcValues (2,6)
  real,    intent(inout) :: poisFact
!  
!
!    ...Check the grid boundaries. Abort, if not isolated
!
!
  if (any (bcTypes (1:2*NDIM) /= GRID_PDE_BND_ISOLATED) ) then
      call Driver_abortFlash ("FATAL: Multipole Poisson solver requires isolated boundaries")
  end if
!  
!
!    ...Start timer.
!
!
  call Timers_start ("Multipole Solver")
!  
!
!    ...Preliminary chores for setting up the run.
!
!
  call Grid_getListOfBlocks         (LEAF, gr_mpoleBlockList, gr_mpoleBlockCount)
  call gr_mpoleCenterOfExpansion    (iSrc)
  call gr_mpoleRadialSampling       ()
  call gr_mpoleAllocateRadialArrays ()
  call gr_mpoleSetRadialBinData     ()
!
!
!     ...Print radial info if requested by user.
!
!
  if (gr_mpoleRadialInfoPrint) then
      call gr_mpolePrintRadialInfo ()
  end if

        call Timers_start             ("gr_mpoleMoments")
        call gr_mpoleMoments          (iSrc)
        call Timers_stop              ("gr_mpoleMoments")

        call Timers_start             ("gr_mpoleCollectMoments")   ! see 'Note' above
          call gr_mpoleCollectMoments ()
        call Timers_stop              ("gr_mpoleCollectMoments")   ! see 'Note' above

        call Timers_start             ("gr_mpolePotentials")
        call gr_mpolePotentials       (iSoln, poisFact)
        call Timers_stop              ("gr_mpolePotentials")

!
!
!    ...Dump the moments if requested by the user.
!
!
  if (gr_mpoleMomentsDump) then
      call gr_mpoleDumpMoments ()
  end if
!
!
!    ...Final chores.
!
!
  call gr_mpoleDeallocateRadialArrays ()
!  
!
!    ...End timer.
!
!
  call Timers_stop ("Multipole Solver")
!
!
!    ...Check the potentials obtained. This routine is only internal
!       and should only be activated by the author of the code for
!       debugging and accuracy check purposes.
!
!
!  call gr_mpolePotential_exact   (iSrc,iSoln, poisFact)
!  
!
!    ...Ready!
!
!
  return
end subroutine Grid_solvePoisson
