!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreInit
!!
!! NAME
!!
!!  gr_hypreInit
!!
!!
!! SYNOPSIS
!!
!!  call gr_hypreInit()
!!
!! Description
!!
!!  Initializes local data for Unit HYPRE defined in Module gr_hypreData.
!!  All the variables here are initialized by calling the
!!  RuntimeParameters_get subroutine. These data variables are for
!!  Unit Scope ->  HYPRE.
!!
!! ARGUMENTS
!!
!!  none  
!!***

subroutine gr_hypreInit()
  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                          RuntimeParameters_mapStrToInt
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Logfile_interface,           ONLY : Logfile_stamp
  use gr_hypreData, ONLY : gr_asol, gr_speedlt, gr_hypreSolverType,            &
                           gr_hypreNStep,                                      & 
                           gr_hyprePcType, gr_hypreRelTol, gr_hypreMaxIter,    & 
                           gr_hypreGridIsSetUp, gr_hyprePrintSolveInfo,        &
                           gr_hypreInfoLevel, gr_hypreUseFloor, gr_hypreFloor, &
                           gr_hypreUse2Norm, gr_hypreRecomputeResidualP,       &
                           gr_hypreAbsTol, gr_hypreRecomputeResidual,          &
                           gr_hypreCfTol, gr_hypreRelChange,                   &
                           gr_hypreSolverAutoAbsTolFact, gr_hypreSolverAbsTolEff,&
                           gr_hypreAnisoDiffusion
  
  implicit none 

#include "constants.h"
#include "gr_hypre_fortran_message.h"
 
  character(len=MAX_STRING_LENGTH) :: solver
  character(len=MAX_STRING_LENGTH) :: pc
  character(len=22) :: msgBuf(3,4)
  character(len=HYPRE_MSG_LEN) :: hypreMsg
  integer :: i,j
  
  real :: ssol
  
  call PhysicalConstants_get("speed of light",gr_speedlt)
  call PhysicalConstants_get("Stefan-Boltzmann",ssol)
  
  gr_asol    = 4.0 * ssol / gr_speedlt
  
  gr_hypreSolverType   = HYPRE_PCG
  gr_hyprePcType       = HYPRE_AMG
  gr_hypreRelTol       = 1.0E-12
  gr_hypreMaxIter      = 10000
  
  call RuntimeParameters_get("gr_hypreSolverType", solver)
  call RuntimeParameters_get("gr_hyprePcType", pc)
  
  call RuntimeParameters_mapStrToInt(solver,gr_hypreSolverType)
  call RuntimeParameters_mapStrToInt(pc,gr_hyprePcType)
  
  call RuntimeParameters_get("gr_hypreRelTol",  gr_hypreRelTol)
  call RuntimeParameters_get("gr_hypreAbsTol",  gr_hypreAbsTol)
  gr_hypreSolverAbsTolEff = gr_hypreAbsTol !may be changed by solver implementation
  call RuntimeParameters_get("gr_hypreSolverAutoAbsTolFact", gr_hypreSolverAutoAbsTolFact)
  if (gr_hypreSolverAutoAbsTolFact < 0.0) then
     gr_hypreSolverAutoAbsTolFact = abs(gr_hypreSolverAutoAbsTolFact) * gr_hypreRelTol
  end if
  call RuntimeParameters_get("gr_hypreCfTol",   gr_hypreCfTol)
  call RuntimeParameters_get("gr_hypreMaxIter", gr_hypreMaxIter)

  call RuntimeParameters_get("gr_hyprePrintSolveInfo", gr_hyprePrintSolveInfo)
  
  call RuntimeParameters_get("gr_hypreInfoLevel", gr_hypreInfoLevel)
  
  call RuntimeParameters_get("gr_hypreUseFloor", gr_hypreUseFloor)
  call RuntimeParameters_get("gr_hypreFloor",    gr_hypreFloor)

  call RuntimeParameters_get("gr_hypreUse2Norm",    gr_hypreUse2Norm)
  call RuntimeParameters_get("gr_hypreRecomputeResidual",  gr_hypreRecomputeResidual)
  call RuntimeParameters_get("gr_hypreRecomputeResidualP", gr_hypreRecomputeResidualP)
  call RuntimeParameters_get("gr_hypreRelChange",   gr_hypreRelChange)
  
  ! Write version info about the HYPRE library used to the log file.
  call gr_hypre_fortran_message(hypreMsg, HYPRE_MSG_LEN)
  call Logfile_stamp(hypreMsg , '[gr_hypreInit]')

  gr_hypreGridIsSetUp = .FALSE.   
  call RuntimeParameters_get("nbegin", gr_hypreNStep)
  gr_hypreNStep = gr_hypreNStep - 1

  ! Anisotropic diffusion, probably to represent anisotropic thermal conduction
  call RuntimeParameters_get("gr_hypreAnisoDiffusion",   gr_hypreAnisoDiffusion)
  
end subroutine gr_hypreInit
