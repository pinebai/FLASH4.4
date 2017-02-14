!!****h* source/physics/utilities/PlasmaState/PlasmaState_interface
!!
!! NAME
!!
!!  PlasmaState_interface
!!
!! SYNOPSIS
!!
!!   use PlasmaState_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Energy Deposition unit that defines its
!!  public interfaces.
!!
!!***

Module PlasmaState_interface

  interface
     subroutine PlasmaState (blockCount, blockList, timeStep, timeSimulation, passSplitDriver)
       integer, intent (in)           :: blockCount
       integer, intent (in)           :: blockList (1:blockCount)
       real,    intent (in)           :: timeStep
       real,    intent (in)           :: timeSimulation
       integer, intent (in), optional :: passSplitDriver
     end subroutine PlasmaState
     subroutine PlasmaState_getComposition(vecZ, vecA, vecFrac, nFrac, solnState)
       implicit none

       real, intent(OUT), dimension(:) :: vecZ, vecA, vecFrac
       integer, intent(OUT)            :: nFrac
       real, intent(IN),  dimension(:) :: solnState

     end subroutine PlasmaState_getComposition
  end interface

  interface
     subroutine PlasmaState_finalize ()    
     end subroutine PlasmaState_finalize
  end interface

  interface
     subroutine PlasmaState_init ()
     end subroutine PlasmaState_init
  end interface

end Module PlasmaState_interface
