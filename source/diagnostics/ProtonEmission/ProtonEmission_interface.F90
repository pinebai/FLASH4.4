!!****h* source/diagnostics/ProtonEmission/ProtonEmission_interface
!!
!! NAME
!!
!!  ProtonEmission_interface
!!
!! SYNOPSIS
!!
!!   use ProtonEmission_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Proton Emission unit that defines its
!!  public interfaces.
!!
!!***

Module ProtonEmission_interface

  interface
     subroutine ProtonEmission (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ProtonEmission
  end interface

  interface
     subroutine ProtonEmission_finalize ()    
     end subroutine ProtonEmission_finalize
  end interface

  interface
     subroutine ProtonEmission_init ()
     end subroutine ProtonEmission_init
  end interface

end Module ProtonEmission_interface
