!!****h* source/physics/sourceTerms/EnergyDeposition/EnergyDeposition_interface
!!
!! NAME
!!
!!  EnergyDeposition_interface
!!
!! SYNOPSIS
!!
!!   use EnergyDeposition_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Energy Deposition unit that defines its
!!  public interfaces.
!!
!!***

Module EnergyDeposition_interface

  interface
     subroutine EnergyDeposition (blockCount, blockList, timeStep, timeSimulation, passSplitDriver)
       integer, intent (in)           :: blockCount
       integer, intent (in)           :: blockList (1:blockCount)
       real,    intent (in)           :: timeStep
       real,    intent (in)           :: timeSimulation
       integer, intent (in), optional :: passSplitDriver
     end subroutine EnergyDeposition
  end interface

  interface
     subroutine EnergyDeposition_finalize ()    
     end subroutine EnergyDeposition_finalize
  end interface

  interface
     subroutine EnergyDeposition_init ()
     end subroutine EnergyDeposition_init
  end interface

end Module EnergyDeposition_interface
