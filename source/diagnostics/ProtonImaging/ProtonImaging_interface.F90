!!****h* source/diagnostics/ProtonImaging/ProtonImaging_interface
!!
!! NAME
!!
!!  ProtonImaging_interface
!!
!! SYNOPSIS
!!
!!   use ProtonImaging_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Proton Imaging unit that defines its
!!  public interfaces.
!!
!!***

Module ProtonImaging_interface

  interface
     subroutine ProtonImaging (blockCount, blockList, timeStep, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeStep
       real,    intent (in) :: timeSimulation
     end subroutine ProtonImaging
  end interface

  interface
     subroutine ProtonImaging_finalize ()    
     end subroutine ProtonImaging_finalize
  end interface

  interface
     subroutine ProtonImaging_init ()
     end subroutine ProtonImaging_init
  end interface

end Module ProtonImaging_interface
