!!****h* source/diagnostics/ThomsonScattering/ThomsonScattering_interface
!!
!! NAME
!!
!!  ThomsonScattering_interface
!!
!! SYNOPSIS
!!
!!   use ThomsonScattering_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Thomson Scattering iagnostic unit that defines its
!!  public interfaces.
!!
!!***

Module ThomsonScattering_interface

  interface
     subroutine ThomsonScattering (blockCount, blockList, timeSimulation)
       integer, intent (in) :: blockCount
       integer, intent (in) :: blockList (1:blockCount)
       real,    intent (in) :: timeSimulation
     end subroutine ThomsonScattering
  end interface

  interface
     subroutine ThomsonScattering_finalize ()    
     end subroutine ThomsonScattering_finalize
  end interface

  interface
     subroutine ThomsonScattering_init ()
     end subroutine ThomsonScattering_init
  end interface

end Module ThomsonScattering_interface
