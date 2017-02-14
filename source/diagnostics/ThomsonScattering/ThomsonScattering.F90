!!****f* source/diagnostics/ThomsonScattering/ThomsonScattering
!!
!! NAME
!!
!!  ThomsonScattering
!!
!! SYNOPSIS
!!
!!  call ThomsonScattering (integer, intent (in) :: blockCount, 
!!                      integer, intent (in) :: blockList (:), 
!!                      real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Simulates Thomson scattering diagnostics.
!!
!! ARGUMENTS
!!
!!  blockCount      : Number of blocks on current processor
!!  blockList       : All block ID numbers
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!
!!***

subroutine ThomsonScattering (blockCount, blockList, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation

  return
end subroutine ThomsonScattering
