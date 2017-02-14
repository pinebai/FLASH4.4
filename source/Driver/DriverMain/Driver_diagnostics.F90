!!****if* source/Driver/DriverMain/Driver_diagnostics
!!
!! NAME
!!
!!  Driver_diagnostics
!!
!! SYNOPSIS
!!
!!  Driver_diagnostics (integer(IN):: blockCount,
!!                      integer(IN):: blockList(blockCount),
!!                      real   (IN):: dt)
!!
!! DESCRIPTION
!!
!!  Driver for diagnostics. Instead of calling all diagnostic routines 
!!  from Driver_evolveFlash we call Driver_diagnostics which then
!!  makes the calls to ProtonImaging, etc.  If a unit is not included
!!  in the simulation, the routine will be a stub and return without
!!  doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the diagnostics
!!  dt           : the current timestep
!!
!!***



subroutine Driver_diagnostics (blockCount, blockList, dt)

  use Driver_data,             ONLY : dr_simTime
  use ProtonImaging_interface, ONLY : ProtonImaging
  use ProtonEmission_interface, ONLY : ProtonEmission
  use ThomsonScattering_interface, ONLY : ThomsonScattering

  implicit none

  real,    intent(IN) :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList

  call ProtonImaging (blockCount, blockList, dt, dr_simTime)
  call ProtonEmission (blockCount, blockList, dt, dr_simTime)
  call ThomsonScattering (blockCount, blockList, dr_simTime)

  return
end subroutine Driver_diagnostics
