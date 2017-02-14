!!****if* source/Simulation/SimulationMain/Cool_Test/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt) 

  use Driver_data, ONLY: dr_simTime
  use Stir_interface, ONLY : Stir
  use Heat_interface, ONLY : Heat
  use Burn_interface, ONLY : Burn
  use Cool_interface, ONLY : Cool
  use Ionize_interface, ONLY : Ionize
  use PrimordialChemistry_interface, ONLY : PrimordialChemistry
  implicit none

  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList

  call Stir(blockCount, blockList, dt) 
  call Burn(blockCount, blockList, dt) 
  call Heat(blockCount, blockList, dt, dr_simTime) 
!  call Cool(blockCount, blockList,dt, dr_simTime)
  call PrimordialChemistry(blockCount, blockList, dt)
!  call Cool(blockCount, blockList, dt, dr_simTime)
  call Ionize(blockCount, blockList, dt, dr_simTime)

  return
end subroutine Driver_sourceTerms
