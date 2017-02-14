!!****if* source/diagnostics/ProtonImaging/localAPI/pi_transportDiskProtons
!!
!! NAME
!!
!!  pi_transportDiskProtons
!!
!! SYNOPSIS
!!
!!  call pi_transportDiskProtons (integer, intent (in) :: blockCount, 
!!                                integer, intent (in) :: blockList (:), 
!!                                real,    intent (in) :: timeStep,
!!                                real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Processes old batches of disk protons residing on disk from a previous time step.
!!  All protons that make it through the domain during the current time step are
!!  recorded on the detector screen(s). Any protons that stay in the domain are
!!  dumped back to disk as disk protons. Currently the protons do not interact with
!!  the domain, hence no domain update is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Read the protons from disk back into the domain
!!         2) Follow the protons motion through the domain
!!         3) Record each proton leaving the domain on the detector screen(s)
!!         4) Store each proton remaining in the domain back to disk
!!         5) Write screen protons to the detector file(s)
!!
!!  The code can handle a larger number of protons than there is memory available
!!  by sending batches of protons and reusing available memory.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : The current time step duration
!!  timeSimulation : current simulation time
!!
!! NOTES
!!          
!!  The proton paths are calculated using classical Newton mechanics and are deflected
!!  due to average electrical and magnetic fields in each cell.
!!
!!***

subroutine pi_transportDiskProtons (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine pi_transportDiskProtons
