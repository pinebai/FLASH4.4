!!****if* source/diagnostics/ProtonImaging/localAPI/pi_transportBeamProtons
!!
!! NAME
!!
!!  pi_transportBeamProtons
!!
!! SYNOPSIS
!!
!!  call pi_transportBeamProtons (integer, intent (in) :: blockCount, 
!!                                integer, intent (in) :: blockList (:), 
!!                                real,    intent (in) :: timeStep,
!!                                real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Launches batches of beam protons onto the domain and transports them through
!!  the domain. Two things can happen to each proton: 1) the proton does not
!!  leave the domain during this time step or 2) the proton leaves the domain
!!  and is recorded on the screen. The domain structure and properties are assumed
!!  to remain frozen during the current time step transportation of the protons.
!!  Currently the protons do not interact with the domain, hence no domain update
!!  is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Create the beam protons on the domain surface
!!         2) Follow the protons motion through the domain
!!         3) Record each proton leaving the domain on the detector screen(s)
!!         4) Store each proton remaining in the domain back to disk
!!         5) Write each screen proton to the detector file(s)
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
!!  The current implementation can handle more than one proton beam and more
!!  than one detector screen. The proton paths are calculated using classical
!!  Newton mechanics and are deflected due to average electrical and magnetic
!!  fields in each cell.
!!
!!***

subroutine pi_transportBeamProtons (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine pi_transportBeamProtons
