!!****f* source/diagnostics/ProtonImaging/ProtonImaging
!!
!! NAME
!!
!!  ProtonImaging
!!
!! SYNOPSIS
!!
!!  call ProtonImaging (integer(in) :: blockCount, 
!!                      integer(in),dimension(1:blockCount) :: blockList (:), 
!!                      real(in)    :: timeStep,
!!                      real(in)    :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Launches bundle(s) of protons onto the domain and records them on detector
!!  screen(s). This is the driver routine for proton imaging at a particular
!!  time during the simulation. The domain structure and properties are assumed
!!  to remain frozen during the traversion of the protons.
!!
!!  Currently the protons do not interact with the domain, hence no domain update
!!  is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Create the protons on the domain surface
!!         2) Follow the protons motion through the domain
!!         3) Record each proton on the detector screen(s)
!!         4) Write each screen proton to the detector file(s)
!!
!!  The code can handle a larger number of protons than there is memory available
!!  by sending bags of protons and reusing available memory.
!!
!! ARGUMENTS
!!
!!  blockCount      : Number of blocks on current processor
!!  blockList       : All block ID numbers
!!  timeStep        : The current time step duration
!!  timeSimulation  : current simulation time
!!
!! NOTES
!!          
!!  The current implementation can handle more than one proton beam and more
!!  than one detector screen. The proton paths are calculated using classical
!!  Newton mechanics and are deflected due to average electrical and magnetic
!!  fields in each cell.
!!
!!***

subroutine ProtonImaging (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ProtonImaging
