!!****f* source/diagnostics/ProtonEmission/ProtonEmission
!!
!! NAME
!!
!!  ProtonEmission
!!
!! SYNOPSIS
!!
!!  call ProtonEmission (integer, intent (in) :: blockCount, 
!!                       integer, intent (in) :: blockList (:),
!!                       real,    intent (in) :: timeStep,
!!                       real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Creates protons (emission) inside the domain and records them on detector screen(s).
!!  This is the driver routine for proton emission at a particular time during the
!!  simulation. The domain structure and properties are assumed to remain frozen during
!!  the travel of the protons through the domain.
!!
!!  Currently the protons do not interact with the domain, hence no domain update
!!  is necessary.
!!
!!  The code consists of the following basic steps:
!!
!!         1) Create the protons inside the domain
!!         2) Follow the protons motion through the domain
!!         3) Record protons on the detector screen(s)
!!         4) Write each screen proton to the detector file(s)
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : The current time step duration
!!  timeSimulation : The current simulation time
!!
!! NOTES
!!          
!!  The current implementation can handle more than one detector screen. The proton
!!  paths are calculated using classical Newton mechanics and are deflected due to
!!  average electrical and magnetic fields in each cell.
!!
!!***

subroutine ProtonEmission (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ProtonEmission
