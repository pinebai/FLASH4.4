!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commGetNewRays
!!
!!  NAME     
!!   ed_commGetNewRays
!!
!!  SYNOPSIS
!!   ed_commGetNewRays()
!!
!!  DESCRIPTION 
!!    This subroutine copies newly received rays from the internal
!!    pipeline buffer into ed_rays.
!!
!!  ARGUMENTS
!!
!!  NOTES
!!    Rays must be sorted both before and after this subroutine.
!!
!!  SIDE EFFECTS
!!    ed_rays and ed_rayCount may be updated.
!!
!!***

subroutine ed_commGetNewRays()
  implicit none
end subroutine ed_commGetNewRays
