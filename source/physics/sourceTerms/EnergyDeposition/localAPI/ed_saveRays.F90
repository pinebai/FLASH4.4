!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_saveRays
!!
!! NAME
!!
!!  ed_saveRays
!!
!! SYNOPSIS
!!
!!  call ed_saveRays ()
!!
!! DESCRIPTION
!!
!!  This routine saves information of those rays which exited the domain, if requested by
!!  the user. If the user specified no saving action, the routine simply sets all those
!!  rays that exited the domain to nonexistent rays and returns. Otherwise the info is
!!  copied into the saved rays array.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_saveRays ()

  implicit none

  return
end subroutine ed_saveRays
