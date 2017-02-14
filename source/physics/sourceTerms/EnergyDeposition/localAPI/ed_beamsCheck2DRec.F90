!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsCheck2DRec
!!
!! NAME
!!
!!  ed_beamsCheck2DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck2DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 2D rectangular
!!  grids (cartesian + cylindrical). All checks which depend on domain grid details should go
!!  in here. Currently it contains the following:
!!
!!         1) Check, if all beam line target areas are completely within the domain.
!!         2) Check, if all beam line lens areas are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck2DRec ()

  implicit none

  return
end subroutine ed_beamsCheck2DRec
