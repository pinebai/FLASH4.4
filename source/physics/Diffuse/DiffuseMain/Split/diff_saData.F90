!!****if* source/physics/Diffuse/DiffuseMain/Split/diff_saData
!!
!! NAME
!!
!!  diff_saData
!!
!! SYNOPSIS
!!  use diff_saData
!!
!! DESCRIPTION
!!
!!  Defines and stores local data for the Split implementation of the main subunit of unit Diffuse.
!!  All the variables defined here should be initialized in Diffuse_init() by calling
!!  the RuntimeParameters_get subroutine or similar.
!!  
!!***


Module diff_saData

  use Grid_interface, ONLY: GRID_PDE_BND_PERIODIC

  implicit none

  logical, save :: updateDiffuse = .TRUE.
  integer, save :: diff_boundary = GRID_PDE_BND_PERIODIC

  real, save :: diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime

  integer,save,dimension(6) :: diff_eleDomainBC
 
  real, save :: diff_thetaImplct
 
end Module diff_saData
