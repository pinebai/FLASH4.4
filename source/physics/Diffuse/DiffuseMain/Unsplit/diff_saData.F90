!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_saData
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
!!  Defines and stores local data for the UG implementation of the main subunit of unit Diffuse.
!!  All the variables defined here should be initialized in Diffuse_init() by calling
!!  RuntimeParameters_get subroutine or similar.
!!
!!  
!!***

#include "Flash.h"

Module diff_saData
  
  use Grid_interface, ONLY: GRID_PDE_BND_PERIODIC


  implicit none

#include "constants.h"

  logical, save :: updateDiffuse = .TRUE.
  integer, save :: diff_boundary = GRID_PDE_BND_PERIODIC

  real, save :: diff_scaleFactThermSaTempDiff, diff_scaleFactThermSaTime

  integer, save :: diff_eleDomainBC(6)
  integer, save :: diff_ionDomainBC(6)
 
  real, save :: diff_thetaImplct
  real, save :: diff_ionThetaImplct

  logical, save :: diff_updEint
  
end Module diff_saData
