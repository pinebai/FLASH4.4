!!****if* source/physics/materialProperties/Opacity/localAPI/op_writeSpeciesROTable
!!
!! NAME
!!
!!  op_writeSpeciesROTable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesROTable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Rosseland Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesROTable (fileUnit,species)

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  return
end subroutine op_writeSpeciesROTable
