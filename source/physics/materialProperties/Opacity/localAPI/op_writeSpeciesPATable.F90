!!****if* source/physics/materialProperties/Opacity/localAPI/op_writeSpeciesPATable
!!
!! NAME
!!
!!  op_writeSpeciesPATable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesPATable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Planck Absorption Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesPATable (fileUnit,species)

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  return
end subroutine op_writeSpeciesPATable
