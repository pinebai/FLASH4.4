!!****if* source/physics/materialProperties/Opacity/localAPI/op_writeSpeciesPETable
!!
!! NAME
!!
!!  op_writeSpeciesPETable
!!
!! SYNOPSIS
!!
!!  call op_writeSpeciesPETable (integer (in) :: fileUnit,
!!                               integer (in) :: species)
!!
!! DESCRIPTION
!!
!!  Prints out the tabulated Planck Emission Opacities for the current species.
!!
!! ARGUMENTS
!!
!!  fileUnit    : unit # for the output file
!!  species     : the species index
!!
!!***
subroutine op_writeSpeciesPETable (fileUnit,species)

  implicit none

  integer, intent (in) :: fileUnit
  integer, intent (in) :: species

  return
end subroutine op_writeSpeciesPETable
