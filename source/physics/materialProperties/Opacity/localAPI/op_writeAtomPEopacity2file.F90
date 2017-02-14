!!****if* source/physics/materialProperties/Opacity/localAPI/op_writeAtomPEopacity2file
!!
!! NAME
!!
!!  op_writeAtomPEopacity2file
!!
!! SYNOPSIS
!!
!!  call op_writeAtomPEopacity2file (integer (in) :: Z,
!!                                   real    (in) :: Elower,
!!                                   real    (in) :: Eupper,
!!                                   integer (in) :: nPoints)
!!
!! DESCRIPTION
!!
!!  This routine writes the photoelectron opacities for element Z between the
!!  energy range Elower - Eupper in nPoints points to a specific file.
!!  The routine is meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!  Z       : atomic number
!!  Elower  : lower energy bound (in keV)
!!  Eupper  : upper energy bound (in keV)
!!  nPoints : number of points to be written to file
!!
!!***
subroutine op_writeAtomPEopacity2file (Z, Elower, Eupper, nPoints)

  implicit none

  integer, intent (in) :: nPoints
  integer, intent (in) :: Z
  real,    intent (in) :: Elower
  real,    intent (in) :: Eupper

  return
end subroutine op_writeAtomPEopacity2file
