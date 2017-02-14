!!****if* source/physics/materialProperties/Opacity/localAPI/op_initNumerics
!!
!! NAME
!!
!!  op_initNumerics
!!
!! SYNOPSIS
!!
!!  call op_initNumerics ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the section of the numerics for the opacity unit. It is here where
!!  several constants are set to enable smooth integrations without dangers of
!!  computational over- or underflow. Also several constants appearing all over the
!!  place like one, zero, etc. are set here.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initNumerics ()

  implicit none

  return
end subroutine op_initNumerics
