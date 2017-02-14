!!****if* source/physics/Gravity/GravityMain/Constant/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for constant gravity.
!!
!! PARAMTERS
!!
!!   grv_direc , string : specifies the direction of constant gravity
!!   grv_vector : value of constant gravity
!!
!!***

module Gravity_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save :: grv_const
  character (len=80), save :: grv_direc

  !! *** Module Variables *** !!

  real, dimension(MDIM), save :: grv_vector
  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity

end module Gravity_data
