!!****if* source/physics/Gravity/GravityMain/UserDefined/Gravity_data
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
!!  Stores the local data for user-defined gravity.
!!
!! PARAMTERS
!!
!!
!!***

module Gravity_data

  implicit none

#include "constants.h"

  !! *** Runtime Parameters *** !!

  !! *** Module Variables *** !!

  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity

end module Gravity_data
