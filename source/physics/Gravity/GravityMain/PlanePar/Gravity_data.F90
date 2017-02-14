!!****if* source/physics/Gravity/GravityMain/PlanePar/Gravity_data
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
!!  Stores the local data for PlanePar gravity.
!!
!!
!!
!!
!!***

module Gravity_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: grv_ptxpos, grv_ptmass, grv_newton
  integer, save :: grv_ptdirn
  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity

end module Gravity_data
