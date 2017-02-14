!!****if* source/numericalTools/Roots/RootsMain/Roots_data
!!
!! NAME
!!
!!  Roots_data
!!
!! SYNOPSIS
!!
!!  use Roots_data
!!  
!! DESCRIPTION
!!
!!  Data module for Roots unit
!!  --------------------------
!!   
!!   Legend: (P) means data that is set as (P)arameters
!!           (G) means data that is (G)et from other units (driver, physical constants, etc ...)
!!           (R) means data that is supplied as input through (R)untime parameters
!!           (C) means data that is (C)alculated internally by the Roots unit
!!           (I) means data that is (I)nitialized internally by the Roots unit
!!
!!   rt_macheps    (C,I) : the machine epsilon for the current precision
!!   rt_LPN        (C,I) : the (L)argest (P)ositive (N)umber for the current precision
!!   rt_sqrtLPN    (C,I) : the square root of rt_LPN
!!
!!***

Module Roots_data

  implicit none

  real, save :: rt_macheps
  real, save :: rt_LPN
  real, save :: rt_sqrtLPN

end Module Roots_data
