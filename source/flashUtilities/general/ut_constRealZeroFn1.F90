!!****if* source/flashUtilities/general/ut_constRealZeroFn1
!!
!! NAME
!!
!!  ut_constRealZeroFn1
!!  ut_constRealOneFn1
!!
!! SYNOPSIS
!!
!!  real a0 = ut_constRealZeroFn1 (real(in)  :: x)
!!
!!  real a1 = ut_constRealOneFn1  (real(in)  :: x)
!!
!! DESCRIPTION
!!
!!  Functions of one real argument that always returns
!!  the same value.
!!
!! ARGUMENTS
!!
!!  x         : ignored argument
!!
!! NOTES
!!
!!***

real function ut_constRealZeroFn1 (x)
  implicit none
  real,    intent (in)  :: x

  ut_constRealZeroFn1 = 0.0

end function ut_constRealZeroFn1

real function ut_constRealOneFn1 (x)
  implicit none
  real,    intent (in)  :: x

  ut_constRealOneFn1 = 1.0

end function ut_constRealOneFn1
