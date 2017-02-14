!!****if* source/flashUtilities/rng/ut_randomSeed
!!
!! NAME
!!
!!  ut_randomSeed
!!
!! SYNOPSIS
!!
!!  call ut_randomSeed(integer, optional, intent(OUT)  :: ut_size,
!!                     integer, optional, dimension(:), intent(IN)  :: ut_put,
!!                     integer, optional, dimension(:), intent(OUT)  :: ut_get)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   ut_size : 
!!
!!   ut_put : 
!!
!!   ut_get : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

#ifndef __ABSOFT__


subroutine ut_randomSeed(ut_size,ut_put, ut_get)
  implicit none
  integer, optional, dimension(:), intent(IN) :: ut_put
  integer, optional, dimension(:), intent(OUT) :: ut_get
  integer, optional, intent(OUT) :: ut_size

  if(present(ut_size))call random_seed(size=ut_size)
  if(present(ut_put))call random_seed(put=ut_put)
  if(present(ut_get))call random_seed(get=ut_get)
end subroutine ut_randomSeed



#else

!The original subroutine (above) cannot be compiled with Absoft 64-bit
!Pro Fortran 11.1.4 on code.uchicago.edu which uses Centos-6.0 x86_64.
!I have created a temporary version of this subroutine which can be
!compiled by Absoft compiler.

subroutine ut_randomSeed(ut_size,ut_put, ut_get)
  implicit none
  integer, optional, dimension(:), intent(IN) :: ut_put
  integer, optional, dimension(:), intent(OUT) :: ut_get
  integer, optional, intent(OUT) :: ut_size
  integer, allocatable, dimension(:) :: ut_put_tmp

  if (present(ut_size)) call random_seed(size=ut_size)
  if (present(ut_put)) then
     allocate(ut_put_tmp(lbound(ut_put,1):ubound(ut_put,1)))
     ut_put_tmp = ut_put
     call random_seed(put=ut_put_tmp)
     deallocate(ut_put_tmp)
  end if
  if (present(ut_get)) call random_seed(get=ut_get)
end subroutine ut_randomSeed

#endif
