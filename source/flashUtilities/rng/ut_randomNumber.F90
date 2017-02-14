!!****if* source/flashUtilities/rng/ut_randomNumber
!!
!! NAME
!!
!!  ut_randomNumber
!!
!! SYNOPSIS
!!
!!  call ut_randomNumber(real, intent(OUT)  :: r1)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   r1 : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

  subroutine ut_randomNumber(r1)
    implicit none
    real, intent(OUT) :: r1
    real :: r2
    call random_number(r2)
    r1=r2
  end subroutine ut_randomNumber
