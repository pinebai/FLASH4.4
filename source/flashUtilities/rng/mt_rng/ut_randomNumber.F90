!!****if* source/flashUtilities/rng/mt_rng/ut_randomNumber
!!
!! NAME
!!
!!  ut_randomNumber
!!
!! SYNOPSIS
!!
!!  call ut_randomNumber(real, intent(OUT)  :: x)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   x : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!! generates a random number on [0,1) with 53-bit resolution
  subroutine ut_randomNumber(x)
    implicit none
    real, intent(OUT) :: x
    integer,parameter :: realEightKind = selected_real_kind(15)
    real(kind=realEightKind)::y
    call ut_rand(y)
    x=y
  end subroutine ut_randomNumber
