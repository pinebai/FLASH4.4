!!****if* source/flashUtilities/rng/ut_randomNumberArray
!!
!! NAME
!!
!!  ut_randomNumberArray
!!
!! SYNOPSIS
!!
!!  call ut_randomNumberArray(real, intent(OUT)  :: r1)
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

  subroutine ut_randomNumberArray(r1)
    implicit none
    real, intent(OUT) :: r1(:)
    call random_number(r1)
  end subroutine ut_randomNumberArray
