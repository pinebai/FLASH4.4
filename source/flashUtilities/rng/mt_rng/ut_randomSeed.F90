!!****if* source/flashUtilities/rng/mt_rng/ut_randomSeed
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

  subroutine ut_randomSeed(ut_size,ut_put, ut_get)

    implicit none

    integer, optional, dimension(:), intent(IN) :: ut_put
    integer, optional, dimension(:), intent(OUT) :: ut_get
    integer, optional, intent(OUT) :: ut_size
    logical :: noargs

    noargs = .not.(present(ut_size).or.present(ut_put).or.present(ut_get))
    if(noargs)then
       call ut_rand_init()
    else
       if(present(ut_size))then
          call ut_rand_init()
          call ut_get_size(ut_size)
       end if
       if(present(ut_get))call ut_get_seed(ut_get)
       if(present(ut_put))then
          print*,'calling put seed'
          call ut_put_seed(ut_put)
          print*,'done'
       end if
    end if
  end subroutine ut_randomSeed

