!!****ih* source/flashUtilities/interpolation/oneDim/ut_interpolationInterface
!!
!! NAME
!!
!!  ut_interpolationInterface
!!
!! SYNOPSIS
!!
!!  use ut_interpolationInterface
!!
!! DESCRIPTION
!!
!!  Interface module for some one-dimensional interpolation utilities.
!!
!!***

! Modification history:
!     Created   June 2007  KW

module ut_interpolationInterface

  implicit none
      

  interface
     subroutine ut_fndpos (check_monot, xarr, n, nl, nu, x, ix, ierr)
       logical, intent(IN)  :: check_monot
       integer, intent(IN)  :: n, nl, nu
       real,    intent(IN)  :: x, xarr(n)

       integer, intent(OUT) :: ix, ierr
     end subroutine ut_fndpos
  end interface

  interface
     subroutine ut_hunt(xx,n,x,low) 
       integer, INTENT(in)            :: n
       real, INTENT(in), DIMENSION(n) :: xx
       real, INTENT(in)               :: x
       integer, INTENT(inout)         :: low
     end subroutine ut_hunt
  end interface

  interface
     subroutine ut_polint(xa, ya, n, x, y, dy)
       real, intent(in)    :: x           ! where function is to be eval'd
       integer,intent(in)  :: n           ! number of input location/value pairs
       real, intent(in)    :: ya(*), xa(*)     ! function values and locations
       real, intent(out)   :: y              ! function evaluated at x
       real, intent(inout) :: dy             ! ignored and unchanged
     end subroutine ut_polint
  end interface

  interface
     subroutine ut_quadraticInterpol(x, y, xint, yint)
       real, intent(in)    :: xint           ! where function is to be eval'd
       real, intent(in)    :: y(*), x(*)     ! function values and locations
       real, intent(out)   :: yint           ! function evaluated at xint
     end subroutine ut_quadraticInterpol
  end interface

  interface
     subroutine ut_quadraticCellAverageInterpol(x, xl, y, xint, yint)
       real, intent(in)    :: xint           ! where function is to be eval'd
       real, intent(in)    :: y(*), x(*), xl(*) ! function values and locations
       real, intent(out)   :: yint           ! function evaluated at xint
     end subroutine ut_quadraticCellAverageInterpol
  end interface

end module ut_interpolationInterface
