!!****if* source/flashUtilities/interpolation/oneDim/ut_fndpos
!!
!! NAME
!!
!!  ut_fndpos
!!
!! SYNOPSIS
!!
!!  ut_fndpos (check_monot, xarr, n, nl, nu, x, ix, ierr)
!!
!!  ut_fndpos (logical, real(), integer, integer, integer,
!!                                         real(), integer, integer)
!!
!! DESCRIPTION
!!
!!  This routine returns the index IX of a value X within an array
!!  XARR of length N, such that X lies in the open interval
!!  [IX,IX+1).  A region of indices where to search can be given
!!  by a lower and an upper index NL and NU, respectively. XARR must be
!!  strictly monotonic increasing within the interval [NL,NU].
!!
!!  Monotonicity of XARR will be checked if CHECK_MONOT is .TRUE. An error
!!  code IERR = -1 is returned if either X lies outside of [XARR(NL),XARR(NU)]
!!  or XARR is not strictly monotonically increasing within this interval
!!  (if CHECK_MONOT).
!! 
!! ARGUMENTS
!!
!!  xarr            target array [in]
!!
!!  n               dimension of xarr, xarr(1...n) [in]
!!
!!  nl, nu          lower and upper bounds for search interval
!!                  nl and nu must be inside [1,n] [in]
!!
!!  x               target location [in]
!!
!!  ix              target location in index space [out]
!!
!!  check_monot     monotonicity test flag [in];
!!                  whether this subroutine shall check for monotonicity of xarr.
!!
!!  ierr            error code (-1 indicates failure) [out]
!!
!! USED BY
!!
!!  setups/wd_convect/init_block.F90        in FLASH2
!!  setups/wd_convect/spherical_hse_1d.F90  in FLASH2
!!  physics/Cosmology/CosmologyMain/MatterLambdaKernel/CosmologicalFunctions in FLASH3
!!
!! HISTORY
!!
!!    This was in FLASH2 as source/util/interp/fndpos.F90, 
!!                    added by Tomek Plewa 2003-03-27  (FLASH2.3)
!!    Changed name to ut_fndpos       - KW 2007-06-07
!!    Changed check_monot to logical  - KW 2007-06-07
!!
!!***

subroutine ut_fndpos (check_monot, xarr, n, nl, nu, x, ix, ierr)

implicit none

logical, intent(IN)  :: check_monot
integer, intent(IN)  :: n, nl, nu
real,    intent(IN)  :: x, xarr(n)

integer, intent(OUT) :: ix, ierr

integer :: i, il, iu, im

!-------------------------------------------------------------------------------

        ierr = 0

!       --------------------
!       check selected range

        if ( x.lt.xarr(nl) ) then
           ix   = nl
           ierr = -1
           return
        else if ( x.gt.xarr(nu) ) then
           ix   = nu
           ierr = -1
           return
        end if

        if ( check_monot ) then

           do i = nl,nu-1
              if ( xarr(i).ge.xarr(i+1) ) then
                 ierr = -1
                 return
              end if
           end do

        end if

!       ------
!       search

        if ( x.eq.xarr(nu) ) then

           ix = nu

        else

           il = nl
           iu = nu

20           if ( iu-il.gt.1 ) then
              im = (iu + il) / 2
              if ( x.lt.xarr(im) ) then
                 iu = im
              else
                 il = im
              end if
              go to 20
           endif

           ix = il

        end if

        return

end subroutine ut_fndpos
