!!****if* source/flashUtilities/interpolation/oneDim/ut_hunt
!!
!!  NAME
!!    ut_hunt  
!!  SYNOPSIS 
!!    call ut_hunt(xx, n, x, low)
!!  FUNCTION 
!!    Given an array xx of length n and a value x, this routine returns a value 
!!    low such that x is between xx(low) and xx(low+1). the array xx must be  
!!    monotonic. low=0 or low=n indicates that x is out of range. on input low is a  
!!    guess for the table entry, and should be used as input for the next hunt.  
!!  INPUTS
!!      xx - real, sorted, array of dimension n to search through
!!      n  - size of xx
!!      x  - value to search for
!!      low - largest value in [1,n] s.t. xx[low] < x
!!  RESULT
!!      most errors return low = 1
!!
!!  HISTORY
!!
!!    This was in FLASH2 as setups/nova/non-numrecipies/hunt.F90.
!!    Changed name to ut_hunt         - KW 2007-06-07
!!    Implemented documented use of 'low' value passed in - KW 2007-06-07
!!***
      subroutine ut_hunt(xx,n,x,low) 

      implicit none

      integer, INTENT(in)            :: n
      real, INTENT(in), DIMENSION(n) :: xx
      real, INTENT(in)               :: x
      integer, INTENT(inout)         :: low

      logical                        :: dirn
      integer                        :: high, middle, step

! initialize step size to 1
      step = 1

! dirn is .true. if in ascending order, .false. otherwise
      dirn = xx(n).ge.xx(1)

      if (dirn) then
         if (x.le.xx(1)) then
            low = 0
            return
!!$         else if (x.eq.xx(1)) then
!!$            low = 1
!!$            return
         end if
         if (x.gt.xx(n)) then
            low = n
            return
         end if
      end if

! this should never, ever happen.
      if (xx(n).eq.xx(1)) then
         print *,'ut_hunt: array xx() is constant!'
         low = 1
         return
      endif

! same here
      if (n < 2) then
         print *,'list too small!'
         low = 1
         return
      endif


! check if low is a valid index
      if ((low.lt.1).or.(low.gt.n)) then
          low = 0
          high = n+1

! low is valid so use it
      else
! if ascending  and x >= xx(low) or
!    descending and x <  xx(low), then hunt up index
         if ((x.ge.xx(low)).eqv.(dirn)) then
            high = low + step
! make sure high is a valid index
            if (high.gt.n) then
               high = n+1
            else
               do while ((x.ge.xx(high)).eqv.(dirn)) 
                  low = high
                  step = 2*step ! double step
                  high = low + step
                  if (high.gt.n) then ! stop hunting
                     high = n+1
                     exit
                  endif
               enddo
            endif
         else ! hunt down in index
            high = low
            low = high - step
! make sure low is a valid index
            if (low.lt.1) then
               low = 0
            else
! continue hunting down if x < xx(low) and xx is ascending or
!                       if x >= xx(low) and xx is descending
               do while ((x.lt.xx(low)).eqv.(dirn))
                  high = low
                  step = 2*step ! double step
                  low = high - step
                  if (low.lt.1) then ! stop hunting
                     low = 0
                     exit
                  endif
               enddo
            endif
         endif
      endif
! after hunt, target value is guaranteed to be 
! xx(low) > x >= xx(high) or
! xx(low) <= x < xx(high)

! now do bisection
      do while ((high - low).gt.1)

         middle = (high + low)/2
         if ((x.ge.xx(middle).eqv.dirn)) then
            low = middle
         else
            high = middle
         endif

      enddo

! make sure low is returned correctly for values on end of table
      if (x.eq.xx(n)) low = n-1
      if (x.eq.xx(1)) low = 1
      return

      end subroutine ut_hunt
