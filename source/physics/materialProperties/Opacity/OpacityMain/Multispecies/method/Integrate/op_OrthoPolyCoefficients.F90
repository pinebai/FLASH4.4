!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_OrthoPolyCoefficients
!!
!! NAME
!!
!!  op_OrthoPolyCoefficients
!!
!! SYNOPSIS
!!
!!  call op_OrthoPolyCoefficients (integer (in)    :: n,
!!                                 integer (in)    :: m,
!!                                 integer (in)    :: nRow,
!!                                 real    (in)    :: Mom,
!!                                 real    (in)    :: A,
!!                                 real    (in)    :: B,
!!                                 real    (inout) :: Row1,
!!                                 real    (inout) :: Row2,
!!                                 real    (out)   :: OrthoA,
!!                                 real    (out)   :: OrthoB)
!!
!! DESCRIPTION
!!
!!  Calculates the orthogonal polynomial 3-term recursion coefficients.
!!
!! ARGUMENTS
!!
!!  n        : the order of the polynomial
!!  m        : the number of modified moments
!!  nRow     : the declared dimension for the row vectors
!!  Mom      : the modified moments
!!  A        : the auxilliary polynomial A coefficients
!!  B        : the auxilliary polynomial B coefficients
!!  Row1     : row vector that will contain intermediate quantities
!!  Row2     : row vector that will contain intermediate quantities
!!  OrthoA   : the orthogonal polynomial A coefficients
!!  OrthoB   : the orthogonal polynomial B coefficients
!!
!!***
subroutine op_OrthoPolyCoefficients (n,m,nRow,Mom,A,B,Row1,Row2,OrthoA,OrthoB)

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent (in)    :: n
  integer, intent (in)    :: m
  integer, intent (in)    :: nRow

  real,    intent (in)    :: Mom    (1:m)
  real,    intent (in)    :: A      (1:m)
  real,    intent (in)    :: B      (1:m)
  real,    intent (inout) :: Row1   (1:nRow)
  real,    intent (inout) :: Row2   (1:nRow) 
  real,    intent (out)   :: OrthoA (1:n)
  real,    intent (out)   :: OrthoB (1:n)    ! full dimension given but only 2:n addressed

  integer :: i,j
  integer :: imax,jmax

  real    :: sigma
  real    :: theta
!
!
!   ...Check conformity between order of polynomial and number of moments. Check also,
!      if sufficient memory is available to hold the row vector quantities.
!
!
  if (m /= (n + n - 1)) then
      call Driver_abortFlash ('[op_OrthoPolyCoefficients] ERROR: # of moments / coefficients mismatch')
  end if

  if (m > nRow) then
      call Driver_abortFlash ('[op_OrthoPolyCoefficients] ERROR: Row vector size too small')
  end if
!
!
!   ...Proceed with Sack & Donovan algorithm in Wheeler's form.
!
!
  select case (n)

    case (1)

      OrthoA (1) = Mom (1) + A (1)

    case (2)

      sigma = Mom (1) + A (1)
      theta = (A (2) - sigma) * Mom (1) + Mom (2) + B (2)

      OrthoA (1) = sigma
      OrthoA (2) = (((A (3) - sigma) * Mom (2) + Mom (3) + B (3) * Mom (1)) / theta) - Mom (1) + A (2)
      OrthoB (2) = theta

    case default

      imax = n - 1
      jmax = m

      do j = 1,jmax
         Row1 (j) = Mom (j)
      end do

      sigma = Row1 (1) + A (1)
      OrthoA (1) = sigma
!
!
!   ...Evaluate the 2nd row of the terminal matrix.
!
!
      jmax = jmax - 1

      Row2 (1) = (A (2) - sigma) * Row1 (1) + Row1 (2) + B (2)
      do j = 2,jmax
         Row2 (j) = (A (j+1) - sigma) * Row1 (j) + Row1 (j+1) + B (j+1) * Row1 (j-1)
      end do

      sigma = A (2) - Row1 (1) + (Row2 (2) / Row2 (1))
      theta = Row2 (1)

      OrthoA (2) = sigma
      OrthoB (2) = theta
!
!
!   ...Evaluate higher rows.
!
!
      do i = 2,imax

         jmax = jmax - 1

         if (mod (i,2) == 0) then

             do j = i,jmax
                Row1 (j) = (A (j+1) - sigma) * Row2 (j) + Row2 (j+1) + B (j+1) * Row2 (j-1) - theta * Row1 (j)
             end do

             sigma = A (i+1) - (Row2 (i) / Row2 (i-1)) + (Row1 (i+1) / Row1 (i))
             theta = Row1 (i) / Row2 (i-1)

         else

             do j = i,jmax
                Row2 (j) = (A (j+1) - sigma) * Row1 (j) + Row1 (j+1) + B (j+1) * Row1 (j-1) - theta * Row2 (j)
             end do

             sigma = A (i+1) - (Row1 (i) / Row1 (i-1)) + (Row2 (i+1) / Row2 (i))
             theta = Row2 (i) / Row1 (i-1)

         end if

         OrthoA (i+1) = sigma
         OrthoB (i+1) = theta

      end do

  end select
!
!
!   ...Ready! 
!
!
  return
end subroutine op_OrthoPolyCoefficients
