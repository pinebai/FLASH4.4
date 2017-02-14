!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/hy_interpNoLim
!!
!! NAME
!! 
!!  hy_interpNoLim
!!
!! SYNOPSIS
!!
!!  call hy_interpNoLim(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells,
!!              real(OUT)   :: al(numCells), 
!!              real(IN)    :: a(numCells), 
!!              real(OUT)   :: ar(numCells), 
!!              real(IN)    :: coeff1(numCells), 
!!              real(IN)    :: coeff2(numCells), 
!!              real(IN)    :: coeff3(numCells), 
!!              real(IN)    :: coeff4(numCells), 
!!              real(IN)    :: coeff5(numCells))
!!
!! DESCRIPTION
!!
!!  Interpolate interface values and do not monotonize
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells :
!! al :
!! a :
!! ar :
!! coeff1 :
!! coeff2 :
!! coeff3 :
!! coeff4 :
!! coeff5 :
!!
!!
!! SIDE EFFECTS
!!
!!  Modifies array hy_dela, which is referenced in subroutine detect.
!!
!!***

subroutine hy_interpNoLim (numIntCells, numCells, al, a, ar, &
                   coeff1, coeff2, coeff3, coeff4, coeff5)


  use Hydro_data, ONLY:  hy_dela
                     
  implicit none
  integer, intent(IN) :: numIntCells, numCells
  real, intent(IN), DIMENSION(numCells) ::  a, coeff1, coeff2, coeff3, coeff4, coeff5
  real, intent(OUT), DIMENSION(numCells) ::  al, ar
  
  real, dimension(numCells) :: scrch1, scrch2, scrch3, scrch4
  integer :: i, numIntCells5, numIntCells6, numIntCells8

    numIntCells5 = numIntCells + 5
    numIntCells6 = numIntCells + 6
    numIntCells8 = numIntCells + 8
    

! compute some common factors
! DONGWOOK: HERE, WE NEED TO USE CHARACTERISTIC VARIABLES FOR DIFFERENCING, NOT PRIMITIVE VARIABLES,
!           BY MULTIPLYING LEFT EIGENVECTORS TO PRIMITIVE VARIABLES
    do i = 2, numIntCells8
       scrch1(i) = a(i) - a(i-1)
       scrch2(i) = abs ( scrch1(i) + scrch1(i) )
       scrch4(i) = sign (1.e00, scrch1(i))
    end do


! apply Eq. 1.8 of Colella & Woodward -- guarantee that a(i+1/2) lies
! between a(i) and a(i+1)
    
    do i = 2, numIntCells6+1
       hy_dela(i)   = coeff1(i) * scrch1(i+1) + coeff2(i) * scrch1(i)
       
!!$       if (hy_dela(i) .LT. 0.e0) then
!!$          scrch3(i) = -1.e0
!!$       else
!!$          scrch3(i) = 1.e0
!!$       endif
!!$       
!!$       hy_dela(i)   = min(abs(hy_dela(i)), scrch2(i), scrch2(i+1))* scrch3(i)
!!$
!!$
!!$       if (-scrch4(i)*scrch4(i+1) >= 0.e0) hy_dela(i) = 0.e0
       
    end do


! assemble the left and right interface values using the quartic polynomial,
! Eq. 1.6 from Colella & Woodward
    
    do i = 2, numIntCells6
       ar(i)  = a (i) + coeff5(i) * scrch1(i+1) + coeff3(i) * hy_dela(i+1)
       ar(i)  = ar(i) + coeff4(i) * hy_dela(i)
       al(i+1)= ar(i)
    end do
    
    return
  end subroutine hy_interpNoLim






