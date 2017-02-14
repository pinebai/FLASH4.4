!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/coeff
!!
!! NAME
!! 
!!  coeff
!!
!!
!! SYNOPSIS
!!
!!  call coeff(integer(IN) :: numIntCells, 
!!             integer(IN) :: numCells, 
!!             real(IN)    :: dx(numCells), 
!!             real(OUT)   :: coeff1(numCells), 
!!             real(OUT)   :: coeff2(numCells), 
!!             real(OUT)   :: coeff3(numCells), 
!!             real(OUT)   :: coeff4(numCells), 
!!             real(OUT)   :: coeff5(numCells))
!!
!! 
!! DESCRIPTION
!!  
!!  Calculate the coefficients of quartic interpolation polynomial
!!  (Eq. 1.6 from Colella & Woodward). 
!!
!!  coeff1 is the coefficient of the first term {a(i+1) - a(i)} in
!!  Eq. 1.7 of Colella & Woodward.  
!!
!!  coeff2 is the coefficient of the second term {a(i) - a(i-1)} in 
!!  Eq. 1.7 of Colella & Woodward.
!!
!!  coeff3 is the entire coefficient (including that before the large
!!  curly braces) of the {\delta a(i+1)} term in Eq. 1.6 of Colella
!!  & Woodward.
!!
!!  coeff4 is the entire coefficient (including that before the large
!!  curly braces) of the {\delta a(i)} term in Eq. 1.6 of Colella 
!!  & Woodward.
!!
!!  coeff5 is the entire coefficient in front of all of the 
!!  a(j+1) - a(j)} terms (there are two of them) in Eq. 1.6 of Colella & 
!!  Woodward.  In that Eq. there are two coefficients, since this term is 
!!  split into two.  Here, we merge these into a single term and return 
!!  the total coefficient.
!!
!! ARGUMENTS
!!
!! numIntCells : 
!! numCells :
!! dx : 
!! coeff1 : Coella & Woodward coefficient 1.
!! coeff2 : Coella & Woodward coefficient 2.
!! coeff3 : Coella & Woodward coefficient 3.
!! coeff4 : Coella & Woodward coefficient 4.
!! coeff5 : Coella & Woodward coefficient 5.
!!
!!
!!***

subroutine coeff(numIntCells,numCells, dx,coeff1,coeff2,coeff3,coeff4,coeff5)

  implicit none

  integer, intent(IN) :: numIntCells,numCells
  real, intent(IN), DIMENSION(numCells) :: dx
  real, intent(OUT), DIMENSION(numCells) :: coeff1,coeff2,coeff3,coeff4,coeff5
  
  real, DIMENSION(numCells) :: scrch1(numCells), scrch2(numCells), scrch3(numCells)
  integer :: i, numIntCells6, numIntCells7, numIntCells8
  real :: temporary

!-----------------------------------

    numIntCells6 = numIntCells + 6
    numIntCells7 = numIntCells + 7
    numIntCells8 = numIntCells + 8

    do  i = 2, numIntCells8
       scrch1(i) = dx(i)     + dx(i-1)
       scrch2(i) = scrch1(i) + dx(i)
       scrch3(i) = scrch1(i) + dx(i-1)
    enddo


    do i = 2, numIntCells7
       temporary = dx(i)  /  ( scrch1(i) + dx(i+1) )

       coeff1(i) = temporary * scrch3(i)   / scrch1(i+1)

       coeff2(i) = temporary * scrch2(i+1) / scrch1(i)

    end do

 
    do i = 2, numIntCells6
       temporary = 1. / ( scrch1(i) + scrch1(i+2) )

       coeff3(i) = -temporary * dx(i)   * scrch1(i)   / scrch3(i+1)

       coeff4(i) =  temporary * dx(i+1) * scrch1(i+2) / scrch2(i+1)

       coeff5(i) = dx(i) - 2.e00 * (dx(i+1) * coeff3(i) + dx(i) * coeff4(i))
       coeff5(i) = coeff5(i) / scrch1(i+1)

    end do

    return
end subroutine coeff
