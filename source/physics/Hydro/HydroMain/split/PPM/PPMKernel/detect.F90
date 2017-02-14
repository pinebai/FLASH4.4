!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/detect
!!
!! NAME
!! 
!!  detect
!!
!!
!! SYNOPSIS
!!
!!  call detect(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells,
!!              integer(IN) :: guard,
!!              real(INOUT) :: al(numCells),
!!              real(IN)    :: a(numCells),
!!              real(INOUT) :: ar(numCells),
!!              real(IN)    :: smalla,
!!              real(IN)    :: rho(numCells),
!!              real(IN)    :: p(numCells),
!!              real(IN)    :: game(numCells),
!!              real(IN)    :: dx(numCells),
!!              real(IN)    :: x(numCells))
!!
!! DESCRIPTION
!!  
!!  Search for contact discontinuities in variable "a" and
!!  steepen the zone structure if needed
!!
!! ARGUMENTS
!!
!!  numIntCells :
!!  numCells :
!!  guard :
!!  al :
!!  a :
!!  ar :
!!  smalla :
!!  rho :
!!  p: 
!!  game :
!!  dx :
!!  x :
!!
!!***

subroutine detect(numIntCells, numCells, al,a,ar,smalla,rho,p,game,dx,x)

  use Hydro_data, ONLY: hy_small, hy_dela

  implicit none

  integer, intent(IN) :: numIntCells, numCells
  real, intent(INOUT), DIMENSION(numCells) :: al, ar
  real, intent(IN),    DIMENSION(numCells) :: a, rho, p, game, dx, x    
  real, dimension(numCells) ::  scrch1, scrch2, scrch3
  real, intent(IN) :: smalla
  integer :: i,numIntCells5,numIntCells6,numIntCells7,numIntCells8
  


    ! the following parameters are set as in Colella and
    ! Woodward (JCP, 54 (1984), 174), Eqs. 1.17 and 3.2

  real, PARAMETER :: eta1 = 20.e0, eta2 = 0.05e0, epsln = 0.01e0, ak0 = 0.1e0

  real :: tmp1, tmp2, tmp3

    !------------------------------------------------------------------------------

  numIntCells5 = numIntCells + 5
  numIntCells6 = numIntCells + 6
  numIntCells7 = numIntCells + 7
  numIntCells8 = numIntCells + 8

! compute some common factors      

  do  i = 2, numIntCells7
     scrch1(i) = dx(i) + dx(i-1)
     scrch2(i) = scrch1(i) + dx(i+1)
     scrch1(i) = (a(i) - a(i-1)) / scrch1(i)
  end do
      
    ! fill scrch2 with {\delta^2 a(i)} as given by Eq. 1.17 in Colella & Woodward
    
  do i = 2, numIntCells6
     scrch2(i) = (scrch1(i+1) - scrch1(i)) / scrch2(i)
     scrch1(i) = x(i) - x(i-1)
     scrch1(i) = scrch1(i) * scrch1(i) * scrch1(i)
  end do

! compute {\tilde \eta (i)} as given in the expression at the top of page 181
! in Colella & Woodward
  
  do i = 3, numIntCells5
     scrch3(i) = (scrch2(i-1) - scrch2(i+1)) * (scrch1(i) + scrch1(i+1))
     
     if (a(i+1) - a(i-1) == 0.e0) then
        tmp3 = hy_small * smalla
     else
        tmp3 = a(i+1) - a(i-1)
     endif
     
     scrch3(i) = scrch3(i) / ((x(i+1) - x(i-1)) * tmp3)

! scrch2 and scrch3 now contain finite difference approximations
! to the second and third derivativess of a.

! apply the first constaint on {\tidle \eta (i) as given in Eq. 1.17

     if (scrch2(i-1)*scrch2(i+1) >= 0.e0) scrch3(i) = 0.e0

! apply the second constraint

     tmp3 = epsln * min(a(i+1),a(i-1)) - abs(a(i+1) - a(i-1))
     
     if (tmp3 >= 0.e0) scrch3(i) = 0.e0
     
     scrch3(i) = max(0.e0, min(1.e0, eta1 * (scrch3(i) - eta2) ))
     
! add an addition constraint (Eq. 3.2) to detect contact discontinuities

     tmp1 = abs (p(i+1)   - p(i-1)  ) / min (p(i+1),   p(i-1)  )
     tmp2 = abs (rho(i+1) - rho(i-1)) / min (rho(i+1), rho(i-1))
     
     if (game(i)*ak0*tmp2-tmp1 < 0.e0) scrch3(i) = 0.e0

!  scrch3 now contains the contact steepening coefficient

     tmp1 = a(i-1) + 0.5e0 * hy_dela(i-1)
     tmp2 = a(i+1) - 0.5e0 * hy_dela(i+1)
     
     al(i) = al(i) + (tmp1 - al(i)) * scrch3(i)
     ar(i) = ar(i) + (tmp2 - ar(i)) * scrch3(i)
  end do
  
  return
end subroutine detect
  
