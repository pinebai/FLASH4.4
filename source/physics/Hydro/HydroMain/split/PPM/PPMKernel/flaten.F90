!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/flaten
!!
!! NAME
!! 
!!  flaten
!!
!! SYNOPSIS
!!
!! 
!!  call flaten(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells, 
!!              real(IN)    :: u(numCells), 
!!              real(IN)    :: p(numCells),  
!!              real(OUT)   :: flatn(numCells), 
!!              real(OUT)   :: flatn1(numCells))
!!
!! 
!! DESCRIPTION
!!  
!!  Flaten zone structure in regions where shocks are too thin.
!!  This version of subroutine FLATEN only uses the simplest form 
!!  of dissipation as described in the appendix of Colella and Woodward 
!!  (JCP, 54 (1984), 174). Therefore the only constants required are 
!!  omg1, omg2 and epsiln, which are read in.
!!
!!  The "standard" values of the constants are:                  
!!                                                                     
!!          epsiln = 0.33                                             
!!                                                                     
!!          omg1   = 0.75                                             
!!          omg2   = 10.0                                             
!!          sig1   = 0.50                                             
!!          sig2   = 1.00                                             
!!          ak1    = 2.00                                             
!!          ak2    = 0.01                                             
!!                                                                     
!!          wig1   = 2.00                                             
!!          wig2   = 0.00           for 1-d                           
!!                   0.10           for 2-d                           
!!          wig3   = 0.3333 - wig2                                    
!!                                                                     
!!
!!
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells : 
!! u :
!! p :
!! flatn :
!! flatn1 :
!!
!! 
!!
!!***

subroutine flaten(numIntCells,numCells, &
                  u, p,  flatn, flatn1)

  use Hydro_data, ONLY: hy_epsiln, hy_omg1, hy_omg2,&
                        hy_igodu,hy_dp,hy_du,hy_smallu
                        
  
  implicit none

  integer, INTENT(IN):: numIntCells, numCells
  real,  INTENT(IN), DIMENSION(numCells) :: u, p
  real,  INTENT(OUT), DIMENSION(numCells) :: flatn, flatn1


  real,dimension(numCells) :: scrch1,scrch2,scrch3, scrch4
  integer :: i, numIntCells5, numIntCells6, numIntCells7, numIntCells8
  real ::  utest, dutest, dp2, dptest, ptest, dpp, ftilde_up
  
  numIntCells5 = numIntCells + 5
  numIntCells6 = numIntCells + 6
  numIntCells7 = numIntCells + 7
  numIntCells8 = numIntCells + 8
  
  do i = 1, numIntCells8
     flatn (i) = 0.e0
     flatn1(i) = 1.e0
  end do

  do i = 2, numIntCells7

! compute the w_j parameter in Eq. A.1 in Colella & Woodward.  w_j is equal to
! 1 if the jth zone is inside a pressure and velocity jump in the sweep direction,
! in a manner consistent with a shock; storage for this in hy_shockd removed - KW

     hy_dp(i)      = p(i+1) - p(i-1)
     hy_du(i)      = u(i+1) - u(i-1)
     scrch1(i)  = hy_epsiln * min (p(i+1), p(i-1)) - abs( hy_dp(i) )
     utest      = hy_smallu - abs (hy_du(i))

     if (utest .LT. 0.e0) then
        dutest = hy_du(i)
     else
        dutest = 0.e0
     endif
     
     if (scrch1(i) .LT. 0.e0) then 
        scrch1(i) = 1.e0
     else
        scrch1(i) = 0.e0
     endif
     
     if (hy_du(i) .GE. 0.e0) scrch1(i) = 0.e0
     
     if (dutest .EQ. 0.e0) scrch1(i) = 0.e0
     
  end do
  
  do i = 3, numIntCells6

       ! hy_shockd removed from here, use Hydro_detectShock instead - KW


       ! compute ftilde, using Eq. A.2 in Colella & Woodward

     dp2 = p(i+2) - p(i-2)
     
     if ( abs(dp2) .GT. 0.e0 ) then
        dpp = hy_dp(i) / dp2 - hy_omg1
     else
        dpp = 0.e0
     end if
     
     scrch3(i) = scrch1(i) * max (0.e0, dpp * hy_omg2)
  end do

! select upstream value

  do i = 4, numIntCells5
     if ( hy_dp(i) .LT. 0.e0 ) then
        ftilde_up = scrch3(i+1)
     else if ( hy_dp(i) .EQ. 0.e0 ) then
        ftilde_up = scrch3(i)
     else
        ftilde_up = scrch3(i-1)
     endif

! select the maximum flattening

     flatn(i) = max (scrch3(i), ftilde_up)
     flatn(i) = max (0.e0, min (1.e0, flatn(i)))
     
     ! select Godunov method, if desired
     
     flatn (i) = flatn(i) * (1.e0 - hy_igodu) + hy_igodu
     
     flatn1(i) = 1.e0 - flatn(i)
  end do
  
  return
end subroutine flaten

