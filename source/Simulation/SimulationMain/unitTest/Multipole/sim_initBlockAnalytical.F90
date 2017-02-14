!!****if* source/Simulation/SimulationMain/unitTest/Multipole/sim_initBlockAnalytical
!!
!! NAME
!!
!!  sim_initBlockAnalytical
!!
!! SYNOPSIS
!!
!!  sim_initBlockAnalytical (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Provides analytical gravity potential solution per block for the Maclaurin spheroid
!!  problem. Analytical solution given in: Ricker, P. 2006 "A Direct Multigrid Poisson
!!  Solver for Oct-tree adaptive meshes".
!!
!!  References:  Maclaurin, C. 1742, 
!!               Chandrasekhar, S. 1987, Ellipsoidal Figures of Equilibrium
!!
!! ARGUMENTS
!!
!!  blockID : current grid block identification number
!!
!! NOTES
!!
!!  The following geometries can be handled:
!!
!!                     3D cartesian
!!                     3D cylindrical
!!                     2D cylindrical
!!                     2D spherical
!!
!!***

subroutine sim_initBlockAnalytical (blockID)

  use Simulation_data, ONLY:  sim_Newton,             &
                              sim_pi,                 &
                              sim_e,                  &
                              sim_a1,                 &
                              sim_a3,                 &
                              sim_a1inv,              &
                              sim_a3inv,              &
                              sim_xctr,               &
                              sim_yctr,               &
                              sim_zctr,               &
                              sim_density,            &
                              sim_nsubinv,            &
                              sim_nsubzones,          &
                              sim_initGeometry,       &
                              GRID_3DCARTESIAN,       &
                              GRID_3DCYLINDRICAL,     &
                              GRID_2DCYLINDRICAL,     &
                              GRID_1DSPHERICAL,       &
                              GRID_2DSPHERICAL

  use Grid_interface,  ONLY : Grid_getBlkIndexLimits, &
                              Grid_getCellCoords,     &
                              Grid_putPointData,      &
                              Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID

  logical  :: gcell = .true.

  integer  :: sizeX, sizeY, sizeZ
  integer  :: k, kk, i, ii, j, jj
  integer  :: ibeg,iend,jbeg,jend,kbeg,kend

  real     :: a,b,c
  real     :: a1squared, a3squared, a1squaredInv, a3squaredInv
  real     :: AA1ch, AA3ch, Ich, distTerm
  real     :: DD, Mass, phiSum, potentialAnalytical
  real     :: dx, dy, dz, xx, yy, zz, dzz, dyy, dxx
  real     :: xdist, ydist, zdist, rInv, rFunc, r2, z2, Rdist
  real     :: lambda, h, hinv, H0, atan_h, H1, H3
  real     :: theta

  integer, dimension (MDIM)   :: Position
  integer, dimension (2,MDIM) :: blkLimits, blkLimitsGC
  real,    dimension (MDIM)   :: deltas

  real, dimension(:), allocatable :: xLeft, yLeft, zLeft
!
!
!    ...Set some convenience variables.
!
!
  a1squared    = sim_a1 * sim_a1
  a3squared    = sim_a3 * sim_a3
  a1squaredInv = 1.0 / a1squared
  a3squaredInv = 1.0 / a3squared
!
!
!    ...The coefficients: A1                       -> AA1ch
!                         A3                       -> AA3ch
!                         {x,y,z} independent term -> Ich
!                         potential prefactor      -> DD
!                         total mass spheroid      -> Mass
!
!
  if (abs(sim_e) < tiny(0.0)) then                ! sim_e=0, a sphere
     AA1ch = 2.0 / 3.0
     AA3ch = AA1ch
  else if (abs(sim_e - 1.0) < tiny(0.0) ) then    ! sim_e = 1, a disk?
     AA1ch = 0.0
     AA3ch = 2
  else                                            ! general case
     AA1ch = sqrt(1.0-sim_e**2)*asin(sim_e)/sim_e**3 - (1.0-sim_e**2)/sim_e**2
     AA3ch = 2.0/sim_e**2 - 2.0*sqrt(1.0-sim_e**2)*asin(sim_e)/sim_e**3
  end if

  Ich  = 2.0*AA1ch*a1squared + AA3ch*a3squared
  DD   = sim_pi*sim_Newton*sim_density
  Mass = 4.0/3.0*sim_pi*a1squared*sim_a3*sim_density
!
!
!    ...Loop over all cells in current block, get the cell's {x,y,z} and
!       calculate the analytical potential. Break each cell into nsubzones^ndim
!       sub-zones and average the results to get the values for the cell.
!       This prevents a blocky spheroid shell.
!
!
  call Grid_getBlkIndexLimits (blockID,    &
                               blkLimits,  &
                               blkLimitsGC )

  select case (sim_initGeometry)
!
!
!   ...The 3D cartesian case.
!
!         Full consideration of subzones in all 3 dimensions necessary,
!         as square of radius (r2) and square of z coordinate (z2) for
!         evaluating the MacLaurin potential vary with all 3 dimensions.
!
!
  case (GRID_3DCARTESIAN)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1
    sizeZ = blkLimitsGC (HIGH,KAXIS) - blkLimitsGC (LOW,KAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))
    allocate (zLeft (sizeZ))

    call Grid_getCellCoords (KAXIS, blockID, LEFT_EDGE, gcell, zLeft, sizeZ)
    call Grid_getCellCoords (JAXIS, blockID, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, blockID, LEFT_EDGE, gcell, xLeft, sizeX)
    call Grid_getDeltas     (blockID, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)
    dz = deltas(KAXIS)

    kbeg = blkLimitsGC (LOW, KAXIS)
    kend = blkLimitsGC (HIGH,KAXIS)
    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    do k = kbeg,kend
       dzz = dz * sim_nsubinv
       do j = jbeg,jend
          dyy = dy * sim_nsubinv
          do i = ibeg,iend
             dxx = dx * sim_nsubinv
 
             phiSum = 0.0

             do kk = 0,sim_nsubzones-1
                zz    = zLeft (k) + (kk + 0.5) * dzz 
                zdist = zz - sim_zctr
                do jj = 0,sim_nsubzones-1
                   yy    = yLeft (j) + (jj + 0.5) * dyy
                   ydist = yy - sim_yctr
                   do ii = 0,sim_nsubzones-1
                      xx    = xLeft (i) + (ii + 0.5) * dxx
                      xdist = xx - sim_xctr

                      r2    = xdist**2 + ydist**2
                      z2    = zdist**2

                      rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                      if (rFunc <= 1.) then                                   ! inside the spheroid
                          phiSum = phiSum + DD * (Ich - AA1ch*r2 - AA3ch*z2) 
                      else                                                    ! outside the spheroid
                          a      = 1.0
                          b      = a1squared + a3squared - r2 - z2
                          c      = a1squared * a3squared - r2 * a3squared - z2 * a1squared
                          lambda = (-b + sqrt (b*b - 4.0*a*c)) / (2.0*a)
                          Rinv   = 1.0 / sqrt (a3squared + lambda)

                          if (abs(sim_e) < tiny(0.0)) then                    ! case for e=0  outside the sphere
                              distTerm = 1.0 - (1.0/3.0)*Rinv**2*(r2 + z2)
                          else                                                ! general case outside the spheroid
                              h        = sim_a1 * sim_e * Rinv
                              hinv     = 1.0/h                 
                              atan_h   = atan(h)
                              H0       = atan_h * hinv
                              H1       = (atan_h - h/(1.0+h**2))*hinv**3
                              H3       = (h - atan_h)*2.0*hinv**3
                              distTerm = H0 - 0.5d0*Rinv**2*(r2*H1 + z2*H3)
                          endif
                          phiSum = phiSum + 1.50*sim_Newton*Mass*Rinv*distTerm
                      endif

                   enddo
                enddo
             enddo

             potentialAnalytical = -phiSum/(sim_nsubzones**3)   

             Position (KAXIS) = k
             Position (JAXIS) = j
             Position (IAXIS) = i

             call Grid_putPointData (blockID, CENTER, APOT_VAR, EXTERIOR, Position, potentialAnalytical)

          end do
       enddo
    enddo

    deallocate(xLeft)
    deallocate(yLeft)
    deallocate(zLeft)
!
!
!   ...The 3D cylindrical case.
!
!         Subzones need to be considered only for the radial and z dimension.
!         The angular dimension does not need to be split into subzones (square
!         of radius r2 and square of z coordinate (z2) for evaluating the MacLaurin
!         potential vary only with radial and z dimension component of 3D cylindrical
!         coordinate system.
!
!
  case (GRID_3DCYLINDRICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, blockID, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, blockID, LEFT_EDGE, gcell, xLeft, sizeX)
    call Grid_getDeltas     (blockID, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    kbeg = blkLimitsGC (LOW, KAXIS)
    kend = blkLimitsGC (HIGH,KAXIS)
    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    do j = jbeg,jend
       dyy = dy * sim_nsubinv           ! z component
       do i = ibeg,iend
          dxx = dx * sim_nsubinv        ! radial component
 
          phiSum = 0.0

          do jj = 0,sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             zdist = yy - sim_yctr
             do ii = 0,sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr

                r2 = Rdist ** 2
                z2 = zdist ** 2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                                   ! inside the spheroid
                    phiSum = phiSum + DD * (Ich - AA1ch*r2 - AA3ch*z2) 
                else                                                    ! outside the spheroid
                    a      = 1.0
                    b      = a1squared + a3squared - r2 - z2
                    c      = a1squared * a3squared - r2 * a3squared - z2 * a1squared
                    lambda = (-b + sqrt (b*b - 4.0*a*c)) / (2.0*a)
                    Rinv   = 1.0 / sqrt (a3squared + lambda)

                    if (abs(sim_e) < tiny(0.0)) then                    ! case for e=0  outside the sphere
                        distTerm = 1.0 - (1.0/3.0)*Rinv**2*(r2 + z2)
                    else                                                ! general case outside the spheroid
                        h        = sim_a1 * sim_e * Rinv
                        hinv     = 1.0/h                 
                        atan_h   = atan(h)
                        H0       = atan_h * hinv
                        H1       = (atan_h - h/(1.0+h**2))*hinv**3
                        H3       = (h - atan_h)*2.0*hinv**3
                        distTerm = H0 - 0.5d0*Rinv**2*(r2*H1 + z2*H3)
                    endif
                    phiSum = phiSum + 1.50*sim_Newton*Mass*Rinv*distTerm
                endif

             enddo
          enddo

          potentialAnalytical = - phiSum / (sim_nsubzones**2)   

          do k = kbeg,kend           ! angular component -> no subzones -> identical potentials

             Position (KAXIS) = k
             Position (JAXIS) = j
             Position (IAXIS) = i

             call Grid_putPointData (blockID, CENTER, APOT_VAR, EXTERIOR, Position, potentialAnalytical)

          end do

       enddo
    enddo

    deallocate(xLeft)
    deallocate(yLeft)
!
!
!   ...The 2D cylindrical case.
!
!
  case (GRID_2DCYLINDRICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, blockID, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, blockID, LEFT_EDGE, gcell, xLeft, sizeX)
    call Grid_getDeltas     (blockID, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    do j = jbeg,jend
       dyy = dy * sim_nsubinv
       do i = ibeg,iend
          dxx = dx * sim_nsubinv

          phiSum = 0.0

          do jj = 0,sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             zdist = yy - sim_yctr
             do ii = 0,sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr

                r2    = Rdist**2
                z2    = zdist**2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                                     ! inside the spheroid
                    phiSum = phiSum + DD * (Ich - AA1ch*r2 - AA3ch*z2) 
                else                                                      ! outside the spheroid
                    a      = 1.0
                    b      = a1squared + a3squared - r2 - z2
                    c      = a1squared * a3squared - r2 * a3squared - z2 * a1squared
                    lambda = (-b + sqrt (b*b - 4.0*a*c)) / (2.0*a)
                    Rinv   = 1.0 / sqrt (a3squared + lambda)

                    if (abs(sim_e) < tiny(0.0)) then                      ! case for e=0  outside the sphere
                        distTerm = 1.0 - (1.0/3.0)*Rinv**2*(r2 + z2)
                    else                                                  ! general case outside the spheroid
                        h        = sim_a1 * sim_e * Rinv
                        hinv     = 1.0/h                 
                        atan_h   = atan(h)
                        H0       = atan_h * hinv
                        H1       = (atan_h - h/(1.0+h**2))*hinv**3
                        H3       = (h - atan_h)*2.0*hinv**3
                        distTerm = H0 - 0.5d0*Rinv**2*(r2*H1 + z2*H3)
                    endif
                    phiSum = phiSum + 1.50*sim_Newton*Mass*Rinv*distTerm
                endif

             enddo
          enddo

          potentialAnalytical = -phiSum/(sim_nsubzones**2)   

          Position (JAXIS) = j
          Position (IAXIS) = i

          call Grid_putPointData(blockID, CENTER, APOT_VAR, EXTERIOR, Position, potentialAnalytical)

       enddo
    enddo

    deallocate(xLeft)
    deallocate(yLeft)
!
!
!   ...The 2D spherical case.
!
!
  case (GRID_2DSPHERICAL)

    sizeX = blkLimitsGC (HIGH,IAXIS) - blkLimitsGC (LOW,IAXIS) + 1
    sizeY = blkLimitsGC (HIGH,JAXIS) - blkLimitsGC (LOW,JAXIS) + 1

    allocate (xLeft (sizeX))
    allocate (yLeft (sizeY))

    call Grid_getCellCoords (JAXIS, blockID, LEFT_EDGE, gcell, yLeft, sizeY)
    call Grid_getCellCoords (IAXIS, blockID, LEFT_EDGE, gcell, xLeft, sizeX)
    call Grid_getDeltas     (blockID, deltas)

    dx = deltas(IAXIS)
    dy = deltas(JAXIS)

    jbeg = blkLimitsGC (LOW, JAXIS)
    jend = blkLimitsGC (HIGH,JAXIS)
    ibeg = blkLimitsGC (LOW, IAXIS)
    iend = blkLimitsGC (HIGH,IAXIS)

    do j = jbeg,jend

       dyy = dy * sim_nsubinv
       do i = ibeg,iend
          dxx = dx * sim_nsubinv

          phiSum = 0.0

          do jj = 0,sim_nsubzones-1
             yy    = yLeft (j) + (jj + 0.5) * dyy
             theta = yy - sim_yctr
             do ii = 0,sim_nsubzones-1
                xx    = xLeft (i) + (ii + 0.5) * dxx
                Rdist = xx - sim_xctr                                      ! Rdist is x^2 + y^2 + z^2 here!

                zdist = Rdist * cos (theta)

                z2    = zdist**2                                           ! this is z^2
                r2    = Rdist**2  - z2                                     ! this is x^2 + y^2

                rFunc = r2 * sim_a1inv * sim_a1inv  +  z2 * sim_a3inv * sim_a3inv

                if (rFunc <= 1.) then                                      ! inside the spheroid
                    phiSum = phiSum + DD * (Ich - AA1ch*r2 - AA3ch*z2) 
                else                                                       ! outside the spheroid
                    a      = 1.0
                    b      = a1squared + a3squared - r2 - z2
                    c      = a1squared * a3squared - r2 * a3squared - z2 * a1squared
                    lambda = (-b + sqrt (b*b - 4.0*a*c)) / (2.0*a)
                    Rinv   = 1.0 / sqrt (a3squared + lambda)

                    if (abs(sim_e) < tiny(0.0)) then                       ! case for e=0  outside the sphere
                        distTerm = 1.0 - (1.0/3.0)*Rinv**2*(r2 + z2)
                    else                                                   ! general case outside the spheroid
                        h        = sim_a1 * sim_e * Rinv
                        hinv     = 1.0/h                 
                        atan_h   = atan(h)
                        H0       = atan_h * hinv
                        H1       = (atan_h - h/(1.0+h**2))*hinv**3
                        H3       = (h - atan_h)*2.0*hinv**3
                        distTerm = H0 - 0.5d0*Rinv**2*(r2*H1 + z2*H3)
                    endif
                    phiSum = phiSum + 1.50*sim_Newton*Mass*Rinv*distTerm
                endif

             enddo
          enddo

          potentialAnalytical = -phiSum/(sim_nsubzones**2)   

          Position (JAXIS) = j
          Position (IAXIS) = i

          call Grid_putPointData(blockID, CENTER, APOT_VAR, EXTERIOR, Position, potentialAnalytical)

       enddo
    enddo

    deallocate(xLeft)
    deallocate(yLeft)

  case (GRID_1DSPHERICAL)

    write (*,*) ' Analytical MacLaurin for 1D spherical not implemented! '

  end select
!
!
!   ...Ready!
!
!
  return
end subroutine sim_initBlockAnalytical
