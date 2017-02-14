!!****if* source/Simulation/SimulationMain/CCSN/Simulation_initBlock
!!
!! NAME
!! 
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID)
!!                       
!!
!! DESCRIPTION
!!
!! Initial conditions for Core Collapse SN problem
!!
!! ARGUMENTS
!!
!!  blockID - my block number
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine Simulation_initBlock(blockId)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, Grid_getBlkPtr, &
       Grid_releaseBlkPtr, Grid_getGeometry, Grid_renormAbundance
  use Eos_interface, ONLY : Eos_getAbarZbar
  use Multispecies_interface, ONLY : Multispecies_getSumInv

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"

  integer, intent(IN) :: blockID


  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:) :: xCenter, xLeft, xRight
  real, allocatable, dimension(:) :: yCenter, yLeft, yRight
  real, allocatable, dimension(:) :: zCenter, zLeft, zRight
  real, dimension(MDIM) :: delta
  real :: dx, dy, dz

  real, dimension(NSPECIES) :: massFrac

  integer,dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  integer :: iSize, jSize, kSize
  integer :: iSizeGC, jSizeGC, kSizeGC
  integer :: ilo, ihi

  integer :: ivar, meshGeom
  integer  ::  jlo, jhi
  integer  ::  n

  real     ::  dens, temp, pres, He4, C12, N14, O16,velx,vely
  real     ::  Ne20
  real     ::  Fe56, Neut, H1,tot

  real :: angle, rho_wind
  real :: dxx_sub, dyy_sub, dzz_sub
  integer :: i, j, k, ii, jj, kk
  integer :: istat

  real :: vol, sum, suminv
  real :: rcc, r_xy, rcc_sub, v_xy
  real :: var_interp, var_sum, vtot
  real :: xcc_sub, ycc_sub, zcc_sub
  real :: sumY,Ye,Ye0

  real :: temp0,xEner,pres0,xEntr,xdedt,xdpderho,xMuNu,xXp,xXn,xXa,xXh
  real :: dTdp, temp1
  integer :: iter
  real :: sign
  real :: tht, tmp, rad, nodes, veltht
  real :: deltaRho, cs, phi

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  iSizeGC = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSizeGC = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSizeGC = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  iSize = blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS)+1
  jSize = blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS)+1
  kSize = blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS)+1

  ilo = blkLimits(LOW,IAXIS)
  ihi = blkLimits(HIGH,IAXIS)

  !! allocate all needed space
  allocate(xCenter(iSizeGC),STAT=istat)
  allocate(xLeft(iSizeGC),STAT=istat)
  allocate(xRight(iSizeGC),STAT=istat)
  allocate(yCenter(jSizeGC),STAT=istat)
  allocate(yLeft(jSizeGC),STAT=istat)
  allocate(yRight(jSizeGC),STAT=istat)
  allocate(zCenter(kSizeGC),STAT=istat)
  allocate(zLeft(kSizeGC),STAT=istat)
  allocate(zRight(kSizeGC),STAT=istat)

  xCenter(:) = 0.e0
  yCenter(:) = 0.e0
  zCenter(:) = 0.e0

  call Grid_getDeltas(blockId, delta)
  dx = delta(IAXIS)
  dy = delta(JAXIS)
  dz = delta(KAXIS)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getGeometry(meshGeom)

  if (NDIM == 1 .AND. meshGeom == SPHERICAL) then

     call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
     call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE, .true.,xLeft, iSizeGC)
     call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,xRight,iSizeGC)

     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)

              ! the interpolation will be done using the parabolic interpolation 
              ! routine
              if (xCenter(i) <= xzn(n1d_total)) then
                 do ivar = 1, NUNK_VARS
                    call parabolic_interp(xzn, model_1d(:,ivar), & 
                         n1d_total, xLeft(i), xCenter(i), xRight(i), var_interp)
                    solnData(ivar,i,j,k) = var_interp
                 end do
              else 
                 call Driver_abortFlash("Beyond the 1D model data")
              endif
              ! multiply initial radial velocity, if desired
              solnData(VELX_VAR,i,j,k) = sim_velMult*solnData(VELX_VAR,i,j,k) 

              sumY = 1.0

#ifdef FLASH_MULTISPECIES
              sum = 0.e0
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,i,j,k) = & 
                      max(sim_smallx, &
                      min(1.e0,solnData(n,i,j,k)))
                 sum = sum + solnData(n,i,j,k)
              enddo
              suminv = 1.e0 / sum
              do n = SPECIES_BEGIN, SPECIES_END
                 solnData(n,i,j,k) =  & 
                      max(sim_smallx, min(1.e0,suminv*&
                      solnData(n,i,j,k)))
              enddo
              call Multispecies_getSumInv(A,sumY,solnData(SPECIES_BEGIN:SPECIES_END,i,j,k))
#endif
#ifdef SUMY_MSCALAR
              solnData(SUMY_MSCALAR,i,j,k) = sumY
#endif
              
           enddo
        enddo
     enddo

  endif

  if (NDIM == 2) then
     !--------------------------------------------------------------------------
     ! create a circular mapping front 
     !--------------------------------------------------------------------------
     if (meshGeom == CARTESIAN .OR. meshGeom == CYLINDRICAL) then 

        call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
        call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE, .true.,xLeft, iSizeGC)
        call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,xRight,iSizeGC)

        call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)
        call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE, .true.,yLeft, jSizeGC)
        call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,.true.,yRight,jSizeGC)

        ! now fill the master arrays
        do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                 ! compute the distance of the current zone from the origin 
                 rcc = sqrt(xCenter(i)**2 + yCenter(j)**2)
                 ! the interpolation will be done using the parabolic interpolation routine
                 if (rcc <= xzn(n1d_total)) then
                    do ivar = 1, NUNK_VARS
                       ! subsample in each zone to get a more accurate zone average -- note,
                       ! this is geometry dependent, so be careful
                       var_sum = 0.0
                       dxx_sub = dx/float(nsub)
                       dyy_sub = dy/float(nsub)
                       do jj = 1, nsub
                          do ii = 1, nsub
                             xcc_sub = xLeft(i) + (ii - 0.5)*dxx_sub
                             ycc_sub = yLeft(j) + (jj - 0.5)*dyy_sub
                             rcc_sub = sqrt(xcc_sub**2 + ycc_sub**2)
                             ! since it is difficult to do the parabolic interpolation, send
                             ! rcc as the left, center, and right positions -- this will do
                             ! a linear interpolation
                             call parabolic_interp(xzn, model_1d(:,ivar), & 
                                  n1d_total, rcc_sub, rcc_sub, rcc_sub, var_interp)
                             ! add the subzone's contribution to entire zone's total -- taking into
                             ! account the geometrical weighting
                             if (meshGeom == CARTESIAN) then
                                var_sum = var_sum + var_interp
                             elseif (meshGeom == CYLINDRICAL) then
                                ! the 'x' coordinate is the cylindrical radial coord, the 'y' coordinate
                                ! is the cylindrical z coord
                                vol = 2.0*PI*2.0*xcc_sub*dxx_sub*dyy_sub
                                var_sum = var_sum + var_interp*vol
                             endif

                          enddo
                       enddo
                       ! divide by the volume of the entire zone to get the average
                       if (meshGeom == CARTESIAN) then
                          var_sum = var_sum / float(nsub*nsub)
                       else if (meshGeom == CYLINDRICAL) then
                          vol = 2.0*PI*2.0*xCenter(i)*dx*dy
                          var_sum = var_sum / vol
                       endif

                       ! fake the velocities -- assume that v_x in the table is v_tot -- it 
                       ! nearly is.  Then compute the angle from xCenter and yCenter and find the 
                       ! x and y compontents of the velocity -- do all of this when velx is 
                       ! read from the table
                       if (ivar .NE. VELY_VAR .AND. ivar .NE. VELX_VAR) then
                          solnData(ivar,i,j,k) = var_sum
                       else
                          ! do both velocities when ivar eq ivelx
                          if (ivar == VELX_VAR) then
                             vtot = var_sum * sim_velMult
                             if (xCenter(i) .NE. 0) then
                                angle = atan(yCenter(j)/xCenter(i))
                             else
                                angle = PI/2.
                             endif
                             solnData(VELX_VAR,i,j,k) = vtot*cos(angle)
                             solnData(VELY_VAR,i,j,k) = vtot*sin(angle)
                          endif

                       endif
                    enddo
                    
                 else 
                    call Driver_abortFlash("Beyond the 1D model data")
                 end if
                 
                 sumY = 1.0
#if NSPECIES > 0
                 sum = 0.e0
                 do n = SPECIES_BEGIN,SPECIES_END
                    solnData(n,i,j,k) = & 
                         max(sim_smallx, &
                         min(1.e0,solnData(n,i,j,k)))
                    sum = sum + solnData(n,i,j,k)
                 enddo
                 suminv = 1.e0 / sum
                 do n = SPECIES_BEGIN, SPECIES_END
                    solnData(n,i,j,k) =  & 
                         max(sim_smallx, min(1.e0,suminv*&
                         solnData(n,i,j,k)))
                 enddo
                 call Eos_getAbarZbar(solnData(:,i,j,k),sumY=sumY)
#endif
#ifdef SUMY_MSCALAR
                 solnData(SUMY_MSCALAR,i,j,k) = sumY
#endif
              enddo
           enddo
        enddo

     else ! Here we may add 2D spherical geometry
        call Driver_abortFlash("incorrect geometry in Simulation_initBlock")
     end if

  else if (NDIM == 3 .and. meshGeom == CARTESIAN) then
     !------------------------------------------------------------------------------
     ! create a spherical mapPIng
     !------------------------------------------------------------------------------

     call Grid_getCellCoords(IAXIS,blockID, CENTER, .true.,xCenter,iSizeGC)
     call Grid_getCellCoords(IAXIS,blockID,LEFT_EDGE, .true.,xLeft, iSizeGC)
     call Grid_getCellCoords(IAXIS,blockID,RIGHT_EDGE,.true.,xRight,iSizeGC)

     call Grid_getCellCoords(JAXIS,blockID, CENTER, .true.,yCenter,jSizeGC)
     call Grid_getCellCoords(JAXIS,blockID,LEFT_EDGE, .true.,yLeft, jSizeGC)
     call Grid_getCellCoords(JAXIS,blockID,RIGHT_EDGE,.true.,yRight,jSizeGC)

     call Grid_getCellCoords(KAXIS,blockID, CENTER, .true.,zCenter,kSizeGC)
     call Grid_getCellCoords(KAXIS,blockID,LEFT_EDGE, .true.,zLeft, kSizeGC)
     call Grid_getCellCoords(KAXIS,blockID,RIGHT_EDGE,.true.,zRight,kSizeGC)

     ! now fill the master arrays
     do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              ! compute the distance of the current zone from the origin 
              rcc = sqrt(xCenter(i)**2 + yCenter(j)**2 + zCenter(k)**2)
              ! the interpolation will be done using the parabolic interpolation 
              ! routine
              do ivar = 1, NUNK_VARS
                 ! subsample in each zone to get a more accurate zone average -- note,
                 ! this is geometry dependent, so be careful
                 var_sum = 0.0
                 dxx_sub = dx/float(nsub)
                 dyy_sub = dy/float(nsub)
                 dzz_sub = dz/float(nsub)

                 do kk = 1, nsub
                    do jj = 1, nsub
                       do ii = 1, nsub
                          xcc_sub = xLeft(i) + (ii - 0.5)*dxx_sub
                          ycc_sub = yLeft(j) + (jj - 0.5)*dyy_sub
                          zcc_sub = zLeft(k) + (kk - 0.5)*dzz_sub

                          rcc_sub = sqrt(xcc_sub**2 + &
                               ycc_sub**2 + &
                               zcc_sub**2)
                          ! since it is difficult to do the parabolic interpolation, send
                          ! rcc as the left, center, and right positions -- this will do
                          ! a linear interpolation
                          call parabolic_interp(xzn, model_1d(:,ivar), & 
                               n1d_total, rcc_sub, rcc_sub, rcc_sub, &
                               var_interp)
                          ! add the subzone's contribution to entire zone's total -- taking into
                          ! account the geometrical weighting
                          var_sum = var_sum + var_interp
                       enddo
                    enddo
                 enddo
                 ! divide by the volume of the entire zone to get the average
                 var_sum = var_sum / float(nsub*nsub*nsub)
                 ! fake the velocities -- assume that v_x in the table is v_tot -- it 
                 ! nearly is.  Then compute the angle from xCenter, yCenter, and zCenter and find the 
                 ! x, y, and z compontents of the velocity -- do all of this one velx is 
                 ! read from the table
                 if (ivar /= VELX_VAR .AND. &
                      ivar /= VELY_VAR .AND. &
                      ivar /= VELZ_VAR) then
                    solnData(ivar,i,j,k) = var_sum
                 else
                    ! do both velocities when ivar eq VELX_VAR
                    if (ivar == VELX_VAR) then
                       vtot = var_sum
                       ! first decompose the velocity into a z component and an 'xy' component
                       r_xy = sqrt(xCenter(i)**2 + yCenter(j)**2)
                       if (r_xy /= 0.0) then
                          angle = atan(zCenter(k)/r_xy)
                       else
                          angle = PI/2.
                       endif

                       solnData(VELZ_VAR,i,j,k) = vtot*sin(angle)*sim_velMult

                       v_xy = vtot*cos(angle)
                       if (xCenter(i) /= 0.0) then
                          angle = atan(yCenter(j)/xCenter(i))
                       else
                          angle = PI/2.0
                       endif
                       sign = xCenter(i)/abs(xCenter(i))
                       solnData(VELX_VAR,i,j,k) = sign*v_xy*cos(angle)*sim_velMult
                       solnData(VELY_VAR,i,j,k) = sign*v_xy*sin(angle)*sim_velMult
                    endif

                 endif
              enddo
              
              sumY = 1.0
#ifdef FLASH_MULTISPECIES
              sum = 0.e0
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,i,j,k) = & 
                      max(sim_smallx, &
                      min(1.e0,solnData(n,i,j,k)))
                 sum = sum + solnData(n,i,j,k)
              enddo
              suminv = 1.e0 / sum
              do n = SPECIES_BEGIN, SPECIES_END
                 solnData(n,i,j,k) =  & 
                      max(sim_smallx, min(1.e0,suminv*&
                      solnData(n,i,j,k)))
              enddo
              call Eos_getAbarZbar(solnData(:,i,j,k),sumY=sumY)
#endif
#ifdef SUMY_MSCALAR
              solnData(SUMY_MSCALAR,i,j,k) = sumY
#endif

           enddo
        enddo
     enddo

  end if

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yLeft)
  deallocate(yRight)
  deallocate(yCenter)
  deallocate(zLeft)
  deallocate(zRight)
  deallocate(zCenter)

  return
end subroutine Simulation_initBlock


subroutine parabolic_interp(x,var,n,y_l,y_c,y_r,var_interp)
!
! Given a vector of coordinates, x, the size, n, and associated function
! values, var, take the zone edges and center (y_l, y_c, y_r), and return
! the parabolic interpolation to this.
!
! x(n)        coordinate values 
! var(n)      function values at x(n)
! n           size of x and var
!
! y_l         coordinate of left edge of the zone
! y_c         coordinate of center of the zone
! y_r         coordinate of right edge of the zone
! 
! var_interp  zone average value of the function in that zone, using 
!             a parabolic structure
!
  implicit none
      
  integer :: n
      
  real :: x(n), var(n)
  real :: y_l, y_c, y_r
      
  real :: var_interp

  real :: var_l, var_c, var_r, err_int

  integer, PARAMETER :: op = 2

  real, PARAMETER :: sixth = 1.e0/6.e0

  integer kat, tmp

! get the function value at the left edge of the zone
real :: entropy, dst, dsd

kat = 0

  if (y_l < x(1)) then  

! the x array is monotonic -- if we are below the minimum in x, then
! just set the index to the first value
     kat = 1
  else
     call ut_hunt(x,n,y_l,kat)
  endif

  kat = max(1, min(kat - op/2 + 1, n - op + 1))
  call ut_polint(x(kat),var(kat),op,y_l,var_l,err_int)

! get the function value at the center of the zone
  call ut_hunt(x,n,y_c,kat)
  kat = max(1, min(kat - op/2 + 1, n - op + 1))
  call ut_polint(x(kat),var(kat),op,y_c,var_c,err_int)

! get the function value at the right edge of the zone
  call ut_hunt(x,n,y_r,kat)
  kat = max(1, min(kat - op/2 + 1, n - op + 1))
  call ut_polint(x(kat),var(kat),op,y_r,var_r,err_int)

! construct the zone averaged value
  var_interp = (var_l + 4.e0*var_c + var_r)*sixth

  return
end subroutine parabolic_interp


subroutine ye_ofRho_init(xDens,outYe)
  
  use Simulation_data, ONLY : sim_rhoOne, sim_rhoTwo, &
       sim_yOne, sim_yTwo, sim_yc, &
       sim_yThree, sim_rhoThree

  implicit none
    
  real, intent(IN) :: xDens 
  real, intent(OUT) :: outYe
  
  real :: xofrho, xofrho2

  xofrho = 2.0*log10(xDens) - log10(sim_rhoTwo) - log10(sim_rhoOne)
  xofrho = xofrho / (log10(sim_rhoTwo) - log10(sim_rhoOne))

  xofrho = max(-1.0,min(1.0,xofrho))

  xofrho2 = log10(xDens) - log10(sim_rhoThree)
  xofrho2 = xofrho2 / (log10(sim_rhoThree) - log10(sim_rhoTwo))
  xofrho2 = max(-1.,min(1.,xofrho2))

  outYe = 0.5*(sim_yTwo + sim_yOne) + 0.5*xofrho*(sim_yTwo - sim_yOne) &
          + sim_yc*(1. - abs(xofrho) &
          + 4.*abs(xofrho)*(abs(xofrho)-0.5)*(abs(xofrho) - 1.)) &
          + max(0.,(sim_yThree-sim_yTwo)*xofrho2)

  return
end subroutine ye_ofRho_init
