!!****if* source/Simulation/SimulationMain/SodSpherical/Simulation_initBlock
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
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up a
!!  Sod-like problem in spherical coordinates to test whether a planar
!!  shock in spherical coordinates stays planar.  This effectively tests the 
!!  fictitious forces in force().  If the forces are setup right, then the planar
!!  shock should stay planar.
!!
!!  Right now, this is setup to do the problem in 2-d spherical coordinates.  
!!  sim_idir = 1 is the x-direction.  sim_idir = 2 is the z-direction.  
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density on the the left side of the interface
!!  sim_rhoRight   Density on the right side of the interface
!!
!!  sim_pLeft      Pressure on the left side of the interface
!!  sim_pRight     Pressure on the righ side of the interface
!!
!!  sim_shockpos   Point of intersection between the shock plane and the x-axis
!!  sim_idir       the direction along which to propagate the shock.
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_shockpos, sim_idir, &    
     &  sim_rhoLeft,  sim_pLeft, sim_rhoRight, sim_pRight, &
     &  sim_smallX, sim_gamma, sim_hydroCfl
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData, &
    Grid_getBlkBC
  use Eos_interface, ONLY: Eos_wrapped


  implicit none

#include "constants.h"
#include "Flash.h"

  integer :: ii, jj, kk

  
  integer, intent(in) :: blockID
  

  real, allocatable, dimension(:) :: yl, yr

  integer, parameter :: nsub = 5
  integer :: nsubj

  logical :: high_state
      
  real :: x_zone, y_zone, z_zone

  real :: rl_zone, r_zone, rr_zone
  real :: thetal_zone, theta_zone, thetar_zone
  real :: phi_zone
  real :: rl_sub, rr_sub, rcc_sub, drr_sub
  real :: tl_sub, tr_sub, tcc_sub, dtt_sub

  real :: vol_zone, vol_sub

  real :: rho_sum, p_sum

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(2,MDIM) :: bcs
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       enerZone, ekinZone
  real :: cfl
  
  logical :: gcell = .true.

  
  
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(yl(sizeY))
  allocate(yr(sizeY))
  allocate(zCoord(sizeZ))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  yl(:) = 0.0
  yr(:) = 0.0
  
  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) then
      call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
      call Grid_getCellCoords&
                      (JAXIS, blockId, LEFT_EDGE,gcell, yl, sizeY)
      call Grid_getCellCoords&
                      (JAXIS, blockId, RIGHT_EDGE,gcell, yr, sizeY)
      nsubj = nsub
  else
     blkLimits(LOW,JAXIS) = 1
     blkLimits(HIGH,JAXIS) = 1
     yCoord(1) = PI / 2
     yl(1) = yCoord(1) - 1.0e-10
     yr(1) = yCoord(1) + 1.0e-10
     nsubj = 1
  end if

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

  if (NDIM > 1) call Grid_getBlkBC(blockID,bcs)


!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

!-----------------------------------------------------------------------------
! loop over all the zones in the current block and set the temperature,
! density, and thermodynamics variables.
!-----------------------------------------------------------------------------

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     axis(KAXIS) = k
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
! z is the phi coordinate in spherical geometry
     phi_zone = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     lPosn0 = sim_shockpos
!     lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
    
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS) = j

! y is the theta coordinate in spherical geometry
        thetal_zone = yl(j)
        theta_zone = yCoord(j)
        thetar_zone = yr(j)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        ! The position of the shock in the current yz-row.
!        lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        lPosn = lPosn0
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS) = i
           
           ! get the cell center, left, and right positions in x
          
! x is the r coordinate in spherical geometry
           rl_zone = xLeft(i)
           r_zone = xCenter(i)
           rr_zone = xRight(i)

! subzone
           rho_sum = 0.0
           p_sum = 0.0

           vol_zone = 0.0

           drr_sub = (rr_zone - rl_zone)/float(nsub)
           dtt_sub = (thetar_zone - thetal_zone)/float(nsub)

           do jj = 0, nsubj-1
              do ii = 0, nsub-1

                 rl_sub = rl_zone + ii*drr_sub
                 rr_sub = rl_zone + (ii + 1.0)*drr_sub
                 rcc_sub = 0.5*(rl_sub + rr_sub)
                 
                 tl_sub = thetal_zone + jj*dtt_sub
                 tr_sub = thetal_zone + (jj + 1.0)*dtt_sub
                 tcc_sub = 0.5*(tl_sub + tr_sub)

! compute the volume of the subzone -- ignore that 2/3 pi prefactor
                 vol_sub = (rr_sub - rl_sub)* &
                      (rr_sub**2 + rl_sub*rr_sub + rl_sub**2)* &
                      (cos(tl_sub) - cos(tr_sub))
                 
! convert the (r, theta) of the subzone to x and z in the Cartesian geometry
                 x_zone = rcc_sub*sin(tcc_sub)
                 z_zone = rcc_sub*cos(tcc_sub)

! stick the discontinuity in the middle of the box in whatever direction is
! specified by sim_idir
                 high_state = .false.

                 select case (sim_idir)

                 case (1)
                    if (x_zone < lPosn) high_state = .true. 
                    
                 case (2)
                    if (z_zone < lPosn) high_state = .true.

                 case default
                    print *, 'ERROR: sim_idir not valid'
                    call Driver_abortFlash('ERROR: sim_idir not valid')

                 end select

                 if (high_state) then
                    rho_sum = rho_sum + sim_rhoLeft*vol_sub
                    p_sum = p_sum + sim_pLeft*vol_sub
                 else
                    rho_sum = rho_sum + sim_rhoRight*vol_sub
                    p_sum = p_sum + sim_pRight*vol_sub
                 endif

                 vol_zone = vol_zone + vol_sub
              enddo
           enddo

! divide through by the volume of the zone
           rhoZone = rho_sum/vol_zone
           presZone = p_sum/vol_zone


           ! data is put stored one cell at a time with these calls to Grid_putPointData
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone) !used in Eos
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone) !used in Eos
           velxZone = 0.0
           velyZone = 0.0
           velzZone = 0.0

           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, 0.0) !dummy use in Eos_wrapped
!           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, enerZone) !ignored in Eos
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone) !used in Eos_wrapped
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone) !used in Eos_wrapped
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone) !used in Eos_wrapped

           !put in value of default species. Ignored in Gamma Eos.
           if (NSPECIES > 0) then 
             call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                        axis, 1.0e0-(NSPECIES-1)*sim_smallX)
             !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if 



#ifdef UNDEFINED
           call Grid_putPointData(blockId, CENTER, VOLU_VAR, EXTERIOR, axis, 0.0)
#endif
#ifdef TEMP_VAR
!           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, temp_zone) !set in Eos
#endif

#ifdef ENER_VAR
!           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone) !set in Eos
#endif
#ifdef GAME_VAR          
!           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma) !set in Eos
#endif
#ifdef GAMC_VAR
!           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma) !set in Eos
#endif

#ifdef CFL_VAR
           cfl = sim_hydroCfl
#if NDIM > 1
           if (bcs(LOW,IAXIS) == REFLECTING .OR. bcs(LOW,IAXIS) == AXISYMMETRIC) then
              if (i==blkLimits(LOW,IAXIS)) then
                 cfl = min(cfl, 0.275)
              end if
           end if
#endif
           call Grid_putPointData(blockId, CENTER, CFL_VAR, EXTERIOR, axis, cfl)
#endif
        enddo
     enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimits, blockID)

! Cleanup
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(yl)
  deallocate(yr)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock










