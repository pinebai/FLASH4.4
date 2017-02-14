!!****if* source/Simulation/SimulationMain/RHD_Sod/Simulation_initBlock
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
!!  a specified block.  This version sets up the RHD Sod shock-tube
!!  problem.
!!
!!  Reference: "Assessment of a high-resolution central scheme for
!!              the solution of the relativistic hydrodynamics equations"
!!
!!              Lucas-Serrano, A.; Font, J. A.; Ibanez, J. M.; Marti, J. M.
!!              A&A, v. 428, p. 703 (2004)
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos, &    
                             sim_rhoLeft,  sim_pLeft,  sim_uLeft, &
                             sim_rhoRight, sim_pRight, sim_uRight, &
                             sim_vLeft,    sim_vRight, sim_wLeft,sim_wRight,&
                             sim_smallX, sim_gamma, sim_smallP

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords, &
                             Grid_putPointData


  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k
  real :: xx, yy,  zz, xxL, xxR
  real :: lPosn0, lPosn
  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
          enerZone, enthalpyZone, lorentzZone, eintZone
  
  logical :: gcell = .true.
  integer :: istat 

  
  ! dump some output to stdout listing the paramters
!!$  if (sim_meshME == MASTER_PE) then
!!$     
!!$     
!!$1    format (1X, 1P, 4(A7, E13.7, :, 1X))
!!$2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
!!$     
!!$  endif
  
  
  ! get the integer index information for the current block
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX),stat=istat)
  allocate(xRight(sizeX),stat=istat)
  allocate(xCenter(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        ! The position of the shock in the current yz-row.
        lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
           xxL = xLeft(i)
           xxR = xRight(i)
          
           
           ! initialize cells to the left of the initial shock.
           if (xxR <= lPosn) then
              
              rhoZone  = sim_rhoLeft
              presZone = sim_pLeft

              ! Note: In order to simulate rotated shock tube problem
              !       uncomment out shock direction angles
              velxZone = sim_uLeft !* sim_xCos
              velyZone = sim_vLeft !* sim_yCos
              velzZone = sim_wLeft !* sim_zCos
              
              ! initialize cells which straddle the shock.  Treat them as though 1/2 of 
              ! the cell lay to the left and 1/2 lay to the right.
           elseif ((xxL < lPosn) .and. (xxR > lPosn)) then
              
              rhoZone  = 0.5 * (sim_rhoLeft+sim_rhoRight)
              presZone = 0.5 * (sim_pLeft+sim_pRight)
              
              velxZone = 0.5 *(sim_uLeft+sim_uRight) !* sim_xCos
              velyZone = 0.5 *(sim_vLeft+sim_uRight) !* sim_yCos
              velzZone = 0.5 *(sim_wLeft+sim_uRight) !* sim_zCos
              
              ! initialize cells to the right of the initial shock.
           else
              
              rhoZone  = sim_rhoRight
              presZone = sim_pRight

              velxZone = sim_uRight !* sim_xCos
              velyZone = sim_vRight !* sim_yCos
              velzZone = sim_wRight !* sim_zCos
           endif

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

#if NSPECIES > 0
           !put in value of default species
           call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                axis, 1.0e0-(NSPECIES-1)*sim_smallX)

           !if there is only 1 species, this loop will not execute
           do n = SPECIES_BEGIN+1,SPECIES_END
              call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                   axis, sim_smallX)
           enddo
#endif

           ! Compute Lorentz factor
           lorentzZone = 1./(1.-(velxZone**2.0 + velyZone**2.0 + velzZone**2.0))
           lorentzZone = sqrt(lorentzZone)

           enthalpyZone = 1.0 + sim_gamma/(sim_gamma-1.0)*presZone/rhoZone

           enerZone = rhoZone*enthalpyZone*lorentzZone**2.0 - presZone

!!$           enerZone = presZone/rhoZone/(sim_gamma - 1.0)

           eintZone = presZone/rhoZone/(sim_gamma - 1.0)

           ! Store the variables in the current zone via the database put methods
           ! data is put stored one cell at a time with this call to Grid_putData           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)
#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif
        enddo
     enddo
  enddo

  !! Cleanup!  Must deallocate arrays
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock










