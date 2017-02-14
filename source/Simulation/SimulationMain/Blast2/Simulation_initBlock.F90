!!****if* source/Simulation/SimulationMain/Blast2/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer (IN) :: blockId)
!!                       
!!
!!
!! DESCRIPTION
!!  Initial conditions for setting up the Woodward Collela 
!!  two blast-wave problem.
!!  Reference:   Woodward, P. & Colella, P. 1984, J. Comp. Phys., 54, 115
!!
!!
!! ARGUMENTS
!!
!!  blockId -           the number of the block to update
!!
!!
!!***


subroutine Simulation_initBlock(blockId)

  use Simulation_data, ONLY: sim_xCos, sim_yCos, sim_zCos, &
     & sim_posnL, sim_posnR, sim_smallX, sim_gamma, sim_smallP, &
       sim_smallE, sim_smallRho, &
     & sim_rhoLeft, sim_rhoMid, sim_rhoRight, &
     & sim_pLeft, sim_pMid, sim_pRight, &
     & sim_uLeft, sim_uMid, sim_uRight, &
     sim_eMassInUAmu
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData

#include "constants.h"
#include "Flash.h"


  implicit none      


  integer, INTENT(in) :: blockId

! compute the maximum length of a vector in each coordinate direction 
! (including guardcells)

  integer :: i, j, k, n
  
  real :: xx, yy, zz, xxL, xxR
  real :: lPosn0L, lPosnL, lPosn0R, lPosnR
  real :: rhoZone, velXZone, velYZone, velzZone, presZone, & 
          eintZone, enerZone, ekinZone
  real :: Zfree, relA, YeFree
  real :: ionEnerFrac, eleEnerFrac
  real :: ionMassFrac, eleMassFrac, radMassFrac
  
  real, allocatable,dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  logical :: gcell = .true.

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  

  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
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


! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

     zz = zCoord(k) ! get cell center coordinates in  z-direction
     axis(KAXIS) = k
     ! Where along the x-axis the shock intersects the xz-plane at the current z
     lPosn0L = sim_posnL - zz*sim_zCos/sim_xCos
     lPosn0R = sim_posnR - zz*sim_zCos/sim_xCos
     
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        axis(JAXIS) = j

        yy = yCoord(j)! get cell center coordinates in the y-direction

        ! The position of the shock in the current yz-row.

        lPosnL = lPosn0L - yy*sim_yCos/sim_xCos
        lPosnR = lPosn0R - yy*sim_yCos/sim_xCos


        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           axis(IAXIS) = i
           
           xx  = xCenter(i) ! get cell center, left, and right positions in x
           
           xxL = xLeft(i)
           xxR = xRight(i)
           
           ! Initialize cells in the left part of the grid.
           
           if (xxR < lPosnL) then
              
              rhoZone = sim_rhoLeft
              presZone = sim_pLeft
              
              velXZone = sim_uLeft * sim_xCos
              velYZone = sim_uLeft * sim_yCos
              velzZone = sim_uLeft * sim_zCos
              
              ! Initialize cells which straddle the left shock. 
              ! Treat them as though 
              ! 1/2 of the cell lay to the left and 1/2 lay to the right.
              
           elseif ((xxL < lPosnL) .and. (xxR > lPosnL)) then
              
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoMid)
              presZone = 0.5 * (sim_pLeft+sim_pMid)
              
              velXZone = 0.5 * (sim_uLeft+sim_uMid) * sim_xCos
              velYZone = 0.5 * (sim_uLeft+sim_uMid) * sim_yCos
              velzZone = 0.5 * (sim_uLeft+sim_uMid) * sim_zCos
              
              ! Initialize cells in the middle part of the grid.
              
           elseif ((xxL > lPosnL) .and. (xxR < lPosnR)) then
              
              rhoZone = sim_rhoMid
              presZone = sim_pMid
              
              velXZone = sim_uMid * sim_xCos
              velYZone = sim_uMid * sim_yCos
              velzZone = sim_uMid * sim_zCos
              
              ! Initialize cells which straddle the right shock.
              
           elseif ((xxL < lPosnR) .and. (xxR > lPosnR)) then
              
              rhoZone = 0.5 * (sim_rhoRight+sim_rhoMid)
              presZone = 0.5 * (sim_pRight+sim_pMid)
              
              velXZone = 0.5 * (sim_uRight+sim_uMid) * sim_xCos
              velYZone = 0.5 * (sim_uRight+sim_uMid) * sim_yCos
              velzZone = 0.5 * (sim_uRight+sim_uMid) * sim_zCos
              
              ! Initialize cells in the right part of the grid.
              
           else
              
              rhoZone = sim_rhoRight
              presZone = sim_pRight
              
              velXZone = sim_uRight * sim_xCos
              velYZone = sim_uRight * sim_yCos
              velzZone = sim_uRight * sim_zCos
              
           endif
           

#if NSPECIES > 0
           !put default species
           call Grid_putPointData(blockID, CENTER, SPECIES_END, EXTERIOR, &
                axis, 1.0e0-(NSPECIES-1)*sim_smallX)


           !if there is only 1 species this loop will not execute
           do n = SPECIES_BEGIN, SPECIES_END-1

              call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                   axis, sim_smallX)



           enddo
#endif
           
           ekinZone = 0.5 * (velXZone**2 + velYZone**2 + velzZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           
!!$           print*,'pres,eint,enerZone:',presZone,eintZone,enerZone

           Zfree = 6
           relA = 12
           YeFree = 6/relA
           ionEnerFrac = 1.0   / (Zfree + 1.0)
           eleEnerFrac = Zfree / (Zfree + 1.0)

           eleMassFrac = Zfree*sim_eMassInUAmu / relA
           ionMassFrac = 1.0 - eleMassFrac

           radMassFrac = 0.0

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
#ifdef PION_VAR
           call Grid_putPointData(blockId, CENTER, PION_VAR, EXTERIOR, axis, presZone*ionEnerFrac)   
#endif
#ifdef PELE_VAR
           call Grid_putPointData(blockId, CENTER, PELE_VAR, EXTERIOR, axis, presZone*eleEnerFrac)   
#endif
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
#endif
#ifdef E1_VAR
           call Grid_putPointData(blockId, CENTER, E1_VAR, EXTERIOR, axis, ekinZone*ionMassFrac+eintZone*ionEnerFrac)
#endif
#ifdef E2_VAR
           call Grid_putPointData(blockId, CENTER, E2_VAR, EXTERIOR, axis, ekinZone*eleMassFrac+eintZone*eleEnerFrac)
#endif
#ifdef E3_VAR
           call Grid_putPointData(blockId, CENTER, E3_VAR, EXTERIOR, axis, sim_smallE)   
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)
#endif
#ifdef EION_VAR
           call Grid_putPointData(blockId, CENTER, EION_VAR, EXTERIOR, axis, eintZone*ionEnerFrac)   
#endif
#ifdef EELE_VAR
           call Grid_putPointData(blockId, CENTER, EELE_VAR, EXTERIOR, axis, eintZone*eleEnerFrac)
#endif
#ifdef ERAD_VAR
           call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, sim_smallE)   
#endif
#ifdef GAME_VAR          
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma)
#endif
#ifdef TEMP_VAR
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, 1e-6)
#endif
#ifdef TION_VAR
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, 1e-6)
#endif
#ifdef TELE_VAR
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, 1e-6)
#endif
#ifdef TRAD_VAR
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, 1e-6)
#endif
#ifdef YE_MSCALAR
!           call Grid_putPointData(blockId, CENTER, YE_MSCALAR, EXTERIOR, axis, 0.5)
           call Grid_putPointData(blockId, CENTER, YE_MSCALAR, EXTERIOR, axis, YeFree)
#endif
#ifdef SUMY_MSCALAR
!           call Grid_putPointData(blockId, CENTER, SUMY_MSCALAR, EXTERIOR, axis, 7.292E-02)
           call Grid_putPointData(blockId, CENTER, SUMY_MSCALAR, EXTERIOR, axis, 1.0/relA)
#endif

        enddo
     enddo
  enddo
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)
  
  return
end subroutine Simulation_initBlock

