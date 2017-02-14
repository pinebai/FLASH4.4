!!****if* source/Simulation/SimulationMain/ShafranovShock/Simulation_initBlock
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
!!  a specified block.
!!
!!  Reference: Shafranov.The structure of shock waves in a plasma, Soviet Phys. JETP, 1957.
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY: sim_x, sim_velx, sim_tele, sim_rho, &
                             sim_tion, sim_pele, sim_pion,       &
                             sim_sele, sim_DataPoints, sim_smallX,sim_meshMe 
  use Simulation_interface, ONLY : Simulation_computeAnalytical
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos_wrapped
  use ut_interpolationInterface, ONLY: ut_hunt, ut_polint


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax,vecLen

  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis

  
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone
  
  logical :: gcell = .true.
  
  ! For hunting
  integer ::  op, kat
  real    ::  dist, dyi

  real    :: TionZone, TEleZone, PeleZone, PionZone

  real, dimension(:), allocatable :: eosData

  vecLen = 1
  allocate(eosData(vecLen * EOS_NUM))
  
  ! dump some output to stdout listing the paramters
!!$  if (sim_meshMe == MASTER_PE) then
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

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.
  
  op  = 2
  kat = 1
  
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        ! The position of the shock in the current yz-row.
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           ! get the cell center, left, and right positions in x         
           
           xx = xCenter(i)
 
           !! Hunt for xx from given data
           call ut_hunt(sim_x,sim_DataPoints,xx,kat)

           
           if (kat == sim_DataPoints) then
              xx = sim_x(kat)
              
              kat = sim_DataPoints - 1
              
           end if

           if (.TRUE.) then
              !! density
              call ut_polint(sim_x(kat),sim_rho(kat),op,xx,rhoZone,dyi)                     
              
              
              
           !! velocity X
              call ut_polint(sim_x(kat),sim_velx(kat),op,xx,velxZone,dyi)
              
              !! ele pressure
              call ut_polint(sim_x(kat),sim_pele(kat),op,xx,PeleZone,dyi)
           
              !! ion pressure
              call ut_polint(sim_x(kat),sim_pion(kat),op,xx,PionZone,dyi)
           
              !! total pressure         
              presZone = PeleZone + PionZone
           
              velyZone = 0.
              velzZone = 0.
           
              !! ion temp
              call ut_polint(sim_x(kat),sim_tion(kat),op,xx,TionZone,dyi)
           
              !! ele pressure
              call ut_polint(sim_x(kat),sim_tele(kat),op,xx,TEleZone,dyi)          
              
           endif
           
           
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k      
           
           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, TEleZone)
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, TionZone)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, TEleZone)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, 1.0e-20)
           
           !put in value of default species
           if (NSPECIES > 0) then
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)
              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if
           
        enddo
     enddo
  enddo

  
#ifdef EELE_VAR
  call Eos_wrapped(MODE_DENS_TEMP_GATHER,blkLimits,blockId)
#endif

!!$  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!$     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!!$        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!!$           axis(IAXIS) = i
!!$           axis(JAXIS) = j
!!$           axis(KAXIS) = k
!!$#ifdef ERAD_VAR
!!$           call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, 0.0  )   
!!$#endif
!!$#ifdef E3_VAR
!!$           call Grid_putPointData(blockId, CENTER, E3_VAR,   EXTERIOR, axis, 0.0  )   
!!$#endif
!!$
!!$#ifdef PRAD_VAR
!!$           call Grid_putPointData(blockId, CENTER, PRAD_VAR, EXTERIOR, axis, 0.0  )   
!!$#endif
!!$#ifdef TRAD_VAR
!!$           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, 0.0  )   
!!$#endif
!!$        enddo
!!$     enddo
!!$  enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)
  deallocate(eosData)

  call Simulation_computeAnalytical(blockId=blockID,tcurr=0.0) 
 
  return
end subroutine Simulation_initBlock










