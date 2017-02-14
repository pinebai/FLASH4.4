!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!  Simulation_initBlock(  integer(in) :: blockID,
!!                         
!!
!! DESCRIPTION   
!!    Initializes fluid data (density, pressure, velocity, etc.) for
!!               a specified block.  This version sets up the Huang & Greengard
!!               Poisson solver test problem (their example 4.3).
!!
!!    Reference:   Huang, J., & Greengard, L. 2000, SIAM J. Sci. Comput.,
!!               21, 1551
!!
!!
!!
!! ARGUMENTS
!!      blockID:     integer(in)      the current block number to be filled
!!      
!!
!!
!!***



subroutine Simulation_initBlock(blockID)

  use Simulation_data, ONLY : sim_sigma,sim_xctr,sim_yctr,sim_zctr,&
                              sim_smlrho,Npeak, sim_subSample
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putRowData, Grid_getDeltas

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(IN)  :: blockID
  
  real,allocatable, dimension(:) :: xlCoord, ylCoord, zlCoord, &
       xrCoord, yrCoord, zrCoord
  
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real,dimension(MDIM) :: deltas
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: startingPos

  logical :: gcell = .true.

  integer       :: i, j, k, n, istat
  integer       :: Nint, ii, jj, kk, h
  real          :: xx, yy, zz
  real          :: xdist, ydist, zdist, dist
  real          :: Nintinv, sum_rho, exponent, Nintinv1
  
 
  real,allocatable,dimension(:) :: gam, rho  

 
  !===============================================================================


  !get the size of the block and allocate the arrays holding the coordinate info
  !this is done on a blk by block basis for compatibility with block sizes that might
  !not be fixed at compile time
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xlCoord(sizeX))
  allocate(xrCoord(sizeX))

  allocate(ylCoord(sizeY))
  allocate(yrCoord(sizeY))

  allocate(zlCoord(sizeZ))
  allocate(zrCoord(sizeZ))

  !initialize coord arrays
  xlCoord = 0.
  xrCoord = 0.
  ylCoord = 0.
  yrCoord = 0.
  zlCoord = 0.
  zrCoord = 0.
  
  
  allocate(rho(sizeX),stat=istat)
  allocate(gam(sizeX),stat=istat)

  gam = 1.5

  call Grid_getDeltas(blockID,deltas)
  if (NDIM > 2) then 
     call Grid_getCellCoords(KAXIS, blockID, LEFT_EDGE, gcell, zlCoord, sizeZ)
     call Grid_getCellCoords(KAXIS, blockID, RIGHT_EDGE, gcell, zrCoord, sizeZ)
  endif
  
  if (NDIM > 1) then
     call Grid_getCellCoords(JAXIS, blockID, LEFT_EDGE, gcell, ylCoord, sizeY)
     call Grid_getCellCoords(JAXIS, blockID, RIGHT_EDGE, gcell, yrCoord, sizeY)
  endif
  
  call Grid_getCellCoords(IAXIS, blockID, LEFT_EDGE, gcell, xlCoord, sizeX)
  call Grid_getCellCoords(IAXIS, blockID, RIGHT_EDGE, gcell, xrCoord, sizeX)
  
  !  Loop over cells in the block.  For each, compute the physical
  !  position of its left and right edge and its center as well as
  !  its physical width.  Then decide whether it is inside the
  !  initial radius or outside and initialize the hydro variables
  !  appropriately.
    
  
  Nint    = sim_subSample  ! Default is 7
  Nintinv = 1./float(Nint)
  Nintinv1= 1./(float(Nint)-1.)
  do k = 1, sizeZ
     do j = 1, sizeY
        do i = 1, sizeX
           
           sum_rho = 0.
           
           do kk = 0, (Nint-1)*K3D
              zz    = zlCoord(k) + kk*Nintinv1*deltas(KAXIS) !!(zrCoord(k)-zlCoord(k))
              do jj = 0, (Nint-1)*K2D
                 yy    = ylCoord(j) + jj*Nintinv1*deltas(JAXIS) !!(yrCoord(j)-ylCoord(j))
                 do ii = 0, Nint-1
                    xx    = xlCoord(i) + ii*Nintinv1*deltas(IAXIS) !! (xrCoord(i)-xlCoord(i))
                    
                    do h = 1, Npeak
                       xdist = (xx - sim_xctr(h))**2
                       ydist = (yy - sim_yctr(h))**2 * K2D
                       zdist = (zz - sim_zctr(h))**2 * K3D
                       dist  = xdist + ydist + zdist
                       exponent = sim_sigma(h) * dist
                       !       Avoid underflow by only computing exponential
                       !       for "small" exponent magnitudes
                       if (exponent < 80.)sum_rho = sum_rho + exp(-exponent)
                    enddo
                    
                 enddo
              enddo
           enddo
           
           rho(i) = max(sum_rho * Nintinv**NDIM, sim_smlrho)
           
        enddo
        
        
        startingPos(IAXIS) = blkLimitsGC (LOW, IAXIS)
        startingPos(JAXIS) = j
        startingPos(KAXIS) = k
        
        !!call Grid_putRowData(blockID, CENTER, GAME_VAR, EXTERIOR, IAXIS, startingPos, gam, sizeX)
        !!call Grid_putRowData(blockID, CENTER, GAMC_VAR, EXTERIOR, IAXIS, startingPos, gam, sizeX)
        call Grid_putRowData(blockID, CENTER, DENS_VAR, EXTERIOR, IAXIS, startingPos, rho, sizeX)
        
     enddo
  enddo

  deallocate(xlCoord)
  deallocate(xrCoord)

  deallocate(ylCoord)
  deallocate(yrCoord)

  deallocate(zlCoord)
  deallocate(zrCoord)
  
  deallocate(rho)
  deallocate(gam)
  
  return

end subroutine Simulation_initBlock
