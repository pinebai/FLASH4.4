!!****if* source/Simulation/SimulationMain/ShafranovShock/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution.
!!
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize.
!!  tcurr   -        current simulation time.
!!
!! SIDE EFFECTS
!!
!!  The analytical solution is computed and stored in an appropriate slot
!!  (or slots) in the solution vector, UNK.
!!
!!
!!  Reference: Shafranov.The structure of shock waves in a plasma, Soviet Phys. JETP, 1957.
!!
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  use Simulation_data, ONLY: sim_x, sim_velx, sim_tele, sim_rho, &
                             sim_tion, sim_pele, sim_pion,       &
                             sim_sele, sim_DataPoints, sim_ShockSpeed, sim_meshMe
  use Grid_interface,  ONLY: Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_putPointData
  use Eos_interface,   ONLY: Eos_wrapped
  use ut_interpolationInterface, ONLY: ut_hunt, ut_polint


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)

  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

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
  if (sim_meshMe == MASTER_PE) then
     
     
1    format (1X, 1P, 4(A7, E13.7, :, 1X))
2    format (1X, 1P, 2(A7, E13.7, 1X), A7, I13)
     
  endif

  !Do some test hunting
  !.. set the number of points in the interpolation
  !.. 2 points = linear, 3 points = quadratic, 4 points = cubic and so on
  !.. kat is a first guess at where in the table we are

  op  = 2
  kat = 1
  
  
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


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        
        
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           xx = xCenter(i) - (sim_ShockSpeed*tcurr)
  
           rhoZone = sim_rho(1)

              !! Hunt for xx from given data
           call ut_hunt(sim_x,sim_DataPoints,xx,kat)

           if (kat .ge. 1) then 
              if (kat /= sim_DataPoints) then 
                  call ut_polint(sim_x(kat),sim_rho(kat),op,xx,rhoZone,dyi)                      
              else 
                  rhoZone = sim_rho(sim_DataPoints)
              endif
           end if

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           call Grid_putPointData(blockId, CENTER, RHOA_VAR, EXTERIOR, axis, rhoZone)

        enddo
     enddo
  enddo 


!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

  deallocate(eosData)

 
  return
end subroutine Simulation_computeAnalytical










