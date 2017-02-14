!!****if* source/Simulation/SimulationMain/unitTest/Poisson/XYdir_2D/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID, 
!!                       integer(in) :: myPE)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference:
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!  myPE   -           my processor number
!!
!! 
!!
!!***

subroutine Simulation_initBlock(blockID, myPE)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkBoundBox,    &
                             Grid_getDeltas,         &
                             Grid_getBlkCenterCoords

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID, myPE
  !!$ ---------------------------------
 
  integer :: i, j, k
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) ::  blIndSize,blIndSizeGC

  real, dimension(MDIM)  :: coord,bsize
  real ::  boundBox(2,MDIM)
  real, pointer, dimension(:,:,:,:) :: solnData

  real :: Lx, Ly, xi, yi, Phi_ijk, F_ijk

  real :: del(3)

  !----------------------------------------------------------------------
  
  Lx = sim_xMax - sim_xMin
  Ly = sim_yMax - sim_yMin  

  ! Get blocks dx, dy ,dz:
  call Grid_getDeltas(blockID,del)

  ! Get Coord and Bsize for the block:
  ! Bounding box:
  call Grid_getBlkBoundBox(blockId,boundBox)
  bsize(:) = boundBox(2,:) - boundBox(1,:)

  call Grid_getBlkCenterCoords(blockId,coord)

  ! Point to Blocks centered variables:
  call Grid_getBlkPtr(blockID,solnData,CENTER)

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,CENTER)


  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

           xi = coord(IAXIS) - 0.5*bsize(IAXIS) + &
                real(i - NGUARD - 1)*del(IAXIS) + 0.5*del(IAXIS)

           yi = coord(JAXIS) - 0.5*bsize(JAXIS) + &
                real(j - NGUARD - 1)*del(JAXIS) + 0.5*del(JAXIS)

           Phi_ijk = sin(2.*PI*xi*pfb_waven_x/Lx + pfb_alpha_x) * &
                     sin(2.*PI*yi*pfb_waven_y/Ly + pfb_alpha_y) 
           
           F_ijk  = -4.*PI**2 * ( (pfb_waven_x/Lx)**2. + (pfb_waven_y/Ly)**2.) * Phi_ijk    
           
           solnData(VANL_VAR,i,j,k) = Phi_ijk

           solnData(VSRC_VAR,i,j,k) = F_ijk

        enddo
     enddo
  enddo

  ! set values for u,v velocities and pressure
  solnData(VERR_VAR,:,:,:) = 0.0
  solnData(VPHI_VAR,:,:,:) = 0.0

  ! Release pointer
  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  return

end subroutine Simulation_initBlock
