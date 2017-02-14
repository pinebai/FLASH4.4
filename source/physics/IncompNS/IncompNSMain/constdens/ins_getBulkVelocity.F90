!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_getBulkVelocity
!!
!! NAME
!!
!!  ins_getBulkVelocity
!!
!! SYNOPSIS
!!
!!  call ins_getBulkVelocity(real(out) :: velb,
!!                           integer(in) :: myaxis)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   velb : 
!!
!!   myaxis : 
!!
!! AUTOGENROBODOC
!!
!!
!!***





subroutine ins_getBulkVelocity(velB,myaxis)

  use Grid_interface, only : Grid_getBlkBoundBox,     &
                             Grid_getBlkCenterCoords, &
                             Grid_getDeltas,          &
                             Grid_getBlkIndexLimits,  &
                             Grid_getListOfBlocks,    &
                             Grid_getBlkPtr,          &
                             Grid_releaseBlkPtr
  use Driver_interface, ONLY : Driver_abortFlash,     &
                               Driver_getSimTime,     &
                               Driver_getDt

  use IncompNS_data, only : ins_area_solids, ins_meshComm, ins_meshMe
  use IncompNS_data, only : ins_globalDomain

  implicit none
#include "Flash.h"
#include "constants.h"
#include "IncompNS.h"
  include "Flash_mpi.h"

  !! ---- Argument List ----------------------------------
  real, intent(out) :: velB
  integer, intent(in) :: myaxis
  !! -----------------------------------------------------

  ! Local Variables:
  real, parameter :: eps =1.e-12

  integer :: blockCount
  integer, dimension(MAXBLOCKS) :: blockList

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(2,MDIM) :: boundBox
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  integer :: lb,blockID,ierr,i,j,k

  real :: coord(MDIM),bsize(MDIM),del(MDIM),dx,dy,dz,dxdy,dydz,dxdz

  real :: area_tot
  real :: velB_proc,area_proc
  real :: globalKmin
  real  :: simTime,dt

  real, parameter :: tau = 0.025

  area_tot  = 0.
  velB      = 0.
  velB_proc = 0.
  area_proc = 0.

  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Add velocitites and cell areas in blocks of each processor
  select case(myaxis)
  case (IAXIS)

     call Driver_abortFlash("ins_getBulkVelocity: Bulk Velocity in X dir not coded.")

  case (JAXIS)

     call Driver_abortFlash("ins_getBulkVelocity: Bulk Velocity in Y dir not coded.")

#if NDIM == MDIM
  case (KAXIS) 
    do lb = 1,blockCount

       blockID = blockList(lb)

       ! Get blocks coord and bsize
       ! Bounding box:
       call Grid_getBlkBoundBox(blockId,boundBox)
       bsize(1:NDIM) = boundBox(2,1:NDIM) - boundBox(1,1:NDIM)
       call Grid_getBlkCenterCoords(blockId,coord)

       ! Get blocks dx, dy ,dz:
       call Grid_getDeltas(blockID,del)

       globalKmin = ins_globalDomain(LOW,KAXIS)
       if (abs(coord(KAXIS)-0.5*bsize(KAXIS)-globalKmin) .lt. eps*del(KAXIS)) then

          dx = del(IAXIS)
          dy = del(JAXIS)
          !dz = del(KAXIS)
          dxdy = dx*dy
          !dydz = dy*dz
          !dxdz = dx*dz

          ! Get Blocks internal limits indexes:
          call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

          ! Get blocks BCs:
          ! call Grid_getBlkBC(blockID,faces,onBoundary)

          ! KAXIS:
          call Grid_getBlkPtr(blockID,facezData,FACEZ)
          ! Low Z Boundary:
          ! Do The Sum:
          do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
             do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
                velB_proc = velB_proc + facezData(VELC_FACE_VAR,i,j,NGUARD+1)*dxdy
                area_proc = area_proc + dxdy
             enddo
          enddo
          call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

       endif

    enddo
#endif

  end select

  ! Reduce areas and fluxes
  ! Areas
  call mpi_allreduce ( area_proc, area_tot, 1, FLASH_REAL, &
                       MPI_SUM, ins_meshComm, ierr )
  call mpi_allreduce ( velB_proc, velB, 1, FLASH_REAL, &
                       MPI_SUM, ins_meshComm, ierr )

#ifdef EXPONENTIAL_WBREF_RAMP
  call Driver_getSimTime(simTime)
  call Driver_getDt(dt)
  area_tot = area_tot - (1.-exp(-(simTime-dt)/tau))*ins_area_solids
#else
  area_tot = area_tot - ins_area_solids
#endif

  if (ins_meshMe .eq. MASTER_PE) write(*,*) 'Area_tot=',area_tot
 
  ! Bulk velocity
  velB = velB/area_tot

  return

end subroutine ins_getBulkVelocity

