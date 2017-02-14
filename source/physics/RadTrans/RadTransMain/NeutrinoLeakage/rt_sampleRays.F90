!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/rt_sampleRays
!!
!!  NAME 
!!
!!  rt_sampleRays
!!
!!  SYNOPSIS
!!
!!  call rt_sampleRays( integer(IN) :: nblk,
!!                      integer(IN) :: blklst(nblk),
!!
!!  DESCRIPTION 
!!
!!   This routine maps dens, temp, and Ye from the hydro grid to
!!   the leakage rays.  The Allreduce at the end distributes the 
!!   data to all processes.
!!
!!  ARGUMENTS
!!
!!   nblk   : The number of blocks in the list
!!   blklst : The list of blocks on which the solution must be updated
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***

subroutine rt_sampleRays(nblk, blklst)

  use Grid_interface, ONLY : Grid_releaseBlkPtr, Grid_getBlkPtr, &
       Grid_getBlkBoundBox, Grid_getBlkIndexLimits, Grid_getCellCoords
  use RadTrans_data, ONLY : rt_meshMe
  use rt_data, ONLY : rt_meshComm, rt_leakNumRad, rt_leakNumTht, rt_leakNumPhi,rt_leakRadii, &
       rt_leakX, rt_leakY, rt_leakZ, rt_leakArr, rt_arraySmall,rt_leakTheta,rt_leakPhi,&
       rt_arraySize, rt_minTht, rt_maxTht, rt_threadBlockList, rt_threadWithinBlock
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "Flash_mpi.h"
#include "Flash.h"
#include "constants.h"


  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)

  real, pointer :: solnData(:,:,:,:)
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer :: error
  logical :: doBlock, doRay
  integer :: i,j,k,n,blockID
  real :: delX, delY, delZ
  logical :: gcell = .true.
  integer :: indX, indY, indZ, indRad, indTheta, indTheta2
  logical :: threadBlockList
  integer :: nray


#if NDIM > 1
  call rt_remapRays(nblk,blklst)
#endif

  rt_minTht = rt_leakNumTht
  rt_maxTht = 1

  !$omp parallel if(rt_threadWithinBlock) &
  !$omp default(none) &
  !$omp private(n,i,j,k,blockID,bndBox,doBlock,solnData,blkLimitsGC,blkLimits,dimSize,xCenter,&
  !$omp         yCenter,zCenter,delX,delY,delZ,doRay,indX,indY,indZ,nray) &
#ifndef LEAK_STATIC
  !$OMP shared(rt_leakNumTht,rt_leakNumRad,rt_leakNumPhi) &
#endif
  !$omp shared(nblk,blkLst,rt_leakArr,rt_leakX,&
  !$omp        rt_leakY,rt_minTht,rt_maxTht,gcell,rt_leakZ,rt_leakPhi)

  ! First zero-out the leakage arrays
  !$omp workshare
  rt_leakArr = 0.
  !$omp end workshare

  !$omp do schedule(static) 
  do n=1,nblk
     blockID = blklst(n)

     call Grid_getBlkBoundBox(blockID,bndBox)

     ! First check if any ray point falls within this block
     doBlock = .FALSE.
     do k=1,rt_leakNumPhi
        do j=1,rt_leakNumTht
           do i=1,rt_leakNumRad
              if ( (rt_leakX(i,j,k) >= bndBox(LOW,IAXIS) .AND. rt_leakX(i,j,k) < bndBox(HIGH,IAXIS)) &
#if NDIM > 1
                   .AND. (rt_leakY(i,j,k) >= bndBox(LOW,JAXIS) .AND. rt_leakY(i,j,k) < bndBox(HIGH,JAXIS)) &
#if NDIM == 3
                   .AND. (rt_leakZ(i,j,k) >= bndBox(LOW,KAXIS) .AND. rt_leakZ(i,j,k) < bndBox(HIGH,KAXIS)) &
#endif
#endif
                   ) then
                 doBlock = .TRUE.
                 exit
              end if
           end do
        end do
     end do

     if (doBlock) then
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID,solnData)
        dimSize(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
        if (NDIM > 2)then
           allocate(zCenter(dimSize(KAXIS)))
           call Grid_getCellCoords(KAXIS,blockID,&
                CENTER,gcell,zCenter,dimSize(KAXIS))
           delZ = zCenter(2) - zCenter(1)
        end if
        if (NDIM > 1)then
           allocate(yCenter(dimSize(JAXIS)))
           call Grid_getCellCoords(JAXIS,blockID,&
                CENTER,gcell,yCenter,dimSize(JAXIS))
           delY = yCenter(2) - yCenter(1)
        end if
        allocate(xCenter(dimSize(IAXIS)))
        call Grid_getCellCoords(IAXIS,blockID,&
             CENTER,gcell,xCenter,dimSize(IAXIS))
        delX = xCenter(2) - xCenter(1)

        do k=1,rt_leakNumPhi
           do j=1,rt_leakNumTht
              do i=1,rt_leakNumRad
                 nray = rt_leakNumTht*(k-1) + j
                 ! Check if this ray point falls in the current block
                 doRay = .FALSE.
                 if ( (rt_leakX(i,j,k) >= bndBox(LOW,IAXIS) .AND. rt_leakX(i,j,k) < bndBox(HIGH,IAXIS)) &
#if NDIM > 1
                      .AND. (rt_leakY(i,j,k) >= bndBox(LOW,JAXIS) .AND. rt_leakY(i,j,k) < bndBox(HIGH,JAXIS)) &
#if NDIM == 3
                      .AND. (rt_leakZ(i,j,k) >= bndBox(LOW,KAXIS) .AND. rt_leakZ(i,j,k) < bndBox(HIGH,KAXIS)) &
#endif
#endif
                      ) then
                    doRay = .TRUE.
                 end if
                 if (doRay) then
                    indX = floor((rt_leakX(i,j,k) - bndBox(LOW,IAXIS))/delX + blkLimits(LOW,IAXIS))
                    indY = blkLimits(LOW,JAXIS)
                    indZ = blkLimits(LOW,KAXIS)
#if NDIM > 1
                    indY = floor(rt_leakY(i,j,k) - bndBox(LOW,JAXIS))/delY + blkLimits(LOW,JAXIS)
#if NDIM == 3
                    indZ = floor(rt_leakZ(i,j,k) - bndBox(LOW,KAXIS))/delZ + blkLimits(LOW,KAXIS)
#endif
#endif
                    !Fill leakage arrays.  Direct injection, no interpolation or averaging.
                    ! TODO: add interpolation of these variables for best consistency
                    rt_leakArr(1,i,nray) = solnData(DENS_VAR,indX,indY,indZ)
                    rt_leakArr(2,i,nray) = solnData(TEMP_VAR,indX,indY,indZ)
                    rt_leakArr(3,i,nray) = solnData(YE_MSCALAR,indX,indY,indZ)
                 end if
              end do
           end do
        end do

        call Grid_releaseBlkPtr(blockID,solnData)
        deallocate(xCenter)
#if NDIM >1
        deallocate(yCenter)
#if NDIM >2
        deallocate(zCenter)
#endif
#endif
     end if
  end do !blkLst
  !$omp enddo
  !$omp end parallel

  call Timers_start("communication")
  ! Now get all the leakage array data from other ranks
  call MPI_Allreduce (MPI_IN_PLACE, rt_leakArr(:,:,:), rt_arraySmall, FLASH_REAL, MPI_SUM, & 
       &                rt_meshComm, error)
  call Timers_stop("communication")

end subroutine rt_sampleRays
