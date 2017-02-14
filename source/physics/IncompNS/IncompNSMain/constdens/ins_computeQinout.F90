!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_computeQinout
!!
!!
!! NAME
!!
!!  ins_computeQinout
!!
!!
!! SYNOPSIS
!!
!!  ins_computeQinout(integer(IN) :: blockCount,
!!                    integer(IN) :: blockList(blockCount)
!!                    logical(IN) :: inou_flg
!!                    real(OUT)   :: Qinout)
!!
!!
!! DESCRIPTION
!!
!! Computes total flow into or out of domain boundaries. This routine is part 
!! of a global mass balance strategy.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  inou_flg   - Flag to set Between Volume Flow in (.true.) or Flow out (.False).
!!               Flow out given by NEUMANN_INS, OUTFLOW_INS (+ sign out of domain) 
!!               the rest of BC are computed as Qin (+ sign into the domain).
!!  Qinout     - Flow in or out of the domain.
!!
!!***

subroutine ins_computeQinout( blockCount, blockList, inou_flg, Qinout)

#include "Flash.h"
  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use IncompNS_data, ONLY : ins_meshComm, ins_domainBC

  implicit none

#include "constants.h"
  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
  logical, INTENT(IN) :: inou_flg
  real,    INTENT(OUT) :: Qinout
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  integer :: lb,blockID,ierr,i,j,k

  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer :: nxc,nyc,nzc

  real :: Qinaux, del(MDIM), dx,dy,dz,dxdy,dydz,dxdz

!=============================================================================

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1

  Qinout = 0.
  Qinaux = 0.


  ! Detect if the problem is an outflow problem to proceed with flow computation.
  if (any(ins_domainBC(LOW:HIGH,1:NDIM) .eq. NEUMANN_INS) .or. &
      any(ins_domainBC(LOW:HIGH,1:NDIM) .eq. OUTFLOW_INS)) then                                              

  if(inou_flg) then  ! Compute Qin

  ! Compute Mass flow from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy
#if NDIM == 2
      dz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
#endif
      dydz = dy*dz
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! IAXIS:
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      ! Low X Boundary:
      if ((faces(LOW,IAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(LOW,IAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(LOW,IAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              Qinaux = Qinaux + facexData(VELC_FACE_VAR,NGUARD+1,j,k)*dydz
           enddo
        enddo
      end if
      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(HIGH,IAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(HIGH,IAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              Qinaux = Qinaux - facexData(VELC_FACE_VAR,nxc,j,k)*dydz !Into domain with - sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)

      ! JAXIS:
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
      ! Low Y Boundary:
      if ((faces(LOW,JAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(LOW,JAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(LOW,JAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux + faceyData(VELC_FACE_VAR,i,NGUARD+1,k)*dxdz
           enddo
        enddo
      end if
      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(HIGH,JAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(HIGH,JAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux - faceyData(VELC_FACE_VAR,i,nyc,k)*dxdz !Into domain with - sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
      ! KAXIS:
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
      ! Low Z Boundary:
      if ((faces(LOW,KAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(LOW,KAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(LOW,KAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux + facezData(VELC_FACE_VAR,i,j,NGUARD+1)*dxdy
           enddo
        enddo
      end if
      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .ne. NOT_BOUNDARY) .and. &
          (faces(HIGH,KAXIS) .ne. NEUMANN_INS)  .and. &
          (faces(HIGH,KAXIS) .ne. OUTFLOW_INS)) then
        ! Do The Sum:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux - facezData(VELC_FACE_VAR,i,j,nzc)*dxdy !Into domain with - sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(Qinaux,Qinout,1,FLASH_REAL,    &
                      FLASH_SUM, ins_meshComm, ierr)


   else ! Compute Qout of the domain

   ! Compute Mass flow from boundaries:
   do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy
#if NDIM == 2
      dz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
#endif
      dydz = dy*dz      
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! IAXIS:
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      ! Low X Boundary:
      if ((faces(LOW,IAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(LOW,IAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              Qinaux = Qinaux - facexData(VELC_FACE_VAR,NGUARD+1,j,k)*dydz  ! sign changed
           enddo
        enddo
      end if
      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              Qinaux = Qinaux + facexData(VELC_FACE_VAR,nxc,j,k)*dydz !Out of domain with + sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,facexData,FACEX)

      ! JAXIS:
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
      ! Low Y Boundary:
      if ((faces(LOW,JAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(LOW,JAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux - faceyData(VELC_FACE_VAR,i,NGUARD+1,k)*dxdz   ! Sign Changed
           enddo
        enddo
      end if
      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux + faceyData(VELC_FACE_VAR,i,nyc,k)*dxdz !Out of domain with + sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)

#if NDIM == 3
      ! KAXIS:
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
      ! Low Z Boundary:
      if ((faces(LOW,KAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(LOW,KAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux - facezData(VELC_FACE_VAR,i,j,NGUARD+1)*dxdy ! Sign Changed
           enddo
        enddo
      end if
      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              Qinaux = Qinaux + facezData(VELC_FACE_VAR,i,j,nzc)*dxdy !Out of domain with + sign.
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(Qinaux,Qinout,1,FLASH_REAL,    &
                      FLASH_SUM, ins_meshComm, ierr)


   end if

 endif  ! Test if there is an OUTFLOW or NEUMANN INS BC

 return
 end subroutine ins_computeQinout
