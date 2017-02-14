!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_rescaleVelout
!!
!!
!! NAME
!!
!!  ins_rescaleVelout
!!
!!
!! SYNOPSIS
!!
!!  ins_rescalVelout (integer(IN) :: blockCount,
!!                    integer(IN) :: blockList(blockCount)
!!                    real(IN   ) :: Qin
!!                    real(IN)    :: Qout)
!!
!!
!! DESCRIPTION
!!
!! Rescales outflow velocities of domain boundaries. This routine is part 
!! of a global mass balance strategy.
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  Qin        - Flow into the domain.
!!  Qout       - Flow out of the domain.
!!
!!***

subroutine ins_rescaleVelout( blockCount, blockList, Qin, Qout)

#include "Flash.h"

  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

   use IncompNS_data, ONLY : ins_meshMe, ins_domainBC

  implicit none

#include "constants.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
  real,    INTENT(IN) :: Qin, Qout
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  integer :: lb,blockID,ierr,i,j,k

  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer :: nxc,nyc,nzc

  real :: Qindout

  real, parameter :: eps = 1.e-11

!=============================================================================

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
  nzc = NZB + NGUARD + 1


  ! Detect if the problem is an outflow problem to proceed with flow computation.
  if (any(ins_domainBC(LOW:HIGH,1:NDIM) .eq. NEUMANN_INS) .or. &
      any(ins_domainBC(LOW:HIGH,1:NDIM) .eq. OUTFLOW_INS)) then 

  ! Check if Qin and Qout are zero
  if ((Qin .lt. eps) .and. (ins_meshMe .eq. MASTER_PE))  &
     write(*,*) 'ins_rescaleVelout: Warning Qin=',Qin,' <',eps  
  if ((Qout .lt. eps) .and. (ins_meshMe .eq. MASTER_PE)) &
     write(*,*) 'ins_rescaleVelout: Warning Qout=',Qout,' <',eps  

  Qindout = 1.
  if (abs(Qout) .gt. eps*abs(Qin)) Qindout = Qin/Qout

  if(ins_meshMe .eq. MASTER_PE) write(*,*) 'Qin, Qout, Qin/Qout=',Qin,Qout,Qindout

  ! Rescale Velocities from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)

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
              facexData(VELC_FACE_VAR,NGUARD+1,j,k) = facexData(VELC_FACE_VAR,NGUARD+1,j,k)*Qindout
           enddo
        enddo
      end if
      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              facexData(VELC_FACE_VAR,nxc,j,k) = facexData(VELC_FACE_VAR,nxc,j,k)*Qindout
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
              faceyData(VELC_FACE_VAR,i,NGUARD+1,k) = faceyData(VELC_FACE_VAR,i,NGUARD+1,k)*Qindout
           enddo
        enddo
      end if
      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              faceyData(VELC_FACE_VAR,i,nyc,k) = faceyData(VELC_FACE_VAR,i,nyc,k)*Qindout
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
              facezData(VELC_FACE_VAR,i,j,NGUARD+1) = facezData(VELC_FACE_VAR,i,j,NGUARD+1)*Qindout
           enddo
        enddo
      end if
      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. NEUMANN_INS)  .or. &
          (faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! Do The Sum:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              facezData(VELC_FACE_VAR,i,j,nzc) = facezData(VELC_FACE_VAR,i,j,nzc)*Qindout
           enddo
        enddo
      end if
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif

   enddo


   
 endif ! Test if there is an OUTFLOW or NEUMANN INS BC

 return
 end subroutine ins_rescaleVelout
