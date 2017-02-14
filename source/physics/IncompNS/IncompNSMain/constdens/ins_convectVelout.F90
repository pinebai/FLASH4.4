!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_convectVelout
!!
!!
!! NAME
!!
!!  ins_convectVelout
!!
!!
!! SYNOPSIS
!!
!!  ins_convectVelout(integer(IN) :: blockCount,
!!                    integer(IN) :: blockList(blockCount)
!!                    real  (OUT) :: convvel(LOW:HIGH,MDIM))
!!
!!
!! DESCRIPTION
!!
!! Computes convective velocities out of OUTFLOW_INS domain boundaries. This routine is part 
!! of a global mass balance strategy. FOR THE MOMENT OUTFLOW_INS is defined only on HIGH x,y or z. 
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  convvel    - Velocities out of domain for OUTFLOW_INS Boundary conditions.
!!
!!***

subroutine ins_convectVelout( blockCount, blockList, convvel)

#include "Flash.h"

  use Grid_interface, only : Grid_getDeltas,         &
                             Grid_getBlkBC,          &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr,     &
                             Grid_getBlkIndexLimits, &
                             Grid_solvePoisson, Grid_getBlkBoundBox, Grid_getBlkCenterCoords

  use IncompNS_data, only : ins_domainBC, ins_globalDomain

#if NDIM == 2
  use IncompNS_data, only : uvel_x,vvel_x,wvel_x,uvel_y,vvel_y,wvel_y
#elif NDIM == 3  
  use IncompNS_data, only : uvel_x,vvel_x,wvel_x,uvel_y,vvel_y,wvel_y, &
                            uvel_z,vvel_z,wvel_z
#endif
  use IncompNS_data, ONLY : ins_meshComm

  implicit none

#include "constants.h"
  include "Flash_mpi.h"


  !! ---- Argument List ----------------------------------
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(MAXBLOCKS) :: blockList 
  real,    INTENT(OUT) :: convvel(LOW:HIGH,MDIM)
  !! -----------------------------------------------------

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
            
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData

  integer :: lb,blockID,ierr,i,j,k

  integer :: faces(2,MDIM),onBoundary(2,MDIM)

  integer :: nxc,nyc,nzc

  real :: convveli(LOW:HIGH,MDIM),del(MDIM),dx,dy,dz,dxdy,dydz,dxdz,Lx,Ly,Lz

  real :: factorarea

!=============================================================================

  nxc = NXB + NGUARD + 1
  nyc = NYB + NGUARD + 1
#if NDIM == 2
  nzc = -1
#elif NDIM == 3
  nzc = NZB + NGUARD + 1
#endif
  convveli = 0.
  convvel  = 0.

  ! Detect if the problem is an outflow problem to proceed with flow computation.
  ! ONLY for HIGH x,y,z.
  if (any(ins_domainBC(HIGH,1:NDIM) .eq. OUTFLOW_INS)) then    

  ! Compute Mass flow from boundaries:
  do lb = 1,blockCount
      blockID = blockList(lb)
      ! Get blocks dx, dy ,dz:
      call Grid_getDeltas(blockID,del)
      dx = del(IAXIS)
      dy = del(JAXIS)
      dxdy = dx*dy

      Lx =(ins_globalDomain(HIGH,IAXIS)-ins_globalDomain(LOW,IAXIS))
      Ly =(ins_globalDomain(HIGH,JAXIS)-ins_globalDomain(LOW,JAXIS))

#if NDIM == 2
      dz = 1.
      Lz = 1.
#elif NDIM == 3
      dz = del(KAXIS)
      Lz = (ins_globalDomain(HIGH,KAXIS)-ins_globalDomain(LOW,KAXIS))
#endif
      dydz = dy*dz      
      dxdz = dx*dz
      ! Get Blocks internal limits indexes:
      call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC) 
      ! Get blocks BCs:
      call Grid_getBlkBC(blockID,faces,onBoundary)
 
      ! AXIS:
      call Grid_getBlkPtr(blockID,facexData,FACEX)
      call Grid_getBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_getBlkPtr(blockID,facezData,FACEZ)
#endif

      ! High X Boundary:
      if ((faces(HIGH,IAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Unxc-1
              uvel_x(NGUARD,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc-1,j,k)
              !Unxc
              uvel_x(NGUARD+1,j,k,HIGH,blockID) = facexData(VELC_FACE_VAR,nxc,j,k)
        
            enddo
         enddo

        factorarea = dy*dz/(Ly*Lz)
        ! Compute the Mean Convective velocity across X High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
              
              convveli(HIGH,IAXIS) = convveli(HIGH,IAXIS) + & 
                                     facexData(VELC_FACE_VAR,nxc,j,k)*factorarea
            enddo
         enddo

        ! redistribute velocities for V and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
              !Vnxc-1
              vvel_x(NGUARD-1,j,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,nxc-1,j,k)
              !Vnxc
              vvel_x(NGUARD,j,k,HIGH,blockID)   = faceyData(VELC_FACE_VAR,nxc,j,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              !Wnxc-1
              wvel_x(NGUARD-1,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc-1,j,k)
              !Wnxc
              wvel_x(NGUARD,j,k,HIGH,blockID) =  facezData(VELC_FACE_VAR,nxc,j,k)

           enddo
        enddo
#endif

      endif

      ! High Y Boundary:
      if ((faces(HIGH,JAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnyc-1
              !faceyData(VELC_FACE_VAR,i,nyc+1,k) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              vvel_y(NGUARD,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc-1,k)
              !Vnyc
              !faceyData(VELC_FACE_VAR,i,nyc+2,k) = faceyData(VELC_FACE_VAR,i,nyc,k) 
              vvel_y(NGUARD+1,i,k,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,nyc,k)
            enddo
         enddo


        factorarea = dx*dz/(Lx*Lz)
        ! Compute the Mean Convective velocity across Y High Boundary:
        do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,JAXIS) = convveli(HIGH,JAXIS) + & 
                                     faceyData(VELC_FACE_VAR,i,nyc,k)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and W:
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unyc-1
              !facexData(VELC_FACE_VAR,i,nyc+1,k) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              uvel_y(NGUARD-1,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc-1,k)
              !Unyc
              uvel_y(NGUARD,i,k,HIGH,blockID) = facexData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo

#if NDIM == 3
        do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnyc-1
              !facezData(VELC_FACE_VAR,i,nyc+1,k) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              wvel_y(NGUARD-1,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc-1,k)
              !Wnyc
              wvel_y(NGUARD,i,k,HIGH,blockID) = facezData(VELC_FACE_VAR,i,nyc,k)
           enddo
        enddo
#endif


      endif


#if NDIM == 3

      ! High Z Boundary:
      if ((faces(HIGH,KAXIS) .eq. OUTFLOW_INS)) then
        ! redistribute velocities to guardcells:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Wnxc-1
              !facezData(VELC_FACE_VAR,i,j,nzc+1) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              wvel_z(NGUARD,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc-1)
              !Wnxc
              !facezData(VELC_FACE_VAR,i,j,nzc+2) = facezData(VELC_FACE_VAR,i,j,nzc) 
              wvel_z(NGUARD+1,i,j,HIGH,blockID) = facezData(VELC_FACE_VAR,i,j,nzc)
            enddo
         enddo


        factorarea = dx*dy/(Lx*Ly)
        ! Compute the Mean Convective velocity across Z High Boundary:
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              
              convveli(HIGH,KAXIS) = convveli(HIGH,KAXIS) + & 
                                     facezData(VELC_FACE_VAR,i,j,nzc)*factorarea
            enddo
         enddo
        

        ! redistribute velocities for U and V:
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)+1
              !Unxc-1
              !facexData(VELC_FACE_VAR,i,j,nzc+1) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              uvel_z(NGUARD-1,i,j,HIGH,blockID) = facexData(VELC_FACE_VAR,i,j,nzc-1)
              !Unxc
              uvel_z(NGUARD,i,j,HIGH,blockID) =  facexData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)+1
           do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              !Vnxc-1
              !faceyData(VELC_FACE_VAR,i,j,nzc+1) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              vvel_z(NGUARD-1,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc-1)
              !Vnxc
              vvel_z(NGUARD,i,j,HIGH,blockID) = faceyData(VELC_FACE_VAR,i,j,nzc)
           enddo
        enddo

      endif
#endif

      call Grid_releaseBlkPtr(blockID,facexData,FACEX)
      call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
#if NDIM == 3
      call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
#endif


   enddo

   ! Gather total inflow volume flow ratio:
   call MPI_Allreduce(convveli,convvel,(HIGH-LOW+1)*MDIM,FLASH_REAL,    &
                      FLASH_SUM, ins_meshComm, ierr)


 endif ! Test if there is an OUTFLOW or NEUMANN INS BC

 return
 end subroutine ins_convectVelout
