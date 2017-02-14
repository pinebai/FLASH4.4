!!****if* source/physics/IncompNS/IncompNSStats/IncompNS_statsIOExport
!!
!! NAME
!!
!!  IncompNS_statsIOExport
!!
!! SYNOPSIS
!!
!!  call IncompNS_statsIOExport(logical(in) :: expt_flag)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   expt_flag : 
!!
!! AUTOGENROBODOC
!!
!!
!!***



subroutine IncompNS_statsIOExport(expt_flag)

  use Grid_interface, ONLY : Grid_getListOfBlocks,Grid_getBlkPtr,Grid_releaseBlkPtr,&
                             Grid_getDeltas,Grid_getBlkBoundBox,Grid_getBlkCenterCoords

  use IncompNS_data, only : ins_meshMe, ins_useIncompNS

  implicit none
#include "constants.h"
#include "Flash.h"

  logical, intent(in) :: expt_flag

  ! Local Variables

  integer :: blockCount
  integer :: blockList(MAXBLOCKS)

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData

  integer :: lb,blockID,i,j,jj,k,iu
  real :: rnpts
  real, dimension(MDIM)  :: coord,bsize,del
  real ::  boundBox(2,MDIM)

  character(20) :: filename  

  real :: uavg_y(NYB,MAXBLOCKS), vavg_face_y(NYB+1,MAXBLOCKS),vavg_y(NYB,MAXBLOCKS),wavg_y(NYB,MAXBLOCKS)
  real :: uuavg_y(NYB,MAXBLOCKS),uvavg_y(NYB,MAXBLOCKS),uwavg_y(NYB,MAXBLOCKS)
  real :: vvavg_y(NYB,MAXBLOCKS),vwavg_y(NYB,MAXBLOCKS),wwavg_y(NYB,MAXBLOCKS)
  real :: ycell(NYB,MAXBLOCKS)
  real :: auxavg_y(NYB,MAXBLOCKS)

  ! External helper function

  integer :: ut_getFreeFileUnit

#ifdef FLASH_GRID_UG 
  if (.NOT. ins_useIncompNS) RETURN

  if (expt_flag) then ! Export

  ! Get List of leaf Blocks:
  call Grid_getListOfBlocks(LEAF,blockList,blockCount)

  ! Loop on Blocks
  do lb = 1,blockCount
     blockID = blockList(lb)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)


     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,solnData,CENTER)
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,facezData,FACEZ)

     ! Plane Averages:
     rnpts = real(NXB*NZB)
     do j = GRID_JLO,GRID_JHI
       jj = j-GRID_JLO+1   

       ycell(jj,blockID) = coord(JAXIS)-0.5*bsize(JAXIS)+real(j-NGUARD-1)*del(JAXIS)+0.5*del(JAXIS) 

       ! Velocities:
       uavg_y(jj,blockID) = sum(facexData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       vavg_y(jj,blockID) = 0.5/rnpts*(sum(faceyData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,j  ,GRID_KLO:GRID_KHI))+&
                                       sum(faceyData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,j+1,GRID_KLO:GRID_KHI)))

       wavg_y(jj,blockID) = sum(facezData(VAVG_FACE_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       ! Reynolds Stresses:
       uuavg_y(jj,blockID)= sum(solnData(UUAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       uvavg_y(jj,blockID)= sum(solnData(UVAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       uwavg_y(jj,blockID)= sum(solnData(UWAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       vvavg_y(jj,blockID)= sum(solnData(VVAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       vwavg_y(jj,blockID)= sum(solnData(VWAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

       wwavg_y(jj,blockID)= sum(solnData(WWAV_VAR,GRID_ILO:GRID_IHI,j,GRID_KLO:GRID_KHI))/rnpts

     enddo


     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  enddo

  ! for now each proc writes a file:
  iu = ut_getFreeFileUnit()
  write(filename, '("IOData/Ystats.", i5.5)') ins_meshMe
  open(unit=iu, file=trim(filename), status='unknown')
  
  do lb = 1,blockCount
     blockID = blockList(lb)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     write(iu,66) blockID,coord(IAXIS:KAXIS),del(IAXIS:KAXIS)
     do jj=1,NYB
       write(iu,67) ycell(jj,blockID),uavg_y(jj,blockID), vavg_y(jj,blockID), wavg_y(jj,blockID),&
                                     uuavg_y(jj,blockID),uvavg_y(jj,blockID),uwavg_y(jj,blockID),&
                                     vvavg_y(jj,blockID),vwavg_y(jj,blockID),wwavg_y(jj,blockID)
     enddo
  enddo   

  close(iu)

  endif
#endif

  return

66   format(i4.4,6g18.10)
67   format(10g18.10)

end subroutine

