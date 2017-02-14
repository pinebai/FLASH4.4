!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_computeDt
!!
!! NAME
!!
!!  IncompNS_computeDt
!!
!!
!! SYNOPSIS
!!
!!  
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***


subroutine IncompNS_computeDt(ins_mindt,ins_minloc)

  use Grid_interface, ONLY : Grid_getListOfBlocks,  &
                             Grid_getBlkIndexLimits,&
                             Grid_getDeltas,        &
                             Grid_getBlkPtr,        &
                             Grid_releaseBlkPtr

  use ins_interface, ONLY : ins_computeDtLocal

  use IncompNS_data, ONLY : ins_useIncompNS

  implicit none
#include "constants.h"
#include "Flash.h"
  real, intent(INOUT) :: ins_mindt
  integer, intent(INOUT) :: ins_minloc(5)

  ! Local Variables:
  real, PARAMETER :: MAX_TSTEP = huge(1.0)
  real    :: dtLocal
  integer :: lminloc(5)
  integer, dimension(MAXBLOCKS) :: blockList
  integer :: numLeafBlocks

  !!prepatory data structures for passing coords to timestep routines
  real, dimension(MDIM) :: del


  !!arrays which hold the starting and ending indices of a block
  integer,dimension(2,MDIM)::blkLimits,blkLimitsGC

  !!coordinate infomration to be passed into physics  
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  integer :: isize,jsize,ksize

  integer :: i, blockID


  if (.NOT. ins_useIncompNS) RETURN

  !!Initialize all timestep variables.
  dtLocal    = MAX_TSTEP
  lminloc(:) = 0     
 
  !! Loop over all local leaf-node blocks
  call Grid_getListOfBlocks(LEAF,blockList, numLeafBlocks)

  do i = 1, numLeafBlocks

     !!Get the coordinate information for all the
     call Grid_getBlkIndexLimits(blockList(i),blkLimits,blkLimitsGC)
     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

     blockID=blockList(i)

     call Grid_getDeltas(blockID, del)

     call Grid_getBlkPtr(blockID, facexData,FACEX)
     call Grid_getBlkPtr(blockID, faceyData,FACEY)
#if NDIM == 3
     call Grid_getBlkPtr(blockID, facezData,FACEZ)
#endif



#ifdef DEBUG_DRIVER
     print*,'going to call INS timestep'
#endif

     ! ins_computeDtLocal computed de blockID min_dt
     ! which if is .lt. dtLocal, is assigned to dtLocal
     call ins_computeDtLocal (blockID,         &
                         isize, jsize, ksize,  &
              del(IAXIS),del(JAXIS),del(KAXIS),&
                         blkLimits,blkLimitsGC,&
                         facexData,faceyData,  &
                         facezData,            &
                         dtLocal,lminloc)



     call Grid_releaseBlkPtr(blockID, facexData, FACEX)
     call Grid_releaseBlkPtr(blockID, faceyData, FACEY)
#if NDIM == 3
     call Grid_releaseBlkPtr(blockID, facezData, FACEZ)
#endif

#ifdef DEBUG_DRIVER
     print*,'returned from INS timestep'
#endif

  enddo

  ! Assign to  ins_mindt
  ins_mindt = dtLocal
  ins_minloc= lminloc
  return

end subroutine IncompNS_computeDt
