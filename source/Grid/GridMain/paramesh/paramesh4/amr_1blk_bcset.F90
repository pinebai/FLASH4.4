!!****if* source/Grid/GridMain/paramesh/paramesh4/amr_1blk_bcset
!!
!! NAME
!!
!!   amr_1blk_bcset
!!
!! SYNOPSIS
!!
!!   call amr_1blk_bcset(integer(IN)    :: mype,
!!                       integer(IN)    :: ibc,
!!                       integer(IN)    :: lb,
!!                       integer(IN)    :: pe,
!!                       integer(IN)    :: idest,
!!                       integer(IN)    :: iopt,
!!                       integer(IN)    :: ibnd,
!!                       integer(IN)    :: jbnd,
!!                       integer(IN)    :: kbnd,
!!                       integer(IN)    :: surrblks(:,:,:,:))
!!
!! DESCRIPTION
!!  
!!   This routine sets guard cell values at external boundaries for a
!!   single guard cell region of a single block that is having its
!!   guard cells filled.
!!
!!   The current block is identified by the (pe,lb) pair.
!!   Note that the following possibilities can occur:
!!   (1) pe==mype and lb <= lnblocks
!!       This is the normal case, lb is the block ID of a block local
!!       to the current PE.
!!   (2) pe==mype and lb > lnblocks
!!       In this case lb is an index that indicates a slot where information
!!       on a remote block is cached on the local PE.
!!   (3) pe/=mype
!!       In this case a remote block is referenced directly, it is the block
!!       that is known as lb on the remote PE pe. Such a call is only valid
!!       if information on that block is locally available so that such a
!!       (pe,lb) reference can be replaced by one of form (2).
!!   When this subroutine is invoked, valid block data will have been loaded
!!   into the appropriate slot (selected by idest) of either UNK1 (and, if
!!   applicable, FACEVARX1 etc.) or WORK1 temporary block storage (depending
!!   on iopt). This is true whether (pe,lb) is local or not.  The cases where
!!   (pe,lb) identifies a remote block should only occur when idest==2.
!!   !!DEV: except in the case of prolongation?
!!
!!   It can be assumed in writing this routine that guard cell regions of the
!!   current block that are of a "less diagonal" character than the current
!!   region have already been properly filled.  A guard cell region is less
!!   diagonal than another one if the set of all its point that touch the
!!   block's region of inner cells is of higher dimension.
!!
!!   It should NOT be assumed that all guard cells for this block which are
!!   not across an external boundary have already been properly filled.
!!  
!!  
!! ARGUMENTS
!!        mype   -         local processor
!!        ibc    -         the integer specifying the particular boundary
!!                          condition to be imposed, as chosen by PARAMESH from
!!                          information available to it in the surr_blks array.
!!                          The dummy argument ibc contain a negative integer,
!!                          either one of the values defined in constants.h or a
!!                          value that encodes a combination of such values.
!!                          The implementation of amr_1blk_bcset need not make
!!                          use of this value but can determine a boundary
!!                          condition to apply by looking at surr_blks
!!                          information provided in the surrblks dummy argument;
!!                          and this implementation does just that. This should
!!                          only matter for guard cells in diagonal directions,
!!                          i.e., in edge and corner regions.
!!        lb     -         block number of selected block
!!        pe     -          processor on which block lb is located
!!        idest  -         selects the storage space in data_1blk.fh which is to
!!                          be used in this call. If the leaf node is having its
!!                          guardcells filled then set this to 1, if its parent
!!                          is being filled set it to 2.
!!        iopt   -         selects which data structure to operate on;
!!                          1 for UNK, 2 for WORK
!!        ibnd   -         a selector setting designating whether the guardcells
!!                          to be set are at the left, center or right section
!!                          of the i index range:
!!                             ibnd = -1      left end
!!                                  =  0      middle
!!                                  = +1      right.  More specifically, if
!!                          ibnd is -1, the i index applied when filling unk will run
!!                          from 1:nguard, if ibnd is 0 from 1+nguard:nxb+nguard,
!!                          and if ibnd is +1 from nxb+nguard+1:nxb+2*nguard.
!!        jbnd   -         a selector setting designating whether the guarcells
!!                          to be set are at the left, center or right section
!!                          of the j index range.
!!        kbnd   -         a selector setting designating whether the guardcells
!!                          to be set are at the left, center or right section
!!                          of the k index range.
!!        surrblks -       information from the surr_blks array maintained by
!!                          PARAMESH for this block.  This array contains
!!                          information on the block's neighbors in all directions
!!                          including diagonal. For directions in which the block
!!                          is at the domain boundary, surrblks contain a negative
!!                          integer indicating the boundary condition type, either
!!                          one of the values defined in constants.h or a value
!!                          that encodes a combination of such values.
!!
!! NOTES
!!
!!   This subroutine is called from PARAMESH subroutine amr_1blk_guardcell_srl.
!!
!! HISTORY
!!
!!   Written       Peter MacNeice          August 1998
!!   Modified      Peter MacNeice          January 2001
!!   Implementation rewritten for FLASH,
!!   logic for determining BC type from surrblks:
!!                 Klaus Weide et al           .. 2007
!!***
  
subroutine amr_1blk_bcset(mype,ibc,lb,pe,idest,&
     iopt,ibnd,jbnd,kbnd,surrblks)
  
!!$  use paramesh_dimensions
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use Grid_data, ONLY : gr_bndOrder,gr_numDataStruct,gr_gridDataStruct,&
       gr_gridDataStructSize
  use gr_bcInterface, ONLY : gr_bcApplyToOneFace
  implicit none
  
#include "constants.h"      
#include "Flash.h"

  integer, intent(in) :: mype,ibc,lb,pe
  integer, intent(in) :: idest,iopt,ibnd,jbnd,kbnd
  integer, intent(in) :: surrblks(:,:,:,:)
  
  integer :: bcDir,bcType
  integer :: leftOrRight, dirToApply
  integer, dimension(MDIM) :: regionType,bnd
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC

  integer :: blockHandle
  integer :: i, struct, blockID,varCount,gridDataStruct

  bnd(1)=ibnd
  bnd(2)=jbnd
  bnd(3)=kbnd

  regionType=NO_VEC
  do i= 1,NDIM
     if(bnd(i)<0)then
        regionType(i)=LEFT_EDGE
     elseif(bnd(i)==0)then
        regionType(i)=CENTER
     else
        regionType(i)=RIGHT_EDGE
     end if
  end do
!!  print*,'in amr bcSet',regionType,ibnd,jbnd,kbnd


  ! Determine which boundary condition type, for which direction, at which face
  ! we are actually going to apply to the given guard cell region, based on the
  ! priorities expressed in bndPriorityOne and bndPriorityTwo. - KW
  call prioselector(bcType,dirToApply,leftOrRight, ibc)
  bcDir = dirToApply

#if(0)
! If we trusted Paramesh's choice:
!  bcType = ibc
! If we only had external boundaries:
!  bcType = gr_domainBC(face,bcDir)
#endif


#ifdef DEBUG_GRID
  print*, 'amr_1blk_bcset '
  print*,mype, ibc, bcDir
  print*, 'boundarytype ',bcType
  select case(bcType)
  case (PERIODIC)
     print*,'periodic'
  case (REFLECTING)
     print*,'reflecting'
  case (OUTFLOW)
     print*,'outflow'
  case (HYDROSTATIC)
     print*,'hydrostatic'
  end select
#endif


!! This function calls (possibly indirectly) the Grid_applyBCEdge with 2*nguard cells if
!! it is filling guardcells in cell centered data structures. To fill guard cells for
!! face centered data structures, for FACEX it calls IAXIS with 2*nguard+1 cells, while
!! the other two directions are called with 2*nguard again. The same logic applies
!! to FACEY and FACEZ data structures. The name of the gridDataStructure is also
!! passed to Grid_applyBCEdge, so it can fill the cells appropriately.
!
!        guardcell guardcell guardcell guardcell  interior interior interior interior
!       |---------|---------|---------|---------|---------|--------|--------|--------|      
!       |         |         |         |         |         |        |        |        |      
!       |         |         |         |         |         |        |        |        |      
!       |---------|---------|---------|---------|---------|--------|--------|--------|      
!       |         |         |         |         |         |        |        |        |      
! facevarx1(1)   (2)       (3)       (4)       (5)       (6)      (7)       (8)     (9)


  ! Get a block handle, we may have been called for a remote block; in that case we assume that
  ! the current processor has valid cached information on that block that will be accessible
  ! via the block handle. - KW
  call internal_getBlkHandle(lb, pe, blockHandle)

  if (regionType(bcDir).NE.leftOrRight) then
     print*,'Internal Error in amr_1blk_bcset',regionType(bcDir),leftOrRight
     call Driver_abortFlash('Internal Error in amr_1blk_bcset')
  end if
  regionType(bcDir)=leftOrRight
  blockID=blockHandle           !only used for Grid_getBlkIndexLimits, which does not really care - KW
  if(iopt==2) then
     varCount=1
     gridDataStruct=WORK
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,&
          CENTER)
     !! Paramesh supports Work for only cell centered data, so CENTER
     !! is the appropriate gridDataStruct for this call

     !! Now that all preparatory work is done, do the calculation.
     !! The three different loops for the three directions are there
     !! because the order of the indices in the array for access changes
     !! and therefore the loops change.
     call gr_bcApplyToOneFace(bcDir,bcType,&
          gridDataStruct,varCount,regionType,&
          blkLimits,blkLimitsGC,blockHandle,idest)
  else
     do struct=1,gr_numDataStruct
        gridDataStruct=gr_gridDataStruct(struct)
        varCount=gr_gridDataStructSize(struct)
     
        !! The next statements gest block index information to 
        !! prepare for the calculation.
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC,&
             gridDataStruct)
        !! Now that all preparatory work is done, call the routine
        !! that will repackage relevant parts of the block data
        !! to pass on the boundary condition routines that do actual
        !! calculation.
        !! The boundary condition routines will see data arrays
        !! that may be equivalent to some transpose of the block
        !! data, depending on the direction given by axis.
!!        print*,bcDir,bcType,gridDataStruct,varCount,regionType(1:NDIM)
        call gr_bcApplyToOneFace(bcDir,bcType,&
             gridDataStruct,varCount,regionType,&
             blkLimits,blkLimitsGC,blockHandle,idest)
     end do
  end if

  return

contains

  subroutine internal_getBlkHandle(blockID, pe, accessBlkOut)
    use tree, ONLY: laddress, strt_buffer, last_buffer
    integer, intent(in) :: blockID, pe
    integer, intent(OUT) :: accessBlkOut

    integer,save :: accessBlk = 0       ! what we return, remember last value between calls - KW
    integer :: iblk
    logical :: lfound

    lfound = .FALSE.
    ! The following logic copied from Tests/amr_1blk_bcset.F90 in Paramesh4 distribution. - KW
    ! This is heavily dependent on Paramesh internals.
    if (pe.EQ.mype) then
       accessBlk = blockID
       lfound = .TRUE.
    else if(accessBlk .gt. 0 .and. &
         accessBlk .ge. strt_buffer .AND. accessBlk .le. last_buffer) then
       if(blockID.eq.laddress(1,accessBlk).and.pe.eq.laddress(2,accessBlk) ) &
          lfound = .true.          ! found cached accessBlk from previous call still valid
    end if
    if (.NOT. lfound) then
       lfound = .false.
       iblk = strt_buffer
       do while(.not.lfound.and.iblk.le.last_buffer)
          if(blockID.eq.laddress(1,iblk).and.pe.eq.laddress(2,iblk) ) then
             accessBlk = iblk
             lfound = .true.
          endif
          iblk = iblk+1
       enddo
       if(.not.lfound) then
          call Driver_abortFlash('Paramesh error: amr_1blk_bcset : '// & 
               &      ' remote block is not in list received on this pe')
       endif
    endif

    accessBlkOut = accessBlk
  end subroutine internal_getBlkHandle


    logical function isBCNumber(bc)
      integer,intent(IN) :: bc
      isBCNumber = (bc <= PARAMESH_PHYSICAL_BOUNDARY .AND. bc .NE. PERIODIC)
    end function isBCNumber

    subroutine setLowOrHigh(loh, axis, dirs)
      integer,intent(OUT) :: loh
      integer,intent(IN)  :: axis,dirs(MDIM)
      loh = LOW
      if (dirs(axis) > 0) loh = HIGH
    end subroutine setLowOrHigh


    subroutine prioselector(bcTypeToApply,directionToApply,leftOrRight,bcTypeFromParamesh)
      use Driver_interface, ONLY : Driver_abortFlash
      integer,intent(OUT) :: bcTypeToApply,directionToApply,leftOrRight
      integer,intent(IN) :: bcTypeFromParamesh !could be used for debugging - KW
      integer :: celery(2,2,2) !for faster(?) access to relevant surrblks info - KW
#ifdef DEBUG_GRID
      integer :: celery0(2,2,2) !for debugging
#endif
      integer :: dirSet(MDIM), dirSetToOrder(MDIM), orderedAxes(MDIM)
      integer :: lowOrHigh
      integer :: dimendiff, nf,directionAux
      integer :: i,j,k

      celery = 0
      
#ifdef DEBUG_GRID
      celery0 = 0
      do k=0,kbnd,sign(1,kbnd)
         do j=0,jbnd,sign(1,jbnd)
            do i=0,ibnd,sign(1,ibnd)
               celery0(1+abs(i),1+abs(j),1+abs(k)) = surrblks(1,2+i,2+j,2+k)
            end do
         end do
      end do
#endif

      do k=0,abs(kbnd)
         do j=0,abs(jbnd)
            do i=0,abs(ibnd)
               celery(1+i,1+j,1+k) = surrblks(1,2+i*ibnd,2+j*jbnd,2+k*kbnd)
            end do
         end do
      end do

#ifdef DEBUG_GRID
      if (ANY(celery .NE. celery0)) then
         call Driver_abortFlash('Something is wrong with initialization of celery array!')
      end if
#endif

      dimenDiff = abs(ibnd) + abs(jbnd) + abs(kbnd)
      dirSet(1) = ibnd
      dirSet(2) = jbnd
      dirSet(3) = kbnd
      select case (dimendiff)
      case(0)
         call Driver_abortFlash('amr_1blk_bcset should not have been called with all of i/j/kbnd zero!');
      case(1)
         if (isBCNumber(celery(2,1,1))) then
            directionToApply = IAXIS
         else if (isBCNumber(celery(1,2,1))) then
            directionToApply = JAXIS
         else
            directionToApply = KAXIS
         end if
         call setLowOrHigh(lowOrHigh, directionToApply, dirSet)
         bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSet,celery,surrblks)
      case(2)
#define NODIR 0
         dirSetToOrder = dirSet
         if (isBCNumber(celery(2,1,1))) dirSetToOrder(1) = 0
         if (isBCNumber(celery(1,2,1))) dirSetToOrder(2) = 0
#if NDIM > 2
         if (isBCNumber(celery(1,1,2))) dirSetToOrder(3) = 0
#endif
         if (ANY(dirSetToOrder(1:NDIM) .NE. 0)) then
            call orderByPriority(dirSetToOrder, orderedAxes, 2)
!!              print*,'case 2A0:',dirSetToOrder,orderedAxes
            if (orderedAxes(2) .EQ. NODIR) then
               dirSetToOrder = dirSet
               dirSetToOrder(orderedAxes(1)) = 0
               call orderByPriority(dirSetToOrder, orderedAxes, 1)
               directionToApply = orderedAxes(1)
            else
!!                  print*,'case 2Aprime:',dirSetToOrder,orderedAxes(1)
!!!!               directionToApply = orderedAxes(2)
               directionToApply = orderedAxes(1)
            end if
!!              print*,'case 2A:',dirSetToOrder,directionToApply
            call setLowOrHigh(lowOrHigh, directionToApply, dirSet)
            bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSet,celery,surrblks)
         else
            call orderByPriority(dirSet, orderedAxes, 1)
            directionToApply = orderedAxes(1)
            dirSetToOrder(directionToApply) = dirSet(directionToApply)
!            print*,'case 2B',dirSetToOrder,directionToApply
            call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
            bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
         end if
!!#if NDIM > 2
      case(3)
         if (isBCNumber(celery(2,2,1)) .AND. isBCNumber(celery(2,1,2)) .AND. isBCNumber(celery(1,2,2))) then
            nf = 0
            if (isBCNumber(celery(2,1,1))) nf = 1
            if (isBCNumber(celery(1,2,1))) nf = nf + 1
            if (isBCNumber(celery(1,1,2))) nf = nf + 1
            select case (nf)
            case(0)
               directionToApply = gr_bndOrder(1)
               dirSetToOrder = dirSet
               dirSetToOrder(gr_bndOrder(2)) = 0
               call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
               bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
            case(1)
               directionToApply = gr_bndOrder(1)
               dirSetToOrder = dirSet
               dirSetToOrder(directionToApply) = 0
               if (isBCNumber(celery(1+abs(dirSetToOrder(1)),1+abs(dirSetToOrder(2)),1+dirSetToOrder(3)))) then
                  dirSetToOrder = dirSet
                  dirSetToOrder(gr_bndOrder(2)) = 0
                  call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
                  bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
               else
                  if (isBCNumber(celery(2,1,1))) then
                     directionAux = IAXIS
                  else if (isBCNumber(celery(1,2,1))) then
                     directionAux = JAXIS
                  else
                     directionAux = KAXIS
                  end if
                  dirSetToOrder = dirSet
                  dirSetToOrder(directionAux) = 0
                  call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
                  bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
               end if
            case(2)
               if (.NOT. isBCNumber(celery(2,1,1))) then
                  directionAux = IAXIS
               else if (.NOT. isBCNumber(celery(1,2,1))) then
                  directionAux = JAXIS
               else
                  directionAux = KAXIS
               end if
               dirSetToOrder = dirSet
               dirSetToOrder(directionAux) = 0
               call orderByPriority(dirSetToOrder,orderedAxes,2)
               directionToApply = orderedAxes(1)
               dirSetToOrder = dirSet
               dirSetToOrder(orderedAxes(2)) = 0
               call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
               bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
            case(3)
               directionToApply = gr_bndOrder(1)
               dirSetToOrder = 0
               dirSetToOrder(directionToApply) = dirSet(directionToApply)
               call setLowOrHigh(lowOrHigh, directionToApply, dirSetToOrder)
               bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSetToOrder,celery,surrblks)
            case default
               call Driver_abortFlash('Impossible nf in amr_1blk_bcset!')
            end select
         else
            dirSetToOrder = 0
            if (.NOT. isBCNumber(celery(2,2,1))) dirSetToOrder(KAXIS) = dirSet(KAXIS)
            if (.NOT. isBCNumber(celery(2,1,2))) dirSetToOrder(JAXIS) = dirSet(JAXIS)
            if (.NOT. isBCNumber(celery(1,2,2))) dirSetToOrder(IAXIS) = dirSet(IAXIS)
            call orderByPriority(dirSetToOrder, orderedAxes, 1)
            directionToApply = orderedAxes(1)
            call setLowOrHigh(lowOrHigh, directionToApply, dirSet)
            bcTypeToApply = bbbc(directionToApply,lowOrHigh,dirSet,celery,surrblks)
         end if
!!#endif
      case default
         call Driver_abortFlash('Impossible dimendiff in amr_1blk_bcset!');
      end select

      leftOrRight = LEFT_EDGE
      if (lowOrHigh==HIGH) leftOrRight = RIGHT_EDGE
    end subroutine prioselector

    integer function bbbc(directionToApply,lowOrHigh,dirset,celery,surrblks)
      integer,intent(IN) :: directionToApply,lowOrHigh
      integer,intent(IN) :: dirset(MDIM)
      integer,intent(IN) :: celery(2,2,2),surrblks(:,:,:,:)
      bbbc = extractBCForDirection(directionToApply,3-lowOrHigh,dirset(1),dirset(2),dirset(3),surrblks)
    end function bbbc

    subroutine orderByPriority(inSet, outSeq, outSeqLen)
      integer,intent(IN)  :: inSet(MDIM)
      integer,intent(OUT) :: outSeq(MDIM)
      integer,intent(IN)  :: outSeqLen
      integer :: i,n
      outSeq = NODIR
      n = 0
!!      print*,'order:',inSet
      do i=1,NDIM
         if (inSet(gr_bndOrder(i)) .NE. 0) then
            n = n + 1
!!            print*,i,'.) adding', gr_bndOrder(i),' as ',n 
            outSeq(n) = gr_bndOrder(i)
            if (n .GE. outSeqLen) return
         end if
      end do
    end subroutine orderByPriority

    integer function extractBCForDirection(idir,lowOrHigh,ibnd,jbnd,kbnd,surrblks)
      use gr_interface, ONLY: gr_extractBCForDirection
      integer,intent(IN) :: idir,ibnd,jbnd,kbnd,lowOrHigh
      integer,intent(IN) :: surrblks(:,:,:,:)
      extractBCForDirection = gr_extractBCForDirection(surrblks(1,2+ibnd,2+jbnd,2+kbnd), idir,lowOrHigh)
    end function extractBCForDirection
end subroutine amr_1blk_bcset
