!!****if* source/Grid/GridMain/gr_GCTransferOneBlk
!!
!! NAME
!!  gr_GCTransferOneBlk
!!
!! SYNOPSIS
!!
!!  gr_GCTransferOneBlk(logical(IN)  :: mode,
!!                      integer(IN)  :: indCnt,
!!                      integer(IN)  :: indList(indCnt),
!!                      integer(IN)  :: offset,
!!                      integer(IN)  :: blkLimits(LOW:HIGH,MDIM),
!!                      integer(IN)  :: blkLimitsGC(LOW:HIGH,MDIM),
!!                      real,pointer :: flatArray(:),
!!                      real,pointer :: blkArray(:,:,:,:))
!!  
!! DESCRIPTION
!!  
!!  This routine can transfers data from a flat array (used for storing saved
!!  guard cell values) to a block array (which will have a full block structure
!!  but invalid interior values; valid values only in the guard cells specified
!!  at the time of their storage), or it can extract data from the mesh arrays one
!!  at a time and store it in the flat array. Normally users are advised to use the
!!  Grid_GCputScratch interface for transferring the data to the flat array, and use 
!!  the current interface for fetching the data from the flat array.
!!
!! ARGUMENTS
!!            
!!  mode     : determines whether data is being transferred from the flatArray 
!!             (mode = .false.), or into the flatArray (mode = .true.)
!!  indCnt   : count of variable indices in the mesh data structure to be 
!!             transferred
!!  indList  : list of the indices, this argument need valid values for mode=.true.
!!             and is ignored for mode = .false. This is because when data are 
!!             being fetched from the flatarray, it is assumed that temporary 
!!             block storage is provided by the calling routine.
!!  offset   : the offset of this block's guardcells in the flatArray
!!  blkLimits: the LOW and HIGH indices of the interior of the block along
!!             each dimension
!!  blkLimitsGC: the LOW and HIGH indices of the block including the guardcells
!!             Note that the values contained in the array here may differ from the
!!             values returned by a call to Grid_getBlkIndexLimits because the
!!             calling routine may only be interested in saving a subset of 
!!             guardcells. The indices here reflect the smaller block if fewer 
!!             than existing guardcells are saved
!!  flatArray : pointer to a single dimensional array that stores all the guardcells
!!             for the grid data structure of interest.
!!  blkArray : a four dimensional real array, that can either be a block from the 
!!             mesh solution data or a temporary
!!
!!  NOTES
!!   variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the fortran module Grid_data. Variables are not
!!   starting with gr_ are local variables or arguments passed to the 
!!   routine.
!!
!!  SEE ALSO
!!   
!!   Grid_GCTransferOneBlk and gr_GCAllocScratch
!!
!!***

#include "constants.h"
#include "Flash.h"

subroutine gr_GCTransferOneBlk(mode,indCnt,indList,offset,&
     blkLimits,blkLimitsGC,&
     flatArray,blkArray)

  implicit none

  integer, intent(IN) :: indCnt,offset
  logical, intent(IN) :: mode
  integer,dimension(indCnt),intent(IN) :: indList
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
  real, pointer, dimension(:) :: flatArray
  real, pointer, dimension(:,:,:,:) :: blkArray

  integer :: i,j,k, n, ind
  logical :: notIntK,notIntJ,notIntI

  notIntK=(NDIM<3)
  n= offset
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     notIntK=notIntK.or.&
          ((k<blkLimits(LOW,KAXIS)).or.(k>blkLimits(HIGH,KAXIS)))
     if(notIntK) then
        notIntJ=(NDIM<2)
        do j=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           notIntJ=notIntJ.or.&
                ((j<blkLimits(LOW,JAXIS)).or.(j>blkLimits(HIGH,JAXIS)))
           if(notIntJ) then
              do i=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
                 notIntI= (i<blkLimits(LOW,IAXIS)).or.(i>blkLimits(HIGH,IAXIS))
                 if(notIntI) then
                    do ind = 1,indCnt
                       n=n+1
                       if(mode) then
                          flatArray(n)=blkArray(indList(ind),i,j,k)
                       else
                          blkArray(ind,i,j,k)=flatArray(n)
                       end if
                    end do
                 end if
              end do
           end if
        end do
     end if
  end do
  
end subroutine gr_GCTransferOneBlk
