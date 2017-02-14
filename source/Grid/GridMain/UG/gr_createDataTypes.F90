!!****if* source/Grid/GridMain/UG/gr_createDataTypes
!!
!! NAME
!!  gr_createDataTypes
!!
!! SYNOPSIS
!!
!!  call gr_createDataTypes()
!!  
!! DESCRIPTION 
!!  create the MPI derived Datatypes in the uniform grid for exchanging
!!  guardcells.  It is more efficient to make a datatype than to do multiple
!!  mpi send/recvs
!!  
!!  
!! ARGUMENTS 
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

subroutine gr_createDataTypes()

  use Grid_data, ONLY : gr_exch,&
       gr_numDataStruct,gr_gridDataStruct, gr_gridDataStructSize
  use Grid_interface, ONLY : Grid_getBlkIndexLimits

implicit none

#include "constants.h"
#include "Flash.h"

  include "Flash_mpi.h"

  integer :: currentDataTypes

  integer :: ierr, flashRealExtent, size, stride, flashCont,i
  integer :: exch1, exch2
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: blockID = 1
  integer, dimension(MDIM) :: guard, blkExtent
  ! x exchange
  ! the datatype for the first dimension is contiguous in memory 
  ! mpi_type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype newtype)
  ! left and right, xdir
  

  do i = 1,gr_numDataStruct
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC, gr_gridDataStruct(i))
     guard(:) = blkLimits(LOW,:)-blkLimitsGC(LOW,:)
     blkExtent(:)=blkLimitsGC(HIGH,:)-blkLimitsGC(LOW,:)+1
     
     if(NDIM==1) then
        call MPI_TYPE_CONTIGUOUS(gr_gridDataStructSize(i)*guard(IAXIS), &
                                 FLASH_REAL, gr_exch(i,IAXIS), ierr)
        call MPI_TYPE_COMMIT(gr_exch(i,IAXIS), ierr)
     else
        
        ! y exchange
        ! for up and down, ydir
        call MPI_TYPE_CONTIGUOUS(guard(JAXIS)*gr_gridDataStructSize(i)*&
             blkExtent(IAXIS), & 
             FLASH_REAL, exch2, ierr)
        call MPI_TYPE_COMMIT(exch2,ierr)
        
        
        ! for left and right, xdir
        ! vector is needed to simplify exchange
        ! call MPI_TYPE_VECTOR(count, blocklength, stride, oldtype, newtype)
        call MPI_TYPE_VECTOR(blkExtent(JAXIS), gr_gridDataStructSize(i)*guard(IAXIS), &
             gr_gridDataStructSize(i)*(blkExtent(IAXIS)), &
             FLASH_REAL, exch1, ierr)
        
        call MPI_TYPE_COMMIT(exch1, ierr)
        
        if(NDIM == 2) then
           gr_exch(i,IAXIS)=exch1
           gr_exch(i,JAXIS)=exch2
        else
           
           ! for z exchange
           ! for left and right, xdir
           ! using an hvector 
           call MPI_TYPE_EXTENT(FLASH_REAL, flashRealExtent, ierr)

           call MPI_TYPE_HVECTOR(blkExtent(KAXIS), 1,  &
                blkExtent(IAXIS)*blkExtent(JAXIS)*gr_gridDataStructSize(i)*&
                flashRealExtent, & 
                exch1, gr_exch(i,IAXIS), ierr)
           
           call MPI_TYPE_COMMIT(gr_exch(i,IAXIS),ierr)
           
           
           !for up and down, ydir
           call MPI_TYPE_HVECTOR(blkExtent(KAXIS), 1,  &
                blkExtent(IAXIS)*blkExtent(JAXIS)*gr_gridDataStructSize(i)*&
                flashRealExtent, & 
                exch2, gr_exch(i,JAXIS), ierr)
           
           call MPI_TYPE_COMMIT(gr_exch(i,JAXIS),ierr)
           
           
           ! for back and forth, zdir
           ! contiguous in memory
           size = guard(KAXIS) * blkExtent(IAXIS) * blkExtent(JAXIS) * &
                gr_gridDataStructSize(i)
           
           call MPI_TYPE_CONTIGUOUS(size, FLASH_REAL, gr_exch(i,KAXIS), ierr)
           
           call MPI_TYPE_COMMIT(gr_exch(i,KAXIS),ierr)
           
        end if
     end if
  end do
  
end subroutine gr_createDataTypes


