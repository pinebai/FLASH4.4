!!****if* source/Grid/GridMain/UG/UGReordered/Grid_dump
!!
!! NAME
!!  Grid_dump
!!
!! SYNOPSIS
!!
!!  call Grid_dump(integer(IN) :: var(num),
!!                 integer(IN) :: num,
!!                 integer(IN) :: blockID,
!!                 logical(IN) :: gcell)
!!  
!! DESCRIPTION 
!!  
!! Dumps the variables specified in "var" to a file. Can be done from 
!! anywhere in the code, and is useful for diagnostic purposes
!! This function can only be used on a single block at a time
!!  
!! ARGUMENTS 
!!
!!  var :: array containing the indices of the variables to be dumped
!!  num :: number of variables being dumped.
!!  blockID :: local number of block to be dumped
!!  gcell :: indicates whether to include guardcells in the dump.
!!             
!! EXAMPLE
!!  
!!  num = 3  !dumping 3 variables
!!  var(1) = DENS_VAR
!!  var(2) = PRES_VAR
!!  var(3) = TEMP_VAR
!!  blockID = 1  ! local block number
!!  gcell = .false.
!!
!!  call Grid_dump(var, num, blockID, gcell)
!!  
!!  will dump the interior cell values of density, pressure and temperature
!!  for local block number 1.
!!
!! NOTES
!!  DENS_VAR, PRES_VAR, TEMP_VAR etc are #defined values in Flash.h
!!  indicating the index in the physical data array.
!!  The routine calling Grid_dump will need to include Flash.h 
!!
!!***


subroutine Grid_dump(var,num,blockID,gcell)
  
  use physicalData, ONLY : unk
  use Grid_data, ONLY : gr_lIndexSize,gr_axisMe,gr_axisNumProcs,gr_guard, gr_meshComm
  use Grid_data, ONLY : gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi, &
       gr_iloGc, gr_ihiGc, gr_jloGc, gr_jhiGc, gr_kloGc, gr_khiGc



  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(IN) :: num, blockID
  integer, dimension(num), intent(IN) :: var
  logical, intent(IN) :: gcell

  character(len=80) :: ff1
  integer,dimension(4), save :: filecount = 0

  integer :: i, count, ierr, siz, fh
  
  integer,dimension(MPI_STATUS_SIZE) :: status
  
  integer(kind=MPI_OFFSET_KIND)::disp, offset
  
  integer :: FLASH_ARR
  integer, dimension(MDIM) :: gsize,lsize,starts

  
  ff1 = "FL3"//char(48+filecount(4))//char(48+filecount(3))//&
       char(48+filecount(2))//char(48+filecount(1))
  filecount(1) = filecount(1) + 1
  do i = 1,4
     if(filecount(i)==10)then
        filecount(i) = 0
        filecount(i+1)=filecount(i+1)+1
     end if
  end do
  
  if((filecount(4) == 0))then !.and.filecount(2) == 0) then
  
!!! create the data structure for outputting data
     if(gcell) then
        lsize = gr_lIndexSize/gr_axisNumProcs+gr_guard*2
     else
        lsize = gr_lIndexSize/gr_axisNumProcs
     end if
     gsize = lsize*gr_axisNumProcs
     starts = gr_axisMe*lsize
     
     call MPI_TYPE_CREATE_SUBARRAY(NDIM,gsize,lsize,starts,&
          MPI_ORDER_FORTRAN,FLASH_REAL,FLASH_ARR,ierr)
     call MPI_TYPE_COMMIT(FLASH_ARR,ierr)
  
  
     call MPI_Type_size(FLASH_ARR,count,ierr)
     call MPI_Type_size(FLASH_REAL,siz,ierr)
  
     count = count/siz
     disp = 0
     offset = 0
  
     call MPI_File_open(gr_meshComm,ff1,&
          MPI_MODE_CREATE+MPI_MODE_RDWR,MPI_INFO_NULL,fh,ierr)
     call MPI_File_set_view(fh,disp,MPI_BYTE,FLASH_ARR,"native",&
          MPI_INFO_NULL,ierr)
     if(gcell) then
        do i = 1,num
!!           gcellarray(:,:,:) = unk(var(i),:,:,:,1)
           call MPI_File_write_all(fh,unk(:,:,:,var(i),1),&
                count,FLASH_REAL,status,ierr)
        end do
     else
        do i = 1,num
!!           nogcellarray(:,:,:) = 
           !!unk(var(i),gr_ilo:gr_ihi,gr_jlo:gr_jhi,gr_klo:gr_khi,1)
           call MPI_File_write_all&
                (fh,unk(gr_ilo:gr_ihi,gr_jlo:gr_jhi,gr_klo:gr_khi,var(i),1),&
                count,FLASH_REAL,status,ierr)
        end do
     end if
     call MPI_File_close(fh,ierr)
     
  end if
  return
end subroutine Grid_dump
