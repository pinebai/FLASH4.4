!!****if* source/Grid/GridMain/UG/Grid_dump
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
!! Binary dumps the physical variables specified in "var". Can be
!! done from anywhere in the code, and is useful for diagnostic
!! purposes. It works only with single block per processor mapping.
!!
!! ARGUMENTS
!!
!!  var :: 1D integer array containing the names of the variables to
!!         be dumped (as defined in Flash.h)
!!  num :: number of variables being dumped.
!!  blockID :: number of block to dump.  In UG always send in 1
!!  gcell :: indicated whether to include guardcells in the dump.
!!           if gcell = .true. guardcells are included in the dump
!!                      which means that at the edges, the values are
!!                      duplicated.
!!           if gcell = .false. the domain is dumped as though it was
!!                       a single block
!!
!! EXAMPLES
!!  if num = 3, and var(1)=DENS_VAR, var(2)=PRES_VAR and var(3)=TEMP_VAR
!!  then the current values of density, pressure and temperature will be
!!  dumped.
!!  To explain the use of gcells, consider a global domain with 8x8 points
!!  mapped on 2x2 processors. each processor has blocks of size
!!  4x4.If there are 2 guard cells along each dimension, then the
!!  block size including guardcells is 8x8 and the distribution on
!!  four processors is as shown below ("*"
!!  are the interior points and the )"o" are guard cells.
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oo****oo        oo****oo
!!             oooooooo        oooooooo
!!             oooooooo        oooooooo
!!
!!  If gcell is true then dump is of size 16x16 and looks like
!!
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oo****oooo****oo
!!                oooooooooooooooo
!!                oooooooooooooooo
!!
!!  and if gcell is false the dump is of size 8x8 and looks like
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!                ********
!!    
!!
!!***



subroutine Grid_dump(var,num,blockID,gcell)
  
  use physicalData, ONLY : unk
  use Grid_data, ONLY : gr_ilo,gr_ihi,gr_jlo,gr_jhi,gr_klo,gr_khi
  use Grid_data, ONLY : gr_iloGc,gr_ihiGc,gr_jloGc,gr_jhiGc,gr_kloGc,gr_khiGc
  use Grid_data, ONLY : gr_axisNumProcs,gr_axisMe, gr_meshComm

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
  
  integer, parameter   :: i16=SELECTED_INT_KIND(16)
  integer(kind=i16)    :: disp, offset
  
  integer :: FLASH_ARR
  integer,dimension(MDIM) :: gsize,lsize,starts

  ff1 = "FL3"//char(48+filecount(4))//char(48+filecount(3))//&
       char(48+filecount(2))//char(48+filecount(1))
  !print*,'filecount',filecount
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
        lsize(1)=gr_ihiGC-gr_iloGc+1
        lsize(2)=gr_jhiGC-gr_jloGc+1
        lsize(3)=gr_khiGC-gr_kloGc+1
     else
        lsize(1)=gr_ihi-gr_ilo+1
        lsize(2)=gr_jhi-gr_jlo+1
        lsize(3)=gr_khi-gr_klo+1
     end if
 
     gsize=lsize*gr_axisNumProcs
     starts = gr_axisMe*lsize

     call MPI_TYPE_CREATE_SUBARRAY(MDIM,gsize,lsize,starts,&
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
           call MPI_File_write_all(fh,unk(var(i),:,:,:,1),&
                count,FLASH_REAL,status,ierr)
        end do
     else
        do i = 1,num
           call MPI_File_write_all&
                (fh,unk(var(i),gr_ilo:gr_ihi,gr_jlo:gr_jhi,gr_klo:gr_khi,1),&
                count,FLASH_REAL,status,ierr)
        end do
     end if
     call MPI_File_close(fh,ierr)
     
  end if
  return
end subroutine Grid_dump
