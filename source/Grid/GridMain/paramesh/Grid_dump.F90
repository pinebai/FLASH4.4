!!****if* source/Grid/GridMain/paramesh/Grid_dump
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
!! With paramesh this function doesn not work in parallel, but works
!! only with a single block
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
!!  This routine works with Paramesh only under very special circumstances
!!***

subroutine Grid_dump(var,num, blockID, gcell)

  use Grid_data, ONLY : gr_ilo, gr_ihi, gr_jlo, gr_jhi, &
       gr_klo, gr_khi, gr_iloGC, gr_ihiGC, gr_jloGC, gr_jhiGC, &
       gr_kloGC, gr_khiGC

  use physicaldata, ONLY :unk

!! uncomment next two lines when using with FLASH 2
!!$  use dBase, ONLY : dBaseKeyNumber
#include "Flash.h"

  implicit none

  integer, intent(IN) :: num, blockID
  integer, dimension(num), intent(IN) :: var
  logical, intent(IN) :: gcell
  

  character(len=80) :: ff1
  integer,dimension(4), save :: filecount = 0
  integer :: i,count

  integer,parameter :: bxn=GRID_IHI_GC-GRID_ILO_GC+1
  integer,parameter :: byn=GRID_JHI_GC-GRID_JLO_GC+1 
  integer,parameter :: bzn=GRID_KHI_GC-GRID_KLO_GC+1 
  real,dimension(1,bxn,byn,bzn,1)::scratch

  if(.not. gcell) then
     print '(8F11.6)', unk(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi, blockID)

  else
     print '(16F7.3)', unk(var(1), gr_iloGc:gr_ihiGc, gr_jloGc:gr_jhiGc, gr_kloGc:gr_khiGc, blockID)
  end if



!!*** and  UNCOMMENT THIS 

!!$  var(TEMP_VAR) =dBaseKeyNumber("temp")
!!$  var(GAME_VAR) =dBaseKeyNumber("game")
!!$  var(PRES_VAR) =dBaseKeyNumber("pres")
!!$  var(EINT_VAR) =dBaseKeyNumber("eint")
!!$  var(VELZ_VAR) =dBaseKeyNumber("velz")
!!$  var(VELY_VAR) =dBaseKeyNumber("vely")
!!$  var(VELX_VAR) =dBaseKeyNumber("velx")
!!$  var(DENS_VAR) =dBaseKeyNumber("dens")
!!$  var(GAMC_VAR) =dBaseKeyNumber("gamc")
!!$  var(ENER_VAR) =dBaseKeyNumber("ener")
 

  !print '(8F11.6)', unk(var(1), gr_ilo:gr_ihi, gr_jlo:gr_jhi, gr_klo:gr_khi, blockID) 
 
  ff1 = "FL2"//char(48+filecount(4))//char(48+filecount(3))//&
       char(48+filecount(2))//char(48+filecount(1))
  !print*,'filecount',filecount
  filecount(1) = filecount(1) + 1
  do i = 1,3
     if(filecount(i)==10)then
        filecount(i) = 0
        filecount(i+1)=filecount(i+1)+1
     end if
  end do
  count = 8*bxn*byn*bzn
  if((filecount(3) == 0).and.filecount(2) == 0) then
  
!!! create the data structure for outputting data

     open(41,file=ff1,access='direct',form='unformatted',recl=count)
     do i = 1,num
        !scratch(1,:,:,:,1)=unk(var(i),:,:,:,blockID)
        !write(41,rec=i)scratch
     end do
     close(41)
  end if
  return
end subroutine Grid_dump
