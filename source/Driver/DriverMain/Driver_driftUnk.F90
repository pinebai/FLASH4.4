!!****if* source/Driver/DriverMain/Driver_driftUnk
!!
!! NAME
!!  Driver_driftUnk
!!
!! DESCRIPTION
!!
!!  Compute one hash per unk variable over all registered block hashes.  Vars with
!!  hashes that have changed since the previous call will be logged.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  flags: bitmask of options (for now there is only one flag)
!!    - DRIFT_NO_PARENTS: if present then only leaf blocks are included
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Driver_driftUnk(src_file, src_line, flags)
  use Driver_data, only: dr_globalMe, dr_nstep, dr_driftFd, &
    dr_driftUnk, dr_driftBlk, dr_driftInst
#ifdef FLASH_GRID_PARAMESH
  use tree, only: nodetype
#endif
  implicit none
  
  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in), optional :: flags

#if DRIFT_ENABLE
  integer :: flgs
  integer(kind=selected_int_kind(18)) :: h
  integer :: b, v
  character(len=30) :: logfile
  character(len=4) :: vname
  character(len=10) :: num1, num2
  character(len=50) :: srcloc
  
  integer, save :: last_inst = 1
  
  flgs = 0
  if(present(flags)) flgs = flags
  
  write(num1,'(I10)') dr_globalMe
  write(logfile,'("drift.",A,".log")') num1(bigdig(dr_globalMe,10):10)
  open(unit=dr_driftFd,file=logfile,position='APPEND')
  
  write(num1,'(I10)') dr_nstep
  write(dr_driftFd,'("step=",A)') num1(bigdig(dr_nstep,10):10)
  
  write(num1,'(I10)') src_line
  srcloc = trim(src_file)//':'//num1(bigdig(src_line,10):10)
  write(dr_driftFd,'("from=",A)') trim(srcloc)
  write(num1,'(I10)') last_inst
  write(num2,'(I10)') dr_driftInst-1
  write(dr_driftFd,'("unks inst=",A," to ",A)') &
    num1(bigdig(last_inst,10):10), num2(bigdig(dr_driftInst-1,10):10)
  last_inst = dr_driftInst
  
  do v=UNK_VARS_BEGIN, UNK_VARS_END
    h = 0
    do b=1, MAXBLOCKS
#ifdef FLASH_GRID_PARAMESH
      if(iand(flgs,DRIFT_NO_PARENTS)==0 .or. nodetype(b)==LEAF) then
#else
      if(.true.) then
#endif
        h = h + dr_driftBlk(v,b)
        h = h + ishft(h,10)
        h = ieor(h, ishft(h,-6))
      end if
    end do
    if(h /= dr_driftUnk(v)) then
      call Simulation_mapIntToStr(v,vname,MAPBLOCK_UNK)
      write(dr_driftFd,'(" ",A4," ",z16.16)') vname, h
      dr_driftUnk(v) = h
    end if
  end do
  
  write(dr_driftFd,'(A)') ''
  close(dr_driftFd)
#endif
contains
  elemental function bigdig(x,w) result(n)
    integer, intent(in) :: x, w
    integer :: n
    n = w - int(log(real(max(1,x)))/log(10.)*(1.0 + 2*epsilon(1.0)))
  end function
end subroutine
