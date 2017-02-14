!!****if* source/Driver/DriverMain/Driver_driftBlock
!!
!! NAME
!!  Driver_driftBlock
!!
!! DESCRIPTION
!!
!!  Register the content hash of the provided block data, if it is different
!!  from the last time this block was hashed then log it to file.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  blk: processor local block number
!!  ptr: block data to hash, all values will be hashed so this should not include guard cells
!!  gds: grid data struct type of data
!!
!!***
#include "Flash.h"
#include "constants.h"

subroutine Driver_driftBlock(src_file, src_line, blk, ptr, gds)
  use Driver_interface, only: Driver_dbgBreak
  use Driver_data, only: dr_globalMe, dr_nstep, &
    dr_driftFd, dr_driftBlk, dr_driftTruncMantissa, &
    dr_driftTuples, dr_driftInst, dr_driftVerboseInst
#ifdef FLASH_GRID_PARAMESH
  use tree, only: nodetype
#endif
  implicit none

  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in) :: blk
  real, intent(in) :: ptr(:,:,:,:)
  integer, intent(in) :: gds

#if DRIFT_ENABLE
  integer :: i,j,k,v,c
  integer(kind=selected_int_kind(18)) :: h(UNK_VARS_BEGIN:UNK_VARS_END)
  integer(kind=selected_int_kind(2)) :: buf(8) ! bytes to hold floating point bits
  integer(kind=selected_int_kind(2)) :: mm(8) ! mantissa mask (bitmask used for rounding)
  
  character(len=30) :: logfile
  character(len=4) :: vname
  character(len=10) :: num1, num2, num3
  character(len=50) :: srcloc
  integer :: mapb
  
  if(gds == CENTER) then
    mapb = MAPBLOCK_UNK
  else
    return
  end if
  
  ! behold...
  mm = not(ieor(transfer(1.0, mm), transfer(1.0 + epsilon(1.0)*(2**dr_driftTruncMantissa-1), mm)))
  
  h(:) = 0
  do k=lbound(ptr,4), ubound(ptr,4)
    do j=lbound(ptr,3), ubound(ptr,3)
      do i=lbound(ptr,2), ubound(ptr,2)
        do v=UNK_VARS_BEGIN, UNK_VARS_END
          buf = transfer(ptr(v,i,j,k),buf)
          do c=1, 8
            h(v) = h(v) + iand(buf(c), mm(c))
            h(v) = h(v) + ishft(h(v),10)
            h(v) = ieor(h(v), ishft(h(v),-6))
          end do
        end do
      end do
    end do
  end do
  
  if(dr_driftVerboseInst > 0 .and. dr_driftInst >= dr_driftVerboseInst) then
    if(any(h(:) /= dr_driftBlk(:,blk))) then
      write(num1,'(I10)') src_line
      srcloc = trim(src_file)//':'//num1(bigdig(src_line,10):10)
      
      write(num1,'(I10)') dr_globalMe
      if(dr_driftTuples) then
        write(logfile,'("drift.tups.",A,".log")') num1(bigdig(dr_globalMe,10):10)
      else
        write(logfile,'("drift.",A,".log")') num1(bigdig(dr_globalMe,10):10)
      end if
      open(unit=dr_driftFd,file=logfile,position='APPEND')
      
      write(num1,'(I10)') dr_driftInst
      
      if(.not. dr_driftTuples) then
        write(dr_driftFd,'("inst=",A)') num1(bigdig(dr_driftInst,10):10)
        write(num2,'(I10)') dr_nstep
        write(dr_driftFd,'("step=",A)') num2(bigdig(dr_nstep,10):10)
        write(dr_driftFd,'("src=",A)') trim(srcloc)
        write(num2,'(I10)') blk
        write(dr_driftFd,'("blk=",A)') num2(bigdig(blk,10):10)
      end if
      
      do v=UNK_VARS_BEGIN, UNK_VARS_END
        if(h(v) /= dr_driftBlk(v,blk)) then
          call Simulation_mapIntToStr(v, vname, mapb)
          if(dr_driftTuples) then
            write(num2,'(I10)') dr_nstep
            write(num3,'(I10)') blk
            write(dr_driftFd,'(A,",",A,",''",A,"'',",A,",",I1,",''",A,"'',0x",z16.16)') &
              num1(bigdig(dr_driftInst,10):10), &
              num2(bigdig(dr_nstep,10):10), &
              trim(srcloc), &
              num3(bigdig(blk,10):10), &
#ifdef FLASH_GRID_PARAMESH
              count((/nodetype(blk)==LEAF/)), &
#else
              1, &
#endif
              trim(vname), &
              h(v)
          else
            write(dr_driftFd,'(" ",A4," ",z16.16)') vname, h(v)
          end if
          dr_driftBlk(v,blk) = h(v)
        end if
      end do
      
      if(.not. dr_driftTuples) write(dr_driftFd,'(A)') ''
      close(dr_driftFd)
    end if
  else
    dr_driftBlk(:,blk) = h(:)
  end if
  
  dr_driftInst = dr_driftInst + 1
#endif
contains
  elemental function bigdig(x,w) result(n)
    integer, intent(in) :: x, w
    integer :: n
    n = w - int(log(real(max(1,x)))/log(10.)*(1.0 + 2*epsilon(1.0)))
  end function
end subroutine Driver_driftBlock
