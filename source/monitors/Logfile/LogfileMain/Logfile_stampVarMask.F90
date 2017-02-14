!!****if* source/monitors/Logfile/LogfileMain/Logfile_stampVarMask
!!
!! NAME
!!
!!  Logfile_stampVarMask
!!
!! SYNOPSIS
!!
!!  Logfile_stampVarMask(logical(IN),dimension(:) :: unkvarmask(NUNK_VARS),
!!                       logical(IN)              :: willcalleos,
!!                       character(len=*)(IN)  :: tag,
!!                       character(len=*)(IN)  :: masktag)
!!
!! DESCRIPTION
!!
!!  Generates a line in the log file that displays an unk mask (as used for
!!  masking of variables for Grid_fillGuardCells calls) and the state of the
!!  willcalleos flag.
!!
!! ARGUMENTS
!!
!!   unkvarmask : A logicall array. For each true element, the logged string has the
!!                first letter of the corresponding variable. For each false element,
!!                the logged string has the character '-'.
!!
!!   willcalleos : A flag. if true, the line logged has the characters '+Eos' appended.
!!
!!   tag : A string that will appear at the beginning of the log line.
!!
!!   masktag : A string that will appear right before the logged mask string, which
!!             follows it separated by '='.
!!
!! NOTES
!!
!!  This subroutine is mostly useful for debugging purposes.
!!
!!***

#include "Flash.h"
#include "constants.h"

subroutine Logfile_stampVarMask(unkVarMask, willCallEos, tag, maskTag)

  use Logfile_interface, ONLY : Logfile_stamp
  use Logfile_data, ONLY :  log_enableGcMaskLogging
  implicit none

  logical, intent(in) :: unkVarMask(:)
  logical, intent(in) :: willCallEos
  character(len=*),intent(in)    :: tag, maskTag

  character(len=size(unkVarMask)+5) :: string
  character(len=1) :: name


  if (.NOT. log_enableGcMaskLogging) return

  call log_formatMaskString(unkVarMask, string, willCallEos)
  call Logfile_stamp( maskTag//'='//string, tag)

contains

  subroutine log_formatMaskString(unkVarMask, string, willCallEos)
    use Simulation_interface, ONLY : Simulation_mapIntToStr
    implicit none

    logical, intent(in) :: unkVarMask(:)
    character(len=size(unkVarMask)+5), intent(out) :: string
    logical, intent(in) :: willCallEos
    integer :: i,j

    string = ' '
    do i=1, NUNK_VARS
       if (unkVarMask(i)) then
          call Simulation_mapIntToStr(i, name, MAPBLOCK_UNK)
          string(i:i) = name
       else
          string(i:i) = '-'
       end if
    end do
#if NFACE_VARS > 0
    if (size(unkVarMask)>NUNK_VARS) then
       string(i:i) = '|'
       i = i+1
    do i=NUNK_VARS+2,1+min(size(unkVarMask),NUNK_VARS+MDIM*NFACE_VARS)
       j = mod(i - 2 - NUNK_VARS,NFACE_VARS) + 1
       if (unkVarMask(i-1)) then
          call Simulation_mapIntToStr(j, name, MAPBLOCK_FACES)
          name(1:1) = char(ieor(32,ichar(name(1:1))))
          string(i:i) = name
       else
          string(i:i) = '-'
       end if
    end do
    end if
#endif
    if (willCallEos) then
       string(i:i+3) = '+Eos'
    else
       string(i:i+3) = '    '
    end if
  end subroutine log_formatMaskString

end subroutine Logfile_stampVarMask
