!!****if* source/Grid/GridMain/Grid_parseNonRep
!!
!! NAME
!!
!!  Grid_parseNonRep
!!
!!
!! SYNOPSIS
!!
!!  Grid_parseNonRep(character(len=*)(IN) :: strlwr,
!!                   integer(OUT) :: nonrep,
!!                   integer(OUT) :: idx)
!!
!!
!! DESCRIPTION
!!
!!  Given the string name of a nonreplicated variable array element, returns the array id and index.
!!
!!
!! ARGUMENTS
!!  
!!   strlwr - the all lowercase string name of the array element, MUST BE ALL lowercase
!!   nonrep - receives the array id
!!   idx - receives the array index
!!
!!***
subroutine Grid_parseNonRep(strlwr, nonrep, idx)
    implicit none
#include "Flash.h"
    character(*), intent(in) :: strlwr
    integer, intent(out) :: nonrep, idx

    character(120) :: namef_flat = NONREP_NAMEF_FLAT_LWR ! this should be a parameter, but it was causing nag to choke
    integer, parameter :: namef_start(1:NONREP_COUNT+1) = NONREP_NAMEF_START

    integer :: i, nf_off, nf_len, c
    character(1) :: c1
    logical :: match
    
    do nonrep = 1, NONREP_COUNT
        nf_off = namef_start(nonrep) - 1
        nf_len = namef_start(nonrep+1) - nf_off - 1
        match = .false.
        idx = 0
        if(nf_len .eq. len(strlwr)) then
            match = .true.
            do i = 1, nf_len
                c1 = strlwr(i:i)
                c = ichar(c1)
                c1 = namef_flat(nf_off+i:nf_off+i)
                if(ichar(c1) .eq. c) cycle
                if(c1 .eq. "?" .and. ichar("0") .le. c .and. c .le. ichar("9")) then
                    idx = 10*idx + c - ichar("0")
                    cycle
                end if
                match = .false.
                exit
            end do
        end if
        if(match) return
    end do
    ! signal no match
    nonrep = 0
    idx = 0
end subroutine
