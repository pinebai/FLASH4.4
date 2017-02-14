!!****f* source/Grid/Grid_formatNonRep
!!
!! NAME
!!
!!  Grid_formatNonRep
!!
!!
!! SYNOPSIS
!!
!!  Grid_formatNonRep(integer(IN) :: nonrep,
!!                    integer(IN) :: idx,
!!                    character(out) :: str(*))
!!
!!
!! DESCRIPTION
!!
!!  Given a nonreplicated variable array id and index into that array, returns a string name suitable for IO
!!
!!
!! ARGUMENTS
!!  
!!   nonrep - array id
!!   idx - index into array
!!   str - receives string name
!!
!!***
subroutine Grid_formatNonRep(nonrep, idx, str)
    implicit none
    integer, intent(in) :: nonrep, idx
    character(*), intent(out) :: str
    str = ''
end subroutine
