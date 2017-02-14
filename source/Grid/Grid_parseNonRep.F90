!!****f* source/Grid/Grid_parseNonRep
!!
!! NAME
!!
!!  Grid_parseNonRep
!!
!!
!! SYNOPSIS
!!
!!  Grid_parseNonRep(character(IN) :: strlwr(*),
!!                     integer(OUT):: nonrep,
!!                     integer(OUT):: idx)
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
    character(*), intent(in) :: strlwr
    integer, intent(out) :: nonrep, idx
    ! signal no match
    nonrep = 0
    idx = 0
end subroutine
