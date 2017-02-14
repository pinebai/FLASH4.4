!!****f* source/Grid/Grid_getVarNonRep
!!
!! NAME
!!
!!  Grid_getVarNonRep
!!
!! SYNOPSIS
!!
!!  call Grid_getVarNonRep(integer(in) :: mapblock, 
!!                         integer(in) :: var, 
!!                         integer(out) :: nonrep,
!!                         integer(out) :: locidx) 
!!
!! DESCRIPTION
!!
!!  Given the index of a variable, this function will determine the id of the non-mesh-replicated
!!  variable array the given unk variable belongs to.
!!
!!
!! ARGUMENTS
!! 
!!  mapblock - integer specifier
!!  var - integer index of variable in mapblock
!!  nonrep - 0 if this variable does not belong to a nonreplicated array, otherwise the 1-based index of that array
!!  locidx - 0 if this variable does not belong to a nonreplicated array, otherwise the 1-based local index within
!!           the array this variable corresponds to
!!
!! NOTES
!!
!!***

subroutine Grid_getVarNonRep(mapblock, var, nonrep, locidx)
   implicit none
   integer, intent(in) :: mapblock, var
   integer, intent(out) :: nonrep
   integer, intent(out), optional :: locidx

   nonrep = 0
   if(present(locidx)) locidx = 0
end subroutine
