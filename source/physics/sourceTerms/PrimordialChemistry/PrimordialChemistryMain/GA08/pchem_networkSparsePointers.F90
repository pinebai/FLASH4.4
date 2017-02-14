!!****ih* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_networkSparsePointers
!!
!! NAME
!!
!!  pchem_networkSparsePointers
!!
!! SYNOPSIS
!!
!!  call pchem_networkSparsePointers(integer, intent(INOUT)    :: iloc(np),
!!                                  integer, intent(INOUT)    :: jloc(np),
!!                                  integer, intent(OUT)      :: nzo,
!!                                  integer, intent(IN)       :: np)
!!
!!
!! DESCRIPTION
!!
!!  this routine is supposed to build the nonzero locations for networkSparseJakob
!!
!!  input is the integer arrays iloc and jloc, both of dimension np, that on
!!  output contain nzo matrix element locations
!!
!!  But is not, just a dummy because our chemistry is not a sparse matrix
!!
!!  called externally from pchem_burner as routine 'bjackob'
!!
!!  ARGUMENTS
!!
!!  iloc  -- integer, array of size np
!!  jloc  -- integer, array of size np
!!  np    -- integer, size of arrays
!!  nzo   -- integer, number of nonzero matrix element locations
!!
!!***


subroutine pchem_networkSparsePointers(iloc,jloc,nzo,np)

  implicit none
  integer, intent(IN)  :: iloc(*),jloc(*),np
  integer, intent(OUT) :: nzo

  return
end subroutine pchem_networkSparsePointers

