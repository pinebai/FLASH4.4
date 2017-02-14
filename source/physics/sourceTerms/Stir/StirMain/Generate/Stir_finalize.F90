!!****if* source/physics/sourceTerms/Stir/StirMain/Generate/Stir_finalize
!!
!! NAME
!!
!!  Stir_finalize
!!
!! SYNOPSIS
!!
!!  call Stir_finalize()
!!
!! ARGUMENTS
!!
!!  none
!!
!! DESCRIPTION
!!
!!  Clean up the Stir unit
!!
!!***

subroutine Stir_finalize()
  use Stir_data, ONLY : st_randseed, st_reproducible, &
       st_saveReproducible, st_randomSaveUnit
  implicit none
  deallocate(st_randseed)
  if(st_reproducible.or.st_saveReproducible)&
       close(unit=st_randomSaveUnit)
  return
end subroutine Stir_finalize
