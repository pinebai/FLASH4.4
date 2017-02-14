!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_updatePlasma
!!
!! NAME
!!
!!  ed_updatePlasma
!!
!! SYNOPSIS
!!
!!  call ed_updatePlasma (integer (in) :: blockCount, 
!!                        integer (in) :: blockList (:))
!!
!! DESCRIPTION
!!
!!  Updates the physical variables of all cells on all blocks owned by the
!!  current processor. This is done by calling Eos.
!!
!! ARGUMENTS
!!
!!  blockCount : number of blocks on current processor
!!  blockList  : all block ID numbers
!!
!! NOTES
!!
!!***

subroutine ed_updatePlasma (blockCount,blockList,scaleFact)

  implicit none
  
  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,OPTIONAL, intent (in) :: scaleFact

  return
end subroutine ed_updatePlasma
