!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_blockData1DRec
!!
!! NAME
!!
!!  ed_blockData1DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData1DRec (integer (in) :: blockID,
!!                          integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: blockData (:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally
!!  of 1D rectangular grids (cartesian + spherical). The block is specified through
!!  its number ID and the blockData array, which contains the needed data for the
!!  block. The following is computed and stored into specific arrays:
!!
!!     1) the cell Densities
!!     2) the cell Volumes
!!     3) the cell Zbar values
!!     4) the cell center Nele (electron number density) values
!!     5) the cell center Tele (electron temperature) values
!!     6) the cell center gradients of Nele, using adjacent cell center Nele info
!!     7) the cell center gradients of Tele, using adjacent cell center Tele info
!!
!!  The necessary arrays for storage must have been allocated before calling this routine.
!!  No checks are done on the passed cell index limits. This is done before calling this
!!  routine.
!!
!! ARGUMENTS
!!
!!  blockID   : the block ID number
!!  iminBlock : minimum cell i-index limit defining the interior block
!!  imaxBlock : maximum cell i-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  deltaInvI : inverse of the cell's x-dimension
!!  blockData : two-dimensional array containing the block data
!!
!!***

subroutine ed_blockData1DRec (blockID,              &
                              iminBlock, imaxBlock, &
                              iminData , imaxData,  &
                              iminDerv , imaxDerv,  &
                              deltaInvI,            &
                              blockData             )

  implicit none

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: iminDerv , imaxDerv
  real,    intent (in) :: deltaInvI
  real,    intent (in) :: blockData (:,:)

  return
end subroutine ed_blockData1DRec
