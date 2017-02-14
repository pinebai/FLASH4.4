!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_blockData2DRec
!!
!! NAME
!!
!!  ed_blockData2DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData2DRec (integer (in) :: blockID,
!!                          integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: jminBlock,
!!                          integer (in) :: jmaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: jminData,
!!                          integer (in) :: jmaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          integer (in) :: jminDerv,
!!                          integer (in) :: jmaxDerv,
!!                          real    (in) :: deltaI,
!!                          real    (in) :: deltaJ,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: deltaInvJ,
!!                          real    (in) :: blockData (:,:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally
!!  of 2D rectangular grids (cartesian + cylindrical). The block is specified through
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
!!  jminBlock : minimum cell j-index limit defining the interior block
!!  jmaxBlock : maximum cell j-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating Nele and Tele values
!!  jminData  : minimum cell j-index limit needed for evaluating Nele and Tele values
!!  jmaxData  : maximum cell j-index limit needed for evaluating Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  jminDerv  : minimum cell j-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  jmaxDerv  : maximum cell j-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  deltaI    : the cell's x-dimension
!!  deltaJ    : the cell's y-dimension
!!  deltaInvI : inverse of the cell's x-dimension
!!  deltaInvJ : inverse of the cell's y-dimension
!!  blockData : three-dimensional array containing the block data
!!
!!***

subroutine ed_blockData2DRec (blockID,              &
                              iminBlock, imaxBlock, &
                              jminBlock, jmaxBlock, &
                              iminData , imaxData,  &
                              jminData , jmaxData,  &
                              iminDerv , imaxDerv,  &
                              jminDerv , jmaxDerv,  &
                              deltaI   , deltaJ,    &
                              deltaInvI, deltaInvJ, &
                              blockData             )

  implicit none

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: jminData , jmaxData
  integer, intent (in) :: iminDerv , imaxDerv
  integer, intent (in) :: jminDerv , jmaxDerv
  real,    intent (in) :: deltaI   , deltaJ
  real,    intent (in) :: deltaInvI, deltaInvJ
  real,    intent (in) :: blockData (:,:,:)

  return
end subroutine ed_blockData2DRec
