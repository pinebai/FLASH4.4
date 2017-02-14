!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_blockData3DRec
!!
!! NAME
!!
!!  ed_blockData3DRec
!!
!! SYNOPSIS
!!
!!  call ed_blockData3DRec (integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: jminBlock,
!!                          integer (in) :: jmaxBlock,
!!                          integer (in) :: kminBlock,
!!                          integer (in) :: kmaxBlock,
!!                          integer (in) :: iminData,
!!                          integer (in) :: imaxData,
!!                          integer (in) :: jminData,
!!                          integer (in) :: jmaxData,
!!                          integer (in) :: kminData,
!!                          integer (in) :: kmaxData,
!!                          integer (in) :: iminDerv,
!!                          integer (in) :: imaxDerv,
!!                          integer (in) :: jminDerv,
!!                          integer (in) :: jmaxDerv,
!!                          integer (in) :: kminDerv,
!!                          integer (in) :: kmaxDerv,
!!                          real    (in) :: deltaI,
!!                          real    (in) :: deltaJ,
!!                          real    (in) :: deltaK,
!!                          real    (in) :: deltaInvI,
!!                          real    (in) :: deltaInvJ,
!!                          real    (in) :: deltaInvK,
!!                          real    (in) :: blockData (:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally of
!!  3D rectangular grids (cartesian). The block is specified through its number ID and the
!!  blockData array, which contains the needed data for the block. The following is computed
!!  and stored into specific arrays:
!!
!!     1) the cell Densities
!!     2) the cell Zbar values
!!     3) the cell center Nele (electron number density) values
!!     4) the cell center Tele (electron temperature) values
!!     5) the cell center gradients of Nele, using adjacent cell center Nele info
!!     6) the cell center gradients of Tele, using adjacent cell center Tele info
!!
!!  The necessary arrays for storage must have been allocated before calling this routine.
!!  No checks are done on the passed cell index limits. This is done before calling this
!!  routine.
!!
!! ARGUMENTS
!!
!!  iminBlock : minimum cell i-index limit defining the interior block
!!  imaxBlock : maximum cell i-index limit defining the interior block
!!  jminBlock : minimum cell j-index limit defining the interior block
!!  jmaxBlock : maximum cell j-index limit defining the interior block
!!  kminBlock : minimum cell k-index limit defining the interior block
!!  kmaxBlock : maximum cell k-index limit defining the interior block
!!  iminData  : minimum cell i-index limit needed for evaluating Nele and Tele values
!!  imaxData  : maximum cell i-index limit needed for evaluating Nele and Tele values
!!  jminData  : minimum cell j-index limit needed for evaluating Nele and Tele values
!!  jmaxData  : maximum cell j-index limit needed for evaluating Nele and Tele values
!!  kminData  : minimum cell k-index limit needed for evaluating Nele and Tele values
!!  kmaxData  : maximum cell k-index limit needed for evaluating Nele and Tele values
!!  iminDerv  : minimum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  imaxDerv  : maximum cell i-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  jminDerv  : minimum cell j-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  jmaxDerv  : maximum cell j-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  kminDerv  : minimum cell k-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  kmaxDerv  : maximum cell k-index limit needed for evaluating gradients (derivatives) of Nele and Tele
!!  deltaI    : the cell's x-dimension
!!  deltaJ    : the cell's y-dimension
!!  deltaK    : the cell's z-dimension
!!  deltaInvI : inverse of the cell's x-dimension
!!  deltaInvJ : inverse of the cell's y-dimension
!!  deltaInvK : inverse of the cell's z-dimension
!!  blockData : four-dimensional array containing the block data
!!
!!***

subroutine ed_blockData3DRec (iminBlock, imaxBlock,            &
                              jminBlock, jmaxBlock,            &
                              kminBlock, kmaxBlock,            &
                              iminData , imaxData,             &
                              jminData , jmaxData,             &
                              kminData , kmaxData,             &
                              iminDerv , imaxDerv,             &
                              jminDerv , jmaxDerv,             &
                              kminDerv , kmaxDerv,             &
                              deltaI   , deltaJ   , deltaK,    &
                              deltaInvI, deltaInvJ, deltaInvK, &
                              blockData                        )

  implicit none

  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: jminData , jmaxData
  integer, intent (in) :: kminData , kmaxData
  integer, intent (in) :: iminDerv , imaxDerv
  integer, intent (in) :: jminDerv , jmaxDerv
  integer, intent (in) :: kminDerv , kmaxDerv
  real,    intent (in) :: deltaI   , deltaJ   , deltaK
  real,    intent (in) :: deltaInvI, deltaInvJ, deltaInvK
  real,    intent (in) :: blockData (:,:,:,:)

  return
end subroutine ed_blockData3DRec
