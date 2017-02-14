!!****if* source/diagnostics/ProtonImaging/localAPI/pi_blockData3DRec
!!
!! NAME
!!
!!  pi_blockData3DRec
!!
!! SYNOPSIS
!!
!!  call pi_blockData3DRec (integer (in) :: iminBlock,
!!                          integer (in) :: imaxBlock,
!!                          integer (in) :: jminBlock,
!!                          integer (in) :: jmaxBlock,
!!                          integer (in) :: kminBlock,
!!                          integer (in) :: kmaxBlock,
!!                          real    (in) :: deltaInvX,
!!                          real    (in) :: deltaInvY,
!!                          real    (in) :: deltaInvZ,
!!                          real    (in) :: blockData (:,:,:,:))
!!
!! DESCRIPTION
!!
!!  Computes cell data for one specific block for those geometries consisting formally of
!!  3D rectangular grids (cartesian). The block is specified through its number ID and the
!!  blockData array, which contains the needed data for the block. The following is computed
!!  and stored into specific arrays:
!!
!!     1) the cell center electric field components (Ex,Ey,Ex) in Gauss
!!     2) the cell center magnetic flux density components (Bx,By,Bx) in Gauss
!!     3) the cell center curl of the magnetic flux density components (CurlBx,CurlBy,CurlBx)
!!        in units of Gauss / cm
!!     4) the cell boundary indicator
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
!!  deltaInvX : inverse of the cell's x-dimension
!!  deltaInvY : inverse of the cell's y-dimension
!!  deltaInvZ : inverse of the cell's z-dimension
!!  blockData : four-dimensional array containing the block data
!!
!!***

subroutine pi_blockData3DRec (iminBlock, imaxBlock,            &
                              jminBlock, jmaxBlock,            &
                              kminBlock, kmaxBlock,            &
                              deltaInvX, deltaInvY, deltaInvZ, &
                              blockData                        )

  implicit none

  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
  real,    intent (in) :: blockData (:,:,:,:)

  return
end subroutine pi_blockData3DRec
