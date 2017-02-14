!!****f* source/Grid/Grid_subcellGeometry
!!
!! NAME
!!
!!  Grid_subcellGeometry
!!
!! SYNOPSIS
!!
!!  call Grid_subcellGeometry(:: nsubi,
!!                             :: nsubj,
!!                             :: nsubk,
!!                            real(in) :: dvcell,
!!                            real, intent(OUT)  :: dvsub(nsubI,nsubJ),
!!                            real,OPTIONAL(in) :: xl,
!!                            real,OPTIONAL(in) :: xr,
!!                            real,OPTIONAL(in) :: yl,
!!                            real,OPTIONAL(in) :: yr,
!!                            integer,OPTIONAL(in) :: pos(*),
!!                            integer,OPTIONAL(in) :: blockid)
!!
!! DESCRIPTION
!!
!!  Geometrically correct computation of the volumnes of subcells.
!!
!!  Currently only implemented for CARTESIAN and CYLINDRICAL geometries.
!!
!! ARGUMENTS
!!
!!   nsubi :     Number of subcell lengths per cell length in X direction.
!!
!!   nsubj :     Number of subcell lengths per cell length in Y direction.
!!
!!   nsubk :     Number of subcell lengths per cell length in Z direction.
!!
!!   dvcell :    Volume of the whole cell.
!!
!!   dvsub :     Volumes of subcells. Note that this is a 2-dimensional array.
!!               For the geometries that FLASH supports, subcells volumes
!!               nrver depend on the third coordinate.
!!
!!   xl :        X-coordinate of left cell face.
!!
!!   xr :        X-coordinate of right cell face.
!!
!!   yl :        Y-coordinate of lower cell face.
!!
!!   yr :        Y-coordinate of upper cell face.
!!
!!   pos :       currently unused.
!!
!!   blockid : ID of block in current processor
!!
!!
!!***

#include "FortranLangFeatures.fh"

subroutine Grid_subcellGeometry(nsubI,nsubJ,nsubK, &
     dvCell, dvSub, xL,xR, yL,yR, pos, blockID)
  implicit none
  integer, VALUE_INTENT(IN) :: nsubI, nsubJ, nsubK
  real,    intent(in)  :: dvCell
  real,    intent(OUT) :: dvSub(nsubI, nsubJ)
  real,OPTIONAL,intent(in) :: xL, xR
  real,OPTIONAL,intent(in) :: yL, yR
  integer,OPTIONAL, intent(in) :: blockID
  integer,OPTIONAL, intent(in) :: pos(*)



end subroutine Grid_subcellGeometry
