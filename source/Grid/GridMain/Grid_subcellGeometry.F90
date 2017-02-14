!!****if* source/Grid/GridMain/Grid_subcellGeometry
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
#include "constants.h"
#include "Flash.h"

subroutine Grid_subcellGeometry(nsubI,nsubJ,nsubK, &
     dvCell, dvSub, xL,xR, yL,yR, pos, blockID)
  use Grid_data, ONLY : gr_geometry
  implicit none
  integer, VALUE_INTENT(IN) :: nsubI, nsubJ, nsubK
  real,    intent(in)  :: dvCell
  real,    intent(OUT) :: dvSub(nsubI, nsubJ)
  real,OPTIONAL,intent(in) :: xL, xR
  real,OPTIONAL,intent(in) :: yL, yR
  integer,OPTIONAL, intent(in) :: blockID
  integer,OPTIONAL, intent(in) :: pos(*)

  real                 :: dvs
  real                 :: xccN2inv, xcsN
  integer              :: i, effGeometry

  effGeometry = gr_geometry
  if (nsubI == 1) effGeometry = CARTESIAN

  select case (effGeometry)
     case(CARTESIAN)
        dvs          = dvCell
        dvs          = dvs   / real(nsubI)
        dvs          = dvs   / real(nsubJ)
        dvSub(:,:)   = dvs   / real(nsubK)
     case(CYLINDRICAL)
!        xcc          =           (xL+xR) * 0.5
!        xccN         =           (xL+xR) * 0.5 * real(nsubI)
!        xccNinv      = 2.0 / (real(nsubI) * (xL+xR))
        xccN2inv     = 2.0 / (real(nsubI)*real(nsubI) * (xL+xR))
!        dxiInv = real(nsubI) / (xR - xL)
!        dxsInv = real(nsubI) / (xR - xL)
        do i=1,nsubI
!           xLs =    ( xL * (nsubI+1-i)  +  xR * (i-1) ) * dxsInv
!           xRs =    ( xL * (nsubI  -i)  +  xR * (i  ) ) * dxsInv
!!$           xcs =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) ) / real(nsubI)
!!$           dvSub(i,1) =dvCell *  xcs / xcc / real(nsubI)
!           xcs =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) ) / real(nsubI)
           xcsN =    ( xL * (real(nsubI-i)+0.5)  +  xR * (real(i)-0.5) )
           dvSub(i,1) =dvCell *  xcsN * xccN2inv
           dvSub(i,1) = dvSub(i,1) / real(nsubJ)
           dvSub(i,:) = dvSub(i,1) / real(nsubK)
        end do
     case default
        call Driver_abortFlash('[Grid_subcellGeometry] Not implemented for current geometry')
  end select


end subroutine Grid_subcellGeometry
