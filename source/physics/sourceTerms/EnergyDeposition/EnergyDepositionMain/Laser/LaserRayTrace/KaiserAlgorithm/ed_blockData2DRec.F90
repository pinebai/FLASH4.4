!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_blockData2DRec
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

  use Driver_interface,         ONLY : Driver_abortFlash

  use Eos_interface,            ONLY : Eos_getAbarZbar

  use Grid_interface,           ONLY : Grid_getSingleCellVol

  use ed_slopeLimiters,         ONLY : ed_mc

  use EnergyDeposition_data,    ONLY : ed_Avogadro,            &
                                       ed_cellDensity,         &
                                       ed_cellGradNele,        &
                                       ed_cellGradTele,        &
                                       ed_cellNele,            &
                                       ed_cellTele,            &
                                       ed_cellVolume,          &
                                       ed_cellZbar,            &
                                       ed_cellWallThickness,   &
                                       ed_computeGradNeleX,    &
                                       ed_computeGradNeleY,    &
                                       ed_enforcePositiveNele, &
                                       ed_enforcePositiveTele, &
                                       ed_gradOrder

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: jminData , jmaxData
  integer, intent (in) :: iminDerv , imaxDerv
  integer, intent (in) :: jminDerv , jmaxDerv
  real,    intent (in) :: deltaI,    deltaJ
  real,    intent (in) :: deltaInvI, deltaInvJ
  real,    intent (in) :: blockData (:,:,:)

  logical :: inBlock

  integer :: i,j

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: DelLeftI,  DelLeftJ
  real    :: DelRightI, DelRightJ
  real    :: gradNeleX, gradNeleY
  real    :: gradTeleX, gradTeleY
  real    :: halfDeltaI, halfDeltaJ
  real    :: Nele, Tele
  real    :: maxDiff
  real    :: scalingFactor

  integer :: cellIndex (1:3)
!
!
!     ...Compute all needed cell quantities for the # of electrons.
!
!             In block cells: Density, Volume, Zbar and Nele.
!             In guard cells: Nele.
!
!
  do j = jminData,jmaxData
     do i = iminData,imaxData

        cellDensity = blockData (DENS_VAR,i,j)
        call Eos_getAbarZbar (blockData (:,i,j),    abar, zbar)

        ed_cellNele (i,j,1) = cellDensity * ed_Avogadro * zbar / abar

        inBlock =       (i >= iminBlock) &
                  .and. (i <= imaxBlock) &
                  .and. (j >= jminBlock) &
                  .and. (j <= jmaxBlock)

        if (inBlock) then

            cellIndex (1) = i
            cellIndex (2) = j
            cellIndex (3) = 1

            call Grid_getSingleCellVol (blockID, EXTERIOR, cellIndex, cellVolume)

            ed_cellDensity (i,j,1) = cellDensity
            ed_cellVolume  (i,j,1) = cellVolume
            ed_cellZbar    (i,j,1) = zbar

        end if

     enddo
  enddo
!
!
!     ...Compute all needed cell quantities for the electron temperature.
!
!             In guard cells: Tele.
!
!
  do j = jminData,jmaxData
     do i = iminData,imaxData
        ed_cellTele (i,j,1) = blockData (TELE_VAR,i,j)
     enddo
  enddo
!
!
!     ...Next calculate the gradients. We must set the gradients for at least one extra layer
!        of guard cells (the rays are entering the block through a specific guard cell next
!        to the interior cell).
!
!
  if (ed_gradOrder == 1) then

      ed_cellGradNele (:,iminDerv:imaxDerv,jminDerv:jmaxDerv,1) = 0.0
      ed_cellGradTele (:,iminDerv:imaxDerv,jminDerv:jmaxDerv,1) = 0.0

  else if (ed_gradOrder == 2) then
!
!
!     ...2nd order x-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleX) then
          do j = jminDerv,jmaxDerv
             do i = iminDerv,imaxDerv
                DelLeftI  = ed_cellNele (i  ,j,1) - ed_cellNele (i-1,j,1)
                DelRightI = ed_cellNele (i+1,j,1) - ed_cellNele (i  ,j,1)
                ed_cellGradNele (IAXIS,i,j,1) = ed_mc (DelLeftI,DelRightI) * deltaInvI   ! ed_mc = slope limiter
             enddo
          enddo
      else
          ed_cellGradNele (IAXIS,iminDerv:imaxDerv,jminDerv:jmaxDerv,1) = 0.0
      end if
!
!
!     ...2nd order y-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleY) then
          do j = jminDerv,jmaxDerv
             do i = iminDerv,imaxDerv
                DelLeftJ  = ed_cellNele (i,j  ,1) - ed_cellNele (i,j-1,1)
                DelRightJ = ed_cellNele (i,j+1,1) - ed_cellNele (i,j  ,1)
                ed_cellGradNele (JAXIS,i,j,1) = ed_mc (DelLeftJ,DelRightJ) * deltaInvJ
             enddo
          enddo
      else
          ed_cellGradNele (JAXIS,iminDerv:imaxDerv,jminDerv:jmaxDerv,1) = 0.0
      end if
!
!
!     ...2nd order electron temperature gradient (all components).
!
!
      do j = jminDerv,jmaxDerv
         do i = iminDerv,imaxDerv
            DelLeftI  = ed_cellTele (i  ,j  ,1) - ed_cellTele (i-1,j  ,1)
            DelLeftJ  = ed_cellTele (i  ,j  ,1) - ed_cellTele (i  ,j-1,1)
            DelRightI = ed_cellTele (i+1,j  ,1) - ed_cellTele (i  ,j  ,1)
            DelRightJ = ed_cellTele (i  ,j+1,1) - ed_cellTele (i  ,j  ,1)
            ed_cellGradTele (IAXIS,i,j,1) = ed_mc (DelLeftI,DelRightI) * deltaInvI
            ed_cellGradTele (JAXIS,i,j,1) = ed_mc (DelLeftJ,DelRightJ) * deltaInvJ
         enddo
      enddo
!
!
!     ...If positive Nele needs to be enforced, do so now.
!
!
      if (ed_enforcePositiveNele) then

          halfDeltaI = 0.5 * (deltaI + ed_cellWallThickness)      ! to exclude tiny -ve Nele
          halfDeltaJ = 0.5 * (deltaJ + ed_cellWallThickness)      ! numbers due to rounding errors

          do j = jminDerv,jmaxDerv
             do i = iminDerv,imaxDerv
                Nele      = ed_cellNele           (i,j,1)
                gradNeleX = ed_cellGradNele (IAXIS,i,j,1)
                gradNeleY = ed_cellGradNele (JAXIS,i,j,1)
                maxDiff   = abs (gradNeleX) * halfDeltaI + abs (gradNeleY) * halfDeltaJ

                if (maxDiff > Nele) then
                    scalingFactor = Nele / maxDiff
                    ed_cellGradNele (IAXIS,i,j,1) = scalingFactor * ed_cellGradNele (IAXIS,i,j,1)
                    ed_cellGradNele (JAXIS,i,j,1) = scalingFactor * ed_cellGradNele (JAXIS,i,j,1)
                end if

             enddo
          enddo
      end if
!
!
!     ...Likewise, if positive Tele needs to be enforced.
!
!
      if (ed_enforcePositiveTele) then

          halfDeltaI = 0.5 * (deltaI + ed_cellWallThickness)      ! to exclude tiny -ve Tele
          halfDeltaJ = 0.5 * (deltaJ + ed_cellWallThickness)      ! numbers due to rounding errors

          do j = jminDerv,jmaxDerv
             do i = iminDerv,imaxDerv
                Tele      = ed_cellTele           (i,j,1)
                gradTeleX = ed_cellGradTele (IAXIS,i,j,1)
                gradTeleY = ed_cellGradTele (JAXIS,i,j,1)
                maxDiff   = abs (gradTeleX) * halfDeltaI + abs (gradTeleY) * halfDeltaJ

                if (maxDiff > Tele) then
                    scalingFactor = Tele / maxDiff
                    ed_cellGradTele (IAXIS,i,j,1) = scalingFactor * ed_cellGradTele (IAXIS,i,j,1)
                    ed_cellGradTele (JAXIS,i,j,1) = scalingFactor * ed_cellGradTele (JAXIS,i,j,1)
                end if

             enddo
          enddo
      end if

  else
      call Driver_abortFlash ("ed_blockData2DRec: No code for gradient order > 2")
  endif
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData2DRec
