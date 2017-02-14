!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_blockData1DRec
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
                                       ed_enforcePositiveNele, &
                                       ed_enforcePositiveTele, &
                                       ed_gradOrder

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: iminData , imaxData
  integer, intent (in) :: iminDerv , imaxDerv
  real,    intent (in) :: deltaInvI
  real,    intent (in) :: blockData (:,:)

  logical :: inBlock

  integer :: i

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: DelLeftI,  DelRightI

  integer :: cellIndex (1:3)
!
!
!     ...Compute all needed cell quantities for the # of electrons.
!
!             In block cells: Density, Volume, Zbar and Nele.
!             In guard cells: Nele.
!
!
  do i = iminData, imaxData

     cellDensity = blockData (DENS_VAR,i)
     call Eos_getAbarZbar (blockData (:,i),    abar, zbar)

     ed_cellNele (i,1,1) = cellDensity * ed_Avogadro * zbar / abar

     inBlock = (i >= iminBlock) .and. (i <= imaxBlock)

     if (inBlock) then

         cellIndex (1) = i
         cellIndex (2) = 1
         cellIndex (3) = 1

         call Grid_getSingleCellVol (blockID, EXTERIOR, cellIndex, cellVolume)

         ed_cellDensity (i,1,1) = cellDensity
         ed_cellVolume  (i,1,1) = cellVolume
         ed_cellZbar    (i,1,1) = zbar

     end if

  enddo
!
!
!     ...Compute all needed cell quantities for the electron temperature.
!
!             In guard cells: Tele.
!
!
  do i = iminData, imaxData
     ed_cellTele (i,1,1) = blockData (TELE_VAR,i)
  enddo
!
!
!     ...Next calculate the gradients. We must set the gradients for at least one extra layer
!        of guard cells (the rays are entering the block through a specific guard cell next
!        to the interior cell).
!
!
  if (ed_gradOrder == 1) then

      ed_cellGradNele (:,iminDerv:imaxDerv,1,1) = 0.0
      ed_cellGradTele (:,iminDerv:imaxDerv,1,1) = 0.0

  else if (ed_gradOrder == 2) then
!
!
!     ...2nd order x-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleX) then
          do i = iminDerv,imaxDerv
             DelLeftI  = ed_cellNele (i  ,1,1) - ed_cellNele (i-1,1,1)
             DelRightI = ed_cellNele (i+1,1,1) - ed_cellNele (i  ,1,1)
             ed_cellGradNele (IAXIS,i,1,1) = ed_mc (DelLeftI,DelRightI) * deltaInvI   ! ed_mc = slope limiter
          enddo
      else
          ed_cellGradNele (IAXIS,iminDerv:imaxDerv,1,1) = 0.0
      end if
!
!
!     ...2nd order electron temperature gradient (all components).
!
!
      do i = iminDerv,imaxDerv
         DelLeftI  = ed_cellTele (i  ,1,1) - ed_cellTele (i-1,1,1)
         DelRightI = ed_cellTele (i+1,1,1) - ed_cellTele (i  ,1,1)
         ed_cellGradTele (IAXIS,i,1,1) = ed_mc (DelLeftI,DelRightI) * deltaInvI
      enddo
!
!
!     ...Higher order gradients do not exist at the moment.
!
!
  else
      call Driver_abortFlash ("ed_blockData1DRec: No code for gradient order > 2")
  endif
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData1DRec
