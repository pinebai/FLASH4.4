!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData1DRec
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
!!  its number ID and the blockData array, which contains the needed data for the block.
!!  The following is computed and/or stored into specific arrays:
!!
!!     1) the cell Densities
!!     2) the cell Volumes
!!     3) the cell Zbar values
!!     4) the cell center Nele (electron number density) values
!!     5) the cell center Tele (electron temperature) values
!!     6) the cell vertex Nele values
!!     7) the cell vertex Tele values
!!     8) the cell vertex d/dx derivative Nele values
!!     9) the cell vertex d/dx derivative Tele values
!!    10) the cell monocubic expansion coefficients for the vertex Nele grid
!!    11) the cell monocubic expansion coefficients for the vertex Tele grid
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
!!  iminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  iminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvI : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  blockData : two-dimensional array containing the block data
!!
!!***

subroutine ed_blockData1DRec (blockID,              &
                              iminBlock, imaxBlock, &
                              iminData , imaxData,  &
                              iminDerv , imaxDerv,  &
                              deltaInvI,            &
                              blockData             )

  use Driver_interface,       ONLY : Driver_abortFlash

  use Eos_interface,          ONLY : Eos_getAbarZbar

  use Grid_interface,         ONLY : Grid_getSingleCellVol

  use Interpolate_interface,  ONLY : Interpolate_cubic1Dcoeffs,   &
                                     Interpolate_cubic1DmonoDerv

  use EnergyDeposition_data,  ONLY : ed_Avogadro,                   &
                                     ed_cellDensity,                &
                                     ed_cellCubicNele,              &
                                     ed_cellCubicTele,              &
                                     ed_cellVolume,                 &
                                     ed_cellZbar,                   &
                                     ed_cellWallThickness,          &
                                     ed_cubicInterpolationZeroDerv

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: iminData , imaxData   ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: iminDerv , imaxDerv   ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: deltaInvI             ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: blockData (:,:)

  logical :: inBlock

  integer :: i, ip
  integer :: iminCenter, imaxCenter
  integer :: numberOfCells
  integer :: vertexRegion

  integer :: cellIndex (1:3)

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: centerNele, centerTele

  real, allocatable :: vertex (:)
  real, allocatable :: dx     (:)
!
!
!     ...Compute total number of interior cells in block.
!
!
  numberOfCells = imaxBlock - iminBlock + 1
!
!
!     ...Allocate the necessary intermediate arrays.
!
!        For the cubic interpolation strategy we need to construct for each cell
!        the surrounding vertex information. For each cell (i) we compute its 2
!        vertex Nele values using linear interpolation of the 2 surrounding center
!        Nele values:
!
!               vertex Nele = 1/2 * sum of two surrounding center Nele
!
!        Note, that if all center Nele values are positive, so will be
!        all the vertex Nele values. The storage of the resulting vertex Nele
!        grid is such, that redundant values are avoided. This is achieved by
!        creating a 'vertex' array, which in its (i) position will contain
!        the lower vertex Nele value of the cell (i). Thus for each cell (i)
!        we have the complete set of 2 vertex Nele values stored as follows:
!
!                  (lower vertex of cell)      vertex (i, )
!                  (upper vertex of cell)      vertex (i+1)
!
!        In case we need to generate derivatives for smooth monotonic surfaces,
!        we need extra 2 layers of cell center values. This is necessary, because
!        each vertex calculation needs 1 outer layer, each 1st order derivative
!        needs an extra outer layer of vertices:
!
!
!
!                                ------ V ----- d1 -----
!                                   C   |       |   C
!
!                                  i-2     i-1      i    <-- cell index
!
!
!               center index ranges:    iminBlock - 2  (iminCenter)
!                                       imaxBlock + 2  (imaxCenter)
!
!               vertex index ranges:    iminBlock - 1
!                                       imaxBlock + 2
!
!                   dx index ranges:    iminBlock
!                                       imaxBlock + 1
!
!
!        If only zero derivatives for bumpy monotonic surfaces is all that we need,
!        then we need only 1 extra layer of cell center values.
!
!
  if (ed_cubicInterpolationZeroDerv) then

      iminCenter = iminBlock - 1
      imaxCenter = imaxBlock + 1

      allocate (vertex (iminBlock : imaxBlock + 1))
      allocate (dx     (iminBlock : imaxBlock + 1))

  else

      iminCenter = iminBlock - 2
      imaxCenter = imaxBlock + 2

      allocate (vertex (iminBlock - 1 : imaxBlock + 2))
      allocate (dx     (iminBlock     : imaxBlock + 1))

  end if
!
!
!     ...Compute all needed center cell quantities for the # of electrons.
!
!             In block cells: Density, Volume, Zbar and Nele.
!             In guard cells: Nele.
!
!        The center Nele values are not stored, but are immediately processed
!        to compute their contribution to the vertex Nele values. Care must be
!        taken not to address vertex array elements out of bounds. This is done
!        below by careful selection of the appropriate vertex i-index region.
!        The following picture gives the location of each region:
!
!
!                -----xxxxxxxxxxxxxxxxxx-----      x = vertices addressed
!                 -1  |        0       |  1
!
!
  vertex (:) = 0.0

  do i = iminCenter, imaxCenter
     ip = i + 1

     cellDensity = blockData (DENS_VAR,i)
     call Eos_getAbarZbar (blockData (:,i),    abar, zbar)

     centerNele = cellDensity * ed_Avogadro * zbar / abar

     vertexRegion = 0

     if (i == iminCenter) then
         vertexRegion = vertexRegion - 1
     else if (i == imaxCenter) then
         vertexRegion = vertexRegion + 1
     end if

     select case (vertexRegion)

     case (-1)
       vertex (ip) = vertex (ip) + centerNele   ! this is left border
     case (0)

       vertex (i ) = vertex (i ) + centerNele   ! this is main body
       vertex (ip) = vertex (ip) + centerNele   ! -> full range of vertices

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

     case (1)
       vertex (i) = vertex (i) + centerNele   ! this is right border
     case default
       call Driver_abortFlash ('[ed_blockData1DRec] ERROR: bad vertex region selected!')
 
     end select

  enddo

  vertex (:) = 0.5 * vertex (:)     ! handle overcounting of vertices contributions
!
!
!     ...If needed, compute the derivatives from the vertices, such that a non-trivial and
!        more smooth monotone cubic interpolation surface will result. Otherwise, set the
!        derivatives to zero, which results in a more bumpy monotone cubic interpolation
!        surface.
!
!
  if (ed_cubicInterpolationZeroDerv) then
      dx (:) = 0.0
  else
      call Interpolate_cubic1DmonoDerv (numberOfCells + 1, &
                                        vertex,            &
                                                        dx )
  end if
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        monocubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the monocubic expansion
!        coefficients:
!
!                                1) lower vertex of cell
!                                2) upper vertex of cell
!
!
!                                  |------------| ---- i
!                                  1            2
!
!
  do i = iminBlock, imaxBlock
     ip = i + 1

     ed_cellCubicNele (1,i,1,1) = vertex (i )
     ed_cellCubicNele (2,i,1,1) = vertex (ip)
     ed_cellCubicNele (3,i,1,1) =     dx (i )
     ed_cellCubicNele (4,i,1,1) =     dx (ip)

  enddo
!
!
!     ...Calculate all Nele monocubic expansion coefficients in one shot.
!
!
  call Interpolate_cubic1Dcoeffs (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values. Reuse the zero derivatives set
!        previously, in case only zero derivatives are wanted.
!
!
  vertex (:) = 0.0

  do i = iminCenter, imaxCenter
     ip = i + 1

     centerTele = blockData (TELE_VAR,i)

     vertexRegion = 0

     if (i == iminCenter) then
         vertexRegion = vertexRegion - 1
     else if (i == imaxCenter) then
         vertexRegion = vertexRegion + 1
     end if

     select case (vertexRegion)

     case (-1)
       vertex (ip) = vertex (ip) + centerTele
     case (0)
       vertex (i ) = vertex (i ) + centerTele
       vertex (ip) = vertex (ip) + centerTele
     case (1)
       vertex (i ) = vertex (i ) + centerTele
     case default
       call Driver_abortFlash ('[ed_blockData1DRec] ERROR: bad vertex region selected!')
 
     end select

  enddo

  vertex (:) = 0.5 * vertex (:)

  if (.not. ed_cubicInterpolationZeroDerv) then
       call Interpolate_cubic1DmonoDerv (numberOfCells + 1, &
                                         vertex,            &
                                                         dx )
  end if

  do i = iminBlock, imaxBlock
     ip = i + 1

     ed_cellCubicTele (1,i,1,1) = vertex (i )
     ed_cellCubicTele (2,i,1,1) = vertex (ip)
     ed_cellCubicTele (3,i,1,1) =     dx (i )
     ed_cellCubicTele (4,i,1,1) =     dx (ip)

  enddo

  call Interpolate_cubic1Dcoeffs (numberOfCells, ed_cellCubicTele)
!
!
!     ...Deallocate the intermediate arrays.
!
!
  deallocate (vertex)
  deallocate (dx)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData1DRec
