!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData2DRec
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
!!  its number ID and the blockData array, which contains the needed data for the block.
!!  The following is computed and/or stored into specific arrays:
!!
!!     1) the cell Densities                                                  (stored)
!!     2) the cell Zbar values                                                (stored)
!!     3) the cell center Nele (electron number density) values               (computed)
!!     4) the cell center Tele (electron temperature) values                  (computed)
!!     5) the cell vertex Nele values                                         (computed)
!!     6) the cell vertex Tele values                                         (computed)
!!     7) the cell vertex (mixed, up to 2nd order) derivative Nele values     (computed)
!!     8) the cell vertex (mixed, up to 2nd order) derivative Tele values     (computed)
!!     9) the cell bicubic expansion coefficients for the vertex Nele grid    (stored)
!!    10) the cell bicubic expansion coefficients for the vertex Tele grid    (stored)
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
!!  iminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jmaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  iminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jmaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaI    : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaJ    : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvI : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvJ : obsolte, not needed here (but must be present for Kaiser routine compatibility)
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

  use Driver_interface,       ONLY : Driver_abortFlash

  use Eos_interface,          ONLY : Eos_getAbarZbar

  use Grid_interface,         ONLY : Grid_getSingleCellVol

  use Interpolate_interface,  ONLY : Interpolate_cubic2Dcoeffs,   &
                                     Interpolate_cubic2DmonoDerv

  use EnergyDeposition_data,  ONLY : ed_Avogadro,                   &
                                     ed_cellCubicNele,              &
                                     ed_cellCubicTele,              &
                                     ed_cellDensity,                &
                                     ed_cellVolume,                 &
                                     ed_cellZbar,                   &
                                     ed_cubicInterpolationZeroDerv

  use ed_interface,           ONLY : ed_printMatrix

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: iminData , imaxData   ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: jminData , jmaxData   ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: iminDerv , imaxDerv   ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: jminDerv , jmaxDerv   ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: deltaI   , deltaJ     ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: deltaInvI, deltaInvJ  ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: blockData (:,:,:)

  logical :: inBlock

  integer :: i, ip
  integer :: iminCenter, imaxCenter
  integer :: j, jp
  integer :: jminCenter, jmaxCenter
  integer :: numberOfCells, numberOfCellsX, numberOfCellsY
  integer :: vertexRegion

  integer :: cellIndex (1:3)

  real    :: abar, zbar
  real    :: cellDensity
  real    :: cellVolume
  real    :: centerNele, centerTele

  real, allocatable :: vertex (:,:)
  real, allocatable :: dx     (:,:)
  real, allocatable :: dy     (:,:)
  real, allocatable :: dxdy   (:,:)
!
!
!     ...Compute total number of interior cells in block as well as the number of cells in
!        each direction.
!
!
  numberOfCellsX = imaxBlock - iminBlock + 1
  numberOfCellsY = jmaxBlock - jminBlock + 1

  numberOfCells  = numberOfCellsX * numberOfCellsY
!
!
!     ...Allocate the necessary intermediate arrays.
!
!        For the cubic interpolation strategy we need to construct for each cell
!        the surrounding vertex information.  For each cell (i,j) we compute its
!        4 vertex values using bilinear interpolation of the 4 surrounding center
!        values:
!
!               vertex value = 1/4 * sum of all surrounding center values
!
!        Note, that if all center values are positive, the vertex values will also
!        be positive. The storage of the resulting vertex grid is such, that redundant
!        values are avoided. This is achieved by creating a 'vertex' array, which
!        in its (i,j) position will contain the lower left corner vertex value of the
!        cell (i,j). Thus for each cell (i,j) we have the complete set of 4 vertex
!        values stored as follows:
!
!                0 0  (lower left vertex of cell)      vertex (i,  j, )
!                1 0  (vertex on i-axis of cell)       vertex (i+1,j, )
!                0 1  (vertex on j-axis of cell)       vertex (i,  j+1)
!                1 1  (vertex on ij-plane of cell)     vertex (i+1,j+1)
!
!        In case we need to generate derivatives for smooth monotonic surfaces,
!        we need an outer layer of 3 extra cell center values. This is necessary,
!        because each vertex calculation needs 1 outer layer, each 1st order derivative
!        needs an extra outer layer of vertices and each 2nd order mixed derivative
!        needs an extra outer layer of 1st order derivatives. The following picture
!        shows the scenario:
!
!                                                 -------
!                                                |       |
!                                                |   C   |
!                                                |       |
!                                       d1 ----- d2 -----
!                                        |       |
!                                        |       |
!                                        |       |
!                                V ---- d1 ----- d1
!                                |       |
!                                |       |
!                                |       |
!                         ------ V ----- V
!                        |       |
!                        |   C   |
!                        |       |
!                         -------
!
!                           i-3     i-2     i-1      i    <-- cell index
!
!
!               center index ranges:    iminBlock - 3  (iminCenter)
!                                       imaxBlock + 3  (imaxCenter)
!                                       jminBlock - 3  (jminCenter)
!                                       jmaxBlock + 3  (jmaxCenter)
!
!               vertex index ranges:    iminBlock - 2
!                                       imaxBlock + 3
!                                       jminBlock - 2
!                                       jmaxBlock + 3
!
!               dx, dy index ranges:    iminBlock - 1
!                                       imaxBlock + 2
!                                       jminBlock - 1
!                                       jmaxBlock + 2
!
!                 dxdy index ranges:    iminBlock
!                                       imaxBlock + 1
!                                       jminBlock
!                                       jmaxBlock + 1
!
!
!        If only zero derivatives for bumpy monotonic surfaces is all that we need, then
!        we need only 1 extra layer of cell center values.
!
!
!
  if (ed_cubicInterpolationZeroDerv) then

      iminCenter = iminBlock - 1
      imaxCenter = imaxBlock + 1
      jminCenter = jminBlock - 1
      jmaxCenter = jmaxBlock + 1

      allocate (vertex (iminBlock : imaxBlock + 1 , jminBlock : jmaxBlock + 1))
      allocate (dx     (iminBlock : imaxBlock + 1 , jminBlock : jmaxBlock + 1))
      allocate (dy     (iminBlock : imaxBlock + 1 , jminBlock : jmaxBlock + 1))
      allocate (dxdy   (iminBlock : imaxBlock + 1 , jminBlock : jmaxBlock + 1))

  else

      iminCenter = iminBlock - 3
      imaxCenter = imaxBlock + 3
      jminCenter = jminBlock - 3
      jmaxCenter = jmaxBlock + 3

      allocate (vertex (iminBlock - 2 : imaxBlock + 3 , jminBlock - 2 : jmaxBlock + 3))
      allocate (dx     (iminBlock - 1 : imaxBlock + 2 , jminBlock - 1 : jmaxBlock + 2))
      allocate (dy     (iminBlock - 1 : imaxBlock + 2 , jminBlock - 1 : jmaxBlock + 2))
      allocate (dxdy   (iminBlock     : imaxBlock + 1 , jminBlock     : jmaxBlock + 1))

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
!        below by careful selection of the appropriate vertex i- and j-index
!        region. The following picture gives the location of each region:
!
!
!                ----------------------------
!               |     |                |     |
!               |  2  |        3       |  4  |
!               |     |                |     |
!                -----xxxxxxxxxxxxxxxxxx-----       x = vertices addressed
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxx   xxxxxxx     |
!               | -1  xxxxxxxx 0 xxxxxxx  1  |
!               |     xxxxxxxx   xxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     |                |     |
!               | -4  |       -3       | -2  |
!               |     |                |     |
!                ----------------------------
!
!
  vertex (:,:) = 0.0

  do j = jminCenter, jmaxCenter
     jp = j + 1
     do i = iminCenter, imaxCenter
        ip = i + 1

        cellDensity = blockData (DENS_VAR,i,j)
        call Eos_getAbarZbar (blockData (:,i,j),    abar, zbar)

        centerNele = cellDensity * ed_Avogadro * zbar / abar

        vertexRegion = 0

        if (i == iminCenter) then
            vertexRegion = vertexRegion - 1
        else if (i == imaxCenter) then
            vertexRegion = vertexRegion + 1
        end if

        if (j == jminCenter) then
            vertexRegion = vertexRegion - 3
        else if (j == jmaxCenter) then
            vertexRegion = vertexRegion + 3
        end if

        select case (vertexRegion)

        case (-4)
          vertex (ip, jp) = vertex (ip, jp) + centerNele   ! this is lower left (i,j) corner
        case (-3)
          vertex (i , jp) = vertex (i , jp) + centerNele   ! this is lower j-border
          vertex (ip, jp) = vertex (ip, jp) + centerNele
        case (-2)
          vertex (i , jp) = vertex (i , jp) + centerNele   ! this is lower right (i,j) corner
        case (-1)
          vertex (ip, j ) = vertex (ip, j ) + centerNele   ! this is left i-border
          vertex (ip, jp) = vertex (ip, jp) + centerNele
        case (0)

          vertex (i , j ) = vertex (i , j ) + centerNele   ! this is main body
          vertex (ip, j ) = vertex (ip, j ) + centerNele   ! -> full range of vertices
          vertex (i , jp) = vertex (i , jp) + centerNele
          vertex (ip, jp) = vertex (ip, jp) + centerNele

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

        case (1)
          vertex (i , j ) = vertex (i , j ) + centerNele   ! this is right i-border
          vertex (i , jp) = vertex (i , jp) + centerNele
        case (2)
          vertex (ip, j ) = vertex (ip, j ) + centerNele   ! this is upper left (i,j) corner
        case (3)
          vertex (i , j ) = vertex (i , j ) + centerNele   ! this is upper j-border
          vertex (ip, j ) = vertex (ip, j ) + centerNele
        case (4)
          vertex (i , j ) = vertex (i , j ) + centerNele   ! this is upper right (i,j) corner
        case default
          call Driver_abortFlash ('[ed_blockData2DRec] ERROR: bad vertex region selected!')
 
        end select

     enddo
  enddo

  vertex (:,:) = 0.25 * vertex (:,:)     ! handle overcounting of vertices contributions
!
!
!     ...If needed, compute the derivatives from the vertices, such that a non-trivial and
!        more smooth monotone cubic interpolation surface will result. Otherwise, set the
!        derivatives to zero, which results in a more bumpy monotone cubic interpolation
!        surface.
!
!
  if (ed_cubicInterpolationZeroDerv) then
      dx   (:,:) = 0.0
      dy   (:,:) = 0.0
      dxdy (:,:) = 0.0
  else
      call Interpolate_cubic2DmonoDerv (numberOfCellsX + 1, &
                                        numberOfCellsY + 1, &
                                        vertex,             &
                                                       dx,  &
                                                       dy,  &
                                                       dxdy )
  end if
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        bicubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the bicubic expansion
!        coefficients:
!
!
!                     1)    0 0  (lower left vertex of cell)
!                     2)    1 0  (vertex on i-axis of cell)
!                     3)    0 1  (vertex on j-axis of cell)
!                     4)    1 1  (vertex on ij-plane of cell)
!
!
!                                   j
! 
!                                   |
!                                   |
!                                   |
!                                  3 ----------- 4
!                                   |           |
!                                   |           |
!                                   |           |
!                                   |           |
!                                   |           |
!                                    ----------- ---- i
!                                  1             2
!
!
  do j = jminBlock, jmaxBlock
     jp = j + 1
     do i = iminBlock, imaxBlock
        ip = i + 1

        ed_cellCubicNele ( 1,i,j,1) = vertex (i , j )
        ed_cellCubicNele ( 2,i,j,1) = vertex (ip, j )
        ed_cellCubicNele ( 3,i,j,1) = vertex (i , jp)
        ed_cellCubicNele ( 4,i,j,1) = vertex (ip, jp)

        ed_cellCubicNele ( 5,i,j,1) =     dx (i , j )
        ed_cellCubicNele ( 6,i,j,1) =     dx (ip, j )
        ed_cellCubicNele ( 7,i,j,1) =     dx (i , jp)
        ed_cellCubicNele ( 8,i,j,1) =     dx (ip, jp)

        ed_cellCubicNele ( 9,i,j,1) =     dy (i , j )
        ed_cellCubicNele (10,i,j,1) =     dy (ip, j )
        ed_cellCubicNele (11,i,j,1) =     dy (i , jp)
        ed_cellCubicNele (12,i,j,1) =     dy (ip, jp)

        ed_cellCubicNele (13,i,j,1) =   dxdy (i , j )
        ed_cellCubicNele (14,i,j,1) =   dxdy (ip, j )
        ed_cellCubicNele (15,i,j,1) =   dxdy (i , jp)
        ed_cellCubicNele (16,i,j,1) =   dxdy (ip, jp)

     enddo
  enddo
!
!
!     ...Calculate all Nele bicubic expansion coefficients in one shot.
!
!
  call Interpolate_cubic2Dcoeffs (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values. Reuse the zero derivatives set
!        previously, in case only zero derivatives are wanted.
!
!
  vertex (:,:) = 0.0

  do j = jminCenter, jmaxCenter
     jp = j + 1
     do i = iminCenter, imaxCenter
        ip = i + 1

        centerTele = blockData (TELE_VAR,i,j)

        vertexRegion = 0

        if (i == iminCenter) then
            vertexRegion = vertexRegion - 1
        else if (i == imaxCenter) then
            vertexRegion = vertexRegion + 1
        end if

        if (j == jminCenter) then
            vertexRegion = vertexRegion - 3
        else if (j == jmaxCenter) then
            vertexRegion = vertexRegion + 3
        end if

        select case (vertexRegion)

        case (-4)
          vertex (ip, jp) = vertex (ip, jp) + centerTele
        case (-3)
          vertex (i , jp) = vertex (i , jp) + centerTele
          vertex (ip, jp) = vertex (ip, jp) + centerTele
        case (-2)
          vertex (i , jp) = vertex (i , jp) + centerTele
        case (-1)
          vertex (ip, j ) = vertex (ip, j ) + centerTele
          vertex (ip, jp) = vertex (ip, jp) + centerTele
        case (0)
          vertex (i , j ) = vertex (i , j ) + centerTele
          vertex (ip, j ) = vertex (ip, j ) + centerTele
          vertex (i , jp) = vertex (i , jp) + centerTele
          vertex (ip, jp) = vertex (ip, jp) + centerTele
        case (1)
          vertex (i , j ) = vertex (i , j ) + centerTele
          vertex (i , jp) = vertex (i , jp) + centerTele
        case (2)
          vertex (ip, j ) = vertex (ip, j ) + centerTele
        case (3)
          vertex (i , j ) = vertex (i , j ) + centerTele
          vertex (ip, j ) = vertex (ip, j ) + centerTele
        case (4)
          vertex (i , j ) = vertex (i , j ) + centerTele
        case default
          call Driver_abortFlash ('[ed_blockData2DRec] ERROR: bad vertex region selected!')
 
        end select

     enddo
  enddo

  vertex (:,:) = 0.25 * vertex (:,:)

  if (.not. ed_cubicInterpolationZeroDerv) then

       call Interpolate_cubic2DmonoDerv (numberOfCellsX + 1, &
                                         numberOfCellsY + 1, &
                                         vertex,             &
                                                        dx,  &
                                                        dy,  &
                                                        dxdy )
  end if

  do j = jminBlock, jmaxBlock
     jp = j + 1
     do i = iminBlock, imaxBlock
        ip = i + 1

        ed_cellCubicTele ( 1,i,j,1) = vertex (i , j )
        ed_cellCubicTele ( 2,i,j,1) = vertex (ip, j )
        ed_cellCubicTele ( 3,i,j,1) = vertex (i , jp)
        ed_cellCubicTele ( 4,i,j,1) = vertex (ip, jp)

        ed_cellCubicTele ( 5,i,j,1) =     dx (i , j )
        ed_cellCubicTele ( 6,i,j,1) =     dx (ip, j )
        ed_cellCubicTele ( 7,i,j,1) =     dx (i , jp)
        ed_cellCubicTele ( 8,i,j,1) =     dx (ip, jp)

        ed_cellCubicTele ( 9,i,j,1) =     dy (i , j )
        ed_cellCubicTele (10,i,j,1) =     dy (ip, j )
        ed_cellCubicTele (11,i,j,1) =     dy (i , jp)
        ed_cellCubicTele (12,i,j,1) =     dy (ip, jp)

        ed_cellCubicTele (13,i,j,1) =   dxdy (i , j )
        ed_cellCubicTele (14,i,j,1) =   dxdy (ip, j )
        ed_cellCubicTele (15,i,j,1) =   dxdy (i , jp)
        ed_cellCubicTele (16,i,j,1) =   dxdy (ip, jp)

     enddo
  enddo

  call Interpolate_cubic2Dcoeffs (numberOfCells, ed_cellCubicTele)
!
!
!     ...Deallocate the intermediate arrays.
!
!
  deallocate (vertex)
  deallocate (dx)
  deallocate (dy)
  deallocate (dxdy)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData2DRec






!  if (blockID == 13) then
!  if (1 == 0) then
!      call ed_printMatrix (6,                                          &
!                           ' Centers ',                                &
!                           iminBlock, imaxBlock, jminBlock, jmaxBlock, &
!                           iminBlock, imaxBlock, jminBlock, jmaxBlock, &
!                           centers                                     )
!
!
!      call ed_printMatrix (6,                                          &
!                           ' Vertices ',                               &
!                           iminData, imaxData+1, jminData, jmaxData+1, &
!                           iminData+1, imaxData, jminData+1, jmaxData, &
!                           vertices                                    )
!
!      call ed_printMatrix (6,                                          &
!                           ' d/dx ',                                   &
!                           iminDerv, imaxDerv, jminDerv, jmaxDerv,     &
!                           iminDerv, imaxDerv, jminDerv, jmaxDerv,     &
!  end if
