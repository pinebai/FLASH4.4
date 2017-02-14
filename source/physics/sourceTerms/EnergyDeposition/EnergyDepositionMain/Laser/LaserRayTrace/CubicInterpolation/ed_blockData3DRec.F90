!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/CubicInterpolation/ed_blockData3DRec
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
!!  Computes cell data for one specific block for those geometries consisting formally
!!  of 3D rectangular grids (cartesian). The block is specified through its number ID
!!  and the blockData array, which contains the needed data for the block. The following
!!  is computed and/or stored into specific arrays:
!!
!!     1) the cell Densities                                                  (stored)
!!     2) the cell Zbar values                                                (stored)
!!     3) the cell center Nele (electron number density) values               (computed)
!!     4) the cell center Tele (electron temperature) values                  (computed)
!!     5) the cell vertex Nele values                                         (computed)
!!     6) the cell vertex Tele values                                         (computed)
!!     7) the cell vertex (mixed, up to 3rd order) derivative Nele values     (computed)
!!     8) the cell vertex (mixed, up to 3rd order) derivative Tele values     (computed)
!!     9) the cell tricubic expansion coefficients for the vertex Nele grid   (stored)
!!    10) the cell tricubic expansion coefficients for the vertex Tele grid   (stored)
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
!!  iminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jmaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  kminData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  kmaxData  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  iminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  imaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  jmaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  kminDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  kmaxDerv  : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaI    : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaJ    : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaK    : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvI : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvJ : obsolte, not needed here (but must be present for Kaiser routine compatibility)
!!  deltaInvK : obsolte, not needed here (but must be present for Kaiser routine compatibility)
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

  use Driver_interface,       ONLY : Driver_abortFlash

  use Eos_interface,          ONLY : Eos_getAbarZbar

  use Interpolate_interface,  ONLY : Interpolate_cubic3Dcoeffs,   &
                                     Interpolate_cubic3DmonoDerv

  use ed_interface,           ONLY : ed_printMatrix

  use EnergyDeposition_data,  ONLY : ed_Avogadro,                   &
                                     ed_cellCubicNele,              &
                                     ed_cellCubicTele,              &
                                     ed_cellDensity,                &
                                     ed_cellZbar,                   &
                                     ed_cubicInterpolationZeroDerv

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  integer, intent (in) :: iminData , imaxData             ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: jminData , jmaxData             ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: kminData , kmaxData             ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: iminDerv , imaxDerv             ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: jminDerv , jmaxDerv             ! not used, but needed for Kaiser routine compatibility
  integer, intent (in) :: kminDerv , kmaxDerv             ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: deltaI   , deltaJ   , deltaK    ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: deltaInvI, deltaInvJ, deltaInvK ! not used, but needed for Kaiser routine compatibility
  real,    intent (in) :: blockData (:,:,:,:)

  logical :: inBlock

  integer :: i, ip
  integer :: iminCenter, imaxCenter
  integer :: j, jp
  integer :: jminCenter, jmaxCenter
  integer :: k, kp
  integer :: kminCenter, kmaxCenter
  integer :: numberOfCells, numberOfCellsX, numberOfCellsY, numberOfCellsZ
  integer :: vertexRegion

  real    :: abar, zbar
  real    :: cellDensity
  real    :: centerNele, centerTele

  real, allocatable :: vertex (:,:,:)
  real, allocatable :: dx     (:,:,:)
  real, allocatable :: dy     (:,:,:)
  real, allocatable :: dz     (:,:,:)
  real, allocatable :: dxdy   (:,:,:)
  real, allocatable :: dxdz   (:,:,:)
  real, allocatable :: dydz   (:,:,:)
  real, allocatable :: dxdydz (:,:,:)
!
!
!     ...Compute total number of interior cells in block as well as the number of cells in
!        each direction.
!
!
  numberOfCellsX = imaxBlock - iminBlock + 1
  numberOfCellsY = jmaxBlock - jminBlock + 1
  numberOfCellsZ = kmaxBlock - kminBlock + 1

  numberOfCells  = numberOfCellsX * numberOfCellsY * numberOfCellsZ
!
!
!     ...Allocate the necessary intermediate arrays.
!
!        For the cubic interpolation strategy we need to construct for each cell
!        the surrounding vertex information. For each cell (i,j,k) we compute its
!        8 vertex values using trilinear interpolation of the 8 surrounding center
!        values:
!
!               vertex value = 1/8 * sum of all surrounding center values
!
!        Note, that if all center values are positive, the vertex values will also
!        be positive. The storage of the resulting vertex  grid is such, that redundant
!        values are avoided. This is achieved by creating a 'vertex' array, which in
!        its (i,j,k) position will contain the lower left corner vertex value of the
!        cell (i,j,k). Thus for each cell (i,j,k) we have the complete set of 8 vertex
!        values stored as follows:
!
!                0 0 0  (lower left vertex of cell)      vertex (i,  j,  k  )
!                1 0 0  (vertex on i-axis of cell)       vertex (i+1,j,  k  )
!                0 1 0  (vertex on j-axis of cell)       vertex (i,  j+1,k  )
!                0 0 1  (vertex on k-axis of cell)       vertex (i,  j,  k+1)
!                1 1 0  (vertex on ij-plane of cell)     vertex (i+1,j+1,k  )
!                1 0 1  (vertex on ik-plane of cell)     vertex (i+1,j,  k+1)
!                0 1 1  (vertex on jk-plane of cell)     vertex (i,  j+1,k+1)
!                1 1 1  (upper right vertex of cell)     vertex (i+1,j+1,k+1)
!
!        In case we need to generate derivatives for smooth monotonic surfaces,
!        we need an outer layer of 4 extra cell center values. This is necessary,
!        because each vertex calculation needs 1 outer layer, each 1st order derivative
!        needs an extra outer layer of vertices, each 2nd order mixed derivative
!        needs an extra outer layer of 1st order derivatives and each 3rd order mixed
!        derivative needs an extra outer layer of 2nd order mixed derivatives. The
!        following picture shows the scenario:
!
!
!                                                         -------
!                                                        |       |
!                                                        |   C   |
!                                                        |       |
!                                                d2 ---- d3 -----
!                                                |       |
!                                                |       |
!                                                |       |
!                                        d1 ---- d2 ---- d2
!                                        |       |
!                                        |       |
!                                        |       |
!                                V ----- d1 ---- d1
!                                |       |
!                                |       |
!                                |       |
!                         ------ V ----- V
!                        |       |
!                        |   C   |
!                        |       |
!                         -------
!
!                           i-4     i-3     i-2      i-1     i    <-- cell index
!
!
!                    center index ranges:    iminBlock - 4  (iminCenter)
!                                            imaxBlock + 4  (imaxCenter)
!                                            jminBlock - 4  (jminCenter)
!                                            jmaxBlock + 4  (jmaxCenter)
!                                            kminBlock - 4  (kminCenter)
!                                            kmaxBlock + 4  (kmaxCenter)
!
!                    vertex index ranges:    iminBlock - 3
!                                            imaxBlock + 4
!                                            jminBlock - 3
!                                            jmaxBlock + 4
!                                            kminBlock - 3
!                                            kmaxBlock + 4
!
!                dx, dy, dz index ranges:    iminBlock - 2
!                                            imaxBlock + 3
!                                            jminBlock - 2
!                                            jmaxBlock + 3
!                                            kminBlock - 2
!                                            kmaxBlock + 3
!
!          dxdy, dxdz, dydz index ranges:    iminBlock - 1
!                                            imaxBlock + 2
!                                            jminBlock - 1
!                                            jmaxBlock + 2
!                                            kminBlock - 1
!                                            kmaxBlock + 2
!
!                    dxdydz index ranges:    iminBlock
!                                            imaxBlock + 1
!                                            jminBlock
!                                            jmaxBlock + 1
!                                            kminBlock
!                                            kmaxBlock + 1
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
      kminCenter = kminBlock - 1
      kmaxCenter = kmaxBlock + 1

      allocate (vertex (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dx     (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dy     (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dz     (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dxdy   (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dxdz   (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dydz   (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))
      allocate (dxdydz (iminBlock:imaxBlock+1 , jminBlock:jmaxBlock+1 , kminBlock:kmaxBlock+1))

  else

      iminCenter = iminBlock - 4
      imaxCenter = imaxBlock + 4
      jminCenter = jminBlock - 4
      jmaxCenter = jmaxBlock + 4
      kminCenter = kminBlock - 4
      kmaxCenter = kmaxBlock + 4

      allocate (vertex (iminBlock-3:imaxBlock+4 , jminBlock-3:jmaxBlock+4 , kminBlock-3:kmaxBlock+4))
      allocate (dx     (iminBlock-2:imaxBlock+3 , jminBlock-2:jmaxBlock+3 , kminBlock-2:kmaxBlock+3))
      allocate (dy     (iminBlock-2:imaxBlock+3 , jminBlock-2:jmaxBlock+3 , kminBlock-2:kmaxBlock+3))
      allocate (dz     (iminBlock-2:imaxBlock+3 , jminBlock-2:jmaxBlock+3 , kminBlock-2:kmaxBlock+3))
      allocate (dxdy   (iminBlock-1:imaxBlock+2 , jminBlock-1:jmaxBlock+2 , kminBlock-1:kmaxBlock+2))
      allocate (dxdz   (iminBlock-1:imaxBlock+2 , jminBlock-1:jmaxBlock+2 , kminBlock-1:kmaxBlock+2))
      allocate (dydz   (iminBlock-1:imaxBlock+2 , jminBlock-1:jmaxBlock+2 , kminBlock-1:kmaxBlock+2))
      allocate (dxdydz (iminBlock  :imaxBlock+1 , jminBlock  :jmaxBlock+1 , kminBlock  :kmaxBlock+1))

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
!        below by careful selection of the appropriate vertex i-,j- and k-index
!        region. The following picture gives the location of each region
!        (x = vertices addressed):
!
!
!         Regions: k = kminCenter
!
!                ----------------------------
!               |     |                |     |
!               |  -7 |       -6       |  -5 |    j = jmaxCenter
!               |     |                |     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxx   xxxxxxx     |
!               | -10 xxxxxxx -9  xxxxxx  -8 |    jminCenter < j < jmaxCenter
!               |     xxxxxxxx   xxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     |                |     |
!               | -13 |      -12       | -11 |    j = jminCenter
!               |     |                |     |
!                ----------------------------
!
!         Regions: kminCenter < k < kmaxCenter
!
!                ----------------------------
!               |     |                |     |
!               |  2  |        3       |  4  |    j = jmaxCenter
!               |     |                |     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxx   xxxxxxx     |
!               |  -1 xxxxxxxx 0 xxxxxxx  1  |    jminCenter < j < jmaxCenter
!               |     xxxxxxxx   xxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     |                |     |
!               |  -4 |       -3       |  -2 |    j = jminCenter
!               |     |                |     |
!                ----------------------------
!
!         Regions: k = kmaxCenter
!
!                ----------------------------
!               |     |                |     |
!               |  11 |       12       |  13 |    j = jmaxCenter
!               |     |                |     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxx   xxxxxxx     |
!               |  8  xxxxxxxx 9 xxxxxxx  10 |    jminCenter < j < jmaxCenter
!               |     xxxxxxxx   xxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!               |     xxxxxxxxxxxxxxxxxx     |
!                -----xxxxxxxxxxxxxxxxxx-----
!               |     |                |     |
!               |  5  |       6        |  7  |    j = jminCenter
!               |     |                |     |
!                ----------------------------
!
!                  |          |           |
!                  |          |           |
!                             |
!           i = iminCenter    |    i = imaxCenter
!                             |
!                iminCenter < i < imaxCenter
!
!
  vertex (:,:,:) = 0.0

  do k = kminCenter, kmaxCenter
     kp = k + 1
     do j = jminCenter, jmaxCenter
        jp = j + 1
        do i = iminCenter, imaxCenter
           ip = i + 1

           cellDensity = blockData (DENS_VAR,i,j,k)
           call Eos_getAbarZbar (blockData (:,i,j,k),    abar, zbar)

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

           if (k == kminCenter) then
               vertexRegion = vertexRegion - 9
           else if (k == kmaxCenter) then
               vertexRegion = vertexRegion + 9
           end if

           select case (vertexRegion)

           case (-13)
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-12)
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-11)
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
           case (-10)
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-9)
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-8)
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
           case (-7)
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
           case (-6)
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
           case (-5)
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
           case (-4)
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-3)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (-2)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
           case (-1)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele
           case (0)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerNele

             inBlock =       (i >= iminBlock) &
                       .and. (i <= imaxBlock) &
                       .and. (j >= jminBlock) &
                       .and. (j <= jmaxBlock) &
                       .and. (k >= kminBlock) &
                       .and. (k <= kmaxBlock)

             if (inBlock) then
                 ed_cellDensity (i,j,k) = cellDensity
                 ed_cellZbar    (i,j,k) = zbar
             end if

           case (1)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerNele
           case (2)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
           case (3)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerNele
           case (4)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (i , j , kp) = vertex (i , j , kp) + centerNele
           case (5)
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
           case (6)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
           case (7)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
           case (8)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
           case (9)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerNele
           case (10)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerNele
           case (11)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
           case (12)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerNele
           case (13)
             vertex (i , j , k ) = vertex (i , j , k ) + centerNele
           case default
             call Driver_abortFlash ('[ed_blockData3DRec] ERROR: bad vertex region selected!')
 
           end select

        enddo
     enddo
  enddo

  vertex (:,:,:) = 0.125 * vertex (:,:,:)     ! handle overcounting of vertices contributions
!
!
!     ...If needed, compute the derivatives from the vertices, such that a non-trivial and
!        more smooth monotone cubic interpolation surface will result. Otherwise, set the
!        derivatives to zero, which results in a more bumpy monotone cubic interpolation surface.
!
!
  if (ed_cubicInterpolationZeroDerv) then
      dx     (:,:,:) = 0.0
      dy     (:,:,:) = 0.0
      dz     (:,:,:) = 0.0
      dxdy   (:,:,:) = 0.0
      dxdz   (:,:,:) = 0.0
      dydz   (:,:,:) = 0.0
      dxdydz (:,:,:) = 0.0
  else
      call Interpolate_cubic3DmonoDerv (numberOfCellsX + 1,              &
                                        numberOfCellsY + 1,              &
                                        numberOfCellsZ + 1,              &
                                        vertex,                          &
                                                       dx,   dy,   dz,   &
                                                       dxdy, dxdz, dydz, &
                                                       dxdydz            )
  end if
!
!
!     ...Scatter the vertices/derivative data into the appropriate cell places for obtaining the
!        tricubic expansion coefficients. The order of the vertices for each cell will be such
!        that it corresponds to the order needed for evaluating the tricubic expansion
!        coefficients:
!
!
!                   1)       0 0 0  (lower left vertex of cell)
!                   2)       1 0 0  (vertex on i-axis of cell)
!                   3)       0 1 0  (vertex on j-axis of cell)
!                   4)       0 0 1  (vertex on k-axis of cell)
!                   5)       1 1 0  (vertex on ij-plane of cell)
!                   6)       1 0 1  (vertex on ik-plane of cell)
!                   7)       0 1 1  (vertex on jk-plane of cell)
!                   8)       1 1 1  (upper right vertex of cell)
!
!
!                             k
!                             
!                             |
!                             | 7  -----------  8
!                             |  /|          /|
!                               / |         / |
!                              /  |        /  |
!                           4  -----------  6 |
!                             |   | j     |   |
!                             | 3  -----------  5
!                             |  /        |  /
!                             | /         | /
!                             |/          |/
!                              -----------    ---- i
!                           1              2
!
!
  do k = kminBlock, kmaxBlock
     kp = k + 1
     do j = jminBlock, jmaxBlock
        jp = j + 1
        do i = iminBlock, imaxBlock
           ip = i + 1

           ed_cellCubicNele ( 1,i,j,k) = vertex (i , j , k )
           ed_cellCubicNele ( 2,i,j,k) = vertex (ip, j  ,k )
           ed_cellCubicNele ( 3,i,j,k) = vertex (i , jp ,k )
           ed_cellCubicNele ( 4,i,j,k) = vertex (i , j  ,kp)
           ed_cellCubicNele ( 5,i,j,k) = vertex (ip, jp ,k )
           ed_cellCubicNele ( 6,i,j,k) = vertex (ip, j  ,kp)
           ed_cellCubicNele ( 7,i,j,k) = vertex (i , jp ,kp)
           ed_cellCubicNele ( 8,i,j,k) = vertex (ip, jp ,kp)

           ed_cellCubicNele ( 9,i,j,k) =     dx (i , j , k )
           ed_cellCubicNele (10,i,j,k) =     dx (ip, j  ,k )
           ed_cellCubicNele (11,i,j,k) =     dx (i , jp ,k )
           ed_cellCubicNele (12,i,j,k) =     dx (i , j  ,kp)
           ed_cellCubicNele (13,i,j,k) =     dx (ip, jp ,k )
           ed_cellCubicNele (14,i,j,k) =     dx (ip, j  ,kp)
           ed_cellCubicNele (15,i,j,k) =     dx (i , jp ,kp)
           ed_cellCubicNele (16,i,j,k) =     dx (ip, jp ,kp)

           ed_cellCubicNele (17,i,j,k) =     dy (i , j , k )
           ed_cellCubicNele (18,i,j,k) =     dy (ip, j  ,k )
           ed_cellCubicNele (19,i,j,k) =     dy (i , jp ,k )
           ed_cellCubicNele (20,i,j,k) =     dy (i , j  ,kp)
           ed_cellCubicNele (21,i,j,k) =     dy (ip, jp ,k )
           ed_cellCubicNele (22,i,j,k) =     dy (ip, j  ,kp)
           ed_cellCubicNele (23,i,j,k) =     dy (i , jp ,kp)
           ed_cellCubicNele (24,i,j,k) =     dy (ip, jp ,kp)

           ed_cellCubicNele (25,i,j,k) =     dz (i , j , k )
           ed_cellCubicNele (26,i,j,k) =     dz (ip, j  ,k )
           ed_cellCubicNele (27,i,j,k) =     dz (i , jp ,k )
           ed_cellCubicNele (28,i,j,k) =     dz (i , j  ,kp)
           ed_cellCubicNele (29,i,j,k) =     dz (ip, jp ,k )
           ed_cellCubicNele (30,i,j,k) =     dz (ip, j  ,kp)
           ed_cellCubicNele (31,i,j,k) =     dz (i , jp ,kp)
           ed_cellCubicNele (32,i,j,k) =     dz (ip, jp ,kp)

           ed_cellCubicNele (33,i,j,k) =   dxdy (i , j , k )
           ed_cellCubicNele (34,i,j,k) =   dxdy (ip, j  ,k )
           ed_cellCubicNele (35,i,j,k) =   dxdy (i , jp ,k )
           ed_cellCubicNele (36,i,j,k) =   dxdy (i , j  ,kp)
           ed_cellCubicNele (37,i,j,k) =   dxdy (ip, jp ,k )
           ed_cellCubicNele (38,i,j,k) =   dxdy (ip, j  ,kp)
           ed_cellCubicNele (39,i,j,k) =   dxdy (i , jp ,kp)
           ed_cellCubicNele (40,i,j,k) =   dxdy (ip, jp ,kp)

           ed_cellCubicNele (41,i,j,k) =   dxdz (i , j , k )
           ed_cellCubicNele (42,i,j,k) =   dxdz (ip, j  ,k )
           ed_cellCubicNele (43,i,j,k) =   dxdz (i , jp ,k )
           ed_cellCubicNele (44,i,j,k) =   dxdz (i , j  ,kp)
           ed_cellCubicNele (45,i,j,k) =   dxdz (ip, jp ,k )
           ed_cellCubicNele (46,i,j,k) =   dxdz (ip, j  ,kp)
           ed_cellCubicNele (47,i,j,k) =   dxdz (i , jp ,kp)
           ed_cellCubicNele (48,i,j,k) =   dxdz (ip, jp ,kp)

           ed_cellCubicNele (49,i,j,k) =   dydz (i , j , k )
           ed_cellCubicNele (50,i,j,k) =   dydz (ip, j  ,k )
           ed_cellCubicNele (51,i,j,k) =   dydz (i , jp ,k )
           ed_cellCubicNele (52,i,j,k) =   dydz (i , j  ,kp)
           ed_cellCubicNele (53,i,j,k) =   dydz (ip, jp ,k )
           ed_cellCubicNele (54,i,j,k) =   dydz (ip, j  ,kp)
           ed_cellCubicNele (55,i,j,k) =   dydz (i , jp ,kp)
           ed_cellCubicNele (56,i,j,k) =   dydz (ip, jp ,kp)

           ed_cellCubicNele (57,i,j,k) = dxdydz (i , j , k )
           ed_cellCubicNele (58,i,j,k) = dxdydz (ip, j  ,k )
           ed_cellCubicNele (59,i,j,k) = dxdydz (i , jp ,k )
           ed_cellCubicNele (60,i,j,k) = dxdydz (i , j  ,kp)
           ed_cellCubicNele (61,i,j,k) = dxdydz (ip, jp ,k )
           ed_cellCubicNele (62,i,j,k) = dxdydz (ip, j  ,kp)
           ed_cellCubicNele (63,i,j,k) = dxdydz (i , jp ,kp)
           ed_cellCubicNele (64,i,j,k) = dxdydz (ip, jp ,kp)

        enddo
     enddo
  enddo
!
!
!     ...Calculate all Nele tricubic expansion coefficients in one shot.
!
!
  call Interpolate_cubic3Dcoeffs (numberOfCells, ed_cellCubicNele)
!
!
!     ...Repeat the same procedure for the Tele values. Reuse the zero derivatives set
!        previously, in case only zero derivatives are wanted.
!
!
  vertex (:,:,:) = 0.0

  do k = kminCenter, kmaxCenter
     kp = k + 1
     do j = jminCenter, jmaxCenter
        jp = j + 1
        do i = iminCenter, imaxCenter
           ip = i + 1

           centerTele = blockData (TELE_VAR,i,j,k)

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

           if (k == kminCenter) then
               vertexRegion = vertexRegion - 9
           else if (k == kmaxCenter) then
               vertexRegion = vertexRegion + 9
           end if

           select case (vertexRegion)

           case (-13)
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-12)
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-11)
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
           case (-10)
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-9)
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-8)
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
           case (-7)
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
           case (-6)
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
           case (-5)
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
           case (-4)
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-3)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (-2)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
           case (-1)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (0)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
             vertex (ip, jp, kp) = vertex (ip, jp, kp) + centerTele
           case (1)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (i , jp, kp) = vertex (i , jp, kp) + centerTele
           case (2)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
           case (3)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
             vertex (ip, j , kp) = vertex (ip, j , kp) + centerTele
           case (4)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (i , j , kp) = vertex (i , j , kp) + centerTele
           case (5)
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
           case (6)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
           case (7)
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
           case (8)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
           case (9)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
             vertex (ip, jp, k ) = vertex (ip, jp, k ) + centerTele
           case (10)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (i , jp, k ) = vertex (i , jp, k ) + centerTele
           case (11)
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
           case (12)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
             vertex (ip, j , k ) = vertex (ip, j , k ) + centerTele
           case (13)
             vertex (i , j , k ) = vertex (i , j , k ) + centerTele
           case default
             call Driver_abortFlash ('[ed_blockData3DRec] ERROR: bad vertex region selected!')
 
           end select

        enddo
     enddo
  enddo

  vertex (:,:,:) = 0.125 * vertex (:,:,:)

  if (.not. ed_cubicInterpolationZeroDerv) then

       call Interpolate_cubic3DmonoDerv (numberOfCellsX + 1,              &
                                         numberOfCellsY + 1,              &
                                         numberOfCellsZ + 1,              &
                                         vertex,                          &
                                                        dx,   dy,   dz,   &
                                                        dxdy, dxdz, dydz, &
                                                        dxdydz            )
  end if

  do k = kminBlock, kmaxBlock
     kp = k + 1
     do j = jminBlock, jmaxBlock
        jp = j + 1
        do i = iminBlock, imaxBlock
           ip = i + 1

           ed_cellCubicTele ( 1,i,j,k) = vertex (i , j , k )
           ed_cellCubicTele ( 2,i,j,k) = vertex (ip, j  ,k )
           ed_cellCubicTele ( 3,i,j,k) = vertex (i , jp ,k )
           ed_cellCubicTele ( 4,i,j,k) = vertex (i , j  ,kp)
           ed_cellCubicTele ( 5,i,j,k) = vertex (ip, jp ,k )
           ed_cellCubicTele ( 6,i,j,k) = vertex (ip, j  ,kp)
           ed_cellCubicTele ( 7,i,j,k) = vertex (i , jp ,kp)
           ed_cellCubicTele ( 8,i,j,k) = vertex (ip, jp ,kp)

           ed_cellCubicTele ( 9,i,j,k) =     dx (i , j , k )
           ed_cellCubicTele (10,i,j,k) =     dx (ip, j  ,k )
           ed_cellCubicTele (11,i,j,k) =     dx (i , jp ,k )
           ed_cellCubicTele (12,i,j,k) =     dx (i , j  ,kp)
           ed_cellCubicTele (13,i,j,k) =     dx (ip, jp ,k )
           ed_cellCubicTele (14,i,j,k) =     dx (ip, j  ,kp)
           ed_cellCubicTele (15,i,j,k) =     dx (i , jp ,kp)
           ed_cellCubicTele (16,i,j,k) =     dx (ip, jp ,kp)

           ed_cellCubicTele (17,i,j,k) =     dy (i , j , k )
           ed_cellCubicTele (18,i,j,k) =     dy (ip, j  ,k )
           ed_cellCubicTele (19,i,j,k) =     dy (i , jp ,k )
           ed_cellCubicTele (20,i,j,k) =     dy (i , j  ,kp)
           ed_cellCubicTele (21,i,j,k) =     dy (ip, jp ,k )
           ed_cellCubicTele (22,i,j,k) =     dy (ip, j  ,kp)
           ed_cellCubicTele (23,i,j,k) =     dy (i , jp ,kp)
           ed_cellCubicTele (24,i,j,k) =     dy (ip, jp ,kp)

           ed_cellCubicTele (25,i,j,k) =     dz (i , j , k )
           ed_cellCubicTele (26,i,j,k) =     dz (ip, j  ,k )
           ed_cellCubicTele (27,i,j,k) =     dz (i , jp ,k )
           ed_cellCubicTele (28,i,j,k) =     dz (i , j  ,kp)
           ed_cellCubicTele (29,i,j,k) =     dz (ip, jp ,k )
           ed_cellCubicTele (30,i,j,k) =     dz (ip, j  ,kp)
           ed_cellCubicTele (31,i,j,k) =     dz (i , jp ,kp)
           ed_cellCubicTele (32,i,j,k) =     dz (ip, jp ,kp)

           ed_cellCubicTele (33,i,j,k) =   dxdy (i , j , k )
           ed_cellCubicTele (34,i,j,k) =   dxdy (ip, j  ,k )
           ed_cellCubicTele (35,i,j,k) =   dxdy (i , jp ,k )
           ed_cellCubicTele (36,i,j,k) =   dxdy (i , j  ,kp)
           ed_cellCubicTele (37,i,j,k) =   dxdy (ip, jp ,k )
           ed_cellCubicTele (38,i,j,k) =   dxdy (ip, j  ,kp)
           ed_cellCubicTele (39,i,j,k) =   dxdy (i , jp ,kp)
           ed_cellCubicTele (40,i,j,k) =   dxdy (ip, jp ,kp)

           ed_cellCubicTele (41,i,j,k) =   dxdz (i , j , k )
           ed_cellCubicTele (42,i,j,k) =   dxdz (ip, j  ,k )
           ed_cellCubicTele (43,i,j,k) =   dxdz (i , jp ,k )
           ed_cellCubicTele (44,i,j,k) =   dxdz (i , j  ,kp)
           ed_cellCubicTele (45,i,j,k) =   dxdz (ip, jp ,k )
           ed_cellCubicTele (46,i,j,k) =   dxdz (ip, j  ,kp)
           ed_cellCubicTele (47,i,j,k) =   dxdz (i , jp ,kp)
           ed_cellCubicTele (48,i,j,k) =   dxdz (ip, jp ,kp)

           ed_cellCubicTele (49,i,j,k) =   dydz (i , j , k )
           ed_cellCubicTele (50,i,j,k) =   dydz (ip, j  ,k )
           ed_cellCubicTele (51,i,j,k) =   dydz (i , jp ,k )
           ed_cellCubicTele (52,i,j,k) =   dydz (i , j  ,kp)
           ed_cellCubicTele (53,i,j,k) =   dydz (ip, jp ,k )
           ed_cellCubicTele (54,i,j,k) =   dydz (ip, j  ,kp)
           ed_cellCubicTele (55,i,j,k) =   dydz (i , jp ,kp)
           ed_cellCubicTele (56,i,j,k) =   dydz (ip, jp ,kp)

           ed_cellCubicTele (57,i,j,k) = dxdydz (i , j , k )
           ed_cellCubicTele (58,i,j,k) = dxdydz (ip, j  ,k )
           ed_cellCubicTele (59,i,j,k) = dxdydz (i , jp ,k )
           ed_cellCubicTele (60,i,j,k) = dxdydz (i , j  ,kp)
           ed_cellCubicTele (61,i,j,k) = dxdydz (ip, jp ,k )
           ed_cellCubicTele (62,i,j,k) = dxdydz (ip, j  ,kp)
           ed_cellCubicTele (63,i,j,k) = dxdydz (i , jp ,kp)
           ed_cellCubicTele (64,i,j,k) = dxdydz (ip, jp ,kp)

        enddo
     enddo
  enddo

  call Interpolate_cubic3Dcoeffs (numberOfCells, ed_cellCubicTele)
!
!
!     ...Deallocate the intermediate arrays.
!
!
  deallocate (vertex)
  deallocate (dx)
  deallocate (dy)
  deallocate (dz)
  deallocate (dxdy)
  deallocate (dxdz)
  deallocate (dydz)
  deallocate (dxdydz)
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData3DRec





!  call ed_printMatrix (6,                      &
!                       ' Vertex for k = 5 ',   &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       vertex (:,:,kminBlock) )
!
!  call ed_printMatrix (6,                      &
!                       ' Vertex for k = 6 ',   &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       vertex (:,:,kminBlock+1) )
!
!  call ed_printMatrix (6,                      &
!                       ' Vertex for k = 7 ',   &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       vertex (:,:,kminBlock+2) )
!
!  call ed_printMatrix (6,                      &
!                       ' Vertex for k = 8 ',   &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       vertex (:,:,kminBlock+3) )
!
!  call ed_printMatrix (6,                      &
!                       ' Vertex for k = 9 ',   &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       iminBlock, imaxBlock+1, &
!                       jminBlock, jmaxBlock+1, &
!                       vertex (:,:,kminBlock+4) )
