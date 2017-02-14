!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/KaiserAlgorithm/ed_blockData3DRec
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

  use Driver_interface,         ONLY : Driver_abortFlash

  use Eos_interface,            ONLY : Eos_getAbarZbar

  use ed_slopeLimiters,         ONLY : ed_mc

  use EnergyDeposition_data,    ONLY : ed_Avogadro,            &
                                       ed_cellDensity,         &
                                       ed_cellGradNele,        &
                                       ed_cellGradTele,        &
                                       ed_cellNele,            &
                                       ed_cellTele,            &
                                       ed_cellZbar,            &
                                       ed_cellWallThickness,   &
                                       ed_computeGradNeleX,    &
                                       ed_computeGradNeleY,    &
                                       ed_computeGradNeleZ,    &
                                       ed_enforcePositiveNele, &
                                       ed_enforcePositiveTele, &
                                       ed_gradOrder

  implicit none

#include "constants.h"
#include "Flash.h"

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

  logical :: inBlock

  integer :: i,j,k

  real    :: abar, zbar
  real    :: cellDensity
  real    :: DelLeftI,  DelLeftJ,  DelLeftK
  real    :: DelRightI, DelRightJ, DelRightK
  real    :: gradNeleX, gradNeleY, gradNeleZ
  real    :: gradTeleX, gradTeleY, gradTeleZ
  real    :: halfDeltaI, halfDeltaJ, halfDeltaK
  real    :: Nele, Tele
  real    :: maxDiff
  real    :: scalingFactor
!
!
!     ...Compute all needed cell quantities for the # of electrons.
!
!             In block cells: Density, Zbar and Nele.
!             In guard cells: Nele.
!
!
  do k = kminData,kmaxData
     do j = jminData,jmaxData
        do i = iminData,imaxData

           cellDensity = blockData (DENS_VAR,i,j,k)
           call Eos_getAbarZbar (blockData (:,i,j,k),    abar, zbar)

           ed_cellNele (i,j,k) = cellDensity * ed_Avogadro * zbar / abar

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

        enddo
     enddo
  enddo
!
!
!     ...Compute all needed cell quantities for the electron temperature.
!
!             In guard cells: Tele.
!
!
  do k = kminData,kmaxData
     do j = jminData,jmaxData
        do i = iminData,imaxData
           ed_cellTele (i,j,k) = blockData (TELE_VAR,i,j,k)
        enddo
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

      ed_cellGradNele (:,iminDerv:imaxDerv,jminDerv:jmaxDerv,kminDerv:kmaxDerv) = 0.0
      ed_cellGradTele (:,iminDerv:imaxDerv,jminDerv:jmaxDerv,kminDerv:kmaxDerv) = 0.0

  else if (ed_gradOrder == 2) then
!
!
!     ...2nd order x-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleX) then
          do k = kminDerv,kmaxDerv
             do j = jminDerv,jmaxDerv
                do i = iminDerv,imaxDerv
                   DelLeftI  = ed_cellNele (i  ,j,k) - ed_cellNele (i-1,j,k)
                   DelRightI = ed_cellNele (i+1,j,k) - ed_cellNele (i  ,j,k)
                   ed_cellGradNele (IAXIS,i,j,k) = ed_mc (DelLeftI,DelRightI) * deltaInvI  ! ed_mc = slope limiter
                enddo
             enddo
          enddo
      else
          ed_cellGradNele (IAXIS,iminDerv:imaxDerv,jminDerv:jmaxDerv,kminDerv:kmaxDerv) = 0.0
      end if
!
!
!     ...2nd order y-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleY) then
          do k = kminDerv,kmaxDerv
             do j = jminDerv,jmaxDerv
                do i = iminDerv,imaxDerv
                   DelLeftJ  = ed_cellNele (i,j  ,k) - ed_cellNele (i,j-1,k)
                   DelRightJ = ed_cellNele (i,j+1,k) - ed_cellNele (i,j  ,k)
                   ed_cellGradNele (JAXIS,i,j,k) = ed_mc (DelLeftJ,DelRightJ) * deltaInvJ
                enddo
             enddo
          enddo
      else
          ed_cellGradNele (JAXIS,iminDerv:imaxDerv,jminDerv:jmaxDerv,kminDerv:kmaxDerv) = 0.0
      end if
!
!
!     ...2nd order z-component # of electrons gradient (if requested).
!
!
      if (ed_computeGradNeleZ) then
          do k = kminDerv,kmaxDerv
             do j = jminDerv,jmaxDerv
                do i = iminDerv,imaxDerv
                   DelLeftK  = ed_cellNele (i,j,k  ) - ed_cellNele (i,j,k-1)
                   DelRightK = ed_cellNele (i,j,k+1) - ed_cellNele (i,j,k  )
                   ed_cellGradNele (KAXIS,i,j,k) = ed_mc (DelLeftK,DelRightK) * deltaInvK
                enddo
             enddo
          enddo
      else
          ed_cellGradNele (KAXIS,iminDerv:imaxDerv,jminDerv:jmaxDerv,kminDerv:kmaxDerv) = 0.0
      end if
!
!
!     ...2nd order electron temperature gradient (all components).
!
!
      do k = kminDerv,kmaxDerv
         do j = jminDerv,jmaxDerv
            do i = iminDerv,imaxDerv

               DelLeftI  = ed_cellTele (i  ,j  ,k  ) - ed_cellTele (i-1,j  ,k  )
               DelLeftJ  = ed_cellTele (i  ,j  ,k  ) - ed_cellTele (i  ,j-1,k  )
               DelLeftK  = ed_cellTele (i  ,j  ,k  ) - ed_cellTele (i  ,j  ,k-1)
               DelRightI = ed_cellTele (i+1,j  ,k  ) - ed_cellTele (i  ,j  ,k  )
               DelRightJ = ed_cellTele (i  ,j+1,k  ) - ed_cellTele (i  ,j  ,k  )
               DelRightK = ed_cellTele (i  ,j  ,k+1) - ed_cellTele (i  ,j  ,k  )
               ed_cellGradTele (IAXIS,i,j,k) = ed_mc (DelLeftI,DelRightI) * deltaInvI
               ed_cellGradTele (JAXIS,i,j,k) = ed_mc (DelLeftJ,DelRightJ) * deltaInvJ
               ed_cellGradTele (KAXIS,i,j,k) = ed_mc (DelLeftK,DelRightK) * deltaInvK
            enddo
         enddo
      enddo
!
!
!     ...If positive Nele needs to be enforced, do so now.
!
!
      if (ed_enforcePositiveNele) then

          halfDeltaI = 0.5 * (deltaI + ed_cellWallThickness)      ! to exclude tiny -ve Nele
          halfDeltaJ = 0.5 * (deltaJ + ed_cellWallThickness)      ! numbers due to rounding
          halfDeltaK = 0.5 * (deltaK + ed_cellWallThickness)      ! errors

          do k = kminDerv,kmaxDerv
             do j = jminDerv,jmaxDerv
                do i = iminDerv,imaxDerv
                   Nele      = ed_cellNele           (i,j,k)
                   gradNeleX = ed_cellGradNele (IAXIS,i,j,k)
                   gradNeleY = ed_cellGradNele (JAXIS,i,j,k)
                   gradNeleZ = ed_cellGradNele (KAXIS,i,j,k)
                   maxDiff   =   abs (gradNeleX) * halfDeltaI &
                               + abs (gradNeleY) * halfDeltaJ &
                               + abs (gradNeleZ) * halfDeltaK

                   if (maxDiff > Nele) then
                       scalingFactor = Nele / maxDiff
                       ed_cellGradNele (IAXIS,i,j,k) = scalingFactor * ed_cellGradNele (IAXIS,i,j,k)
                       ed_cellGradNele (JAXIS,i,j,k) = scalingFactor * ed_cellGradNele (JAXIS,i,j,k)
                       ed_cellGradNele (KAXIS,i,j,k) = scalingFactor * ed_cellGradNele (KAXIS,i,j,k)
                   end if

                enddo
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
          halfDeltaJ = 0.5 * (deltaJ + ed_cellWallThickness)      ! numbers due to rounding
          halfDeltaK = 0.5 * (deltaK + ed_cellWallThickness)      ! errors

          do k = kminDerv,kmaxDerv
             do j = jminDerv,jmaxDerv
                do i = iminDerv,imaxDerv
                   Tele      = ed_cellTele           (i,j,k)
                   gradTeleX = ed_cellGradTele (IAXIS,i,j,k)
                   gradTeleY = ed_cellGradTele (JAXIS,i,j,k)
                   gradTeleZ = ed_cellGradTele (KAXIS,i,j,k)
                   maxDiff   =   abs (gradTeleX) * halfDeltaI &
                               + abs (gradTeleY) * halfDeltaJ &
                               + abs (gradTeleZ) * halfDeltaK

                   if (maxDiff > Tele) then
                       scalingFactor = Tele / maxDiff
                       ed_cellGradTele (IAXIS,i,j,k) = scalingFactor * ed_cellGradTele (IAXIS,i,j,k)
                       ed_cellGradTele (JAXIS,i,j,k) = scalingFactor * ed_cellGradTele (JAXIS,i,j,k)
                       ed_cellGradTele (KAXIS,i,j,k) = scalingFactor * ed_cellGradTele (KAXIS,i,j,k)
                   end if

                enddo
             enddo
          enddo
      end if

  else
      call Driver_abortFlash ("ed_blockData3DRec: No code for gradient order > 2")
  endif
!
!
!     ...Ready!
!
!
  return
end subroutine ed_blockData3DRec
