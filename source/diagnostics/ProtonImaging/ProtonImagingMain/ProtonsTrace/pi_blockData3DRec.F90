!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonsTrace/pi_blockData3DRec
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

  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_cellBfield,              &
                                  pi_cellBoundary,            &
                                  pi_cellCurlBfield,          &
                                  pi_cellEfield,              &
                                  pi_ignoreElectricalField,   &
                                  pi_screenProtonDiagnostics, &
                                  pi_speedOfLightInv,         &
                                  pi_squareRoot4Pi

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: iminBlock, imaxBlock
  integer, intent (in) :: jminBlock, jmaxBlock
  integer, intent (in) :: kminBlock, kmaxBlock
  real,    intent (in) :: deltaInvX, deltaInvY, deltaInvZ
  real,    intent (in) :: blockData (:,:,:,:)

  integer :: i,j,k

  real    :: factorE, factorB, factorX, factorY, factorZ
!
!
!     ...Get the boundary indicator.
!
!
  do k = kminBlock,kmaxBlock
  do j = jminBlock,jmaxBlock
  do i = iminBlock,imaxBlock
     pi_cellBoundary (i,j,k) = blockData (BDRY_VAR,i,j,k)
  enddo
  enddo
  enddo
!
!
!     ...Get the electric field. The electric field components need to be multiplied
!        by sqrt(4pi)/c to convert to Gauss units.
!
!
  if (pi_ignoreElectricalField) then

      do k = kminBlock,kmaxBlock
      do j = jminBlock,jmaxBlock
      do i = iminBlock,imaxBlock
         pi_cellEfield (1,i,j,k) = 0.0
         pi_cellEfield (2,i,j,k) = 0.0
         pi_cellEfield (3,i,j,k) = 0.0
      enddo
      enddo
      enddo

  else

      factorE = pi_squareRoot4Pi * pi_speedOfLightInv    ! electric field -> in Gauss

      do k = kminBlock,kmaxBlock
      do j = jminBlock,jmaxBlock
      do i = iminBlock,imaxBlock
         pi_cellEfield (1,i,j,k) = factorE * blockData (ELEX_VAR,i,j,k)
         pi_cellEfield (2,i,j,k) = factorE * blockData (ELEY_VAR,i,j,k)
         pi_cellEfield (3,i,j,k) = factorE * blockData (ELEZ_VAR,i,j,k)
      enddo
      enddo
      enddo

  end if
!
!
!     ...Get the magnetic data. Note the different units in case diagnostics are to be calculated.
!        If diagnostics are needed, we need not only B/c but also the bare B fields.
!
!
  if (pi_screenProtonDiagnostics) then

      factorB = pi_squareRoot4Pi                      ! magnetic flux density -> Gauss
      factorX = 0.5 * deltaInvX * pi_squareRoot4Pi    !
      factorY = 0.5 * deltaInvY * pi_squareRoot4Pi    ! curl of magnetic flux density -> Gauss / cm
      factorZ = 0.5 * deltaInvZ * pi_squareRoot4Pi    !

      do k = kminBlock,kmaxBlock
      do j = jminBlock,jmaxBlock
      do i = iminBlock,imaxBlock

         pi_cellBfield     (1,i,j,k) = factorB *  blockData (MAGX_VAR,i,j,k)
         pi_cellBfield     (2,i,j,k) = factorB *  blockData (MAGY_VAR,i,j,k)
         pi_cellBfield     (3,i,j,k) = factorB *  blockData (MAGZ_VAR,i,j,k)

         pi_cellCurlBfield (1,i,j,k) = factorY * (blockData (MAGZ_VAR,i,j+1,k) - blockData (MAGZ_VAR,i,j-1,k)) &
                                     - factorZ * (blockData (MAGY_VAR,i,j,k+1) - blockData (MAGY_VAR,i,j,k-1))
         pi_cellCurlBfield (2,i,j,k) = factorZ * (blockData (MAGX_VAR,i,j,k+1) - blockData (MAGX_VAR,i,j,k-1)) &
                                     - factorX * (blockData (MAGZ_VAR,i+1,j,k) - blockData (MAGZ_VAR,i-1,j,k))
         pi_cellCurlBfield (3,i,j,k) = factorX * (blockData (MAGY_VAR,i+1,j,k) - blockData (MAGY_VAR,i-1,j,k)) &
                                     - factorY * (blockData (MAGX_VAR,i,j+1,k) - blockData (MAGX_VAR,i,j-1,k))
      enddo
      enddo
      enddo

  else

      factorB = pi_squareRoot4Pi * pi_speedOfLightInv    ! magnetic flux density -> Gauss / c

      do k = kminBlock,kmaxBlock
      do j = jminBlock,jmaxBlock
      do i = iminBlock,imaxBlock
         pi_cellBfield (1,i,j,k) = factorB * blockData (MAGX_VAR,i,j,k)
         pi_cellBfield (2,i,j,k) = factorB * blockData (MAGY_VAR,i,j,k)
         pi_cellBfield (3,i,j,k) = factorB * blockData (MAGZ_VAR,i,j,k)
      enddo
      enddo
      enddo

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine pi_blockData3DRec
