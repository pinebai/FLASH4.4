!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_initIntegrate
!!
!! NAME
!!
!!  op_initIntegrate
!!
!! SYNOPSIS
!!
!!  call op_initIntegrate ()
!!
!! DESCRIPTION
!!
!!  Inititalizes the integration section of the opacity unit.
!!  Several arrays are allocated here according to specific needs.
!!
!! ARGUMENTS
!!
!!***
subroutine op_initIntegrate ()

  use Driver_interface,            ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use op_integrateData,            ONLY : op_initializedIntegrate, &
                                          op_useQuadrature,        &
                                          op_useRomberg,           &
                                          op_maxRoots,             &
                                          op_maxWork,              &
                                          op_AuxPolynomialA,       &
                                          op_AuxPolynomialB,       &
                                          op_Moments,              &
                                          op_JmatrixDiagonals,     &
                                          op_JmatrixOffdiagonals,  &
                                          op_OrthoPolynomialA,     &
                                          op_OrthoPolynomialB,     &
                                          op_work1,                &
                                          op_work2,                &
                                          op_work3,                &
                                          op_printQuadratureData,  &
                                          op_maxRombergSteps,      &
                                          op_RombergAccuracy,      &
                                          op_RombergRow,           &
                                          op_RombergIntegral

  implicit none

# include "Opacity.h"

  integer :: maxMoments
  integer :: status
!
!
!    ...Do everything for the quadrature integration (if needed).
!       Query the maximum number of quadrature roots.
!       Allocate the arrays of size the number of roots.
!
!
  call RuntimeParameters_get ("opacity_useQuadrature",  op_useQuadrature)

  if (op_useQuadrature) then

      call RuntimeParameters_get ("opacity_printQuadratureData",   op_printQuadratureData)
      call RuntimeParameters_get ("opacity_maxQuadratureRoots",    op_maxRoots)

      allocate (op_JmatrixDiagonals (1:op_maxRoots), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_JmatrixDiagonals allocation failed')
      end if

      allocate (op_JmatrixOffdiagonals (1:op_maxRoots), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_JmatrixOffdiagonals allocation failed')
      end if

      allocate (op_OrthoPolynomialA (1:op_maxRoots), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_OrthoPolynomialA allocation failed')
      end if

      allocate (op_OrthoPolynomialB (1:op_maxRoots), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_OrthoPolynomialB allocation failed')
      end if
!
!
!    ...Determine the maximum number of moments needed.
!       Allocate the arrays of size the number of moments.
!
!
      maxMoments = 2 * op_maxRoots - 1

      allocate (op_AuxPolynomialA (1:maxMoments), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_AuxPolynomialA allocation failed')
      end if

      allocate (op_AuxPolynomialB (1:maxMoments), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_AuxPolynomialB allocation failed')
      end if

      allocate (op_Moments (1:maxMoments), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_Moments allocation failed')
      end if
!
!
!   ...Allocate the working arrays of size op_maxWork. Their size will be set here and will
!      be tentatively 3 times the maximum number of moments expected. As these working arrays
!      will be used to evaluate moments by backwards recursion, they must accomodate more
!      moments than the ones finally wanted.
!
!
      op_maxWork = 3 * maxMoments

      allocate (op_work1 (1:op_maxWork), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_work1 allocation failed')
      end if

      allocate (op_work2 (1:op_maxWork), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_work2 allocation failed')
      end if

      allocate (op_work3 (1:op_maxWork), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_work3 allocation failed')
      end if

  end if
!
!
!    ...Do everything for the Romberg integration (if needed).
!
!
  call RuntimeParameters_get ("opacity_useRomberg",  op_useRomberg)

  if (op_useRomberg) then

      call RuntimeParameters_get ("opacity_RombergAccuracy",    op_RombergAccuracy)

      allocate (op_RombergIntegral (0:op_maxRombergSteps), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_RombergIntegral alocation failed')
      end if

      allocate (op_RombergRow (0:op_maxRombergSteps,OLD:NEW), stat = status)

      if (status > 0) then
          call Driver_abortFlash ('[op_initIntegrate] ERROR: op_RombergRow allocation failed')
      end if

  end if
!
!
!    ...Set initialization status.
!
!
  op_initializedIntegrate = .true.
!
!
!   ...Ready! 
!
!
  return
end subroutine op_initIntegrate
