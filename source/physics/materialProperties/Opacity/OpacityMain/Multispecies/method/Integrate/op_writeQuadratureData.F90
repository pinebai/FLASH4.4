!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_writeQuadratureData
!!
!! NAME
!!
!!  op_writeQuadratureData
!!
!! SYNOPSIS
!!
!!  call op_writeQuadratureData (integer (in) :: nRoots,
!!                               integer (in) :: nMoments,
!!                               logical (in) :: Ttiny,
!!                               logical (in) :: Tsmall,
!!                               logical (in) :: Tlarge,
!!                               logical (in) :: Thuge,
!!                               real    (in) :: T,
!!                               real    (in) :: p,
!!                               real    (in) :: q,
!!                               real    (in) :: beta,
!!                               real    (in) :: MomZero,
!!                               real    (in) :: Roots   (1:nRoots),
!!                               real    (in) :: Weights (1:nRoots))
!!
!! DESCRIPTION
!!
!!  Utility routine, which prints detailed info regarding the presently calculated
!!  quadrature rule to a text file. The information is written out to a file
!!  named basenm_opacityQuadratureData.txt, where basenm is the runtime parameter for
!!  output file names. The file is appended at each time for each quadrature.
!!  The routine is mainly meant for checking purpose only.
!!
!! ARGUMENTS
!!
!!  nRoots   : the number of roots and weights calculated
!!  nMoments : the number of moments calculated
!!  Ttiny    : indicator, if the integration limit T was considered tiny
!!  Tsmall   : indicator, if the integration limit T was considered small
!!  Tlarge   : indicator, if the integration limit T was considered large
!!  Thuge    : indicator, if the integration limit T was considered huge
!!  T        : the upper integration limit
!!  p        : the p parameter for the shifted Jacobi polynomials G(p,q,x)
!!  q        : the q parameter for the shifted Jacobi polynomials G(p,q,x)
!!  beta     : the exponent in x^beta for the Laguerre-type weight
!!  MomZero  : the 0-th moment
!!  Roots    : the roots for i = 1,2,...,nRoots
!!  Weights  : the weights for i = 1,2,...,nRoots
!!
!!***
subroutine op_writeQuadratureData (nRoots, nMoments,                  &
                                   Ttiny, Tsmall, Tlarge, Thuge, T,   &
                                   p, q, beta, MomZero, Roots, Weights)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use op_integrateData,            ONLY : op_AuxPolynomialA,     &
                                          op_AuxPolynomialB,     &
                                          op_Moments,            &
                                          op_OrthoPolynomialA,   &
                                          op_OrthoPolynomialB
  use Driver_interface,            ONLY : Driver_abortFlash
  use op_numericsData,             ONLY : one

  implicit none

# include "constants.h"

  integer, intent (in) :: nRoots
  integer, intent (in) :: nMoments
  logical, intent (in) :: Ttiny
  logical, intent (in) :: Tsmall
  logical, intent (in) :: Tlarge
  logical, intent (in) :: Thuge
  real,    intent (in) :: T
  real,    intent (in) :: p,q
  real,    intent (in) :: beta
  real,    intent (in) :: MomZero

  real,    intent (in) :: Roots   (1:nRoots)
  real,    intent (in) :: Weights (1:nRoots)

  logical, save :: firstCall = .true.

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: i
  integer :: posBlank

  character (len=MAX_STRING_LENGTH), save :: baseName
  character (len=MAX_STRING_LENGTH), save :: fileName
!
!
!   ...Do the printout only on the master processor.
!
!
!  if (op_myPE == MASTER_PE) then

      fileUnit = ut_getFreeFileUnit ()

      if (firstCall) then
          call RuntimeParameters_get ("basenm",baseName)
          posBlank = index (baseName,' ')
          fileName = baseName (:posBlank-1) // 'QuadratureData.txt'
          open  (fileUnit, file=fileName)
          firstCall = .false.
      else
          open (fileUnit, file=fileName, position='APPEND')
      end if
!
!
!     ...Printout the quadrature info.
!
!
      write (fileUnit,*)
      write (fileUnit,*) '     QUADRATURE INFO ( GENERALIZED LAGUERRE WEIGHT:  X^BETA * EXP (-X) )'
      write (fileUnit,*)
      write (fileUnit,'(A33,ES14.6)') ' Laguerre-type weight exponent = ',beta
      write (fileUnit,'(A33,ES14.6)') ' Upper integration limit T     = ',T
      write (fileUnit,'(A33,I3)')     ' Number of roots and weights   = ',nRoots
      write (fileUnit,'(A33,I3)')     ' Number of modified moments    = ',nMoments

      if (Ttiny) then

      else if (Tsmall) then

           write (fileUnit,'(A33,A23)')    ' Auxilliary polynomial type    = ','shifted Jacobi G(p,q,x)'
           write (fileUnit,'(A33,ES14.6)') ' shJacobi polynomial p         = ',p
           write (fileUnit,'(A33,ES14.6)') ' shJacobi polynomial q         = ',q

      else if (Tlarge) then

           write (fileUnit,'(A33,A31)')    ' Auxilliary polynomial type    = ','generalized Laguerre L(alpha,x)'
           write (fileUnit,'(A33,ES14.6)') ' Laguerre polynomial alpha     = ',beta

      else if (Thuge) then

      end if
!
!
!     ...Modified moment values (if any).
!
!
      if (Tsmall .or. Tlarge) then

          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,*)
          write (fileUnit,'(A35)') '             MODIFIED MOMENTS      '
          write (fileUnit,*)
          write (fileUnit,'(A35)') '   i              Moments          '
          write (fileUnit,'(A35)') ' ----------------------------------'

          if (Tsmall) then
              write (fileUnit,'(I4, 10X, ES20.12,A24)') 0,MomZero,' (Integration limit = 1)'
          else if (Tlarge) then
              write (fileUnit,'(I4, 10X, ES20.12,A24)') 0,MomZero,' (Integration limit = T)'
          end if

          do i = 1,nMoments
             write (fileUnit,'(I4,10X,ES20.12)') i, op_Moments (i)
          end do

      end if
!
!
!     ...Auxilliary and/or orthogonal polynomial 3-term recursion coefficients.
!
!
      if (Ttiny .or. Thuge) then

        write (fileUnit,*)
        write (fileUnit,*)
        write (fileUnit,*)
        write (fileUnit,'(A45)') ' ORTHOGONAL POLYNOMIAL RECURSION COEFFICIENTS'
        write (fileUnit,*)
        write (fileUnit,'(A45)') '   i     Orthogonal (A)      Orthogonal (B)  '
        write (fileUnit,'(A45)') ' --------------------------------------------'

        write (fileUnit,'(I4, ES20.12, 20X)') 1, op_OrthoPolynomialA (1)
        do i = 2,nRoots
           write (fileUnit,'(I4, ES20.12, ES20.12)') i, op_OrthoPolynomialA (i), op_OrthoPolynomialB (i)
        end do

      else if (Tsmall .or. Tlarge) then

        write (fileUnit,*)
        write (fileUnit,*)
        write (fileUnit,*)
        write (fileUnit,'(A85)') '              AUXILLIARY / ORTHOGONAL POLYNOMIAL RECURSION COEFFICIENTS              '
        write (fileUnit,*)
        write (fileUnit,'(A85)') '   i    Auxilliary (A)      Auxilliary (B)       Orthogonal (A)      Orthogonal (B)  '
        write (fileUnit,'(A85)') ' ------------------------------------------------------------------------------------'

        write (fileUnit,'(I4, ES20.12, 20X, ES20.12, 20X)') &
                           1, op_AuxPolynomialA (1), op_OrthoPolynomialA (1)
        do i = 2,nRoots
           write (fileUnit,'(I4, ES20.12, ES20.12, ES20.12, ES20.12)')             &
                              i, op_AuxPolynomialA   (i), op_AuxPolynomialB   (i), &
                                 op_OrthoPolynomialA (i), op_OrthoPolynomialB (i)
        end do

      end if
!
!
!     ...The roots and weights.
!
!
      write (fileUnit,*)
      write (fileUnit,*)
      write (fileUnit,*)
      write (fileUnit,'(A45)') '             ROOTS AND WEIGHTS               '
      write (fileUnit,*)
      write (fileUnit,'(A45)') '   i         Roots              Weights      '
      write (fileUnit,'(A45)') ' --------------------------------------------'

      do i = 1,nRoots
         write (fileUnit,'(I4,   ES20.12,    ES20.12)') &
                            i, Roots (i), Weights (i)
      end do

      write (fileUnit,*)
      write (fileUnit,'(A38)') '  ### finished present quadrature ###  '
      close (fileUnit)

!  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_writeQuadratureData
