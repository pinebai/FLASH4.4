!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMomBins2Dspherical
!!
!! NAME
!!
!!  gr_mpoleMomBins2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMomBins2Dspherical (integer (in) :: maxQtype)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in 2D spherical
!!  geometry. On output, the gr_mpoleMomentR and gr_mpoleMomentI arrays contain
!!  the mass moments over the regular and irregular solid harmonics. The moments
!!  are evaluated by an outer (threaded) loop over all radial bin types, ensuring
!!  separate evaluation of different radial bin types for all threads.
!!
!! ARGUMENTS
!!
!!  maxQtype : the total number of different radial bin types
!!
!!***

subroutine gr_mpoleMomBins2Dspherical (maxQtype)

  use gr_mpoleData,      ONLY : gr_mpoleSymmetryPlane2D,        &
                                gr_mpoleNumberInv,              &
                                gr_mpoleMaxL,                   &
                                gr_mpoleMaxM,                   &
                                gr_mpoleMaxLM,                  &
                                gr_mpoleMaxQ,                   &
                                gr_mpoleQ,                      &
                                gr_mpoleQnumberOfCells,         &
                                gr_mpoleQdataCells2D

  use gr_mpoleData,      ONLY : gr_mpoleQDampingR,              &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: maxQtype

  integer :: L
  integer :: nC,nQ
  integer :: nCells
  integer :: Q

  real    :: cellMass
  real    :: f,g,h
  real    :: Ic0, Ic1, Ic2, IcL
  real    :: r, rR, rI, rsqR, rinvI, rsqinvI
  real    :: Rc0, Rc1, Rc2, RcL
  real    :: Rdamping, Idamping
  real    :: z, zR, zI
!
!
!        ...Only one thread initializes the shared moment arrays.
!
!
!$omp single
  gr_mpoleMomentR (:,1:gr_mpoleMaxQ) = ZERO
  gr_mpoleMomentI (:,1:gr_mpoleMaxQ) = ZERO
!$omp end single
!
!
!        ...Outer (threaded) loop over all different radial bin types.
!           The 2D spherical case. Only a subset of moments for:
!
!                          0 =< L =< gr_mpoleMaxL
!                               M = 0 (cosine part)
!
!             are calculated and summed into the appropriate radial bins.
!             The reason for using only M=0 moments is that the M quantum
!             number controls variations in the azimuthal phi angle in the
!             xy-plane and around the z-axis. But in the 2D spherical case
!             there is no variation around the z-axis. The corresponding
!             M=0 regular and irregular solid harmonics are calculated by
!             using only the M=0 recursion relation.
!
!
!$omp do schedule (dynamic)
  do nQ = 1,maxQtype

     Q        = gr_mpoleQ              (nQ)
     nCells   = gr_mpoleQnumberOfCells (nQ)
     Idamping = gr_mpoleQDampingI      ( Q)
     Rdamping = gr_mpoleQDampingR      ( Q)

     do nC = 1,nCells

        z        = gr_mpoleQdataCells2D (nC,nQ) % coord1
        cellMass = gr_mpoleQdataCells2D (nC,nQ) % cellMass
        r        = gr_mpoleQdataCells2D (nC,nQ) % radius
!
!
!        ...Calculate the (damped) moments.
!
!
        rI = r * Idamping
        rinvI = ONE / rI

        RcL = cellMass                    !  regular  solid harmonic for cos 00
        IcL = cellMass * rinvI            ! irregular solid harmonic for cos 00

        gr_mpoleMomentR (1,Q) = gr_mpoleMomentR (1,Q) + RcL
        gr_mpoleMomentI (1,Q) = gr_mpoleMomentI (1,Q) + IcL
!
!
!        ...Proceed, if we go beyond monopoles.
!
!
        if (gr_mpoleMaxL > 0) then

         zR = z * Rdamping
         rR = r * Rdamping
         zI = z * Idamping

         rsqR = rR * rR
         rsqinvI = rinvI * rinvI

         Rc1 = RcL * zR                    !  regular  solid harmonic for cos 10
         Ic1 = IcL * zI * rsqinvI          ! irregular solid harmonic for cos 10

         gr_mpoleMomentR (2,Q) = gr_mpoleMomentR (2,Q) + Rc1
         gr_mpoleMomentI (2,Q) = gr_mpoleMomentI (2,Q) + Ic1
!
!
!        ...Proceed, if we go beyond dipoles.
!
!
         if (gr_mpoleMaxL > 1) then

          Rc0 = RcL
          Ic0 = IcL

          do L = 2,gr_mpoleMaxL

             h = real (L + L - 1)
             g = real ((L - 1) * (L - 1))
             f = gr_mpoleNumberInv (L) * gr_mpoleNumberInv (L)

             Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f        !  regular  solid harmonic for cos L0
             Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI  ! irregular solid harmonic for cos L0

             gr_mpoleMomentR (L+1,Q) = gr_mpoleMomentR (L+1,Q) + Rc2
             gr_mpoleMomentI (L+1,Q) = gr_mpoleMomentI (L+1,Q) + Ic2

             Rc0 = Rc1                                    !
             Rc1 = Rc2                                    ! swap, to accumulate next L regular
             Ic0 = Ic1                                    ! and irregular solid harmonic
             Ic1 = Ic2                                    !

          end do

         end if    ! gr_mpoleMaxL > 1 condition
        end if     ! gr_mpoleMaxL > 0 condition

     end do   ! loop over cells in radial bin type
  end do      ! loop over radial bin types (threaded)
!$omp end do
!
!
!    ...If a symmetry plane is assumed along the theta = 90 degrees (pi/2) part of
!       the domain, the overall moment values including the omitted symmetry part
!       can be calculated from the existing moments by the following rule:
!
!                      odd L values -> set to zero
!                     even L values -> multiply by 2
!
!       The reason for this rule is that the symmetry plane effectively changes
!       the z-coordinate value of the moments to -z. Since all odd L moments have
!       the structure z*(even function of z) their +z and -z sum cancels. Likewise
!       the even L moments are all even functions of z and thus the +z and -z moments
!       are of eqaul sign and magnitude.
!
!
!$omp single
  if (gr_mpoleSymmetryPlane2D) then

   gr_mpoleMomentR (1:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ) = TWO * gr_mpoleMomentR (1:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ)
   gr_mpoleMomentR (2:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ) = ZERO
   gr_mpoleMomentI (1:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ) = TWO * gr_mpoleMomentI (1:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ)
   gr_mpoleMomentI (2:gr_mpoleMaxLM:2,1:gr_mpoleMaxQ) = ZERO

  end if
!$omp end single
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleMomBins2Dspherical
