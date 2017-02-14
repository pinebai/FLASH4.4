!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMomBins1Dspherical
!!
!! NAME
!!
!!  gr_mpoleMomBins1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMomBins1Dspherical (integer (in) :: maxQtype)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in 1D spherical
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

subroutine gr_mpoleMomBins1Dspherical (maxQtype)

  use gr_mpoleData,      ONLY : gr_mpoleMaxQ,            &
                                gr_mpoleQDampingI,       &
                                gr_mpoleMomentR,         &
                                gr_mpoleMomentI,         &
                                gr_mpoleQ,               &
                                gr_mpoleQnumberOfCells,  &
                                gr_mpoleQdataCells1D

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"
  
  integer, intent (in) :: maxQtype

  integer :: nC,nQ
  integer :: nCells
  integer :: Q

  real    :: cellMass
  real    :: Idamping
  real    :: r
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
!           The 1D spherical case uses only L = 0 moments.
!
!
!$omp do schedule (dynamic)
  do nQ = 1,maxQtype

     Q        = gr_mpoleQ              (nQ)
     nCells   = gr_mpoleQnumberOfCells (nQ)
     Idamping = gr_mpoleQDampingI      ( Q)

     do nC = 1,nCells

        cellMass = gr_mpoleQdataCells1D (nC,nQ) % cellMass
        r        = gr_mpoleQdataCells1D (nC,nQ) % radius
!
!
!        ...Calculate the (damped) moments.
!
!
        gr_mpoleMomentR (1,Q) = gr_mpoleMomentR (1,Q) + cellMass
        gr_mpoleMomentI (1,Q) = gr_mpoleMomentI (1,Q) + cellMass / (r * Idamping)

     end do   ! loop over cells in radial bin type
  end do      ! loop over radial bin types (threaded)
!$omp end do
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleMomBins1Dspherical
