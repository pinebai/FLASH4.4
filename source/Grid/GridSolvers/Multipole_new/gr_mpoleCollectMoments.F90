!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCollectMoments
!!
!! NAME
!!
!!  gr_mpoleCollectMoments
!!
!! SYNOPSIS
!!
!!  gr_mpoleCollectMoments ()
!!
!! DESCRIPTION
!!
!!  The individual moments in each radial bin are collected and distributed to each
!!  processor. Each cell should contribute its regular moment to all cells with radius
!!  greater than it and its irregular moment to all cells with radius less than it.
!!  After exiting this routine, all processors have a copy of the entire moments array.
!!
!! ARGUMENTS
!!
!!  none
!!
!!***

subroutine gr_mpoleCollectMoments ()

  use Grid_data,         ONLY : gr_meshComm

  use gr_mpoleData,      ONLY : gr_mpoleMaxL,                 &
                                gr_mpoleMaxM,                 &
                                gr_mpoleMaxLM,                &
                                gr_mpoleMaxQ,                 &
                                gr_mpoleTotalNrCosineMoments, &
                                gr_mpoleMomentR,              &
                                gr_mpoleMomentI,              &
                                gr_mpoleQDampingR,            &
                                gr_mpoleQDampingI,            &
                                gr_mpoleScratch

  implicit none
  
#include "Flash.h"
#include "gr_mpole.h"

  include "Flash_mpi.h"

  integer :: allMoments
  integer :: c,s
  integer :: error
  integer :: L,M,n
  integer :: Q

  real    :: dampingQuotient
  real    :: dampI, Idamp
  real    :: dampR, Rdamp
!
!
!    ...Handle trivial case of gr_mpoleMaxL = 0 first. For comments about what is
!       done, see below.
!
!
  if (gr_mpoleMaxL == 0) then

        do Q = 2,gr_mpoleMaxQ
           gr_mpoleMomentR (1,Q) = gr_mpoleMomentR (1,Q) + gr_mpoleMomentR (1,Q-1)
        end do

        do Q = gr_mpoleMaxQ-1,1,-1
           dampingQuotient = gr_mpoleQDampingI (Q+1) / gr_mpoleQDampingI (Q)
           gr_mpoleMomentI (1,Q) = gr_mpoleMomentI (1,Q) + dampingQuotient * gr_mpoleMomentI (1,Q+1)
        end do
!
!
!    ...The case when gr_mpoleMaxL > 0. Each cell should contribute its regular
!       moment to all cells with radius greater than it. Rescale the added moments.
!
!
  else

        do Q = 2,gr_mpoleMaxQ

           dampingQuotient = gr_mpoleQDampingR (Q) / gr_mpoleQDampingR (Q-1)

           Rdamp = ONE     ! save for accumulation during M loop below
           dampR = ONE     ! local accumulation over L values

           do L = 0,gr_mpoleMaxL
              gr_mpoleMomentR (L+1,Q) = gr_mpoleMomentR (L+1,Q) + dampR * gr_mpoleMomentR (L+1,Q-1)
              dampR = dampR * dampingQuotient
           end do

           if (gr_mpoleMaxM /= 0) then

               c = gr_mpoleMaxL + 1
               s = gr_mpoleTotalNrCosineMoments

               n = 0
               do M = 1,gr_mpoleMaxM
                  Rdamp = Rdamp * dampingQuotient
                  dampR = Rdamp
                  do L = M,gr_mpoleMaxL
                     n = n + 1
                     gr_mpoleMomentR (c+n,Q) = gr_mpoleMomentR (c+n,Q) + dampR * gr_mpoleMomentR (c+n,Q-1)
                     gr_mpoleMomentR (s+n,Q) = gr_mpoleMomentR (s+n,Q) + dampR * gr_mpoleMomentR (s+n,Q-1)
                     dampR = dampR * dampingQuotient
                  end do
               end do

           end if
        end do
!
!
!    ...Each cell should contribute its irregular moment to all cells with
!       radius less than it. Rescale the added moments.
!
!
        do Q = gr_mpoleMaxQ-1,1,-1

           dampingQuotient = gr_mpoleQDampingI (Q+1) / gr_mpoleQDampingI (Q)

           Idamp = dampingQuotient     ! save for accumulation during M loop below
           dampI = dampingQuotient     ! local accumulation over L values

           do L = 0,gr_mpoleMaxL
              gr_mpoleMomentI (L+1,Q) = gr_mpoleMomentI (L+1,Q) + dampI * gr_mpoleMomentI (L+1,Q+1)
              dampI = dampI * dampingQuotient
           end do

           if (gr_mpoleMaxM /= 0) then

               c = gr_mpoleMaxL + 1
               s = gr_mpoleTotalNrCosineMoments

               n = 0
               do M = 1,gr_mpoleMaxM
                  Idamp = Idamp * dampingQuotient
                  dampI = Idamp
                  do L = M,gr_mpoleMaxL
                     n = n + 1
                     gr_mpoleMomentI (c+n,Q) = gr_mpoleMomentI (c+n,Q) + dampI * gr_mpoleMomentI (c+n,Q+1)
                     gr_mpoleMomentI (s+n,Q) = gr_mpoleMomentI (s+n,Q) + dampI * gr_mpoleMomentI (s+n,Q+1)
                     dampI = dampI * dampingQuotient
                  end do
               end do

           end if
        end do

  end if
!
!
!    ...Make global sum and give all processors a copy of it. Don't
!       forget to zero the extra outer and inner radial bin of the
!       moments.
!
!
  allMoments = gr_mpoleMaxLM * gr_mpoleMaxQ

  call MPI_AllReduce (gr_mpoleMomentR (:,1:gr_mpoleMaxQ), &
                      gr_mpoleScratch,                    &
                      allMoments,                         &
                      FLASH_REAL,                         &
                      MPI_SUM,                            &
                      gr_meshComm,                        &
                      error                               )

  gr_mpoleMomentR (:,0) = ZERO
  gr_mpoleMomentR (:,1:gr_mpoleMaxQ) = gr_mpoleScratch (:,1:gr_mpoleMaxQ)

  call MPI_AllReduce (gr_mpoleMomentI (:,1:gr_mpoleMaxQ), &
                      gr_mpoleScratch,                    &
                      allMoments,                         &
                      FLASH_REAL,                         &
                      MPI_SUM,                            &
                      gr_meshComm,                        &
                      error                               )

  gr_mpoleMomentI (:,1:gr_mpoleMaxQ) = gr_mpoleScratch (:,1:gr_mpoleMaxQ)
  gr_mpoleMomentI (:,gr_mpoleMaxQ+1) = ZERO
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleCollectMoments
