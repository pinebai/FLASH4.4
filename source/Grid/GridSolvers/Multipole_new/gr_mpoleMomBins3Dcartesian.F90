!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMomBins3Dcartesian
!!
!! NAME
!!
!!  gr_mpoleMomBins3Dcartesian
!!
!! SYNOPSIS
!!
!!  gr_mpoleMomBins3Dcartesian (integer (in) :: maxQtype)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in 3D cartesian
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

subroutine gr_mpoleMomBins3Dcartesian (maxQtype)

  use gr_mpoleData,      ONLY : gr_mpoleSymmetryAxis3D,         &
                                gr_mpoleNumberInv,              &
                                gr_mpoleTotalNrCosineMoments,   &
                                gr_mpoleMaxL,                   &
                                gr_mpoleMax2L,                  &
                                gr_mpoleMaxM,                   &
                                gr_mpoleMaxLM,                  &
                                gr_mpoleMaxQ,                   &
                                gr_mpoleQ,                      &
                                gr_mpoleQnumberOfCells,         &
                                gr_mpoleQdataCells3D

  use gr_mpoleData,      ONLY : gr_mpoleQDampingR,              &
                                gr_mpoleQDampingI,              &
                                gr_mpoleMomentR,                &
                                gr_mpoleMomentI

  implicit none
  
#include "Flash.h"
#include "constants.h"
#include "gr_mpole.h"

  integer, intent (in) :: maxQtype

  integer :: c,s
  integer :: M,MM,L
  integer :: nC, nQ
  integer :: nCells
  integer :: Q

  real    :: cellMass
  real    :: f,g,h
  real    :: Ic0, Is0, Ic1, Is1, Ic2, Is2, IcL, IsL
  real    :: r, rR, rI, rsqR, rinvI, rsqinvI
  real    :: Rc0, Rs0, Rc1, Rs1, Rc2, Rs2, RcL, RsL
  real    :: Rdamping, Idamping
  real    :: x,y,z
  real    :: xI, yI, zI
  real    :: xR, yR, zR
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
!           The 3D cartesian case. The full set of moments for:
!
!                          0 =< L =< gr_mpoleMaxL
!                          0 =< M =< L (cosine part)
!                          1 =< M =< L (sine part)
!
!           are calculated and summed into the appropriate radial bins.
!
!           If an axisymmetric symmetry has been specified, one enforces
!           computationally an axial symmetry on the problem, i.e. the problem
!           is treated as if axisymmetry is present (which in the real simulation
!           is never the case due to the finite resolution of the grid).
!           Axisymmetry means rotational invariance around the z-axis and hence
!           we evaluate only the M = 0 cosine components:
!
!                             0 =< L =< gr_mpoleMaxL
!                                  M = 0 (cosine part)
!
!
!$omp do schedule (dynamic)
  do nQ = 1,maxQtype

     Q        = gr_mpoleQ              (nQ)
     nCells   = gr_mpoleQnumberOfCells (nQ)
     Idamping = gr_mpoleQDampingI      ( Q)
     Rdamping = gr_mpoleQDampingR      ( Q)

     do nC = 1,nCells

        x        = gr_mpoleQdataCells3D (nC,nQ) % coord1
        y        = gr_mpoleQdataCells3D (nC,nQ) % coord2
        z        = gr_mpoleQdataCells3D (nC,nQ) % coord3
        cellMass = gr_mpoleQdataCells3D (nC,nQ) % cellMass
        r        = gr_mpoleQdataCells3D (nC,nQ) % radius
!
!
!        ...Calculate the (damped) moments.
!
!
        rI = r * Idamping
        rinvI = ONE / rI

        RcL = cellMass
        IcL = cellMass * rinvI

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

         Rc0 = RcL                         !  regular  solid harmonic for cos 00
         Ic0 = IcL                         ! irregular solid harmonic for cos 00
         Rc1 = RcL * zR                    !  regular  solid harmonic for cos 10
         Ic1 = IcL * zI * rsqinvI          ! irregular solid harmonic for cos 10

         gr_mpoleMomentR (2,Q) = gr_mpoleMomentR (2,Q) + Rc1
         gr_mpoleMomentI (2,Q) = gr_mpoleMomentI (2,Q) + Ic1

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

         if (.not. gr_mpoleSymmetryAxis3D) then

          xR = x * Rdamping
          yR = y * Rdamping
          xI = x * Idamping
          yI = y * Idamping

          c = gr_mpoleMaxL + 2
          s = gr_mpoleTotalNrCosineMoments + 1

          Rc0 = - RcL * xR * HALF                 !  regular  solid harmonic for cos 11
          Rs0 = - RcL * yR * HALF                 !  regular  solid harmonic for sin 11
          Ic0 = - IcL * xI * rsqinvI              ! irregular solid harmonic for cos 11
          Is0 = - IcL * yI * rsqinvI              ! irregular solid harmonic for sin 11

          gr_mpoleMomentR (c,Q) = gr_mpoleMomentR (c,Q) + Rc0
          gr_mpoleMomentR (s,Q) = gr_mpoleMomentR (s,Q) + Rs0
          gr_mpoleMomentI (c,Q) = gr_mpoleMomentI (c,Q) + Ic0
          gr_mpoleMomentI (s,Q) = gr_mpoleMomentI (s,Q) + Is0
!
!
!        ...Proceed, if we go beyond dipoles.
!
!
          if (gr_mpoleMaxL > 1) then

           RcL = Rc0                               ! store cos 11 for later (22 case)
           RsL = Rs0                               ! store sin 11 for later (22 case)
           IcL = Ic0                               ! store cos 11 for later (22 case)
           IsL = Is0                               ! store sin 11 for later (22 case)

           Rc1 =         zR * RcL                  !  regular  solid harmonic for cos 21
           Rs1 =         zR * RsL                  !  regular  solid harmonic for sin 21
           Ic1 = THREE * zI * IcL * rsqinvI        ! irregular solid harmonic for cos 21
           Is1 = THREE * zI * IsL * rsqinvI        ! irregular solid harmonic for sin 21

           gr_mpoleMomentR (c+1,Q) = gr_mpoleMomentR (c+1,Q) + Rc1
           gr_mpoleMomentR (s+1,Q) = gr_mpoleMomentR (s+1,Q) + Rs1
           gr_mpoleMomentI (c+1,Q) = gr_mpoleMomentI (c+1,Q) + Ic1
           gr_mpoleMomentI (s+1,Q) = gr_mpoleMomentI (s+1,Q) + Is1

           do L = 2,gr_mpoleMaxL-1                 ! M=1 and L=3,gr_mpoleMaxL (shifted L-loop!)

              h = real (L + L + 1)
              g = real ((L + 1) * (L - 1))
              f = gr_mpoleNumberInv (L + 2) * gr_mpoleNumberInv (L)

              Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f        !  regular  solid harmonic for cos L1
              Rs2 = (h * zR * Rs1 - rsqR * Rs0) * f        !  regular  solid harmonic for sin L1
              Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI  ! irregular solid harmonic for cos L1
              Is2 = (h * zI * Is1 -    g * Is0) * rsqinvI  ! irregular solid harmonic for sin L1

              gr_mpoleMomentR (c+L,Q) = gr_mpoleMomentR (c+L,Q) + Rc2
              gr_mpoleMomentR (s+L,Q) = gr_mpoleMomentR (s+L,Q) + Rs2
              gr_mpoleMomentI (c+L,Q) = gr_mpoleMomentI (c+L,Q) + Ic2
              gr_mpoleMomentI (s+L,Q) = gr_mpoleMomentI (s+L,Q) + Is2

              Rc0 = Rc1                                    !
              Rs0 = Rs1                                    !
              Rc1 = Rc2                                    !
              Rs1 = Rs2                                    ! swap, to accumulate next regular
              Ic0 = Ic1                                    ! and irregular solid harmonic
              Is0 = Is1                                    !
              Ic1 = Ic2                                    !
              Is1 = Is2                                    !

           end do

           c = c + gr_mpoleMaxL
           s = s + gr_mpoleMaxL

           do M = 2,gr_mpoleMaxL-1

              MM = M + M
              h = real (MM + 1)
              g = real (MM - 1)
              f = gr_mpoleNumberInv (MM)

              Rc0 =       (yR * RsL - xR * RcL) * f            !  regular  solid harmonic for cos MM
              Rs0 =     - (yR * RcL + xR * RsL) * f            !  regular  solid harmonic for sin MM
              Ic0 =   g * (yI * IsL - xI * IcL) * rsqinvI      ! irregular solid harmonic for cos MM
              Is0 = - g * (yI * IcL + xI * IsL) * rsqinvI      ! irregular solid harmonic for sin MM

              gr_mpoleMomentR (c,Q) = gr_mpoleMomentR (c,Q) + Rc0
              gr_mpoleMomentR (s,Q) = gr_mpoleMomentR (s,Q) + Rs0
              gr_mpoleMomentI (c,Q) = gr_mpoleMomentI (c,Q) + Ic0
              gr_mpoleMomentI (s,Q) = gr_mpoleMomentI (s,Q) + Is0

              RcL = Rc0                                        ! store cos MM for M+1,M+1
              RsL = Rs0                                        ! store sin MM for M+1,M+1
              IcL = Ic0                                        ! store cos MM for M+1,M+1
              IsL = Is0                                        ! store sin MM for M+1,M+1

              Rc1 =     zR * RcL                               !  regular  solid harmonic for cos L=M+1,M
              Rs1 =     zR * RsL                               !  regular  solid harmonic for sin L=M+1,M
              Ic1 = h * zI * IcL * rsqinvI                     ! irregular solid harmonic for cos L=M+1,M
              Is1 = h * zI * IsL * rsqinvI                     ! irregular solid harmonic for sin L=M+1,M

              gr_mpoleMomentR (c+1,Q) = gr_mpoleMomentR (c+1,Q) + Rc1
              gr_mpoleMomentR (s+1,Q) = gr_mpoleMomentR (s+1,Q) + Rs1
              gr_mpoleMomentI (c+1,Q) = gr_mpoleMomentI (c+1,Q) + Ic1
              gr_mpoleMomentI (s+1,Q) = gr_mpoleMomentI (s+1,Q) + Is1

              do L = 2,gr_mpoleMaxL-M                          ! L=M+2,gr_mpoleMaxL (shifted L-loop!)

                 h = real (L + L + MM - 1)
                 g = real (L + MM - 1) * (L - 1)
                 f = gr_mpoleNumberInv (L + MM) * gr_mpoleNumberInv (L)

                 Rc2 = (h * zR * Rc1 - rsqR * Rc0) * f         !  regular  solid harmonic for cos LM
                 Rs2 = (h * zR * Rs1 - rsqR * Rs0) * f         !  regular  solid harmonic for sin LM
                 Ic2 = (h * zI * Ic1 -    g * Ic0) * rsqinvI   ! irregular solid harmonic for cos LM
                 Is2 = (h * zI * Is1 -    g * Is0) * rsqinvI   ! irregular solid harmonic for sin LM

                 gr_mpoleMomentR (c+L,Q) = gr_mpoleMomentR (c+L,Q) + Rc2
                 gr_mpoleMomentR (s+L,Q) = gr_mpoleMomentR (s+L,Q) + Rs2
                 gr_mpoleMomentI (c+L,Q) = gr_mpoleMomentI (c+L,Q) + Ic2
                 gr_mpoleMomentI (s+L,Q) = gr_mpoleMomentI (s+L,Q) + Is2

                 Rc0 = Rc1                                     !
                 Rs0 = Rs1                                     !
                 Rc1 = Rc2                                     !
                 Rs1 = Rs2                                     ! swap, to accumulate next regular
                 Ic0 = Ic1                                     ! and irregular solid harmonic
                 Is0 = Is1                                     !
                 Ic1 = Ic2                                     !
                 Is1 = Is2                                     !

              end do

              c = c + gr_mpoleMaxL - M + 1
              s = s + gr_mpoleMaxL - M + 1

           end do

           g = real (gr_mpoleMax2L - 1)
           f = gr_mpoleNumberInv (gr_mpoleMax2L)

           Rc0 =       (yR * RsL - xR * RcL) * f               !  regular  solid harmonic for cos MaxL,MaxL
           Rs0 =     - (yR * RcL + xR * RsL) * f               !  regular  solid harmonic for sin MaxL,MaxL
           Ic0 =   g * (yI * IsL - xI * IcL) * rsqinvI         ! irregular solid harmonic for cos MaxL,MaxL
           Is0 = - g * (yI * IcL + xI * IsL) * rsqinvI         ! irregular solid harmonic for sin MaxL,MaxL

           gr_mpoleMomentR (c,Q) = gr_mpoleMomentR (c,Q) + Rc0
           gr_mpoleMomentR (s,Q) = gr_mpoleMomentR (s,Q) + Rs0
           gr_mpoleMomentI (c,Q) = gr_mpoleMomentI (c,Q) + Ic0
           gr_mpoleMomentI (s,Q) = gr_mpoleMomentI (s,Q) + Is0

          end if   ! gr_mpoleMaxL > 1       condition
         end if    ! gr_mpoleSymmetryAxis3D condition
        end if     ! gr_mpoleMaxL > 0       condition

     end do   ! loop over cells in radial bin type
  end do      ! loop over radial bin types (threaded)
!$omp end do
!
!
!        ...Ready!
!
!
  return
end subroutine gr_mpoleMomBins3Dcartesian
