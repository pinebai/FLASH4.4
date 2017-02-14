!*******************************************************************************

! Routine:      amr_prolong_work_fun

! Description:

! This routine takes data from the array recv1, originally extracted
! from the workspace array work, and performs a prolongation
! operation on it, between the bounds ranges ia to ib, ja to jb, and
! ka to kb. The data in recv1 is from a parent block and the
! result of the prolongation operation is written directly into block
! isg, which is one of its children. The position of the child within the
! parent block is specified by the ioff, joff and koff arguments.
! The argument jface allows the call to limit its effect to a specific
! face of the block if required. If jface is set to a value between 1 and
! 6, then guard cells for that face are set. If jface is not between 1 to
! 6 then the prolongation operation is applied to the whole block.

! Written:     Peter MacNeice          January 1997
! Revised:     Paul Ricker             March 2001


      subroutine amr_prolong_work_fun & 
     &       (ia, ib, ja, jb, ka, kb, nlayers, isg, & 
     &        ioff, joff, koff, jface, mype)

!===============================================================================

use physicaldata
      use tree
      use workspace
      implicit none
      include 'mpif.h'


      integer :: ia, ib, ja, jb, ka, kb, nlayers
      integer :: isg, ioff, joff, koff, jface, mype

      integer :: icl, icu, jcl, jcu, kcl, kcu

!===============================================================================

! Set the bounds on the loop controlling the interpolation.  If jface is not
! in the range 1-6, interpolate on all faces, otherwise interpolate only the
! face indicated by jface.

      if ( (jface < 1) .or. (jface > 6)) then      ! all faces

        icl = ia + nguard
        icu = ib - nguard
        jcl = ja + nguard
        jcu = jb - nguard
        kcl = ka + nguard
        kcu = kb - nguard

      else                                         ! one face only

        icl = 1+nguard_work
        icu = nxb+nguard_work
        jcl = 1+nguard_work*k2d
        jcu = nyb+nguard_work*k2d
        kcl = 1+nguard_work*k3d
        kcu = nzb+nguard_work*k3d

        select case (jface)
          case (1)
            icl = nguard_work+1-nlayers
            icu = nguard_work
          case (2)
            icl = nxb+nguard_work+1
            icu = nxb+nguard_work+nlayers
          case (3)
            jcl = nguard_work+1-nlayers
            jcu = nguard_work
          case (4)
            jcl = nyb+nguard_work+1
            jcu = nyb+nguard_work+nlayers
          case (5)
            kcl = nguard_work+1-nlayers
            kcu = nguard_work
          case (6)
            kcl = nzb+nguard_work+1
            kcu = nzb+nguard_work+nlayers
        end select

      endif

!===============================================================================

! Call the generic interpolation function.

      call amr_prolong_gen_work_fun (icl, icu, jcl, jcu, & 
     &                               kcl, kcu, isg, ioff, joff, & 
     &                               koff, mype)

!===============================================================================

      return
      end

