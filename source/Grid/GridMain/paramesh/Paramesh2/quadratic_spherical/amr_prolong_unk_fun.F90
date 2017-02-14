!*******************************************************************************

! Routine:      amr_prolong_unk_fun

! Description:

! This routine takes data from the array recv, originally extracted
! from the solution array unk, and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written directly into block
! isg, which is one of its children. The position of the child within the
! parent block is specified by the ioff, joff and koff arguments.
! The argument jface allows the call to limit its effect to a specific
! face of the block if required. If jface is set to a value between 1 and
! 6, then guard cells for that face are set. If jface is not between 1 to
! 6 then the prolongation operation is applied to the whole block.

! Written:     Peter MacNeice          January 1997
! Revised:     Paul Ricker             March 2001


subroutine amr_prolong_unk_fun & 
     &     (recv, isg, ioff, joff, koff, jface, mype, isrc)
  
  !===============================================================================
  
  use physicaldata
  use tree
  use paramesh_interfaces, ONLY : amr_prolong_gen_unk_fun
  implicit none
  include 'mpif.h'
  
  
  integer :: isg, ioff, joff, koff, jface, mype
  real    :: recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
  integer,intent(IN),OPTIONAL :: isrc
  
  integer :: icl, icu, jcl, jcu, kcl, kcu
  
!===============================================================================
  
  ! Set the bounds on the loop controlling the interpolation.  If jface is not
  ! in the range 1-7, interpolate on all faces - really, on the whole block
  ! including guardcells.  If jface is 7, interpolate in interior cells (i.e.,
  ! non- guardcells) only.  Otherwise interpolate only the face indicated by
  ! jface.
  
  if ( (jface < 1) .or. (jface > 7)) then      ! all faces
   
     icl = il_bnd
     icu = iu_bnd
     jcl = jl_bnd
     jcu = ju_bnd
     kcl = kl_bnd
     kcu = ku_bnd
     
  else
     
     icl = 1+nguard
     icu = nxb+nguard
     jcl = 1+nguard*k2d
     jcu = nyb+nguard*k2d
     kcl = 1+nguard*k3d
     kcu = nzb+nguard*k3d
     
     select case (jface)                         ! one face only for cases 1-6
     case (1)
        icl = il_bnd
        icu = il_bnd + (nguard-1)
     case (2)
        icl = iu_bnd - (nguard-1)
        icu = iu_bnd
     case (3)
        jcl = jl_bnd
        jcu = jl_bnd + (nguard-1)
     case (4)
        jcl = ju_bnd - (nguard-1)
        jcu = ju_bnd
     case (5)
        kcl = kl_bnd
        kcu = kl_bnd + (nguard-1)
     case (6)
        kcl = ku_bnd - (nguard-1)
        kcu = ku_bnd
     end select
     
  endif
      
!===============================================================================

! Call the generic interpolation function.

      call amr_prolong_gen_unk_fun (recv, icl, icu, jcl, jcu, & 
     &                              kcl, kcu, isg, ioff, joff, & 
     &                              koff, mype)

!===============================================================================
      
      return
end subroutine amr_prolong_unk_fun

