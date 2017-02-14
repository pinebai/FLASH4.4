!!****if* source/Grid/GridParticles/GridParticlesMapToMesh/Paramesh/gr_ptProlongSmear
!!
!! NAME
!!  
!!
!! SYNOPSIS
!!
!!  gr_ptProlongSmear(integer,intent(IN) :: ioff, &
!!                      integer,intent(IN) :: joff, &
!!                      integer,intent(IN) :: koff, &
!!                      integer,dimension(1:2,1:2,1:2),intent(OUT) :: prolongedSection)
!!
!!
!! DESCRIPTION
!!
!! This routine is a wrapper around the FLASH2 prolongation subroutine. 
!! It uses a single internal cell from the source block, along with its low / high 
!! neighbors in each dimension, to approximate the value that should exist in 
!! the cells that occupy the same spatial location as the single internal cell 
!! but on a finer mesh.  This technique is known as prolongation (interpolation) 
!!
!!
!! ARGUMENTS
!!               ioff:  The IAXIS coordinate of the internal cell.
!!               joff:  The JAXIS coordinate of the internal cell.
!!               koff:  The KAXIS coordinate of the internal cell.
!!               prolongedSection: An array containing the interpolated data.
!!               (This is the data that needs to be accumulated on the finer mesh).
!!
!!***

subroutine gr_ptProlongSmear(ioff,joff,koff,prolongedSection)

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  use gr_ptData, ONLY : gr_ptBuf

  implicit none

  integer, INTENT(in) :: ioff, joff, koff
  real,  dimension(1:2, 1:2, 1:2), INTENT(out) :: prolongedSection

  interface 
     subroutine prolong_temp_patch (src, ioff, joff, koff, &
          ia, ib, ja, jb, ka, kb, recv)

       use Grid_data, ONLY : gr_iloGc, gr_ihiGc, gr_jloGc, gr_jhiGc, gr_kloGc, gr_khiGc
       implicit none
       real, dimension(gr_iLoGc:gr_iHiGc,gr_jLoGc:gr_jHiGc,gr_kLoGc:gr_kHiGc), intent(IN) :: src
       integer, intent(IN) :: ioff, joff, koff, ia, ib, ja, jb, ka, kb
       real, intent(OUT)    :: recv(ia:ib, ja:jb, ka:kb)

     end subroutine prolong_temp_patch
  end interface


  call prolong_temp_patch (gr_ptBuf, ioff, joff, koff, & 
       1, 2, 1, 2, 1, 2, prolongedSection)

  return

end subroutine gr_ptProlongSmear




!===============================================================================

! Routine:      prolong_temp_patch

! Description:  Prolongates (interpolates) a section of one array (src) into a
!               section of another (recv).  The source section begins with
!               indices (ioff,joff,koff), and the receiving section runs from
!               (ia,ja,ka) to (ib,jb,kb).  The routine currently implements
!               conservative quadratic interpolation, assuming a refinement
!               factor (REF_FAC) of 2.


subroutine prolong_temp_patch (src, ioff, joff, koff, &
                               ia, ib, ja, jb, ka, kb, recv)


!-------------------------------------------------------------------------------

  use Grid_data, ONLY : gr_iloGc, gr_ihiGc, gr_jloGc, gr_jhiGc, gr_kloGc, gr_khiGc

  implicit none

#include "constants.h"
#include "Flash.h"
#include "gr_ptMapToMesh.h"

  real, dimension(gr_iLoGc:gr_iHiGc,gr_jLoGc:gr_jHiGc,gr_kLoGc:gr_kHiGc), intent(IN) :: src
  integer, intent(IN) :: ioff, joff, koff, ia, ib, ja, jb, ka, kb
  real, intent(OUT)    :: recv(ia:ib, ja:jb, ka:kb)  !Edit to the original FLASH2 subroutine.

  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1
  integer :: i, j, k, ico, jco, kco
  integer :: ic, jc, kc, icmin, icmax, jcmin, jcmax, kcmin, kcmax
  real    :: cmin, cmax, fmin, fmax
  real    :: PC(-1:1,0:1) = reshape( (/-0.125, 1.,  0.125, & 
  &                                     0.125, 1., -0.125 /), & 
  &                                   shape = (/ 3, 2 /) )

!-------------------------------------------------------------------------------

! Interpolation loop:  loop over zones in the target (fine, child) block.
! (i,j,k) are the indices of the fine zones.  (ip,jp,kp) are the indices
! of the coarse zone enclosing each fine zone.  ip1 and im1 refer to coarse
! zones offset by one to the right and left, respectively; likewise for
! jp1, jm1, kp1, km1.  ico, jco, and kco indicate which set of coefficients
! to use.

  do k = ka, kb

    kp  = (k-ka)/2 + koff
    kp1 = kp + K3D
    km1 = kp - K3D
    kco = mod(k, 2)

    do j = ja, jb

      jp  = (j-ja)/2 + joff
      jp1 = jp + K2D
      jm1 = jp - K2D
      jco = mod(j, 2)

      do i = ia, ib

        ip  = (i-ia)/2 + ioff
        ip1 = ip + 1
        im1 = ip - 1
        ico = mod(i, 2)

#if NDIM == 1
        recv(i,j,k) = &
     &     PC(0,ico) * PC(0,jco) * PC(0,kco) * src(ip,jp,kp) & 
     &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * src(im1,jp,kp) & 
     &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * src(ip1,jp,kp)
#endif
#if NDIM == 2 
        recv(i,j,k) = &
     &     PC(0,ico) * PC(0,jco) * PC(0,kco) * src(ip,jp,kp) & 
     &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * src(im1,jp,kp) & 
     &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * src(ip1,jp,kp) &
     &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * src(ip,jm1,kp) & 
     &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * src(ip,jp1,kp) & 
     &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * src(im1,jm1,kp) & 
     &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * src(ip1,jm1,kp) & 
     &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * src(im1,jp1,kp) & 
     &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * src(ip1,jp1,kp)
#endif
#if NDIM == 3  
        recv(i,j,k) = &
     &     PC(0,ico) * PC(0,jco) * PC(0,kco) * src(ip,jp,kp) & 
     &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * src(im1,jp,kp) & 
     &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * src(ip1,jp,kp) &
     &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * src(ip,jm1,kp) & 
     &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * src(ip,jp1,kp) & 
     &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * src(im1,jm1,kp) & 
     &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * src(ip1,jm1,kp) & 
     &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * src(im1,jp1,kp) & 
     &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * src(ip1,jp1,kp) &
     &   + PC(0,ico) * PC(0,jco) * PC(-1,kco)* src(ip,jp,km1) & 
     &   + PC(-1,ico)* PC(0,jco) * PC(-1,kco)* src(im1,jp,km1) & 
     &   + PC(1,ico) * PC(0,jco) * PC(-1,kco)* src(ip1,jp,km1) & 
     &   + PC(0,ico) * PC(-1,jco)* PC(-1,kco)* src(ip,jm1,km1) & 
     &   + PC(0,ico) * PC(1,jco) * PC(-1,kco)* src(ip,jp1,km1) & 
     &   + PC(-1,ico)* PC(-1,jco)* PC(-1,kco)* src(im1,jm1,km1) & 
     &   + PC(1,ico) * PC(-1,jco)* PC(-1,kco)* src(ip1,jm1,km1) & 
     &   + PC(-1,ico)* PC(1,jco) * PC(-1,kco)* src(im1,jp1,km1) & 
     &   + PC(1,ico) * PC(1,jco) * PC(-1,kco)* src(ip1,jp1,km1) & 
     &   + PC(0,ico) * PC(0,jco) * PC(1,kco) * src(ip,jp,kp1) & 
     &   + PC(-1,ico)* PC(0,jco) * PC(1,kco) * src(im1,jp,kp1) & 
     &   + PC(1,ico) * PC(0,jco) * PC(1,kco) * src(ip1,jp,kp1) & 
     &   + PC(0,ico) * PC(-1,jco)* PC(1,kco) * src(ip,jm1,kp1) & 
     &   + PC(0,ico) * PC(1,jco) * PC(1,kco) * src(ip,jp1,kp1) & 
     &   + PC(-1,ico)* PC(-1,jco)* PC(1,kco) * src(im1,jm1,kp1) & 
     &   + PC(1,ico) * PC(-1,jco)* PC(1,kco) * src(ip1,jm1,kp1) & 
     &   + PC(-1,ico)* PC(1,jco) * PC(1,kco) * src(im1,jp1,kp1) & 
     &   + PC(1,ico) * PC(1,jco) * PC(1,kco) * src(ip1,jp1,kp1)
#endif

      enddo
    enddo
  enddo

!-------------------------------------------------------------------------------

! Monotonicity checking:  force interpolants to be flat if fine-zone averages
! fall outside the range of values of coarse-zone averages used to constrain
! them.

  icmin = ioff
  icmax = (ib-ia)/2 + ioff
  jcmin = joff
  jcmax = (jb-ja)/2 + joff
  kcmin = koff
  kcmax = (kb-ka)/2 + koff

  k = ka - 2
  do kc = kcmin, kcmax
    k = k + 2
    j = ja - 2
    do jc = jcmin, jcmax
      j = j + 2
      i = ia - 2
      do ic = icmin, icmax
        i = i + 2

        cmin = minval(src(ic-1:ic+1,jc-K2D:jc+K2D,kc-K3D:kc+K3D))
        cmax = maxval(src(ic-1:ic+1,jc-K2D:jc+K2D,kc-K3D:kc+K3D))
        fmin = minval(recv(i:i+1,j:j+K2D,k:k+K3D))
        fmax = maxval(recv(i:i+1,j:j+K2D,k:k+K3D))

        if ( (fmin < cmin) .or. (fmax > cmax) ) then
          do kp = k, k+K3D
            do jp = j, j+K2D
              do ip = i, i+1
                recv(ip,jp,kp) = src(ic,jc,kc)
              enddo
            enddo
          enddo
        endif

      enddo
    enddo
  enddo

  return

!-------------------------------------------------------------------------------

end subroutine prolong_temp_patch
