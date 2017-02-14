!*******************************************************************************

! Routine:      amr_prolong_gen_work_fun

! Description:

! This routine takes data from the array recv1, originally extracted
! from the workspace array work on some block, and performs a prolongation
! operation on it, between the bounds ranges ia to ib, ja to jb, and ka to kb.
! The data in recv1 is from a parent block and the result of the prolongation
! operation is returned in outp.  The position of the child within the parent
! block is specified by the ioff, joff and koff arguments.

! This particular prolongation is conservative, cell-averaged, triquadratic
! interpolation. It can be used for blocks with an even or odd number of grid
! cells.


subroutine amr_prolong_gen_work_fun & 
     &       (ia, ib, ja, jb, ka, kb, lb, ioff, joff, koff, mype)
  
  !===============================================================================
  
  use Grid_data, ONLY: gr_monotone
  use physicaldata
  use tree
  use workspace
  implicit none
  include 'mpif.h'
  

  
  integer :: ia, ib, ja, jb, ka, kb, lb
  integer :: ioff, joff, koff, mype
  
  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1
  integer :: i, j, k, ico, jco, kco
  integer :: ic, jc, kc, icmin, icmax, jcmin, jcmax
  integer :: kcmin, kcmax
  real    :: cmin, cmax, fmin, fmax
  real    :: PC(-1:1,0:1) = reshape( (/-0.125, 1.,  0.125, & 
       &                                      0.125, 1., -0.125 /), & 
       &                                   shape = (/ 3, 2 /) )

!===============================================================================

  ! Interpolation loop:  loop over zones in the target (fine, child) block.
  ! (i,j,k) are the indices of the fine zones.  (ip,jp,kp) are the indices
  ! of the coarse zone enclosing each fine zone.  ip1 and im1 refer to coarse
  ! zones offset by one to the right and left, respectively; likewise for
  ! jp1, jm1, kp1, km1.  ico, jco, and kco indicate which set of coefficients
  ! to use.
  
  do k = ka, kb
     
     kp  = ((k-1)/2) + 1 + (nguard_work/2)*k3d + koff
     kp1 = kp + k3d
     km1 = kp - k3d
     kco = mod(k, 2)
     
     do j = ja, jb
        
        jp  = ((j-1)/2) + 1 + (nguard_work/2)*k2d + joff
        jp1 = jp + k2d
        jm1 = jp - k2d
        jco = mod(j, 2)
        
        do i = ia, ib
           
           ip  = ((i-1)/2) + 1 + (nguard_work/2) + ioff
           ip1 = ip + 1
           im1 = ip - 1
           ico = mod(i, 2)
           
#if N_DIM == 1
           work(i,j,k,lb,1) = &
                &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv1(ip,jp,kp) & 
                &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv1(im1,jp,kp) & 
                &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv1(ip1,jp,kp)
#endif
#if N_DIM == 2 
           work(i,j,k,lb,1) = &
                &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv1(ip,jp,kp) & 
                &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv1(im1,jp,kp) & 
                &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv1(ip1,jp,kp) &
                &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * recv1(ip,jm1,kp) & 
                &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * recv1(ip,jp1,kp) & 
                &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * recv1(im1,jm1,kp) & 
                &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * recv1(ip1,jm1,kp) & 
                &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * recv1(im1,jp1,kp) & 
                &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * recv1(ip1,jp1,kp)
#endif
#if N_DIM == 3  
           work(i,j,k,lb,1) = &
                &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv1(ip,jp,kp) & 
                &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv1(im1,jp,kp) & 
                &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv1(ip1,jp,kp) &
                &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * recv1(ip,jm1,kp) & 
                &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * recv1(ip,jp1,kp) & 
                &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * recv1(im1,jm1,kp) & 
                &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * recv1(ip1,jm1,kp) & 
                &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * recv1(im1,jp1,kp) & 
                &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * recv1(ip1,jp1,kp) &
                &   + PC(0,ico) * PC(0,jco) * PC(-1,kco)* recv1(ip,jp,km1) & 
                &   + PC(-1,ico)* PC(0,jco) * PC(-1,kco)* recv1(im1,jp,km1) & 
                &   + PC(1,ico) * PC(0,jco) * PC(-1,kco)* recv1(ip1,jp,km1) & 
                &   + PC(0,ico) * PC(-1,jco)* PC(-1,kco)* recv1(ip,jm1,km1) & 
                &   + PC(0,ico) * PC(1,jco) * PC(-1,kco)* recv1(ip,jp1,km1) & 
                &   + PC(-1,ico)* PC(-1,jco)* PC(-1,kco)* recv1(im1,jm1,km1) & 
                &   + PC(1,ico) * PC(-1,jco)* PC(-1,kco)* recv1(ip1,jm1,km1) & 
                &   + PC(-1,ico)* PC(1,jco) * PC(-1,kco)* recv1(im1,jp1,km1) & 
                &   + PC(1,ico) * PC(1,jco) * PC(-1,kco)* recv1(ip1,jp1,km1) & 
                &   + PC(0,ico) * PC(0,jco) * PC(1,kco) * recv1(ip,jp,kp1) & 
                &   + PC(-1,ico)* PC(0,jco) * PC(1,kco) * recv1(im1,jp,kp1) & 
                &   + PC(1,ico) * PC(0,jco) * PC(1,kco) * recv1(ip1,jp,kp1) & 
                &   + PC(0,ico) * PC(-1,jco)* PC(1,kco) * recv1(ip,jm1,kp1) & 
                &   + PC(0,ico) * PC(1,jco) * PC(1,kco) * recv1(ip,jp1,kp1) & 
                &   + PC(-1,ico)* PC(-1,jco)* PC(1,kco) * recv1(im1,jm1,kp1) & 
                &   + PC(1,ico) * PC(-1,jco)* PC(1,kco) * recv1(ip1,jm1,kp1) & 
                &   + PC(-1,ico)* PC(1,jco) * PC(1,kco) * recv1(im1,jp1,kp1) & 
                &   + PC(1,ico) * PC(1,jco) * PC(1,kco) * recv1(ip1,jp1,kp1)
#endif
           
        enddo
     enddo
  enddo
  
!===============================================================================
  
  if (gr_monotone) then
     
     ! Monotonicity checking:  force interpolants to be flat if fine-zone averages
     ! fall outside the range of values of coarse-zone averages used to constrain
     ! them.
     
     icmin = ((ia-1)/2) + 1 + (nguard/2) + ioff
     icmax = ((ib-1)/2) + 1 + (nguard/2) + ioff
     jcmin = ((ja-1)/2) + 1 + (nguard/2)*k2d + joff
     jcmax = ((jb-1)/2) + 1 + (nguard/2)*k2d + joff
     kcmin = ((ka-1)/2) + 1 + (nguard/2)*k3d + koff
     kcmax = ((kb-1)/2) + 1 + (nguard/2)*k3d + koff
     
     k = ka - 2
     do kc = kcmin, kcmax
        k = k + 2
        j = ja - 2
        do jc = jcmin, jcmax
           j = j + 2
           i = ia - 2
           do ic = icmin, icmax
              i = i + 2
              
              cmin = minval(recv1(ic-1:ic+1,jc-k2d:jc+k2d, & 
                   &                            kc-k3d:kc+k3d))
              cmax = maxval(recv1(ic-1:ic+1,jc-k2d:jc+k2d, & 
                   &                            kc-k3d:kc+k3d))
              fmin = minval(work(i:i+1,j:j+k2d,k:k+k3d,lb,1))
              fmax = maxval(work(i:i+1,j:j+k2d,k:k+k3d,lb,1))
              
              if ( (fmin < cmin) .or. (fmax > cmax) ) then
                 do kp = k, k+k3d
                    do jp = j, j+k2d
                       do ip = i, i+1
                          work(ip,jp,kp,lb,1) = recv1(ic,jc,kc)
                       enddo
                    enddo
                 enddo
              endif
              
           enddo
        enddo
     enddo
     
  endif
  
  !===============================================================================
  
  return
end subroutine amr_prolong_gen_work_fun
