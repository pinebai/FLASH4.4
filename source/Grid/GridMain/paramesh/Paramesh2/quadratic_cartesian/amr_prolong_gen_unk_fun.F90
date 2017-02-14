!*******************************************************************************
!
! Routine:      amr_prolong_gen_unk_fun
!
! Description:
!
! This routine takes data from the array recv, originally extracted 
! from the solution array unk, and performs a prolongation operation 
! on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
! The data in recv is from a parent block and the
! result of the prolongation operation is written directly into block
! isg, which is one of its children. The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! This particular prolongation is conservative, cell-averaged, triquadratic
! interpolation. It can be used for blocks with an even or odd number of grid
! cells.
!
! For 1-d prolongation, the following is done:
!
!
!   +---------+---------+---------+     (coarse level/parent)
!
!      im1        ip        ip1         coarse zone indices
!
!
!   +----+----+xxxx+----+----+----+     (children)
!
!
! we want to find the value of each zone-average variable in the children
! of zone ip (we denote one in the figure above with 'xxxx').  We use a
! second order zone-average preserving polynomial reconstruction of the
! coarse zones around our target fine zone (ip is the parent, ip1 and 1m1
! are the zones to the right and left respectively).
!
! The interpolating polynomial is:
!
!  f(z) = 1/(2 dz**2) (<f>_1 - 2 <f>_0 + <f>_-1) z**2 +
!         1/(2 dz) (<f>_1 - <f>_-1) z +
!         1/24 (-<f>_1 + 26 <f>_0 - <f>_-1)
!
! (see Laney, Computational Gasdynamics, for example).  This satistifies the
! constraint that is properly reproduces the zone average in each of the 
! coarse zones
!
! We now integrate this reconstruction over one of the fine zones that is a
! child of the ip coarse zone, and divide by the size of the fine zone to
! yield the zone average value of variable f in the fine zone:
!
!  <f>_fine, L = <f>_0 + (1/8) <f>_1 - (1/8) <f>_-1
!
! The procedure is similar for higher dimensions.
!
!

subroutine amr_prolong_gen_unk_fun & 
     &       (recv, ia, ib, ja, jb, ka, kb, isg, & 
     &        ioff, joff, koff, mype, isrc)
  
!===============================================================================
  
  use Grid_data, ONLY: gr_monotone
  use physicaldata
  use tree
  implicit none
  include 'mpif.h'
  


  real,intent(IN)    :: recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
  integer,intent(IN) :: ia, ib, ja, jb, ka, kb
  integer,intent(IN) :: isg, ioff, joff, koff, mype
  integer,intent(IN),OPTIONAL :: isrc
  
  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1
  integer :: ic, jc, kc, icmin, icmax, jcmin, jcmax
  integer :: kcmin, kcmax
  integer :: ivar, i, j, k, ico, jco, kco
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
     
     kp  = ((k-1)/2) + 1 + (nguard/2)*k3d + koff
     kp1 = kp + k3d
     km1 = kp - k3d
     kco = mod(k, 2)
     
     do j = ja, jb
        
        jp  = ((j-1)/2) + 1 + (nguard/2)*k2d + joff
        jp1 = jp + k2d
        jm1 = jp - k2d
        jco = mod(j, 2)
        
        do i = ia, ib
           
           ip  = ((i-1)/2) + 1 + (nguard/2) + ioff
           ip1 = ip + 1
           im1 = ip - 1
           ico = mod(i, 2)
           
           do ivar = 1, nvar
              
#if N_DIM == 1
              unk(ivar,i,j,k,isg) = &
                   &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip,jp,kp) & 
                   &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv(ivar,im1,jp,kp) & 
                   &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip1,jp,kp)
#endif
#if N_DIM == 2 
              unk(ivar,i,j,k,isg) = &
                   &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip,jp,kp) & 
                   &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv(ivar,im1,jp,kp) & 
                   &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip1,jp,kp) &
                   &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * recv(ivar,ip,jm1,kp) & 
                   &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * recv(ivar,ip,jp1,kp) & 
                   &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * recv(ivar,im1,jm1,kp) & 
                   &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * recv(ivar,ip1,jm1,kp) & 
                   &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * recv(ivar,im1,jp1,kp) & 
                   &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * recv(ivar,ip1,jp1,kp)
#endif
#if N_DIM == 3  
              unk(ivar,i,j,k,isg) = &
                   &     PC(0,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip,jp,kp) & 
                   &   + PC(-1,ico)* PC(0,jco) * PC(0,kco) * recv(ivar,im1,jp,kp) & 
                   &   + PC(1,ico) * PC(0,jco) * PC(0,kco) * recv(ivar,ip1,jp,kp) &
                   &   + PC(0,ico) * PC(-1,jco)* PC(0,kco) * recv(ivar,ip,jm1,kp) & 
                   &   + PC(0,ico) * PC(1,jco) * PC(0,kco) * recv(ivar,ip,jp1,kp) & 
                   &   + PC(-1,ico)* PC(-1,jco)* PC(0,kco) * recv(ivar,im1,jm1,kp) & 
                   &   + PC(1,ico) * PC(-1,jco)* PC(0,kco) * recv(ivar,ip1,jm1,kp) & 
                   &   + PC(-1,ico)* PC(1,jco) * PC(0,kco) * recv(ivar,im1,jp1,kp) & 
                   &   + PC(1,ico) * PC(1,jco) * PC(0,kco) * recv(ivar,ip1,jp1,kp) &
                   &   + PC(0,ico) * PC(0,jco) * PC(-1,kco)* recv(ivar,ip,jp,km1) & 
                   &   + PC(-1,ico)* PC(0,jco) * PC(-1,kco)* recv(ivar,im1,jp,km1) & 
                   &   + PC(1,ico) * PC(0,jco) * PC(-1,kco)* recv(ivar,ip1,jp,km1) & 
                   &   + PC(0,ico) * PC(-1,jco)* PC(-1,kco)* recv(ivar,ip,jm1,km1) & 
                   &   + PC(0,ico) * PC(1,jco) * PC(-1,kco)* recv(ivar,ip,jp1,km1) & 
                   &   + PC(-1,ico)* PC(-1,jco)* PC(-1,kco)* recv(ivar,im1,jm1,km1) & 
                   &   + PC(1,ico) * PC(-1,jco)* PC(-1,kco)* recv(ivar,ip1,jm1,km1) & 
                   &   + PC(-1,ico)* PC(1,jco) * PC(-1,kco)* recv(ivar,im1,jp1,km1) & 
                   &   + PC(1,ico) * PC(1,jco) * PC(-1,kco)* recv(ivar,ip1,jp1,km1) & 
                   &   + PC(0,ico) * PC(0,jco) * PC(1,kco) * recv(ivar,ip,jp,kp1) & 
                   &   + PC(-1,ico)* PC(0,jco) * PC(1,kco) * recv(ivar,im1,jp,kp1) & 
                   &   + PC(1,ico) * PC(0,jco) * PC(1,kco) * recv(ivar,ip1,jp,kp1) & 
                   &   + PC(0,ico) * PC(-1,jco)* PC(1,kco) * recv(ivar,ip,jm1,kp1) & 
                   &   + PC(0,ico) * PC(1,jco) * PC(1,kco) * recv(ivar,ip,jp1,kp1) & 
                   &   + PC(-1,ico)* PC(-1,jco)* PC(1,kco) * recv(ivar,im1,jm1,kp1) & 
                   &   + PC(1,ico) * PC(-1,jco)* PC(1,kco) * recv(ivar,ip1,jm1,kp1) & 
                   &   + PC(-1,ico)* PC(1,jco) * PC(1,kco) * recv(ivar,im1,jp1,kp1) & 
                   &   + PC(1,ico) * PC(1,jco) * PC(1,kco) * recv(ivar,ip1,jp1,kp1)
#endif
              
           enddo
           
        enddo
     enddo
  enddo
  
!===============================================================================
  
  ! Monotonicity checking:  force interpolants to be flat if fine-zone averages
  ! fall outside the range of values of coarse-zone averages used to constrain
  ! them.

  if (gr_monotone) then
     
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
              
              do ivar = 1, nvar
                 
                 cmin = minval(recv(ivar,ic-1:ic+1,jc-k2d:jc+k2d, & 
                      &                             kc-k3d:kc+k3d))
                 cmax = maxval(recv(ivar,ic-1:ic+1,jc-k2d:jc+k2d, & 
                      &                             kc-k3d:kc+k3d))
                 fmin = minval(unk(ivar,i:i+1,j:j+k2d,k:k+k3d,isg))
                 fmax = maxval(unk(ivar,i:i+1,j:j+k2d,k:k+k3d,isg))

                 if ( (fmin < cmin) .or. (fmax > cmax) ) then
                    do kp = k, k+k3d
                       do jp = j, j+k2d
                          do ip = i, i+1
                             unk(ivar,ip,jp,kp,isg) = recv(ivar,ic,jc,kc)
                          enddo
                       enddo
                    enddo
                 endif
                 
              enddo
              
           enddo
        enddo
     enddo
     
  endif
  
  !===============================================================================
  return
end subroutine amr_prolong_gen_unk_fun
