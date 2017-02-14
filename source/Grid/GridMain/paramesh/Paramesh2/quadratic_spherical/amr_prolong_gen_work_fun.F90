#define MONOTONIC


!!****if* source/mesh/amr/paramesh2.0/quadratic_cylindrical/amr_prolong_gen_work_fun
!!
!! NAME
!!
!!  amr_prolong_gen_work_fun
!!
!!
!! SYNOPSIS
!!
!!  amr_prolong_gen_work_fun(ia, ib, ja, jb, ka, kb, isg, &
!!                          ioff, joff, koff, mype)
!!
!!
!!  amr_prolong_gen_work_fun(integer, integer, integer, integer, &
!!                          integer, integer, integer, &
!!                          integer, integer, integer, integer)
!!
!!
!! DESCRIPTION
!!
!!
!!  This routine takes data from the array recv1, originally extracted 
!!  from the solution array work, and performs a prolongation operation 
!!  on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!  The data in recv1 is from a parent block (and may have come from off
!!  processor) and the result of the prolongation operation is written 
!!  directly into block isg, which is one of its children (and on the 
!!  local processor).  The position of the child within the parent block 
!!  is specified by the ioff, joff, and koff arguments (given in terms 
!!  of parent block zones, excluding guardcells).
!!
!!  This particular prolongation is conservative, cell-averaged, quadratic
!!  interpolation.  Geometrical effects are accounted for 2-d cylindrical
!!  coords.  
!!
!!  In 1-d prolongation, the following is done:
!!
!!
!!   +---------+---------+---------+     (coarse level/parent)
!!
!!      im1        ip        ip1         coarse zone indices
!!
!!
!!   +----+----+xxxx+----+----+----+     (children)
!!
!!
!!  We want to find the value of each zone-average variable in the children
!!  of zone ip (we denote one in the figure above with 'xxxx').  We use a
!!  (bi-)quadratic zone-average preserving polynomial reconstruction of the
!!  parent zones around our target fine zone (ip is the parent, ip1 and 1m1
!!  are the zones to the right and left respectively).  This is third order
!!  accurate.
!!
!!  For 2-d cylindrical coordinates, the polynomial is of the form
!!
!!   f(r,z) = a_1 r**2 + a_2 r + a_3 + a_4 z**2 + a_5 z
!!
!!  (note, we neglect the cross-terms, because this was a whole lot 
!!   easier to derive on a piece of paper, and appears to be sufficient
!!   enough).
!!
!!  5 parent zones are required to do the prolongation:
!!
!!
!!                    +-------+
!!                    |       |
!!                    |       |
!!                    |       |
!!            +-------+---+---+-------+
!!            |       |   |   |       |
!!            |       +---+---+       |
!!            |       |   |   |       |
!!            +-------+---+---+-------+
!!                    |       |
!!                    |       |
!!                    |       |
!!                    +-------+
!!
!!
!!  subject to the constraints that
!!
!!  1/V_i,j int_{z_{j-1/2}}^{z_{j+1/2}} dz 
!!          int_{r_{i-1/2}}^{r_{i+1/2}} r dr f(r,z) = <f>_{i,j}
!!
!!
!!
!!  There is no guarantee that this is monotonic, so we check to see
!!  whether the newly prolonged data falls outside of the bounds set
!!  by the parents, and if so, we replace the data with a direct copy
!!  from the immediate parent (for all children of that parent)
!!
!!  Furthermore, when you use a quadratic reconstruction polynomial 
!!  for your prolongation, or use a linear one in a non-Cartesian 
!!  coordinate system, there is no guarantee that you will preserve 
!!  the sum of the abundances equals 1 in the children.
!!
!!  if they do equal 1 to some tolerance by chance, there may still be
!!  oscillations.  Left to itself, the monotonicity constraint treats
!!  each variable separately.  Therefore, it may force one species to be 
!!  monotonic by just doing a direct copy from the parent, but not the 
!!  other species, breaking the sum of the abundances equal to 1.
!!
!!  so now, we test to see if the child obeys the sum of abundances = 1, 
!!  and if not, I replace all of them with the parent values.  I don't 
!!  think osciallations in the abundances is an issue, but if so, I 
!!  can also force monotonicity on all of them if one abundance fails 
!!  the test.  The trick is that they must be treated the same.
!!
!!
!! ARGUMENTS
!!
!!  ia, ib        the range of *child* zones to fill via prolongation. 
!!  ja, jb        (sometimes we are only doing the guardcells).  This is set 
!!  ka, kb        by amr_prolong_work_fun, the wrapper for the generic 
!!                prolongation.
!!
!!  isg           the block number of the child block whose data we are 
!!                filling.
!!
!!  ioff          the offsets of the child block into the parent (in terms of 
!!  joff          parent zones, excluding guardcells).
!!  koff
!!
!!  mype          the local processor number (why is this here?)
!!
!!
!! NOTES
!!
!!  To check for conservation, define CONSERVATION_CHECK.  This will 
!!  restrict the newly created child data up to the parent, and 
!!  compute the relative error and yell at you if it exceeds some
!!  tolerance.  The errors should be solely due to roundoff error.
!!
!!
!!***

subroutine amr_prolong_gen_work_fun(ia, ib, ja, jb, ka, kb, isg, & 
     ioff, joff, koff, mype)

!  use dBase, ONLY : nvar, ndim

  use workspace, ONLY : work, recv1, nguard_work
  use Driver_interface, ONLY : Driver_abortFlash

  use Grid_data, ONLY : gr_geometry, gr_oneBlock

  implicit none
#include "Flash.h"
#include "constants.h"

  integer :: ia, ib, ja, jb, ka, kb

  integer :: isg
  integer :: ioff, joff, koff

  integer :: mype

  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1

  integer :: i, j, k

  integer :: ico, jco, kco

  real :: a1, a2, a3, a4, a5

  real :: D3x, D3xp1, D4x, D4xp1
  real :: D3x_child, D4x_child
  real :: D3_D2, D3_D2p1, D3_D2m1, D4_D2, D4_D2p1, D4_D2m1, D4_D3, D4_D3p1
  real :: D4_D3_child, D4_D2_child, D3_D2_child
  real :: rratio, rratiop1, rratiom1, rratio_child

  real :: factor

  real :: cvol(4), pvol

  real :: error

  real, PARAMETER :: one_third = 1.0/3.0, &
                     two_thirds = 2.0/3.0, &
                     four_thirds = 4.0/3.0, &
                     one_twelvth = 1.0/12.0


! compute the maximum length of a vector in each coordinate direction 
! (including guardcells)
!!!  real, dimension(1:2*nguard+NXB) :: x, xl, xr
!!!  real, dimension(1:2*nguard*K2D+NYB) :: y, yl, yr 
!!!  real, dimension(1:2*nguard*K3D+NZB) :: z, zl, zr

  real :: xlc, xcc, xrc, ylc, ycc, yrc, zlc, zcc, zrc
  real :: xlcp1, xcp1, xrcp1, xlcm1, xcm1, xrcm1
  real :: ylcp1, ycp1, yrcp1, ylcm1, ycm1, yrcm1
  real :: zlcp1, zcp1, zrcp1, zlcm1, zcm1, zrcm1

  real :: dxc, dyc, dzc, dy_child

  integer :: ipmin, ipmax, jpmin, jpmax, kpmin, kpmax
  integer :: ichild1, jchild1, kchild1

  real :: parent_min, parent_max, child_min, child_max

  integer :: ii, jj, kk



! start by getting the coordinates of the child block that we are prolonging to
! -- remember we cannot get the parent coordinate information as easily, as it
! may live on another processor, therefore, we need to derive it from the child
! information and the fact that blocks differ by 2x in resolution between 
! levels of refinement.
!
! also compute the parent zone width (dxc, dyc, dzc) from the fine zones, 
! assuming a factor of 2 between levels

!!!  call dBaseGetCoords(izn,  iZcoord, isg, z)
!!!  call dBaseGetCoords(iznl, iZcoord, isg, zl)
!!!  call dBaseGetCoords(iznr, iZcoord, isg, zr)
#define ZCOORD(I) (gr_oneBlock(isg)%thirdAxisCoords(CENTER,I))
#define ZCOORD_LEFT(I) (gr_oneBlock(isg)%thirdAxisCoords(LEFT_EDGE,I))
#define ZCOORD_RIGHT(I) (gr_oneBlock(isg)%thirdAxisCoords(RIGHT_EDGE,I))

#if N_DIM == 3
  dzc = ZCOORD(3) - ZCOORD(1)
#else
  dzc = 0.0
#endif

!!!  call dBaseGetCoords(izn,  iYcoord, isg, y)
!!!  call dBaseGetCoords(iznl, iYcoord, isg, yl)
!!!  call dBaseGetCoords(iznr, iYcoord, isg, yr)
#define YCOORD(I) (gr_oneBlock(isg)%secondAxisCoords(CENTER,I))
#define YCOORD_LEFT(I) (gr_oneBlock(isg)%secondAxisCoords(LEFT_EDGE,I))
#define YCOORD_RIGHT(I) (gr_oneBlock(isg)%secondAxisCoords(RIGHT_EDGE,I))

#if N_DIM >= 2
  dyc = YCOORD(3) - YCOORD(1)
#else
  dyc = 0.0
#endif

!!!  call dBaseGetCoords(izn,  iXcoord, isg, x)
!!!  call dBaseGetCoords(iznl, iXcoord, isg, xl)
!!!  call dBaseGetCoords(iznr, iXcoord, isg, xr)
#define XCOORD(I) (gr_oneBlock(isg)%firstAxisCoords(CENTER,I))
#define XCOORD_LEFT(I) (gr_oneBlock(isg)%firstAxisCoords(LEFT_EDGE,I))
#define XCOORD_RIGHT(I) (gr_oneBlock(isg)%firstAxisCoords(RIGHT_EDGE,I))
  dxc = XCOORD(3) - XCOORD(1)

! loop over zones in the target (fine/child) block, and prolong the parent
! data to each child zone.  
!
! (i,j,k) are the indices of the fine zones.  
!
! (ip,jp,kp) are the indices of the parent zone enclosing each fine zone.  
!
! ip1 and im1 refer to parent zones offset by one to the right and left, 
! respectively; likewise for jp1, jm1, kp1, km1.  
!
! ico, jco, and kco indicate the location of the child in the parent zone in
! each coordinate direction.

! .         .         .         .         .
! +---------+---------+---------+---------+
! .    1    .    2    .    3    .    4    .     parent zone index 
! .         .         .         .         .
! +----+----+----+----+----+----+----+----+
! . 1    2  . 3    4  . 5    6  . 7    8  .     fine zone index
! .         .         .         .         .
!
!
!  here, fine zones 1 and 2 share the same parent/coarse zone 1
!
!  ico = 1 for fine zone 1  (assuming that this is the z-direction)
!  ico = 0 for fine zone 2 
!
!  we know the fine zone coordinates, XCOORD_LEFT, x, XCOORD_RIGHT, from the dBase,
!  and we use them to get the parent zone coordinates.  
!
!  The parent zone coordinates for parent zone 1 are:
!
!     xlc = XCOORD_LEFT(1)
!     xrc = XCOORD_RIGHT(2)
!     xcc = 0.5*(xlc + xrc)
!

  do k = ka, kb

     kp  = ((k-1)/2) + 1 + (nguard_work/2)*K3D + koff
     kp1 = kp + K3D
     km1 = kp - K3D

! kco tells us where the child falls in the parent zone.  kco = 1 means it 
! shares the lower z interface, kco = 0 means it shares the upper z interface 
! of the parent zone
     kco = mod(k, 2)

! find the z-coords of the parent cell containing the present child
     zlc = ZCOORD_LEFT(k-(1-kco)*K3D)
     zrc = ZCOORD_RIGHT(k+kco*K3D)
     zcc = 0.5*(zlc + zrc)

! we also need the adjacent cell coords.  XXcp1 is the k+1 parent zone, 
! XXcm1 is the k-1 parent zone.
     zlcp1 = zrc
     zrcp1 = zlcp1 + dzc
     zcp1 = 0.5*(zlcp1 + zrcp1)

     zrcm1 = zlc
     zlcm1 = zrcm1 - dzc
     zcm1 = 0.5*(zlcm1 + zrcm1) 

     do j = ja, jb

        jp  = ((j-1)/2) + 1 + (nguard_work/2)*K2D + joff
        jp1 = jp + K2D
        jm1 = jp - K2D

        jco = mod(j, 2)

! find the y-coords of the parent cell containing the present child
        ylc = YCOORD_LEFT(j-(1-jco)*K2D)
        yrc = YCOORD_RIGHT(j+jco*K2D)
        ycc = 0.5*(ylc + yrc)

! we also need the adjacent cell coords.  XXcp1 is the j+1 parent zone, 
! XXcm1 is the j-1 parent zone.
        ylcp1 = yrc
        yrcp1 = ylcp1 + dyc
        ycp1 = 0.5*(ylcp1 + yrcp1)

        yrcm1 = ylc
        ylcm1 = yrcm1 - dyc
        ycm1 = 0.5*(ylcm1 + yrcm1) 

        do i = ia, ib

           ip  = ((i-1)/2) + 1 + (nguard_work/2) + ioff
           ip1 = ip + 1
           im1 = ip - 1

           ico = mod(i, 2)

! find the x-coords of the parent cell containing the present child
           xlc = XCOORD_LEFT(i-(1-ico))
           xrc = XCOORD_RIGHT(i+ico)
           xcc = 0.5*(xlc + xrc)

! we also need the adjacent cell coords.  XXp1 is the i+1 zone, XXm1 is the i-1
! zone.
           xlcp1 = xrc
           xrcp1 = xlcp1 + dxc
           xcp1 = 0.5*(xlcp1 + xrcp1)

           xrcm1 = xlc
           xlcm1 = xrcm1 - dxc
           xcm1 = 0.5*(xlcm1 + xrcm1) 


! do the right prolongation, depending on the mesh geometry
           select case (gr_geometry)

           case (CARTESIAN)

              call Driver_abortFlash &
                   ("ERROR: I was too lazy to write the cartesian prolongation")

           case (CYLINDRICAL)

! remember that in FLASH, x is the cylindrical radial coordinate, and y is 
! the cylindrical z coordinate.

! compute some common factors
!
! DNx   = xrc**N - xlc**N
! DNxp1 = xrcp1**N - xlcp1**N
! DNxm1 = xrcm1**N - xlcm1**N
!
! DN_DM = DNx/DMx (note, where possible, we simplify the resulting 
!                  expression)
!
              if (xrc /= 0.0) then
                 rratio = xlc/xrc
                 D3_D2 = xrc*(1.0 + rratio + rratio**2)/(1.0 + rratio)
              else
                 D3_D2 = (xrc**3 - xlc**3)/(xrc**2 - xlc**2)
              endif

              D4_D2 = (xrc**2 + xlc**2)


              if (xrcp1 /= 0.0) then
                 rratiop1 = xlcp1/xrcp1
                 D3_D2p1 = xrcp1*(1.0 + rratiop1 + rratiop1**2)/(1.0 + rratiop1)
              else
                 D3_D2p1 = (xrcp1**3 - xlcp1**3)/(xrcp1**2 - xlcp1**2)
              endif

              D4_D2p1 = (xrcp1**2 + xlcp1**2)


              if (xrcm1 /= 0.0) then
                 rratiom1 = xlcm1/xrcm1
                 D3_D2m1 = xrcm1*(1.0 + rratiom1 + rratiom1**2)/(1.0 + rratiom1)
              else
                 D3_D2m1 = (xrcm1**3 - xlcm1**3)/(xrcm1**2 - xlcm1**2)
              endif

              D4_D2m1 = (xrcm1**2 + xlcm1**2)
              

              if (XCOORD_RIGHT(i) /= 0.e0) then
                 rratio_child = XCOORD_LEFT(i)/XCOORD_RIGHT(i)
                 D3_D2_child = XCOORD_RIGHT(i)*(1.0 + rratio_child + rratio_child**2)/ &
                      (1.0 + rratio_child)
              else
                 D3_D2_child = (XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)/(XCOORD_RIGHT(i)**2 - XCOORD_LEFT(i)**2)
              endif

              D4_D2_child = (XCOORD_RIGHT(i)**2 + XCOORD_LEFT(i)**2)

              dy_child = YCOORD_RIGHT(j) - YCOORD_LEFT(j)


! the dy_child term has a sign in front of it depending on where the
! child falls in the parent in the z-direction.  jco = 1 means it is 
! at the lower z interface (and the sign is -).  jco = 0 means it is
! at the upper z interface (and the sign is +).
              factor = float(1 - 2*jco)

!              print *, 'i,j,k,isg = ', i,j,k,isg
!              print *, 'ip, jp, kp = ', ip, jp, kp
!              print *, 'D4_D2, D4_D2m1, D4_D2p1 = ', D4_D2, D4_D2m1, D4_D2p1
!              print *, 'D3_D2, D3_D2m1, D3_D2p1 = ', D3_D2, D3_D2m1, D3_D2p1
!              print *, 'D4_D2_child, D3_D2_child = ', D4_D2_child, D3_D2_child
!              print *, 'xlcm1, xlc, xlcp1 = ', xlcm1, xlc, xlcp1
!              print *, 'xrcm1, xrc, xrcp1 = ', xrcm1, xrc, xrcp1
!              print *, 'dyc, dy_child = ', dyc, dy_child

! compute the coefficients of the reconstruction polynomial.  For simplicity,
! we express some of them in terms of others.
              a4 = 0.5*(recv1(ip,jp1,kp) - 2.0*recv1(ip,jp,kp) + &
                        recv1(ip,jm1,kp))/dyc**2

              a5 = (recv1(ip,jp1,kp) - recv1(ip,jm1,kp))/ &
                   (2.0*dyc)


! if we are right at the left boundary, r=0, and some terms can blow up --
! treat these specially.
              if (D4_D2 /= D4_D2p1) then

                 a2 = &
                      1.5*( (D4_D2 - D4_D2p1)* &
                             (recv1(ip,jp,kp) - recv1(im1,jp,kp)) - &
                            (D4_D2 - D4_D2m1)* &
                             (recv1(ip,jp,kp) - recv1(ip1,jp,kp)) )/ &
                          ( (D4_D2 - D4_D2p1)*(D3_D2 - D3_D2m1) - &
                            (D4_D2 - D4_D2m1)*(D3_D2 - D3_D2p1))

                 a1 = (2.0/(D4_D2 - D4_D2p1))* &
                      (recv1(ip,jp,kp) - recv1(ip1,jp,kp) - &
                      two_thirds*a2*(D3_D2 - D3_D2p1))

              else

                 a2 = &
                      1.5*( (D4_D2 - D4_D2m1)* &
                             (recv1(ip,jp,kp) - recv1(ip1,jp,kp)) - &
                            (D4_D2 - D4_D2p1)* &
                             (recv1(ip,jp,kp) - recv1(im1,jp,kp)) )/ &
                          ( (D4_D2 - D4_D2m1)*(D3_D2 - D3_D2p1) - &
                            (D4_D2 - D4_D2p1)*(D3_D2 - D3_D2m1))

                 a1 = (2.0/(D4_D2 - D4_D2m1))* &
                      (recv1(ip,jp,kp) - recv1(im1,jp,kp) - &
                      two_thirds*a2*(D3_D2 - D3_D2m1))


              endif

              a3 = recv1(ip,jp,kp) - 0.5*a1*D4_D2 - &
                   two_thirds*a2*D3_D2 - one_twelvth*a4*dyc**2


! the interpolating polynomial is integrated over the child cell volume to
! get the prolonged cell-average value for that child.  Compute that now
              work(i,j,k,isg,1) =  &
                   0.5*a1*D4_D2_child + &
                   two_thirds*a2*D3_D2_child + &
                   a3 + &
                   one_third*a4*dy_child**2 + &
                   0.5*factor*a5*dy_child
                 
!                 print *, work(i,j,k,isg,1), recv1(ip,jp,kp)
!                 print *, 'a1, a2, a3 = ', a1, a2, a3
!                 print *, 'a4, a5 = ', a4, a5

!              print *, '---------------------------------------'
!              print *, ' '

           case (SPHERICAL)

              call Driver_abortFlash &
                   ("ERROR: I was too lazy to write the spherical prolongation")

           end select
            
        enddo
     enddo
  enddo



#ifdef MONOTONIC
!-----------------------------------------------------------------------------
! monotonicity hack -- with the geometry factors (r dr in cylindrical, 
! r**2 dr in spherical), we can oscillate -- this is a bad thing.  Amongst
! other things, the abundances can go negative -- that's really bad.  To 
! fix this up, loop over all the newly prolonged child zones and, if we've
! gone outside the bounds of the parent zone data, then do just a direct
! copy for the prolongation.  Also, if the abundances don't sum to 1,
! do a direct copy.
!-----------------------------------------------------------------------------

! start by figuring out the coarse cells are the parents of the children
! we just filled
  ipmin  = ((ia-1)/2) + 1 + (nguard_work/2) + ioff
  ipmax  = ((ib-1)/2) + 1 + (nguard_work/2) + ioff

  jpmin  = ((ja-1)/2) + 1 + (nguard_work/2)*K2D + joff
  jpmax  = ((jb-1)/2) + 1 + (nguard_work/2)*K2D + joff

  kpmin  = ((ka-1)/2) + 1 + (nguard_work/2)*K3D + koff
  kpmax  = ((kb-1)/2) + 1 + (nguard_work/2)*K3D + koff


! loop over the parent cells and make sure that the zone-average values
! of the children call between the limits set by the parents that 
! are contained in the prolongation stencil
  kchild1 = ka
  do k = kpmin, kpmax

     jchild1 = ja
     do j = jpmin, jpmax

        ichild1 = ia
        do i = ipmin, ipmax

! for this cylindrical prolongation stencil, 5 parents were involved
           parent_min = min(recv1(i,  j,  k), &
                            recv1(i+1,j,  k), &
                            recv1(i-1,j,  k), &
                            recv1(i,  j+1,k), &
                            recv1(i,  j-1,k))

           parent_max = max(recv1(i,  j,  k), &
                            recv1(i+1,j,  k), &
                            recv1(i-1,j,  k), &
                            recv1(i,  j+1,k), &
                            recv1(i,  j-1,k))

! now find the extrema of the children 
           child_min = min(work(ichild1,  jchild1,    kchild1,    isg,1),&
                           work(ichild1+1,jchild1,    kchild1,    isg,1),&
                           work(ichild1,  jchild1+K2D,kchild1,    isg,1),&
                           work(ichild1+1,jchild1+K2D,kchild1,    isg,1),&
                           work(ichild1,  jchild1,    kchild1+K3D,isg,1),&
                           work(ichild1+1,jchild1,    kchild1+K3D,isg,1),&
                           work(ichild1,  jchild1+K2D,kchild1+K3D,isg,1),&
                           work(ichild1+1,jchild1+K2D,kchild1+K3D,isg,1))

           child_max = max(work(ichild1,  jchild1,    kchild1,    isg,1),&
                           work(ichild1+1,jchild1,    kchild1,    isg,1),&
                           work(ichild1,  jchild1+K2D,kchild1,    isg,1),&
                           work(ichild1+1,jchild1+K2D,kchild1,    isg,1),&
                           work(ichild1,  jchild1,    kchild1+K3D,isg,1),&
                           work(ichild1+1,jchild1,    kchild1+K3D,isg,1),&
                           work(ichild1,  jchild1+K2D,kchild1+K3D,isg,1),&
                           work(ichild1+1,jchild1+K2D,kchild1+K3D,isg,1))

! check to see if the extrema of the children go out of the bounds of the 
! parents
           if (child_min < parent_min .OR. &
               child_max > parent_max) then

! reset the children
              work(ichild1,  jchild1,    kchild1,    isg,1) = recv1(i,j,k)
              work(ichild1+1,jchild1,    kchild1,    isg,1) = recv1(i,j,k)
              work(ichild1,  jchild1+K2D,kchild1,    isg,1) = recv1(i,j,k)
              work(ichild1+1,jchild1+K2D,kchild1,    isg,1) = recv1(i,j,k)
              work(ichild1,  jchild1,    kchild1+K3D,isg,1) = recv1(i,j,k)
              work(ichild1+1,jchild1,    kchild1+K3D,isg,1) = recv1(i,j,k)
              work(ichild1,  jchild1+K2D,kchild1+K3D,isg,1) = recv1(i,j,k)
              work(ichild1+1,jchild1+K2D,kchild1+K3D,isg,1) = recv1(i,j,k)
                 
           endif
                             
           ichild1 = ichild1 + 2
        enddo
     
        jchild1 = jchild1 + 2*K2D
     enddo
  
     kchild1 = kchild1 + 2*K3D
  enddo

#endif


#ifdef CONSERVATION_CHECK
!-----------------------------------------------------------------------------
! check if we are conservative by performing the restriction operation
! on the children and comaring to the parent
!-----------------------------------------------------------------------------

! start by figuring out the coarse cells are the parents of the children
! we just filled
  ipmin  = ((ia-1)/2) + 1 + (nguard_work/2) + ioff
  ipmax  = ((ib-1)/2) + 1 + (nguard_work/2) + ioff

  jpmin  = ((ja-1)/2) + 1 + (nguard_work/2)*K2D + joff
  jpmax  = ((jb-1)/2) + 1 + (nguard_work/2)*K2D + joff

  kpmin  = ((ka-1)/2) + 1 + (nguard_work/2)*K3D + koff
  kpmax  = ((kb-1)/2) + 1 + (nguard_work/2)*K3D + koff


! loop over the parent cells and make sure that the zone-average values
! of the children call between the limits set by the parents that 
! are contained in the prolongation stencil
  kchild1 = ka
  do k = kpmin, kpmax

     jchild1 = ja
     do j = jpmin, jpmax

        ichild1 = ia
        do i = ipmin, ipmax

! compute the volume of the children
           select case (gr_geometry)

           case (CARTESIAN)

              call Driver_abortFlash( &
                   'ERROR: cartesian not implemented in amr_prolong_gen_work_fun.F90')

           case (CYLINDRICAL)
              cvol(1) = (XCOORD_RIGHT(ichild1) + XCOORD_LEFT(ichild1))* &
                        (XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                        (YCOORD_RIGHT(jchild1) - YCOORD_LEFT(jchild1))

              cvol(2) = (XCOORD_RIGHT(ichild1+1) + XCOORD_LEFT(ichild1+1))* &
                        (XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                        (YCOORD_RIGHT(jchild1) - YCOORD_LEFT(jchild1))

              cvol(3) = (XCOORD_RIGHT(ichild1) + XCOORD_LEFT(ichild1))* &
                        (XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                        (YCOORD_RIGHT(jchild1+1) - YCOORD_LEFT(jchild1+1))

              cvol(4) = (XCOORD_RIGHT(ichild1+1) + XCOORD_LEFT(ichild1+1))* &
                        (XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                        (YCOORD_RIGHT(jchild1+1) - YCOORD_LEFT(jchild1+1))

              pvol = cvol(1) + cvol(2) + cvol(3) + cvol(4)
              
              
           case (SPHERICAL)

              call Driver_abortFlash( &
                   'ERROR: spherical not implemented in amr_prolong_gen_work_fun.F90')
                 
           end select


           select case (gr_geometry)

           case (CARTESIAN)
              
              call Driver_abortFlash( &
                   'ERROR: cartesian not implemented in amr_prolong_gen_work_fun.F90')
              
           case (CYLINDRICAL)

              error = (recv1(i,j,k) - &
                   (cvol(1)*work(ichild1,  jchild1,  kchild1,isg,1) + &
                   cvol(2)*work(ichild1+1,jchild1,  kchild1,isg,1) + &
                   cvol(3)*work(ichild1,  jchild1+1,kchild1,isg,1) + &
                   cvol(4)*work(ichild1+1,jchild1+1,kchild1,isg,1))/pvol)/ &
                   recv1(i,j,k)

           case (SPHERICAL)

              call Driver_abortFlash( &
                   'ERROR: spherical not implemented in amr_prolong_gen_work_fun.F90')


           end select

           if (error > 1.e-10 .AND. recv1(i,j,k) /= 0.0) then
              print *, 'error = ', error, recv1(i,j,k)
              print *, a1, a2, a3, a4, a5
           endif
              
           ichild1 = ichild1 + 2
        enddo

        jchild1 = jchild1 + 2*K2D
     enddo

     kchild1 = kchild1 + 2*K3D
  enddo
  
#endif

        
  return
end subroutine amr_prolong_gen_work_fun
