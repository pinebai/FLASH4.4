!#define ABUNDANCE_CHECK_INTERIORS
#define MONOTONIC
!#define CONSERVATION_CHECK

!!****if* source/mesh/amr/paramesh2.0/quadratic_spherical/amr_prolong_gen_unk_fun
!!
!! NAME
!!
!!  amr_prolong_gen_unk_fun
!!
!!
!! SYNOPSIS
!!
!!  amr_prolong_gen_unk_fun(recv, ia, ib, ja, jb, ka, kb, isg, &
!!                          ioff, joff, koff, mype)
!!
!!
!!  amr_prolong_gen_unk_fun(real(), integer, integer, integer, integer, &
!!                          integer, integer, integer, &
!!                          integer, integer, integer, integer)
!!
!!
!! DESCRIPTION
!!
!!
!!  This routine takes data from the array recv, originally extracted 
!!  from the solution array unk, and performs a prolongation operation 
!!  on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
!!  The data in recv is from a parent block (and may have come from off
!!  processor) and the result of the prolongation operation is written 
!!  directly into block isg, which is one of its children (and on the 
!!  local processor).  The position of the child within the parent block 
!!  is specified by the ioff, joff, and koff arguments (given in terms 
!!  of parent block zones, excluding guardcells).
!!
!!  This particular prolongation is conservative, cell-averaged, quadratic
!!  interpolation.  Geometrical effects are accounted for 1-d spherical
!!  coords (r only), and 2-d axisymmetric spherical coords (r,theta)
!!
!!-----------------------------------------------------------------------------
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
!!  For 1-d spherical coordinates, the polynomial is of the form
!!
!!   f(r) = a_1 r**2 + a_2 r + a_3 
!!
!!  (note, we neglect the cross-terms, because this was a whole lot 
!!   easier to derive on a piece of paper, and appears to be sufficient
!!   enough).
!!
!!  3 parent zones are required to do the prolongation:
!!
!!  subject to the constraints that
!!
!!  4 pi/V_i int_{r_{i-1/2}}^{r_{i+1/2}} r**2 dr f(r,z) = <f>_{i,j}
!!
!!
!!-----------------------------------------------------------------------------
!!  In 2-d axisymmetric spherical coordinates (r, theta), the volume of 
!!  a zone is:
!!
!!    V_{i,j} = (2/3) pi (r_{i+1/2}^3 - r_{i-1/2}^3) *
!!                       (cos t_{j-1/2} - cos t{j+1/2})
!!
!!  (we use t for theta)
!!
!!  The reconstruction polynomial is 
!!
!!   f(r,t) = a_1 r^2 + a_2 r + a_3 + a_4 t^2 + a_5 t
!!
!!
!!  subject to the constraints that
!!
!!  2 pi / V_{i,j} \int_{t_{j-1/2}}^{t_{j+1/2}} sin t dt 
!!                 \int_{r_{i-1/2}}^{r_{i+1/2}} r^2 dr f(r,t) = <f>_{i,j}
!!
!!
!!  This leads to five equations for five unknowns (a_1 to a_5), which
!!  can be solved.
!!
!!-----------------------------------------------------------------------------
!!  In 3-d spherical coordinates (r, theta, phi), the volume of 
!!  a zone is:
!!
!!    V_{i,j} = (1/3) pi (r_{i+1/2}^3 - r_{i-1/2}^3) *
!!                       (cos t_{j-1/2} - cos t{j+1/2}) *
!!                       (p_{k+1/2} - p_{k-1/2})
!!
!!  (we use t for theta, and p for phi)
!!
!!  The reconstruction polynomial is 
!!
!!   f(r,t) = a_1 r^2 + a_2 r + a_3 + a_4 t^2 + a_5 t + a_6 p^2 + a_7 p
!!
!!
!!  subject to the constraints that
!!
!!  1 / V_{i,j} \int_{p_{j-1/2}}^{p_{j+1/2}} dp 
!!              \int_{t_{j-1/2}}^{t_{j+1/2}} sin t dt 
!!              \int_{r_{i-1/2}}^{r_{i+1/2}} r^2 dr f(r,t,p) = <f>_{i,j,k}
!!
!!
!!  This leads to seven equations for seven unknowns (a_1 to a_7), which
!!  can be solved.
!!
!!  We note, that since we do not include the cross-terms, a_1 through a_5 in
!!  the 3-d case are identical to those in the 2-d case.
!!
!!-----------------------------------------------------------------------------
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
!!  recv          the coarse/parent data that will be used to initialize the
!!                child block.  An entire block (including guardcells) of 
!!                coarse data is provided here
!!
!!  ia, ib        the range of *child* zones to fill via prolongation. 
!!  ja, jb        (sometimes we are only doing the guardcells).  This is set 
!!  ka, kb        by amr_prolong_unk_fun, the wrapper for the generic 
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
!!  The expressions here are at present unpublished.
!!
!!  M. Zingale 12-2002
!!
!!  Occasionally, interpolation may yield negative values for abundances.
!!  This usually happens inside unresolved density gradients. In this
!!  case the fine data is set to coarse data. The value triggering this
!!  mechanism is xnuc_neg (negative small epsilon).
!!
!!  T. Plewa, 09-2003
!!
!!***

subroutine amr_prolong_gen_unk_fun(recv, ia, ib, ja, jb, ka, kb, isg, & 
     ioff, joff, koff, mype, isrc)

!  use dBase, ONLY :iLo_gc, iHi_gc, jLo_gc, jHi_gc, kLo_gc, kHi_gc

  use physicaldata, ONLY: unk
  use Driver_interface, ONLY : Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBC
  use Grid_data, ONLY : gr_geometry, gr_oneBlock, gr_convertToConsvdForMeshCalls
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get 

  implicit none
#include "constants.h"
#include "Flash.h"

  real, intent(IN) :: recv(NUNK_VARS,GRID_ILO_GC:GRID_IHI_GC,GRID_JLO_GC:GRID_JHI_GC,GRID_KLO_GC:GRID_KHI_GC)

  integer,intent(IN) :: ia, ib, ja, jb, ka, kb

  integer,intent(IN) :: isg
  integer,intent(IN) :: ioff, joff, koff

  integer,intent(IN) :: mype
  integer,intent(IN),OPTIONAL :: isrc

  integer :: ip, ip1, im1, jp, jp1, jm1, kp, kp1, km1

  integer :: stencil
  integer, parameter :: LEFT = 1, NORMAL = 2, RIGHT = 3

  integer :: ivar
  integer :: i, j, k, n

  integer :: ico, jco, kco

  real :: a1, a2, a3, a4, a5, a6, a7

  real :: D4_D3, D4_D3p1, D4_D3m1
  real :: D5_D3, D5_D3p1, D5_D3m1

  real :: D4_D3_child, D5_D3_child
  real :: rratio, rratiop1, rratiom1, rratio_child

  real :: A_j, A_jp1, A_jm1, A_jchild
  real :: B_j, B_jp1, B_jm1, B_jchild

  real :: cos_t_jl, cos_t_jr, sin_t_jl, sin_t_jr
  real :: cos_t_jlp1, cos_t_jrp1, sin_t_jlp1, sin_t_jrp1
  real :: cos_t_jlm1, cos_t_jrm1, sin_t_jlm1, sin_t_jrm1
  real :: cos_t_jlchild, cos_t_jrchild, sin_t_jlchild, sin_t_jrchild

  real :: dcos_t_j, dcos_t_jp1, dcos_t_jm1, dcos_t_jchild
  real :: dsin_t_j, dsin_t_jp1, dsin_t_jm1, dsin_t_jchild

  real :: Dphi3_Dphi, Dphi3_Dphip1, Dphi3_Dphim1, Dphi3_Dphi_child


  real :: cvol(8), pvol

  integer, save :: inuc_begin, inuc_end

  real, save :: smallx

  real :: error, sum, scale, smallrho

  real, PARAMETER :: one_third = 1.0/3.0, &
                     two_thirds = 2.0/3.0, &
                     four_thirds = 4.0/3.0, &
                     five_thirds = 5.0/3.0, &
                     one_twelvth = 1.0/12.0

!!!  real, DIMENSION(:,:,:,:), POINTER :: unk

! compute the maximum length of a vector in each coordinate direction 
! (including guardcells)
!!!  real, dimension(1:2*NGUARD+NXB) :: x, xl, xr
!!!  real, dimension(1:2*NGUARD*K2D+NYB) :: y, yl, yr
  real, dimension(1:2*NGUARD*K2D+NYB) :: dy 
!!!  real, dimension(1:2*NGUARD*K3D+NZB) :: z, zl, zr
  real, dimension(MDIM) :: del
  integer, dimension(2,MDIM) :: bcs

  real :: ypmin, ypmax, dely
  logical :: ypmin_bc_reflect, ypmax_bc_reflect
  real, dimension(2*NGUARD*K2D+NYB) :: yp, ypl, ypr, dyp
 
  real :: xlc, xcc, xrc, ylc, ycc, yrc, zlc, zcc, zrc
  real :: xlcp1, xcp1, xrcp1, xlcm1, xcm1, xrcm1
  real :: ylcp1, ycp1, yrcp1, ylcm1, ycm1, yrcm1
  real :: zlcp1, zcp1, zrcp1, zlcm1, zcm1, zrcm1

  real :: dxc, dyc, dzc

  integer :: ipmin, ipmax, jpmin, jpmax, kpmin, kpmax
  integer :: ichild1, jchild1, kchild1

  real :: parent_min, parent_max, child_min, child_max

  integer :: ii, jj, kk
  logical :: monotonize_abundances
  logical :: monotonize(NUNK_VARS)

  real, parameter :: xnuc_neg = -1.e-15
  real            :: xnuc, xnuc_min

  ! storage for the variable names in FLASH
!!  character (len=4), save :: flashVars(NUNK_VARS)

  logical, save :: firstCall = .true.


  if (firstCall) then

     inuc_begin = SPECIES_BEGIN
     inuc_end   = SPECIES_END


     call RuntimeParameters_get('smallx', smallx)

! get the names of the fluids being followed -- useful for debugging
!!     call getvarlabels(flashVars)

     firstCall = .false.

  endif

!!!  unk => dBaseGetDataPtrSingleBlock(isg, GC)

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
!!!  call dBaseGetCoords(idz,  iYcoord, isg, dy)
  call Grid_getDeltas(isg,del)
  dy = del(JAXIS)
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


! we need to treat the angular coordinates specially if we have a reflecting
! boundary condition -- mainly, we want to use a different stencil in the
! reconstruction of the parent data when we fall on a physical boundary.

! Here we get the min and max y-coordinate for the parent block, by looking
! at the child block data, and check if we are reflecting.  We use the fact
! that the angular coordinate is always spaced uniformly.

  call Grid_getBlkBC(isg,bcs)
  if (joff == 0) then

! we are dealing with the child block that is inside parent zones
! NGUARD+1, NGUARD+2, ..., NGUARD+NYB/2 in the y-direction

! this means that YCOORD_LEFT(NGUARD+1) is the left boundary of the parent block,
! and YCOORD_LEFT(NGUARD+1)+2*NYB*dy(NGUARD+1) is the right boundary
     ypmin = YCOORD_LEFT(NGUARD*K2D+1)
     ypmax = YCOORD_LEFT(NGUARD*K2D+1) + 2*NYB*dy(NGUARD*K2D+1)

! we also need to know whether we are on a reflecting boundary.  Since
! we are dealing with the minimum portion of the parent block in the
! y-direction, we can only get this information for the minimum y edge
! (and that is all we will need it at)
     if (bcs(LOW,JAXIS) == REFLECTING) then
        ypmin_bc_reflect = .true.
     else
        ypmin_bc_reflect = .false.
     endif

     ypmax_bc_reflect = .false.


  elseif (joff == NYB/2) then

! we are dealing with the child block that is inside parent zones
! NGUARD+NYB/2+1, ... NGUARD+NYB in the y-direction

! this means that YCOORD_RIGHT(NGUARD+NYB) is the right boundary of the parent block,
! and YCOORD_RIGHT(NGUARD+NYB)-2*NYB*dy(NGUARD+NYB) is the left boundary
     ypmax = YCOORD_RIGHT(NGUARD*K2D+NYB)
     ypmin = YCOORD_RIGHT(NGUARD*K2D+NYB) - 2*NYB*dy(NGUARD*K2D+NYB)

! now figure out if we are reflecting -- we can only know this about the
! maximum y-boundary, because of where our child block falls inside the
! parent.
     if (bcs(HIGH,JAXIS) == REFLECTING) then
        ypmax_bc_reflect = .true.
     else
        ypmax_bc_reflect = .false.
     endif

     ypmin_bc_reflect = .false.


  else
     call Driver_abortFlash("ERROR: joff does not make sense in amr_prolong_gen_unk_fun")
  endif

  dely = (ypmax - ypmin)/float(NYB)
  
  do j = 1, 2*NGUARD*K2D+NYB

     ypl(j) = ypmin + float(j-1-NGUARD)*dely
     ypr(j) = ypmax + float(j-NGUARD)*dely

     yp(j) = 0.5*(ypl(j) + ypr(j))

     dyp(j) = dely

  enddo

  do k = ka, kb

     kp  = ((k-1)/2) + 1 + (NGUARD/2)*K3D + koff
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

        jp  = ((j-1)/2) + 1 + (NGUARD/2)*K2D + joff
        jp1 = jp + K2D
        jm1 = jp - K2D

        if (ypmin_bc_reflect .AND. jp == NGUARD+1) then
           stencil = LEFT
        elseif (ypmax_bc_reflect .AND. jp == NGUARD+NYB) then
           stencil = RIGHT
        else
           stencil = NORMAL
        endif

        jco = mod(j, 2)

! find the y-coords of the parent cell containing the present child
        ylc = ypl(jp)
        yrc = ypr(jp)
        ycc = 0.5*(ylc + yrc)

! we also need the adjacent cell coords.  XXcp1 is the j+1 parent zone, 
! XXcm1 is the j-1 parent zone.
        ylcp1 = ypl(jp1)
        yrcp1 = ypr(jp1)
        ycp1 = 0.5*(ylcp1 + yrcp1)

        ylcm1 = ypl(jm1)
        yrcm1 = ypr(jm1)
        ycm1 = 0.5*(ylcm1 + yrcm1) 

        do i = ia, ib

           ip  = ((i-1)/2) + 1 + (NGUARD/2) + ioff
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


!------------------------------------------------------------------------------
! do the right prolongation, depending on the mesh geometry
!------------------------------------------------------------------------------
           select case (gr_geometry)

           case (CARTESIAN)

              call Driver_abortFlash &
                   ("ERROR: these prolongation routines are for spherical coords only")

           case (CYLINDRICAL)


              call Driver_abortFlash &
                   ("ERROR: these prolongation routines are for spherical coords only")


           case (SPHERICAL)

#if N_DIM == 1
!------------------------------------------------------------------------------
! 1-d spherical coords
!------------------------------------------------------------------------------

! remember that in FLASH, x is the spherical radial coordinate

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

                 D5_D3 = xrc**2 * &
                  (1.0 + rratio + rratio**2 + rratio**3 + rratio**4)/ &
                  (1.0 + rratio + rratio**2)
                 
                 D4_D3 = xrc * &
                  (1.0 + rratio)*(1.0 + rratio**2)/ &
                  (1.0 + rratio + rratio**2)
                 
              else
                 D5_D3 = (xrc**5 - xlc**5)/(xrc**3 - xlc**3)
                 D4_D3 = (xrc**4 - xlc**4)/(xrc**3 - xlc**3)
              endif



              if (xrcp1 /= 0.0) then
                 rratiop1 = xlcp1/xrcp1

                 D5_D3p1 = xrcp1**2 * &
                  (1.0 + rratiop1 + rratiop1**2 + rratiop1**3 + rratiop1**4)/ &
                  (1.0 + rratiop1 + rratiop1**2)

                 D4_D3p1 = xrcp1 * &
                  (1.0 + rratiop1)*(1.0 + rratiop1**2)/ &
                  (1.0 + rratiop1 + rratiop1**2)
                                  
              else
                 D5_D3p1 = (xrcp1**5 - xlcp1**5)/(xrcp1**3 - xlcp1**3)
                 D4_D3p1 = (xrcp1**4 - xlcp1**4)/(xrcp1**3 - xlcp1**3)
              endif



              if (xrcm1 /= 0.0) then
                 rratiom1 = xlcm1/xrcm1

                 D5_D3m1 = xrcm1**2 * &
                  (1.0 + rratiom1 + rratiom1**2 + rratiom1**3 + rratiom1**4)/ &
                  (1.0 + rratiom1 + rratiom1**2)

                 D4_D3m1 = xrcm1 * &
                  (1.0 + rratiom1)*(1.0 + rratiom1**2)/ &
                  (1.0 + rratiom1 + rratiom1**2)
                                  
              else
                 D5_D3m1 = (xrcm1**5 - xlcm1**5)/(xrcm1**3 - xlcm1**3)
                 D4_D3m1 = (xrcm1**4 - xlcm1**4)/(xrcm1**3 - xlcm1**3)
              endif


              if (XCOORD_RIGHT(i) /= 0.e0) then
                 rratio_child = XCOORD_LEFT(i)/XCOORD_RIGHT(i)

                 D5_D3_child = XCOORD_RIGHT(i)**2 * &
                      (1.0 + rratio_child + rratio_child**2 + &
                       rratio_child**3 + rratio_child**4)/ &
                      (1.0 + rratio_child + rratio_child**2)

                 D4_D3_child = XCOORD_RIGHT(i) * &
                      (1.0 + rratio_child)*(1.0 + rratio_child**2)/ &
                      (1.0 + rratio_child + rratio_child**2)

              else
                 D5_D3_child = (XCOORD_RIGHT(i)**5 - XCOORD_LEFT(i)**5)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
                 D4_D3_child = (XCOORD_RIGHT(i)**4 - XCOORD_LEFT(i)**4)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
              endif


              do ivar = 1, NUNK_VARS

! compute the coefficients of the reconstruction polynomial.  For simplicity,
! we express some of them in terms of others.
                 a2 = four_thirds * &
                      ( (D5_D3 - D5_D3p1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp)) - &
                        (D5_D3 - D5_D3m1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp)) ) / &
                      ( (D4_D3 - D4_D3m1)*(D5_D3 - D5_D3p1) - &
                        (D4_D3 - D4_D3p1)*(D5_D3 - D5_D3m1))
                         

! if we are right at the left boundary, r=0, and some terms can blow up --
! treat these specially.
                 if (D5_D3 /= D5_D3p1) then

                    a1 = (five_thirds/(D5_D3 - D5_D3p1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3p1))

                 else

                    a1 = (five_thirds/(D5_D3 - D5_D3m1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3m1))

                 endif

                 a3 = recv(ivar,ip,jp,kp) - 0.6*a1*D5_D3 - 0.75*a2*D4_D3


! the interpolating polynomial is integrated over the child cell volume to
! get the prolonged cell-average value for that child.  Compute that now
                 unk(ivar,i,j,k,isg) = &
                      0.6*a1*D5_D3_child + 0.75*a2*D4_D3_child + a3
                 
              enddo

#elif N_DIM == 2

!------------------------------------------------------------------------------
! 2-d spherical coords
!------------------------------------------------------------------------------

! in 2-d spherical coordinates, x is the radial coordinate, and y is the 
! theta angle.


! compute some common factors.  
!
! A_j is the factor in front of a_4 when the reconstruction poly. is integrated
!  over the volume.
!
! B_j is the factor in front of a_5
!
! A_jp1, A_jm1 are the same factors for the j +/- 1 coarse zone
!

              cos_t_jl = cos(ylc)
              cos_t_jr = cos(yrc)

              sin_t_jl = sin(ylc)
              sin_t_jr = sin(yrc)

              dcos_t_j = cos_t_jr - cos_t_jl
              dsin_t_j = sin_t_jr - sin_t_jl

              A_j = -1.0/dcos_t_j * &
                   (yrc*(2.0*sin_t_jr - yrc*cos_t_jr) - &
                    ylc*(2.0*sin_t_jl - ylc*cos_t_jl) + &
                    2.0*dcos_t_j)

              B_j = -1.0/dcos_t_j * &
                   (dsin_t_j - yrc*cos_t_jr + ylc*cos_t_jl)



              cos_t_jlp1 = cos(ylcp1)
              cos_t_jrp1 = cos(yrcp1)

              sin_t_jlp1 = sin(ylcp1)
              sin_t_jrp1 = sin(yrcp1)

              dcos_t_jp1 = cos_t_jrp1 - cos_t_jlp1
              dsin_t_jp1 = sin_t_jrp1 - sin_t_jlp1

              A_jp1 = -1.0/dcos_t_jp1 * &
                   (yrcp1*(2.0*sin_t_jrp1 - yrcp1*cos_t_jrp1) - &
                    ylcp1*(2.0*sin_t_jlp1 - ylcp1*cos_t_jlp1) + &
                    2.0*dcos_t_jp1)

              B_jp1 = -1.0/dcos_t_jp1 * &
                   (dsin_t_jp1 - yrcp1*cos_t_jrp1 + ylcp1*cos_t_jlp1)



              cos_t_jlm1 = cos(ylcm1)
              cos_t_jrm1 = cos(yrcm1)

              sin_t_jlm1 = sin(ylcm1)
              sin_t_jrm1 = sin(yrcm1)

              dcos_t_jm1 = cos_t_jrm1 - cos_t_jlm1
              dsin_t_jm1 = sin_t_jrm1 - sin_t_jlm1

              A_jm1 = -1.0/dcos_t_jm1 * &
                   (yrcm1*(2.0*sin_t_jrm1 - yrcm1*cos_t_jrm1) - &
                    ylcm1*(2.0*sin_t_jlm1 - ylcm1*cos_t_jlm1) + &
                    2.0*dcos_t_jm1)

              B_jm1 = -1.0/dcos_t_jm1 * &
                   (dsin_t_jm1 - yrcm1*cos_t_jrm1 + ylcm1*cos_t_jlm1)



              cos_t_jlchild = cos(YCOORD_LEFT(j))
              cos_t_jrchild = cos(YCOORD_RIGHT(j))

              sin_t_jlchild = sin(YCOORD_LEFT(j))
              sin_t_jrchild = sin(YCOORD_RIGHT(j))

              dcos_t_jchild = cos_t_jrchild - cos_t_jlchild
              dsin_t_jchild = sin_t_jrchild - sin_t_jlchild

              A_jchild = -1.0/dcos_t_jchild * &
                   (YCOORD_RIGHT(j)*(2.0*sin_t_jrchild - YCOORD_RIGHT(j)*cos_t_jrchild) - &
                    YCOORD_LEFT(j)*(2.0*sin_t_jlchild - YCOORD_LEFT(j)*cos_t_jlchild) + &
                    2.0*dcos_t_jchild)

              B_jchild = -1.0/dcos_t_jchild * &
                   (dsin_t_jchild - YCOORD_RIGHT(j)*cos_t_jrchild + YCOORD_LEFT(j)*cos_t_jlchild)


! compute the radial factors.  These all have the form of 
!
! DN_DM = DNx/DMx (note, where possible, we simplify the resulting 
!                  expression)
! 
! where
!   DNx   = xrc**N - xlc**N
!   DNxp1 = xrcp1**N - xlcp1**N
!   DNxm1 = xrcm1**N - xlcm1**N
!

              if (xrc /= 0.0) then
                 rratio = xlc/xrc

                 D5_D3 = xrc**2 * &
                  (1.0 + rratio + rratio**2 + rratio**3 + rratio**4)/ &
                  (1.0 + rratio + rratio**2)
                 
                 D4_D3 = xrc * &
                  (1.0 + rratio)*(1.0 + rratio**2)/ &
                  (1.0 + rratio + rratio**2)
                 
              else
                 D5_D3 = (xrc**5 - xlc**5)/(xrc**3 - xlc**3)
                 D4_D3 = (xrc**4 - xlc**4)/(xrc**3 - xlc**3)
              endif



              if (xrcp1 /= 0.0) then
                 rratiop1 = xlcp1/xrcp1

                 D5_D3p1 = xrcp1**2 * &
                  (1.0 + rratiop1 + rratiop1**2 + rratiop1**3 + rratiop1**4)/ &
                  (1.0 + rratiop1 + rratiop1**2)

                 D4_D3p1 = xrcp1 * &
                  (1.0 + rratiop1)*(1.0 + rratiop1**2)/ &
                  (1.0 + rratiop1 + rratiop1**2)
                                  
              else
                 D5_D3p1 = (xrcp1**5 - xlcp1**5)/(xrcp1**3 - xlcp1**3)
                 D4_D3p1 = (xrcp1**4 - xlcp1**4)/(xrcp1**3 - xlcp1**3)
              endif



              if (xrcm1 /= 0.0) then
                 rratiom1 = xlcm1/xrcm1

                 D5_D3m1 = xrcm1**2 * &
                  (1.0 + rratiom1 + rratiom1**2 + rratiom1**3 + rratiom1**4)/ &
                  (1.0 + rratiom1 + rratiom1**2)

                 D4_D3m1 = xrcm1 * &
                  (1.0 + rratiom1)*(1.0 + rratiom1**2)/ &
                  (1.0 + rratiom1 + rratiom1**2)
                                  
              else
                 D5_D3m1 = (xrcm1**5 - xlcm1**5)/(xrcm1**3 - xlcm1**3)
                 D4_D3m1 = (xrcm1**4 - xlcm1**4)/(xrcm1**3 - xlcm1**3)
              endif


              if (XCOORD_RIGHT(i) /= 0.e0) then
                 rratio_child = XCOORD_LEFT(i)/XCOORD_RIGHT(i)

                 D5_D3_child = XCOORD_RIGHT(i)**2 * &
                      (1.0 + rratio_child + rratio_child**2 + &
                       rratio_child**3 + rratio_child**4)/ &
                      (1.0 + rratio_child + rratio_child**2)

                 D4_D3_child = XCOORD_RIGHT(i) * &
                      (1.0 + rratio_child)*(1.0 + rratio_child**2)/ &
                      (1.0 + rratio_child + rratio_child**2)

              else
                 D5_D3_child = (XCOORD_RIGHT(i)**5 - XCOORD_LEFT(i)**5)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
                 D4_D3_child = (XCOORD_RIGHT(i)**4 - XCOORD_LEFT(i)**4)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
              endif


! go through variable by variable and compute the coefficients, a_1
! through a_5.  These are found by integrating the reconstruction
! polynomial over each of the parent cells in our stencil and solving
! the resulting eqs. for the unknowns.


              do ivar = 1, NUNK_VARS

! compute the coefficients of the reconstruction polynomial.  For simplicity,
! we express some of them in terms of others.

! the a_1 and a_2 coefficients are the same as in the 1-d case, since we are
! not including cross terms
                 a2 = four_thirds * &
                      ( (D5_D3 - D5_D3p1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp)) - &
                        (D5_D3 - D5_D3m1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp)) ) / &
                      ( (D4_D3 - D4_D3m1)*(D5_D3 - D5_D3p1) - &
                        (D4_D3 - D4_D3p1)*(D5_D3 - D5_D3m1))
                         

! if we are right at the left boundary, r=0, and some terms can blow up --
! treat these specially.
                 if (D5_D3 /= D5_D3p1) then

                    a1 = (five_thirds/(D5_D3 - D5_D3p1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3p1))

                 else

                    a1 = (five_thirds/(D5_D3 - D5_D3m1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3m1))

                 endif


                 if (stencil == NORMAL) then
                    a5 = ( (A_jm1 - A_j)* &
                         (recv(ivar,ip,jp1,kp) - recv(ivar,ip,jp,kp)) - &
                         (A_jp1 - A_j)* &
                         (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp)) ) / &
                         ( (B_jp1 - B_j)*(A_jm1 - A_j) - &
                         (B_jm1 - B_j)*(A_jp1 - A_j))

                    a4 = (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp) - &
                         a5*(B_jm1 - B_j))/(A_jm1 - A_j)

                 else if (stencil == LEFT) then
                    a5 = (recv(ivar,ip,jp1,kp) - recv(ivar,ip,jp,kp)) / &
                         (B_jp1 - B_j)

                    a4 = 0.0

                 else if (stencil == RIGHT) then
                    a5 = (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp)) / &
                         (B_jm1 - B_j)

                    a4 = 0.0
                 endif

                 a3 = recv(ivar,ip,jp,kp) - &
                      0.6*a1*D5_D3 - &
                      0.75*a2*D4_D3 - &
                      a4*A_j -&
                      a5*B_j


! the interpolating polynomial is integrated over the child cell volume to
! get the prolonged cell-average value for that child.  Compute that now
                 unk(ivar,i,j,k,isg) = 0.6*a1*D5_D3_child + &
                                       0.75*a2*D4_D3_child + &
                                       a3 + &
                                       a4*A_jchild + &
                                       a5*B_jchild
                 
              enddo


#else

!------------------------------------------------------------------------------
! 3-d spherical coords
!------------------------------------------------------------------------------

! in 3-d spherical coordinates, x is the radial coordinate, y is the theta 
! angle, and z is the phi angle.


! compute some common factors.  
!
! A_j is the factor in front of a_4 when the reconstruction poly. is integrated
!  over the volume.
!
! B_j is the factor in front of a_5
!
! A_jp1, A_jm1 are the same factors for the j +/- 1 coarse zone
!

              cos_t_jl = cos(ylc)
              cos_t_jr = cos(yrc)

              sin_t_jl = sin(ylc)
              sin_t_jr = sin(yrc)

              dcos_t_j = cos_t_jr - cos_t_jl
              dsin_t_j = sin_t_jr - sin_t_jl

              A_j = -1.0/dcos_t_j * &
                   (yrc*(2.0*sin_t_jr - yrc*cos_t_jr) - &
                    ylc*(2.0*sin_t_jl - ylc*cos_t_jl) + &
                    2.0*dcos_t_j)

              B_j = -1.0/dcos_t_j * &
                   (dsin_t_j - yrc*cos_t_jr + ylc*cos_t_jl)



              cos_t_jlp1 = cos(ylcp1)
              cos_t_jrp1 = cos(yrcp1)

              sin_t_jlp1 = sin(ylcp1)
              sin_t_jrp1 = sin(yrcp1)

              dcos_t_jp1 = cos_t_jrp1 - cos_t_jlp1
              dsin_t_jp1 = sin_t_jrp1 - sin_t_jlp1

              A_jp1 = -1.0/dcos_t_jp1 * &
                   (yrcp1*(2.0*sin_t_jrp1 - yrcp1*cos_t_jrp1) - &
                    ylcp1*(2.0*sin_t_jlp1 - ylcp1*cos_t_jlp1) + &
                    2.0*dcos_t_jp1)

              B_jp1 = -1.0/dcos_t_jp1 * &
                   (dsin_t_jp1 - yrcp1*cos_t_jrp1 + ylcp1*cos_t_jlp1)



              cos_t_jlm1 = cos(ylcm1)
              cos_t_jrm1 = cos(yrcm1)

              sin_t_jlm1 = sin(ylcm1)
              sin_t_jrm1 = sin(yrcm1)

              dcos_t_jm1 = cos_t_jrm1 - cos_t_jlm1
              dsin_t_jm1 = sin_t_jrm1 - sin_t_jlm1

              A_jm1 = -1.0/dcos_t_jm1 * &
                   (yrcm1*(2.0*sin_t_jrm1 - yrcm1*cos_t_jrm1) - &
                    ylcm1*(2.0*sin_t_jlm1 - ylcm1*cos_t_jlm1) + &
                    2.0*dcos_t_jm1)

              B_jm1 = -1.0/dcos_t_jm1 * &
                   (dsin_t_jm1 - yrcm1*cos_t_jrm1 + ylcm1*cos_t_jlm1)



              cos_t_jlchild = cos(YCOORD_LEFT(j))
              cos_t_jrchild = cos(YCOORD_RIGHT(j))

              sin_t_jlchild = sin(YCOORD_LEFT(j))
              sin_t_jrchild = sin(YCOORD_RIGHT(j))

              dcos_t_jchild = cos_t_jrchild - cos_t_jlchild
              dsin_t_jchild = sin_t_jrchild - sin_t_jlchild

              A_jchild = -1.0/dcos_t_jchild * &
                   (YCOORD_RIGHT(j)*(2.0*sin_t_jrchild - YCOORD_RIGHT(j)*cos_t_jrchild) - &
                    YCOORD_LEFT(j)*(2.0*sin_t_jlchild - YCOORD_LEFT(j)*cos_t_jlchild) + &
                    2.0*dcos_t_jchild)

              B_jchild = -1.0/dcos_t_jchild * &
                   (dsin_t_jchild - YCOORD_RIGHT(j)*cos_t_jrchild + YCOORD_LEFT(j)*cos_t_jlchild)


! compute the radial factors.  These all have the form of 
!
! DN_DM = DNx/DMx (note, where possible, we simplify the resulting 
!                  expression)
! 
! where
!   DNx   = xrc**N - xlc**N
!   DNxp1 = xrcp1**N - xlcp1**N
!   DNxm1 = xrcm1**N - xlcm1**N
!

              if (xrc /= 0.0) then
                 rratio = xlc/xrc

                 D5_D3 = xrc**2 * &
                  (1.0 + rratio + rratio**2 + rratio**3 + rratio**4)/ &
                  (1.0 + rratio + rratio**2)
                 
                 D4_D3 = xrc * &
                  (1.0 + rratio)*(1.0 + rratio**2)/ &
                  (1.0 + rratio + rratio**2)
                 
              else
                 D5_D3 = (xrc**5 - xlc**5)/(xrc**3 - xlc**3)
                 D4_D3 = (xrc**4 - xlc**4)/(xrc**3 - xlc**3)
              endif



              if (xrcp1 /= 0.0) then
                 rratiop1 = xlcp1/xrcp1

                 D5_D3p1 = xrcp1**2 * &
                  (1.0 + rratiop1 + rratiop1**2 + rratiop1**3 + rratiop1**4)/ &
                  (1.0 + rratiop1 + rratiop1**2)

                 D4_D3p1 = xrcp1 * &
                  (1.0 + rratiop1)*(1.0 + rratiop1**2)/ &
                  (1.0 + rratiop1 + rratiop1**2)
                                  
              else
                 D5_D3p1 = (xrcp1**5 - xlcp1**5)/(xrcp1**3 - xlcp1**3)
                 D4_D3p1 = (xrcp1**4 - xlcp1**4)/(xrcp1**3 - xlcp1**3)
              endif



              if (xrcm1 /= 0.0) then
                 rratiom1 = xlcm1/xrcm1

                 D5_D3m1 = xrcm1**2 * &
                  (1.0 + rratiom1 + rratiom1**2 + rratiom1**3 + rratiom1**4)/ &
                  (1.0 + rratiom1 + rratiom1**2)

                 D4_D3m1 = xrcm1 * &
                  (1.0 + rratiom1)*(1.0 + rratiom1**2)/ &
                  (1.0 + rratiom1 + rratiom1**2)
                                  
              else
                 D5_D3m1 = (xrcm1**5 - xlcm1**5)/(xrcm1**3 - xlcm1**3)
                 D4_D3m1 = (xrcm1**4 - xlcm1**4)/(xrcm1**3 - xlcm1**3)
              endif


              if (XCOORD_RIGHT(i) /= 0.e0) then
                 rratio_child = XCOORD_LEFT(i)/XCOORD_RIGHT(i)

                 D5_D3_child = XCOORD_RIGHT(i)**2 * &
                      (1.0 + rratio_child + rratio_child**2 + &
                       rratio_child**3 + rratio_child**4)/ &
                      (1.0 + rratio_child + rratio_child**2)

                 D4_D3_child = XCOORD_RIGHT(i) * &
                      (1.0 + rratio_child)*(1.0 + rratio_child**2)/ &
                      (1.0 + rratio_child + rratio_child**2)

              else
                 D5_D3_child = (XCOORD_RIGHT(i)**5 - XCOORD_LEFT(i)**5)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
                 D4_D3_child = (XCOORD_RIGHT(i)**4 - XCOORD_LEFT(i)**4)/(XCOORD_RIGHT(i)**3 - XCOORD_LEFT(i)**3)
              endif

! Dphi3_Dphi = (phi_{k+1/2}^3 - phi_{k-1/2}^3)/(phi_{k+1/2} - phi_{k-1/2}).  This 
! appear in front of the a_6 term, and has k+1, k-1, and child versions
              Dphi3_Dphi = (zrc**2 + zrc*zlc + zlc**2)
              Dphi3_Dphip1 = (zrcp1**2 + zrcp1*zlcp1 + zlcp1**2)
              Dphi3_Dphim1 = (zrcm1**2 + zrcm1*zlcm1 + zlcm1**2)

              Dphi3_Dphi = (ZCOORD_RIGHT(k)**2 + ZCOORD_RIGHT(k)*ZCOORD_LEFT(k) + ZCOORD_LEFT(k)**2)


! go through variable by variable and compute the coefficients, a_1
! through a_7.  These are found by integrating the reconstruction
! polynomial over each of the parent cells in our stencil and solving
! the resulting eqs. for the unknowns.


              do ivar = 1, NUNK_VARS

! compute the coefficients of the reconstruction polynomial.  For simplicity,
! we express some of them in terms of others.

! the a_1 and a_2 coefficients are the same as in the 1-d case, since we are
! not including cross terms
                 a2 = four_thirds * &
                      ( (D5_D3 - D5_D3p1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp)) - &
                        (D5_D3 - D5_D3m1)* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp)) ) / &
                      ( (D4_D3 - D4_D3m1)*(D5_D3 - D5_D3p1) - &
                        (D4_D3 - D4_D3p1)*(D5_D3 - D5_D3m1))
                         

! if we are right at the left boundary, r=0, and some terms can blow up --
! treat these specially.
                 if (D5_D3 /= D5_D3p1) then

                    a1 = (five_thirds/(D5_D3 - D5_D3p1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,ip1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3p1))

                 else

                    a1 = (five_thirds/(D5_D3 - D5_D3m1))* &
                         (recv(ivar,ip,jp,kp) - recv(ivar,im1,jp,kp) - &
                         0.75*a2*(D4_D3 - D4_D3m1))

                 endif


                 if (stencil == NORMAL) then
                    a5 = ( (A_jm1 - A_j)* &
                         (recv(ivar,ip,jp1,kp) - recv(ivar,ip,jp,kp)) - &
                         (A_jp1 - A_j)* &
                         (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp)) ) / &
                         ( (B_jp1 - B_j)*(A_jm1 - A_j) - &
                         (B_jm1 - B_j)*(A_jp1 - A_j))

                    a4 = (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp) - &
                         a5*(B_jm1 - B_j))/(A_jm1 - A_j)

                 else if (stencil == LEFT) then
                    a5 = (recv(ivar,ip,jp1,kp) - recv(ivar,ip,jp,kp)) / &
                         (B_jp1 - B_j)

                    a4 = 0.0

                 else if (stencil == RIGHT) then
                    a5 = (recv(ivar,ip,jm1,kp) - recv(ivar,ip,jp,kp)) / &
                         (B_jm1 - B_j)

                    a4 = 0.0
                 endif

                 a7 = ((Dphi3_Dphip1 - Dphi3_Dphi)* &
                        (recv(ivar,ip,jp,km1) - recv(ivar,ip,jp,kp)) - &
                       (Dphi3_Dphim1 - Dphi3_Dphi)* &
                        (recv(ivar,ip,jp,kp1) - recv(ivar,ip,jp,kp))) / &
                      ((zcm1 - zcc)*(Dphi3_Dphip1 - Dphi3_Dphi) - &
                       (zcp1 - zcc)*(Dphi3_Dphim1 - Dphi3_Dphi))

                 a6 = (3.0/(Dphi3_Dphip1 - Dphi3_Dphi))* &
                      (recv(ivar,ip,jp,kp1) - recv(ivar,ip,jp,kp) - a7*(zcp1 - zcc))

                 a3 = recv(ivar,ip,jp,kp) - &
                      0.6*a1*D5_D3 - &
                      0.75*a2*D4_D3 - &
                      a4*A_j -&
                      a5*B_j


! the interpolating polynomial is integrated over the child cell volume to
! get the prolonged cell-average value for that child.  Compute that now
                 unk(ivar,i,j,k,isg) = 0.6*a1*D5_D3_child + &
                                       0.75*a2*D4_D3_child + &
                                       a3 + &
                                       a4*A_jchild + &
                                       a5*B_jchild + &
                                       a6*Dphi3_Dphi_child + &
                                       a7*ZCOORD(k)
                 
              enddo

#endif

           end select
            
        enddo
     enddo
  enddo


#ifdef ABUNDANCE_CHECK_INTERIORS
!-----------------------------------------------------------------------------
! compute the error in the nuclear abundances in the block *interiors*
! to see if we've screwed up.  
!-----------------------------------------------------------------------------

  do k = NGUARD*K3D+1, NGUARD*K3D+NZB
     do j = NGUARD*K2D+1, NGUARD*K2D+NYB
        do i = NGUARD+1, NGUARD+NXB

           error = -1.0
           do n = 1, NSPECIES
              error = error + min(1.0,max(smallx,unk(inuc_begin-1+n,i,j,k,isg)))
           enddo
           
           if (error > 1.e-5) print *, '*****  ERROR before monotonicity constraint = ', error, &
                ' in amr_prolong_gen_unk_fun.F90', i,j,k,isg
           
        enddo
     enddo
  enddo

#endif


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
  ipmin  = ((ia-1)/2) + 1 + (NGUARD/2) + ioff
  ipmax  = ((ib-1)/2) + 1 + (NGUARD/2) + ioff

  jpmin  = ((ja-1)/2) + 1 + (NGUARD/2)*K2D + joff
  jpmax  = ((jb-1)/2) + 1 + (NGUARD/2)*K2D + joff

  kpmin  = ((ka-1)/2) + 1 + (NGUARD/2)*K3D + koff
  kpmax  = ((kb-1)/2) + 1 + (NGUARD/2)*K3D + koff


! loop over the parent cells and make sure that the zone-average values
! of the children call between the limits set by the parents that 
! are contained in the prolongation stencil
  kchild1 = ka
  do k = kpmin, kpmax

     jchild1 = ja
     do j = jpmin, jpmax

        ichild1 = ia
        do i = ipmin, ipmax

! check the abundance conservation of the children of the current zone.
! doing linear prolongation does not always preserve sum{X_i} = 1 in
! the children

           monotonize_abundances = .false.

           error    = 0.e0
           xnuc_min = 1.e0

           if (gr_convertToConsvdForMeshCalls) then

! if we are dealing with conserved variables, then we are storing 
! partial densities now, not mass fractions.

              do kk = kchild1, kchild1+K3D
                 do jj = jchild1, jchild1+K2D
                    do ii = ichild1, ichild1+1

                       sum = -1.e0

                       do n = 1, NSPECIES
                          xnuc = unk(inuc_begin-1+n,ii,jj,kk,isg) &
                                /unk(DENS_VAR,ii,jj,kk,isg)
                          sum = sum + max(smallx,min(1.e0, xnuc ))
                          xnuc_min = min(xnuc_min, xnuc)
                       enddo

                       error = max(error, abs(sum))

                    enddo
                 enddo
              enddo


           else

              do kk = kchild1, kchild1+K3D
                 do jj = jchild1, jchild1+K2D
                    do ii = ichild1, ichild1+1

                       sum = -1.e0

                       do n = 1, NSPECIES
                          xnuc = unk(inuc_begin-1+n,ii,jj,kk,isg)
                          sum = sum + max(smallx,min(1.e0, xnuc ))
                          xnuc_min = min(xnuc_min, xnuc)
                       enddo

                       error = max(error, abs(sum))
                       
                    enddo
                 enddo
              enddo
              
           endif

           monotonize_abundances = error > 1.e-8 .OR. xnuc_min < xnuc_neg

!!$           if ( abs(error) > 1.e-8 ) then
!!$              print *, 'need to monotonize child: ',ii,jj,kk,isg
!!$           else
!!$              print *, 'need to reset child: ',ii,jj,kk,isg
!!$           end if


! check to see what variables fall outside the limits set by the parent
! stencil.  The logical variable monotonize(ivar) is true if variable
! ivar needs to be monotonized
           monotonize(:) = .false.

! don't monotonize the abundances unless they don't sum to 1.  If we do 
! want to monotonize the abundances, then we must make sure that if we
! do one, that we do *ALL* of them, otherwise, we will break the summing
! to 1.
           do ivar = 1, inuc_begin-1 

#if N_DIM == 1 
! for this 1-d spherical prolongation stencil, 3 parents were involved
              parent_min = min(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k))

              parent_max = max(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k))

! now find the extrema of the children 
              child_min = min(unk(ivar,ichild1,  jchild1,    kchild1,isg),&
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg))

              child_max = max(unk(ivar,ichild1,  jchild1,    kchild1,isg),&
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg))

#elif N_DIM == 2
! for this 2-d spherical prolongation stencil, 5 parents were involved
              parent_min = min(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k), &
                               recv(ivar,i,  j+1,k), &
                               recv(ivar,i,  j-1,k))

              parent_max = max(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k), &
                               recv(ivar,i,  j+1,k), &
                               recv(ivar,i,  j-1,k))

! now find the extrema of the children 
              child_min = min(unk(ivar,ichild1,  jchild1,    kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg))

              child_max = max(unk(ivar,ichild1,  jchild1,    kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg))
#else

! for this 3-d spherical prolongation stencil, 7 parents were involved
              parent_min = min(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k), &
                               recv(ivar,i,  j+1,k), &
                               recv(ivar,i,  j-1,k), &
                               recv(ivar,i,  j,  k+1), &
                               recv(ivar,i,  j,  k-1))

              parent_max = max(recv(ivar,i,  j,  k), &
                               recv(ivar,i+1,j,  k), &
                               recv(ivar,i-1,j,  k), &
                               recv(ivar,i,  j+1,k), &
                               recv(ivar,i,  j-1,k), &
                               recv(ivar,i,  j,  k+1), &
                               recv(ivar,i,  j,  k-1))

! now find the extrema of the children 
              child_min = min(unk(ivar,ichild1,  jchild1,    kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1,  jchild1,    kchild1+K3D,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1+K3D,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1+K3D,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1+K3D,isg))

              child_max = max(unk(ivar,ichild1,  jchild1,    kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg), &
                              unk(ivar,ichild1,  jchild1,    kchild1+K3D,isg), &
                              unk(ivar,ichild1+1,jchild1,    kchild1+K3D,isg), &
                              unk(ivar,ichild1,  jchild1+K2D,kchild1+K3D,isg), &
                              unk(ivar,ichild1+1,jchild1+K2D,kchild1+K3D,isg))

#endif

! check to see if the extrema of the children go out of the bounds of the 
! parents
              if (child_min < parent_min .OR. &
                  child_max > parent_max) then
                 monotonize(ivar) = .true.
              endif

                 
           enddo

           if (monotonize_abundances) monotonize(inuc_begin:inuc_end) = .true.

           do ivar = 1, NUNK_VARS

              if (monotonize(ivar)) then


! reset the children
#if N_DIM == 1
                 unk(ivar,ichild1,  jchild1,    kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1,    kchild1,isg) = recv(ivar,i,j,k)
#elif N_DIM == 2
                 unk(ivar,ichild1,  jchild1,    kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1,    kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1,  jchild1+K2D,kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg) = recv(ivar,i,j,k)
#else
                 unk(ivar,ichild1,  jchild1,    kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1,    kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1,  jchild1+K2D,kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1+K2D,kchild1,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1,  jchild1,    kchild1+K3D,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1,    kchild1+K3D,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1,  jchild1+K2D,kchild1+K3D,isg) = recv(ivar,i,j,k)
                 unk(ivar,ichild1+1,jchild1+K2D,kchild1+K3D,isg) = recv(ivar,i,j,k)
#endif
                 
              endif
                             
           enddo

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
  ipmin  = ((ia-1)/2) + 1 + (NGUARD/2) + ioff
  ipmax  = ((ib-1)/2) + 1 + (NGUARD/2) + ioff

  jpmin  = ((ja-1)/2) + 1 + (NGUARD/2)*K2D + joff
  jpmax  = ((jb-1)/2) + 1 + (NGUARD/2)*K2D + joff

  kpmin  = ((ka-1)/2) + 1 + (NGUARD/2)*K3D + koff
  kpmax  = ((kb-1)/2) + 1 + (NGUARD/2)*K3D + koff


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
                   'ERROR: cartesian not implemented in amr_prolong_gen_unk_fun.F90')

           case (CYLINDRICAL)
              
              call Driver_abortFlash( &
                   'ERROR: cylindrical not implemented in amr_prolong_gen_unk_fun.F90')
              
           case (SPHERICAL)

#if N_DIM == 1
              cvol(1) = (XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + &
                    XCOORD_LEFT(ichild1)**2)

              cvol(2) = (XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + &
                    XCOORD_LEFT(ichild1+1)**2)

              pvol = cvol(1) + cvol(2)

#elif N_DIM == 2
              
              cvol(1) = two_thirds*PI*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))

              cvol(2) = two_thirds*PI*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))

              cvol(3) = two_thirds*PI*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))

              cvol(4) = two_thirds*PI*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))


! force conservation
              pvol = cvol(1) + cvol(2) + cvol(3) + cvol(4)

#else
              cvol(1) = one_third*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))* &
                   (ZCOORD_RIGHT(kchild1) - ZCOORD_LEFT(kchild1))

              cvol(2) = one_third*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))* &
                   (ZCOORD_RIGHT(kchild1) - ZCOORD_LEFT(kchild1))

              cvol(3) = one_third*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))* &
                   (ZCOORD_RIGHT(kchild1) - ZCOORD_LEFT(kchild1))

              cvol(4) = one_third*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))* &
                   (ZCOORD_RIGHT(kchild1) - ZCOORD_LEFT(kchild1))

              cvol(5) = one_third*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))* &
                   (ZCOORD_RIGHT(kchild1+1) - ZCOORD_LEFT(kchild1+1))

              cvol(6) = one_third*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1)) - cos(YCOORD_RIGHT(jchild1)))* &
                   (ZCOORD_RIGHT(kchild1+1) - ZCOORD_LEFT(kchild1+1))

              cvol(7) = one_third*(XCOORD_RIGHT(ichild1) - XCOORD_LEFT(ichild1))* &
                   (XCOORD_RIGHT(ichild1)**2 + XCOORD_RIGHT(ichild1)*XCOORD_LEFT(ichild1) + XCOORD_LEFT(ichild1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))* &
                   (ZCOORD_RIGHT(kchild1+1) - ZCOORD_LEFT(kchild1+1))

              cvol(8) = one_third*(XCOORD_RIGHT(ichild1+1) - XCOORD_LEFT(ichild1+1))* &
                   (XCOORD_RIGHT(ichild1+1)**2 + XCOORD_RIGHT(ichild1+1)*XCOORD_LEFT(ichild1+1) + XCOORD_LEFT(ichild1+1)**2)* &
                   (cos(YCOORD_LEFT(jchild1+1)) - cos(YCOORD_RIGHT(jchild1+1)))* &
                   (ZCOORD_RIGHT(kchild1+1) - ZCOORD_LEFT(kchild1+1))


! force conservation
              pvol = cvol(1) + cvol(2) + cvol(3) + cvol(4) + &
                     cvol(5) + cvol(6) + cvol(7) + cvol(8)

#endif


           end select

           do ivar = 1, NUNK_VARS

              select case (gr_geometry)

              case (CARTESIAN)

                 call Driver_abortFlash( &
                      'ERROR: cartesian not implemented in amr_prolong_gen_unk_fun.F90')
                 
              case (CYLINDRICAL)

                 call Driver_abortFlash( &
                      'ERROR: cylindrical not implemented in amr_prolong_gen_unk_fun.F90')

              case (SPHERICAL)

#if N_DIM == 1
                 error = (recv(ivar,i,j,k) - &
                      (cvol(1)*unk(ivar,ichild1,  jchild1,  kchild1,isg) + &
                       cvol(2)*unk(ivar,ichild1+1,jchild1,  kchild1,isg))/pvol)/&
                      recv(ivar,i,j,k)
              
#elif N_DIM == 2
                 error = (recv(ivar,i,j,k) - &
                      (cvol(1)*unk(ivar,ichild1,  jchild1,  kchild1,isg) + &
                       cvol(2)*unk(ivar,ichild1+1,jchild1,  kchild1,isg) + &
                       cvol(3)*unk(ivar,ichild1  ,jchild1+1,kchild1,isg) + &
                       cvol(4)*unk(ivar,ichild1+1,jchild1+1,kchild1,isg))/pvol)/&
                      recv(ivar,i,j,k)
#else
                 error = (recv(ivar,i,j,k) - &
                      (cvol(1)*unk(ivar,ichild1,  jchild1,  kchild1,isg) + &
                       cvol(2)*unk(ivar,ichild1+1,jchild1,  kchild1,isg) + &
                       cvol(3)*unk(ivar,ichild1  ,jchild1+1,kchild1,isg) + &
                       cvol(4)*unk(ivar,ichild1+1,jchild1+1,kchild1,isg) + &
                       cvol(5)*unk(ivar,ichild1,  jchild1,  kchild1+1,isg) + &
                       cvol(6)*unk(ivar,ichild1+1,jchild1,  kchild1+1,isg) + &
                       cvol(7)*unk(ivar,ichild1  ,jchild1+1,kchild1+1,isg) + &
                       cvol(8)*unk(ivar,ichild1+1,jchild1+1,kchild1+1,isg))/pvol)/&
                      recv(ivar,i,j,k)

#endif   

              end select

              if (error > 1.e-10 .AND. recv(ivar,i,j,k) /= 0.0) then
!!                 print *, 'error = ', error, flashVars(ivar), recv(ivar,i,j,k)
                 print *, 'error = ', error, ivar, recv(ivar,i,j,k)
                 print *, 'x, y, z = ', XCOORD(ichild1), YCOORD(jchild1), ZCOORD(kchild1)
!!                 print *, flashVars
              endif
              
           enddo

           ichild1 = ichild1 + 2
        enddo

        jchild1 = jchild1 + 2*K2D
     enddo

     kchild1 = kchild1 + 2*K3D
  enddo
  
#endif


#ifdef ABUNDANCE_CHECK_INTERIORS
!-----------------------------------------------------------------------------
! compute the error in the nuclear abundances in the block *interiors*
! to see if we've screwed up
!-----------------------------------------------------------------------------

  do k = NGUARD*K3D+1, NGUARD*K3D+NZB
     do j = NGUARD*K2D+1, NGUARD*K2D+NYB
        do i = NGUARD+1, NGUARD+NXB

           error = -1.0
           do n = 1, NSPECIES
              error = error + min(1.0,max(smallx,unk(inuc_begin-1+n,i,j,k,isg)))
           enddo
           
           if (error > 1.e-5) then
              print *, '*****  ERROR after prolongation = ', error, ' in amr_prolong_gen_unk_fun.F90', i,j,k,isg
!              print *, 'abundances = '
!              do n = 1, NSPECIES
!                 print *, unk(inuc_begin-1+n,i,j,k,isg)
!              enddo
           endif

        enddo
     enddo
  enddo

#endif

  return
end subroutine amr_prolong_gen_unk_fun

