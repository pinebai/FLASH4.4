!!****if* source/Grid/GridMain/paramesh/Paramesh2/amr_restrict_red
!!
!! NAME
!!
!!  amr_restrict_red
!!
!!
!! SYNOPSIS
!!
!!  amr_restrict_red(icoord, block_parent, ioff, joff, koff)
!!
!!  amr_restrict_red(integer, integer, integer, integer, integer)
!!
!!
!! DESCRIPTION
!!
!!  Perform an area weighted flux averaging at block boundaries where
!!  there is a jump in refinement.  
!!
!!  In amr_restrict_bnd_data, we move the fluxes in the children of a
!!  block which is at a jump in refinement to the processor that their
!!  parent is stored kon.  These child fluxes are stored in recvar{x,y,z}.
!!
!!  At a coarse-fine interface, we want to take the child fluxes as the
!!  more accurate ones, and find the area-weighted average of those fluxes
!!  and pass it to the neighboring coarser block.  This is accomplished by
!!  first averaging the child fluxes and passing the result to their 
!!  parent.  This is then passed to the neighboring block of the children/
!!  parent.
!!
!!  The averaging should take into account the area of the face that the
!!  fluxes pass through.
!!
!!  In this routine, the parent data is on the local processor, so we can
!!  only get the coordinate information for that.  We will need to derive
!!  the child coordinate data from this.  Note, to ensure conservation,
!!  the sum of the child areas must equal the area of the parent.
!!
!!  Note that this does not update guard cell elements of bndtempx(y)(z).
!!
!!  This particular version is only appropriate for 2nd order schemes 
!!  using linear interpolation with even number of mesh points along 
!!  each block axis.
!!
!!***

subroutine amr_restrict_red(icoord, block_parent, ioff, joff, koff)
  
  ! the flux information is passed through common blocks in physicaldata
  use physicaldata
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_data, ONLY : gr_geometry
#include "constants.h"
#include "Flash.h"

  implicit none
  
  integer :: icoord
  
  integer :: ioff, joff, koff
  integer :: block_parent
  
  integer :: i, j, k, ivar
  integer :: ipar, jpar, kpar
  
  real, dimension(2*nguard+nxb) :: x, xl, xr
  real, dimension(2*nguard*k2d+nyb) :: y, yl, yr
  real, dimension(2*nguard*k3d+nzb) :: z, zl, zr

  real, dimension(4) :: carea
  real :: parea

  logical :: validGeom, gcell=.true.
  integer :: csize

!------------------------------------------------------------------------------
#ifdef DEBUG
  validGeom =  ((gr_geometry == CARTESIAN)                      .OR. &
       ((gr_geometry == CYLINDRICAL) .AND. (ndim == 2)) .OR. &
       ((gr_geometry == SPHERICAL))                     .OR. &
       ((gr_geometry == POLAR)       .AND. (ndim == 2)))
  
  if (.NOT. validGeom) then
     call Driver_abortFlash("[amr_restrict_red] ERROR: invalid geometry.")
  endif
  
#endif

! start by getting the coordinates of the parent block -- remember
! we are only filling a fraction of the parent here.  ioff, joff, koff
! are the offsets (in terms of cells) of the child into the parent
! block
  csize = 2*nguard+nxb
  call Grid_getCellCoords(IAXIS, block_parent, CENTER, gcell,   x,csize)
  call Grid_getCellCoords(IAXIS, block_parent, LEFT_EDGE, gcell,  xl,csize)
  call Grid_getCellCoords(IAXIS, block_parent, RIGHT_EDGE, gcell,  xr,csize)
  csize = 2*nguard*k2d+nyb
  call Grid_getCellCoords(JAXIS, block_parent, CENTER, gcell,   y,csize)
  call Grid_getCellCoords(JAXIS, block_parent, LEFT_EDGE, gcell,  yl,csize)
  call Grid_getCellCoords(JAXIS, block_parent, RIGHT_EDGE, gcell,  yr,csize)
  csize = 2*nguard*k3d+nzb
  call Grid_getCellCoords(KAXIS, block_parent, CENTER, gcell, z,csize)
  call Grid_getCellCoords(KAXIS, block_parent, LEFT_EDGE, gcell, zl,csize)
  call Grid_getCellCoords(KAXIS, block_parent, RIGHT_EDGE, gcell, zr,csize)


!-----------------------------------------------------------------------------
! flux going through the x-face
!-----------------------------------------------------------------------------
  if(icoord == 1) then                       

     select case (gr_geometry)

     case (CARTESIAN)

! cartesian is easy, the parent is twice as big as the child -- just do
! straight averaging.  Note, we write this in a dimension-independent
! fashion

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = 0.25*(        & 
                         recvarx1(ivar,i,j,    k)     +   & 
                         recvarx1(ivar,i,j+k2d,k)     +   & 
                         recvarx1(ivar,i,j,    k+k3d) +   & 
                         recvarx1(ivar,i,j+k2d,k+k3d)) 

                 enddo

              enddo
           enddo
        enddo


     case (CYLINDRICAL)

! this is 2-d cylindrical geometry only.  The area of the x-face is 2*pi*r*dz,
! but that 'r' is common to the coarse and fine cells that share the 
! interface, so we basically reduce down to a Cartesian average here.

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = 0.5*(      & 
                         recvarx1(ivar,i,j,    k)  +   & 
                         recvarx1(ivar,i,j+k2d,k))

                 enddo

              enddo
           enddo
        enddo


     case (SPHERICAL)

#if N_DIM == 1 
! this is 1-d spherical only, so things are very, very simple.  The area of 
! the x-face is 4*pi*r**2, but, as with the cylindrical case, the r is shared 
! by both the fine and coarse cells that share that interface, so the 
! averaging is just the same as the Cartesian.

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = recvarx1(ivar,i,j,k)

                 enddo

              enddo
           enddo
        enddo

#elif N_DIM == 2
! here we are going through the r-face, which is just an arc

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 2-d only here (r,theta)
! We always store the fluxes at the minimum face for a zone.  
!
! The area is just the integral over \phi (which yields 2\pi) and over \theta
! (which must include the sin \theta term from the \phi integral):
!
!  A = 2 \pi \int_{\theta_{j-1/2}}^{\theta_{j+1/2}} r_{i-1/2}^2 sin \theta d\theta
!
!    = 2 \pi r_{i-1/2}^2 { cos(\theta_{j-1/2}) - cos(\theta_{j+1/2}) }
!
! the '2\pi r^2' term just cancels out, since it is the same for all faces

                 carea(1) = (cos(yl(jpar)) - cos(y(jpar)))
                 carea(2) = (cos(y(jpar)) - cos(yr(jpar)))

                 parea = carea(1) + carea(2)

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = (     & 
                         carea(1)*recvarx1(ivar,i,j,k)     +  & 
                         carea(2)*recvarx1(ivar,i,j+1,k))/parea

                 enddo

              enddo
           enddo
        enddo

        

#else
! here we are going through the r-face, which is just an arc

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 3-d here (r,theta,phi)
! We always store the fluxes at the minimum face for a zone.  
!
! The area is just the integral over \phi and over \theta
! (which must include the sin \theta term from the \phi integral):
!
!  A = \int_{\theta_{j-1/2}}^{\theta_{j+1/2}} 
!      \int_{\phi_{k-1/2}}^{\phi_{k+1/2}} r_{i-1/2}^2 sin \theta d\theta d\phi
!
!    = r_{i-1/2}^2 ( cos(\theta_{j-1/2}) - cos(\theta_{j+1/2}) )
!      (\phi_{k+1/2} - \phi_{k-1/2})
!
! the 'r^2' term just cancels out, since it is the same for all faces

                 carea(1) = (cos(yl(jpar)) - cos(y(jpar)))*(z(kpar) - zl(kpar))
                 carea(2) = (cos(y(jpar)) - cos(yr(jpar)))*(z(kpar) - zl(kpar))
                 carea(3) = (cos(yl(jpar)) - cos(y(jpar)))*(zr(kpar) - z(kpar))
                 carea(4) = (cos(y(jpar)) - cos(yr(jpar)))*(zr(kpar) - z(kpar))

                 parea = carea(1) + carea(2) + carea(3) + carea(4)

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = (     & 
                         carea(1)*recvarx1(ivar,i,j,  k)   +  & 
                         carea(2)*recvarx1(ivar,i,j+1,k)   +  &
                         carea(3)*recvarx1(ivar,i,j,  k+1) +  &
                         carea(4)*recvarx1(ivar,i,j+1,k+1))/parea

                 enddo

              enddo
           enddo
        enddo

#endif

      case (POLAR)
! this is 2-d polar geometry only.  The area of the x-face is r*dphi,
! but that 'r' is common to the coarse and fine cells that share the 
! interface, so we basically reduce down to a Cartesian average here.

#if N_DIM == 2   

         do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
             do i=1,2
               do ivar = 1,nfluxes

                 bndtempx1(ivar,i,j,k) = 0.5*(      & 
                      recvarx1(ivar,i,j,    k)  +   & 
                      recvarx1(ivar,i,j+k2d,k))

               end do
             end do
           end do
         end do

#else
! here we are going through the r-face, which is just an arc

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1+nguard*k2d,nyb+nguard*k2d,2
              do i=1,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 3-d here (r,theta,phi)
! We always store the fluxes at the minimum face for a zone.  
!
! The area is just the integral over \phi and over \theta
! (which must include the sin \theta term from the \phi integral):
!
!  A = \int_{\theta_{j-1/2}}^{\theta_{j+1/2}} 
!      \int_{\phi_{k-1/2}}^{\phi_{k+1/2}} r_{i-1/2}^2 sin \theta d\theta d\phi
!
!    = r_{i-1/2}^2 ( cos(\theta_{j-1/2}) - cos(\theta_{j+1/2}) )
!      (\phi_{k+1/2} - \phi_{k-1/2})
!
! the 'r^2' term just cancels out, since it is the same for all faces

                 carea(1) = (cos(yl(jpar)) - cos(y(jpar)))*(z(kpar) - zl(kpar))
                 carea(2) = (cos(y(jpar)) - cos(yr(jpar)))*(z(kpar) - zl(kpar))
                 carea(3) = (cos(yl(jpar)) - cos(y(jpar)))*(zr(kpar) - z(kpar))
                 carea(4) = (cos(y(jpar)) - cos(yr(jpar)))*(zr(kpar) - z(kpar))

                 parea = carea(1) + carea(2) + carea(3) + carea(4)

                 do ivar=1,nfluxes

                    bndtempx1(ivar,i,j,k) = (     & 
                         carea(1)*recvarx1(ivar,i,j,  k)   +  & 
                         carea(2)*recvarx1(ivar,i,j+1,k)   +  &
                         carea(3)*recvarx1(ivar,i,j,  k+1) +  &
                         carea(4)*recvarx1(ivar,i,j+1,k+1))/parea

                 enddo

              enddo
           enddo
        enddo

#endif

     end select


!-----------------------------------------------------------------------------
! flux going through the y-face
!-----------------------------------------------------------------------------
  elseif (icoord == 2) then                    

     select case (gr_geometry) 

     case (CARTESIAN)

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1,2
              do i=1+nguard,nxb+nguard,2

                 do ivar=1,nfluxes

                    bndtempy1(ivar,i,j,k) = 0.25*(     & 
                         recvary1(ivar,i,  j,k)     +  & 
                         recvary1(ivar,i+1,j,k)     +  & 
                         recvary1(ivar,i,  j,k+k3d) +  & 
                         recvary1(ivar,i+1,j,k+k3d)) 

                 enddo

              enddo
           enddo
        enddo


     case (CYLINDRICAL)

!
! in 2-d cylindrical coordinates, we are going through an annular face, so
! we need to compute the area of the child faces and the parent, and find
! the weighted average.
!
!
!   +-------------+
!   |             |
!   |             |
!   |       F     |            ^ z
!   |      ^      |            |
!   |   f1 |   f2 |            |
!   |   ^  |   ^  |            +----->  r
!   |   |  |   |  |
!   +------+------+
!   |      .      |
!   |      .      |
!   |      .      |     f1 and f2 are the fluxes from the children.  We
!   +......+......+     want to find F, which is the flux from the parent
!   |      .      |     of those children.  The area that f1 and f2 go 
!   |      .      |     through is actually an annulus, since we are in 
!   |      .      |     cylindrical geometry.
!   +------+------+
!

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1,2
              do i=1+nguard,nxb+nguard,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 2-d only here (r,z)
! this is just 2*pi*(r_r**2 - r_l**2) = 2*pi*(r_r + r_l)*(r_r - r_l).
                 carea(1) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))
                 carea(2) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))

                 parea = carea(1) + carea(2)

                 do ivar=1,nfluxes

                    bndtempy1(ivar,i,j,k) = (     & 
                         carea(1)*recvary1(ivar,i,  j,k)     +  & 
                         carea(2)*recvary1(ivar,i+1,j,k))/parea

                 enddo

              enddo
           enddo
        enddo


     case (SPHERICAL)

#if N_DIM == 2 

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1,2
              do i=1+nguard,nxb+nguard,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 2-d only here (r,theta)
!
! the area is just the integral over \phi and r,
!
! A = 2 \pi \int_{r_{i-1/2}}^{r_{i+1/2}} r sin \theta_{j-1/2} dr
!
!   = 2 \pi r_i \Delta r sin \theta_{j-1/2}
!
! The 2 \pi sin(\theta_{j-1/2}) term is common for all faces, and cancels out

                 carea(1) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))
                 carea(2) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))

                 parea = carea(1) + carea(2)

                 do ivar=1,nfluxes

                    bndtempy1(ivar,i,j,k) = (     & 
                         carea(1)*recvary1(ivar,i,  j,k)     +  & 
                         carea(2)*recvary1(ivar,i+1,j,k))/parea

                 enddo

              enddo
           enddo
        enddo

#else

        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1,2
              do i=1+nguard,nxb+nguard,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 3-d here (r,theta,phi)
!
! the area is just the integral over \phi and r,
!
! A = \int_{r_{i-1/2}}^{r_{i+1/2}} 
!     \int_{\phi_{k-1/2}}^{\phi_{k+1/2}} r sin \theta_{j-1/2} dr d\phi
!
!   = r_i \Delta r sin \theta_{j-1/2} \Delta \phi
!
! The sin(\theta_{j-1/2}) term is common for all faces, and cancels out

                 carea(1) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))* &
                      (z(kpar) - zl(kpar))
                 carea(2) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))* &
                      (z(kpar) - zl(kpar))
                 carea(3) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))* &
                      (zr(kpar) - z(kpar))
                 carea(4) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))* &
                      (zr(kpar) - z(kpar))

                 parea = carea(1) + carea(2) + carea(3) + carea(4)

                 do ivar=1,nfluxes

                    bndtempy1(ivar,i,j,k) = (     & 
                         carea(1)*recvary1(ivar,i,  j,k)    + & 
                         carea(2)*recvary1(ivar,i+1,j,k)    + &
                         carea(3)*recvary1(ivar,i,  j,k+1)  + & 
                         carea(4)*recvary1(ivar,i+1,j,k+1)  &
                         )/parea

                 enddo

              enddo
           enddo
        enddo

#endif
        
      case (POLAR)
        
        do k=1+nguard*k3d,nzb+nguard*k3d,2
           do j=1,2
              do i=1+nguard,nxb+nguard,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2

! start by computing the child areas -- remember, we are 2-d only here (r,phi)
! this is just (r_r - r_l).
                 carea(1) = (x(ipar) - xl(ipar))
                 carea(2) = (xr(ipar) - x(ipar))

                 parea = carea(1) + carea(2)

                 do ivar=1,nfluxes

                    bndtempy1(ivar,i,j,k) = (     & 
                         carea(1)*recvary1(ivar,i,  j,k)     +  & 
                         carea(2)*recvary1(ivar,i+1,j,k))/parea

                 enddo

              enddo
           enddo
        enddo

     end select


!-----------------------------------------------------------------------------
! flux going through the z-face
!-----------------------------------------------------------------------------
  elseif(icoord == 3) then                    

     select case(gr_geometry)

     case (CARTESIAN)

        do k=1,2
           do j=1+nguard,nyb+nguard,2
              do i=1+nguard,nxb+nguard,2

                 do ivar=1,nfluxes

                    bndtempz1(ivar,i,j,k) = 0.25*(  & 
                         recvarz1(ivar,i,  j,  k) + & 
                         recvarz1(ivar,i+1,j,  k) + & 
                         recvarz1(ivar,i,  j+1,k) + & 
                         recvarz1(ivar,i+1,j+1,k))

                 enddo

              enddo
           enddo
        enddo


     case (CYLINDRICAL)

        call Driver_abortFlash("[amr_restrict_red] ERROR: invalid geometry.")
     

     case (SPHERICAL)
   
        do k=1,2
           do j=1+nguard,nyb+nguard,2
              do i=1+nguard,nxb+nguard,2

! compute the location of the parent current cell in the parent block
                 ipar = nguard + 1 + ioff + (i - nguard)/2
                 jpar = nguard*k2d + 1 + joff + (j - nguard*k2d)/2
                 kpar = nguard*k3d + 1 + koff + (k - nguard*k3d)/2
           

! start by computing the child areas -- remember, we are 3-d here (r,theta,phi)
!
! the area is just the integral over \theta and r,
!
! A = \int_{r_{i-1/2}}^{r_{i+1/2}} 
!     \int_{\theta_{j-1/2}}^{\theta_{j+1/2}} r dr d\theta
!
!   = r_i \Delta r \Delta \theta
!

                 carea(1) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))* &
                      (y(jpar) - yl(jpar))
                 carea(2) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))* &
                      (y(jpar) - yl(jpar))
                 carea(3) = (xl(ipar) + x(ipar))*(x(ipar) - xl(ipar))* &
                      (yr(jpar) - y(jpar))
                 carea(4) = (x(ipar) + xr(ipar))*(xr(ipar) - x(ipar))* &
                      (yr(jpar) - y(jpar))

                 parea = carea(1) + carea(2) + carea(3) + carea(4)

                 do ivar=1,nfluxes

                    bndtempz1(ivar,i,j,k) = (     & 
                         carea(1)*recvarz1(ivar,i,  j,k)    + & 
                         carea(2)*recvarz1(ivar,i+1,j,k)    + &
                         carea(3)*recvarz1(ivar,i,  j+1,k)  + & 
                         carea(4)*recvarz1(ivar,i+1,j+1,k)  &
                         )/parea

                 enddo

              enddo
           enddo
        enddo


     end select
  
  endif

  return
end subroutine amr_restrict_red
