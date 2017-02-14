!!****if* source/Grid/GridMain/paramesh/Paramesh2/amr_restrict_work_fun
!!
!! NAME
!!
!!  amr_restrict_work_fun
!!
!! 
!! SYNOPSIS
!!
!!  amr_restrict_work_fun(datain, dataout, block_parent, ioff, joff, koff)
!!
!!  amr_restrict_work_fun(real(), real(), integer, integer, integer, integer)
!!
!!
!! DESCRIPTION
!!
!!  Restrict the data from the children up to the parent block.  This is a 
!!  modified version of the paramesh restriction routine that takes into 
!!  account the geometry of the mesh by using a volume-weighted average.
!!
!!  It is always the case that the parent is on the local processor (see 
!!  amr_restrict_cc), so we can use the dBase to get the coordinate 
!!  information for the parent block.  We can derive the child coordinate
!!  information from this (we know there is a 2x change in resolution),
!!  and the offsets of the child block into the parent (ioff, joff, koff)
!!  We cannot simply grab the dBase coordinate information for the children,
!!  since they may be on different processors.
!!
!!  
!! ARGUMENTS
!!
!!  datain        the data for the entire block of one of the children of 
!!                the parent work array.
!!
!!  dataout       the parent block buffer, this is where we are placing the
!!                restricted data.  This is for the work array, not unk
!!
!!  block_parent  the block ID (on the local processor) of the parent
!!
!!  ioff          the offsets (in parent zones, excluding guardcells) of the
!!  joff          current child block into the parent block, in each 
!!  koff          coordinate direction.
!!
!!
!! NOTES
!!
!!  For now, we support cartesian, 2-d cylindrical, spherical, 2-d polar
!!  geometries only.
!!
!!  The interface for this routine has been modified from the standard 
!!  paramesh verion, and the corresponding modifications were made to the
!!  amr_restrict_cc routine.
!!
!!  The dataout array (parent) has twice as many slots (in each direction) for 
!!  data than the restricted child block (datain) will provide.  This is 
!!  because the parent and child differ in resolution by a factor of 2, but
!!  are logically identical in layout.  We fill every other slot in dataout
!!  when doing restriction, and amr_restrict_cc handles the mapping.
!!
!!***


subroutine amr_restrict_work_fun(datain, dataout, block_parent, ioff, joff, koff)

  use physicaldata, ONLY :nxb, nyb, nzb, k2d, k3d, nguard, nvar,ndim
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getCellCoords
  use Grid_data, ONLY : gr_geometry
  use tree, ONLY : mchild
  use workspace

  implicit none
#include "constants.h"

  real datain(ilw:iuw,jlw:juw,klw:kuw)
  real dataout(ilw:iuw,jlw:juw,klw:kuw)

  integer i, j, k

  integer :: block_parent

  integer :: ioff, joff, koff

  integer :: ipar, jpar, kpar

  real :: pvol_i
  real, DIMENSION(mchild) :: cvol
  integer :: csize

! compute the maximum length of a vector in each coordinate direction 
! (including guardcells)
  real, dimension(1:2*nguard+nxb) :: x, xl, xr
  real, dimension(1:2*nguard*k2d+nyb) :: y, yl, yr
  real, dimension(1:2*nguard*k3d+nzb) :: z, zl, zr

  real, PARAMETER :: one_third = 1.0/3.0, &
                     two_thirds = 2.0/3.0

  real, PARAMETER :: pi = 3.141592653589793238E0

  logical :: validGeom, gcell = .true.


!------------------------------------------------------------------------------

! check for validity of the geometry
  validGeom =  ((gr_geometry == CARTESIAN)                      .OR. &
       ((gr_geometry == CYLINDRICAL) .AND. (ndim == 2)) .OR. &
       ((gr_geometry == SPHERICAL))                     .OR. &
       ((gr_geometry == POLAR)       .AND. (ndim == 2)))
  
  if (.NOT. validGeom) then
     call Driver_abortFlash("[amr_restrict_work_fun] ERROR: invalid geometry.")
  endif

  !call PhysicalConstants_getConstant("pi", pi)

! loop over all of the cells in the parent that this child
! can fill (note, dataout is dimensioned the same as datain,
! so we skip every other value).  Use a volume weighted average
! to fill the parent given the child information.

! start by getting the coordinates of the parent block -- remember
! we are only filling a fraction of the parent here.  ioff, joff, koff
! are the offsets (in terms of cells) of the child into the parent
! block
  csize = 2*nguard+nxb
  call Grid_getCellCoords(KAXIS, block_parent, CENTER, gcell, z,csize)
  call Grid_getCellCoords(KAXIS, block_parent, LEFT_EDGE, gcell, zl,csize)
  call Grid_getCellCoords(KAXIS, block_parent, RIGHT_EDGE, gcell, zr,csize)
  csize = 2*nguard*k2d+nyb
  call Grid_getCellCoords(JAXIS, block_parent, CENTER, gcell, y,csize)
  call Grid_getCellCoords(JAXIS, block_parent, LEFT_EDGE, gcell, yl,csize)
  call Grid_getCellCoords(JAXIS, block_parent, RIGHT_EDGE, gcell, yr,csize)
  csize = 2*nguard*k3d+nzb
  call Grid_getCellCoords(IAXIS, block_parent, CENTER, gcell, x,csize)
  call Grid_getCellCoords(IAXIS, block_parent, LEFT_EDGE, gcell, xl,csize)
  call Grid_getCellCoords(IAXIS, block_parent, RIGHT_EDGE, gcell, xr,csize)
  
! fill the dataout array with the restricted info from the children -- note
! we only use indicies nguard+1, nguard+3, nguard+5, ... in the dataout array, 
! REGARDLESS of the location of the child in the parent block

  do k = 1+nguard_work*k3d, nzb+nguard_work*k3d, 2
     do j = 1+nguard_work*k2d, nyb+nguard_work*k2d, 2
        do i = 1+nguard_work, nxb+nguard_work, 2

! compute the location of the parent current cell in the parent block
! note, the position must be with respect to unk, not work, since that
! is what the coordinates are defined with -- this means use nguard
! as the offset
           ipar = nguard + 1 + ioff + (i - nguard_work)/2
           jpar = nguard*k2d + 1 + joff + (j - nguard_work*k2d)/2
           kpar = nguard*k3d + 1 + koff + (k - nguard_work*k3d)/2
           
           select case (gr_geometry)

!----------------------------------------------------------------------------
! 1-, 2-, and 3-d cartesian
!----------------------------------------------------------------------------
           case (CARTESIAN)

              dataout(i,j,k) = ( & 
                   datain(i,   j,     k    ) + & 
                   datain(i+1, j,     k    ) + & 
                   datain(i,   j+k2d, k    ) + & 
                   datain(i+1, j+k2d, k    ) + & 
                   datain(i,   j,     k+k3d) + & 
                   datain(i+1, j,     k+k3d) + & 
                   datain(i,   j+k2d, k+k3d) + & 
                   datain(i+1, j+k2d, k+k3d))*.125
              
!----------------------------------------------------------------------------
! 2-d cylindrical
!----------------------------------------------------------------------------
           case (CYLINDRICAL)

! compute the volume of the parent cell 
!             pvol = (xr(ipar)**2 - xl(ipar)**2)*(yr(jpar) - yl(jpar))

! compute the volume of the children -- note the center coords (x, y, z)
! are the right coord of some children and the left coord of others
!              cvol(1) = (x(ipar)**2 - xl(ipar)**2)*(y(jpar) - yl(jpar))
!              cvol(2) = (xr(ipar)**2 - x(ipar)**2)*(y(jpar) - yl(jpar))
!              cvol(3) = (x(ipar)**2 - xl(ipar)**2)*(yr(jpar) - y(jpar))
!              cvol(4) = (xr(ipar)**2 - x(ipar)**2)*(yr(jpar) - y(jpar))

              cvol(1) = (x(ipar) + xl(ipar))*(x(ipar) - xl(ipar))* &
                   (y(jpar) - yl(jpar))
              cvol(2) = (xr(ipar) + x(ipar))*(xr(ipar) - x(ipar))* &
                   (y(jpar) - yl(jpar))
              cvol(3) = (x(ipar) + xl(ipar))*(x(ipar) - xl(ipar))* &
                   (yr(jpar) - y(jpar))
              cvol(4) = (xr(ipar) + x(ipar))*(xr(ipar) - x(ipar))* &
                   (yr(jpar) - y(jpar))

! force conservation
              pvol_i = 1.e0/(cvol(1) + cvol(2) + cvol(3) + cvol(4))

! do the volume weighted averaging
              dataout(i,j,k) = ( & 
                   cvol(1)*datain(i,   j,     k) + & 
                   cvol(2)*datain(i+1, j,     k) + & 
                   cvol(3)*datain(i,   j+k2d, k) + & 
                   cvol(4)*datain(i+1, j+k2d, k) )*pvol_i

           case (SPHERICAL)

#if N_DIM == 1
!----------------------------------------------------------------------------
! 1-d spherical
!----------------------------------------------------------------------------
              
              cvol(1) = (x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)

              cvol(2) = (xr(ipar) - x(ipar))* &
                   (xr(ipar)**2 + xr(ipar)*x(ipar) + x(ipar)**2)

              pvol_i = 1.e0/(cvol(1) + cvol(2))

! do the volume weighted averaging
              dataout(i,j,k) = ( & 
                   cvol(1)*datain(i,   j, k) + & 
                   cvol(2)*datain(i+1, j, k) )*pvol_i
              
#elif N_DIM == 2
!----------------------------------------------------------------------------
! 2-d spherical (r, theta)
!----------------------------------------------------------------------------

! here, x() is the radial coordinate, and y() is the theta angle (in radians)
!
! the volume is -(2/3) pi (r_r**3 - r_l**3) (cos theta_r - cos theta_l)
!
              cvol(1) = two_thirds*pi*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))

              cvol(2) = two_thirds*pi*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))

              cvol(3) = two_thirds*pi*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))

              cvol(4) = two_thirds*pi*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))

! force conservation
              pvol_i = 1.e0/(cvol(1) + cvol(2) + cvol(3) + cvol(4))

! do the volume weighted averaging
              dataout(i,j,k) = ( & 
                   cvol(1)*datain(i,   j,     k) + & 
                   cvol(2)*datain(i+1, j,     k) + & 
                   cvol(3)*datain(i,   j+k2d, k) + & 
                   cvol(4)*datain(i+1, j+k2d, k) )*pvol_i
#else
!----------------------------------------------------------------------------
! 3-d spherical (r, theta, phi)
!----------------------------------------------------------------------------

! here, x() is the radial coordinate, and y() is the theta angle (in radians),
! and z() is the phi coordinate (in radians)
!
! the volume is -(1/3) pi (r_r**3 - r_l**3) (cos theta_r - cos theta_l) 
!                         (phi_r - phi_l)
!
              cvol(1) = one_third*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))* &
                   (z(kpar) - zl(kpar))

              cvol(2) = one_third*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))* &
                   (z(kpar) - zl(kpar))

              cvol(3) = one_third*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))* &
                   (z(kpar) - zl(kpar))

              cvol(4) = one_third*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))* &
                   (z(kpar) - zl(kpar))

              cvol(5) = one_third*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))* &
                   (zr(kpar) - z(kpar))

              cvol(6) = one_third*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(yl(jpar)) - cos(y(jpar)))* &
                   (zr(kpar) - z(kpar))

              cvol(7) = one_third*(x(ipar) - xl(ipar))* &
                   (x(ipar)**2 + x(ipar)*xl(ipar) + xl(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))* &
                   (zr(kpar) - z(kpar))

              cvol(8) = one_third*(xr(ipar) - x(ipar))* &
                   (x(ipar)**2 + x(ipar)*xr(ipar) + xr(ipar)**2)* &
                   (cos(y(jpar)) - cos(yr(jpar)))* &
                   (zr(kpar) - z(kpar))

! force conservation
              pvol_i = 1.e0/(cvol(1) + cvol(2) + cvol(3) + cvol(4) + &
                             cvol(5) + cvol(6) + cvol(7) + cvol(8) )

! do the volume weighted averaging
                 dataout(i,j,k) = ( & 
                      cvol(1)*datain(i,   j,     k    ) + & 
                      cvol(2)*datain(i+1, j,     k    ) + & 
                      cvol(3)*datain(i,   j+k2d, k    ) + & 
                      cvol(4)*datain(i+1, j+k2d, k    ) + &
                      cvol(5)*datain(i,   j,     k+k3d) + & 
                      cvol(6)*datain(i+1, j,     k+k3d) + & 
                      cvol(7)*datain(i,   j+k2d, k+k3d) + & 
                      cvol(8)*datain(i+1, j+k2d, k+k3d))*pvol_i
#endif

!----------------------------------------------------------------------------
! 2-d polar (r, phi)
!----------------------------------------------------------------------------
           case (POLAR)

! compute the volume of the parent cell: dphi*(r_r^2-r_l^2)
!              pvol = (xr(ipar)**2 - xl(ipar)**2)*(yr(jpar) - yl(jpar))

! compute the volume of the children -- note the center coords (x, y, z)
! are the right coord of some children and the left coord of others
!
! Note: these do not need to be real volumes, rather their relative values
!       should be correct.

              cvol(1) = (x (ipar) + xl(ipar)) &
                       *(x (ipar) - xl(ipar)) &
                       *(y (jpar) - yl(jpar))
              cvol(2) = (xr(ipar) + x (ipar)) &
                       *(xr(ipar) - x (ipar)) &
                       *(y (jpar) - yl(jpar))
              cvol(3) = (x (ipar) + xl(ipar)) &
                       *(x (ipar) - xl(ipar)) &
                       *(yr(jpar) - y (jpar))
              cvol(4) = (xr(ipar) + x (ipar)) &
                       *(xr(ipar) - x (ipar)) &
                       *(yr(jpar) - y (jpar))

! force conservation
              pvol_i = 1.e0/(cvol(1) + cvol(2) + cvol(3) + cvol(4))

! do the volume weighted averaging
              dataout(i,j,k) = ( & 
                   cvol(1)*datain(i,   j,     k) + & 
                   cvol(2)*datain(i+1, j,     k) + & 
                   cvol(3)*datain(i,   j+k2d, k) + & 
                   cvol(4)*datain(i+1, j+k2d, k) )*pvol_i

           end select

        enddo
     enddo
  enddo

  return
end subroutine amr_restrict_work_fun
