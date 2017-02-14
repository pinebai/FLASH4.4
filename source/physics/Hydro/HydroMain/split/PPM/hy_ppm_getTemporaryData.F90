!!****if* source/physics/Hydro/HydroMain/split/PPM/hy_ppm_getTemporaryData
!!
!! NAME
!!  hy_ppm_getTemporaryData
!!
!! SYNOPSIS
!!
!!  hy_ppm_getTemporaryData(integer, intent(in) :: axis,
!!                         integer, intent(in) :: blockID, 
!!                         integer, intent(in) :: size(3),
!!                         real, intent(out) :: area(:,:,:),
!!                         real, intent(out) :: dtdx(:,:,:),
!!                         real, intent(out) :: grav(:,:,:),
!!                         real, intent(out) :: ngrav(:,:,:),
!!                         real, intent(out) :: fict(:,:,:),
!!                         real, intent(out) :: areaLeft(:,:,:))
!!
!! DESCRIPTION 
!!  This routine is needed because of a peculiarity of the PPM Hydro
!!  unit in combination with AMR.  The update of the solution is
!!  done in two steps.  The first step updates all the interior cells,
!!  while the temporary quantities needed for updating boundary cells
!!  are saved. Then there is a flux conservation stage to take care of the
!!  fine/coarse boundaries.  The second step updates the solution in
!!  boundary cells.  Boundary cells are block cells (not guardcells!)
!!  that have a face that coincides with a part of the boundary of the block
!!  and that is normal to the direction of the sweep.  Interior cells are
!!  the remaining block cells, i.e., cells that do not face a guardcell in
!!  the sweep direction.
!!  This routine fetches the saved temporary quantities for performing the
!!  boundary cell update in the second step mentioned above.
!!
!!  ARGUMENTS
!!
!!  axis - the sweep direction 
!!  blockID - the local identification of the block
!!  size - determines the shape of the remaining arguments; each of
!!         those arrays has dimension(size(1),size(2),size(3)).
!!  area - temporary quantity needed for hy_ppm_updateSolution   
!!  dtdx - temporary quantity needed for hy_ppm_updateSolution   
!!  grav - temporary quantity needed for hy_ppm_updateSolution   
!!  ngrav - temporary quantity needed for hy_ppm_updateSolution   
!!  fict - temporary quantity needed for hy_ppm_updateSolution   
!!  areaLeft - temporary quantity needed for flux correction
!! 
!! SEE ALSO
!!
!!  hy_ppm_putTemporaryData
!!***

subroutine hy_ppm_getTemporaryData(axis, blockID, size, &
                               area, dtdx, grav, ngrav, fict, areaLeft)



  use Hydro_data, ONLY: hy_xarea,hy_xdtdx,hy_xgrav,hy_xngrav,hy_xfict,hy_xareaAtFaces
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_data, ONLY: hy_yarea,hy_ydtdy,hy_ygrav,hy_yngrav,hy_yfict,hy_yareaAtFaces
  use Hydro_data, ONLY: hy_zarea,hy_zdtdz,hy_zgrav,hy_zngrav,hy_zfict,hy_zareaAtFaces

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: axis, blockID
  integer, intent(in), dimension(3) :: size
  real, intent(out), dimension(size(1),size(2),size(3)) :: &
                    area,dtdx,grav,ngrav,fict,areaLeft
  integer :: sx,ex,sy,ey,sz,ez

#ifdef DEBUG
  if((blockid<1).or.(blockid>MAXBLOCKS)) then
     call Driver_abortFlash("invalid blockid for getting coordinates ")
  end if
#endif

  sx = NGUARD+1
  sy = NGUARD*K2D+1
  sz = NGUARD*K3D+1
  ex = size(1)-NGUARD
  ey = size(2)-NGUARD*K2D
  ez = size(3)-NGUARD*K3D

  select case(axis)
  case(IAXIS)
     area(sx,sy:ey,sz:ez) = hy_xarea(1,:,:,blockID) 
     area(ex,sy:ey,sz:ez) = hy_xarea(2,:,:,blockID) 
     dtdx(sx,sy:ey,sz:ez) = hy_xdtdx(1,:,:,blockID)
     dtdx(ex,sy:ey,sz:ez) = hy_xdtdx(2,:,:,blockID)
     grav(sx,sy:ey,sz:ez) = hy_xgrav(1,:,:,blockID) 
     grav(ex,sy:ey,sz:ez) = hy_xgrav(2,:,:,blockID)
     ngrav(sx,sy:ey,sz:ez) = hy_xngrav(1,:,:,blockID)
     ngrav(ex,sy:ey,sz:ez) = hy_xngrav(2,:,:,blockID) 
     fict(sx,sy:ey,sz:ez) = hy_xfict(1,:,:,blockID) 
     fict(ex,sy:ey,sz:ez) = hy_xfict(2,:,:,blockID)
     areaLeft(sx,sy:ey,sz:ez)   = hy_xareaAtFaces(1,1,:,:,blockID) 
     areaLeft(sx+1,sy:ey,sz:ez) = hy_xareaAtFaces(2,1,:,:,blockID) 
     areaLeft(ex,sy:ey,sz:ez)   = hy_xareaAtFaces(1,2,:,:,blockID) 
     areaLeft(ex+1,sy:ey,sz:ez) = hy_xareaAtFaces(2,2,:,:,blockID) 

     
  case(JAXIS)

     area(sx:ex,sy,sz:ez) = hy_yarea(:,1,:,blockID) 
     area(sx:ex,ey,sz:ez) =hy_yarea(:,2,:,blockID) 
     dtdx(sx:ex,sy,sz:ez) = hy_ydtdy(:,1,:,blockID) 
     dtdx(sx:ex,ey,sz:ez) = hy_ydtdy(:,2,:,blockID) 
     grav(sx:ex,sy,sz:ez) = hy_ygrav(:,1,:,blockID) 
     grav(sx:ex,ey,sz:ez) = hy_ygrav(:,2,:,blockID) 
     ngrav(sx:ex,sy,sz:ez) = hy_yngrav(:,1,:,blockID)
     ngrav(sx:ex,ey,sz:ez) = hy_yngrav(:,2,:,blockID)
     fict(sx:ex,sy,sz:ez) = hy_yfict(:,1,:,blockID) 
     fict(sx:ex,ey,sz:ez) = hy_yfict(:,2,:,blockID)
     areaLeft(sx:ex,sy,sz:ez)   = hy_yareaAtFaces(1,:,1,:,blockID) 
     areaLeft(sx:ex,sy+1,sz:ez) = hy_yareaAtFaces(2,:,1,:,blockID) 
     areaLeft(sx:ex,ey,sz:ez)   = hy_yareaAtFaces(1,:,2,:,blockID) 
     areaLeft(sx:ex,ey+1,sz:ez) = hy_yareaAtFaces(2,:,2,:,blockID) 
     
  case(KAXIS)
     area(sx:ex,sy:ey,sz) = hy_zarea(:,:,1,blockID)
     area(sx:ex,sy:ey,ez) = hy_zarea(:,:,2,blockID) 
     dtdx(sx:ex,sy:ey,sz) = hy_zdtdz(:,:,1,blockID) 
     dtdx(sx:ex,sy:ey,ez) = hy_zdtdz(:,:,2,blockID)
     grav(sx:ex,sy:ey,sz) = hy_zgrav(:,:,1,blockID)
     grav(sx:ex,sy:ey,ez) = hy_zgrav(:,:,2,blockID)
     ngrav(sx:ex,sy:ey,sz) = hy_zngrav(:,:,1,blockID) 
     ngrav(sx:ex,sy:ey,ez) = hy_zngrav(:,:,2,blockID) 
     fict(sx:ex,sy:ey,sz) = hy_zfict(:,:,1,blockID) 
     fict(sx:ex,sy:ey,ez) = hy_zfict(:,:,2,blockID) 
     areaLeft(sx:ex,sy:ey,sz)   = hy_zareaAtFaces(1,:,:,1,blockID)
     areaLeft(sx:ex,sy:ey,sz+1) = hy_zareaAtFaces(2,:,:,1,blockID)
     areaLeft(sx:ex,sy:ey,ez)   = hy_zareaAtFaces(1,:,:,2,blockID) 
     areaLeft(sx:ex,sy:ey,ez+1) = hy_zareaAtFaces(2,:,:,2,blockID) 
  end select

  return
end subroutine hy_ppm_getTemporaryData





