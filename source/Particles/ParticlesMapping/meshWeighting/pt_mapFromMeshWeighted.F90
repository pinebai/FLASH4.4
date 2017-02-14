!!****if* source/Particles/ParticlesMapping/meshWeighting/pt_mapFromMeshWeighted
!!
!! NAME
!!
!!  pt_mapFromMeshWeighted
!!
!! SYNOPSIS
!!
!!  pt_mapFromMeshWeighted(integer, INTENT(in)    :: numAttrib,
!!                        integer, INTENT(in)    :: attrib(2,numAttrib)
!!                           real, INTENT(in)    :: pos(MDIM)
!!                           real, INTENT(in)    :: bndBox(2,MDIM)
!!                           real, INTENT(in)    :: delta(MDIM)
!!                           real, pointer       :: solnVec(:,:,:,:)
!!                           real, INTENT(OUT)   :: partAttribVec(numAttrib)
!!
!! DESCRIPTION
!!
!!  Routine to map a list of  quantities defined on the mesh onto a single
!!  particle position.
!!
!!  Currently volume elements are assumed to be Cartesian or
!!  2D axisymmetric (r-z).
!!
!!  This version does quadratic interpolation onto particle
!!  positions (rather than the standard linear particle-mesh
!!  mapping).
!!  The routines assumes that the input arguement attrib contains the 
!!  list of particle attributes that need updating, and the corresponding
!!  mesh variables. The following constants are define in "Particles.h"
!!  which must be included in the routine
!!
!!     PART_DS_IND  : index into particle data structure   
!!     GRID_DS_IND  : index into the grid data structure
!!  As and example, if it is desired to map
!!  all three directional velocities, and the density then
!!  numAttrib is 4, and the array attrib will have the following values
!!     attrib(PART_DS_IND,1)=VELX_PART_PROP, attrib(GRID_DS_IND,1)=VELX_VAR
!!     attrib(PART_DS_IND,2)=VELY_PART_PROP, attrib(GRID_DS_IND,2)=VELY_VAR
!!     attrib(PART_DS_IND,3)=VELZ_PART_PROP, attrib(GRID_DS_IND,3)=VELZ_VAR
!!     attrib(PART_DS_IND,4)=DENS_PART_PROP, attrib(GRID_DS_IND,4)=DENS_VAR
!!  Here the order of the attributes is completely unimportant, only the 
!!  map between the particle property and the grid variable should be on the
!!  same index
!!
!! ARGUMENTS
!!
!!   numAttrib : number of attributes to update
!!   attrib    : list containing the attributes and the corresponding grid
!!               data structure index
!!   pos       : The physical coordinates of the particle
!!   bndBox    : the bounding box of the block containing the particle
!!   delta     : the dx,dy,dz of the block
!!   solnVec   : Pointer to the solution data block in which particle is found
!!   partAttribVec  : calculated attribute values
!!
!!
!!***

!=======================================================================
!!REORDER(4):solnVec


subroutine pt_mapFromMeshWeighted (numAttrib, attrib, pos, bndBox,&
     deltaCell,solnVec, partAttribVec)
  
  use Particles_data, ONLY : pt_geometry, pt_str_geometry
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: numAttrib
  integer, dimension(2, numAttrib),intent(IN) :: attrib
  real,dimension(MDIM), INTENT(in)    :: pos,deltaCell
  real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox
  real, pointer       :: solnVec(:,:,:,:)
  real,dimension(numAttrib), intent(OUT) :: partAttribVec

  real                :: dx_block_i, dy_block_i, dz_block_i, xp, yp, zp
  real                :: hx, hy, hz, result
  integer             :: ip, jp, kp, ivar,i, iface,jface,kface



  real,dimension(MDIM) :: h
  real,dimension(LEFT_EDGE:RIGHT_EDGE):: wx,wy,wz
  real,dimension(LEFT_EDGE:RIGHT_EDGE,LEFT_EDGE:RIGHT_EDGE,LEFT_EDGE:RIGHT_EDGE) :: wt
  integer, parameter :: L=LEFT_EDGE, C=CENTER, R=RIGHT_EDGE

!-------------------------------------------------------------------------------

! Get a pointer to the grid data for the specified block, and obtain block
! size and coordinate information.


! x-coord
! John ZuHone correctly pointed out that dx_block_i is the INVERSE of dx
!original  dx_block_i = NXB / size(1)  !! dx of this block
  dx_block_i = 1.0/deltaCell(1)
  xp = (pos(IAXIS) - bndBox(LOW,IAXIS)) * dx_block_i  !! offset of particle from block edge
  ip = floor(xp) + 1 + NGUARD       !! actual cell index (including guards)
  h(IAXIS) = modulo(xp,1.) - 5.e-1        !! remainder of ??
     
! y-coord
  if (NDIM >= 2) then
    dy_block_i = 1.0/deltaCell(2)
    yp = (pos(JAXIS) - bndBox(LOW,JAXIS)) * dy_block_i
    jp = floor(yp) + 1 + NGUARD
    h(JAXIS) = modulo(yp,1.) - 5.e-1
  else
    jp = 1
    h(JAXIS) = 0.
  endif

! z-coord
  if (NDIM == 3) then
    dz_block_i = 1.0/deltaCell(3)
    zp = (pos(KAXIS) - bndBox(LOW,KAXIS)) * dz_block_i
    kp = floor(zp) + 1 + NGUARD
    h(KAXIS) = modulo(zp,1.) - 5.e-1
  else
    kp = 1
    h(KAXIS) = 0.
  endif


! Pass .false. temporarily to use the traditional CIC scheme.
  call pt_assignWeights (.false., h, wx, wy, wz)

  do kface=LEFT_EDGE,RIGHT_EDGE
     do jface=LEFT_EDGE,RIGHT_EDGE
        do iface=LEFT_EDGE,RIGHT_EDGE
           wt(iface,jface,kface)=wx(iface)*wy(jface)*wz(kface)
        end do
     end do
  end do

  do i=1,numAttrib
     ivar=attrib(GRID_DS_IND,i)
  
! Compute interpolation coefficients of some sort


     result =  solnVec(ivar,ip  ,jp  ,kp  ) * wt(C,C,C)  + & ! center
               solnVec(ivar,ip+1,jp  ,kp  ) * wt(R,C,C)  + & ! right
               solnVec(ivar,ip-1,jp  ,kp  ) * wt(L,C,C)      ! left  

     if (NDIM >= 2) then
        result  =  result + &
               solnVec(ivar,ip-1,jp-1,kp  ) * wt(L,L,C)  + & 
               solnVec(ivar,ip  ,jp-1,kp  ) * wt(C,L,C)  + & 
               solnVec(ivar,ip+1,jp-1,kp  ) * wt(R,L,C)  + & 
               solnVec(ivar,ip+1,jp+1,kp  ) * wt(R,R,C)  + & 
               solnVec(ivar,ip  ,jp+1,kp  ) * wt(C,R,C)  + & 
               solnVec(ivar,ip-1,jp+1,kp  ) * wt(L,R,C)
     endif

     if (NDIM == 3) then
        result  =  result + &
               solnVec(ivar,ip  ,jp  ,kp-1) * wt(C,C,L)  + & 
               solnVec(ivar,ip-1,jp  ,kp-1) * wt(L,C,L)  + &  
               solnVec(ivar,ip-1,jp-1,kp-1) * wt(L,L,L)  + & 
               solnVec(ivar,ip  ,jp-1,kp-1) * wt(C,L,L)  + & 
               solnVec(ivar,ip+1,jp-1,kp-1) * wt(R,L,L)  + & 
               solnVec(ivar,ip+1,jp  ,kp-1) * wt(R,C,L)  + & 
               solnVec(ivar,ip+1,jp+1,kp-1) * wt(R,R,L)  + & 
               solnVec(ivar,ip  ,jp+1,kp-1) * wt(C,R,L)  + & 
               solnVec(ivar,ip-1,jp+1,kp-1) * wt(L,R,L)  + &  
                                                  
               solnVec(ivar,ip  ,jp  ,kp+1) * wt(C,C,R)  + & 
               solnVec(ivar,ip-1,jp  ,kp+1) * wt(L,C,R)  + & 
               solnVec(ivar,ip-1,jp-1,kp+1) * wt(L,L,R)  + & 
               solnVec(ivar,ip  ,jp-1,kp+1) * wt(C,L,R)  + & 
               solnVec(ivar,ip+1,jp-1,kp+1) * wt(R,L,R)  + & 
               solnVec(ivar,ip+1,jp  ,kp+1) * wt(R,C,R)  + & 
               solnVec(ivar,ip+1,jp+1,kp+1) * wt(R,R,R)  + & 
               solnVec(ivar,ip  ,jp+1,kp+1) * wt(C,R,R)  + & 
               solnVec(ivar,ip-1,jp+1,kp+1) * wt(L,R,R) 
     endif
     
     partAttribVec(i) = result
  end do

  return

end subroutine pt_mapFromMeshWeighted

