!!****f* source/Particles/Particles_mapFromMesh
!!
!! NAME
!!
!!  Particles_mapFromMesh
!!
!! SYNOPSIS
!!
!!  Particles_mapFromMesh(integer, INTENT(in)    :: mapType,
!!                        integer, INTENT(in)    :: numAttrib,
!!                        integer, INTENT(in)    :: attrib(2,numAttrib)
!!                           real, INTENT(in)    :: pos(MDIM)
!!                           real, INTENT(in)    :: bndBox(2,MDIM)
!!                           real, INTENT(in)    :: delta(MDIM)
!!                           real, pointer       :: solnVec(:,:,:,:)
!!                           real, INTENT(INOUT) :: partAttribVec(numAttrib)
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
!!  mesh variables. The following constants are defined in "Particles.h"
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
!!   mapType   : method used to map grid quantities to particles
!!   numAttrib : number of attributes to update
!!   attrib    : list containing the attributes and the corresponding grid
!!               data structure index
!!   pos       : The physical coordinates of the particle
!!   bndBox    : the bounding box of the block containing the particle
!!   delta     : the dx,dy,dz of the block
!!   solnVec   : Pointer to the solution data block in which particle is found
!!   partAttribVec  : calculated particle attribute values 
!!
!!
!!
!!
!!***

subroutine Particles_mapFromMesh (mapType,numAttrib, attrib, pos, bndBox,&
     deltaCell,solnVec, partAttribVec)
  

  implicit none

#include "constants.h"

  integer, INTENT(in) :: mapType,numAttrib
  integer, dimension(2, numAttrib), INTENT(in) :: attrib
  real,dimension(MDIM), INTENT(in)    :: pos,deltaCell
  real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox
  real, pointer       :: solnVec(:,:,:,:)
  real,dimension(numAttrib), intent(OUT) :: partAttribVec
  partAttribVec = 0.0
  return
end subroutine Particles_mapFromMesh

!===============================================================================

