!!****f* source/Grid/Grid_mapMeshToParticles
!!
!! NAME
!!
!!  Grid_mapMeshToParticles
!!
!! SYNOPSIS
!!
!!  Grid_mapMeshToParticles(real(inout) :: particles(part_props,numParticles),
!!                          integer(in) :: part_props,
!!                          integer(in) :: part_blkID,
!!                          integer(in) :: numParticles,
!!                          integer(in) :: posAttrib(POS),
!!                          integer(in) :: numAttrib,
!!                          integer(in) :: attrib(:,numAttrib),
!!                          integer(in) :: mapType,
!!               optional,  integer(in) :: gridDataStruct)
!!
!! DESCRIPTION
!!
!!  Routine to map a quantity defined on the mesh onto the particle positions.
!!
!! ARGUMENTS
!!
!!     particles:     Data structure containing particles information
!!     part_props : number of particle attributes
!!     part_blkID : the index of particle attributes that carries the block number
!!     numParticles : the number of particles on my proc
!!     posAttrib           : particles data structure indices pointing
!!                           to particle positions
!!     numAttrib           : number of attributes that need to be mapped
!!     attrib              : list of attributes and their corresponding
!!                           mesh data structure indices
!!                           processor
!!     mapType : method for mapping grid quantities to particles
!!     gridDataStruct: The Grid data structure that varGrid refers to; one of
!!                    CENTER, FACEX, FACEY, FACEZ, GRIDVAR.  
!!                    If this argument is not present, the default is CENTER, 
!!                    so that attrib indicates a list of variable in UNK.
!!
!!  
!!***
!*******************************************************************************

subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                    numParticles,posAttrib,&
                                    numAttrib, attrib,&
                                    mapType,gridDataStruct)

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Particles.h"

  integer, INTENT(in) :: part_props, numParticles, numAttrib, part_blkID
  real, INTENT(inout),dimension(part_props,numParticles) :: particles
  integer,dimension(MDIM), intent(IN) :: posAttrib
  integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
  integer, INTENT(IN) :: mapType
  integer, optional, intent(IN) :: gridDataStruct

  return
end subroutine Grid_mapMeshToParticles

