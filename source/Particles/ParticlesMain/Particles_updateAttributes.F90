!!****if* source/Particles/ParticlesMain/Particles_updateAttributes
!!
!! NAME
!!
!!  Particles_updateAttributes
!!
!! SYNOPSIS
!!
!!  Particles_updateAttributes()
!!
!! DESCRIPTION
!!
!!  Particle attribute advancement routines.  Use this routine to 
!!  update user-define particle attributes beyond the usual
!!  particle attributes of position, velocity, block, tag, and mass.
!!  It is usually called indirectly from IO_output. It makes sure that guard
!!  cells are current and, if necessary, that Eos has been applied to them
!!  before particles attributes are calculated by interpolation from the grid.
!!   
!!  The particle attributes to map from the grid can be specified at runtime
!!  using runtime parameters particle_attribute_1, particle_attribute_2,
!!  etc.
!!
!! ARGUMENTS
!!  
!! NOTES
!!
!! The map between particle property and a mesh variable on which the
!! property is dependent has to be specified in the Config file of 
!! the Simulation directory for this routine to work right. Please
!! see the Config file of IsentropicVortex setup for an example, and
!! also see the Setup chapter of the User's Guide.
!! 
!!
!!
!!***


subroutine Particles_updateAttributes()

  use Particles_data, ONLY : particles,pt_numLocal,pt_maxPerProc,&
       pt_gcMaskForWrite,pt_gcMaskSizeForWrite,pt_attributes,&
       pt_meshVar,pt_numAttributes, pt_meshMe, useParticles, pt_posAttrib,&
       pt_typeInfo, pt_numLost, pt_keepLostParticles, pt_reduceGcellFills
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_fillGuardCells,Grid_mapMeshToParticles,&
       Grid_sortParticles
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer,parameter :: GRID_STRUCT=1,PARTICLE_MAP=2
  integer :: i,j,k,mapType
  integer, dimension(PART_ATTR_DS_SIZE,pt_numAttributes) :: attrib
  integer :: part_props, numAttrib
  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  integer       ::p_begin,p_end, p_count,numGridDataStruct, pbak, pfor, lostNow 

  integer,dimension(GRID_STRUCT:PARTICLE_MAP,MAX_GRID_DATA_STRUCT) :: gridStruct
  integer :: savePtNumLocal
  !-------------------------------------------------------------------------

  if (.NOT. useParticles) return
  if (pt_numAttributes==0) return

  savePtNumLocal=pt_numLocal
  if(pt_keepLostParticles) then
     !! if the last particle has a block number of LOST, then the lost particles
     !! are included in the count. the count needs to adjust temporarily to avoid
     !! mapping particles that have left the domain
     if(particles(BLK_PART_PROP,pt_numLocal)==LOST)then
        pt_numLocal=pt_numLocal-pt_numLost
        if(particles(BLK_PART_PROP,pt_numLocal+1)/=LOST)then
           print*,'the value is ', particles(BLK_PART_PROP,pt_numLocal+1),pt_numLocal,pt_numLost,savePtNumLocal
           call Driver_abortFlash("at update of particles counting is wrong")
        end if
        if(particles(BLK_PART_PROP,pt_numLocal)==LOST)call &
             Driver_abortFlash("at update too many particles lost")
     end if
  end if

  part_props=NPART_PROPS
  call Timers_start ("updateAttributes map from mesh")

  !! This section of the code fills the data strucure "gridStruct"
  !! This data structure is a 2D array that stores the constants
  !! representing the gridDataStucture and its corresponding 
  !! particle_map data structure needed by the setup script generated
  !! routine that finds out which grid variable in which data structure
  !! corresponds to the given particle attribute. This information is
  !! stored for every grid data structure included in the simulation.

  numGridDataStruct=0

  if(NUNK_VARS > 0) then
     numGridDataStruct=numGridDataStruct+1
     gridStruct(GRID_STRUCT,numGridDataStruct)=CENTER
     gridStruct(PARTICLE_MAP,numGridDataStruct)=PARTICLEMAP_UNK
  end if

  if(NSCRATCH_GRID_VARS>0) then
     numGridDataStruct=numGridDataStruct+1  
     gridStruct(GRID_STRUCT,numGridDataStruct)=SCRATCH
     gridStruct(PARTICLE_MAP,numGridDataStruct)=PARTICLEMAP_SCRATCH
  end if

  if(NFACE_VARS>0) then
     numGridDataStruct=numGridDataStruct+1  
     gridStruct(GRID_STRUCT,numGridDataStruct)=FACEX
     gridStruct(PARTICLE_MAP,numGridDataStruct)=PARTICLEMAP_FACEX
     if(NDIM>1) then
        numGridDataStruct=numGridDataStruct+1  
        gridStruct(GRID_STRUCT,numGridDataStruct)=FACEY
        gridStruct(PARTICLE_MAP,numGridDataStruct)=PARTICLEMAP_FACEY
     end if
     if(NDIM>2) then
        numGridDataStruct=numGridDataStruct+1  
        gridStruct(GRID_STRUCT,numGridDataStruct)=FACEZ
        gridStruct(PARTICLE_MAP,numGridDataStruct)=PARTICLEMAP_FACEZ
     end if
  end if

  !! If no valid data structures were found, something is very wrong
  if(numGridDataStruct==0)&
       call Driver_abortFlash("Particles_updateAttributes: no data structures")

  !! We reach here only if at least one valid data structure was found.
  !! The number of attributes to be output is greater than 0,
  !! so fill guardcell in preparation for updating the attributes.
  if (ANY(pt_gcMaskForWrite)) then
     if (pt_reduceGcellFills) then
        call Grid_fillGuardCells(CENTER_FACES,ALLDIR, doEos=.true.,&
             unitReadsMeshDataOnly=.true.)
     else
        call Grid_fillGuardCells(CENTER_FACES,ALLDIR, doEos=.true.,&
             maskSize=pt_gcMaskSizeForWrite,mask=pt_gcMaskForWrite,&
             makeMaskConsistent=.true.)
     end if
  end if

  !! Now sort the particles based upon their type and block numbers.
  !! The sort routine first separates various particles types, and then
  !! within each type, sorts the particles based on their associated
  !! block number

#ifdef TYPE_PART_PROP
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
  call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
       pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif

  if(pt_keepLostParticles) then
     pfor=pt_numLocal
     do while(particles(BLK_PART_PROP,pt_numLocal)==LOST)
        pt_numLocal=pt_numLocal-1
     end do
     lostNow=pfor-pt_numLocal
     pbak=pt_maxPerProc-pt_numLost
     if(pbak<pt_numLocal)call Driver_abortFlash("no more space for lost particles")
     particles(:,pbak-lostNow:pbak)=particles(:,pt_numLocal+1:pt_numLocal+lostNow)
     pt_numLost=pt_numLost+lostNow
  end if

  !! After running a sort, it is necessary to update the information about
  !! starting point and count for each particle type. If only one type is 
  !! included in the simulation, then starting point is always 1,
  !! and the count is pt_numLocal

  call pt_updateTypeDS(particlesPerBlk)
  
  !! Now loop over all included data structures to see if any attribute
  !! is associated with it, and if so, update the values of that attribute
  !! by mapping them from the corresponding Grid variable.
  do j=1,numGridDataStruct

     numAttrib=0
     !! This loop finds the number of attributes associated with the
     !! current grid data structure.
     do i = 1,pt_numAttributes
        if(pt_meshVar(PT_MAP,i)==gridStruct(PARTICLE_MAP,j)) then
           numAttrib=numAttrib+1
           attrib(PART_DS_IND,numAttrib)=pt_attributes(i)
           attrib(GRID_DS_IND,numAttrib)=pt_meshVar(PT_VAR,i)
        end if
     end do
     
     !! if at least one attribute was found to be associated with 
     !! the current data structure then update that attribute
     !! appropriately for each particle type. Here different
     !! particle types may be using different methods to get the
     !! value from the grid variable.
     if(numAttrib>0) then
        do i = 1,NPART_TYPES
           p_count=pt_typeInfo(PART_LOCAL,i)
           p_begin=pt_typeInfo(PART_TYPE_BEGIN,i)
           p_end=p_begin+p_count-1
           mapType=pt_typeInfo(PART_MAPMETHOD,i)
           call Grid_mapMeshToParticles(particles(:,p_begin:p_end),&
                part_props,BLK_PART_PROP,p_count,&
                pt_posAttrib,numAttrib,attrib,mapType,&
                gridDataStruct=gridStruct(GRID_STRUCT,j))
        end do
     end if
  end do
  call Timers_stop ("updateAttributes map from mesh")
  pt_numLocal=savePtNumLocal
  
  !-----------------------------------------------------------------------
  return
  
end subroutine Particles_updateAttributes



